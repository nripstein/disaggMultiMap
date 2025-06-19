#' @description
#' Internal helper for 'predict.disag_model_mmap_aghq()'.
#' Builds per-time design matrices (with intercept), the SPDE projection matrix,
#' and the coordinate table for raster reconstruction.
#'
#' @param data A 'disag_data_mmap' object (from 'prepare_data_mmap()').
#' @param new_data Optional new covariate data:
#'   - a single 'SpatRaster' (recycled across all times), or
#'   - a list of length 'length(data$time_points)' of 'SpatRaster' objects.
#' @return A list with elements:
#'   - 'X_list': list of design matrices (each n_cells × p, with "Intercept").
#'   - 'A': SPDE projection matrix (n_cells × n_knots).
#'   - 'coords': data.frame of x/y coordinates for each cell.
#' @keywords internal
get_predict_matrices <- function(data, new_data = NULL) {
  #-- Validate inputs --
  stopifnot(inherits(data, "disag_data_mmap"))
  n_times <- length(data$time_points)

  #-- Build covariate list --
  # If the user didn’t supply new_data, use the original rasters.
  # Otherwise accept a single SpatRaster (recycled) or a list matching time points.
  cov_list <- if (is.null(new_data)) {
    data$covariate_rasters_list
  } else if (inherits(new_data, "SpatRaster")) {
    rep(list(new_data), n_times)
  } else if (is.list(new_data) && length(new_data) == n_times) {
    new_data
  } else {
    stop("`new_data` must be NULL, a SpatRaster, or a list of length ", n_times, ".")
  }

  #-- Compute SPDE projection matrix & coords --
  # coords: data.frame with x,y for each prediction cell.
  coords  <- data$coords_for_prediction
  # A matrix projects from mesh nodes to cells.
  Amatrix <- fmesher::fm_evaluator(data$mesh, loc = as.matrix(coords))$proj$A

  #-- Construct design matrices for each time slice --
  X_list <- vector("list", n_times)

  for (i in seq_len(n_times)) {
    cov_i <- cov_list[[i]]

    if (is.null(cov_i)) {
      # No covariates → intercept-only design: column of 1’s
      X_list[[i]] <- matrix(
        1,
        nrow = nrow(Amatrix),
        ncol = 1,
        dimnames = list(NULL, "Intercept")
      )
    } else {
      # Extract raster values as matrix: rows=cells, cols=variables
      X_cov <- terra::values(cov_i, mat = TRUE)
      colnames(X_cov) <- names(cov_list[[i]])
      # Prepend intercept column
      X_list[[i]] <- cbind(Intercept = 1, X_cov)
    }

    # Ensure design & projection dims align
    if (nrow(X_list[[i]]) != nrow(Amatrix)) {
      stop("Dim mismatch at time ", i,
           ": nrow(X) = ", nrow(X_list[[i]]),
           " vs nrow(A) = ", nrow(Amatrix), ".")
    }
  }

  list(
    X_list = X_list,
    A      = Amatrix,
    coords = coords
  )
}

#' Predict mean & credible intervals for AGHQ-fitted disaggregation model
#'
#' @description
#' Given a 'disag_model_mmap_aghq' object, draws from the AGHQ marginal, builds
#' per-cell posterior samples, and returns means and credible-interval rasters.
#'
#' @param object A 'disag_model_mmap_aghq' fit (from 'disag_model_mmap_aghq()').
#' @param new_data Optional covariates for prediction (see helper).
#' @param predict_iid Currently not implemented; must be FALSE.
#' @param N Number of marginal draws to sample (default 1000).
#' @param CI Credible-interval level in (0,1) (default 0.95).
#' @param verbose If TRUE, prints runtime in minutes.
#' @param ... Unused.
#'
#' @return An object of class 'disag_prediction_mmap_aghq' containing:
#'   - 'mean_prediction': list of SpatRasters ('prediction', 'field', 'covariates').
#'   - 'uncertainty_prediction': list with 'predictions_ci$lower' & 'upper'.
#' @method predict disag_model_mmap_aghq
#' @export
predict.disag_model_mmap_aghq <- function(object,
                                          new_data    = NULL,
                                          predict_iid = FALSE,
                                          N           = 1e3,
                                          CI          = 0.95,
                                          verbose     = FALSE,
                                          ...) {
  #-- Input validation --
  stopifnot(inherits(object, "disag_model_mmap_aghq"))
  if (predict_iid) stop("`predict_iid = TRUE` is not yet supported.")
  if (!is.numeric(N) || length(N) != 1 || N < 1) {
    stop("`N` must be a single positive integer.")
  }
  if (!is.numeric(CI) || CI <= 0 || CI >= 1) {
    stop("`CI` must be a number strictly between 0 and 1.")
  }

  start_time <- Sys.time()

  #-- Extract design & projection matrices --
  mats    <- get_predict_matrices(object$data, new_data)
  X_list  <- mats$X_list      # list of design matrices per time
  Amatrix <- mats$A           # projection matrix
  coords  <- mats$coords      # coords for raster building

  #-- Draw from posterior marginal via AGHQ --
  samps  <- aghq::sample_marginal(object$aghq_model, N)
  W      <- samps$samps       # matrix: (n_params × N_draws)
  nd_idx <- which(rownames(W) == "nodemean")

  n_times     <- length(X_list)
  layer_names <- paste0("time_", object$data$time_points)

  # Prepare lists to collect outputs
  mean_preds  <- vector("list", n_times)
  mean_fields <- vector("list", n_times)
  mean_covs   <- vector("list", n_times)
  ci_low      <- vector("list", n_times)
  ci_high     <- vector("list", n_times)

  #-- Setup link function --
  link_fn <- switch(
    object$model_setup$link,
    log      = exp,
    identity = identity,
    logit    = function(x) 1/(1 + exp(-x)),
    stop("Unsupported link: ", object$model_setup$link)
  )

  #-- Loop over time points --
  for (i in seq_len(n_times)) {
    # 1. Extract design matrix
    X        <- X_list[[i]]
    ncov     <- ncol(X) - 1
    beta_idx <- seq_len(ncov + 1)

    # 2. Subset posterior draws for betas & field nodes
    W_beta   <- W[beta_idx, , drop = FALSE]  # (p × N)
    field_w  <- W[nd_idx, , drop = FALSE]  # (1 × N)

    # 3. Compute linear predictors
    lin_cov  <- X %*% W_beta  # (n_cells × N)
    field_mat<- Amatrix %*% field_w   # (n_cells × N)
    lin_all  <- lin_cov + field_mat

    # 4. Apply link to get lambda draws
    lam_mat  <- link_fn(lin_all)

    # 5. Summarize posterior mean
    mean_vals  <- rowMeans(lam_mat)
    field_vals <- rowMeans(field_mat)

    # 6. Build SpatRasters for mean & components
    mean_preds[[i]]  <- terra::rast(cbind(coords, y = mean_vals),  type = "xyz")
    mean_fields[[i]] <- terra::rast(cbind(coords, y = field_vals), type = "xyz")
    mean_covs[[i]]   <- mean_preds[[i]] - mean_fields[[i]]

    # 7. Compute cellwise CIs
    probs     <- c((1 - CI)/2, 1 - (1 - CI)/2)
    ci_mat    <- apply(lam_mat, 1, quantile, probs = probs, na.rm = TRUE)
    ci_low[[i]]  <- terra::rast(cbind(coords, y = ci_mat[1, ]), type = "xyz")
    ci_high[[i]] <- terra::rast(cbind(coords, y = ci_mat[2, ]), type = "xyz")
  }

  #-- Name layers consistently --
  names(mean_preds) <- layer_names
  names(mean_fields) <- layer_names
  names(mean_covs) <- layer_names
  names(ci_low) <- layer_names
  names(ci_high) <- layer_names

  #-- Assemble final output --
  mean_prediction <- list(
    prediction = do.call(c, mean_preds),
    field      = do.call(c, mean_fields),
    iid        = NULL,
    covariates = do.call(c, mean_covs)
  )

  uncertainty_prediction <- list(
    realisations   = NULL,
    predictions_ci = list(
      lower = do.call(c, ci_low),
      upper = do.call(c, ci_high)
    )
  )

  out <- list(
    mean_prediction        = mean_prediction,
    uncertainty_prediction = uncertainty_prediction
  )
  class(out) <- c("disag_prediction_mmap_aghq", "list")

  #-- Optionally report runtime --
  if (verbose) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    message(sprintf("predict.disag_model_mmap_aghq() runtime: %.2f mins", as.numeric(elapsed)))
  }

  invisible(out)
}
