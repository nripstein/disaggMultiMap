#' Construct design and projection matrices for prediction
#'
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
#' @param expected_cov_names Character vector of training covariate names (order matters).
#'   If length 0, predictions are intercept-only regardless of provided rasters.
#' @param time_varying_betas Logical; used for clearer error messages when aligning layers.
#' @keywords internal
get_predict_matrices <- function(data,
                                 new_data = NULL,
                                 expected_cov_names = NULL,
                                 time_varying_betas = FALSE) {
  #-- Validate inputs --
  stopifnot(inherits(data, "disag_data_mmap"))
  n_times <- length(data$time_points)

  #-- Build covariate list --
  # If the user didn’t supply new_data, use the original rasters.
  # Otherwise accept a single SpatRaster (recycled) or a list matching time points.
  training_rasters <- is.null(new_data)
  cov_list <- if (training_rasters) {
    data$covariate_rasters_list
  } else if (inherits(new_data, "SpatRaster")) {
    rep(list(new_data), n_times)
  } else if (is.list(new_data) && length(new_data) == n_times) {
    # ensure each element is a SpatRaster if provided
    if (length(new_data) && any(!vapply(new_data, inherits, logical(1), what = "SpatRaster"))) {
      stop("Each element of `new_data` must be a SpatRaster.")
    }
    new_data
  } else {
    stop("`new_data` must be NULL, a SpatRaster, or a list of length ", n_times, ".")
  }

  #-- Compute SPDE projection matrix & coords --
  # coords: data.frame with x,y for each prediction cell.
  coords <- data$coords_for_prediction
  # A matrix projects from mesh nodes to cells.
  Amatrix <- fmesher::fm_evaluator(data$mesh, loc = as.matrix(coords))$proj$A
  categorical_schema <- tryCatch(data$categorical_covariate_schema, error = function(e) NULL)

  #-- Construct design matrices for each time slice --
  X_list <- vector("list", n_times)

  # Expected covariate names (training order); default to character(0)
  if (is.null(expected_cov_names)) expected_cov_names <- character(0)

  for (i in seq_len(n_times)) {
    cov_i <- cov_list[[i]]
    if (!is.null(cov_i) && inherits(cov_i, "SpatRaster")) {
      cov_i <- encode_categorical_raster_stack(
        covariates = cov_i,
        categorical_schema = categorical_schema,
        time_index = i,
        context = if (isTRUE(training_rasters)) {
          "predict.disag_model_mmap_aghq() using training covariates"
        } else {
          "predict.disag_model_mmap_aghq() using new_data"
        }
      )
    }

    # Intercept-only model (no covariates in training): ignore provided covariates
    if (length(expected_cov_names) == 0L || is.null(cov_i)) {
      # No covariates → intercept-only design: column of 1’s
      X_list[[i]] <- matrix(
        1,
        nrow = nrow(Amatrix),
        ncol = 1,
        dimnames = list(NULL, "Intercept")
      )
    } else {
      # Align covariate layer order to training order
      cur_names <- names(cov_i)
      if (!identical(cur_names, expected_cov_names)) {
        if (!setequal(cur_names, expected_cov_names)) {
          stop(paste0(
            "Covariate layer names for prediction do not match training.\n",
            if (isTRUE(time_varying_betas)) "When time_varying_betas=TRUE, names and order must match across time.\n" else "",
            "Expected: ", paste(expected_cov_names, collapse = ", "), "\n",
            "Got     : ", paste(cur_names, collapse = ", ")
          ))
        }
        # Same set, different order -> reorder to training order
        cov_i <- cov_i[[expected_cov_names]]
      }

      # Extract raster values as matrix: rows=cells, cols=variables (training order)
      X_cov <- terra::values(cov_i, mat = TRUE)
      colnames(X_cov) <- names(cov_i)
      # Prepend intercept column
      X_list[[i]] <- cbind(Intercept = 1, X_cov)
    }

    # Ensure design & projection dims align
    if (nrow(X_list[[i]]) != nrow(Amatrix)) {
      stop(
        "Dim mismatch at time ", i,
        ": nrow(X) = ", nrow(X_list[[i]]),
        " vs nrow(A) = ", nrow(Amatrix), "."
      )
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
                                          new_data = NULL,
                                          predict_iid = FALSE,
                                          N = 1e3,
                                          CI = 0.95,
                                          verbose = FALSE,
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
  mats <- get_predict_matrices(
    data = object$data,
    new_data = new_data,
    expected_cov_names = tryCatch(object$model_setup$coef_meta$cov_names, error = function(e) NULL),
    time_varying_betas = isTRUE(object$model_setup$time_varying_betas)
  )
  X_list <- mats$X_list # list of design matrices per time
  Amatrix <- mats$A # projection matrix
  coords <- mats$coords # coords for raster building

  #-- Draw from posterior marginal via AGHQ (theta only; no random effects) --
  samps <- aghq::sample_marginal(object$aghq_model, N)
  W <- samps$samps # matrix: (n_theta × N_draws)
  # Canonicalize draw names to match training normalization
  if (is.null(rownames(W))) stop("AGHQ draws lack row names; cannot align to model parameters.")
  coef_meta <- tryCatch(object$model_setup$coef_meta, error = function(e) NULL)
  tv_flag   <- isTRUE(object$model_setup$time_varying_betas)
  rownames(W) <- canonicalize_draw_names(rownames(W), coef_meta, tv_flag)
  theta_order <- object$model_setup$theta_order
  if (is.null(theta_order)) {
    stop("AGHQ prediction requires parameter name order (theta_order). Please refit the model with a recent version.")
  }

  beta_map <- object$model_setup$beta_index_map
  if (is.null(beta_map)) {
    stop("AGHQ prediction requires beta_index_map. Please refit the model with a recent version.")
  }
  tv   <- isTRUE(beta_map$tv)
  p    <- as.integer(beta_map$p)
  Tn   <- as.integer(beta_map$Tn)

  n_times <- length(X_list)
  layer_names <- paste0("time_", object$data$time_points)

  # Prepare lists to collect outputs
  mean_preds <- vector("list", n_times)
  mean_fields <- vector("list", n_times)
  mean_covs <- vector("list", n_times)
  ci_low <- vector("list", n_times)
  ci_high <- vector("list", n_times)

  #-- Setup link function --
  link_fn <- switch(object$model_setup$link,
                    log      = exp,
                    identity = identity,
                    logit    = function(x) 1 / (1 + exp(-x)),
                    stop("Unsupported link: ", object$model_setup$link)
  )

  #-- Compute field at conditional mode (theta mode) if available --
  field_vec <- NULL
  if (isTRUE(object$model_setup$field)) {
    # Try to extract nodemean at mode from AGHQ object first
    mode_sources <- list(
      tryCatch(object$aghq_model$optresults$mode, error = function(e) NULL),
      tryCatch(object$aghq_model$modesandhessians$mode, error = function(e) NULL),
      tryCatch(object$aghq_model$mode, error = function(e) NULL)
    )
    nodemean_from_mode <- NULL
    for (mv in mode_sources) {
      if (!is.null(mv)) {
        nm <- names(mv)
        if (!is.null(nm)) {
          idx <- grepl("^nodemean", nm)
          if (any(idx)) {
            nodemean_from_mode <- as.numeric(mv[idx])
            break
          }
        }
      }
    }

    if (!is.null(nodemean_from_mode) && length(nodemean_from_mode) == ncol(Amatrix)) {
      field_vec <- as.numeric((Amatrix %*% matrix(nodemean_from_mode, ncol = 1))[, 1])
    } else if (!is.null(object$obj)) {
      # Fallback to TMB objective's last parameters
      raw_pars <- tryCatch(object$obj$env$last.par.best, error = function(e) NULL)
      if (!is.null(raw_pars)) {
        pl <- tryCatch(slice_params_tmb(raw_pars, object), error = function(e) NULL)
        if (!is.null(pl) && length(pl$nodemean) == ncol(Amatrix)) {
          field_vec <- as.numeric((Amatrix %*% matrix(pl$nodemean, ncol = 1))[, 1])
        }
      }
    }
  }

  #-- Loop over time points --
  for (i in seq_len(n_times)) {
    # 1. Extract design matrix
    X <- X_list[[i]]

    # 2. Subset posterior draws for betas (aligned to training order)
    if (!tv) {
      rows_t <- c(beta_map$intercept_idx, if (p > 0L) beta_map$slope_idx)
    } else {
      # Guard on time index sizes
      if (n_times != Tn) stop("Prediction time points do not match training time points.")
      if (p > 0L) {
        rows_t <- c(beta_map$intercept_idx[i], beta_map$slope_idx[, i])
      } else {
        rows_t <- beta_map$intercept_idx[i]
      }
    }
    # Map required theta names to rows of W (only for betas we need)
    theta_req <- theta_order[rows_t]
    W_rows <- match(theta_req, rownames(W))
    if (anyNA(W_rows)) {
      # Fallback: if draw rows follow theta_order positionally, use positions
      idx_in_theta <- match(theta_req, theta_order)
      if (!anyNA(idx_in_theta) && nrow(W) >= max(idx_in_theta)) {
        W_rows <- idx_in_theta
      } else {
        miss <- theta_req[is.na(W_rows)]
        stop("AGHQ draws are missing required beta parameters: ", paste(miss, collapse = ", "))
      }
    }
    W_beta <- W[W_rows, , drop = FALSE] # ((1+p) × N)

    # 3. Compute linear predictors (covariates + optional field at mode)
    lin_cov <- X %*% W_beta # (n_cells × N)
    if (!is.null(field_vec)) {
      lin_cov <- sweep(lin_cov, 1, field_vec, `+`)
    }

    # 4. Apply link to get lambda draws
    lam_mat <- link_fn(lin_cov)

    # 5. Summarize posterior mean
    mean_vals <- Matrix::rowMeans(lam_mat)

    # 6. Build SpatRasters for mean & components (field from mode if available)
    mean_preds[[i]] <- terra::rast(cbind(coords, y = mean_vals), type = "xyz")
    # covariate component on linear scale (mean across draws before link)
    cov_lin_mean <- Matrix::rowMeans(X %*% W_beta)
    mean_covs[[i]]  <- terra::rast(cbind(coords, y = cov_lin_mean), type = "xyz")
    if (!is.null(field_vec)) {
      mean_fields[[i]] <- terra::rast(cbind(coords, y = field_vec), type = "xyz")
    } else {
      mean_fields[[i]] <- NULL
    }

    # 7. Compute cellwise CIs
    probs <- c((1 - CI) / 2, 1 - (1 - CI) / 2)
    ci_mat <- apply(lam_mat, 1, stats::quantile, probs = probs, na.rm = TRUE)
    ci_low[[i]] <- terra::rast(cbind(coords, y = ci_mat[1, ]), type = "xyz")
    ci_high[[i]] <- terra::rast(cbind(coords, y = ci_mat[2, ]), type = "xyz")
  }

  #-- Name layers consistently --
  names(mean_preds) <- layer_names
  # Field layers may be identical across time; if present, replicate names
  if (!all(vapply(mean_fields, is.null, logical(1)))) names(mean_fields) <- layer_names
  names(mean_covs) <- layer_names
  names(ci_low) <- layer_names
  names(ci_high) <- layer_names

  #-- Assemble final output --
  pred_stack <- do.call(c, mean_preds)
  cov_stack  <- do.call(c, mean_covs)
  # Stack field if available
  field_stack <- if (!all(vapply(mean_fields, is.null, logical(1)))) do.call(c, mean_fields) else NULL
  mean_prediction <- list(
    prediction = pred_stack,
    field      = field_stack,
    iid        = NULL,
    covariates = cov_stack
  )

  uncertainty_prediction <- list(
    realisations = NULL,
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
    msg <- sprintf("predict.disag_model_mmap_aghq() runtime: %.2f mins", as.numeric(elapsed))
    if (isTRUE(object$model_setup$field) && is.null(field_vec)) {
      msg <- paste0(msg, "; note: field contribution omitted (no mode available)")
    } else if (isTRUE(object$model_setup$field)) {
      msg <- paste0(msg, "; field added at mode; uncertainty excludes field")
    }
    message(msg)
  }

  invisible(out)
}
