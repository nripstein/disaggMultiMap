#' Predict for Multi-Map Disaggregation Model fit with MCMC
#'
#' @param object A fitted \code{disag_model_mmap_mcmc} object.
#' @param new_data Optionally, a new SpatRaster (or list of them) for prediction.
#' @param predict_iid Logical. If TRUE, include the polygon IID effect in predictions.
#' @param N Number of posterior draws to use for uncertainty estimation.
#' @param CI Credible interval level (default 0.95).
#' @param verbose Logical; if TRUE, print runtime diagnostics.
#' @param ... Further arguments.
#'
#' @export
#' @method predict disag_model_mmap_mcmc
predict.disag_model_mmap_mcmc <- function(object,
                                          new_data = NULL,
                                          predict_iid = FALSE,
                                          N = 100,
                                          CI = 0.95,
                                          verbose = FALSE,
                                          ...) {
  stopifnot(inherits(object, "disag_model_mmap_mcmc"))
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N < 1) {
    stop("`N` must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(CI) || length(CI) != 1L || CI <= 0 || CI >= 1) {
    stop("`CI` must be a number strictly between 0 and 1.", call. = FALSE)
  }

  start_time <- Sys.time()
  N <- as.integer(round(N))

  draws <- extract_mcmc_draws_matrix_mmap(object$mcmc_fit)
  n_available <- nrow(draws)
  if (!n_available) {
    stop("No posterior draws available from `mcmc_fit`.", call. = FALSE)
  }

  if (N > n_available) {
    warning(
      paste0(
        "`N` (", N, ") exceeds available posterior draws (", n_available,
        "). Using all available draws."
      ),
      call. = FALSE
    )
    N_use <- n_available
  } else {
    N_use <- N
  }

  draw_idx <- if (N_use < n_available) {
    sample.int(n_available, size = N_use, replace = FALSE)
  } else {
    seq_len(n_available)
  }
  draw_mat <- draws[draw_idx, , drop = FALSE]

  param_draws <- vector("list", nrow(draw_mat))
  for (rr in seq_len(nrow(draw_mat))) {
    draw <- draw_mat[rr, ]
    names(draw) <- colnames(draw_mat)
    param_draws[[rr]] <- draw_to_param_list_mcmc_mmap(draw, object)
  }

  times <- object$data$time_points
  n_times <- length(times)
  layer_names <- paste0("time_", times)

  mean_preds <- vector("list", n_times)
  mean_fields <- vector("list", n_times)
  mean_iids <- vector("list", n_times)
  mean_covs <- vector("list", n_times)
  realizations_list <- vector("list", n_times)
  ci_lower <- vector("list", n_times)
  ci_upper <- vector("list", n_times)
  probs <- c((1 - CI) / 2, 1 - (1 - CI) / 2)

  for (ii in seq_along(times)) {
    if (is.null(new_data)) {
      new_data_i <- NULL
      use_training <- TRUE
    } else {
      if (is.list(new_data) && length(new_data) == n_times) {
        new_data_i <- new_data[[ii]]
      } else {
        new_data_i <- new_data
      }
      use_training <- FALSE
    }

    objs <- setup_objects_mmap(
      model_output = object,
      new_data = new_data_i,
      predict_iid = predict_iid,
      time_index = ii,
      use_training = use_training
    )

    draw_preds <- vector("list", N_use)
    draw_covs <- vector("list", N_use)
    draw_fields <- vector("list", N_use)
    draw_iids <- vector("list", N_use)

    for (rr in seq_len(N_use)) {
      out <- predict_single_raster_mmap(
        model_parameters = param_draws[[rr]],
        objects = objs,
        link_function = object$model_setup$link
      )
      draw_preds[[rr]] <- out$prediction
      draw_covs[[rr]] <- out$covariates
      draw_fields[[rr]] <- out$field
      draw_iids[[rr]] <- out$iid
    }

    pred_stack <- do.call(c, draw_preds)
    realizations_list[[ii]] <- pred_stack
    mean_preds[[ii]] <- terra::app(pred_stack, mean, na.rm = TRUE)
    ci_stack <- terra::app(pred_stack, function(x) stats::quantile(x, probs = probs, na.rm = TRUE))
    ci_lower[[ii]] <- ci_stack[[1]]
    ci_upper[[ii]] <- ci_stack[[2]]

    cov_stack <- do.call(c, draw_covs)
    mean_covs[[ii]] <- terra::app(cov_stack, mean, na.rm = TRUE)

    if (!all(vapply(draw_fields, is.null, logical(1)))) {
      field_stack <- do.call(c, draw_fields)
      mean_fields[[ii]] <- terra::app(field_stack, mean, na.rm = TRUE)
    } else {
      mean_fields[[ii]] <- NULL
    }

    if (!all(vapply(draw_iids, is.null, logical(1)))) {
      iid_stack <- do.call(c, draw_iids)
      mean_iids[[ii]] <- terra::app(iid_stack, mean, na.rm = TRUE)
    } else {
      mean_iids[[ii]] <- NULL
    }
  }

  names(realizations_list) <- layer_names

  pred_stack <- do.call(c, mean_preds)
  names(pred_stack) <- layer_names
  cov_stack <- do.call(c, mean_covs)
  names(cov_stack) <- layer_names
  lower_stack <- do.call(c, ci_lower)
  names(lower_stack) <- layer_names
  upper_stack <- do.call(c, ci_upper)
  names(upper_stack) <- layer_names

  field_stack <- if (!all(vapply(mean_fields, is.null, logical(1)))) do.call(c, mean_fields) else NULL
  if (!is.null(field_stack)) names(field_stack) <- layer_names
  iid_stack <- if (!all(vapply(mean_iids, is.null, logical(1)))) do.call(c, mean_iids) else NULL
  if (!is.null(iid_stack)) names(iid_stack) <- layer_names

  out <- list(
    mean_prediction = list(
      prediction = pred_stack,
      field = field_stack,
      iid = iid_stack,
      covariates = cov_stack
    ),
    uncertainty_prediction = list(
      realisations = realizations_list,
      predictions_ci = list(
        lower = lower_stack,
        upper = upper_stack
      )
    )
  )
  class(out) <- c("disag_prediction_mmap_mcmc", "disag_prediction_mmap", "list")

  if (verbose) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    message(sprintf("predict.disag_model_mmap_mcmc() runtime: %.2f mins", elapsed))
  }

  out
}


extract_mcmc_draws_matrix_mmap <- function(mcmc_fit) {
  draws <- tryCatch(as.matrix(mcmc_fit), error = function(e) NULL)
  if (is.null(draws) || !is.matrix(draws) || nrow(draws) < 1L || ncol(draws) < 1L) {
    stop("Could not extract posterior draws from `mcmc_fit` as a matrix.", call. = FALSE)
  }
  if (is.null(colnames(draws)) || !length(colnames(draws))) {
    stop("Posterior draws are missing column names; cannot align parameters.", call. = FALSE)
  }

  keep <- !grepl("__$", colnames(draws))
  draws <- draws[, keep, drop = FALSE]
  if (!ncol(draws)) {
    stop("Posterior draw matrix has no model parameter columns after filtering diagnostics.", call. = FALSE)
  }

  draws
}


draw_to_param_list_mcmc_mmap <- function(draw, model_output) {
  ms <- model_output$model_setup
  dat <- model_output$data
  p <- tryCatch(as.integer(ms$coef_meta$p), error = function(e) NA_integer_)
  if (is.na(p)) {
    p <- if (!is.null(dat$covariate_rasters_list) && !is.null(dat$covariate_rasters_list[[1]])) {
      terra::nlyr(dat$covariate_rasters_list[[1]])
    } else {
      0L
    }
  }

  Tn <- length(dat$time_points)
  n_poly <- nrow(dat$polygon_data)
  n_s <- if (isTRUE(ms$field)) nrow(dat$mesh$loc) else 0L
  tv <- isTRUE(ms$time_varying_betas)
  fam <- ms$family
  iidf <- isTRUE(ms$iid)
  fld <- isTRUE(ms$field)

  par_vec <- numeric(0)

  if (tv) {
    par_vec <- c(par_vec, extract_draw_vector_mmap(draw, "intercept_t", Tn))
    par_vec <- c(par_vec, extract_draw_vector_mmap(draw, "slope_t", p * Tn))
  } else {
    par_vec <- c(par_vec, extract_draw_scalar_mmap(draw, "intercept"))
    par_vec <- c(par_vec, extract_draw_vector_mmap(draw, "slope", p))
  }

  if (identical(fam, "gaussian")) {
    par_vec <- c(par_vec, extract_draw_scalar_mmap(draw, "log_tau_gaussian"))
  }

  if (identical(fam, "negbinomial")) {
    par_vec <- c(par_vec, extract_draw_scalar_mmap(draw, "iideffect_log_tau"))
  } else if (iidf) {
    par_vec <- c(par_vec, extract_draw_scalar_mmap(draw, "iideffect_log_tau"))
    par_vec <- c(par_vec, extract_draw_vector_mmap(draw, "iideffect", n_poly))
  }

  if (fld) {
    par_vec <- c(par_vec, extract_draw_scalar_mmap(draw, "log_sigma"))
    par_vec <- c(par_vec, extract_draw_scalar_mmap(draw, "log_rho"))
    par_vec <- c(par_vec, extract_draw_vector_mmap(draw, "nodemean", n_s))
  }

  slice_params_tmb(par_vec, model_output)
}


extract_draw_scalar_mmap <- function(draw, key) {
  nm <- names(draw)
  if (key %in% nm) {
    return(as.numeric(draw[[key]]))
  }

  idx <- grep(paste0("^", key, "\\[[0-9]+\\]$"), nm, value = TRUE)
  if (length(idx) == 1L) {
    index_num <- as.integer(sub(paste0("^", key, "\\[([0-9]+)\\]$"), "\\1", idx))
    if (!is.na(index_num) && index_num == 1L) {
      return(as.numeric(draw[[idx]]))
    }
  }

  stop("Posterior draws are missing required scalar parameter `", key, "`.", call. = FALSE)
}


extract_draw_vector_mmap <- function(draw, key, expected_len) {
  expected_len <- as.integer(expected_len)
  if (expected_len == 0L) {
    return(numeric(0))
  }

  nm <- names(draw)
  idx_names <- grep(paste0("^", key, "\\[[0-9]+\\]$"), nm, value = TRUE)
  if (length(idx_names)) {
    idx_values <- as.integer(sub(paste0("^", key, "\\[([0-9]+)\\]$"), "\\1", idx_names))
    ord <- order(idx_values)
    idx_values <- idx_values[ord]
    idx_names <- idx_names[ord]

    if (anyDuplicated(idx_values)) {
      stop("Posterior draws contain duplicated indices for `", key, "`.", call. = FALSE)
    }
    if (!identical(idx_values, seq_len(expected_len))) {
      stop(
        "Posterior draws for `", key, "` have indices ",
        paste(idx_values, collapse = ", "),
        " but expected 1:", expected_len, ".",
        call. = FALSE
      )
    }
    return(as.numeric(draw[idx_names]))
  }

  if (expected_len == 1L && key %in% nm) {
    return(as.numeric(draw[[key]]))
  }

  stop(
    "Posterior draws are missing required vector parameter `",
    key,
    "` of length ",
    expected_len,
    ".",
    call. = FALSE
  )
}
