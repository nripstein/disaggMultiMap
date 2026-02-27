#' Predict for Multi-Map Disaggregation Model fit with TMB
#'
#' @param object A fitted disag_model_mmap_tmb object.
#' @param new_data Optionally, a new SpatRaster (or list of them) for prediction.
#' @param predict_iid Logical. If TRUE, include the polygon iid effect in predictions.
#' @param N Number of Monte Carlo draws for uncertainty estimation.
#' @param CI Credible interval level (default 0.95).
#' @param ... Further arguments.
#'
#' @export
#' @method predict disag_model_mmap_tmb
predict.disag_model_mmap_tmb <- function(object, new_data = NULL,
                                         predict_iid = FALSE, N = 100, CI = 0.95, ...) {
  times   <- object$data$time_points
  n_times <- length(times)

  preds     <- vector("list", n_times)
  fields    <- vector("list", n_times)
  iids      <- vector("list", n_times)
  cov_terms <- vector("list", n_times)

  raw_pars <- object$obj$env$last.par.best
  pars     <- parvec_to_param_list(object, raw_pars)

  for (i in seq_along(times)) {
    cov_i <- object$data$covariate_rasters_list[[i]]
    object$data$covariate_rasters <- cov_i

    objs <- setup_objects_mmap(object,
                               new_data    = NULL,
                               predict_iid = predict_iid,
                               time_index  = i,
                               use_training = TRUE)

    out <- predict_single_raster_mmap(pars,
                                      objs,
                                      link_function = object$model_setup$link)

    preds[[i]]     <- out$prediction
    fields[[i]]    <- out$field
    iids[[i]]      <- out$iid
    cov_terms[[i]] <- out$covariates
  }

  mean_prediction <- list(
    prediction = do.call(c, preds),
    field      = if (object$model_setup$field)   do.call(c, fields)   else NULL,
    iid        = if (object$model_setup$iid)     do.call(c, iids)     else NULL,
    covariates = do.call(c, cov_terms)
  )
  layer_names <- paste0("time_", times)

  names(mean_prediction$prediction) <- layer_names



  if (!is.null(mean_prediction$field)) names(mean_prediction$field) <- layer_names
  if (!is.null(mean_prediction$iid))   names(mean_prediction$iid)   <- layer_names
  names(mean_prediction$covariates)    <- layer_names

  uncertainty_prediction <- predict_uncertainty_mmap(
    object, new_data = new_data, predict_iid = predict_iid,
    N = N, CI = CI
  )

  prediction <- list(
    mean_prediction        = mean_prediction,
    uncertainty_prediction = uncertainty_prediction
  )
  class(prediction) <- c("disag_prediction_mmap", "list")
  prediction
}


#' Predict mean for multi-map disaggregation (TMB)
#'
#' @description
#' Given a fitted TMB model object and optional new covariate data,
#' compute the mean‐only prediction (no uncertainty) for one raster.
#'
#' @param model_output A 'disag_model_mmap_tmb' model fit.
#' @param new_data Optional SpatRaster (or list) of new covariates.
#' @param predict_iid Logical; if TRUE, include the IID polygon effect.
#' @return A list with components:
#'   - prediction: SpatRaster of the mean prediction on the response scale.
#'   - field: SpatRaster of the spatial field component (or NULL).
#'   - iid: SpatRaster of the IID effect (or NULL).
#'   - covariates: SpatRaster of the linear predictor from covariates only.
#' @keywords internal
predict_model_mmap <- function(model_output, new_data = NULL, predict_iid = FALSE, newdata = NULL) {
  if (!is.null(newdata) && missing(new_data)) {
    new_data <- newdata
    message("newdata is deprecated and will be removed in a future version - please use new_data instead")
  }

  objects_for_prediction <- setup_objects_mmap(model_output, new_data = new_data, predict_iid = predict_iid)

  # Get the last parameter estimates
  pars <- model_output$obj$env$last.par.best
  pars <- parvec_to_param_list(model_output, pars)

  prediction <- predict_single_raster_mmap(pars,
                                           objects_for_prediction,
                                           link_function = model_output$model_setup$link)

  return(prediction)
}

#' Estimate uncertainty via Monte Carlo for multi-map disaggregation (TMB)
#'
#' @description
#' Draw Monte Carlo samples of model parameters, propagate through the
#' prediction function, and compute credible intervals at each cell.
#'
#' @param model_output A 'disag_model_mmap_tmb' model fit.
#' @param new_data Optional SpatRaster (or list) of new covariates.
#' @param predict_iid Logical; if TRUE, include the IID polygon effect.
#' @param N Integer; number of Monte Carlo draws (default 100).
#' @param CI Numeric in (0,1); credible‐interval level (default 0.95).
#' @return A list with components:
#'   - realisations: list of SpatRasters of each draw.
#'   - predictions_ci: list with 'lower' and 'upper' SpatRaster stacks.
#' @keywords internal
predict_uncertainty_mmap <- function(model_output, new_data = NULL,
                                     predict_iid = FALSE, N = 100, CI = 0.95, newdata = NULL) {
  if (!is.null(newdata) && missing(new_data)) {
    new_data <- newdata
    message("'newdata' is deprecated; please use 'new_data' instead")
  }

  times   <- model_output$data$time_points
  n_times <- length(times)

  realizations_list <- vector("list", n_times)
  ci_lower <- vector("list", n_times)
  ci_upper <- vector("list", n_times)

  parameters <- model_output$obj$env$last.par.best
  if (model_output$model_setup$iid || model_output$model_setup$field) {
    ch <- Matrix::Cholesky(model_output$sd_out$jointPrecision)
    par_draws <- sparseMVN::rmvn.sparse(N, parameters, ch, prec = TRUE)
  } else {
    covm <- Matrix::Matrix(model_output$sd_out$cov.fixed, sparse = TRUE)
    ch   <- Matrix::Cholesky(covm)
    par_draws <- sparseMVN::rmvn.sparse(N, parameters, ch, prec = FALSE)
  }

  for (i in seq_along(times)) {
    cov_i <- model_output$data$covariate_rasters_list[[i]]
    model_output$data$covariate_rasters <- cov_i

    objs <- setup_objects_mmap(model_output,
                               new_data    = NULL,
                               predict_iid = predict_iid,
                               time_index  = i,
                               use_training = TRUE)

    draw_preds <- vector("list", N)
    for (r in seq_len(N)) {
      pdraw <- slice_params_tmb(par_draws[r, ], model_output)
      draw_preds[[r]] <- predict_single_raster_mmap(
        pdraw, objs, link_function = model_output$model_setup$link
      )$prediction
    }
    realizations_list[[i]] <- do.call(c, draw_preds)
    names(realizations_list) <- paste0("time_", times)

    probs <- c((1 - CI)/2, 1 - (1 - CI)/2)
    ci_stack <- terra::app(realizations_list[[i]], function(x) stats::quantile(x, probs, na.rm = TRUE))
    ci_lower[[i]] <- ci_stack[[1]]
    ci_upper[[i]] <- ci_stack[[2]]
  }

  layer_names <- paste0("time_", times)  # harmonize with mean stack
  lower_ci_stack <- do.call(c, ci_lower); names(lower_ci_stack) <- layer_names
  upper_ci_stack <- do.call(c, ci_upper); names(upper_ci_stack) <- layer_names

  list(
    realisations   = realizations_list,
    predictions_ci = list(lower = lower_ci_stack, upper = upper_ci_stack)
  )
}

#' Prepare prediction objects for multi-map disaggregation (TMB)
#'
#' @description
#' Constructs the covariate rasters, field projector, and IID shapefile
#' objects needed by the single‐raster prediction routines.
#'
#' @param model_output A 'disag_model_mmap_tmb' model fit.
#' @param new_data Optional SpatRaster (or list) of new covariates.
#' @param predict_iid Logical; if TRUE, include the IID polygon effect.
#'
#' @return A list with elements:
#'   - covariates: SpatRaster of covariate layers.
#'   - field_objects: list with 'coords' matrix and 'Amatrix' projector (or NULL).
#'   - iid_objects: list with 'shapefile' and 'template' (or NULL).
#' @keywords internal
setup_objects_mmap <- function(model_output, new_data = NULL, predict_iid = FALSE, time_index = NULL, use_training = FALSE) {
  new_data <- disagg_check_new_data(new_data, model_output)

  data <- model_output$data

  # Decide time index early
  ti <- if (is.null(time_index)) length(data$time_points) else as.integer(time_index)
  if (ti < 1L) stop("time_index must be a positive integer.", call. = FALSE)

  # Choose covariates
  if (isTRUE(use_training)) {
    if (!is.null(data$covariate_rasters_list)) {
      covariates <- data$covariate_rasters_list[[ti]]
    } else {
      covariates <- NULL
    }
  } else {
    if (is.null(new_data)) {
      if (!is.null(data$covariate_rasters_list)) {
        covariates <- data$covariate_rasters_list[[length(data$covariate_rasters_list)]]
      } else {
        covariates <- NULL
      }
    } else {
      covariates <- new_data
    }
  }

  if (is.null(covariates)) {
    template_agg <- data$aggregation_rasters_list[[ti]]
    zero_vals <- rep(0, terra::ncell(template_agg))
    covariates <- template_agg
    terra::values(covariates) <- zero_vals
    names(covariates) <- "intercept_only"
  }

  categorical_schema <- tryCatch(model_output$data$categorical_covariate_schema, error = function(e) NULL)
  if (inherits(covariates, "SpatRaster") &&
      !(terra::nlyr(covariates) == 1L && identical(names(covariates)[1], "intercept_only"))) {
    covariates <- encode_categorical_raster_stack(
      covariates = covariates,
      categorical_schema = categorical_schema,
      time_index = ti,
      context = if (isTRUE(use_training)) {
        "predict.disag_model_mmap_tmb() using training covariates"
      } else {
        "predict.disag_model_mmap_tmb() using new_data"
      }
    )
  }

  data$covariate_rasters <- covariates

  if (model_output$model_setup$field) {
    coords <- if (is.null(new_data)) data$coords_for_prediction else disagg_get_coords(data)
    Amatrix <- fmesher::fm_evaluator(data$mesh, loc = as.matrix(coords))$proj$A
    # Sanity check projector against mesh nodes
    if (ncol(Amatrix) != nrow(data$mesh$loc)) {
      stop(sprintf(
        paste0(
          "Field projector columns (", ncol(Amatrix), ") do not match mesh nodes (", nrow(data$mesh$loc), ").\n",
          "This suggests a mismatch in the mesh used for prediction vs training."
        )
      ), call. = FALSE)
    }
    field_objects <- list(coords = coords, Amatrix = Amatrix)
  } else {
    field_objects <- NULL
  }

  if (predict_iid && model_output$model_setup$iid && !identical(model_output$model_setup$family, "negbinomial")) {
    iid_objects <- list(shapefile = model_output$data$polygon_shapefile_list[[1]],
                        template  = data$covariate_rasters)
  } else {
    iid_objects <- NULL
  }

  # --- NEW: time plumbing + guards ------------------------------------------
  tv <- isTRUE(model_output$model_setup$time_varying_betas)
  if (is.null(model_output$model_setup$time_varying_betas)) {
    # Infer tv flag from presence of intercept_t in the fitted parameter vector
    pl_try <- try(slice_params_tmb(model_output$obj$env$last.par.best, model_output), silent = TRUE)
    if (!inherits(pl_try, "try-error")) {
      tv <- !is.null(pl_try$intercept_t) && length(pl_try$intercept_t) > 0
    }
  }

  # If real covariates are present (not the intercept_only placeholder),
  # enforce count and names against training metadata.
  has_real_covs <- !(terra::nlyr(covariates) == 1 && names(covariates)[1] == "intercept_only")
  if (has_real_covs) {
    meta <- tryCatch(model_output$model_setup$coef_meta, error = function(e) NULL)
    p_expected <- if (is.null(meta)) NA_integer_ else meta$p
    train_names <- if (is.null(meta)) NULL else meta$cov_names

    if (!is.na(p_expected) && terra::nlyr(covariates) != p_expected) {
      stop(sprintf("Mismatch in covariate layer count for prediction (got %d, expected %d).",
                   terra::nlyr(covariates), p_expected), call. = FALSE)
    }

    if (!is.null(train_names) && length(train_names)) {
      cur_names <- names(covariates)
      if (setequal(cur_names, train_names) && !identical(cur_names, train_names)) {
        covariates <- covariates[[train_names]]
      } else if (!setequal(cur_names, train_names)) {
        stop(sprintf(
          paste0(
            "Covariate layer names do not match training.\n",
            "Expected: %s\nGot     : %s\n",
            "If names are the same set but different order, ensure they can be reordered."
          ),
          paste(train_names, collapse = ", "),
          paste(cur_names, collapse = ", ")
        ), call. = FALSE)
      }
    }
  }

  # Dimension guard: projector rows must match prediction cells
  if (model_output$model_setup$field) {
    n_cells <- terra::ncell(covariates)
    if (nrow(field_objects$Amatrix) != n_cells) {
      stop(sprintf(
        paste0(
          "Field projector dimension mismatch.\n",
          "Rows(A) = %d; ncell(covariates) = %d.\n",
          "Check that prediction rasters share extent/resolution with the training grid."
        ),
        nrow(field_objects$Amatrix), n_cells
      ), call. = FALSE)
    }
  }
  # ---------------------------------------------------------------------------

  list(
    covariates          = covariates,
    field_objects       = field_objects,
    iid_objects         = iid_objects,
    time_index          = ti,  # <-- NEW
    time_varying_betas  = tv   # <-- NEW (read in predict_single_raster_mmap)
  )
}

#' Predict on a single raster for multi-map disaggregation (TMB)
#'
#' Apply the linear predictor, add spatial field and IID components,
#' then transform via the link to return a SpatRaster prediction.
#'
#' @param model_parameters Named list of parameter vectors (split by name).
#' @param objects List from 'setup_objects_mmap' containing data and projectors.
#' @param link_function Character; one of 'identity', 'log', or 'logit'.
#'
#' @return A list with components:
#'   - prediction: SpatRaster on the response scale.
#'   - field: SpatRaster of field contribution (or NULL).
#'   - iid: SpatRaster of IID contribution (or NULL).
#'   - covariates: SpatRaster of the covariate linear predictor.
#' @keywords internal
predict_single_raster_mmap <- function(model_parameters, objects, link_function) {
  # How many real covariate layers?
  n_layers <- 0L
  if (!is.null(objects$covariates)) {
    n_layers <- terra::nlyr(objects$covariates)
    if (n_layers == 1L && identical(names(objects$covariates)[1], "intercept_only")) {
      n_layers <- 0L
    }
  }

  # Time-varying flag and index (plumbed via setup_objects_mmap)
  tv <- isTRUE(objects$time_varying_betas)
  ti <- if (!is.null(objects$time_index)) as.integer(objects$time_index) else 1L
  if (ti < 1L) stop("time_index must be a positive integer.", call. = FALSE)

  # ---------- Covariate contribution on linear scale ----------
  covs_by_betas <- list()
  if (!tv) {
    # Shared betas
    if (n_layers == 0L) {
      cov_contribution <- objects$covariates * 0 + model_parameters$intercept
    } else {
      for (j in seq_len(n_layers)) {
        covs_by_betas[[j]] <- model_parameters$slope[j] * objects$covariates[[j]]
      }
      cov_by_betas <- terra::rast(covs_by_betas)
      sum_cov_by_betas <- if (terra::nlyr(cov_by_betas) > 1L) terra::app(cov_by_betas, sum) else cov_by_betas
      cov_contribution <- sum_cov_by_betas + model_parameters$intercept
    }
  } else {
    # Time-varying betas
    if (is.null(model_parameters$intercept_t) || ti > length(model_parameters$intercept_t)) {
      stop(sprintf("intercept_t does not contain entry for time_index %d.", ti), call. = FALSE)
    }
    if (n_layers == 0L) {
      cov_contribution <- objects$covariates * 0 + model_parameters$intercept_t[ti]
    } else {
      p <- n_layers
      if (is.null(model_parameters$slope_t))
        stop("slope_t is missing in model_parameters for time-varying prediction.", call. = FALSE)
      if ((length(model_parameters$slope_t) %% p) != 0L)
        stop(sprintf("Length of slope_t (%d) is not a multiple of p (%d).",
                     length(model_parameters$slope_t), p), call. = FALSE)
      start <- (ti - 1L) * p + 1L
      end   <- start + p - 1L
      if (end > length(model_parameters$slope_t))
        stop(sprintf("Requested slopes for time_index %d exceed length of slope_t.", ti), call. = FALSE)
      slopes_this_time <- model_parameters$slope_t[start:end]

      for (j in seq_len(p)) {
        covs_by_betas[[j]] <- slopes_this_time[j] * objects$covariates[[j]]
      }
      cov_by_betas <- terra::rast(covs_by_betas)
      sum_cov_by_betas <- if (terra::nlyr(cov_by_betas) > 1L) terra::app(cov_by_betas, sum) else cov_by_betas
      cov_contribution <- sum_cov_by_betas + model_parameters$intercept_t[ti]
    }
  }

  # Start with covariate contribution as linear predictor
  linear_pred <- cov_contribution
  lin_template <- terra::rast(linear_pred)  # canonical grid for everything that follows

  # ---------- Spatial field (align to template, keep 1 layer) ----------
  if (!is.null(objects$field_objects)) {
    # Project to prediction points
    A <- objects$field_objects$Amatrix
    nm <- model_parameters$nodemean
    if (ncol(A) != length(nm)) {
      stop(sprintf(
        paste0(
          "Field vector length mismatch.\n",
          "ncol(A) = %d, length(nodemean) = %d.\n",
          "This indicates a parameter slicing/order issue or a mesh mismatch."
        ), ncol(A), length(nm)
      ), call. = FALSE)
    }
    field_vec <- as.numeric((A %*% nm)[, 1])

    # Safer mapping: write field values to the exact grid using cell indices
    idx <- terra::cellFromXY(lin_template, objects$field_objects$coords)
    if (anyNA(idx)) stop("Could not map field coordinates onto prediction grid.", call. = FALSE)

    field_ras <- terra::rast(lin_template)   # one-layer, correct geom
    field_ras[] <- NA_real_
    field_ras[idx] <- field_vec

    linear_pred <- linear_pred + field_ras
  } else {
    field_ras <- NULL
  }

  # ---------- IID polygon effect (align to template, keep 1 layer) ----------
  if (!is.null(objects$iid_objects)) { # not used for NB likelihood
    # Stamp polygon IID values onto the template grid directly
    objects$iid_objects$shapefile$iid <- model_parameters$iideffect
    iid_ras <- terra::rasterize(
      objects$iid_objects$shapefile,
      lin_template,
      field = "iid"
    )

    # Fill missing pixels with draws from the IID prior SD
    iideffect_sd <- 1 / sqrt(exp(model_parameters$iideffect_log_tau))
    na_pixels <- which(is.na(terra::values(iid_ras, mat = FALSE)))
    if (length(na_pixels))
      iid_ras[na_pixels] <- stats::rnorm(length(na_pixels), 0, iideffect_sd)

    linear_pred <- linear_pred + iid_ras
  } else {
    iid_ras <- NULL
  }

  # ---------- Link to response scale ----------
  if (identical(link_function, "logit")) {
    prediction_ras <- 1 / (1 + exp(-linear_pred))
  } else if (identical(link_function, "log")) {
    prediction_ras <- exp(linear_pred)
  } else if (identical(link_function, "identity")) {
    prediction_ras <- linear_pred
  } else {
    stop("Link function not implemented.", call. = FALSE)
  }

  # Single-layer invariant (catch regressions early)
  if (terra::nlyr(prediction_ras) != 1L) {
    stop(sprintf("Internal: expected a single-layer prediction, got %d layers.",
                 terra::nlyr(prediction_ras)), call. = FALSE)
  }

  list(
    prediction = prediction_ras,
    field      = field_ras,
    iid        = iid_ras,
    covariates = cov_contribution
  )
}
