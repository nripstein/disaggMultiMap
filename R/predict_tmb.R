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

  # Extract time points
  times   <- object$data$time_points
  n_times <- length(times)

  # Prepare lists for stacking
  preds     <- vector("list", n_times)
  fields    <- vector("list", n_times)
  iids      <- vector("list", n_times)
  cov_terms <- vector("list", n_times)

  # Split parameters once
  raw_pars <- object$obj$env$last.par.best
  pars     <- split(raw_pars, names(raw_pars))

  # Loop over each time point
  for (i in seq_along(times)) {
    # i-th covariate raster
    cov_i <- object$data$covariate_rasters_list[[i]]

    # Temporarily assign covariate raster to data so check_new_data() passes
    object$data$covariate_rasters <- cov_i

    # Build prediction objects
    objs <- setup_objects_mmap(object,
                               new_data    = cov_i,
                               predict_iid = predict_iid)

    # Generate prediction components
    out <- predict_single_raster_mmap(pars,
                                      objs,
                                      link_function = object$model_setup$link)

    # Store each component
    preds[[i]]     <- out$prediction
    fields[[i]]    <- out$field
    iids[[i]]      <- out$iid
    cov_terms[[i]] <- out$covariates
  }

  # Stack components into multi-layer SpatRasters
  mean_prediction <- list(
    prediction = do.call(c, preds),
    field      = if (object$model_setup$field)   do.call(c, fields)   else NULL,
    iid        = if (object$model_setup$iid)     do.call(c, iids)     else NULL,
    covariates = do.call(c, cov_terms)
  )

  # Name layers by time index
  layer_names <- paste0("time_", times)
  names(mean_prediction$prediction) <- layer_names
  if (!is.null(mean_prediction$field))      names(mean_prediction$field)      <- layer_names
  if (!is.null(mean_prediction$iid))        names(mean_prediction$iid)        <- layer_names
  names(mean_prediction$covariates) <- layer_names

  # Uncertainty
  uncertainty_prediction <- predict_uncertainty_mmap(
    object, new_data = new_data, predict_iid = predict_iid,
    N = N, CI = CI
  )

  # Return combined predictions
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
#' @param ... Unused.
#'
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
  pars <- split(pars, names(pars))

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
#' @param ... Unused.
#'
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

  # Storage lists
  realizations_list <- vector("list", n_times)
  ci_lower   <- vector("list", n_times)
  ci_upper   <- vector("list", n_times)

  # Draw parameters once
  parameters <- model_output$obj$env$last.par.best
  if (model_output$model_setup$iid || model_output$model_setup$field) {
    ch <- Matrix::Cholesky(model_output$sd_out$jointPrecision)
    par_draws <- sparseMVN::rmvn.sparse(N, parameters, ch, prec = TRUE)
  } else {
    covm <- Matrix::Matrix(model_output$sd_out$cov.fixed, sparse = TRUE)
    ch   <- Matrix::Cholesky(covm)
    par_draws <- sparseMVN::rmvn.sparse(N, parameters, ch, prec = FALSE)
  }

  # Loop over time for uncertainty draws
  for (i in seq_along(times)) {
    cov_i <- model_output$data$covariate_rasters_list[[i]]
    model_output$data$covariate_rasters <- cov_i

    objs <- setup_objects_mmap(model_output,
                               new_data    = cov_i,
                               predict_iid = predict_iid)

    # Monte Carlo draws per time
    draw_preds <- vector("list", N)
    for (r in seq_len(N)) {
      pdraw <- split(par_draws[r, ], names(parameters))
      draw_preds[[r]] <- predict_single_raster_mmap(
        pdraw, objs, link_function = model_output$model_setup$link
      )$prediction
    }
    realizations_list[[i]] <- do.call(c, draw_preds)
    names(realizations_list) <- paste0("time_", times)

    # Compute credible intervals
    probs <- c((1 - CI)/2, 1 - (1 - CI)/2)
    ci_stack <- terra::app(realizations_list[[i]], function(x) stats::quantile(x, probs, na.rm = TRUE))
    ci_lower[[i]] <- ci_stack[[1]]
    ci_upper[[i]] <- ci_stack[[2]]
  }

  # Stack CI rasters and name layers
  layer_names <- paste0("t", times)
  lower_ci_stack <- do.call(c, ci_lower)
  upper_ci_stack <- do.call(c, ci_upper)
  names(lower_ci_stack) <- layer_names
  names(upper_ci_stack) <- layer_names

  return(list(
    realisations   = realizations_list,
    predictions_ci = list(lower = lower_ci_stack,
                          upper = upper_ci_stack)
  ))
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
setup_objects_mmap <- function(model_output, new_data = NULL, predict_iid = FALSE) {
  new_data <- disaggregation:::check_new_data(new_data, model_output)

  # Pull out original data from model_output$data (which is the dis_data_mmap object)
  data <- model_output$data

  # If new_data is provided, use it (assuming same covariate names and ordering)
  if(is.null(new_data)){
    if(!is.null(data$covariate_rasters_list)) {
      covariates <- data$covariate_rasters_list[[length(data$covariate_rasters_list)]]  # using last time point by default
    } else {
      covariates <- NULL
    }
  } else {
    covariates <- new_data
  }

  if(is.null(covariates)) {
    # No covariates were provided so build a zero raster from the stored aggregation template
    template_agg <- data$aggregation_rasters_list[[1]]
    zero_vals <- rep(0, terra::ncell(template_agg))
    covariates <- template_agg
    terra::values(covariates) <- zero_vals
    names(covariates) <- "intercept_only"
  }

  # In multi-map setting, we assume prediction is done over a common spatial extent
  data$covariate_rasters <- covariates

  # For field: extract prediction coordinates from the provided covariate raster
  if(model_output$model_setup$field) {
    if(is.null(new_data)) {
      coords <- data$coords_for_prediction
    } else {
      coords <- disaggregation:::getCoords(data)  # getCoords() can be reused as in original
    }
    Amatrix <- fmesher::fm_evaluator(data$mesh, loc = as.matrix(coords))$proj$A
    field_objects <- list(coords = coords, Amatrix = Amatrix)
  } else {
    field_objects <- NULL
  }

  # For iid effect: use the polygon shapefile from data (if available)
  if(predict_iid && model_output$model_setup$iid) {
    iid_objects <- list(shapefile = model_output$data$polygon_shapefile_list[[1]],  # use first time point's shapefile for prediction grid
                        template = data$covariate_rasters)
  } else {
    iid_objects <- NULL
  }

  return(list(covariates = covariates,
              field_objects = field_objects,
              iid_objects = iid_objects))
}

#' Predict on a single raster for multi-map disaggregation (TMB)
#'
#' @description
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

  # Determine how many layers really exist
  n_layers <- 0
  if (!is.null(objects$covariates)) {
    n_layers <- terra::nlyr(objects$covariates)
    # If the lone layer is our “intercept_only” fallback, treat as zero real covariates
    if (n_layers == 1 && names(objects$covariates)[1] == "intercept_only") {
      n_layers <- 0
    }
  }
  covs_by_betas <- list()

  # consider case of no covariates
  if (n_layers == 0) {
    # Build a raster that everywhere equals model_parameters$intercept on the linear scale.
    # objects$covariates is a single‐layer raster of zeros named "intercept_only"
    cov_contribution <- objects$covariates * 0 + model_parameters$intercept
  } else {
    # Standard path: multiply each real covariate layer by its slope
    for (i in seq_len(n_layers)) {
      covs_by_betas[[i]] <- model_parameters$slope[i] * objects$covariates[[i]]
    }
    cov_by_betas <- terra::rast(covs_by_betas)
    if (terra::nlyr(cov_by_betas) > 1) {
      sum_cov_by_betas <- sum(cov_by_betas)
    } else {
      sum_cov_by_betas <- cov_by_betas
    }
    cov_contribution <- sum_cov_by_betas + model_parameters$intercept
  }

  # Start with covariate contribution as linear predictor
  linear_pred <- cov_contribution

  if (!is.null(objects$field_objects)) {
    # Add the spatial field component
    field <- (objects$field_objects$Amatrix %*% model_parameters$nodemean)[,1]
    field_ras <- terra::rast(cbind(objects$field_objects$coords, field),
                             type = 'xyz',
                             crs  = terra::crs(linear_pred))
    linear_pred <- linear_pred + field_ras
  } else {
    field_ras <- NULL
  }

  if (!is.null(objects$iid_objects)) { # this won't happen if using NB likelihood
    # Incorporate iid effect: rasterize the shapefile with iid values
    objects$iid_objects$shapefile$iid <- model_parameters$iideffect
    iid_ras <- terra::rasterize(objects$iid_objects$shapefile,
                                objects$iid_objects$template,
                                field = 'iid')
    iideffect_sd <- 1 / sqrt(exp(model_parameters$iideffect_log_tau))
    na_pixels <- which(is.na(terra::values(iid_ras, mat = FALSE)))
    if (length(na_pixels) != 0) {
      iid_ras[na_pixels] <- stats::rnorm(length(na_pixels), 0, iideffect_sd)
    }
    if (terra::ext(iid_ras) != terra::ext(linear_pred)) {
      raster_new_extent <- linear_pred
      raster_new_extent[1:terra::ncell(raster_new_extent)] <- NA
      iid_ras <- terra::merge(iid_ras, raster_new_extent)
      iid_ras <- terra::crop(iid_ras, raster_new_extent)
      missing_pixels <- which(is.na(terra::values(iid_ras, mat = FALSE)))
      if (length(missing_pixels) != 0) {
        iid_ras[missing_pixels] <- stats::rnorm(length(missing_pixels), 0, iideffect_sd)
      }
    }
    linear_pred <- linear_pred + iid_ras
  } else {
    iid_ras <- NULL
  }

  # Apply link function to convert linear predictor to prediction scale
  if (link_function == 'logit') {
    prediction_ras <- 1 / (1 + exp(-linear_pred))
  } else if (link_function == 'log') {
    prediction_ras <- exp(linear_pred)
  } else if (link_function == 'identity') {
    prediction_ras <- linear_pred
  } else {
    stop("Link function not implemented.")
  }

  predictions <- list(
    prediction = prediction_ras,
    field      = field_ras,
    iid        = iid_ras,
    covariates = cov_contribution
  )
  return(predictions)
}
