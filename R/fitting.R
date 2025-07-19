#' Fit a multi-map disaggregation model via TMB + AGHQ
#'
#' @description
#' Builds the TMB ADFun object for a multi-map disaggregation model, then
#' fits the model via AGHQ with desired number of quadrature points.
#'
#' @param data A 'disag_data_mmap' object (from 'prepare_data_mmap()').
#' @param priors Optional named list of prior specifications (see internal helper).
#' @param family One of "gaussian", "binomial", "poisson", or "negbinomial".
#' @param link One of "identity", "logit", or "log".
#' @param k Integer ≥ 1: number of quadrature nodes for AGHQ ('1' = Laplace).
#' @param field Logical: include the spatial random field?
#' @param iid Logical: include polygon-specific IID effects?
#' @param silent Logical: if TRUE, suppress TMB's console output.
#' @param starting_values Optional named list of starting parameter values.
#' @param verbose Logical: if TRUE, print total runtime.
#' @return An object of class 'disag_model_mmap_aghq' (a list with '$aghq_model', '$data', and '$model_setup').
#' @export
disag_model_mmap_aghq <- function(data,
                                  priors = NULL,
                                  family = "poisson",
                                  link = "log",
                                  k = 3,
                                  field = TRUE,
                                  iid = TRUE,
                                  silent = TRUE,
                                  starting_values = NULL,
                                  verbose = FALSE) {
  start_time <- Sys.time()

  #-- 1. Input validation --
  if (!inherits(data, "disag_data_mmap")) {
    stop("`data` must be a 'disag_data_mmap' object; run prepare_data_mmap() first.")
  }
  if (!is.null(priors) && !is.list(priors)) {
    stop("`priors` must be NULL or a named list of prior values.")
  }

  #-- 2. Build TMB ADFun object --
  obj <- make_model_object_mmap(
    data            = data,
    priors          = priors,
    family          = family,
    link            = link,
    field           = field,
    iid             = iid,
    silent          = silent,
    starting_values = starting_values,
    verbose         = verbose
  )

  #-- 3. Run AGHQ --
  message("Fitting ", family," disaggregation model via AGHQ (k = ", k, ").")
  aghq_model <- aghq::marginal_laplace_tmb(
    obj,
    k             = k,
    startingvalue = obj$par,
    control       = aghq::default_control_tmb()
  )

  #-- 4. Assemble output --
  out <- list(
    aghq_model = aghq_model,
    data = data,
    model_setup = list(
      family = family,
      link   = link,
      field  = field,
      iid    = iid
    )
  )
  class(out) <- c("disag_model_mmap_aghq", "list")

  #-- 5. Runtime message --
  if (verbose) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    message(sprintf("disag_model_mmap_aghq() runtime: %.2f minutes", as.numeric(elapsed)))
  }

  return(out)
}

#' Build the TMB ADFun object for multi-map disaggregation
#'
#' @description
#' Internal helper. Converts data, priors, and model settings into the list
#' of inputs required by 'TMB::MakeADFun()'.
#'
#' @param data A 'disag_data_mmap' object.
#' @param priors NULL or named list overriding default hyperpriors.
#' @param family One of "gaussian", "binomial", "poisson", "negbinomial".
#' @param link One of "identity", "logit", "log".
#' @param field Logical: include spatial field?
#' @param iid Logical: include IID polygon effects?
#' @param silent Logical: pass to 'MakeADFun()' to suppress output.
#' @param starting_values NULL or named list of starting values.
#' @return A 'TMB::ADFun' object ready for 'marginal_laplace_tmb()'.
#' @keywords internal
make_model_object_mmap <- function(data,
                                   priors = NULL,
                                   family = "gaussian",
                                   link = "identity",
                                   field = TRUE,
                                   iid = TRUE,
                                   silent = TRUE,
                                   starting_values = NULL,
                                   verbose = FALSE) {
  #-- 1. Validate data prerequisites --
  if (!inherits(data, "disag_data_mmap")) {
    stop("Internal error: `data` must be 'disag_data_mmap'.")
  }
  if (family == "binomial" && anyNA(data$polygon_data$N)) {
    stop("Binomial models require non-NA sample sizes (`data$polygon_data$N`).")
  }
  if (is.null(data$mesh)) {
    stop("`data` must include a mesh; run prepare_data_mmap(make_mesh = TRUE).")
  }

  #-- 2. Map family name to integer --
  if (family == "gaussian") {
    family_id <- 0
  } else if (family == "binomial") {
    family_id <- 1
  } else if (family == "poisson") {
    family_id <- 2
  } else if (family == "negbinomial") {
    family_id <- 3
  } else {
    stop("`family` must be one of 'gaussian', 'binomial', 'poisson', 'negbinomial'.")
  }

  #-- 3. Map link name to integer --
  if (link == "logit") {
    link_id <- 0
  } else if (link == "log") {
    link_id <- 1
  } else if (link == "identity") {
    link_id <- 2
  } else {
    stop("`link` must be one of 'identity', 'logit', 'log'.")
  }

  #-- 4. Warning for redundant settings --
  if (family == "gaussian" && iid) {
    warning(
      "Gaussian likelihood with IID effect is redundant; consider `iid = FALSE`.",
      call. = FALSE
    )
  }

  #-- 5. Build SPDE operators & projection matrix --
  nu <- 1
  spde <- rSPDE::matern.operators(
    mesh = data$mesh,
    alpha = nu + 1,
    compute_higher_order = TRUE
  )$fem_mesh_matrices
  spde[[4]] <- NULL
  names(spde) <- c("M0", "M1", "M2")

  Apix <- fmesher::fm_evaluator(data$mesh, loc = data$coords_for_fit)$proj$A
  n_s <- nrow(spde$M0)

  #-- 6. Prepare covariate matrix --
  # Select covariate columns after preprocessing, including any one-hot dummies:
  reserved <- c("ID", "cell", "poly_local_id", "time")
  cov_cols <- setdiff(names(data$covariate_data), reserved)
  if (length(cov_cols) == 0) {
    has_covariates <- FALSE
    cov_matrix <- matrix(0.0, nrow = nrow(data$covariate_data), ncol = 0)
  } else {
    has_covariates <- TRUE
    cov_matrix <- as.matrix(data$covariate_data[, cov_cols, drop = FALSE])
    # ensure numeric
    cov_matrix <- apply(cov_matrix, 2, as.numeric)
    colnames(cov_matrix) <- cov_cols
  }


  if (verbose) {
    total_polys <- length(unique(data$covariate_data$poly_local_id))
    for (col in cov_cols) {
      v <- tapply(
        data$covariate_data[[col]],
        data$covariate_data$poly_local_id,
        var,
        na.rm = TRUE
      )
      zero_var_pols <- names(v)[!is.finite(v) | v == 0]
      n_zero <- length(zero_var_pols)
      if (n_zero > 0) {
        msg_list <- paste(zero_var_pols, collapse = ", ")
        if (n_zero > 10) {
          msg_list <- paste0(paste(head(zero_var_pols, 10), collapse = ", "), ", ...")
        }
        message(sprintf(
          "Covariate '%s' has zero variance in %d/%d polygons (%s). Poisson models may be unidentifiable.",
          col, n_zero, total_polys, msg_list
        ))
      }
    }
  }

  #-- 7. Set up default hyperpriors --
  bbox <- sf::st_bbox(data$polygon_shapefile_list[[1]])
  diag_len <- sqrt((bbox$xmax - bbox$xmin)^2 + (bbox$ymax - bbox$ymin)^2)
  prior_rho <- diag_len / 3
  prior_sigma <- sd(data$polygon_data$response / mean(data$polygon_data$response))

  default_priors <- list(
    priormean_intercept     = 0,
    priorsd_intercept       = 10.0,
    priormean_slope         = 0.0,
    priorsd_slope           = 0.5,
    prior_rho_min           = prior_rho,
    prior_rho_prob          = 0.1,
    prior_sigma_max         = prior_sigma,
    prior_sigma_prob        = 0.1,
    prior_iideffect_sd_max  = 0.1,
    prior_iideffect_sd_prob = 0.01
  )

  # Merge user priors
  if (!is.null(priors)) {
    invalid <- setdiff(names(priors), names(default_priors))
    if (length(invalid)) {
      stop("Invalid prior names: ", paste(invalid, collapse = ", "))
    }
    final_priors <- modifyList(default_priors, priors)
  } else {
    final_priors <- default_priors
  }

  #-- 8. Flatten start/end index across time --
  pixel_offset <- 0L
  start_end_index_combined <- NULL
  for (idx in seq_along(data$start_end_index)) {
    cur <- data$start_end_index[[idx]]
    adj <- cur + pixel_offset
    start_end_index_combined <- rbind(start_end_index_combined, adj)
    pixel_offset <- pixel_offset + (max(cur) + 1L)
  }

  #-- 9. Build parameter list & map for TMB --
  default_parameters <- list(
    intercept         = 0,
    slope             = rep(0, ncol(cov_matrix)),
    log_tau_gaussian  = 8,
    iideffect         = rep(0, nrow(data$polygon_data)),
    iideffect_log_tau = 1,
    log_sigma         = 0,
    log_rho           = 4,
    nodemean          = rep(0, n_s)
  )
  parameters <- if (is.null(starting_values)) {
    default_parameters
  } else {
    modifyList(default_parameters, starting_values)
  }

  # rename covariates
  names(default_parameters$slope) <- cov_cols

  # Build data list for TMB
  input_data <- c(
    list(
      x = cov_matrix,
      aggregation_values = data$aggregation_pixels,
      Apixel = Apix,
      spde = spde,
      start_end_index = start_end_index_combined,
      polygon_response_data = data$polygon_data$response,
      response_sample_size = data$polygon_data$N,
      time = data$polygon_data$time,
      family = family_id,
      link = link_id,
      nu = nu,
      field = as.integer(field),
      iid = as.integer(iid)
    ),
    final_priors
  )

  #-- 10. Build map for fixed parameters --
  tmb_map <- list()
  if (!field) {
    tmb_map <- c(tmb_map, list(
      log_sigma = factor(NA),
      log_rho   = factor(NA),
      nodemean  = factor(rep(NA, n_s))
    ))
  }
  if (family_id == 3) { # NB
    tmb_map <- c(tmb_map, list(
      iideffect = factor(rep(NA, nrow(data$polygon_data))) # no iid effect for each polygon, but still include iideffect_log_tau
    ))
  }
  if (!iid) {
    tmb_map <- c(tmb_map, list(
      iideffect_log_tau = factor(NA),
      iideffect         = factor(rep(NA, nrow(data$polygon_data)))
    ))
  }
  if (family_id != 0) { # non‐Gaussian
    tmb_map <- c(tmb_map, list(log_tau_gaussian = factor(NA)))
  }
  if (!has_covariates) {
    tmb_map <- c(tmb_map, list(slope = factor(rep(NA, 0))))
  }

  #-- 11. Identify random effects --
  random_effects <- character(0)
  if (field) random_effects <- c(random_effects, "nodemean")
  if (iid) random_effects <- c(random_effects, "iideffect")

  #-- 12. Make objective function in TMB--
  obj <- TMB::MakeADFun(
    data       = input_data,
    parameters = parameters,
    map        = tmb_map,
    random     = random_effects,
    silent     = silent,
    DLL        = "disaggMultiMap"
  )

  return(obj)
}
