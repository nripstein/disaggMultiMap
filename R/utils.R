#' Get Default Prior Values for Disaggregation Model
#'
#' @description
#' Calculates the default Penalized Complexity (PC) prior parameters and Gaussian
#' priors that will be used by \code{disag_model_mmap()} if the user does not
#' provide overrides.
#'
#' @details
#' The default priors are dynamic and depend on the input data:
#' \itemize{
#'   \item \strong{Range (Rho):} The lower bound \code{prior_rho_min} is set to
#'   1/3 of the diagonal length of the study area's bounding box.
#'   \item \strong{Spatial SD (Sigma):} The upper bound \code{prior_sigma_max}
#'   is set to the coefficient of variation of the polygon response counts.
#' }
#'
#' @param data A \code{disag_data_mmap} object (output from \code{prepare_data_mmap}).
#'
#' @return A named list of prior specifications.
#' @export
#'
#' @examples
#' \dontrun{
#' # Check defaults before fitting
#' defaults <- get_priors(my_data)
#' print(defaults)
#'
#' # Use defaults as a base to modify specific values
#' my_priors <- defaults
#' my_priors$prior_rho_prob <- 0.05 # Stricter probability
#' fit <- disag_model_mmap(my_data, priors = my_priors)
#' }
#' @keywords external
get_priors <- function(data) {

  if (!inherits(data, "disag_data_mmap")) {
    stop("Input 'data' must be a 'disag_data_mmap' object created by prepare_data_mmap().")
  }

  # 1. Calculate dynamic Rho (Range) lower bound
  bbox <- sf::st_bbox(data$polygon_shapefile_list[[1]])
  diag_len <- sqrt((bbox$xmax - bbox$xmin)^2 + (bbox$ymax - bbox$ymin)^2)
  calc_rho_min <- diag_len / 3

  # 2. Calculate dynamic Sigma (SD) upper bound
  resp <- data$polygon_data$response
  # Handle edge case where mean is 0 to avoid Inf
  if (mean(resp, na.rm = TRUE) == 0) {
    calc_sigma_max <- 1.0
  } else {
    calc_sigma_max <- sd(resp / mean(resp, na.rm = TRUE), na.rm = TRUE)
  }

  # 3. Assemble the full list of defaults
  # These match the hardcoded values in make_model_object_mmap
  priors <- list(
    # Intercepts & Slopes
    priormean_intercept     = 0.0,
    priorsd_intercept       = 10.0,
    priormean_slope         = 0.0,
    priorsd_slope           = 0.5,

    # Spatial Field (PC Priors)
    prior_rho_min           = calc_rho_min,
    prior_rho_prob          = 0.1,  # 10% prob range is smaller than rho_min

    prior_sigma_max         = calc_sigma_max,
    prior_sigma_prob        = 0.1,  # 10% prob sigma is larger than sigma_max

    # IID / Negative Binomial Dispersion (PC Priors)
    prior_iideffect_sd_max  = 0.1,  # Default assumes small extra variance
    prior_iideffect_sd_prob = 0.01  # 1% prob SD is larger than 0.1
  )

  return(priors)
}

#' Compute coefficient metadata (number of timepoints, covariates, and their names)
#' @param data A disag_data_mmap object.
#' @keywords internal
compute_coef_meta <- function(data) {
  reserved  <- c("ID","cell","poly_local_id","time")
  cov_names <- setdiff(names(data$covariate_data), reserved)
  if (is.null(cov_names)) cov_names <- character(0)
  list(n_times   = length(data$time_points),
       p         = length(cov_names),
       cov_names = cov_names)
}

#' Normalize fixed-effect parameter names
#'
#' @description
#' Normalize fixed-effect parameter names to consistent labels.
#' @param nm Character vector of parameter names.
#' @param coef_meta List with coefficient metadata from `compute_coef_meta()`.
#' @param time_varying_betas Logical indicating whether time-varying betas are used.
#' @return Character vector of normalized parameter names.
#' @keywords internal
normalize_fixed_names <- function(nm, coef_meta, time_varying_betas) {
  if (is.null(nm)) return(nm)
  p   <- coef_meta$p
  Tn  <- coef_meta$n_times
  cvs <- coef_meta$cov_names

  if (!isTRUE(time_varying_betas)) {
    # Shared case: slope -> cov_names (in order), intercept unchanged
    idx <- which(nm == "slope" | grepl("^slope(\\[|\\.|$)", nm))
    if (length(idx) == p && p > 0L) nm[idx] <- cvs
    return(nm)
  }

  # Time-varying case: make intercept_t1..Tn and cov_tX names (time-major)
  idx_i <- which(grepl("^intercept_t(\\[|\\.|$|[0-9]+$)?", nm))
  if (length(idx_i) == Tn && Tn > 0L) nm[idx_i] <- paste0("intercept_t", seq_len(Tn))

  idx_s <- which(grepl("^slope_t(\\[|\\.|$|[0-9]+$)?", nm))
  if (length(idx_s) == p * Tn && p > 0L) {
    nm[idx_s] <- as.vector(vapply(
      seq_len(Tn),
      function(t) paste0(cvs, "_t", t),
      character(p)
    ))
  }
  nm
}


#' Rename AGHQ model parameter names
#'
#' @description
#' Rename parameter entries in an AGHQ model object (marginals, optimizer
#' outputs, modes, Hessians) using canonicalized names based on coefficient
#' metadata and time-varying structure.
#' @param aghq_model Fitted AGHQ model object.
#' @param coef_meta List with coefficient metadata from `compute_coef_meta()`.
#' @param time_varying_betas Logical indicating whether time-varying betas are used.
#' @return AGHQ model object with normalized parameter names in key components.
#' @keywords internal
rename_aghq_model_names <- function(aghq_model, coef_meta, time_varying_betas) {
  # 1) 1D marginals (this is what aghq::summary() displays)
  if (!is.null(aghq_model$marginals) && !is.null(names(aghq_model$marginals))) {
    names(aghq_model$marginals) <-
      normalize_fixed_names(names(aghq_model$marginals), coef_meta, time_varying_betas)
  }

  # 2) Optimizer results (mode/parameters)
  if (!is.null(aghq_model$optresults)) {
    for (fld in c("par","mode","theta")) {
      if (!is.null(aghq_model$optresults[[fld]]))
        names(aghq_model$optresults[[fld]]) <-
          normalize_fixed_names(names(aghq_model$optresults[[fld]]), coef_meta, time_varying_betas)
    }
  }

  # 3) Modes and Hessians
  if (!is.null(aghq_model$modesandhessians)) {
    if (!is.null(aghq_model$modesandhessians$mode))
      names(aghq_model$modesandhessians$mode) <-
        normalize_fixed_names(names(aghq_model$modesandhessians$mode), coef_meta, time_varying_betas)

    if (!is.null(aghq_model$modesandhessians$H)) {
      rn <- rownames(aghq_model$modesandhessians$H)
      cn <- colnames(aghq_model$modesandhessians$H)
      if (!is.null(rn)) rownames(aghq_model$modesandhessians$H) <-
        normalize_fixed_names(rn, coef_meta, time_varying_betas)
      if (!is.null(cn)) colnames(aghq_model$modesandhessians$H) <-
        normalize_fixed_names(cn, coef_meta, time_varying_betas)
    }
  }

  # 4) Normalized posterior thetanames (ordering for samples)
  if (!is.null(aghq_model$normalized_posterior)) {
    tn <- try(aghq_model$normalized_posterior$thetanames, silent = TRUE)
    if (!inherits(tn, "try-error") && !is.null(tn)) {
      aghq_model$normalized_posterior$thetanames <-
        normalize_fixed_names(tn, coef_meta, time_varying_betas)
    }
  }
  return(aghq_model)
}


#' Canonicalize posterior draw parameter names
#'
#' @description
#' Canonicalize parameter draw names for consistency across shared vs.
#' time-varying intercepts/slopes while leaving random effects and
#' hyperparameters unchanged.
#' @param old_names Character vector of draw names.
#' @param coef_meta List with coefficient metadata from `compute_coef_meta()`.
#' @param time_varying_betas Logical indicating whether time-varying betas are used.
#' @return Character vector of canonicalized draw names.
#' @keywords internal
canonicalize_draw_names <- function(old_names, coef_meta, time_varying_betas) {
  # Defensive defaults
  if (is.null(old_names)) return(old_names)
  new <- old_names

  # Pull meta safely
  p   <- tryCatch(coef_meta$p,       error = function(e) 0L)
  Tn  <- tryCatch(coef_meta$n_times, error = function(e) NA_integer_)
  cvs <- tryCatch(coef_meta$cov_names, error = function(e) character(0))

  if (is.null(p) || is.na(p))  p  <- 0L
  if (is.null(Tn) || is.na(Tn)) {
    # fallback: infer n_times from intercept_t count if present
    Tn <- sum(grepl("^intercept_t(\\[|\\.|$|[0-9]+$)?", old_names))
  }
  if (is.null(cvs)) cvs <- character(0)

  # --- collapse random effects & IID to base names so split() groups correctly ---
  new[grepl("^nodemean(\\[.*\\])?$", old_names)]  <- "nodemean"
  new[grepl("^iideffect(\\[.*\\])?$", old_names)] <- "iideffect"

  # Shared case
  if (!isTRUE(time_varying_betas)) {
    if (p > 0L) {
      idx_s <- which(old_names == "slope" | grepl("^slope(\\[|\\.|$)", old_names))
      # Only rename if lengths match exactly
      if (length(idx_s) == p) new[idx_s] <- cvs
    }
    return(new)
  }

  # Time-varying case
  # Intercepts (match variants: intercept_t, intercept_t1, intercept_t.1, intercept_t[1])
  if (!is.na(Tn) && Tn > 0L) {
    idx_i <- which(grepl("^intercept_t(\\[[0-9]+\\]|\\.[0-9]+|[0-9]+)?$", old_names))
    if (length(idx_i) == Tn) new[idx_i] <- paste0("intercept_t", seq_len(Tn))
  }

  # Slopes (match variants: slope_t, slope_t1, slope_t.1, slope_t[1])
  if (p > 0L && !is.na(Tn) && Tn > 0L) {
    idx_s <- which(grepl("^slope_t(\\[[0-9]+\\]|\\.[0-9]+|[0-9]+)?$", old_names))
    if (length(idx_s) == p * Tn) {
      tv_names <- as.vector(vapply(seq_len(Tn), function(t) paste0(cvs, "_t", t), character(p)))
      new[idx_s] <- tv_names
    }
  }

  new
}


# Slice a raw parameter vector into a minimal list the predictor needs, by index.
# No names are used; only sizes/flags.
slice_params_tmb <- function(par_vec, model_output) {
  # Pull meta
  ms   <- model_output$model_setup
  dat  <- model_output$data
  p    <- if (!is.null(dat$covariate_rasters_list) && !is.null(dat$covariate_rasters_list[[1]]))
    terra::nlyr(dat$covariate_rasters_list[[1]]) else 0L
  Tn   <- length(dat$time_points)
  fam  <- ms$family
  fld  <- isTRUE(ms$field)
  iidf <- isTRUE(ms$iid)
  tv   <- isTRUE(ms$time_varying_betas)

  n_poly <- nrow(dat$polygon_data)

  idx <- 1L

  # Fixed effects
  if (tv) {
    intercept_t <- par_vec[idx:(idx+Tn-1L)]; idx <- idx + Tn
    slope_t_len <- p * Tn
    slope_t <- if (slope_t_len > 0L) par_vec[idx:(idx+slope_t_len-1L)] else numeric(0)
    idx <- idx + slope_t_len
    intercept <- NULL
    slope     <- numeric(0)
  } else {
    intercept <- par_vec[idx]; idx <- idx + 1L
    slope     <- if (p > 0L) par_vec[idx:(idx+p-1L)] else numeric(0)
    idx <- idx + p
    intercept_t <- NULL
    slope_t     <- numeric(0)
  }

  # Gaussian-only noise (mapped out otherwise)
  log_tau_gaussian <- if (identical(fam, "gaussian")) { val <- par_vec[idx]; idx <- idx + 1L; val } else NULL

  # IID / NB overdispersion
  # NB (family == "negbinomial"): keep iideffect_log_tau, no iideffect vector
  # Non-NB with iid: iideffect_log_tau + iideffect vector
  iideffect_log_tau <- NULL
  iideffect         <- numeric(0)
  if (identical(fam, "negbinomial")) {
    iideffect_log_tau <- par_vec[idx]; idx <- idx + 1L
  } else if (iidf) {
    iideffect_log_tau <- par_vec[idx]; idx <- idx + 1L
    iideffect <- par_vec[idx:(idx+n_poly-1L)]; idx <- idx + n_poly
  }

  # Field hyper + nodemean
  log_sigma <- NULL; log_rho <- NULL; nodemean <- numeric(0)
  if (fld) {
    log_sigma <- par_vec[idx]; idx <- idx + 1L
    log_rho   <- par_vec[idx]; idx <- idx + 1L
    # Remaining tail is nodemean
    if (idx <= length(par_vec)) {
      nodemean <- par_vec[idx:length(par_vec)]
      idx <- length(par_vec) + 1L
    }
  }

  list(
    # shared
    intercept = intercept,
    slope     = slope,
    # time-varying
    intercept_t = intercept_t,
    slope_t     = slope_t,
    # random/hyper
    iideffect_log_tau = iideffect_log_tau,
    iideffect         = iideffect,
    log_tau_gaussian  = log_tau_gaussian,
    log_sigma         = log_sigma,
    log_rho           = log_rho,
    nodemean          = nodemean,
    # tiny meta the predictor needs
    p   = p,
    Tn  = Tn,
    tv  = tv
  )
}


#' Convert a TMB parameter vector to a named parameter list using TMB's parList
#' Falls back to slice_params_tmb if parList fails (e.g., shape mismatch)
#' @keywords internal
parvec_to_param_list <- function(model_output, par_vec) {
  # Use robust, order-based slicing aligned with MakeADFun parameter list
  slice_params_tmb(par_vec, model_output)
}
