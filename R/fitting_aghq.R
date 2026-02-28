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
#' @param time_varying_betas Logical; if TRUE, each time point has its own fixed-effect
#' @param aghq_k Integer >= 1: number of quadrature nodes for AGHQ ('1' = Laplace).
#' @param field Logical: include the spatial random field?
#' @param iid Logical: include polygon-specific IID effects?
#' @param silent Logical: if TRUE, suppress TMB's console output.
#' @param starting_values Optional named list of starting parameter values.
#' @param optimizer Optional optimizer name passed to AGHQ control.
#' @param verbose Logical: if TRUE, print total runtime.
#' @return An object of class 'disag_model_mmap_aghq' (a list with '$aghq_model', '$data', and '$model_setup').
#' @export
disag_model_mmap_aghq <- function(data,
                                  priors = NULL,
                                  family = "poisson",
                                  link = "log",
                                  time_varying_betas = FALSE,
                                  aghq_k = 1,
                                  field = TRUE,
                                  iid = TRUE,
                                  silent = TRUE,
                                  starting_values = NULL,
                                  optimizer = NULL,
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
    time_varying_betas = time_varying_betas,
    beta_random_effects = TRUE,
    field           = field,
    iid             = iid,
    silent          = silent,
    starting_values = starting_values,
    verbose         = verbose
  )

  #-- 3. Select optimizer --
  if (is.null(optimizer)) {
    optimizer <- "BFGS"
  }
  if (verbose) {
    message("Using optimizer: ", optimizer)
  }

  #-- 4. Run AGHQ --
  message("Fitting ", family, " disaggregation model via AGHQ (k = ", aghq_k, ").")

  if (optimizer == "nlminb") {
    # nlminb is not supported by CRAN aghq's optimize_theta(), so we run it
    # manually and pass pre-computed results via the optresults parameter.
    # This mirrors the exact computation that the nripstein/aghq fork performs
    # internally: nlminb for optimization, numDeriv::jacobian with Richardson
    # extrapolation for the Hessian.
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr)

    if (opt$convergence != 0) {
      warning("nlminb optimizer did not converge (code ", opt$convergence, "): ", opt$message)
    }

    # Hessian of obj$fn (neg-log-posterior) at mode — positive definite.
    # Uses Richardson extrapolation on obj$gr, matching the fork's approach.
    hess <- numDeriv::jacobian(obj$gr, opt$par, method = "Richardson")

    # Structure mirrors what aghq::optimize_theta() returns internally:
    # - ff$fn/gr/he: log-posterior (negated TMB objective) for quadrature evaluation
    # - hessian: Hessian of neg-log-posterior (positive definite) for quadrature grid
    # - mode: parameter values at the posterior mode
    optresults <- list(
      ff = list(
        fn = function(x) -obj$fn(x),
        gr = function(x) -obj$gr(x),
        he = function(x) -numDeriv::jacobian(obj$gr, x, method = "Richardson")
      ),
      mode        = opt$par,
      hessian     = hess,
      convergence = opt$convergence
    )

    aghq_model <- aghq::marginal_laplace_tmb(
      obj,
      k             = aghq_k,
      startingvalue = opt$par,
      optresults    = optresults
    )
  } else {
    # Standard path: delegate optimization entirely to aghq
    control <- aghq::default_control_tmb(method = optimizer)
    aghq_model <- aghq::marginal_laplace_tmb(
      obj,
      k             = aghq_k,
      startingvalue = obj$par,
      control       = control
    )
  }

  # RENAME
  coef_meta  <- compute_coef_meta(data)
  aghq_model <- rename_aghq_model_names(aghq_model, coef_meta, time_varying_betas)

  #-- 5. Extract random-effect mode via sdreport for later prediction
  sd_out <- NULL
  try({ sd_out <- TMB::sdreport(obj) }, silent = TRUE)

  #-- 6. Assemble output --
  # Build fixed-parameter order and beta index map for robust prediction
  theta_order <- try(names(aghq_model$optresults$mode), silent = TRUE)
  if (inherits(theta_order, "try-error")) theta_order <- NULL

  random_effect_layout <- NULL
  beta_index_map <- NULL
  random_info <- tryCatch(
    build_random_effect_layout_mmap(
      obj = obj,
      aghq_model = aghq_model,
      data = data,
      coef_meta = coef_meta,
      time_varying_betas = time_varying_betas
    ),
    error = function(e) NULL
  )
  if (!is.null(random_info)) {
    random_effect_layout <- random_info$random_effect_layout
    beta_index_map <- random_info$beta_index_map
  }
  # Backward-compatible fallback for older model shapes where betas are in theta.
  if (is.null(beta_index_map) && !is.null(theta_order)) {
    beta_index_map <- build_beta_index_map_from_theta(
      theta_order = theta_order,
      coef_meta = coef_meta,
      time_varying_betas = time_varying_betas,
      Tn = length(data$time_points)
    )
  }

  out <- list(
    aghq_model = aghq_model,
    obj = obj,
    data = data,
    sd_out = sd_out,
    model_setup = list(
      family = family,
      link   = link,
      field  = field,
      iid    = iid,
      time_varying_betas = time_varying_betas,
      coef_meta = coef_meta,
      theta_order = theta_order,
      beta_index_map = beta_index_map,
      random_effect_layout = random_effect_layout
    )
  )
  class(out) <- c("disag_model_mmap_aghq", "disag_model_mmap", "list")

  #-- 7. Runtime message --
  if (verbose) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    message(sprintf("disag_model_mmap_aghq() runtime: %.2f minutes", as.numeric(elapsed)))
  }

  return(out)
}

build_beta_index_map_from_theta <- function(theta_order, coef_meta, time_varying_betas, Tn) {
  p <- coef_meta$p
  if (isTRUE(time_varying_betas)) {
    intercept_idx <- vapply(seq_len(Tn), function(t)
      match(paste0("intercept_t", t), theta_order), integer(1))
    if (p > 0L) {
      slope_idx <- sapply(seq_len(Tn), function(t)
        match(paste0(coef_meta$cov_names, "_t", t), theta_order), simplify = "matrix")
      if (is.null(dim(slope_idx))) slope_idx <- matrix(slope_idx, nrow = p)
    } else {
      slope_idx <- matrix(integer(0), nrow = 0, ncol = Tn)
    }
    if (any(is.na(intercept_idx)) || (p > 0L && any(is.na(slope_idx)))) {
      stop("Internal: could not map time-varying betas to theta order in AGHQ fit.")
    }
    return(list(
      intercept_idx = intercept_idx,
      slope_idx = slope_idx,
      tv = TRUE,
      p = p,
      Tn = Tn,
      space = "theta"
    ))
  }

  intercept_idx <- match("intercept", theta_order)
  if (is.na(intercept_idx)) stop("Internal: 'intercept' not found in theta order.")
  slope_idx <- if (p > 0L) match(coef_meta$cov_names, theta_order) else integer(0)
  if (p > 0L && any(is.na(slope_idx))) {
    missing <- coef_meta$cov_names[is.na(slope_idx)]
    slope_like <- theta_order[grepl("^slope(\\[[0-9]+\\]|\\.[0-9]+|[0-9]+)?$", theta_order)]
    stop(
      "Internal: some slope names not found in theta order.\n",
      "Missing expected slopes: ", paste(missing, collapse = ", "), "\n",
      "Expected slope order: ", paste(coef_meta$cov_names, collapse = ", "), "\n",
      "Observed slope-like theta names: ",
      if (length(slope_like)) paste(slope_like, collapse = ", ") else "<none>"
    )
  }
  list(
    intercept_idx = intercept_idx,
    slope_idx = slope_idx,
    tv = FALSE,
    p = p,
    Tn = Tn,
    space = "theta"
  )
}

build_random_effect_layout_mmap <- function(obj,
                                            aghq_model,
                                            data,
                                            coef_meta,
                                            time_varying_betas) {
  theta_mode <- tryCatch(aghq_model$optresults$mode, error = function(e) NULL)
  random_mode <- tryCatch(aghq_model$modesandhessians$mode[[1]], error = function(e) NULL)
  if (is.null(theta_mode) || is.null(random_mode)) return(NULL)

  random_names <- names(random_mode)
  if (is.null(random_names)) return(NULL)
  random_names <- as.character(random_names)
  get_component_idx <- function(name) as.integer(which(random_names == name))

  random_effect_layout <- list(
    random_length = length(random_mode),
    intercept_idx = get_component_idx("intercept"),
    slope_idx = get_component_idx("slope"),
    intercept_t_idx = get_component_idx("intercept_t"),
    slope_t_idx = get_component_idx("slope_t"),
    nodemean_idx = get_component_idx("nodemean"),
    iideffect_idx = get_component_idx("iideffect")
  )

  Tn <- length(data$time_points)
  p <- coef_meta$p
  if (isTRUE(time_varying_betas)) {
    intercept_idx <- random_effect_layout$intercept_t_idx
    if (length(intercept_idx) != Tn) return(NULL)
    if (p > 0L) {
      slope_flat <- random_effect_layout$slope_t_idx
      if (length(slope_flat) != p * Tn) return(NULL)
      slope_idx <- matrix(slope_flat, nrow = p, ncol = Tn)
    } else {
      slope_idx <- matrix(integer(0), nrow = 0, ncol = Tn)
    }
    beta_index_map <- list(
      intercept_idx = intercept_idx,
      slope_idx = slope_idx,
      tv = TRUE,
      p = p,
      Tn = Tn,
      space = "random_effects"
    )
  } else {
    intercept_idx <- random_effect_layout$intercept_idx
    if (length(intercept_idx) != 1L) return(NULL)
    slope_idx <- if (p > 0L) random_effect_layout$slope_idx else integer(0)
    if (p > 0L && length(slope_idx) != p) return(NULL)
    beta_index_map <- list(
      intercept_idx = intercept_idx,
      slope_idx = slope_idx,
      tv = FALSE,
      p = p,
      Tn = Tn,
      space = "random_effects"
    )
  }

  list(
    random_effect_layout = random_effect_layout,
    beta_index_map = beta_index_map
  )
}
