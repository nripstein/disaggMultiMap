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
#' @param fixed_effect_betas Logical; if TRUE (default), beta coefficients are
#'   in AGHQ outer parameters. If FALSE, active betas are treated as TMB random effects.
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
                                  fixed_effect_betas = TRUE,
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
    fixed_effect_betas = fixed_effect_betas,
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
  # Build fixed/random parameter order and beta index map for robust prediction
  theta_order <- try(names(aghq_model$optresults$mode), silent = TRUE)
  if (inherits(theta_order, "try-error")) theta_order <- NULL
  random_order <- try(names(obj$env$par)[obj$env$random], silent = TRUE)
  if (inherits(random_order, "try-error")) random_order <- NULL

  beta_index_map <- NULL
  Tn <- length(data$time_points)
  p  <- coef_meta$p
  beta_source <- if (isTRUE(fixed_effect_betas)) "theta" else "random"
  beta_names <- if (isTRUE(time_varying_betas)) {
    intercept_names <- paste0("intercept_t", seq_len(Tn))
    slope_names <- if (p > 0L) {
      as.vector(vapply(
        seq_len(Tn),
        function(t) paste0(coef_meta$cov_names, "_t", t),
        character(p)
      ))
    } else {
      character(0)
    }
    c(intercept_names, slope_names)
  } else {
    c("intercept", coef_meta$cov_names)
  }
  source_order <- if (identical(beta_source, "theta")) theta_order else random_order
  if (!is.null(source_order)) {
    source_order_norm <- normalize_fixed_names(source_order, coef_meta, time_varying_betas)
    idx_all <- match(beta_names, source_order_norm)
    if (!anyNA(idx_all)) {
      if (isTRUE(time_varying_betas)) {
        intercept_idx <- idx_all[seq_len(Tn)]
        slope_idx <- if (p > 0L) {
          matrix(
            idx_all[-seq_len(Tn)],
            nrow = p,
            ncol = Tn,
            byrow = FALSE
          )
        } else {
          matrix(integer(0), nrow = 0, ncol = Tn)
        }
      } else {
        intercept_idx <- idx_all[1L]
        slope_idx <- if (p > 0L) idx_all[-1L] else integer(0)
      }
      beta_index_map <- list(
        source = beta_source,
        intercept_idx = intercept_idx,
        slope_idx = slope_idx,
        tv = isTRUE(time_varying_betas),
        p = p,
        Tn = Tn
      )
    } else {
      miss <- beta_names[is.na(idx_all)]
      stop(
        "Internal: could not map beta names in ", beta_source, " order. Missing: ",
        paste(miss, collapse = ", ")
      )
    }
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
      fixed_effect_betas = fixed_effect_betas,
      coef_meta = coef_meta,
      theta_order = theta_order,
      random_order = random_order,
      beta_index_map = beta_index_map
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
