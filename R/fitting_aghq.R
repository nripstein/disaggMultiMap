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
  if (!optimizer %in% c("BFGS", "sparse_trust", "trust")) {
    stop("`optimizer` must be one of 'BFGS', 'sparse_trust', or 'trust'.")
  }
  if (verbose) {
    message("Using optimizer: ", optimizer)
  }

  control <- aghq::default_control_tmb(method = optimizer)

  #-- 4. Run AGHQ --

  message("Fitting ", family," disaggregation model via AGHQ (k = ", aghq_k, ").")
  aghq_model <- aghq::marginal_laplace_tmb(
    obj,
    k             = aghq_k,
    startingvalue = obj$par,
    control       = control
  )

  # RENAME
  coef_meta  <- compute_coef_meta(data)
  aghq_model <- rename_aghq_model_names(aghq_model, coef_meta, time_varying_betas)

  #-- 5. (Optional) Extract random-effect mode via sdreport for later prediction
  sd_out <- NULL
  try({ sd_out <- TMB::sdreport(obj) }, silent = TRUE)

  #-- 6. Assemble output --
  # Build fixed-parameter order and beta index map for robust prediction
  theta_order <- try(names(aghq_model$optresults$mode), silent = TRUE)
  if (inherits(theta_order, "try-error")) theta_order <- NULL

  beta_index_map <- NULL
  if (!is.null(theta_order)) {
    Tn <- length(data$time_points)
    p  <- coef_meta$p
    if (isTRUE(time_varying_betas)) {
      # time-varying: intercept_t1..T, and cov_names with _t1..T
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
      beta_index_map <- list(intercept_idx = intercept_idx, slope_idx = slope_idx,
                             tv = TRUE, p = p, Tn = Tn)
    } else {
      # shared betas: 'intercept' and cov_names
      intercept_idx <- match("intercept", theta_order)
      if (is.na(intercept_idx)) stop("Internal: 'intercept' not found in theta order.")
      slope_idx <- if (p > 0L) match(coef_meta$cov_names, theta_order) else integer(0)
      if (p > 0L && any(is.na(slope_idx))) stop("Internal: some slope names not found in theta order.")
      beta_index_map <- list(intercept_idx = intercept_idx, slope_idx = slope_idx,
                             tv = FALSE, p = p, Tn = Tn)
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
      coef_meta = coef_meta,
      theta_order = theta_order,
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
