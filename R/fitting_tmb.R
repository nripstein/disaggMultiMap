#' Fit a multi-map disaggregation model via TMB
#'
#' @description
#' Builds the TMB ADFun object for a multi-map disaggregation model, then
#' fits the model by maximizing the TMB objective and approximates uncertainty
#' via the optimized Hessian.
#'
#' @param data A 'disag_data_mmap' object (from 'prepare_data_mmap()').
#' @param priors Optional named list of prior specifications (see internal helper).
#' @param family One of 'gaussian', 'binomial', 'poisson', or 'negbinomial'.
#' @param link One of 'identity', 'logit', or 'log'.
#' @param time_varying_betas Logical; if TRUE, each time point has its own fixed-effect
#' @param iterations Integer >= 1: maximum number of optimizer iterations.
#' @param field Logical: include the spatial random field?
#' @param iid Logical: include polygon-specific IID effects?
#' @param hess_control_parscale Optional numeric vector for scaling the Hessian steps.
#' @param hess_control_ndeps Numeric; relative step size for Hessian finite-difference (default 1e-4).
#' @param silent Logical: if TRUE, suppress TMB's console output.
#' @param starting_values Optional named list of starting parameter values.
#' @param verbose Logical: if TRUE, print total runtime.
#'
#' @return An object of class 'disag_model_mmap_tmb' (a list with '$obj', '$opt',
#'   '$sd_out', '$data', and '$model_setup').
#' @export
disag_model_mmap_tmb <- function(data,
                                 priors = NULL,
                                 family = 'poisson',
                                 link = 'log',
                                 time_varying_betas = FALSE,
                                 iterations = 1000,
                                 field = TRUE,
                                 iid = TRUE,
                                 hess_control_parscale = NULL,
                                 hess_control_ndeps = 1e-4,
                                 silent = TRUE,
                                 starting_values = NULL,
                                 verbose = FALSE) {
  start_time <- Sys.time()

  if(!inherits(data, 'disag_data_mmap')) {
    stop("data must be an object of class 'disag_data_mmap'. Use prepare_data_mmap() first.")
  }
  if(!is.null(priors)) stopifnot(inherits(priors, 'list'))
  stopifnot(is.numeric(iterations))

  obj <- make_model_object_mmap(data = data,
                                priors = priors,
                                family = family,
                                link = link,
                                time_varying_betas = time_varying_betas,
                                field = field,
                                iid = iid,
                                silent = silent,
                                starting_values = starting_values)

  message("Fitting ", family," disaggregation model via TMB.")
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                       control = list(iter.max = iterations, trace = 0))

  if(opt$convergence != 0) warning('The model did not converge. Try changing starting_values')

  hess_control <- disagg_setup_hess_control(opt, hess_control_parscale, hess_control_ndeps)
  hess <- stats::optimHess(opt$par, fn = obj$fn, gr = obj$gr, control = hess_control)

  sd_out <- TMB::sdreport(obj, getJointPrecision = TRUE, hessian.fixed = hess)

  # Rename slope parameters to match covariate layer names (assuming all time points share the same names)
  # cov_names <- names(data$covariate_rasters_list[[1]])
  # names(sd_out$par.fixed)[names(sd_out$par.fixed) == "slope"] <- cov_names
  # names(opt$par)[names(opt$par) == "slope"] <- cov_names

  ## UPGRADED RENAMING:
  coef_meta <- compute_coef_meta(data)
  # Normalize names on sd_out$par.fixed and opt$par
  if (!is.null(sd_out$par.fixed)) {
    names(sd_out$par.fixed) <- normalize_fixed_names(names(sd_out$par.fixed), coef_meta, time_varying_betas)
  }
  if (!is.null(opt$par)) {
    names(opt$par) <- normalize_fixed_names(names(opt$par), coef_meta, time_varying_betas)
  }

  model_output <- list(obj = obj,
                       opt = opt,
                       sd_out = sd_out,
                       data = data,
                       model_setup = list(family = family,
                                          link = link,
                                          field = field,
                                          iid = iid,
                                          time_varying_betas = time_varying_betas,
                                          coef_meta = compute_coef_meta(data)))

  class(model_output) <- c("disag_model_mmap_tmb", "disag_model_mmap", "list")

  if (verbose == TRUE) {
    end_time <- Sys.time()
    time_taken <- as.numeric(difftime(end_time, start_time, units = "mins"))
    message(sprintf("disag_model_mmap() runtime: %.2f minutes\n", time_taken))
  }
  return(model_output)
}
