#' Fit a multi-map disaggregation model via TMB + MCMC (tmbstan)
#'
#' @description
#' Builds the TMB ADFun object for a multi-map disaggregation model, then
#' fits the model via \code{tmbstan::tmbstan()}.
#'
#' @param data A \code{disag_data_mmap} object (from \code{prepare_data_mmap()}).
#' @param priors Optional named list of prior specifications (see internal helper).
#' @param family One of \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, or \code{"negbinomial"}.
#' @param link One of \code{"identity"}, \code{"logit"}, or \code{"log"}.
#' @param time_varying_betas Logical; if TRUE, each time point has its own fixed-effect.
#' @param field Logical: include the spatial random field?
#' @param iid Logical: include polygon-specific IID effects?
#' @param silent Logical: if TRUE, suppress TMB/tmbstan output where supported.
#' @param starting_values Optional named list of starting parameter values.
#' @param verbose Logical: if TRUE, print total runtime.
#' @param ... Additional arguments passed through to \code{tmbstan::tmbstan()}
#'   (e.g., \code{chains}, \code{iter}, \code{warmup}, \code{cores}, \code{control}).
#'
#' @return An object of class \code{disag_model_mmap_mcmc} (a list with
#'   \code{$mcmc_fit}, \code{$obj}, \code{$data}, \code{$model_setup}, and
#'   \code{$engine_args_used}).
#' @export
disag_model_mmap_mcmc <- function(data,
                                  priors = NULL,
                                  family = "poisson",
                                  link = "log",
                                  time_varying_betas = FALSE,
                                  field = TRUE,
                                  iid = TRUE,
                                  silent = TRUE,
                                  starting_values = NULL,
                                  verbose = FALSE,
                                  ...) {
  start_time <- Sys.time()

  if (!inherits(data, "disag_data_mmap")) {
    stop("`data` must be a 'disag_data_mmap' object; run prepare_data_mmap() first.", call. = FALSE)
  }
  if (!is.null(priors) && !is.list(priors)) {
    stop("`priors` must be NULL or a named list of prior values.", call. = FALSE)
  }
  if (!requireNamespace("tmbstan", quietly = TRUE)) {
    stop(
      paste0(
        "MCMC engine requires suggested package `tmbstan`.\n",
        "Install it with e.g. `install.packages(\"tmbstan\")` ",
        "or a source install appropriate for your R setup."
      ),
      call. = FALSE
    )
  }

  obj <- make_model_object_mmap(
    data = data,
    priors = priors,
    family = family,
    link = link,
    time_varying_betas = time_varying_betas,
    field = field,
    iid = iid,
    silent = silent,
    starting_values = starting_values,
    verbose = verbose
  )

  dots <- list(...)
  dot_names <- names(dots)
  if (is.null(dot_names)) {
    dot_names <- rep("", length(dots))
  }

  if ("obj" %in% dot_names) {
    warning("Ignoring `obj` in `...`; the TMB object is constructed internally.", call. = FALSE)
    dots <- dots[dot_names != "obj"]
    dot_names <- names(dots)
    if (is.null(dot_names)) {
      dot_names <- rep("", length(dots))
    }
  }

  if ("silent" %in% dot_names) {
    warning(
      "Ignoring `silent` in `...`; use the top-level `silent` argument in `disag_model_mmap_mcmc()`.",
      call. = FALSE
    )
    dots <- dots[dot_names != "silent"]
  }

  tmbstan_args <- c(list(obj = obj, silent = silent), dots)

  message("Fitting ", family, " disaggregation model via MCMC (tmbstan).")
  mcmc_fit <- do.call(tmbstan::tmbstan, tmbstan_args)

  out <- list(
    mcmc_fit = mcmc_fit,
    obj = obj,
    data = data,
    model_setup = list(
      family = family,
      link = link,
      field = field,
      iid = iid,
      time_varying_betas = time_varying_betas,
      coef_meta = compute_coef_meta(data)
    ),
    engine_args_used = tmbstan_args[names(tmbstan_args) != "obj"]
  )
  class(out) <- c("disag_model_mmap_mcmc", "disag_model_mmap", "list")

  if (verbose) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    out$runtime_mins <- elapsed
    message(sprintf("disag_model_mmap_mcmc() runtime: %.2f minutes", elapsed))
  }

  out
}
