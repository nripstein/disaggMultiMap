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
#' @param k Integer >= 1: number of quadrature nodes for AGHQ ('1' = Laplace).
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
                                  k = 1,
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

  message("Fitting ", family," disaggregation model via AGHQ (k = ", k, ").")
  aghq_model <- aghq::marginal_laplace_tmb(
    obj,
    k             = k,
    startingvalue = obj$par,
    control       = control
  )

  #-- 5. Assemble output --
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
  class(out) <- c("disag_model_mmap_aghq", "disag_model_mmap", "list")

  #-- 6. Runtime message --
  if (verbose) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    message(sprintf("disag_model_mmap_aghq() runtime: %.2f minutes", as.numeric(elapsed)))
  }

  return(out)
}
