#' Compute coefficient metadata (number of timepoints, covariates, and their names)
#' @param data A disag_data_mmap object.
#' @keywords internal
compute_coef_meta <- function(data) {
  reserved  <- c("ID","cell","poly_local_id","time")
  cov_names <- setdiff(names(data$covariate_data), reserved)
  list(
    n_times   = length(data$time_points),
    p         = length(cov_names),
    cov_names = cov_names
  )
}

# Derive coef meta from data (no flag duplication)
compute_coef_meta <- function(data) {
  reserved  <- c("ID", "cell", "poly_local_id", "time")
  cov_names <- setdiff(names(data$covariate_data), reserved)
  list(
    n_times   = length(data$time_points),
    p         = length(cov_names),
    cov_names = cov_names
  )
}

# Normalize parameter names in outputs
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

# R/utils_names.R
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
  return(aghq_model)
}
