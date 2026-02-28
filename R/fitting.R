#' Fit a multi-map disaggregation model (via AGHQ or TMB)
#'
#' @description
#' Top-level fitting wrapper with engine dispatch and engine-specific argument
#' handling. Engine-specific controls should be supplied via \code{engine.args}.
#'
#' @param data A \code{disag_data_mmap} object.
#' @param priors Optional named list of prior overrides.
#' @param family One of \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, or \code{"negbinomial"}.
#' @param link One of \code{"identity"}, \code{"logit"}, or \code{"log"}.
#' @param engine Character; either \code{"AGHQ"} or \code{"TMB"}.
#' @param time_varying_betas Logical; if TRUE, each time point has its own fixed-effect.
#' @param engine.args Optional named list of engine-specific options.
#'   Supported keys:
#'   \itemize{
#'     \item AGHQ: \code{aghq_k}, \code{optimizer}
#'     \item TMB: \code{iterations}, \code{hess_control_parscale}, \code{hess_control_ndeps}
#'   }
#' @param aghq_k Deprecated at wrapper level; use \code{engine.args = list(aghq_k = ...)}.
#'   Retained for backward compatibility.
#' @param field Logical; include spatial field?
#' @param iid Logical; include IID polygon effects?
#' @param silent Logical; pass through to engine fit function.
#' @param starting_values Optional named list of starting values.
#' @param optimizer Deprecated at wrapper level; use
#'   \code{engine.args = list(optimizer = ...)}. Retained for backward compatibility.
#' @param verbose Logical; print runtime diagnostics.
#' @param ... Additional arguments. Engine-specific arguments passed via \code{...}
#'   are deprecated in this wrapper and should be moved to \code{engine.args}.
#'
#' @return A fitted model object of class \code{disag_model_mmap_tmb} or
#'   \code{disag_model_mmap_aghq} (both also inherit \code{disag_model_mmap}).
#' @export
disag_model_mmap <- function(data,
                             priors = NULL,
                             family = "poisson",
                             link   = "log",
                             engine = c("AGHQ","TMB"),
                             time_varying_betas = FALSE,
                             engine.args = NULL,
                             aghq_k = 2, # Deprecated wrapper argument; prefer engine.args$aghq_k
                             field           = TRUE,
                             iid             = TRUE,
                             silent          = TRUE,
                             starting_values = NULL,
                             optimizer       = NULL, # Deprecated wrapper argument; prefer engine.args$optimizer
                             verbose         = FALSE,
                             ...) {
  engine <- match.arg(engine)
  engine_specs <- get_engine_specs_mmap()
  engine_spec <- engine_specs[[engine]]

  dots <- list(...)
  dot_names <- names(dots)
  if (is.null(dot_names)) {
    dot_names <- rep("", length(dots))
  }

  engine.args <- validate_engine_args_container(engine.args)

  # Extract named engine-specific arguments supplied via ... (deprecated path).
  idx_engine_dots <- which(!is.na(dot_names) & nzchar(dot_names) & dot_names %in% engine_spec$engine_keys)
  dot_engine_args <- if (length(idx_engine_dots)) dots[idx_engine_dots] else list()
  passthrough_non_engine_dots <- if (length(idx_engine_dots)) dots[-idx_engine_dots] else dots

  if (length(dot_engine_args)) {
    warning(
      paste0(
        "Engine-specific arguments in `...` are deprecated in `disag_model_mmap()` (",
        paste(unique(names(dot_engine_args)), collapse = ", "),
        "). Use `engine.args = list(...)`."
      ),
      call. = FALSE
    )
  }

  # Legacy dedicated arguments retained for backward compatibility.
  legacy_named_args <- list()
  used_aghq_k <- !missing(aghq_k)
  used_optimizer <- !missing(optimizer)

  if (engine == "AGHQ") {
    if (used_aghq_k) {
      warning(
        "`aghq_k` in `disag_model_mmap()` is deprecated; use `engine.args = list(aghq_k = ...)`.",
        call. = FALSE
      )
      legacy_named_args$aghq_k <- aghq_k
    }
    if (used_optimizer) {
      warning(
        "`optimizer` in `disag_model_mmap()` is deprecated; use `engine.args = list(optimizer = ...)`.",
        call. = FALSE
      )
      legacy_named_args$optimizer <- optimizer
    }
  } else {
    if (used_aghq_k) {
      warning(
        paste0("`aghq_k` is AGHQ-specific and was ignored because `engine = \"", engine, "\"`."),
        call. = FALSE
      )
    }
    if (used_optimizer) {
      warning(
        paste0("`optimizer` is AGHQ-specific and was ignored because `engine = \"", engine, "\"`."),
        call. = FALSE
      )
    }
  }

  resolved_engine_args <- resolve_engine_args_mmap(
    engine = engine,
    engine_spec = engine_spec,
    engine_args = engine.args,
    legacy_named_args = legacy_named_args,
    dot_engine_args = dot_engine_args
  )
  resolved_engine_args <- validate_engine_specific_values(
    engine = engine,
    resolved_engine_args = resolved_engine_args
  )

  common_args <- list(
    data               = data,
    priors             = priors,
    family             = family,
    link               = link,
    time_varying_betas = time_varying_betas,
    field              = field,
    iid                = iid,
    silent             = silent,
    starting_values    = starting_values,
    verbose            = verbose
  )

  dispatch_args <- c(common_args, resolved_engine_args, passthrough_non_engine_dots)
  do.call(engine_spec$fit_fun, dispatch_args)
}

get_engine_specs_mmap <- function() {
  list(
    AGHQ = list(
      fit_fun = disag_model_mmap_aghq,
      engine_keys = c("aghq_k", "optimizer"),
      defaults = list(aghq_k = 2L)
    ),
    TMB = list(
      fit_fun = disag_model_mmap_tmb,
      engine_keys = c("iterations", "hess_control_parscale", "hess_control_ndeps"),
      defaults = list()
    )
  )
}

validate_engine_args_container <- function(engine.args) {
  if (is.null(engine.args)) {
    return(list())
  }
  if (!is.list(engine.args)) {
    stop("`engine.args` must be NULL or a named list.", call. = FALSE)
  }
  if (!length(engine.args)) {
    return(list())
  }

  nm <- names(engine.args)
  if (is.null(nm)) {
    stop("`engine.args` must be a named list.", call. = FALSE)
  }
  if (any(is.na(nm) | !nzchar(nm))) {
    stop("`engine.args` contains missing or empty names.", call. = FALSE)
  }
  if (anyDuplicated(nm)) {
    dup <- unique(nm[duplicated(nm)])
    stop(
      "`engine.args` contains duplicated names: ",
      paste(dup, collapse = ", "),
      call. = FALSE
    )
  }

  engine.args
}

resolve_engine_args_mmap <- function(engine,
                                     engine_spec,
                                     engine_args,
                                     legacy_named_args,
                                     dot_engine_args) {
  allowed <- engine_spec$engine_keys

  unknown_engine_args <- setdiff(names(engine_args), allowed)
  if (length(unknown_engine_args)) {
    warning(
      paste0(
        "Ignoring unknown `engine.args` key(s) for engine ",
        engine,
        ": ",
        paste(unknown_engine_args, collapse = ", "),
        ". Allowed keys are: ",
        if (length(allowed)) paste(allowed, collapse = ", ") else "<none>",
        "."
      ),
      call. = FALSE
    )
    keep_idx <- names(engine_args) %in% allowed
    engine_args <- engine_args[keep_idx]
  }

  resolved <- engine_spec$defaults
  source_map <- if (length(resolved)) {
    stats::setNames(rep("default", length(resolved)), names(resolved))
  } else {
    stats::setNames(character(0), character(0))
  }

  apply_source <- function(src_list, source_label) {
    if (!length(src_list)) return(invisible(NULL))
    for (ii in seq_along(src_list)) {
      nm <- names(src_list)[ii]
      if (!nzchar(nm) || !(nm %in% allowed)) next

      if (nm %in% names(source_map) && source_map[[nm]] != "default" && source_map[[nm]] != source_label) {
        warning(
          paste0(
            "Argument conflict for `", nm, "`: using value from ", source_label,
            " and overriding value from ", source_map[[nm]], "."
          ),
          call. = FALSE
        )
      }
      resolved[[nm]] <<- src_list[[ii]]
      source_map[[nm]] <<- source_label
    }
  }

  # Lowest to highest precedence.
  apply_source(dot_engine_args, "`...`")
  apply_source(legacy_named_args, "deprecated top-level arguments")
  apply_source(engine_args, "`engine.args`")

  resolved
}

validate_engine_specific_values <- function(engine, resolved_engine_args) {
  is_scalar_integerish <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && abs(x - round(x)) < 1e-8
  }

  if (engine == "AGHQ") {
    if ("aghq_k" %in% names(resolved_engine_args)) {
      k <- resolved_engine_args$aghq_k
      if (!is_scalar_integerish(k) || k < 1) {
        stop("`aghq_k` must be an integer-like scalar >= 1.", call. = FALSE)
      }
      resolved_engine_args$aghq_k <- as.integer(round(k))
    }
    if ("optimizer" %in% names(resolved_engine_args)) {
      opt <- resolved_engine_args$optimizer
      if (!is.character(opt) || length(opt) != 1L || !nzchar(opt)) {
        stop("`optimizer` must be a non-empty character scalar.", call. = FALSE)
      }
    }
  }

  if (engine == "TMB") {
    if ("iterations" %in% names(resolved_engine_args)) {
      it <- resolved_engine_args$iterations
      if (!is_scalar_integerish(it) || it < 1) {
        stop("`iterations` must be an integer-like scalar >= 1.", call. = FALSE)
      }
      resolved_engine_args$iterations <- as.integer(round(it))
    }
    if ("hess_control_parscale" %in% names(resolved_engine_args)) {
      ps <- resolved_engine_args$hess_control_parscale
      if (!(is.null(ps) || is.numeric(ps))) {
        stop("`hess_control_parscale` must be numeric or NULL.", call. = FALSE)
      }
    }
    if ("hess_control_ndeps" %in% names(resolved_engine_args)) {
      nd <- resolved_engine_args$hess_control_ndeps
      if (!is.numeric(nd) || length(nd) != 1L || !is.finite(nd) || nd <= 0) {
        stop("`hess_control_ndeps` must be a numeric scalar > 0.", call. = FALSE)
      }
    }
  }

  resolved_engine_args
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
#' @param time_varying_betas Logical; if TRUE, each time point has its own fixed-effect
#' @param beta_random_effects Logical; internal switch controlling whether beta
#'   coefficients are treated as random effects in TMB's inner Laplace step.
#' @param field Logical: include spatial field?
#' @param iid Logical: include IID polygon effects?
#' @param silent Logical: pass to 'MakeADFun()' to suppress output.
#' @param starting_values NULL or named list of starting values.
#' @param optimizer Optional; For changing the arguments used in AGHQ.
#' @param verbose Logical: if TRUE, print details throughout including runtime.
#' @return A 'TMB::ADFun' object ready for 'marginal_laplace_tmb()'.
#' @keywords external
make_model_object_mmap <- function(data,
                                   priors = NULL,
                                   family = "gaussian",
                                   link = "identity",
                                   time_varying_betas = FALSE,
                                   beta_random_effects = FALSE,
                                   field = TRUE,
                                   iid = TRUE,
                                   silent = TRUE,
                                   starting_values = NULL,
                                   optimizer = NULL,
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

  #-- Early validation for time-varying betas and covariate consistency --
  validate_timevarying_covariates(
    covariate_rasters_list = data$covariate_rasters_list,
    time_varying_betas     = time_varying_betas,
    where = "make_model_object_mmap()"
  )

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
        stats::var,
        na.rm = TRUE
      )
      zero_var_pols <- names(v)[!is.finite(v) | v == 0]
      n_zero <- length(zero_var_pols)
      if (n_zero > 0) {
        msg_list <- paste(zero_var_pols, collapse = ", ")
        if (n_zero > 10) {
          msg_list <- paste0(paste(utils::head(zero_var_pols, 10), collapse = ", "), ", ...")
        }
        message(sprintf(
          "Covariate '%s' has zero variance in %d/%d polygons (%s). Poisson models may be unidentifiable.",
          col, n_zero, total_polys, msg_list
        ))
      }
    }
  }

  # -- 6b. Time metadata for pixel rows (for time-dependent betas) --
  n_times <- length(data$time_points)
  pixel_time_index <- as.integer(data$covariate_data$time)
  ncov <- as.integer(ncol(cov_matrix))

  #-- 7. Set up default hyperpriors --
  bbox <- sf::st_bbox(data$polygon_shapefile_list[[1]])
  diag_len <- sqrt((bbox$xmax - bbox$xmin)^2 + (bbox$ymax - bbox$ymin)^2)
  prior_rho <- diag_len / 3
  prior_sigma <- stats::sd(data$polygon_data$response / mean(data$polygon_data$response))

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
    final_priors <- utils::modifyList(default_priors, priors)
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
    # Shared betas
    intercept         = -5,
    slope             = rep(0, ncov),

    # Time-varying betas
    intercept_t       = rep(-5, n_times),
    slope_t           = rep(0, ncov * n_times),

    # Hyper/latent parameters
    log_tau_gaussian  = 8,
    iideffect         = rep(0, nrow(data$polygon_data)),
    iideffect_log_tau = 1,
    log_sigma         = 0,
    log_rho           = 4,
    nodemean          = rep(0, n_s)
  )

  # Names for slopes
  # Shared
  if (ncov > 0L) {
    names(default_parameters$slope) <- colnames(cov_matrix)  # or cov_cols
  }

  # Time-varying (time-major order: all p covariates for t=1, then t=2, …)
  if (n_times > 0L) {
    names(default_parameters$intercept_t) <- paste0("intercept_t", seq_len(n_times))
  }
  if (ncov > 0L && n_times > 0L) {
    names(default_parameters$slope_t) <- paste0(
      rep(colnames(cov_matrix), times = n_times),
      "_t",
      rep(seq_len(n_times), each = ncov)
    )
  }


  parameters <- if (is.null(starting_values)) {
    default_parameters
  } else {
    utils::modifyList(default_parameters, starting_values)
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
      iid = as.integer(iid),

      # new
      pixel_time_index   = pixel_time_index,
      n_times            = as.integer(n_times),
      ncov               = as.integer(ncov),
      time_varying_betas = as.integer(time_varying_betas)
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
  # -- Always drop the iideffect vector for NB --
  if (family_id == 3) {
    tmb_map <- c(tmb_map, list(
      iideffect = factor(rep(NA, nrow(data$polygon_data)))
    ))
  }

  # -- Only drop iideffect_log_tau when iid=FALSE and it's not NB --
  if (!iid && family_id != 3) {
    tmb_map <- c(tmb_map, list(
      iideffect_log_tau = factor(NA),
      iideffect         = factor(rep(NA, nrow(data$polygon_data)))
    ))
  }
  if (family_id != 0) { # non‐Gaussian
    tmb_map <- c(tmb_map, list(log_tau_gaussian = factor(NA)))
  }
  # if (!has_covariates) {
  #   tmb_map <- c(tmb_map, list(slope = factor(rep(NA, 0))))
  # }

  if (time_varying_betas) {     # Keep: intercept_t, slope_t
    tmb_map <- c(tmb_map, list(
      intercept = factor(NA),
      slope     = factor(rep(NA, ncov))   # length zero if ncov == 0 (which is ok)
    ))
  } else {     # Keep: intercept, slope
    tmb_map <- c(tmb_map, list(
      intercept_t = factor(rep(NA, n_times)),
      slope_t     = factor(rep(NA, ncov * n_times))  # length zero if p == 0 (which is ok)
    ))
  }

  #-- 11. Identify random effects --
  random_effects <- character(0)
  if (isTRUE(beta_random_effects)) {
    # AGHQ path: integrate only covariance hyperparameters in the outer step.
    if (time_varying_betas) {
      random_effects <- c(random_effects, "intercept_t", "slope_t")
    } else {
      random_effects <- c(random_effects, "intercept", "slope")
    }
  }
  if (field) random_effects <- c(random_effects, "nodemean")
  if (iid && family_id != 3) { # include polygon-specific random‐effect vector when iid and not NB
    random_effects <- c(random_effects, "iideffect")
  }

  #-- 12. Make objective function in TMB--
  obj <- TMB::MakeADFun(
    data       = input_data,
    parameters = parameters,
    map        = tmb_map,
    random     = random_effects,
    silent     = silent,
    DLL        = "disaggMultiMap"
    # inner_control = list(trace = 10, REPORT=1) # for debugging
  )

  return(obj)
}


#' Validate covariate layer consistency across time
#'
#' @param covariate_rasters_list list or NULL; each element is a multilayer raster/brick/spatRaster
#' @param time_varying_betas logical; when TRUE, enforce identical layer names and order across time
#' @param where character; caller label for clearer error messages
#' @keywords internal
validate_timevarying_covariates <- function(covariate_rasters_list,
                                            time_varying_betas,
                                            where = "make_model_object_mmap()") {
  if (!isTRUE(time_varying_betas)) return(TRUE)

  # Intercept-only automatically passes test
  if (is.null(covariate_rasters_list)) return(TRUE)

  if (!is.list(covariate_rasters_list) || length(covariate_rasters_list) < 1L) {
    stop(where, ": `covariate_rasters_list` must be a non-empty list when provided.",
         call. = FALSE)
  }

  base_names <- terra::names(covariate_rasters_list[[1L]])
  base_len <- length(base_names)

  for (ii in seq_along(covariate_rasters_list)) {
    nm  <- terra::names(covariate_rasters_list[[ii]])
    len <- length(nm)
    same_len   <- identical(len, base_len)
    same_names <- identical(nm, base_names)
    if (!(same_len && same_names)) {
      exp_str <- if (base_len == 0L) "<none>" else paste(base_names, collapse = ", ")
      got_str <- if (len      == 0L) "<none>" else paste(nm,         collapse = ", ")
      stop(
        paste0(
          where, ": when `time_varying_betas = TRUE`, all time slices must have ",
          "identical covariate layers (same names and order).\n",
          "Mismatch at time index ", ii, ".\n",
          "Expected: ", exp_str, "\n",
          "Got     : ", got_str
        ),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}
