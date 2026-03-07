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
    calc_sigma_max <- stats::sd(resp / mean(resp, na.rm = TRUE), na.rm = TRUE)
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

# Build an internal level table (code + label) for one categorical raster layer.
extract_categorical_level_table <- function(raster_layer, layer_name, context = "categorical schema") {
  lvl_list <- tryCatch(terra::levels(raster_layer), error = function(e) NULL)
  if (!is.null(lvl_list) &&
      length(lvl_list) > 0L &&
      !is.null(lvl_list[[1]]) &&
      is.data.frame(lvl_list[[1]]) &&
      nrow(lvl_list[[1]]) > 0L) {
    lvl_df <- lvl_list[[1]]
    code <- lvl_df[[1]]
    if (ncol(lvl_df) >= 2L) {
      label <- lvl_df[[2]]
    } else {
      label <- lvl_df[[1]]
    }
  } else {
    vals <- sort(unique(stats::na.omit(terra::values(raster_layer))))
    if (!length(vals)) {
      stop(sprintf(
        "%s: layer '%s' has no non-missing values to define categorical levels.",
        context, layer_name
      ), call. = FALSE)
    }
    code <- vals
    label <- as.character(vals)
  }

  tbl <- data.frame(
    code = code,
    label = as.character(label),
    stringsAsFactors = FALSE
  )
  tbl <- tbl[!is.na(tbl$code) & !is.na(tbl$label), , drop = FALSE]
  tbl <- unique(tbl)
  if (!nrow(tbl)) {
    stop(sprintf(
      "%s: layer '%s' produced an empty categorical level table.",
      context, layer_name
    ), call. = FALSE)
  }
  tbl
}

# Normalize user-provided categorical baseline value to canonical label.
normalize_categorical_baseline_value <- function(level_table, baseline_value, layer_name, context = "categorical baseline") {
  baseline_chr <- as.character(baseline_value)
  labels <- level_table$label
  codes_chr <- as.character(level_table$code)

  idx <- match(baseline_chr, labels)
  if (is.na(idx)) {
    idx <- match(baseline_chr, codes_chr)
  }
  if (is.na(idx)) {
    baseline_num <- suppressWarnings(as.numeric(baseline_chr))
    code_num <- suppressWarnings(as.numeric(codes_chr))
    if (!is.na(baseline_num)) {
      hit <- which(!is.na(code_num) & code_num == baseline_num)
      if (length(hit) == 1L) {
        idx <- hit
      }
    }
  }

  if (is.na(idx)) {
    stop(sprintf(
      paste0(
        "%s: layer '%s' baseline '%s' not found.\n",
        "Available labels: %s\n",
        "Available codes : %s"
      ),
      context,
      layer_name,
      baseline_chr,
      paste(unique(labels), collapse = ", "),
      paste(unique(codes_chr), collapse = ", ")
    ), call. = FALSE)
  }

  labels[[idx]]
}

# Build global categorical schema from all time slices.
build_categorical_schema <- function(covariate_rasters_list, categorical_covariate_baselines) {
  if (is.null(categorical_covariate_baselines) || !length(categorical_covariate_baselines)) {
    return(list(schema = list(), baselines = list()))
  }

  schema <- list()
  normalized_baselines <- list()

  for (layer_name in names(categorical_covariate_baselines)) {
    layer_tables <- lapply(seq_along(covariate_rasters_list), function(i) {
      cov_i <- covariate_rasters_list[[i]]
      if (is.null(cov_i) || !inherits(cov_i, "SpatRaster")) {
        stop(sprintf(
          "prepare_data_mmap(): categorical layer '%s' requested but covariate raster for time %d is missing.",
          layer_name, i
        ), call. = FALSE)
      }
      extract_categorical_level_table(
        cov_i[[layer_name]],
        layer_name = layer_name,
        context = sprintf("prepare_data_mmap() time %d", i)
      )
    })

    merged_tbl <- do.call(rbind, layer_tables)
    merged_tbl$code_chr <- as.character(merged_tbl$code)

    by_label <- split(merged_tbl$code_chr, merged_tbl$label)
    bad_label <- names(by_label)[vapply(by_label, function(x) length(unique(x)) > 1L, logical(1))]
    if (length(bad_label)) {
      stop(sprintf(
        "prepare_data_mmap(): categorical layer '%s' has labels mapped to multiple codes: %s.",
        layer_name,
        paste(bad_label, collapse = ", ")
      ), call. = FALSE)
    }

    by_code <- split(merged_tbl$label, merged_tbl$code_chr)
    bad_code <- names(by_code)[vapply(by_code, function(x) length(unique(x)) > 1L, logical(1))]
    if (length(bad_code)) {
      stop(sprintf(
        "prepare_data_mmap(): categorical layer '%s' has codes mapped to multiple labels: %s.",
        layer_name,
        paste(bad_code, collapse = ", ")
      ), call. = FALSE)
    }

    label_order <- character(0)
    for (tbl in layer_tables) {
      for (lab in tbl$label) {
        if (!(lab %in% label_order)) {
          label_order <- c(label_order, lab)
        }
      }
    }

    first_idx <- match(label_order, merged_tbl$label)
    level_codes <- merged_tbl$code[first_idx]
    level_table <- data.frame(
      code = level_codes,
      label = label_order,
      stringsAsFactors = FALSE
    )

    baseline_label <- normalize_categorical_baseline_value(
      level_table = level_table,
      baseline_value = categorical_covariate_baselines[[layer_name]],
      layer_name = layer_name,
      context = "prepare_data_mmap()"
    )

    dummy_levels <- setdiff(label_order, baseline_label)
    dummy_names <- if (length(dummy_levels)) {
      paste0(layer_name, "_", make.names(dummy_levels))
    } else {
      character(0)
    }

    baseline_idx <- match(baseline_label, level_table$label)
    schema[[layer_name]] <- list(
      layer_name = layer_name,
      level_labels = level_table$label,
      level_codes = level_table$code,
      baseline_label = baseline_label,
      baseline_code = level_table$code[[baseline_idx]],
      dummy_levels = dummy_levels,
      dummy_names = dummy_names
    )
    normalized_baselines[[layer_name]] <- baseline_label
  }

  list(schema = schema, baselines = normalized_baselines)
}

# Map raw categorical values (codes or labels) to canonical labels from schema.
map_values_to_categorical_labels <- function(values, layer_schema, layer_name, context = "categorical mapping") {
  labels_out <- rep(NA_character_, length(values))
  if (!length(values)) return(labels_out)

  not_na <- !is.na(values)
  if (!any(not_na)) return(labels_out)

  vals_chr <- as.character(values[not_na])
  code_chr <- as.character(layer_schema$level_codes)
  labels <- layer_schema$level_labels

  idx <- match(vals_chr, code_chr)

  unmatched <- which(is.na(idx))
  if (length(unmatched)) {
    idx_label <- match(vals_chr[unmatched], labels)
    idx[unmatched[!is.na(idx_label)]] <- idx_label[!is.na(idx_label)]
  }

  unmatched <- which(is.na(idx))
  if (length(unmatched)) {
    vals_num <- suppressWarnings(as.numeric(vals_chr[unmatched]))
    code_num <- suppressWarnings(as.numeric(code_chr))
    if (length(code_num)) {
      for (uu in seq_along(unmatched)) {
        if (!is.na(vals_num[[uu]])) {
          hits <- which(!is.na(code_num) & code_num == vals_num[[uu]])
          if (length(hits) == 1L) {
            idx[unmatched[[uu]]] <- hits[[1]]
          }
        }
      }
    }
  }

  unmatched <- which(is.na(idx))
  if (length(unmatched)) {
    bad_vals <- unique(vals_chr[unmatched])
    stop(sprintf(
      paste0(
        "%s: layer '%s' contains values not present in training categorical schema: %s.\n",
        "Allowed labels: %s\n",
        "Allowed codes : %s"
      ),
      context,
      layer_name,
      paste(bad_vals, collapse = ", "),
      paste(unique(labels), collapse = ", "),
      paste(unique(code_chr), collapse = ", ")
    ), call. = FALSE)
  }

  labels_out[not_na] <- labels[idx]
  labels_out
}

# Encode one categorical vector into deterministic treatment dummies.
encode_categorical_values <- function(values, layer_schema, layer_name, context = "categorical encoding", na_to_baseline = FALSE) {
  labels <- map_values_to_categorical_labels(
    values = values,
    layer_schema = layer_schema,
    layer_name = layer_name,
    context = context
  )

  if (isTRUE(na_to_baseline)) {
    labels[is.na(labels)] <- layer_schema$baseline_label
  }

  n <- length(labels)
  k <- length(layer_schema$dummy_levels)
  if (k == 0L) {
    return(matrix(numeric(0), nrow = n, ncol = 0, dimnames = list(NULL, character(0))))
  }

  dm <- vapply(layer_schema$dummy_levels, function(lvl) {
    out <- rep(NA_real_, n)
    keep <- !is.na(labels)
    out[keep] <- as.numeric(labels[keep] == lvl)
    out
  }, numeric(n))

  if (k == 1L) {
    dm <- matrix(dm, ncol = 1L)
  }
  colnames(dm) <- layer_schema$dummy_names
  dm
}

# Expand raw categorical layers in a SpatRaster into schema-consistent dummy layers.
encode_categorical_raster_stack <- function(covariates,
                                            categorical_schema,
                                            time_index = NULL,
                                            context = "prediction") {
  if (is.null(covariates) || !inherits(covariates, "SpatRaster")) return(covariates)
  if (is.null(categorical_schema) || !length(categorical_schema)) return(covariates)

  layer_names <- names(covariates)
  schema_layers <- names(categorical_schema)
  if (!length(schema_layers)) return(covariates)

  context_tag <- if (is.null(time_index)) context else sprintf("%s (time %s)", context, time_index)

  for (lay in schema_layers) {
    sch <- categorical_schema[[lay]]
    has_raw <- lay %in% layer_names
    has_any_dummy <- length(sch$dummy_names) > 0L && any(sch$dummy_names %in% layer_names)
    has_all_dummies <- length(sch$dummy_names) > 0L && all(sch$dummy_names %in% layer_names)

    if (has_raw && has_any_dummy) {
      stop(sprintf(
        "%s: layer '%s' appears as both raw categorical and encoded dummy names.",
        context_tag, lay
      ), call. = FALSE)
    }
    if (!has_raw && has_any_dummy && !has_all_dummies) {
      stop(sprintf(
        "%s: layer '%s' has a partial dummy set. Expected all of: %s",
        context_tag, lay, paste(sch$dummy_names, collapse = ", ")
      ), call. = FALSE)
    }
  }

  out_layers <- list()
  for (lay in layer_names) {
    if (!(lay %in% schema_layers)) {
      out_layers[[length(out_layers) + 1L]] <- covariates[[lay]]
      next
    }

    sch <- categorical_schema[[lay]]
    vals <- terra::values(covariates[[lay]], mat = FALSE)
    dm <- encode_categorical_values(
      values = vals,
      layer_schema = sch,
      layer_name = lay,
      context = context_tag,
      na_to_baseline = FALSE
    )

    if (ncol(dm) > 0L) {
      idx <- rep(which(layer_names == lay)[1], ncol(dm))
      dummy_stack <- covariates[[idx]]
      names(dummy_stack) <- colnames(dm)
      terra::values(dummy_stack) <- dm
      for (jj in seq_len(ncol(dm))) {
        out_layers[[length(out_layers) + 1L]] <- dummy_stack[[jj]]
      }
    }
  }

  if (!length(out_layers)) {
    template <- covariates[[1]]
    terra::values(template) <- 0
    names(template) <- "intercept_only"
    return(template)
  }

  do.call(c, out_layers)
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
  shared_slope_pat <- "^slope(\\[[0-9]+\\]|\\.[0-9]+|[0-9]+)?$"

  if (!isTRUE(time_varying_betas)) {
    # Shared case: slope/slope1/slope.1/slope[1] -> cov_names (in order)
    idx <- which(grepl(shared_slope_pat, nm))
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
  shared_slope_pat <- "^slope(\\[[0-9]+\\]|\\.[0-9]+|[0-9]+)?$"

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
      idx_s <- which(grepl(shared_slope_pat, old_names))
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


# Build TMB beta index metadata across fixed/random/full parameter orders.
# @keywords internal
build_tmb_beta_metadata <- function(model_output) {
  ms <- model_output$model_setup
  coef_meta <- tryCatch(ms$coef_meta, error = function(e) NULL)
  if (is.null(coef_meta)) {
    stop("Internal: model_setup$coef_meta is required for TMB beta metadata.")
  }
  tv <- isTRUE(ms$time_varying_betas)
  p <- as.integer(coef_meta$p)
  Tn <- as.integer(coef_meta$n_times)
  cov_names <- coef_meta$cov_names

  fixed_order_raw <- tryCatch(names(model_output$sd_out$par.fixed), error = function(e) NULL)
  random_order_raw <- tryCatch(names(model_output$sd_out$par.random), error = function(e) NULL)
  full_order_raw <- tryCatch(names(model_output$obj$env$par), error = function(e) NULL)

  fixed_order <- if (is.null(fixed_order_raw)) {
    character(0)
  } else {
    normalize_fixed_names(fixed_order_raw, coef_meta, tv)
  }
  random_order <- if (is.null(random_order_raw)) {
    character(0)
  } else {
    canonicalize_draw_names(random_order_raw, coef_meta, tv)
  }
  full_order <- if (is.null(full_order_raw)) {
    character(0)
  } else {
    canonicalize_draw_names(full_order_raw, coef_meta, tv)
  }

  beta_source <- if (isTRUE(ms$fixed_effect_betas)) "fixed" else "random"
  beta_names <- if (tv) {
    c(
      paste0("intercept_t", seq_len(Tn)),
      as.vector(vapply(seq_len(Tn), function(t) paste0(cov_names, "_t", t), character(p)))
    )
  } else {
    c("intercept", cov_names)
  }

  source_order <- if (identical(beta_source, "fixed")) fixed_order else random_order
  idx_source_all <- match(beta_names, source_order)
  if (anyNA(idx_source_all)) {
    miss <- beta_names[is.na(idx_source_all)]
    stop(
      "Internal: could not map beta names in ", beta_source, " order. Missing: ",
      paste(miss, collapse = ", ")
    )
  }

  idx_full_all <- match(beta_names, full_order)
  if (anyNA(idx_full_all)) {
    miss <- beta_names[is.na(idx_full_all)]
    stop(
      "Internal: could not map beta names in full parameter order. Missing: ",
      paste(miss, collapse = ", ")
    )
  }

  if (tv) {
    intercept_idx <- idx_source_all[seq_len(Tn)]
    full_intercept_idx <- idx_full_all[seq_len(Tn)]
    slope_idx <- if (p > 0L) {
      matrix(idx_source_all[-seq_len(Tn)], nrow = p, ncol = Tn, byrow = FALSE)
    } else {
      matrix(integer(0), nrow = 0, ncol = Tn)
    }
    full_slope_idx <- if (p > 0L) {
      matrix(idx_full_all[-seq_len(Tn)], nrow = p, ncol = Tn, byrow = FALSE)
    } else {
      matrix(integer(0), nrow = 0, ncol = Tn)
    }
  } else {
    intercept_idx <- idx_source_all[1L]
    full_intercept_idx <- idx_full_all[1L]
    slope_idx <- if (p > 0L) idx_source_all[-1L] else integer(0)
    full_slope_idx <- if (p > 0L) idx_full_all[-1L] else integer(0)
  }

  list(
    theta_order = fixed_order,
    fixed_order = fixed_order,
    random_order = random_order,
    full_order = full_order,
    beta_index_map = list(
      source = beta_source,
      intercept_idx = intercept_idx,
      slope_idx = slope_idx,
      full_intercept_idx = full_intercept_idx,
      full_slope_idx = full_slope_idx,
      tv = tv,
      p = p,
      Tn = Tn
    )
  )
}


# Slice a raw parameter vector into a minimal list the predictor needs, by index.
# No names are used; only sizes/flags.
slice_params_tmb <- function(par_vec, model_output) {
  # Pull meta
  ms   <- model_output$model_setup
  dat  <- model_output$data
  coef_meta <- tryCatch(ms$coef_meta, error = function(e) NULL)
  p    <- if (!is.null(coef_meta) && !is.null(coef_meta$p)) {
    as.integer(coef_meta$p)
  } else if (!is.null(dat$covariate_rasters_list) && !is.null(dat$covariate_rasters_list[[1]])) {
    terra::nlyr(dat$covariate_rasters_list[[1]])
  } else {
    0L
  }
  Tn   <- if (!is.null(coef_meta) && !is.null(coef_meta$n_times)) {
    as.integer(coef_meta$n_times)
  } else {
    length(dat$time_points)
  }
  fam  <- ms$family
  fld  <- isTRUE(ms$field)
  iidf <- isTRUE(ms$iid)
  tv   <- isTRUE(ms$time_varying_betas)

  n_poly <- nrow(dat$polygon_data)

  par_vec <- as.numeric(par_vec)

  # Keep a compatibility fallback for unnamed vectors.
  slice_sequential <- function(vec) {
    idx <- 1L

    if (tv) {
      intercept_t <- vec[idx:(idx + Tn - 1L)]; idx <- idx + Tn
      slope_t_len <- p * Tn
      slope_t <- if (slope_t_len > 0L) vec[idx:(idx + slope_t_len - 1L)] else numeric(0)
      idx <- idx + slope_t_len
      intercept <- NULL
      slope <- numeric(0)
    } else {
      intercept <- vec[idx]; idx <- idx + 1L
      slope <- if (p > 0L) vec[idx:(idx + p - 1L)] else numeric(0)
      idx <- idx + p
      intercept_t <- NULL
      slope_t <- numeric(0)
    }

    log_tau_gaussian <- if (identical(fam, "gaussian")) { val <- vec[idx]; idx <- idx + 1L; val } else NULL

    iideffect_log_tau <- NULL
    iideffect <- numeric(0)
    if (identical(fam, "negbinomial")) {
      iideffect_log_tau <- vec[idx]; idx <- idx + 1L
    } else if (iidf) {
      iideffect_log_tau <- vec[idx]; idx <- idx + 1L
      iideffect <- vec[idx:(idx + n_poly - 1L)]; idx <- idx + n_poly
    }

    log_sigma <- NULL
    log_rho <- NULL
    nodemean <- numeric(0)
    if (fld) {
      log_sigma <- vec[idx]; idx <- idx + 1L
      log_rho <- vec[idx]; idx <- idx + 1L
      if (idx <= length(vec)) nodemean <- vec[idx:length(vec)]
    }

    list(
      intercept = intercept,
      slope = slope,
      intercept_t = intercept_t,
      slope_t = slope_t,
      iideffect_log_tau = iideffect_log_tau,
      iideffect = iideffect,
      log_tau_gaussian = log_tau_gaussian,
      log_sigma = log_sigma,
      log_rho = log_rho,
      nodemean = nodemean,
      p = p,
      Tn = Tn,
      tv = tv
    )
  }

  ref_names <- tryCatch(names(model_output$obj$env$par), error = function(e) NULL)
  par_names <- names(par_vec)
  has_valid_names <- !is.null(par_names) &&
    length(par_names) == length(par_vec) &&
    !anyNA(par_names) &&
    all(nzchar(par_names))

  if (!has_valid_names && !is.null(ref_names) && length(ref_names) == length(par_vec)) {
    names(par_vec) <- ref_names
    par_names <- names(par_vec)
    has_valid_names <- TRUE
  }

  if (!has_valid_names) return(slice_sequential(par_vec))

  pull_idx <- function(key) which(par_names == key)
  pull_vec <- function(key, expected_len = NULL, required = FALSE) {
    idx <- pull_idx(key)
    if (!length(idx)) {
      if (isTRUE(required)) {
        stop("Internal: required parameter block `", key, "` was not found in TMB parameter vector.")
      }
      return(NULL)
    }
    vals <- as.numeric(par_vec[idx])
    if (!is.null(expected_len) && length(vals) != expected_len) {
      stop(
        "Internal: parameter block `", key, "` has length ", length(vals),
        ", expected ", expected_len, "."
      )
    }
    vals
  }
  pull_scalar <- function(key, required = FALSE) {
    vals <- pull_vec(key, expected_len = 1L, required = required)
    if (is.null(vals)) return(NULL)
    vals[[1L]]
  }
  idx_in_range <- function(idx) {
    idx <- as.integer(idx)
    if (!length(idx)) return(TRUE)
    all(idx >= 1L & idx <= length(par_vec))
  }

  beta_map <- tryCatch(ms$beta_index_map, error = function(e) NULL)
  use_full_beta_idx <- !is.null(beta_map) &&
    !is.null(beta_map$full_intercept_idx) &&
    !is.null(beta_map$full_slope_idx) &&
    idx_in_range(beta_map$full_intercept_idx) &&
    idx_in_range(beta_map$full_slope_idx)

  if (isTRUE(tv)) {
    if (use_full_beta_idx) {
      intercept_t <- as.numeric(par_vec[as.integer(beta_map$full_intercept_idx)])
      slope_t <- if (p > 0L) as.numeric(par_vec[as.integer(as.vector(beta_map$full_slope_idx))]) else numeric(0)
    } else {
      intercept_t <- pull_vec("intercept_t", expected_len = Tn, required = TRUE)
      slope_t <- if (p > 0L) pull_vec("slope_t", expected_len = p * Tn, required = TRUE) else numeric(0)
    }
    intercept <- NULL
    slope <- numeric(0)
  } else {
    if (use_full_beta_idx) {
      intercept <- as.numeric(par_vec[as.integer(beta_map$full_intercept_idx[[1L]])])
      slope <- if (p > 0L) as.numeric(par_vec[as.integer(beta_map$full_slope_idx)]) else numeric(0)
    } else {
      intercept <- pull_scalar("intercept", required = TRUE)
      slope <- if (p > 0L) pull_vec("slope", expected_len = p, required = TRUE) else numeric(0)
    }
    intercept_t <- NULL
    slope_t <- numeric(0)
  }

  log_tau_gaussian <- if (identical(fam, "gaussian")) pull_scalar("log_tau_gaussian", required = TRUE) else NULL

  iideffect_log_tau <- NULL
  iideffect <- numeric(0)
  if (identical(fam, "negbinomial")) {
    iideffect_log_tau <- pull_scalar("iideffect_log_tau", required = TRUE)
  } else if (iidf) {
    iideffect_log_tau <- pull_scalar("iideffect_log_tau", required = TRUE)
    iideffect <- pull_vec("iideffect", expected_len = n_poly, required = TRUE)
  }

  log_sigma <- NULL
  log_rho <- NULL
  nodemean <- numeric(0)
  if (fld) {
    n_nodes <- tryCatch(nrow(dat$mesh$loc), error = function(e) NA_integer_)
    log_sigma <- pull_scalar("log_sigma", required = TRUE)
    log_rho <- pull_scalar("log_rho", required = TRUE)
    nodemean <- if (is.na(n_nodes)) {
      pull_vec("nodemean", required = TRUE)
    } else {
      pull_vec("nodemean", expected_len = n_nodes, required = TRUE)
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
