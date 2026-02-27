#' Summary function for disag_data_mmap objects
#'
#' @description
#' Prints counts of time points, polygons, pixels, per-time largest/smallest polygon
#' sizes, number of covariates and their summaries and a mesh summary
#'
#' @param object A 'disag_data_mmap' object (from 'prepare_data_mmap()').
#' @param ... Additional arguments (unused).
#' @return Invisibly returns a list with components:
#'   - 'n_times', 'n_polygons', 'n_pixels'
#'   - 'per_time': data.frame with 'time', 'min_pixels', 'max_pixels'
#'   - 'n_covariates', 'covariate_summaries' (named list of summaries)
#'   - 'mesh_nodes', 'mesh_triangles'
#' @method summary disag_data_mmap
#' @export
summary.disag_data_mmap <- function(object, ...) {
  # 1. Basic counts
  n_times <- length(object$time_points)
  n_polygons <- nrow(object$polygon_data)
  n_pixels <- nrow(object$covariate_data)

  cat("Disaggregation data (multi-map) summary\n")
  cat("========================================\n")
  cat(sprintf("Time points: %d\n", n_times))
  cat(sprintf("Total polygons: %d\n", n_polygons))
  cat(sprintf("Total pixels: %d\n\n", n_pixels))

  # 2. Per-time polygon size
  per_time <- do.call(rbind, lapply(seq_len(n_times), function(t) {
    pix_idx <- object$covariate_data$time == t
    ct <- table(object$covariate_data$poly_local_id[pix_idx])
    data.frame(
      time       = object$time_points[t],
      min_pixels = min(ct),
      max_pixels = max(ct)
    )
  }))
  cat("Per-time polygon sizes (pixels):\n")
  print(per_time)
  cat("\n")

  # 3. Covariate summaries
  cov_rasts <- object$covariate_rasters_list
  if (is.null(cov_rasts) ||
    length(cov_rasts) == 0 ||
    is.null(cov_rasts[[1]]) ||
    length(names(cov_rasts[[1]])) == 0) {
    n_covariates <- 0
    cov_summaries <- list()
    cat("No covariates in data.\n\n")
  } else {
    cov_names <- names(cov_rasts[[1]])
    n_covariates <- length(cov_names)
    cat(sprintf("Covariates: %d\n", n_covariates))
    cov_df <- object$covariate_data[, cov_names, drop = FALSE]
    cov_summaries <- lapply(cov_df, summary)
    cat("\nCovariate summaries:\n")
    print(cov_summaries)
    # cat("\n")
  }

  # 4. Mesh summaries
  mesh <- object$mesh
  mesh_nodes <- if (!is.null(mesh)) nrow(mesh$loc) else NA_integer_
  mesh_tri <- if (!is.null(mesh)) nrow(mesh$graph$tv) else NA_integer_
  cat("Mesh summary:\n")
  cat(sprintf("  Nodes: %d\n", mesh_nodes))
  cat(sprintf("  Triangles: %d\n\n", mesh_tri))


  # Assemble return value
  out <- list(
    n_times = n_times,
    n_polygons = n_polygons,
    n_pixels = n_pixels,
    per_time = per_time,
    n_covariates = n_covariates,
    covariate_summaries = cov_summaries,
    mesh_nodes = mesh_nodes,
    mesh_triangles = mesh_tri
  )
  return(invisible(out))
}


#' Print method for 'disag_data_mmap' objects
#'
#' @description
#' Displays a brief overview of a multi-map disaggregation dataset:
#' number of time points, total polygons, and total pixels.
#'
#' @param x A 'disag_data_mmap' object.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns the original 'disag_data_mmap' object.
#' @method print disag_data_mmap
#' @export
print.disag_data_mmap <- function(x, ...) {
  n_times <- length(x$time_points)
  n_polygons <- nrow(x$polygon_data)
  n_pixels <- nrow(x$covariate_data)

  cat("Disaggregation data (multi-map) info\n")
  cat("=====================================\n")
  cat(sprintf("Time points: %d\n", n_times))
  cat(sprintf("Total polygons: %d\n", n_polygons))
  cat(sprintf("Total pixels: %d\n\n", n_pixels))

  cat("Use `summary(...)` for more details.\n")
  return(invisible(x))
}


#' Print method for 'disag_model_mmap_aghq' objects
#'
#' @description
#' Displays a brief overview of a multi-map disaggregation model:
#' model family, link function, and components included.
#'
#' @param x A 'disag_model_mmap_aghq' object.
#' @param ... Additional arguments (not used).
#' @param max_print Maximum number of random effects details to print.
#'
#' @return Invisibly returns the original 'disag_model_mmap_aghq' object.
#' @method print disag_model_mmap_aghq
#' @export
print.disag_model_mmap_aghq <- function(x, ..., max_print = 30) {
  # Check if x is a valid disag_model_mmap_aghq object
  if (!inherits(x, "disag_model_mmap_aghq")) {
    stop("Object must be of class 'disag_model_mmap_aghq'")
  }

  # Extract model setup information with error handling
  model_info <- tryCatch({
    list(
      family = x$model_setup$family,
      link = x$model_setup$link,
      field = x$model_setup$field,
      iid = x$model_setup$iid
    )
  }, error = function(e) {
    warning("Could not extract model setup information: ", e$message)
    list(
      family = "unknown",
      link = "unknown",
      field = NA,
      iid = NA
    )
  })

  # Print header and basic model information
  cat("Disaggregation model (multi-map) fit with AGHQ\n")
  cat("==============================================\n")
  cat(sprintf("Family: %s\n", model_info$family))
  cat(sprintf("Link function: %s\n", model_info$link))

  # Handle potential NA values in field and iid
  if (is.na(model_info$field)) {
    cat("Spatial field included: Unknown\n")
  } else {
    cat(sprintf("Spatial field included: %s\n", ifelse(model_info$field, "Yes", "No")))
  }

  if (is.na(model_info$iid)) {
    cat("IID effects included: Unknown\n")
  } else {
    cat(sprintf("IID effects included: %s\n", ifelse(model_info$iid, "Yes", "No")))
  }

  # Extract fixed effects information with robust error handling
  fixed_effects <- tryCatch({
    # Check if aghq_model exists
    if (is.null(x$aghq_model)) {
      warning("aghq_model component is NULL")
      return(character(0))
    }

    # Check if mode exists
    if (is.null(x$aghq_model$mode)) {
      warning("aghq_model$mode component is NULL")
      return(character(0))
    }

    # Get parameter names from the mode vector
    mode <- x$aghq_model$mode
    # Filter out random effects (nodemean and iideffect)
    fixed_names <- names(mode)[!grepl("^nodemean|^iideffect", names(mode))]
    fixed_names
  }, error = function(e) {
    warning("Error extracting fixed effects: ", e$message)
    character(0)  # Return empty vector if there's an error
  })

  cat(sprintf("Quadrature Points: %s\n", x$aghq_model$normalized_posterior$grid$level[[1]]))

  # Print fixed effects information
  if (length(fixed_effects) > 0) {
    cat(sprintf("\nFixed effects parameters: %d\n", length(fixed_effects)))
    cat("Parameter names: ", paste(fixed_effects, collapse = ", "), "\n")
  } else {
    cat("\nNo fixed effects parameters found.\n")
  }

  # Extract random effects information with robust error handling
  random_effects <- tryCatch({
    # Check if aghq_model exists
    if (is.null(x$aghq_model)) {
      warning("aghq_model component is NULL")
      return(character(0))
    }

    # Check if mode exists
    if (is.null(x$aghq_model$mode)) {
      warning("aghq_model$mode component is NULL")
      return(character(0))
    }

    # Get parameter names from the mode vector
    mode <- x$aghq_model$mode
    # Filter for random effects (nodemean and iideffect)
    random_names <- names(mode)[grepl("^nodemean|^iideffect", names(mode))]
    random_names
  }, error = function(e) {
    warning("Error extracting random effects: ", e$message)
    character(0)  # Return empty vector if there's an error
  })

  # Print random effects information
  n_random <- length(random_effects)
  if (n_random > 0) {
    cat(sprintf("\nRandom effects: %d\n", n_random))

    # Count by type with error handling
    n_nodemean <- sum(grepl("^nodemean", random_effects))
    n_iideffect <- sum(grepl("^iideffect", random_effects))

    if (n_nodemean > 0) {
      cat(sprintf("  Spatial field nodes: %d\n", n_nodemean))
    }

    if (n_iideffect > 0) {
      cat(sprintf("  IID effects: %d\n", n_iideffect))
    }

    # Print warning if too many random effects to display
    if (n_random > max_print) {
      cat(sprintf("\nThere are %d random effects, but max_print = %d, so not showing details.\n",
                  n_random, max_print))
      cat(sprintf("Set max_print higher than %d if you would like to see random effects details.\n",
                  n_random))
    }
  } else {
    cat("\nNo random effects found.\n")
  }

  # Add message directing to summary
  cat("\nUse `summary(...)` for more detailed information about the model fit.\n")

  return(invisible(x))
}


#' Summary method for 'disag_model_mmap_aghq' objects (direct approach)
#'
#' @description
#' Creates a simplified summary of a multi-map disaggregation model fit with AGHQ,
#' directly using the AGHQ model's summary information.
#'
#' @param object A 'disag_model_mmap_aghq' object.
#' @param ... Additional arguments (not used).
#'
#' @return An object of class 'summary.disag_model_mmap_aghq' containing the summary information.
#' @method summary disag_model_mmap_aghq
#' @export
summary.disag_model_mmap_aghq <- function(object, ...) {
  # Check if object is a valid disag_model_mmap_aghq object
  if (!inherits(object, "disag_model_mmap_aghq")) {
    stop("Object must be of class 'disag_model_mmap_aghq'")
  }

  # Extract model setup information
  model_info <- list(
    family = object$model_setup$family,
    link = object$model_setup$link,
    field = object$model_setup$field,
    iid = object$model_setup$iid
  )

  # Get AGHQ summary directly
  aghq_summary <- NULL
  if (!is.null(object$aghq_model)) {
    # Store the direct summary object for later use
    aghq_summary <- summary(object$aghq_model)
    # Normalize fixed-effect names in AGHQ summary output (especially summarytable rownames)
    # so shared/time-varying slopes reflect training covariate names.
    coef_meta <- tryCatch(object$model_setup$coef_meta, error = function(e) NULL)
    tv_flag <- isTRUE(tryCatch(object$model_setup$time_varying_betas, error = function(e) FALSE))
    if (!is.null(coef_meta)) {
      if (!is.null(aghq_summary$mode) && !is.null(names(aghq_summary$mode))) {
        names(aghq_summary$mode) <- normalize_fixed_names(
          names(aghq_summary$mode),
          coef_meta = coef_meta,
          time_varying_betas = tv_flag
        )
      }
      if (!is.null(aghq_summary$summarytable) && !is.null(rownames(aghq_summary$summarytable))) {
        rownames(aghq_summary$summarytable) <- normalize_fixed_names(
          rownames(aghq_summary$summarytable),
          coef_meta = coef_meta,
          time_varying_betas = tv_flag
        )
      }
    }
    quad_points <- object$aghq_model$normalized_posterior$grid$level[[1]]
  }

  # Create the summary object
  out <- list(
    model_info = model_info,
    aghq_summary = aghq_summary,
    quad_points = quad_points
  )

  # Set the class
  class(out) <- c("summary.disag_model_mmap_aghq", "list")

  return(out)
}

#' Print method for 'summary.disag_model_mmap_aghq' objects (direct approach)
#'
#' @description
#' Displays the summary information for a multi-map disaggregation model
#' in a well-formatted way, directly using the AGHQ model's summary information.
#'
#' @param x A 'summary.disag_model_mmap_aghq' object.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the original summary object.
#' @method print summary.disag_model_mmap_aghq
#' @export
print.summary.disag_model_mmap_aghq <- function(x, ...) {
  # Check if x is a valid summary.disag_model_mmap_aghq object
  if (!inherits(x, "summary.disag_model_mmap_aghq")) {
    stop("Object must be of class 'summary.disag_model_mmap_aghq'")
  }

  # Print header and basic model information
  cat("Summary of disaggregation model (multi-map) fit with AGHQ\n")
  cat("=======================================================\n")
  cat(sprintf("Family: %s\n", x$model_info$family))
  cat(sprintf("Link function: %s\n", x$model_info$link))
  cat(sprintf("Spatial field included: %s\n", ifelse(x$model_info$field, "Yes", "No")))
  cat(sprintf("IID effects included: %s\n", ifelse(x$model_info$iid, "Yes", "No")))
  cat(sprintf("Quadrature Points: %s\n", x$quad_points))


  # Print AGHQ summary if available
  if (!is.null(x$aghq_summary)) {
    cat("\nParameter estimates:\n")
    cat("------------------\n")

    # Print the summary table directly
    print(x$aghq_summary)
  } else {
    cat("\nNo AGHQ summary information available.\n")
    cat("Try using summary(object$aghq_model) directly for more information.\n")
  }

  # Return invisibly
  return(invisible(x))
}
