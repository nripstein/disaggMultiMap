#' Summary function for disag_data_mmap objects
#'
#' @description
#' Prints counts of time points, polygons, pixels, per-time largest/smallest polygon
#' sizes, number of covariates and their summaries and a mesh summary
#'
#' @param object A 'disag_data_mmap' object (from 'prepare_data_mmap()').
#' @return Invisibly returns a list with components:
#'   - 'n_times', 'n_polygons', 'n_pixels'
#'   - 'per_time': data.frame with 'time', 'min_pixels', 'max_pixels'
#'   - 'n_covariates', 'covariate_summaries' (named list of summaries)
#'   - 'mesh_nodes', 'mesh_triangles'
#' @method summary disag_data_mmap
#' @export
summary.disag_data_mmap <- function(object) {
  # 1. Basic counts
  n_times    <- length(object$time_points)
  n_polygons <- nrow(object$polygon_data)
  n_pixels   <- nrow(object$covariate_data)

  cat("Disaggregation data (multi‐map) summary\n")
  cat("========================================\n")
  cat(sprintf("Time points: %d\n", n_times))
  cat(sprintf("Total polygons: %d\n", n_polygons))
  cat(sprintf("Total pixels: %d\n\n", n_pixels))

  # 2. Per‐time polygon size
  per_time <- do.call(rbind, lapply(seq_len(n_times), function(t) {
    pix_idx <- object$covariate_data$time == t
    ct      <- table(object$covariate_data$poly_local_id[pix_idx])
    data.frame(
      time       = object$time_points[t],
      min_pixels = min(ct),
      max_pixels = max(ct)
    )
  }))
  cat("Per‐time polygon sizes (pixels):\n")
  print(per_time)
  cat("\n")

  # 3. Covariate summaries
  cov_rasts <- object$covariate_rasters_list
  if (is.null(cov_rasts) ||
      length(cov_rasts)==0 ||
      is.null(cov_rasts[[1]]) ||
      length(names(cov_rasts[[1]]))==0) {
    n_covariates     <- 0
    cov_summaries    <- list()
    cat("No covariates in data.\n\n")
  } else {
    cov_names        <- names(cov_rasts[[1]])
    n_covariates     <- length(cov_names)
    cat(sprintf("Covariates: %d\n", n_covariates))
    cov_df           <- object$covariate_data[, cov_names, drop = FALSE]
    cov_summaries    <- lapply(cov_df, summary)
    cat("\nCovariate summaries:\n")
    print(cov_summaries)
    # cat("\n")
  }

  # 4. Mesh summaries
  mesh         <- object$mesh
  mesh_nodes   <- if (!is.null(mesh)) nrow(mesh$loc) else NA_integer_
  mesh_tri     <- if (!is.null(mesh)) nrow(mesh$graph$tv) else NA_integer_
  cat("Mesh summary:\n")
  cat(sprintf("  Nodes: %d\n", mesh_nodes))
  cat(sprintf("  Triangles: %d\n\n", mesh_tri))


  # Assemble return value
  out <- list(
    n_times            = n_times,
    n_polygons         = n_polygons,
    n_pixels           = n_pixels,
    per_time           = per_time,
    n_covariates       = n_covariates,
    covariate_summaries= cov_summaries,
    mesh_nodes         = mesh_nodes,
    mesh_triangles     = mesh_tri
  )
  return(invisible(out))
}


#' Print method for 'disag_data_mmap' objects
#'
#' @description
#' Displays a brief overview of a multi-map disaggregation dataset:
#' number of time points, total polygons, and total pixels.
#'
#' @param object A 'disag_data_mmap' object.
#'
#' @return Invisibly returns the original 'disag_data_mmap' object.
#' @method print disag_data_mmap
#' @export
print.disag_data_mmap <- function(object) {
  n_times    <- length(object$time_points)
  n_polygons <- nrow(object$polygon_data)
  n_pixels   <- nrow(object$covariate_data)

  cat("Disaggregation data (multi‐map) info\n")
  cat("====================================\n")
  cat(sprintf("Time points: %d\n", n_times))
  cat(sprintf("Total polygons: %d\n", n_polygons))
  cat(sprintf("Total pixels: %d\n\n", n_pixels))

  cat("Use `summary(...)` for more details.\n")
  return(invisible(object))
}


#' Summarize model fit for multi-map disaggregation model fit with AGHQ using AGHQ's builtin summary
#'
#' @param object
#'
#' @returns
#' @export
print.disag_model_mmap_aghq <- function(object) {
  cat(summary(object$aghq_model)$summarytable, "\n")
  return(invisible(object))
}
