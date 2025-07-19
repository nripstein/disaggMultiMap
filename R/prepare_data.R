#' Prepare multi-map disaggregation data
#'
#' @description
#' Given lists of polygon sf's, covariate rasters, and aggregation rasters,
#' combines them into a single 'disag_data_mmap' object ready for model fitting.
#' @param polygon_shapefile_list List of 'sf' polygon objects, one per time point.
#' @param covariate_rasters_list Optional list of 'SpatRaster' stacks; may be NULL.
#' @param aggregation_rasters_list Optional list of 'SpatRaster'; if NULL, uses uniform counts.
#' @param id_var Name of the polygon ID column in each 'sf'.
#' @param response_var Name of the response column.
#' @param categorical_covariate_baselines Named list; names are raster layers to treat as categorical, values are baseline level to drop.
#' @param sample_size_var Name of the sample-size column (for binomial models); may be NULL.
#' @param mesh_args Passed to 'build_mesh()'.
#' @param na_action Logical; if TRUE, drop or impute NAs instead of stopping.
#' @param make_mesh Logical; if TRUE, build the spatial mesh over all polygons.
#' @param verbose Logical; if TRUE, print timing info.
#' @return An object of class 'disag_data_mmap' with components:
#'   - 'polygon_data', 'covariate_data', 'aggregation_pixels', …
#' @export
prepare_data_mmap <- function(polygon_shapefile_list,
                              covariate_rasters_list = NULL,
                              aggregation_rasters_list = NULL,
                              id_var = "area_id",
                              response_var = "response",
                              categorical_covariate_baselines = NULL,
                              sample_size_var = NULL,
                              mesh_args = NULL,
                              na_action = FALSE,
                              make_mesh = TRUE,
                              verbose = FALSE) {
  start_time <- Sys.time()

  # Validate inputs
  validate_prepare_data_inputs(
    polygon_shapefile_list,
    covariate_rasters_list,
    aggregation_rasters_list,
    id_var,
    response_var,
    sample_size_var,
    make_mesh,
    categorical_covariate_baselines
  )

  n_times <- length(polygon_shapefile_list)

  # Process each time point
  per_time <- lapply(
    seq_along(polygon_shapefile_list),
    function(i) {
      prepare_time_point(
        t               = i,
        poly_sf         = polygon_shapefile_list[[i]],
        cov_rasters     = covariate_rasters_list[[i]],
        agg_raster      = aggregation_rasters_list[[i]],
        id_var          = id_var,
        response_var    = response_var,
        sample_size_var = sample_size_var,
        na_action       = na_action,
        categorical_covariate_baselines = categorical_covariate_baselines
      )
    }
  )

  # Extract each element across time:
  polygon_data_list <- lapply(per_time, `[[`, "poly_data")
  covariate_data_list <- lapply(per_time, `[[`, "cov_data")
  aggregation_pixels_list <- lapply(per_time, `[[`, "agg_pixels")
  coords_for_fit_list <- lapply(per_time, `[[`, "coords_fit")
  coords_for_prediction_list <- lapply(per_time, `[[`, "coords_pred")
  start_end_index_list <- lapply(per_time, `[[`, "start_end_index")

  # Combine data across time points
  polygon_data_combined <- do.call(rbind, polygon_data_list)
  covariate_data_combined <- do.call(rbind, covariate_data_list)
  aggregation_pixels_combined <- unlist(aggregation_pixels_list)
  coords_for_fit_combined <- do.call(rbind, coords_for_fit_list)
  # For prediction, we use the coordinates from the first time point
  coords_for_prediction_combined <- coords_for_prediction_list[[1]]

  # Build a single mesh using the union of all polygon geometries (assumes same CRS)
  mesh <- build_combined_mesh(polygon_shapefile_list, mesh_args, make_mesh)

  disag_data_mmap <- list(
    polygon_shapefile_list = polygon_shapefile_list,
    shapefile_names = list(id_var = id_var, response_var = response_var),
    covariate_rasters_list = covariate_rasters_list,
    aggregation_rasters_list = aggregation_rasters_list,
    polygon_data = polygon_data_combined,
    covariate_data = covariate_data_combined,
    aggregation_pixels = aggregation_pixels_combined,
    coords_for_fit = coords_for_fit_combined,
    coords_for_prediction = coords_for_prediction_combined,
    start_end_index = start_end_index_list, # kept as list by time point
    mesh = mesh,
    time_points = seq_len(n_times)
  )

  class(disag_data_mmap) <- c("disag_data_mmap", "list")

  if (verbose == TRUE) {
    end_time <- Sys.time()
    time_taken <- as.numeric(difftime(end_time, start_time, units = "mins"))
    message(sprintf("prepare_data_mmap() runtime: %.2f minutes\n", time_taken))
  }

  return(disag_data_mmap)
}


#### UTILS
# helper function using local polygon IDs
getStartendindex_mmap <- function(covariates, polygon_data) {
  # Ensure that covariates and polygon_data have the column 'poly_local_id' generated in prepare_data_mmap()
  if (!("poly_local_id" %in% names(covariates)) || !("poly_local_id" %in% names(polygon_data))) {
    stop("Both covariates and polygon_data must contain the 'poly_local_id' column.")
  }

  # For each polygon (by its local id), find the range (start and end index) in the covariate data
  unique_ids <- sort(unique(covariates$poly_local_id))

  startendindex <- lapply(unique_ids, function(pid) {
    idx <- which(covariates$poly_local_id == pid)
    if (length(idx) == 0) stop("No matching covariate pixels found for poly_local_id ", pid)
    range(idx)
  })

  startendindex <- do.call(rbind, startendindex)

  # Match order to polygon_data (by poly_local_id)
  # Assuming polygon_data has one row per polygon and poly_local_id is unique within the time slice.
  order_idx <- match(polygon_data$poly_local_id, unique_ids)
  startendindex <- startendindex[order_idx, , drop = FALSE]

  # Downstream TMB C++ expects 0-indexing, so subtract 1
  startendindex <- startendindex - 1L

  return(startendindex)
}

#' Validate inputs to prepare_data_mmap()
#'
#' @description
#' Check that all list-arguments are the same length (where required),
#' that required arguments have the correct type, and that
#' 'id_var', 'response_var', and 'sample_size_var' are valid strings.
#'
#' @param polygon_shapefile_list A list of 'sf' objects.
#' @param covariate_rasters_list NULL or a list of 'SpatRaster' objects.
#' @param aggregation_rasters_list NULL or a list of 'SpatRaster' objects.
#' @param id_var Character of length 1: name of polygon ID column.
#' @param response_var Character of length 1: name of response column.
#' @param sample_size_var NULL or character of length 1: sample-size column.
#' @param make_mesh Logical flag indicating whether to build a mesh.
#' @param categorical_covariate_baselines passed from prepare_data_mmap
#'
#' @return Invisibly 'TRUE' if all checks pass; otherwise stops with an error.
#' @keywords internal
validate_prepare_data_inputs <- function(polygon_shapefile_list,
                                         covariate_rasters_list,
                                         aggregation_rasters_list,
                                         id_var,
                                         response_var,
                                         sample_size_var,
                                         make_mesh,
                                         categorical_covariate_baselines = NULL) {
  # polygon list
  if (!is.list(polygon_shapefile_list) ||
    length(polygon_shapefile_list) < 1) {
    stop("`polygon_shapefile_list` must be a non-empty list of sf objects.")
  }
  if (any(!vapply(polygon_shapefile_list, inherits, logical(1), what = "sf"))) {
    stop("All elements of `polygon_shapefile_list` must inherit from class 'sf'.")
  }

  n_times <- length(polygon_shapefile_list)

  # covariate list
  if (!is.null(covariate_rasters_list)) {
    if (!is.list(covariate_rasters_list) ||
      length(covariate_rasters_list) != n_times) {
      stop("`covariate_rasters_list` must be NULL or a list of length ", n_times, ".")
    }
    if (any(!vapply(covariate_rasters_list, inherits, logical(1), what = "SpatRaster"))) {
      stop("All elements of `covariate_rasters_list` must be SpatRaster objects.")
    }
  }

  # aggregation list
  if (!is.null(aggregation_rasters_list)) {
    if (!is.list(aggregation_rasters_list) ||
      length(aggregation_rasters_list) != n_times) {
      stop("`aggregation_rasters_list` must be NULL or a list of length ", n_times, ".")
    }
    if (any(!vapply(aggregation_rasters_list, inherits, logical(1), what = "SpatRaster"))) {
      stop("All elements of `aggregation_rasters_list` must be SpatRaster objects.")
    }
  }

  # id_var / response_var / sample_size_var
  if (!is.character(id_var) || length(id_var) != 1) {
    stop("`id_var` must be a single string naming the polygon ID column.")
  }
  if (!is.character(response_var) || length(response_var) != 1) {
    stop("`response_var` must be a single string naming the response column.")
  }
  if (!is.null(sample_size_var) &&
    (!is.character(sample_size_var) || length(sample_size_var) != 1)) {
    stop("`sample_size_var` must be NULL or a single string naming the sample-size column.")
  }

  # make_mesh flag
  if (!is.logical(make_mesh) || length(make_mesh) != 1) {
    stop("`make_mesh` must be a single logical value (TRUE or FALSE).")
  }

  # categorical_baselines must be NULL or named list of character
  if (!is.null(categorical_covariate_baselines)) {
    if (!is.list(categorical_covariate_baselines) || is.null(names(categorical_covariate_baselines))) {
      stop("'categorical_covariate_baselines' must be a named list of baseline values such as 'list(landuse = 'urban', soil_type = 'clay')'.")
    }
    # ensure each name matches a covariate layer
    cov_names <- if (is.null(covariate_rasters_list)) character(0) else names(covariate_rasters_list[[1]])
    bad <- setdiff(names(categorical_covariate_baselines), cov_names)
    if (length(bad)) stop("Unknown categorical layers: ", paste(bad, collapse=", "))
  }

  return(invisible(TRUE))
}

#' Process a single time point for prepare_data_mmap()
#'
#' @description
#' Internal helper called by 'prepare_data_mmap()'.
#' For time index 't', it:
#' 1. Validates the polygon sf and rasters.
#' 2. Handles NAs in the response.
#' 3. Builds or validates the aggregation raster.
#' 4. Extracts and merges covariate + aggregation pixel data.
#' 5. Computes coordinates for mesh fitting and for prediction.
#' 6. Computes the start/end pixel indices per polygon.
#'
#' @param t Integer time-point index (used for messaging).
#' @param poly_sf An 'sf' polygon object for time 't'.
#' @param cov_rasters A 'SpatRaster' of covariates for time 't', or NULL.
#' @param agg_raster A 'SpatRaster' of aggregation weights for time 't', or NULL.
#' @param id_var Name of the polygon ID column.
#' @param response_var Name of the response column.
#' @param sample_size_var Name of the sample-size column, or NULL.
#' @param na_action Logical; if TRUE, drop/impute NAs instead of erroring.
#'
#' @return A list with elements:
#'   - 'poly_data': data.frame of polygon-level info (incl. 'poly_local_id' & 'time').
#'   - 'cov_data': data.frame of pixel-level covariates + 'poly_local_id' + 'time' + 'cell'.
#'   - 'agg_pixels': numeric vector of aggregation weights per pixel.
#'   - 'coords_fit': coords for mesh-building (only used pixels).
#'   - 'coords_pred': coords for full-extent prediction.
#'   - 'start_end_index': integer matrix of 0-indexed start/end for each polygon.
#' @keywords internal
prepare_time_point <- function(t,
                               poly_sf,
                               cov_rasters,
                               agg_raster,
                               id_var,
                               response_var,
                               sample_size_var,
                               na_action,
                               categorical_covariate_baselines) {
  #-- 1. Validate inputs --
  if (!inherits(poly_sf, "sf")) {
    stop("Time ", t, ": `polygon_shapefile_list[[t]]` must be an sf object.")
  }
  if (!is.null(cov_rasters) && !inherits(cov_rasters, "SpatRaster")) {
    stop("Time ", t, ": covariate_rasters_list[[t]] must be a SpatRaster or NULL.")
  }
  if (!is.null(agg_raster) && !inherits(agg_raster, "SpatRaster")) {
    stop("Time ", t, ": aggregation_rasters_list[[t]] must be a SpatRaster or NULL.")
  }

  #-- 2. Handle missing responses --
  resp_vals <- poly_sf[[response_var]]
  na_resp <- is.na(resp_vals)
  if (any(na_resp)) {
    if (na_action) {
      poly_sf <- poly_sf[!na_resp, ]
      message("Time ", t, ": dropped ", sum(na_resp), " polygons with NA response.")
    } else {
      stop("Time ", t, ": found NAs in response; set na_action = TRUE to drop.")
    }
  }

  #-- 3. Build or validate aggregation raster --
  if (is.null(agg_raster)) {
    # No user‐provided aggregation → uniform weights from first covariate
    if (is.null(cov_rasters)) {
      stop("Time ", t, ": no aggregation raster and no covariates to derive one.")
    }
    agg_raster <- cov_rasters[[1]]
    terra::values(agg_raster) <- rep(1, terra::ncell(agg_raster))
  }
  names(agg_raster) <- "aggregation_raster"

  #-- 4. Merge covariate & aggregation pixel data --
  # 4a) Build poly_data with local IDs
  poly_sf$poly_local_id <- seq_len(nrow(poly_sf))
  poly_sf$time <- t
  poly_data <- disaggregation::getPolygonData(
    poly_sf, id_var, response_var, sample_size_var
  )
  poly_data$poly_local_id <- poly_sf$poly_local_id
  poly_data$time <- t

  # 4b) Combine rasters for extraction
  if (is.null(cov_rasters)) {
    rast_comb <- agg_raster
    cov_names <- character(0)
  } else {
    rast_comb <- c(cov_rasters, agg_raster)
    cov_names <- names(cov_rasters)
  }

  # 4c) Extract cell values and merge with polygon info
  raw_cov <- terra::extract(
    rast_comb,
    terra::vect(poly_sf),
    cells = TRUE,
    na.rm = TRUE,
    ID = TRUE
  )
  # Temporarily tag rows to merge back
  poly_data$temp_id <- seq_len(nrow(poly_data))
  cov_data <- merge(
    raw_cov,
    poly_data[, c("temp_id", "poly_local_id", "time")],
    by.x = "ID",
    by.y = "temp_id"
  )
  cov_data$time <- t
  cov_data$temp_id <- NULL

  # 4d) Separate aggregation weights
  agg_filter <- names(cov_data) == "aggregation_raster"
  agg_pixels <- as.numeric(cov_data[, agg_filter])
  cov_data <- cov_data[, !agg_filter, drop = FALSE]

  # 4e) Handle NAs in aggregation weights
  na_agg <- is.na(agg_pixels)
  if (any(na_agg)) {
    if (na_action) {
      agg_pixels[na_agg] <- 0
      message("Time ", t, ": replaced ", sum(na_agg), " NA aggregation weights with 0.")
    } else {
      stop("Time ", t, ": found NAs in aggregation weights; set na_action = TRUE to impute.")
    }
  }

  # 4f) Handle NAs in covariates (if any exist)
  if (!is.null(cov_rasters)) {
    cov_cols <- setdiff(names(cov_data), c("poly_local_id", "cell", "time"))
    na_cov <- unlist(lapply(cov_data[, cov_cols, drop = FALSE], function(x) any(is.na(x))))
    if (any(na_cov)) {
      if (na_action) {
        for (col in cov_cols) {
          col_na <- is.na(cov_data[[col]])
          cov_data[[col]][col_na] <- stats::median(cov_data[[col]], na.rm = TRUE)
        }
        message("Time ", t, ": imputed NAs in covariates with column medians.")
      } else {
        stop("Time ", t, ": found NAs in covariates; set na_action = TRUE to impute.")
      }
    }
  }

  # 4g) Encode categorical rasters (one-hot, drop baseline)
  if (!is.null(cov_rasters) && length(categorical_covariate_baselines) > 0) {
    for (lay in intersect(names(cov_rasters), names(categorical_covariate_baselines))) {
      vals <- cov_data[[lay]]
      # derive levels via terra; fallback to unique raster values
      lvl_list <- terra::levels(cov_rasters[[lay]])
      if (!is.null(lvl_list) && nrow(lvl_list[[1]]) > 0) {
        df <- lvl_list[[1]]
        # assume first column is ID, second column is category label
        if (ncol(df) >= 2) {
          lvl_tab <- as.character(df[[2]])
        } else {
          lvl_tab <- as.character(df[[1]])
        }
      } else {
        lvl_tab <- sort(unique(terra::values(cov_rasters[[lay]])))
      }
      f <- factor(vals, levels = lvl_tab)
      # validate baseline exists
      bl <- categorical_covariate_baselines[[lay]]
      if (!bl %in% lvl_tab) {
        stop(sprintf("Layer '%s': specified baseline '%s' not found among levels: %s",
                     lay, bl, paste(lvl_tab, collapse=", ")))
      }
      # set baseline
      f <- relevel(f, ref = bl)
      # treatment contrasts: intercept + (k-1) dummies
      X <- model.matrix(~ f)
      # drop intercept column
      dm <- X[, -1, drop = FALSE]

      colnames(dm) <- sub("^f", paste0(lay, "_"), colnames(dm))
      # bind into cov_data and remove original
      cov_data <- cbind(cov_data, dm)
      cov_data[[lay]] <- NULL
    }
  }

  #-- 5. Compute coordinates for mesh & prediction --
  raster_for_coords <- if (is.null(cov_rasters)) agg_raster else cov_rasters
  coords_fit <- disaggregation:::extractCoordsForMesh(
    raster_for_coords,
    selectIds = cov_data$cell
  )
  coords_pred <- disaggregation:::extractCoordsForMesh(raster_for_coords)

  #-- 6. Compute start/end indices for each polygon --
  start_end_index <- getStartendindex_mmap(
    covariates   = cov_data,
    polygon_data = poly_data
  )

  return(list(
    poly_data         = poly_data,
    cov_data          = cov_data,
    agg_pixels        = agg_pixels,
    coords_fit        = coords_fit,
    coords_pred       = coords_pred,
    start_end_index   = start_end_index
  ))
}

#' (Internal) Build combined mesh from a list of sf polygons
#'
#' @description
#' Just in case some maps have additional polygons outside the first extent
#' @param polygon_list List of 'sf' objects (same CRS).
#' @param mesh_args Passed to 'build_mesh()'.
#' @return An 'inla.mesh' object or 'NULL' if 'make_mesh = FALSE'.
#' @keywords internal
build_combined_mesh <- function(polygon_list, mesh_args, make_mesh) {
  if (!make_mesh) {
    message("A mesh is not being built. Spatial model will not run without a mesh.")
    return(NULL)
  }
  combined_geoms <- do.call(c, lapply(polygon_list, sf::st_geometry))
  combined_sf <- sf::st_sf(
    geometry = sf::st_combine(combined_geoms),
    crs      = sf::st_crs(polygon_list[[1]])
  )
  disaggregation::build_mesh(combined_sf, mesh_args)
}
