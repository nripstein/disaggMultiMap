#' Detect if a raster layer is categorical
#'
#' @description
#' Determines if a raster layer should be treated as categorical based on multiple criteria:
#' 1. Explicit definition in categorical_covariate_baselines
#' 2. Presence of defined levels in the raster
#' 3. Small number of unique values
#'
#' @param raster_layer A SpatRaster layer to check
#' @param layer_name The name of the layer
#' @param categorical_baselines Named list of categorical baselines from disag_data_mmap
#' @param max_categories Maximum number of unique values to consider categorical (default = 10)
#' @return Logical indicating if the layer should be treated as categorical
#' @keywords internal
is_categorical_layer <- function(raster_layer, layer_name, categorical_baselines = NULL, max_categories = 10) {
  # Check if explicitly defined as categorical
  if (!is.null(categorical_baselines) && layer_name %in% names(categorical_baselines)) {
    return(TRUE)
  }
  
  # Check if raster has defined levels
  lvl_list <- tryCatch({
    terra::levels(raster_layer)
  }, error = function(e) {
    return(NULL)
  })
  
  if (!is.null(lvl_list) && length(lvl_list) > 0) {
    if (!is.null(lvl_list[[1]]) && is.data.frame(lvl_list[[1]]) && nrow(lvl_list[[1]]) > 0) {
      return(TRUE)
    }
  }
  
  # Check if small number of unique values (excluding NA)
  vals <- terra::values(raster_layer)
  unique_vals <- length(unique(stats::na.omit(vals)))
  return(!is.na(unique_vals) && unique_vals <= max_categories)
}

#' Get categorical levels from a raster layer
#'
#' @description
#' Extracts categorical levels from a raster layer, either from defined levels
#' or from unique values in the raster.
#'
#' @param raster_layer A SpatRaster layer
#' @return Character vector of level names
#' @keywords internal
get_categorical_levels <- function(raster_layer) {
  # Try to get levels from raster
  lvl_list <- tryCatch({
    terra::levels(raster_layer)
  }, error = function(e) {
    return(NULL)
  })
  
  if (!is.null(lvl_list) && length(lvl_list) > 0 && !is.null(lvl_list[[1]]) && is.data.frame(lvl_list[[1]]) && nrow(lvl_list[[1]]) > 0) {
    df <- lvl_list[[1]]
    # If there are at least 2 columns, assume second column is label
    if (ncol(df) >= 2) {
      return(as.character(df[[2]]))
    } else {
      return(as.character(df[[1]]))
    }
  }
  
  # Fall back to unique values
  vals <- terra::values(raster_layer)
  return(as.character(sort(unique(stats::na.omit(vals)))))
}