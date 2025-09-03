# Simplified helper data generation for multi-map disaggregation testing
# Follows original disaggregation pattern but handles multi-map complexity

#' Generate basic multi-map test data
#'
#' Creates synthetic spatial data across multiple time points
#' Keeps complexity where needed for multi-map functionality
#'
#' @param n_times Number of time points
#' @param n_polygon_per_side Number of polygons per side (creates n^2 total)
#' @param n_pixels_per_side Number of pixels per side for rasters
#' @param add_categorical Whether to include categorical covariates
#' @param add_nas Whether to inject NAs for testing
generate_mmap_test_data <- function(n_times = 3,
                                    n_polygon_per_side = 5,
                                    n_pixels_per_side = 10,
                                    add_categorical = FALSE,
                                    add_nas = FALSE) {

  n_polygons <- n_polygon_per_side * n_polygon_per_side

  # Create base polygon grid (same across time)
  polygons <- vector("list", n_polygons)
  for(i in seq_len(n_polygons)) {
    row <- ceiling(i / n_polygon_per_side)
    col <- ifelse(i %% n_polygon_per_side != 0, i %% n_polygon_per_side, n_polygon_per_side)
    xmin <- 2 * (col - 1)
    xmax <- 2 * col
    ymin <- 2 * (row - 1)
    ymax <- 2 * row
    polygons[[i]] <- list(cbind(c(xmin, xmax, xmax, xmin, xmin),
                                c(ymax, ymax, ymin, ymin, ymax)))
  }
  polys <- lapply(polygons, sf::st_polygon)

  # Generate data for each time point
  polygon_shapefile_list <- vector("list", n_times)
  covariate_rasters_list <- vector("list", n_times)
  aggregation_rasters_list <- vector("list", n_times)

  for(t in seq_len(n_times)) {
    # Simple time-varying response data
    response_values <- runif(n_polygons, min = t*10, max = t*100)

    # Add NAs if requested
    if(add_nas && t == 2) {
      response_values[1] <- NA
    }

    response_df <- data.frame(
      area_id = seq_len(n_polygons),
      response = response_values
    )
    polygon_shapefile_list[[t]] <- sf::st_sf(response_df, geometry = polys)

    # Create covariate rasters (simple time-varying pattern)
    r1 <- terra::rast(ncol = n_pixels_per_side, nrow = n_pixels_per_side)
    terra::ext(r1) <- terra::ext(polygon_shapefile_list[[t]])
    r1[] <- sapply(seq_len(terra::ncell(r1)), function(x) {
      rnorm(1, ifelse(x %% n_pixels_per_side != 0, x %% n_pixels_per_side, n_pixels_per_side) + t, 3)
    })

    r2 <- terra::rast(ncol = n_pixels_per_side, nrow = n_pixels_per_side)
    terra::ext(r2) <- terra::ext(polygon_shapefile_list[[t]])
    r2[] <- sapply(seq_len(terra::ncell(r2)), function(x) {
      rnorm(1, ceiling(x / n_pixels_per_side) + t * 2, 3)
    })

    cov_stack <- c(r1, r2)
    names(cov_stack) <- c('temp', 'precip')

    # Add categorical covariate if requested
    if(add_categorical) {
      r3 <- terra::rast(ncol = n_pixels_per_side, nrow = n_pixels_per_side)
      terra::ext(r3) <- terra::ext(polygon_shapefile_list[[t]])
      r3[] <- sample(c("urban", "rural", "forest"), terra::ncell(r3), replace = TRUE)
      names(r3) <- "landuse"
      cov_stack <- c(cov_stack, r3)
    }

    covariate_rasters_list[[t]] <- cov_stack

    # Simple aggregation raster (population-like)
    agg_r <- terra::rast(ncol = n_pixels_per_side, nrow = n_pixels_per_side)
    terra::ext(agg_r) <- terra::ext(polygon_shapefile_list[[t]])
    agg_r[] <- rpois(terra::ncell(agg_r), lambda = 10 + t)
    aggregation_rasters_list[[t]] <- agg_r
  }

  return(list(
    polygon_shapefile_list = polygon_shapefile_list,
    covariate_rasters_list = covariate_rasters_list,
    aggregation_rasters_list = aggregation_rasters_list,
    n_times = n_times,
    n_polygons = n_polygons
  ))
}

#' Create binomial test data variant
create_binomial_variant <- function(base_data) {
  result <- base_data

  for(t in seq_len(base_data$n_times)) {
    N <- floor(runif(base_data$n_polygons, min = 10, max = 100))
    max_response <- max(base_data$polygon_shapefile_list[[t]]$response, na.rm = TRUE)
    proportions <- pmin(base_data$polygon_shapefile_list[[t]]$response / max_response, 0.8)

    result$polygon_shapefile_list[[t]]$response <- rbinom(base_data$n_polygons, N, proportions)
    result$polygon_shapefile_list[[t]]$sample_size <- N
  }

  return(result)
}

# Generate standard test datasets (like original disaggregation helper_data.R)
test_data_mmap <- generate_mmap_test_data()
test_data_mmap_categorical <- generate_mmap_test_data(add_categorical = TRUE)
test_data_mmap_nas <- generate_mmap_test_data(add_nas = TRUE)
test_data_mmap_binomial <- create_binomial_variant(test_data_mmap)
