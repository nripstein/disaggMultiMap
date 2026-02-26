# Deterministic synthetic data and lightweight fit fixtures for multi-map tests.

with_temp_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(force(expr))
  }

  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(expr)
}

#' Generate basic multi-map test data
#'
#' @param n_times Number of time points
#' @param n_polygon_per_side Number of polygons per side (creates n^2 total)
#' @param n_pixels_per_side Number of pixels per side for rasters
#' @param add_categorical Whether to include categorical covariates
#' @param add_nas Whether to inject NAs for testing
#' @param seed Integer seed for deterministic generation
generate_mmap_test_data <- function(n_times = 3,
                                    n_polygon_per_side = 5,
                                    n_pixels_per_side = 10,
                                    add_categorical = FALSE,
                                    add_nas = FALSE,
                                    seed = 1L) {
  with_temp_seed(seed, {
    n_polygons <- n_polygon_per_side * n_polygon_per_side

    polygons <- vector("list", n_polygons)
    for (i in seq_len(n_polygons)) {
      row <- ceiling(i / n_polygon_per_side)
      col <- ifelse(i %% n_polygon_per_side != 0, i %% n_polygon_per_side, n_polygon_per_side)
      xmin <- 2 * (col - 1)
      xmax <- 2 * col
      ymin <- 2 * (row - 1)
      ymax <- 2 * row
      polygons[[i]] <- list(cbind(c(xmin, xmax, xmax, xmin, xmin),
                                  c(ymax, ymax, ymin, ymin, ymax)))
    }
    polys <- sf::st_sfc(lapply(polygons, sf::st_polygon), crs = 3857)

    polygon_shapefile_list <- vector("list", n_times)
    covariate_rasters_list <- vector("list", n_times)
    aggregation_rasters_list <- vector("list", n_times)

    for (t in seq_len(n_times)) {
      response_values <- runif(n_polygons, min = t * 10, max = t * 100)
      if (add_nas && t == 2) {
        response_values[1] <- NA
      }

      response_df <- data.frame(
        area_id = seq_len(n_polygons),
        response = response_values
      )
      polygon_shapefile_list[[t]] <- sf::st_sf(response_df, geometry = polys)

      r1 <- terra::rast(ncol = n_pixels_per_side, nrow = n_pixels_per_side)
      terra::ext(r1) <- terra::ext(polygon_shapefile_list[[t]])
      terra::crs(r1) <- "EPSG:3857"
      r1[] <- sapply(seq_len(terra::ncell(r1)), function(x) {
        rnorm(1, ifelse(x %% n_pixels_per_side != 0, x %% n_pixels_per_side, n_pixels_per_side) + t, 3)
      })

      r2 <- terra::rast(ncol = n_pixels_per_side, nrow = n_pixels_per_side)
      terra::ext(r2) <- terra::ext(polygon_shapefile_list[[t]])
      terra::crs(r2) <- "EPSG:3857"
      r2[] <- sapply(seq_len(terra::ncell(r2)), function(x) {
        rnorm(1, ceiling(x / n_pixels_per_side) + t * 2, 3)
      })

      cov_stack <- c(r1, r2)
      names(cov_stack) <- c("temp", "precip")

      if (add_categorical) {
        r3 <- terra::rast(ncol = n_pixels_per_side, nrow = n_pixels_per_side)
        terra::ext(r3) <- terra::ext(polygon_shapefile_list[[t]])
        terra::crs(r3) <- "EPSG:3857"
        r3[] <- sample(c("urban", "rural", "forest"), terra::ncell(r3), replace = TRUE)
        names(r3) <- "landuse"
        cov_stack <- c(cov_stack, r3)
      }

      covariate_rasters_list[[t]] <- cov_stack

      agg_r <- terra::rast(ncol = n_pixels_per_side, nrow = n_pixels_per_side)
      terra::ext(agg_r) <- terra::ext(polygon_shapefile_list[[t]])
      terra::crs(agg_r) <- "EPSG:3857"
      agg_r[] <- rpois(terra::ncell(agg_r), lambda = 10 + t)
      aggregation_rasters_list[[t]] <- agg_r
    }

    list(
      polygon_shapefile_list = polygon_shapefile_list,
      covariate_rasters_list = covariate_rasters_list,
      aggregation_rasters_list = aggregation_rasters_list,
      n_times = n_times,
      n_polygons = n_polygons
    )
  })
}

#' Create binomial test data variant
create_binomial_variant <- function(base_data, seed = 100L) {
  with_temp_seed(seed, {
    result <- base_data

    for (t in seq_len(base_data$n_times)) {
      n_vals <- floor(runif(base_data$n_polygons, min = 10, max = 100))
      max_response <- max(base_data$polygon_shapefile_list[[t]]$response, na.rm = TRUE)
      proportions <- pmin(base_data$polygon_shapefile_list[[t]]$response / max_response, 0.8)

      result$polygon_shapefile_list[[t]]$response <- rbinom(base_data$n_polygons, n_vals, proportions)
      result$polygon_shapefile_list[[t]]$sample_size <- n_vals
    }

    result
  })
}

make_fixture_prepare <- function(seed = 1L,
                                 n_times = 3,
                                 n_polygon_per_side = 4,
                                 n_pixels_per_side = 8,
                                 ...) {
  generate_mmap_test_data(
    seed = seed,
    n_times = n_times,
    n_polygon_per_side = n_polygon_per_side,
    n_pixels_per_side = n_pixels_per_side,
    ...
  )
}

make_fixture_fit_tmb <- function(seed = 10L,
                                 n_times = 3,
                                 n_polygon_per_side = 5,
                                 n_pixels_per_side = 10,
                                 ...) {
  generate_mmap_test_data(
    seed = seed,
    n_times = n_times,
    n_polygon_per_side = n_polygon_per_side,
    n_pixels_per_side = n_pixels_per_side,
    ...
  )
}

.mmap_fit_cache <- new.env(parent = emptyenv())
.mmap_prepare_cache <- new.env(parent = emptyenv())
.mmap_aghq_prepare_cache <- new.env(parent = emptyenv())
.mmap_aghq_fit_cache <- new.env(parent = emptyenv())

get_cached_prepared_data <- function(name = c(
                                       "prep_default_mesh",
                                       "prep_default_no_mesh",
                                       "prep_categorical_mesh",
                                       "prep_na_no_mesh",
                                       "prep_binomial_no_mesh",
                                       "prep_nocov_no_mesh"
                                     )) {
  name <- match.arg(name)
  if (exists(name, envir = .mmap_prepare_cache, inherits = FALSE)) {
    return(get(name, envir = .mmap_prepare_cache, inherits = FALSE))
  }

  prepared <- switch(
    name,
    prep_default_mesh = prepare_data_mmap(
      polygon_shapefile_list = test_data_mmap$polygon_shapefile_list,
      covariate_rasters_list = test_data_mmap$covariate_rasters_list,
      aggregation_rasters_list = test_data_mmap$aggregation_rasters_list,
      make_mesh = TRUE
    ),
    prep_default_no_mesh = prepare_data_mmap(
      polygon_shapefile_list = test_data_mmap$polygon_shapefile_list,
      covariate_rasters_list = test_data_mmap$covariate_rasters_list,
      aggregation_rasters_list = test_data_mmap$aggregation_rasters_list,
      make_mesh = FALSE
    ),
    prep_categorical_mesh = prepare_data_mmap(
      polygon_shapefile_list = test_data_mmap_categorical$polygon_shapefile_list,
      covariate_rasters_list = test_data_mmap_categorical$covariate_rasters_list,
      aggregation_rasters_list = test_data_mmap_categorical$aggregation_rasters_list,
      categorical_covariate_baselines = list(landuse = "urban"),
      make_mesh = TRUE
    ),
    prep_na_no_mesh = prepare_data_mmap(
      polygon_shapefile_list = test_data_mmap_nas$polygon_shapefile_list,
      covariate_rasters_list = test_data_mmap_nas$covariate_rasters_list,
      aggregation_rasters_list = test_data_mmap_nas$aggregation_rasters_list,
      na_action = TRUE,
      make_mesh = FALSE
    ),
    prep_binomial_no_mesh = prepare_data_mmap(
      polygon_shapefile_list = test_data_mmap_binomial$polygon_shapefile_list,
      covariate_rasters_list = test_data_mmap_binomial$covariate_rasters_list,
      aggregation_rasters_list = test_data_mmap_binomial$aggregation_rasters_list,
      sample_size_var = "sample_size",
      make_mesh = FALSE
    ),
    prep_nocov_no_mesh = prepare_data_mmap(
      polygon_shapefile_list = test_data_mmap$polygon_shapefile_list,
      covariate_rasters_list = NULL,
      aggregation_rasters_list = test_data_mmap$aggregation_rasters_list,
      make_mesh = FALSE
    )
  )

  assign(name, prepared, envir = .mmap_prepare_cache)
  prepared
}

get_cached_tmb_fit <- function(name = "default",
                               seed = 10L,
                               iterations = 60,
                               family = "poisson",
                               link = "log",
                               field = TRUE,
                               iid = FALSE,
                               time_varying_betas = FALSE) {
  key <- paste(name, seed, iterations, family, link, field, iid, time_varying_betas, sep = "|")
  if (exists(key, envir = .mmap_fit_cache, inherits = FALSE)) {
    return(get(key, envir = .mmap_fit_cache, inherits = FALSE))
  }

  fixture <- make_fixture_fit_tmb(seed = seed)
  prepared <- prepare_data_mmap(
    polygon_shapefile_list = fixture$polygon_shapefile_list,
    covariate_rasters_list = fixture$covariate_rasters_list,
    aggregation_rasters_list = fixture$aggregation_rasters_list,
    make_mesh = TRUE
  )

  fit <- disag_model_mmap_tmb(
    data = prepared,
    family = family,
    link = link,
    iterations = iterations,
    field = field,
    iid = iid,
    time_varying_betas = time_varying_betas,
    silent = TRUE
  )

  out <- list(fixture = fixture, data = prepared, fit = fit)
  assign(key, out, envir = .mmap_fit_cache)
  out
}

get_core_field_fit <- function() {
  get_cached_tmb_fit(
    name = "fit_core_field",
    seed = 10L,
    iterations = 60,
    family = "poisson",
    link = "log",
    field = TRUE,
    iid = FALSE,
    time_varying_betas = FALSE
  )
}

get_core_nofield_fit <- function() {
  get_cached_tmb_fit(
    name = "fit_core_nofield",
    seed = 10L,
    iterations = 20,
    family = "poisson",
    link = "log",
    field = FALSE,
    iid = FALSE,
    time_varying_betas = FALSE
  )
}

get_cached_aghq_prepared_data <- function(name = c(
                                            "aghq_small_onecov_mesh",
                                            "aghq_small_twocov_mesh"
                                          )) {
  name <- match.arg(name)
  if (exists(name, envir = .mmap_aghq_prepare_cache, inherits = FALSE)) {
    return(get(name, envir = .mmap_aghq_prepare_cache, inherits = FALSE))
  }

  fixture <- make_fixture_fit_tmb(
    seed = 51L,
    n_times = 2,
    n_polygon_per_side = 3,
    n_pixels_per_side = 6
  )

  prepared <- switch(
    name,
    aghq_small_onecov_mesh = {
      cov_one <- lapply(fixture$covariate_rasters_list, function(r) r[[1]])
      prepare_data_mmap(
        polygon_shapefile_list = fixture$polygon_shapefile_list,
        covariate_rasters_list = cov_one,
        aggregation_rasters_list = fixture$aggregation_rasters_list,
        make_mesh = TRUE
      )
    },
    aghq_small_twocov_mesh = prepare_data_mmap(
      polygon_shapefile_list = fixture$polygon_shapefile_list,
      covariate_rasters_list = fixture$covariate_rasters_list,
      aggregation_rasters_list = fixture$aggregation_rasters_list,
      make_mesh = TRUE
    )
  )

  assign(name, prepared, envir = .mmap_aghq_prepare_cache)
  prepared
}

get_cached_aghq_fit <- function(name = "aghq_small_onecov_shared",
                                aghq_k = 1,
                                optimizer = "BFGS",
                                field = TRUE,
                                iid = TRUE,
                                time_varying_betas = FALSE) {
  key <- paste(name, aghq_k, optimizer, field, iid, time_varying_betas, sep = "|")
  if (exists(key, envir = .mmap_aghq_fit_cache, inherits = FALSE)) {
    return(get(key, envir = .mmap_aghq_fit_cache, inherits = FALSE))
  }

  prepared <- switch(
    name,
    aghq_small_onecov_shared = get_cached_aghq_prepared_data("aghq_small_onecov_mesh"),
    stop("Unknown AGHQ fit fixture name: ", name)
  )

  fit <- disag_model_mmap(
    data = prepared,
    engine = "AGHQ",
    family = "poisson",
    link = "log",
    aghq_k = aghq_k,
    field = field,
    iid = iid,
    time_varying_betas = time_varying_betas,
    silent = TRUE,
    optimizer = optimizer
  )

  out <- list(data = prepared, fit = fit)
  assign(key, out, envir = .mmap_aghq_fit_cache)
  out
}

# Baseline fixtures used in prepare-data tests
test_data_mmap <- make_fixture_prepare(seed = 1L)
test_data_mmap_categorical <- make_fixture_prepare(seed = 2L, add_categorical = TRUE)
test_data_mmap_nas <- make_fixture_prepare(seed = 3L, add_nas = TRUE)
test_data_mmap_binomial <- create_binomial_variant(test_data_mmap, seed = 4L)
