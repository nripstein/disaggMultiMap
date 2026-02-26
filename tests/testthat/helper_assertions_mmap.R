# Shared assertion helpers for multi-map tests.

validate_disag_data_mmap_structure <- function(result,
                                               expected_times = NULL,
                                               expected_polygons_total = NULL) {
  expected_names <- c(
    "polygon_shapefile_list",
    "shapefile_names",
    "covariate_rasters_list",
    "aggregation_rasters_list",
    "polygon_data",
    "covariate_data",
    "aggregation_pixels",
    "coords_for_fit",
    "coords_for_prediction",
    "start_end_index",
    "categorical_covariate_baselines",
    "categorical_covariate_schema",
    "mesh",
    "time_points"
  )

  expect_s3_class(result, "disag_data_mmap")
  expect_type(result, "list")
  expect_equal(names(result), expected_names)

  n_times <- length(result$time_points)
  if (!is.null(expected_times)) {
    expect_equal(n_times, expected_times)
  }

  expect_true(is.list(result$polygon_shapefile_list))
  expect_equal(length(result$polygon_shapefile_list), n_times)
  expect_true(all(vapply(result$polygon_shapefile_list, inherits, logical(1), what = "sf")))

  expect_true(is.list(result$covariate_rasters_list))
  expect_equal(length(result$covariate_rasters_list), n_times)

  expect_true(is.list(result$aggregation_rasters_list))
  expect_equal(length(result$aggregation_rasters_list), n_times)

  expect_true(is.list(result$start_end_index))
  expect_equal(length(result$start_end_index), n_times)
  expect_true(all(vapply(result$start_end_index, is.matrix, logical(1))))
  expect_true(all(vapply(result$start_end_index, function(x) ncol(x) == 2L, logical(1))))

  expect_true(is.data.frame(result$polygon_data))
  expect_true(is.data.frame(result$covariate_data))
  expect_true(is.matrix(result$coords_for_fit))
  expect_true(is.matrix(result$coords_for_prediction))
  expect_type(result$aggregation_pixels, "double")
  expect_equal(nrow(result$covariate_data), nrow(result$coords_for_fit))
  expect_equal(length(result$aggregation_pixels), nrow(result$covariate_data))

  expect_true(is.list(result$categorical_covariate_baselines))
  expect_true(is.list(result$categorical_covariate_schema))

  if (!is.null(expected_polygons_total)) {
    expect_equal(nrow(result$polygon_data), expected_polygons_total)
  }
}

validate_temporal_consistency <- function(result) {
  n_times <- length(result$time_points)

  expect_setequal(sort(unique(result$polygon_data$time)), result$time_points)
  expect_setequal(sort(unique(result$covariate_data$time)), result$time_points)

  for (ti in seq_len(n_times)) {
    n_poly_t <- sum(result$polygon_data$time == ti)
    idx_t <- result$start_end_index[[ti]]

    expect_equal(nrow(idx_t), n_poly_t)
    expect_true(all(idx_t[, 1] <= idx_t[, 2]))
    expect_true(all(idx_t >= 0L))

    cov_t <- result$covariate_data[result$covariate_data$time == ti, , drop = FALSE]
    expect_true(nrow(cov_t) > 0L)
  }
}

validate_categorical_encoding <- function(covariate_data,
                                          levels_all,
                                          baseline,
                                          layer_name = "landuse") {
  non_baseline <- setdiff(levels_all, baseline)
  dummy_cols <- paste0(layer_name, "_", make.names(non_baseline))

  expect_false(layer_name %in% names(covariate_data))
  expect_true(all(dummy_cols %in% names(covariate_data)))

  dm <- as.matrix(covariate_data[, dummy_cols, drop = FALSE])
  expect_false(any(is.na(dm)))
  expect_true(all(dm %in% c(0, 1)))

  # Treatment coding with dropped baseline should result in at most one active dummy.
  rs <- rowSums(dm)
  expect_true(all(rs %in% c(0, 1)))
}
