test_that("prepare_data_mmap basic structure works as expected", {
  result <- get_cached_prepared_data("prep_default_mesh")

  validate_disag_data_mmap_structure(
    result,
    expected_times = test_data_mmap$n_times,
    expected_polygons_total = test_data_mmap$n_polygons * test_data_mmap$n_times
  )
  validate_temporal_consistency(result)
  expect_equal(result$time_points, seq_len(test_data_mmap$n_times))
})

test_that("prepare_data_mmap with sample size works as expected", {
  result <- get_cached_prepared_data("prep_binomial_no_mesh")

  validate_disag_data_mmap_structure(
    result,
    expected_times = test_data_mmap_binomial$n_times,
    expected_polygons_total = test_data_mmap_binomial$n_polygons * test_data_mmap_binomial$n_times
  )
  expect_true(is.null(result$mesh))
  expect_equal(sum(is.na(result$polygon_data$N)), 0)
})

test_that("prepare_data_mmap handles NAs as expected", {
  expect_error(
    prepare_data_mmap(
      polygon_shapefile_list = test_data_mmap_nas$polygon_shapefile_list,
      covariate_rasters_list = test_data_mmap_nas$covariate_rasters_list,
      aggregation_rasters_list = test_data_mmap_nas$aggregation_rasters_list
    )
  )

  result <- get_cached_prepared_data("prep_na_no_mesh")

  validate_disag_data_mmap_structure(result, expected_times = test_data_mmap_nas$n_times)
  expect_true(is.null(result$mesh))
  expect_equal(sum(is.na(result$polygon_data$response)), 0)
  expect_equal(sum(is.na(result$covariate_data)), 0)
  expect_equal(sum(is.na(result$aggregation_pixels)), 0)
})

test_that("prepare_data_mmap handles categorical covariates correctly", {
  result <- get_cached_prepared_data("prep_categorical_mesh")

  validate_disag_data_mmap_structure(result, expected_times = test_data_mmap_categorical$n_times)
  expect_equal(result$categorical_covariate_baselines$landuse, "urban")
  expect_true("landuse" %in% names(result$categorical_covariate_schema))
  validate_categorical_encoding(
    covariate_data = result$covariate_data,
    levels_all = c("urban", "rural", "forest"),
    baseline = "urban",
    layer_name = "landuse"
  )
})

test_that("prepare_data_mmap validates multi-map inputs", {
  expect_error(
    prepare_data_mmap(
      polygon_shapefile_list = test_data_mmap$polygon_shapefile_list,
      covariate_rasters_list = test_data_mmap$covariate_rasters_list[1:2],
      aggregation_rasters_list = test_data_mmap$aggregation_rasters_list
    )
  )

  expect_error(
    prepare_data_mmap(
      polygon_shapefile_list = list("not_sf", "objects"),
      covariate_rasters_list = test_data_mmap$covariate_rasters_list,
      aggregation_rasters_list = test_data_mmap$aggregation_rasters_list
    )
  )
})

test_that("prepare_data_mmap works without covariates", {
  result <- get_cached_prepared_data("prep_nocov_no_mesh")

  validate_disag_data_mmap_structure(result, expected_times = test_data_mmap$n_times)
  expect_true(all(vapply(result$covariate_rasters_list, is.null, logical(1))))
  expect_true(is.null(result$mesh))
  expect_setequal(names(result$covariate_data), c("ID", "cell", "poly_local_id", "time"))
})

test_that("getStartendindex_mmap helper works correctly", {
  set.seed(42)
  covariates <- data.frame(
    poly_local_id = c(1, 1, 1, 2, 2, 3, 3, 3, 3),
    cell = 1:9,
    temp = rnorm(9)
  )

  polygon_data <- data.frame(
    poly_local_id = c(1, 2, 3),
    response = c(10, 20, 30)
  )

  result <- getStartendindex_mmap(covariates, polygon_data)

  expect_type(result, "integer")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), nrow(polygon_data))
  expect_equal(ncol(result), 2)
  expect_equal(result, matrix(c(0, 3, 5, 2, 4, 8), nrow = 3, ncol = 2))
})
