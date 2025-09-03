test_that("prepare_data_mmap function works as expected", {

  result <- prepare_data_mmap(
    polygon_shapefile_list = test_data_mmap$polygon_shapefile_list,
    covariate_rasters_list = test_data_mmap$covariate_rasters_list,
    aggregation_rasters_list = test_data_mmap$aggregation_rasters_list
  )

  # Use modern validation function
  validate_disag_data_mmap_structure(
    result,
    expected_times = test_data_mmap$n_times,
    expected_polygons = test_data_mmap$n_polygons
  )

  # Check temporal consistency
  validate_temporal_consistency(result)
})

test_that("prepare_data_mmap function with sample size works as expected", {

  result <- prepare_data_mmap(
    polygon_shapefile_list = test_data_mmap_binomial$polygon_shapefile_list,
    covariate_rasters_list = test_data_mmap_binomial$covariate_rasters_list,
    aggregation_rasters_list = test_data_mmap_binomial$aggregation_rasters_list,
    sample_size_var = 'sample_size',
    make_mesh = FALSE
  )

  validate_disag_data_mmap_structure(result)
  expect_true(is.null(result$mesh))
  expect_equal(sum(is.na(result$polygon_data$N)), 0)
})

test_that("prepare_data_mmap function deals with NAs as expected", {

  # Should error without na_action
  expect_error(
    prepare_data_mmap(
      polygon_shapefile_list = test_data_mmap_nas$polygon_shapefile_list,
      covariate_rasters_list = test_data_mmap_nas$covariate_rasters_list,
      aggregation_rasters_list = test_data_mmap_nas$aggregation_rasters_list
    )
  )

  # Should work with na_action = TRUE
  result <- prepare_data_mmap(
    polygon_shapefile_list = test_data_mmap_nas$polygon_shapefile_list,
    covariate_rasters_list = test_data_mmap_nas$covariate_rasters_list,
    aggregation_rasters_list = test_data_mmap_nas$aggregation_rasters_list,
    na_action = TRUE,
    make_mesh = FALSE
  )

  validate_disag_data_mmap_structure(result)
  expect_equal(sum(is.na(result$polygon_data$response)), 0)
  expect_equal(sum(is.na(result$covariate_data)), 0)
  expect_equal(sum(is.na(result$aggregation_pixels)), 0)
})

test_that("prepare_data_mmap handles categorical covariates correctly", {

  result <- prepare_data_mmap(
    polygon_shapefile_list = test_data_mmap_categorical$polygon_shapefile_list,
    covariate_rasters_list = test_data_mmap_categorical$covariate_rasters_list,
    aggregation_rasters_list = test_data_mmap_categorical$aggregation_rasters_list,
    categorical_covariate_baselines = list(landuse = "urban")
  )

  validate_disag_data_mmap_structure(result)

  # Check categorical encoding
  validate_categorical_encoding(
    result$covariate_data,
    c("urban", "rural", "forest"),
    "urban"
  )
})

test_that("prepare_data_mmap validates multi-map inputs correctly", {

  # Test mismatched list lengths
  expect_error(
    prepare_data_mmap(
      polygon_shapefile_list = test_data_mmap$polygon_shapefile_list,
      covariate_rasters_list = test_data_mmap$covariate_rasters_list[1:2],
      aggregation_rasters_list = test_data_mmap$aggregation_rasters_list
    )
  )

  # Test wrong object types
  expect_error(
    prepare_data_mmap(
      polygon_shapefile_list = list("not_sf", "objects"),
      covariate_rasters_list = test_data_mmap$covariate_rasters_list,
      aggregation_rasters_list = test_data_mmap$aggregation_rasters_list
    )
  )
})

test_that("prepare_data_mmap works without covariates", {

  result <- prepare_data_mmap(
    polygon_shapefile_list = test_data_mmap$polygon_shapefile_list,
    covariate_rasters_list = NULL,
    aggregation_rasters_list = test_data_mmap$aggregation_rasters_list,
    make_mesh = FALSE
  )

  validate_disag_data_mmap_structure(result)
  expect_null(result$covariate_rasters_list)

  # Should only have ID columns and aggregation in covariate_data
  expected_cols <- c("ID", "cell", "poly_local_id", "time")
  expect_setequal(names(result$covariate_data), expected_cols)
})

test_that("getStartendindex_mmap helper function works correctly", {

  # Create simple test data
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

  # Check 0-indexing (C++ style)
  expected_matrix <- matrix(c(0, 3, 5, 2, 4, 8), nrow = 3, ncol = 2)
  expect_equal(result, expected_matrix)
})
