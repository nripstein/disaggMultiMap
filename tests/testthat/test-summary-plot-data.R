test_that("summary.disag_data_mmap returns expected structure", {
  data_obj <- get_cached_prepared_data("prep_default_mesh")

  out <- summary(data_obj)

  expect_true(is.list(out))
  expect_equal(
    names(out),
    c(
      "n_times",
      "n_polygons",
      "n_pixels",
      "per_time",
      "n_covariates",
      "covariate_summaries",
      "mesh_nodes",
      "mesh_triangles"
    )
  )
  expect_type(out$n_times, "integer")
  expect_type(out$n_polygons, "integer")
  expect_type(out$n_pixels, "integer")
  expect_true(is.data.frame(out$per_time))
  expect_true(is.list(out$covariate_summaries))
})

test_that("print.disag_data_mmap returns object invisibly", {
  data_obj <- get_cached_prepared_data("prep_default_no_mesh")

  out <- print(data_obj)
  expect_s3_class(out, "disag_data_mmap")
  expect_identical(out, data_obj)
})

test_that("plot helpers return expected classes when mesh is present", {
  data_obj <- get_cached_prepared_data("prep_default_mesh")

  p_poly <- plot_polygons(data_obj, time = 1)
  p_cov <- plot_covariate_raster(data_obj, covariate = 1, time = 1)
  p_agg <- plot_aggregation_raster(data_obj, time = 1)
  p_mesh <- plot_mesh(data_obj)
  p_sum <- plot_prepare_summary(data_obj, time = 1)
  p_s3 <- plot(data_obj, time = 1)

  expect_s3_class(p_poly, "ggplot")
  expect_s3_class(p_cov, "ggplot")
  expect_s3_class(p_agg, "ggplot")
  expect_s3_class(p_mesh, "ggplot")
  expect_true(inherits(p_sum, "gg"))
  expect_true(inherits(p_s3, "gg"))
})

test_that("plot_mesh errors when mesh is missing", {
  data_no_mesh <- get_cached_prepared_data("prep_default_no_mesh")

  expect_error(plot_mesh(data_no_mesh), "No mesh")
})
