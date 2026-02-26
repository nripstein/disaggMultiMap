test_that("AGHQ smoke fit and predict (gated)", {
  run_aghq <- tolower(Sys.getenv("RUN_AGHQ_TESTS", "false")) %in% c("1", "true", "yes")
  skip_if_not(run_aghq, message = "Set RUN_AGHQ_TESTS=true to run AGHQ tests.")

  fixture <- make_fixture_fit_tmb(seed = 41L)
  data_obj <- prepare_data_mmap(
    polygon_shapefile_list = fixture$polygon_shapefile_list,
    covariate_rasters_list = fixture$covariate_rasters_list,
    aggregation_rasters_list = fixture$aggregation_rasters_list,
    make_mesh = TRUE
  )

  fit <- disag_model_mmap(
    data = data_obj,
    engine = "AGHQ",
    family = "poisson",
    link = "log",
    aghq_k = 1,
    field = TRUE,
    iid = TRUE,
    silent = TRUE,
    optimizer = "BFGS"
  )
  expect_s3_class(fit, "disag_model_mmap_aghq")

  pred <- predict(fit, predict_iid = FALSE, N = 50)
  expect_true(inherits(pred$mean_prediction$prediction, "SpatRaster"))
})
