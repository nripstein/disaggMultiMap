skip_if_mcmc_opted_out <- function() {
  run_mcmc <- !tolower(Sys.getenv("RUN_MCMC_TESTS", "false")) %in% c("0", "false", "no")
  skip_if_not(run_mcmc, message = "Set RUN_MCMC_TESTS=true to run MCMC integration tests.")
  skip_if_not_installed("tmbstan")
}


test_that("disag_model_mmap MCMC returns expected fit contract (gated)", {
  skip_if_mcmc_opted_out()
  bundle <- suppressWarnings(get_cached_mcmc_fit())
  fit <- bundle$fit

  expect_s3_class(fit, "disag_model_mmap_mcmc")
  expect_s3_class(fit, "disag_model_mmap")
  expect_true(all(c("mcmc_fit", "obj", "data", "model_setup", "engine_args_used") %in% names(fit)))
  expect_true(is.list(fit$model_setup))
  expect_equal(fit$model_setup$family, "poisson")
  expect_equal(fit$model_setup$link, "log")
  expect_false(isTRUE(fit$model_setup$field))
  expect_false(isTRUE(fit$model_setup$iid))
})


test_that("predict.disag_model_mmap_mcmc returns expected structure (gated)", {
  skip_if_mcmc_opted_out()
  bundle <- suppressWarnings(get_cached_mcmc_fit())
  fit <- bundle$fit
  n_times <- length(fit$data$time_points)

  pred <- suppressWarnings(predict(fit, predict_iid = FALSE, N = 2, CI = 0.9))

  expect_s3_class(pred, "disag_prediction_mmap_mcmc")
  expect_s3_class(pred, "disag_prediction_mmap")
  expect_true(all(c("mean_prediction", "uncertainty_prediction") %in% names(pred)))

  mean_pred <- pred$mean_prediction$prediction
  expect_true(inherits(mean_pred, "SpatRaster"))
  expect_equal(terra::nlyr(mean_pred), n_times)
  expect_equal(names(mean_pred), paste0("time_", fit$data$time_points))

  ci <- pred$uncertainty_prediction$predictions_ci
  expect_true(all(c("lower", "upper") %in% names(ci)))
  expect_true(inherits(ci$lower, "SpatRaster"))
  expect_true(inherits(ci$upper, "SpatRaster"))
  expect_equal(terra::nlyr(ci$lower), n_times)
  expect_equal(terra::nlyr(ci$upper), n_times)
})
