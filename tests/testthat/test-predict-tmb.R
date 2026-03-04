test_that("predict.disag_model_mmap_tmb returns expected structure", {
  bundle <- get_core_field_fit()
  fit <- bundle$fit
  n_times <- length(fit$data$time_points)

  pred <- suppressMessages(predict(fit, N = 2, predict_iid = FALSE))

  expect_s3_class(pred, "disag_prediction_mmap")
  expect_equal(names(pred), c("mean_prediction", "uncertainty_prediction"))

  expect_equal(names(pred$mean_prediction), c("prediction", "field", "iid", "covariates"))
  expect_true(inherits(pred$mean_prediction$prediction, "SpatRaster"))
  expect_true(inherits(pred$mean_prediction$field, "SpatRaster"))
  expect_true(is.null(pred$mean_prediction$iid))
  expect_true(inherits(pred$mean_prediction$covariates, "SpatRaster"))
  expect_equal(terra::nlyr(pred$mean_prediction$prediction), n_times)
  expect_equal(terra::nlyr(pred$mean_prediction$field), n_times)
  expect_equal(terra::nlyr(pred$mean_prediction$covariates), n_times)

  expect_equal(names(pred$uncertainty_prediction), c("realisations", "predictions_ci"))
  expect_true(is.list(pred$uncertainty_prediction$realisations))
  expect_equal(length(pred$uncertainty_prediction$realisations), n_times)
  for (ii in seq_len(n_times)) {
    expect_true(inherits(pred$uncertainty_prediction$realisations[[ii]], "SpatRaster"))
    expect_equal(terra::nlyr(pred$uncertainty_prediction$realisations[[ii]]), 2)
  }

  expect_equal(names(pred$uncertainty_prediction$predictions_ci), c("lower", "upper"))
  expect_true(inherits(pred$uncertainty_prediction$predictions_ci$lower, "SpatRaster"))
  expect_true(inherits(pred$uncertainty_prediction$predictions_ci$upper, "SpatRaster"))
  expect_equal(terra::nlyr(pred$uncertainty_prediction$predictions_ci$lower), n_times)
  expect_equal(terra::nlyr(pred$uncertainty_prediction$predictions_ci$upper), n_times)
})

test_that("predict.disag_model_mmap_tmb handles no-field/no-iid model", {
  bundle <- suppressWarnings(get_core_nofield_fit())
  fit <- bundle$fit
  n_times <- length(fit$data$time_points)

  pred <- suppressMessages(predict(fit, N = 2, predict_iid = FALSE))

  expect_true(inherits(pred$mean_prediction$prediction, "SpatRaster"))
  expect_true(is.null(pred$mean_prediction$field))
  expect_true(is.null(pred$mean_prediction$iid))
  expect_true(inherits(pred$mean_prediction$covariates, "SpatRaster"))
  expect_equal(terra::nlyr(pred$mean_prediction$prediction), n_times)

  expect_true(is.list(pred$uncertainty_prediction$realisations))
  for (ii in seq_len(n_times)) {
    expect_equal(terra::nlyr(pred$uncertainty_prediction$realisations[[ii]]), 2)
  }
  expect_equal(terra::nlyr(pred$uncertainty_prediction$predictions_ci$lower), n_times)
  expect_equal(terra::nlyr(pred$uncertainty_prediction$predictions_ci$upper), n_times)
})

test_that("predict.disag_model_mmap_tmb works with shared random-betas mode", {
  bundle <- suppressWarnings(get_cached_tmb_fit(
    name = "pred_shared_random_betas",
    seed = 14L,
    iterations = 20,
    family = "poisson",
    link = "log",
    field = FALSE,
    iid = FALSE,
    time_varying_betas = FALSE,
    fixed_effect_betas = FALSE
  ))
  fit <- bundle$fit
  n_times <- length(fit$data$time_points)

  pred <- suppressMessages(predict(fit, N = 2, predict_iid = FALSE))

  expect_s3_class(pred, "disag_prediction_mmap")
  expect_equal(terra::nlyr(pred$mean_prediction$prediction), n_times)
  expect_equal(length(pred$uncertainty_prediction$realisations), n_times)
  expect_equal(terra::nlyr(pred$uncertainty_prediction$predictions_ci$lower), n_times)
  expect_equal(terra::nlyr(pred$uncertainty_prediction$predictions_ci$upper), n_times)
})

test_that("predict.disag_model_mmap_tmb works with field + shared random-betas mode", {
  bundle <- suppressWarnings(get_cached_tmb_fit(
    name = "pred_shared_random_betas_field",
    seed = 17L,
    iterations = 60,
    family = "poisson",
    link = "log",
    field = TRUE,
    iid = FALSE,
    time_varying_betas = FALSE,
    fixed_effect_betas = FALSE
  ))
  fit <- bundle$fit
  n_times <- length(fit$data$time_points)

  pred <- suppressMessages(predict(fit, N = 2, predict_iid = FALSE))

  expect_s3_class(pred, "disag_prediction_mmap")
  expect_equal(terra::nlyr(pred$mean_prediction$prediction), n_times)
  expect_true(inherits(pred$mean_prediction$field, "SpatRaster"))
  expect_equal(terra::nlyr(pred$mean_prediction$field), n_times)
  expect_equal(length(pred$uncertainty_prediction$realisations), n_times)
  expect_equal(terra::nlyr(pred$uncertainty_prediction$predictions_ci$lower), n_times)
  expect_equal(terra::nlyr(pred$uncertainty_prediction$predictions_ci$upper), n_times)
})

test_that("predict.disag_model_mmap_tmb works with time-varying random-betas mode", {
  bundle <- suppressWarnings(get_cached_tmb_fit(
    name = "pred_tv_random_betas",
    seed = 15L,
    iterations = 20,
    family = "poisson",
    link = "log",
    field = FALSE,
    iid = FALSE,
    time_varying_betas = TRUE,
    fixed_effect_betas = FALSE
  ))
  fit <- bundle$fit
  n_times <- length(fit$data$time_points)

  pred <- suppressMessages(predict(fit, N = 2, predict_iid = FALSE))

  expect_s3_class(pred, "disag_prediction_mmap")
  expect_equal(terra::nlyr(pred$mean_prediction$prediction), n_times)
  expect_equal(length(pred$uncertainty_prediction$realisations), n_times)
  expect_equal(terra::nlyr(pred$uncertainty_prediction$predictions_ci$lower), n_times)
  expect_equal(terra::nlyr(pred$uncertainty_prediction$predictions_ci$upper), n_times)
})
