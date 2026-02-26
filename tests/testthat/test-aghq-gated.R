skip_if_aghq_opted_out <- function() {
  run_aghq <- !tolower(Sys.getenv("RUN_AGHQ_TESTS", "true")) %in% c("0", "false", "no")
  skip_if_not(run_aghq, message = "Set RUN_AGHQ_TESTS=false to skip heavy AGHQ tests.")
}

test_that("disag_model_mmap AGHQ returns expected fit contract (gated)", {
  skip_if_aghq_opted_out()
  bundle <- suppressWarnings(get_cached_aghq_fit("aghq_small_onecov_shared"))
  fit <- bundle$fit

  expect_s3_class(fit, "disag_model_mmap_aghq")
  expect_s3_class(fit, "disag_model_mmap")
  expect_true(all(c("aghq_model", "obj", "data", "sd_out", "model_setup") %in% names(fit)))

  expect_true(is.list(fit$model_setup))
  expect_equal(fit$model_setup$family, "poisson")
  expect_equal(fit$model_setup$link, "log")
  expect_true(isTRUE(fit$model_setup$field))
  expect_true(isTRUE(fit$model_setup$iid))
  expect_false(isTRUE(fit$model_setup$time_varying_betas))

  expect_true(is.character(fit$model_setup$theta_order))
  expect_true(length(fit$model_setup$theta_order) > 0L)
  expect_true(is.list(fit$model_setup$beta_index_map))
})

test_that("predict.disag_model_mmap_aghq returns expected structure (gated)", {
  skip_if_aghq_opted_out()
  bundle <- suppressWarnings(get_cached_aghq_fit("aghq_small_onecov_shared"))
  fit <- bundle$fit
  n_times <- length(fit$data$time_points)

  pred <- suppressWarnings(predict(fit, predict_iid = FALSE, N = 12, CI = 0.9))

  expect_s3_class(pred, "disag_prediction_mmap_aghq")
  expect_equal(names(pred), c("mean_prediction", "uncertainty_prediction"))
  expect_equal(names(pred$mean_prediction), c("prediction", "field", "iid", "covariates"))
  expect_null(pred$mean_prediction$iid)

  expect_true(is.list(pred$mean_prediction$prediction))
  expect_true(is.list(pred$mean_prediction$field))
  expect_true(is.list(pred$mean_prediction$covariates))
  expect_equal(length(pred$mean_prediction$prediction), n_times)
  expect_equal(length(pred$mean_prediction$field), n_times)
  expect_equal(length(pred$mean_prediction$covariates), n_times)

  for (ii in seq_len(n_times)) {
    expect_true(inherits(pred$mean_prediction$prediction[[ii]], "SpatRaster"))
    expect_true(inherits(pred$mean_prediction$field[[ii]], "SpatRaster"))
    expect_true(inherits(pred$mean_prediction$covariates[[ii]], "SpatRaster"))
  }

  expect_equal(names(pred$uncertainty_prediction), c("realisations", "predictions_ci"))
  expect_null(pred$uncertainty_prediction$realisations)
  expect_equal(names(pred$uncertainty_prediction$predictions_ci), c("lower", "upper"))
  expect_true(is.list(pred$uncertainty_prediction$predictions_ci$lower))
  expect_true(is.list(pred$uncertainty_prediction$predictions_ci$upper))
  expect_equal(length(pred$uncertainty_prediction$predictions_ci$lower), n_times)
  expect_equal(length(pred$uncertainty_prediction$predictions_ci$upper), n_times)

  for (ii in seq_len(n_times)) {
    expect_true(inherits(pred$uncertainty_prediction$predictions_ci$lower[[ii]], "SpatRaster"))
    expect_true(inherits(pred$uncertainty_prediction$predictions_ci$upper[[ii]], "SpatRaster"))
  }
})

test_that("predict.disag_model_mmap_aghq validates new_data covariate names (gated)", {
  skip_if_aghq_opted_out()
  bundle <- suppressWarnings(get_cached_aghq_fit("aghq_small_onecov_shared"))
  fit <- bundle$fit

  bad_new_data <- lapply(fit$data$covariate_rasters_list, function(r) {
    out <- r
    names(out) <- "wrong_name"
    out
  })

  expect_error(
    predict(fit, new_data = bad_new_data, predict_iid = FALSE, N = 4),
    "Covariate layer names for prediction do not match training"
  )
})

test_that("AGHQ shared-betas fit with two covariates maps slope names and order (gated regression)", {
  skip_if_aghq_opted_out()
  data_obj <- get_cached_aghq_prepared_data("aghq_small_twocov_mesh")

  fit <- suppressWarnings(
    disag_model_mmap(
      data = data_obj,
      engine = "AGHQ",
      family = "poisson",
      link = "log",
      aghq_k = 1,
      field = TRUE,
      iid = TRUE,
      time_varying_betas = FALSE,
      silent = TRUE,
      optimizer = "BFGS"
    )
  )

  expect_s3_class(fit, "disag_model_mmap_aghq")
  expect_false(isTRUE(fit$model_setup$time_varying_betas))

  theta_order <- fit$model_setup$theta_order
  beta_map <- fit$model_setup$beta_index_map
  cov_names <- fit$model_setup$coef_meta$cov_names

  expect_equal(length(cov_names), 2L)
  expect_false(any(is.na(beta_map$slope_idx)))
  expect_equal(theta_order[beta_map$slope_idx], cov_names)
})
