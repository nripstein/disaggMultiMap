test_that("predict.disag_model_mmap_aghq rejects predict_iid = TRUE", {
  fake <- structure(list(), class = c("disag_model_mmap_aghq", "list"))

  expect_error(
    predict(fake, predict_iid = TRUE),
    "`predict_iid = TRUE` is not yet supported\\."
  )
})

test_that("predict.disag_model_mmap_aghq validates N", {
  fake <- structure(list(), class = c("disag_model_mmap_aghq", "list"))

  expect_error(
    predict(fake, N = 0),
    "`N` must be a single positive integer\\."
  )
})

test_that("predict.disag_model_mmap_aghq validates CI", {
  fake <- structure(list(), class = c("disag_model_mmap_aghq", "list"))

  expect_error(
    predict(fake, N = 1, CI = 1),
    "`CI` must be a number strictly between 0 and 1\\."
  )
})

test_that("normalize_fixed_names maps shared slope variants to covariate names in order", {
  coef_meta <- list(p = 2L, n_times = 1L, cov_names = c("temp", "precip"))
  nm <- c("intercept", "slope", "slope1", "log_sigma")

  out <- normalize_fixed_names(nm, coef_meta, time_varying_betas = FALSE)

  expect_equal(out, c("intercept", "temp", "precip", "log_sigma"))
})

test_that("canonicalize_draw_names maps shared slope variants to covariate names in order", {
  coef_meta <- list(p = 2L, n_times = 1L, cov_names = c("temp", "precip"))
  nm <- c("slope", "slope1", "nodemean[1]", "iideffect[2]")

  out <- canonicalize_draw_names(nm, coef_meta, time_varying_betas = FALSE)

  expect_equal(out, c("temp", "precip", "nodemean", "iideffect"))
})

test_that("normalize_fixed_names keeps time-varying slope naming stable", {
  coef_meta <- list(p = 2L, n_times = 2L, cov_names = c("temp", "precip"))
  nm <- c("intercept_t", "intercept_t1", "slope_t", "slope_t1", "slope_t2", "slope_t3")

  out <- normalize_fixed_names(nm, coef_meta, time_varying_betas = TRUE)

  expect_equal(out, c("intercept_t1", "intercept_t2", "temp_t1", "precip_t1", "temp_t2", "precip_t2"))
})
