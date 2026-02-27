test_that("disag_model_mmap wrapper dispatches to TMB", {
  data_obj <- get_cached_prepared_data("prep_default_mesh")

  fit <- suppressMessages(
    disag_model_mmap(
      data = data_obj,
      engine = "TMB",
      family = "poisson",
      link = "log",
      engine.args = list(iterations = 40),
      field = FALSE,
      iid = FALSE,
      silent = TRUE
    )
  )

  expect_s3_class(fit, "disag_model_mmap_tmb")
  expect_s3_class(fit, "disag_model_mmap")
})

test_that("disag_model_mmap_tmb returns expected object contract", {
  bundle <- get_core_field_fit()
  fit <- bundle$fit

  expect_s3_class(fit, "disag_model_mmap_tmb")
  expect_true(all(c("obj", "opt", "sd_out", "data", "model_setup") %in% names(fit)))
  expect_true(is.list(fit$model_setup))
  expect_equal(fit$model_setup$family, "poisson")
  expect_equal(fit$model_setup$link, "log")
  expect_true(isTRUE(fit$model_setup$field))
  expect_false(isTRUE(fit$model_setup$iid))
  expect_false(isTRUE(fit$model_setup$time_varying_betas))
  expect_true(is.numeric(fit$opt$par))
})

test_that("disag_model_mmap_tmb keeps model flags for no-field/no-iid setup", {
  bundle <- suppressWarnings(get_core_nofield_fit())
  fit <- bundle$fit

  expect_false(isTRUE(fit$model_setup$field))
  expect_false(isTRUE(fit$model_setup$iid))
  expect_equal(fit$model_setup$family, "poisson")
  expect_equal(fit$model_setup$link, "log")
})

test_that("disag_model_mmap_tmb validates family and link inputs", {
  data_obj <- get_cached_prepared_data("prep_default_mesh")

  expect_error(
    disag_model_mmap_tmb(
      data = data_obj,
      family = "banana",
      link = "log",
      iterations = 20,
      field = FALSE,
      iid = FALSE,
      silent = TRUE
    ),
    "family"
  )

  expect_error(
    disag_model_mmap_tmb(
      data = data_obj,
      family = "poisson",
      link = "apple",
      iterations = 20,
      field = FALSE,
      iid = FALSE,
      silent = TRUE
    ),
    "link"
  )
})
