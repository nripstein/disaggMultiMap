test_that("slice_params_tmb honors full beta indices when betas are not leading", {
  fake_model <- list(
    model_setup = list(
      family = "poisson",
      field = TRUE,
      iid = FALSE,
      time_varying_betas = FALSE,
      fixed_effect_betas = FALSE,
      coef_meta = list(p = 2L, n_times = 1L, cov_names = c("x1", "x2")),
      beta_index_map = list(
        full_intercept_idx = 3L,
        full_slope_idx = c(4L, 5L)
      )
    ),
    data = list(
      time_points = 1L,
      polygon_data = data.frame(response = numeric(0)),
      mesh = list(loc = matrix(0, nrow = 3L, ncol = 2L))
    ),
    obj = list(
      env = list(
        par = stats::setNames(
          rep(0, 8L),
          c("log_sigma", "log_rho", "intercept", "slope", "slope", "nodemean", "nodemean", "nodemean")
        )
      )
    )
  )

  par_vec <- c(10, 20, 30, 40, 50, 60, 70, 80)
  out <- slice_params_tmb(par_vec, fake_model)

  expect_equal(out$intercept, 30)
  expect_equal(out$slope, c(40, 50))
  expect_equal(out$log_sigma, 10)
  expect_equal(out$log_rho, 20)
  expect_equal(out$nodemean, c(60, 70, 80))
})

test_that("slice_params_tmb matches beta_index_map for field + random-betas TMB fits", {
  bundle <- suppressWarnings(get_cached_tmb_fit(
    name = "slice_field_random_betas",
    seed = 24L,
    iterations = 20,
    family = "poisson",
    link = "log",
    field = TRUE,
    iid = FALSE,
    time_varying_betas = FALSE,
    fixed_effect_betas = FALSE
  ))
  fit <- bundle$fit

  raw <- fit$obj$env$last.par.best
  out <- slice_params_tmb(raw, fit)
  beta_map <- fit$model_setup$beta_index_map

  expect_equal(out$intercept, unname(raw[beta_map$full_intercept_idx]))
  expect_equal(out$slope, unname(raw[as.integer(beta_map$full_slope_idx)]))
  expect_equal(length(out$nodemean), nrow(fit$data$mesh$loc))
})
