skip_if_aghq_opted_out_g12 <- function() {
  run_aghq <- !tolower(Sys.getenv("RUN_AGHQ_TESTS", "true")) %in% c("0", "false", "no")
  skip_if_not(run_aghq, message = "Set RUN_AGHQ_TESTS=false to skip heavy AGHQ tests.")
}

test_that("G12: AGHQ outer dimension and quadrature grid are invariant to beta count (gated)", {
  skip_if_aghq_opted_out_g12()

  family <- get_test_family_mmap()
  optimizer <- get_test_aghq_optimizer_mmap(family)

  fit_onecov <- suppressWarnings(
    get_cached_aghq_fit(
      name = "aghq_small_onecov_shared",
      aghq_k = 2,
      family = family,
      optimizer = optimizer,
      field = TRUE,
      iid = TRUE,
      time_varying_betas = FALSE
    )$fit
  )

  data_twocov <- get_cached_aghq_prepared_data("aghq_small_twocov_mesh")
  fit_twocov <- suppressWarnings(
    disag_model_mmap(
      data = data_twocov,
      engine = "AGHQ",
      family = family,
      link = "log",
      engine.args = list(
        aghq_k = 2,
        optimizer = optimizer
      ),
      field = TRUE,
      iid = TRUE,
      time_varying_betas = FALSE,
      silent = TRUE
    )
  )

  fit_tv <- suppressWarnings(
    get_cached_aghq_fit(
      name = "aghq_small_onecov_shared",
      aghq_k = 2,
      family = family,
      optimizer = optimizer,
      field = TRUE,
      iid = TRUE,
      time_varying_betas = TRUE
    )$fit
  )

  d_one <- length(fit_onecov$obj$par)
  d_two <- length(fit_twocov$obj$par)
  d_tv <- length(fit_tv$obj$par)

  expect_equal(d_one, d_two)
  expect_equal(d_one, d_tv)

  n_one <- nrow(fit_onecov$aghq_model$normalized_posterior$nodesandweights)
  n_two <- nrow(fit_twocov$aghq_model$normalized_posterior$nodesandweights)
  n_tv <- nrow(fit_tv$aghq_model$normalized_posterior$nodesandweights)

  expect_equal(n_one, n_two)
  expect_equal(n_one, n_tv)
})
