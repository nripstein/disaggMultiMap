test_that("disag_model_mmap accepts engine = MCMC and validates engine.args early", {
  data_obj <- get_cached_prepared_data("prep_default_mesh")

  expect_error(
    disag_model_mmap(
      data = data_obj,
      engine = "MCMC",
      engine.args = list(iter = 0)
    ),
    "`iter` must be an integer-like scalar >= 1"
  )

  expect_error(
    disag_model_mmap(
      data = data_obj,
      engine = "MCMC",
      engine.args = list(iter = 10, warmup = 10)
    ),
    "`warmup` must be strictly less than `iter`"
  )
})


test_that("MCMC engine.args allows pass-through keys", {
  spec <- get_engine_specs_mmap()[["MCMC"]]
  resolved <- resolve_engine_args_mmap(
    engine = "MCMC",
    engine_spec = spec,
    engine_args = list(
      iter = 20,
      warmup = 10,
      refresh = 0,
      control = list(adapt_delta = 0.9)
    ),
    legacy_named_args = list(),
    dot_engine_args = list()
  )

  expect_equal(resolved$iter, 20)
  expect_equal(resolved$warmup, 10)
  expect_equal(resolved$refresh, 0)
  expect_true("control" %in% names(resolved))
  expect_equal(resolved$control$adapt_delta, 0.9)
})


test_that("MCMC fit path errors clearly when tmbstan is unavailable", {
  data_obj <- get_cached_prepared_data("prep_default_mesh")
  skip_if(requireNamespace("tmbstan", quietly = TRUE), "tmbstan is installed; missing dependency path not triggered.")

  expect_error(
    disag_model_mmap_mcmc(
      data = data_obj,
      field = FALSE,
      iid = FALSE,
      silent = TRUE
    ),
    "requires suggested package `tmbstan`"
  )
})
