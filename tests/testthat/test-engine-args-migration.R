capture_warnings_mmap <- function(expr) {
  warnings <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

skip_if_aghq_opted_out_engine_args <- function() {
  run_aghq <- !tolower(Sys.getenv("RUN_AGHQ_TESTS", "true")) %in% c("0", "false", "no")
  skip_if_not(run_aghq, message = "Set RUN_AGHQ_TESTS=false to skip heavy AGHQ tests.")
}

test_that("engine.args routes TMB-specific controls", {
  data_obj <- get_cached_prepared_data("prep_default_mesh")

  fit <- suppressMessages(
    disag_model_mmap(
      data = data_obj,
      engine = "TMB",
      family = "poisson",
      link = "log",
      engine.args = list(
        iterations = 20,
        hess_control_ndeps = 1e-4
      ),
      field = FALSE,
      iid = FALSE,
      silent = TRUE
    )
  )

  expect_s3_class(fit, "disag_model_mmap_tmb")
  expect_s3_class(fit, "disag_model_mmap")
})

test_that("engine.args unknown keys warn and are ignored", {
  data_obj <- get_cached_prepared_data("prep_default_mesh")
  out <- capture_warnings_mmap(
    suppressMessages(
      disag_model_mmap(
        data = data_obj,
        engine = "TMB",
        family = "poisson",
        link = "log",
        engine.args = list(iterations = 20, banana = 123),
        field = FALSE,
        iid = FALSE,
        silent = TRUE
      )
    )
  )

  expect_true(any(grepl("Ignoring unknown `engine.args` key", out$warnings)))
  expect_s3_class(out$value, "disag_model_mmap_tmb")
})

test_that("legacy engine-specific args in dots are deprecated but still supported", {
  data_obj <- get_cached_prepared_data("prep_default_mesh")
  out <- capture_warnings_mmap(
    suppressMessages(
      disag_model_mmap(
        data = data_obj,
        engine = "TMB",
        family = "poisson",
        link = "log",
        iterations = 20,
        field = FALSE,
        iid = FALSE,
        silent = TRUE
      )
    )
  )

  expect_true(any(grepl("Engine-specific arguments in `\\.\\.\\.` are deprecated", out$warnings)))
  expect_s3_class(out$value, "disag_model_mmap_tmb")
})

test_that("AGHQ-specific top-level args are warned and ignored under TMB", {
  data_obj <- get_cached_prepared_data("prep_default_mesh")
  out <- capture_warnings_mmap(
    suppressMessages(
      disag_model_mmap(
        data = data_obj,
        engine = "TMB",
        family = "poisson",
        link = "log",
        engine.args = list(iterations = 20),
        aghq_k = 3,
        optimizer = "BFGS",
        field = FALSE,
        iid = FALSE,
        silent = TRUE
      )
    )
  )

  expect_true(any(grepl("`aghq_k` is AGHQ-specific and was ignored", out$warnings)))
  expect_true(any(grepl("`optimizer` is AGHQ-specific and was ignored", out$warnings)))
  expect_s3_class(out$value, "disag_model_mmap_tmb")
})

test_that("invalid engine.args container and names are rejected", {
  data_obj <- get_cached_prepared_data("prep_default_mesh")

  expect_error(
    disag_model_mmap(
      data = data_obj,
      engine = "TMB",
      engine.args = 1
    ),
    "`engine.args` must be NULL or a named list"
  )

  expect_error(
    disag_model_mmap(
      data = data_obj,
      engine = "TMB",
      engine.args = list(1)
    ),
    "`engine.args` must be a named list"
  )

  expect_error(
    disag_model_mmap(
      data = data_obj,
      engine = "TMB",
      engine.args = structure(list(10, 20), names = c("iterations", "iterations"))
    ),
    "duplicated names"
  )
})

test_that("invalid engine-specific values are rejected before fit", {
  data_obj <- get_cached_prepared_data("prep_default_mesh")

  expect_error(
    disag_model_mmap(
      data = data_obj,
      engine = "TMB",
      engine.args = list(iterations = 0)
    ),
    "`iterations` must be an integer-like scalar >= 1"
  )

  expect_error(
    disag_model_mmap(
      data = data_obj,
      engine = "TMB",
      engine.args = list(hess_control_ndeps = 0)
    ),
    "`hess_control_ndeps` must be a numeric scalar > 0"
  )

  expect_error(
    disag_model_mmap(
      data = data_obj,
      engine = "AGHQ",
      engine.args = list(aghq_k = 0)
    ),
    "`aghq_k` must be an integer-like scalar >= 1"
  )

  expect_error(
    disag_model_mmap(
      data = data_obj,
      engine = "AGHQ",
      engine.args = list(optimizer = 1)
    ),
    "`optimizer` must be a non-empty character scalar"
  )
})

test_that("AGHQ engine.args are routed and default k in wrapper remains 2 (gated)", {
  skip_if_aghq_opted_out_engine_args()
  data_obj <- get_cached_aghq_prepared_data("aghq_small_onecov_mesh")

  fit <- suppressWarnings(
    disag_model_mmap(
      data = data_obj,
      engine = "AGHQ",
      family = "poisson",
      link = "log",
      engine.args = list(aghq_k = 1, optimizer = "BFGS"),
      field = TRUE,
      iid = TRUE,
      silent = TRUE
    )
  )
  expect_s3_class(fit, "disag_model_mmap_aghq")
  expect_equal(as.integer(fit$aghq_model$normalized_posterior$grid$level[[1]]), 1L)

  fit_default <- suppressWarnings(
    disag_model_mmap(
      data = data_obj,
      engine = "AGHQ",
      family = "poisson",
      link = "log",
      field = TRUE,
      iid = TRUE,
      silent = TRUE
    )
  )
  expect_s3_class(fit_default, "disag_model_mmap_aghq")
  expect_equal(as.integer(fit_default$aghq_model$normalized_posterior$grid$level[[1]]), 2L)
})

test_that("legacy AGHQ top-level args are deprecated but still supported (gated)", {
  skip_if_aghq_opted_out_engine_args()
  data_obj <- get_cached_aghq_prepared_data("aghq_small_onecov_mesh")
  out <- capture_warnings_mmap(
    suppressMessages(
      disag_model_mmap(
        data = data_obj,
        engine = "AGHQ",
        family = "poisson",
        link = "log",
        aghq_k = 1,
        optimizer = "BFGS",
        field = TRUE,
        iid = TRUE,
        silent = TRUE
      )
    )
  )

  expect_true(any(grepl("`aghq_k` in `disag_model_mmap\\(\\)` is deprecated", out$warnings)))
  expect_true(any(grepl("`optimizer` in `disag_model_mmap\\(\\)` is deprecated", out$warnings)))
  expect_s3_class(out$value, "disag_model_mmap_aghq")
})

test_that("AGHQ helper keeps a legacy fixture path for compatibility checks (gated)", {
  skip_if_aghq_opted_out_engine_args()
  out <- capture_warnings_mmap(
    get_cached_aghq_fit("aghq_small_onecov_shared", use_legacy_args = TRUE)
  )

  expect_true(any(grepl("deprecated", out$warnings)))
  expect_s3_class(out$value$fit, "disag_model_mmap_aghq")
})

test_that("engine.args wins in AGHQ conflicts with legacy top-level args (gated)", {
  skip_if_aghq_opted_out_engine_args()
  data_obj <- get_cached_aghq_prepared_data("aghq_small_onecov_mesh")
  out <- capture_warnings_mmap(
    suppressMessages(
      disag_model_mmap(
        data = data_obj,
        engine = "AGHQ",
        family = "poisson",
        link = "log",
        engine.args = list(aghq_k = 1),
        aghq_k = 2,
        field = TRUE,
        iid = TRUE,
        silent = TRUE
      )
    )
  )

  expect_true(any(grepl("Argument conflict for `aghq_k`", out$warnings)))
  expect_s3_class(out$value, "disag_model_mmap_aghq")
  expect_equal(as.integer(out$value$aghq_model$normalized_posterior$grid$level[[1]]), 1L)
})
