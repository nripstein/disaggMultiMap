with_use_poisson_tests <- function(value, expr) {
  old <- Sys.getenv("USE_POISSON_TESTS", unset = NA_character_)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("USE_POISSON_TESTS")
    } else {
      Sys.setenv(USE_POISSON_TESTS = old)
    }
  }, add = TRUE)

  if (is.null(value)) {
    Sys.unsetenv("USE_POISSON_TESTS")
  } else {
    Sys.setenv(USE_POISSON_TESTS = value)
  }

  force(expr)
}

test_that("get_test_family_mmap defaults to negative binomial", {
  with_use_poisson_tests(NULL, {
    expect_equal(get_test_family_mmap(), "negbinomial")
  })
})

test_that("get_test_family_mmap switches to poisson for true-like values", {
  with_use_poisson_tests("true", {
    expect_equal(get_test_family_mmap(), "poisson")
  })
  with_use_poisson_tests("1", {
    expect_equal(get_test_family_mmap(), "poisson")
  })
  with_use_poisson_tests("yes", {
    expect_equal(get_test_family_mmap(), "poisson")
  })
})

test_that("get_test_family_mmap stays negative binomial for false-like values", {
  with_use_poisson_tests("false", {
    expect_equal(get_test_family_mmap(), "negbinomial")
  })
  with_use_poisson_tests("0", {
    expect_equal(get_test_family_mmap(), "negbinomial")
  })
  with_use_poisson_tests("no", {
    expect_equal(get_test_family_mmap(), "negbinomial")
  })
})

test_that("get_test_aghq_optimizer_mmap maps family to expected optimizer", {
  expect_equal(get_test_aghq_optimizer_mmap("negbinomial"), "nlminb")
  expect_equal(get_test_aghq_optimizer_mmap("poisson"), "BFGS")
  expect_error(
    get_test_aghq_optimizer_mmap("gaussian"),
    "Unsupported family"
  )
})
