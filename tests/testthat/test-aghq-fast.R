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
