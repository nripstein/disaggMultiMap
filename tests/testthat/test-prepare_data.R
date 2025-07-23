library(testthat)
library(sf)
library(terra)

test_that("validate_prepare_data_inputs accepts valid inputs", {
  # ---- create two simple square polygons ----
  # square 1: (0,0)-(1,1)
  poly1 <- st_sf(
    data.frame(id = 1, response = 5),
    geometry = st_sfc(st_polygon(list(rbind(
      c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)
    )))),
    crs = 4326
  )
  # square 2: (1,1)-(2,2)
  poly2 <- st_sf(
    data.frame(id = 2, response = 3),
    geometry = st_sfc(st_polygon(list(rbind(
      c(1,1), c(2,1), c(2,2), c(1,2), c(1,1)
    )))),
    crs = 4326
  )
  polygon_list <- list(poly1, poly2)

  # ---- matching covariate rasters ----
  r1 <- rast(nrows=2, ncols=2, xmin=0, xmax=2, ymin=0, ymax=2)
  values(r1) <- matrix(1:4, 2, 2)
  r2 <- r1 + 4
  cov_list <- list(r1, r2)

  # ---- matching aggregation rasters ----
  a1 <- rast(r1); values(a1) <- rep(1, ncell(a1))
  a2 <- rast(r2); values(a2) <- rep(2, ncell(a2))
  agg_list <- list(a1, a2)

  # ---- call function and expect no error ----
  expect_silent(
    out <- validate_prepare_data_inputs(
      polygon_shapefile_list   = polygon_list,
      covariate_rasters_list   = cov_list,
      aggregation_rasters_list = agg_list,
      id_var       = "id",
      response_var = "response",
      sample_size_var = NULL,
      make_mesh    = TRUE
    )
  )

  # ---- it should return invisibly TRUE ----
  expect_true(is.logical(out) && length(out) == 1 && out)
})
