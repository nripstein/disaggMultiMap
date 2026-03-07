[![R-CMD-check](https://github.com/nripstein/disaggMultiMap/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nripstein/disaggMultiMap/actions/workflows/R-CMD-check.yaml)

## Installation

You can install from Github with:

``` r
# install.packages("devtools")
devtools::install_github("nripstein/disaggMultiMap")
```

or the development branch from
``` r
devtools::install_github("nripstein/disaggMultiMap", ref = "development")
```

## Current Workflow

```r
# 1) Prepare multi-time data
dat <- prepare_data_mmap(
  polygon_shapefile_list = polygon_list,
  covariate_rasters_list = covariate_list,
  aggregation_rasters_list = agg_list
)

# 2) Fit model (TMB or AGHQ)
fit <- disag_model_mmap(dat, engine = "AGHQ")

# 3) Predict
pred <- predict(fit)
```

## Passing Engine-Specific Arguments

Use `engine.args` as the preferred interface for engine-specific controls.

```r
# AGHQ controls
fit_aghq <- disag_model_mmap(
  dat,
  engine = "AGHQ",
  engine.args = list(
    aghq_k = 2,
    optimizer = "BFGS"
  )
)

# TMB controls
fit_tmb <- disag_model_mmap(
  dat,
  engine = "TMB",
  engine.args = list(
    iterations = 1000,
    hess_control_ndeps = 1e-4
  )
)
```

