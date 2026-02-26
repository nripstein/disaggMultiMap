## Installation

You can install from Github with:

``` r
# install.packages("devtools")
devtools::install_github("https://github.com/nripstein/disaggMultiMap")
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
fit <- disag_model_mmap(dat, engine = "TMB")
# fit <- disag_model_mmap(dat, engine = "AGHQ", aghq_k = 2)

# 3) Predict
pred <- predict(fit)
```
