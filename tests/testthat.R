# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(disaggMultiMap)

# Use devtools::load_all() to load all package functions during development
# This ensures all package functions are available for testing
if (!exists("plot_polygons") && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all()
}

test_check("disaggMultiMap")
