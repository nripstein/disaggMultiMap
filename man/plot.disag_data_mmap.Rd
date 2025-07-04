% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot.disag_data_mmap}
\alias{plot.disag_data_mmap}
\title{Visual summary plot of prepared data}
\usage{
\method{plot}{disag_data_mmap}(disag_data, covariate = 1, time = 1)
}
\arguments{
\item{disag_data}{A `disag_data_mmap` object.}

\item{covariate}{Integer or name of the covariate to display (default = 1).}

\item{time}{Integer time‐slice (default = 1).}
}
\value{
A ggdraw object (from cowplot) which can be printed.
}
\description{
Combines polygons, aggregation raster, mesh, and (if present) a covariate
into a 2×2 grid.
}
