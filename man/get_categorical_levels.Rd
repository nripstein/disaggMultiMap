% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_helpers.R
\name{get_categorical_levels}
\alias{get_categorical_levels}
\title{Get categorical levels from a raster layer}
\usage{
get_categorical_levels(raster_layer)
}
\arguments{
\item{raster_layer}{A SpatRaster layer}
}
\value{
Character vector of level names
}
\description{
Extracts categorical levels from a raster layer, either from defined levels
or from unique values in the raster.
}
\keyword{internal}
