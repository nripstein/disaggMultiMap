#' Roxygen commands
#'
#' @useDynLib disaggMultiMap
#'
dummy <- function() {
  return(NULL)
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".data", "x", "y", "z", "xend", "yend", "value"))
}
