.disaggregation_ns <- function(name) {
  utils::getFromNamespace(name, "disaggregation")
}

disagg_check_new_data <- function(new_data, model_output) {
  .disaggregation_ns("check_new_data")(new_data, model_output)
}

disagg_get_coords <- function(data) {
  .disaggregation_ns("getCoords")(data)
}

disagg_extract_coords_for_mesh <- function(raster, selectIds = NULL) {
  fn <- .disaggregation_ns("extractCoordsForMesh")
  if (is.null(selectIds)) {
    fn(raster)
  } else {
    fn(raster, selectIds = selectIds)
  }
}

disagg_setup_hess_control <- function(opt, hess_control_parscale, hess_control_ndeps) {
  .disaggregation_ns("setup_hess_control")(opt, hess_control_parscale, hess_control_ndeps)
}
