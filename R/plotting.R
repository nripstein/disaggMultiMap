#' Plot polygon response data
#'
#' @description
#' Draws the prepared polygons colored by the response variable, with an optional title.
#'
#' @param disag_data A 'disag_data_mmap' object.
#' @param time Integer index of time-slice to plot (default = 1).
#' @param show_title Logical; if TRUE (default), add a title "Response at time X".
#' @return A ggplot2 object.
#' @examples
#' # default, with title
#' p1 <- plot_polygons(my_disag_data, time = 2)
#' # no title
#' p2 <- plot_polygons(my_disag_data, time = 2, show_title = FALSE)
#' @export
plot_polygons <- function(disag_data,
                          time = 1,
                          show_title = TRUE) {
  sf_obj <- disag_data$polygon_shapefile_list[[time]]
  field_name <- disag_data$shapefile_names$response_var
  # Build map
  p <- ggplot2::ggplot(sf_obj) +
    ggplot2::geom_sf(aes(fill = .data[[field_name]])) +
    ggplot2::scale_fill_viridis_c(
      limits = range(sf_obj[[field_name]], na.rm = TRUE)
    ) +
    ggplot2::labs(fill = field_name) +
    ggplot2::theme_minimal()
  # Add title if requested
  if (show_title) {
    tp <- disag_data$time_points[time]
    p <- p + ggplot2::ggtitle(paste0(field_name, " at time ", tp)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  p
}


#' Plot a single covariate raster
#'
#' @description
#' Renders one layer of the covariate raster stack, preserving the raster's CRS,
#' and coloring by value with a Viridis scale.
#'
#' @param disag_data A 'disag_data_mmap' object.
#' @param covariate Integer index or name of the covariate layer.
#' @param time Integer time-slice (default = 1).
#' @return A ggplot2 object.
#' @export
plot_covariate_raster <- function(disag_data,
                                  covariate = 1,
                                  time = 1) {
  # 1. Pull out the SpatRaster for this time
  cov_list <- disag_data$covariate_rasters_list
  if (is.null(cov_list) || length(cov_list) < time) {
    stop("No covariate rasters available for time = ", time)
  }
  rast <- cov_list[[time]]
  lyr_names <- names(rast)
  if (is.null(lyr_names) || length(lyr_names) == 0) {
    stop("Covariate raster at time = ", time, " has no layers")
  }

  # 2. Resolve which layer
  if (is.character(covariate)) {
    if (!covariate %in% lyr_names) {
      stop(
        "Covariate '", covariate, "' not found; available: ",
        paste(lyr_names, collapse = ", ")
      )
    }
    lyr <- covariate
  } else {
    if (covariate < 1 || covariate > length(lyr_names)) {
      stop("Covariate index must be between 1 and ", length(lyr_names))
    }
    lyr <- lyr_names[covariate]
  }

  single_layer <- rast[[lyr]]

  # 3. Convert to data.frame with x,y,value
  df <- terra::as.data.frame(single_layer, xy = TRUE, na.rm = FALSE)
  # rename the layer’s column to “value”
  names(df)[names(df) == lyr] <- "value"

  # 4. Plot
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(name = lyr) +
    ggplot2::coord_sf(crs = sf::st_crs(disag_data$polygon_shapefile_list[[time]])) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal()
}




#' Plot the offset raster
#'
#' @description
#' Draws the aggregation pixel values used in the fit
#'
#' @param disag_data A 'disag_data_mmap' object.
#' @param time Integer time-slice (default = 1).
#' @return A ggplot2 object.
#' @export
plot_aggregation_raster <- function(disag_data,
                                    time = 1) {
  # 1. Extract coords_for_fit (list or matrix)
  cf <- disag_data$coords_for_fit
  coords <- if (is.list(cf)) {
    if (length(cf) < time) {
      stop(
        "coords_for_fit has length ", length(cf),
        " but you requested time = ", time
      )
    }
    cf[[time]]
  } else if (is.matrix(cf) || inherits(cf, "data.frame")) {
    as.matrix(cf)
  } else {
    stop("`coords_for_fit` must be a list or a matrix/data.frame")
  }

  # 2. Extract aggregation_pixels (list or vector)
  ap <- disag_data$aggregation_pixels
  agg_vals <- if (is.list(ap)) {
    if (length(ap) < time) {
      stop(
        "aggregation_pixels has length ", length(ap),
        " but you requested time = ", time
      )
    }
    ap[[time]]
  } else if (is.numeric(ap)) {
    ap
  } else {
    stop("`aggregation_pixels` must be a list or a numeric vector")
  }

  if (nrow(coords) != length(agg_vals)) {
    stop(
      "Number of coords (", nrow(coords),
      ") does not match number of aggregation values (", length(agg_vals), ")"
    )
  }

  # 3. Build data.frame
  df <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    z = agg_vals
  )

  # 4. Plot with geom_tile() + coord_equal()
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = z)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = "aggregation") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal()
}






#' Plot the SPDE mesh with custom outer/inner boundaries
#'
#' @param disag_data A 'disag_data_mmap' object.
#' @param edge_col Colour for internal mesh edges (default = "grey70").
#' @param edge_size Line width for those edges (default = 0.2).
#' @param outer_col Colour for the outer perimeter (default = "black").
#' @param outer_size Line width for the outer perimeter (default = 1).
#' @param inner_col Colour for any inner perimeter (default = "blue").
#' @param inner_size Line width for inner perimeter (default = 1).
#' @param node_col Colour for mesh nodes (default = "black").
#' @param node_size Size for mesh nodes (default = 0.5).
#' @return A ggplot2 object.
#' @export
plot_mesh <- function(disag_data,
                      edge_col = "grey70",
                      edge_size = 0.2,
                      outer_col = "black",
                      outer_size = 1,
                      inner_col = "blue",
                      inner_size = 1,
                      node_col = "black",
                      node_size = 0.5) {
  mesh <- disag_data$mesh
  if (is.null(mesh)) stop("No mesh found in disag_data$mesh")

  # 1. All triangle edges
  tv <- mesh$graph$tv
  idx <- rbind(
    tv[, 1:2, drop = FALSE],
    tv[, 2:3, drop = FALSE],
    tv[, c(3, 1), drop = FALSE]
  )
  edges_df <- data.frame(
    x    = mesh$loc[idx[, 1], 1],
    y    = mesh$loc[idx[, 1], 2],
    xend = mesh$loc[idx[, 2], 1],
    yend = mesh$loc[idx[, 2], 2]
  )

  # 2. Outer perimeter segments
  bnd_idx <- mesh$segm$bnd$idx
  outer_df <- data.frame(
    x    = mesh$loc[bnd_idx[, 1], 1],
    y    = mesh$loc[bnd_idx[, 1], 2],
    xend = mesh$loc[bnd_idx[, 2], 1],
    yend = mesh$loc[bnd_idx[, 2], 2]
  )

  # 3. Inner perimeter segments (if any)
  if (!is.null(mesh$segm$int$idx) && nrow(mesh$segm$int$idx) > 0) {
    int_idx <- mesh$segm$int$idx
    inner_df <- data.frame(
      x    = mesh$loc[int_idx[, 1], 1],
      y    = mesh$loc[int_idx[, 1], 2],
      xend = mesh$loc[int_idx[, 2], 1],
      yend = mesh$loc[int_idx[, 2], 2]
    )
  } else {
    inner_df <- NULL
  }

  # 4. Nodes
  nodes_df <- data.frame(
    x = mesh$loc[, 1],
    y = mesh$loc[, 2]
  )

  # 5. Build plot
  p <- ggplot2::ggplot() +
    # a) light grey mesh edges
    ggplot2::geom_segment(
      data = edges_df,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = edge_col,
      linewidth = edge_size,
      lineend = "round",
      na.rm = TRUE
    ) +
    # b) outer perimeter
    ggplot2::geom_segment(
      data = outer_df,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = outer_col,
      linewidth = outer_size,
      lineend = "round",
      na.rm = TRUE
    ) +
    # c) inner perimeter (if present)
    {
      if (!is.null(inner_df)) {
        ggplot2::geom_segment(
          data = inner_df,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
          color = inner_col,
          linewidth = inner_size,
          lineend = "round",
          na.rm = TRUE
        )
      } else {
        NULL
      }
    } +
    # d) nodes on top
    ggplot2::geom_point(
      data = nodes_df,
      ggplot2::aes(x = x, y = y),
      color = node_col,
      size = node_size,
      na.rm = TRUE
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_void()

  p
}



#' Visual summary plot of prepared data
#'
#' @description
#' Combines polygons, aggregation raster, mesh, and (if present) a covariate
#' into a 2×2 grid.
#'
#' @param disag_data A `disag_data_mmap` object.
#' @param covariate Integer or name of the covariate to display (default = 1).
#' @param time Integer time‐slice (default = 1).
#' @return A ggdraw object (from cowplot) which can be printed.
#' @method plot disag_data_mmap
#' @export
plot_prepare_summary <- function(disag_data,
                                 covariate = 1,
                                 time = 1) {
  # 1) Prepare the three core panels
  p1 <- plot_polygons(disag_data, time = time, show_title = FALSE)
  p2 <- plot_aggregation_raster(disag_data, time = time)
  p4 <- plot_mesh(disag_data) # mesh bottom‐left

  # 2) Decide if we have a covariate to show
  cov_list <- disag_data$covariate_rasters_list
  has_cov <- !is.null(cov_list) &&
    length(cov_list) >= time &&
    inherits(cov_list[[time]], "SpatRaster") &&
    length(names(cov_list[[time]])) > 0

  if (has_cov) {
    # build the covariate panel
    p3 <- plot_covariate_raster(disag_data, covariate = covariate, time = time)
    # 2×2: p1 p2 / p4 p3
    grid <- cowplot::plot_grid(
      p1, p2,
      p4, p3,
      ncol  = 2,
      labels = "AUTO"
    )
  } else {
    # 3 panels only: p1 p2 / p4 [blank]
    empty <- ggplot2::ggplot() +
      ggplot2::theme_void()
    grid <- cowplot::plot_grid(
      p1, p2,
      p4, empty,
      ncol  = 2,
      labels = c("A", "B", "C", "")
    )
  }

  grid
}


#' Visual summary plot of prepared data
#'
#' @description
#' Combines polygons, aggregation raster, mesh, and (if present) a covariate
#' into a 2×2 grid.
#'
#' @param disag_data A `disag_data_mmap` object.
#' @param covariate Integer or name of the covariate to display (default = 1).
#' @param time Integer time‐slice (default = 1).
#' @return A ggdraw object (from cowplot) which can be printed.
#' @export
plot.disag_data_mmap <- function(disag_data,
                                 covariate = 1,
                                 time = 1) {
  plot_prepare_summary(disag_data,
    covariate = 1,
    time = 1
  )
}
