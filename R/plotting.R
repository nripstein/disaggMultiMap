#' Plot polygon response data
#'
#' @description
#' Draws the prepared polygons colored by the response variable, with an optional title.
#'
#' @param disag_data A 'disag_data_mmap' object.
#' @param time Integer index of time-slice to plot (default = 1).
#' @param show_title Logical; if TRUE (default), add a title "Response at time X".
#' @return A ggplot2 object.
#' @export
plot_polygons <- function(disag_data,
                          time = 1,
                          show_title = TRUE) {
  sf_obj <- disag_data$polygon_shapefile_list[[time]]
  field_name <- disag_data$shapefile_names$response_var
  # Build map
  p <- ggplot2::ggplot(sf_obj) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[field_name]])) +
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
#' and coloring by value with a Viridis scale. Automatically detects and handles
#' categorical covariates with appropriate discrete color scales.
#'
#' @param disag_data A 'disag_data_mmap' object.
#' @param covariate Integer index or name of the covariate layer.
#' @param time Integer time-slice (default = 1).
#' @param max_categories Maximum number of unique values to consider categorical (default = 10).
#' @return A ggplot2 object.
#' @export
plot_covariate_raster <- function(disag_data,
                                  covariate = 1,
                                  time = 1,
                                  max_categories = 10) {
  # Validate input object
  if (!inherits(disag_data, "disag_data_mmap")) {
    stop("Input must be a 'disag_data_mmap' object")
  }

  # Validate time parameter
  if (!is.numeric(time) || length(time) != 1 || time < 1) {
    stop("'time' must be a positive integer")
  }

  # 1. Pull out the SpatRaster for this time
  cov_list <- disag_data$covariate_rasters_list

  # Check if covariates exist at all
  if (is.null(cov_list)) {
    stop("No covariate rasters available in this disag_data_mmap object")
  }

  # Check if covariates exist for the requested time
  if (length(cov_list) < time) {
    stop("No covariate rasters available for time = ", time,
         ". Maximum available time is ", length(cov_list))
  }

  # Get the raster for this time
  rast <- cov_list[[time]]
  if (!inherits(rast, "SpatRaster")) {
    stop("Expected a SpatRaster object for time = ", time)
  }

  # Get layer names
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
    if (!is.numeric(covariate) || length(covariate) != 1) {
      stop("'covariate' must be a single integer or character string")
    }
    if (covariate < 1 || covariate > length(lyr_names)) {
      stop("Covariate index must be between 1 and ", length(lyr_names),
           "; got ", covariate)
    }
    lyr <- lyr_names[covariate]
  }

  # Extract the single layer
  single_layer <- rast[[lyr]]

  # Check if all values are NA
  if (all(is.na(terra::values(single_layer)))) {
    warning("Covariate '", lyr, "' contains only NA values")
  }

  # Check if this is a categorical layer
  is_categorical <- is_categorical_layer(
    single_layer,
    lyr,
    disag_data$categorical_covariate_baselines,
    max_categories
  )

  # 3. Convert to data.frame with x,y,value
  df <- terra::as.data.frame(single_layer, xy = TRUE, na.rm = FALSE)
  # rename the layer's column to "value"
  names(df)[names(df) == lyr] <- "value"

  # For categorical data, convert to factor with proper levels
  if (is_categorical) {
    levels <- get_categorical_levels(single_layer)
    df$value <- factor(df$value, levels = levels)
  }

  # 4. Plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::coord_sf(crs = sf::st_crs(disag_data$polygon_shapefile_list[[time]])) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal()

  # Apply appropriate scale based on whether data is categorical
  if (is_categorical) {
    # Use a consistent color palette for categorical variables
    # We use a custom palette that's more distinct than viridis for categories
    cat_colors <- c(
      "#440154", # dark purple
      "#21918c", # teal
      "#fde725", # yellow
      "#35b779", # green
      "#31688e", # blue
      "#443983", # indigo
      "#90d743", # lime
      "#e35932", # orange
      "#d8576b"  # pink
    )

    # If there are more categories than colors, recycle the colors
    n_levels <- length(levels(df$value))
    if (n_levels > length(cat_colors)) {
      cat_colors <- rep_len(cat_colors, n_levels)
    } else {
      cat_colors <- cat_colors[1:n_levels]
    }

    p <- p + ggplot2::scale_fill_manual(
      name = lyr,
      values = cat_colors,
      na.value = "transparent"
    )
  } else {
    p <- p + ggplot2::scale_fill_viridis_c(name = lyr, na.value = "transparent")
  }

  return(p)
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
#' @param max_categories Maximum number of unique values to consider categorical (default = 10).
#' @return A ggdraw object (from cowplot) which can be printed.
#' @keywords internal
plot_prepare_summary <- function(disag_data,
                                 covariate = 1,
                                 time = 1,
                                 max_categories = 10) {
  # Validate input
  if (!inherits(disag_data, "disag_data_mmap")) {
    stop("Input must be a 'disag_data_mmap' object")
  }

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

  # Set panel labels
  panel_labels <- c("Response", "Offset", "INLA Mesh", "Covariates")


  if (has_cov) {
    # Try to build the covariate panel
    tryCatch({
      p3 <- plot_covariate_raster(disag_data, covariate = covariate, time = time,
                                 max_categories = max_categories)

      # For categorical covariates, adjust the legend to be more compact
      if (is.factor(p3$data$value)) {
        p3 <- p3 +
          ggplot2::theme(
            legend.key.size = ggplot2::unit(0.8, "lines"),
            legend.text = ggplot2::element_text(size = 8),
            legend.title = ggplot2::element_text(size = 9),
            legend.margin = ggplot2::margin(0, 0, 0, 0),
            legend.box.margin = ggplot2::margin(-10, 0, 0, 0)
          )
      }

      # 2×2: p1 p2 / p4 p3
      grid <- cowplot::plot_grid(
        p1, p2,
        p4, p3,
        ncol  = 2,
        labels = panel_labels
      )
      # Convert to ggdraw object
      grid <- cowplot::ggdraw(grid)
    }, error = function(e) {
      # If there's an error plotting the covariate, show a message instead
      message("Could not plot covariate: ", e$message)
      empty <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                         label = paste("No covariate available:", e$message),
                         size = 3) +
        ggplot2::theme_void() +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)

      grid <- cowplot::plot_grid(
        p1, p2,
        p4, empty,
        ncol  = 2,
        labels = panel_labels
      )
      # Convert to ggdraw object
      grid <- cowplot::ggdraw(grid)
    })
  } else {
    # Create an informative empty panel
    if (is.null(cov_list)) {
      msg <- "No covariates in this dataset"
    } else if (length(cov_list) < time) {
      msg <- paste("No covariates for time =", time)
    } else {
      msg <- "No covariate layers found"
    }

    empty <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg, size = 3) +
      ggplot2::theme_void() +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)

    grid <- cowplot::plot_grid(
      p1, p2,
      p4, empty,
      ncol  = 2,
      labels = panel_labels
    )
    # Convert to ggdraw object
    grid <- cowplot::ggdraw(grid)
  }

  return(grid)
}


#' Visual summary plot of prepared data
#'
#' @description
#' Combines polygons, aggregation raster, mesh, and (if present) a covariate
#' into a 2×2 grid.
#'
#' @param x A `disag_data_mmap` object.
#' @param y Not used (required for S3 method compatibility).
#' @param covariate Integer or name of the covariate to display (default = 1).
#' @param time Integer time‐slice (default = 1).
#' @param max_categories Maximum number of unique values to consider categorical (default = 10).
#' @param ... Additional arguments passed to plot_prepare_summary.
#' @return A ggdraw object (from cowplot) which can be printed.
#' @export
plot.disag_data_mmap <- function(x,
                                 y = NULL,
                                 ...,
                                 covariate = 1,
                                 time = 1,
                                 max_categories = 10) {
  plot_prepare_summary(x,
    covariate = covariate,
    time = time,
    max_categories = max_categories,
    ...
  )
}
