# Function to plot class-average SIF emission maps
plot_avg_sif_maps <- function(raster_data, class_colors, gamma = 0.4, ncol = 2,
                              titles = NULL, units = NULL, axis_font_size = 1.5,
                              title_font_size = 2) {
  
  # Overview:
  #   Plots class-average SIF emission maps using a common scale and class-specific
  #   black-to-colour gradients with gamma correction.
  
  # Requires:
  #   - raster_data (SpatRaster): Raster containing class-average SIF layers.
  #   - class_colors (character): Colours assigned to each vegetation class.
  #   - gamma (numeric): Gamma correction applied to the colour scaling.
  #   - ncol (integer): Number of map columns.
  #   - titles (character, optional): Layer titles.
  #   - units (character, optional): Units displayed in the legend.
  #   - axis_font_size (numeric): Axis text size.
  #   - title_font_size (numeric): Map title size.
  
  # Effects:
  #   - Creates a common SIF scale across all maps.
  #   - Applies gamma correction and class-specific colour gradients.
  #   - Overlays the valid raster extent as an outline.
  
  # Returns:
  #   - No explicit return; produces the requested map plots.
  
  # Basic checks.
  if (!inherits(raster_data, "SpatRaster"))
    stop("`raster_data` must be a SpatRaster.")
  n_layers <- terra::nlyr(raster_data)
  if (length(class_colors) != n_layers)
    stop("`class_colors` length must match number of layers.")
  if (!is.null(titles) && length(titles) == n_layers)
    names(raster_data) <- titles
  
  # Set the plot layout.
  nrow <- ceiling(n_layers / ncol)
  oldpar <- par(mfrow = c(nrow, ncol), mar = c(2, 2, 2, 5))
  on.exit(par(oldpar), add = TRUE)
  
  # Calculate a common scale across all class-specific maps.
  global_max <- max(terra::global(raster_data, "max", na.rm = TRUE)[[1]], na.rm = TRUE)
  
  # Plot each class-specific SIF layer.
  for (i in seq_len(n_layers)) {
    lyr <- raster_data[[i]]
    study_area <- terra::as.polygons(!is.na(lyr), dissolve = TRUE)
    lyr[lyr == 0] <- NA
    
    max_val <- global_max
    if (is.na(max_val) || max_val == 0) {
      plot(lyr, col = "#000000", main = names(raster_data)[i])
      next
    }
    
    # Normalize, gamma correct, and rescale.
    scaled <- lyr / max_val
    corrected <- scaled ^ gamma
    plotvals <- corrected * max_val
    
    # Create a black-to-class-colour gradient.
    pal <- colorRampPalette(c("#000000", class_colors[i]))(100)
    
    # Plot the class-average SIF map.
    terra::plot(
      plotvals, col = pal, main = names(raster_data)[i], axes = FALSE,
      cex.axis = axis_font_size, cex.main = title_font_size,
      plg = list(title = if (!is.null(units)) paste0("(", units, ")") else "", cex = 1.5),
      range = c(0, max_val + 1e-6), legend = TRUE, box = TRUE
    )
    terra::lines(study_area, col = "grey30", lwd = 0.5)
  }
}