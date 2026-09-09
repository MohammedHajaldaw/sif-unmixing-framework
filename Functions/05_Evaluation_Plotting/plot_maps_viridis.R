# Function to plot raster maps using a Viridis colour scale
plot_maps_viridis <- function(
    raster_data, gamma = 0.4, ncol = 2, titles = NULL, units = NULL,
    axis_font_size = 1.5, title_font_size = 2
) {
  
  # Overview:
  #   Plots raster maps using a Viridis colour scale and
  #   gamma correction applied independently to each raster layer.
  
  # Requires:
  #   - raster_data (SpatRaster): Raster containing multiple layers.
  #   - gamma (numeric): Gamma correction applied to the colour scaling.
  #   - ncol (integer): Number of map columns.
  #   - titles (character, optional): Layer titles.
  #   - units (character, optional): Units displayed in the legend.
  #   - axis_font_size (numeric): Axis text size.
  #   - title_font_size (numeric): Map title size.
  
  # Effects:
  #   - Applies per-layer normalization and gamma correction.
  #   - Produces the requested maps using a Viridis colour scale.
  
  # Returns:
  #   - No explicit return; produces the requested map plots.
  
  if (!inherits(raster_data, "SpatRaster"))
    stop("`raster_data` must be a SpatRaster.")
  
  n_layers <- terra::nlyr(raster_data)
  if (!is.null(titles) && length(titles) == n_layers)
    names(raster_data) <- titles
  
  # Set up the plotting layout.
  nrow <- ceiling(n_layers / ncol)
  oldpar <- par(mfrow = c(nrow, ncol), mar = c(2, 2, 2, 5))
  on.exit(par(oldpar), add = TRUE)
  
  for (i in seq_len(n_layers)) {
    lyr <- raster_data[[i]]
    max_val <- terra::global(lyr, "max", na.rm = TRUE)[[1]]
    
    if (is.na(max_val) || max_val == 0) {
      plot(lyr, col = "#000000", main = names(raster_data)[i])
      next
    }
    
    # Normalize, apply gamma correction, and rescale to the original range.
    scaled <- lyr / max_val
    corrected <- scaled ^ gamma
    plotvals <- corrected * max_val
    
    pal <- viridis::viridis(100)
    
    # Plot using a common Viridis palette for each layer.
    terra::plot(
      plotvals, col = pal, main = names(raster_data)[i], axes = FALSE,
      cex.axis = axis_font_size, cex.main = title_font_size,
      plg = list(
        title = if (!is.null(units)) paste0("(", units, ")") else "",
        cex = 1.5
      ),
      range = c(0, max_val + 1e-6), # Ensure max_val is included in the plot range.
      legend = TRUE, box = TRUE
    )
  }
}