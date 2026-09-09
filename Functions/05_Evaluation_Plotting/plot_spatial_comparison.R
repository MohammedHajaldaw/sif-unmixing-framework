# Function to plot two rasters and their difference
plot_spatial_comparison <- function(raster1, raster2, title1, title2, title_diff,
                                    title_size = 1.5, interval = 0.2) {
  
  # Overview:
  #   Displays two SpatRaster objects side by side followed by their difference
  #   map (raster1 - raster2). The first two rasters use a common Viridis scale,
  #   while the difference map uses a discrete diverging colour scale.
  
  # Requires:
  #   - raster1, raster2 (SpatRaster): Rasters to compare.
  #   - title1, title2 (character): Titles for the input raster plots.
  #   - title_diff (character): Title for the difference map.
  #   - title_size (numeric): Font size multiplier for plot titles.
  #   - interval (numeric): Interval width for the difference map.
  
  # Effects:
  #   - Computes a common range for the two input rasters.
  #   - Calculates the difference raster and symmetric breaks around zero.
  #   - Produces the three plots in a single row.
  
  # Returns:
  #   - No explicit return; produces the requested plots.
  
  if (!inherits(raster1, "SpatRaster") || !inherits(raster2, "SpatRaster"))
    stop("Both inputs must be SpatRaster objects.")
  
  # Compute the common range for the first two rasters.
  min1 <- terra::global(raster1, "min", na.rm = TRUE)$min
  max1 <- terra::global(raster1, "max", na.rm = TRUE)$max
  min2 <- terra::global(raster2, "min", na.rm = TRUE)$min
  max2 <- terra::global(raster2, "max", na.rm = TRUE)$max
  global_min <- min(min1, min2)
  global_max <- max(max1, max2)
  
  # Calculate the difference raster and its range.
  diff_raster <- raster1 - raster2
  diff_min <- round(terra::global(diff_raster, "min", na.rm = TRUE)$min, 2)
  diff_max <- round(terra::global(diff_raster, "max", na.rm = TRUE)$max, 2)
  max_abs <- max(abs(diff_min), abs(diff_max))
  
  # Create symmetric breaks around zero using the specified interval.
  half_int <- interval / 2
  n <- ceiling((max_abs + half_int) / interval)
  breaks_full <- seq(
    -n * interval - half_int, n * interval + half_int, by = interval
  )
  
  breaks <- breaks_full[breaks_full >= diff_min & breaks_full <= diff_max]
  breaks <- round(breaks, 2)
  
  # Ensure the actual data range is included in the breaks.
  if (min(breaks) > diff_min) breaks <- c(diff_min, breaks)
  if (max(breaks) < diff_max) breaks <- c(breaks, diff_max)
  
  diff_palette <- c("#0571b0", "lightgrey", "#f4a582", "#ca0020")
  
  # Generate interval labels automatically.
  labels <- paste0(
    ifelse(head(breaks, -1) >= 0, " ", ""), head(breaks, -1), " to ",
    ifelse(tail(breaks, -1) >= 0, " ", ""), tail(breaks, -1)
  )
  
  # Set up the three-panel layout.
  par(mfrow = c(1, 3))
  
  color_palette <- viridis::viridis(100, option = "D")
  
  # Plot the first raster.
  plot(
    raster1, main = title1, col = color_palette, cex.main = title_size,
    range = c(global_min, global_max), plg = list(cex = 2),
    axes = FALSE, box = TRUE
  )
  
  # Plot the second raster.
  plot(
    raster2, main = title2, col = color_palette, cex.main = title_size,
    range = c(global_min, global_max), plg = list(cex = 2),
    axes = FALSE, box = TRUE
  )
  
  # Plot the difference map with discrete intervals.
  terra::plot(
    diff_raster, main = title_diff, col = diff_palette, breaks = breaks,
    type = "interval", mar = c(3, 3, 3, 3),
    plg = list(cex = 2, x = "topright", legend = labels),
    axes = FALSE, box = FALSE
  )
  
  # Reset the plotting layout.
  par(mfrow = c(1, 1))
}