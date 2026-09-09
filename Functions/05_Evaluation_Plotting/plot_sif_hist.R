# Function to plot a SIF histogram
plot_sif_hist <- function(sif_raster, x_label = expression("SIF (mW" ~ m^{-2} ~ nm^{-1} ~ sr^{-1} ~ ")"),
                          plot_title = NULL, bins = 50, x_interval = 1, fill_palette = viridis::viridis(100, option = "D"),
                          breaks = NULL
) {
  
  # Overview:
  #   Generates a histogram of retrieved SIF values with an optional
  #   discrete colour scheme and a vertical line showing mean SIF.
  
  # Requires:
  #   - sif_raster (SpatRaster): Raster containing SIF values.
  #   - x_label (expression/character): Label for the x-axis.
  #   - plot_title (character, optional): Plot title.
  #   - bins (integer): Number of histogram bins.
  #   - x_interval (numeric): Interval for x-axis ticks when breaks are not provided.
  #   - fill_palette (vector): Colours used to fill the histogram.
  #   - breaks (numeric vector, optional): Breaks for discrete histogram colouring.
  
  # Effects:
  #   - Extracts finite SIF values from the raster.
  #   - Calculates the mean SIF value.
  #   - Generates the histogram and mean SIF line.
  #   - Applies the specified colour scheme and plot formatting.
  
  # Returns:
  #   - A ggplot object containing the SIF histogram.
  
  # Extract finite SIF values.
  sif_values <- values(sif_raster)
  sif_values <- sif_values[is.finite(sif_values)]
  data <- data.frame(SIF = sif_values)
  
  # Calculate the mean SIF value.
  mean_sif <- mean(data$SIF, na.rm = TRUE)
  mean_line <- data.frame(mean_sif = mean_sif, label = "Mean")
  
  # Determine fill mapping and x-axis breaks.
  if (!is.null(breaks)) {
    data$bin <- cut(data$SIF, breaks = breaks, include.lowest = TRUE)
    fill_mapping <- aes(fill = bin)
    fill_scale <- scale_fill_manual(values = fill_palette, guide = "none")
    x_breaks <- breaks
  } else {
    fill_mapping <- aes(fill = ..x..)
    fill_scale <- scale_fill_gradientn(colors = fill_palette, guide = "none")
    x_breaks <- seq(floor(min(data$SIF)),ceiling(max(data$SIF)),by = x_interval)
  }
  
  # Create the histogram.
  p_sif_histogram <- ggplot(data, aes(x = SIF)) +
    geom_histogram(fill_mapping, bins = bins, color = "black", linewidth = 0.1) +
    fill_scale +
    geom_vline(
      aes(xintercept = mean_sif, color = label),
      data = mean_line, linetype = "dashed", size = 1, show.legend = TRUE
    ) +
    scale_color_manual(values = c("Mean" = "red")) +
    scale_x_continuous(
      breaks = sort(unique(c(x_breaks, 0))),
      labels = function(x) ifelse(x == 0, "0", as.character(x))
    ) +
    labs(title = plot_title, x = x_label, y = "Frequency", color = NULL
    ) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 44, face = "bold"),
      axis.title = element_text(size = 38, face = "bold"),
      axis.text = element_text(size = 36),
      axis.line = element_line(color = "black", size = 1),
      panel.grid.major = element_line(color = "grey85", size = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 36)
    )
  
  return(p_sif_histogram)
}