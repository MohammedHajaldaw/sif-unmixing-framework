# Function to plot combined SIF consistency results with metrics
plot_purity_sif_consistency <- function(
    predicted_raster, validation_raster, class_map, pure_class_map,
    plot_title = "Reconstructed vs Retrieved SIF",
    x_label = expression("Retrieved SIF (mW" ~ m^{-2} ~ nm^{-1} ~ sr^{-1} ~ ")"),
    y_label = expression("Reconstructed SIF (mW" ~ m^{-2} ~ nm^{-1} ~ sr^{-1} ~ ")"),
    metrics_pos = c(0.1, 0.9)) {
  
  # Overview:
  #   Creates combined scatter and density plots comparing reconstructed and
  #   retrieved SIF for pure, mixed, and combined pixels, including metrics.
  
  # Requires:
  #   - predicted_raster (SpatRaster): Predicted class-specific SIF contributions.
  #   - validation_raster (SpatRaster): Retrieved SIF used for consistency assessment.
  #   - class_map (SpatRaster): Vegetation classification raster.
  #   - pure_class_map (SpatRaster): Raster identifying pure vegetation pixels.
  #   - plot_title (character): Overall plot title.
  #   - x_label, y_label (expression/character): Scatter plot axis labels.
  #   - metrics_pos (numeric): Relative position of the metrics annotation.
  
  # Effects:
  #   - Identifies pure and mixed pixels.
  #   - Calculates bias, RMSE, NRMSE, and R² for each pixel group.
  #   - Produces combined scatter and density plots.
  
  # Returns:
  #   - A combined patchwork plot.
  
  # Define vegetation classes and colours.
  class_names <- c("Crops", "Mixed Vegetation", "Trees")
  class_colors <- c(
    "Crops" = "#984EA3",
    "Mixed Vegetation" = "#a6d854",
    "Trees" = "#1C7255"
  )
  
  # Identify pure and mixed pixels.
  pure_vs_mixed <- crop(
    classify(
      ifel(is.na(pure_class_map), 0, pure_class_map),
      rcl = matrix(c(0, Inf, 1), ncol = 3, byrow = TRUE)
    ),
    class_map, touches = TRUE, mask = TRUE
  )
  
  # Prepare raster values.
  validation_raster[validation_raster <= 0 | is.na(class_map)] <- NA
  predicted_combined <- if (nlyr(predicted_raster) > 1)
    sum(predicted_raster) else predicted_raster
  
  predicted_values <- values(predicted_combined)
  validation_values <- values(validation_raster)
  class_labels <- values(class_map)
  pm_values <- values(pure_vs_mixed)
  valid_idx <- complete.cases(predicted_values, validation_values, class_labels, pm_values)
  
  df <- data.frame(
    Predicted = predicted_values[valid_idx],
    Retrieved = validation_values[valid_idx],
    Class = factor(class_labels[valid_idx], labels = class_names),
    PixelType = pm_values[valid_idx]
  )
  
  # Create pure, mixed, and combined pixel groups.
  df$PixelGroup <- ifelse(df$PixelType == 1, "Pure", "Mixed")
  df_combined <- df
  df_combined$PixelGroup <- "Combined"
  df_plot <- rbind(df, df_combined)
  df_plot$PixelGroup <- factor(
    df_plot$PixelGroup,
    levels = c("Pure", "Mixed", "Combined")
  )
  
  # Calculate validation metrics.
  metrics_df <- df_plot %>%
    group_by(PixelGroup) %>%
    summarise(
      n = n(),
      bias = mean(Predicted - Retrieved),
      rmse = sqrt(mean((Predicted - Retrieved)^2)),
      nrmse = rmse / diff(range(Retrieved)) * 100,
      r2 = 1 - sum((Retrieved - Predicted)^2) /
        sum((Retrieved - mean(Retrieved))^2),
      x = min(Retrieved) + diff(range(Retrieved)) * metrics_pos[1],
      y = max(Predicted) * metrics_pos[2],
      .groups = "drop"
    )
  
  metrics_df$label <- paste0(
    "R² = ", round(metrics_df$r2, 2), "\n",
    "RMSE = ", round(metrics_df$rmse, 3), "\n",
    "NRMSE = ", round(metrics_df$nrmse, 2), "%\n",
    "Bias = ", round(metrics_df$bias, 3)
  )
  
  metrics_df$x <- min(df_plot$Retrieved) +
    diff(range(df_plot$Retrieved)) * metrics_pos[1]
  metrics_df$y <- max(df_plot$Predicted) * metrics_pos[2]
  
  # Scatter plot.
  p_scatter <- ggplot(df_plot, aes(x = Retrieved, y = Predicted, color = Class)) +
    geom_point(size = 2, alpha = 0.6) +
    scale_color_manual(values = class_colors) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
    geom_abline(
      intercept = 0, slope = 1, color = "gray70",
      linetype = "dashed", linewidth = 0.8
    ) +
    geom_text(
      data = metrics_df, aes(x = x, y = y, label = label),
      inherit.aes = FALSE, hjust = 0, size = 5, color = "gray30"
    ) +
    geom_text(
      data = metrics_df,
      aes(
        x = max(df_plot$Retrieved),
        y = min(df_plot$Predicted),
        label = paste0("n = ", n)
      ),
      inherit.aes = FALSE, hjust = 1, vjust = 0, size = 5, color = "gray30"
    ) +
    facet_wrap(~PixelGroup, ncol = 3) +
    labs(title = plot_title, x = x_label, y = y_label) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 21, face = "bold"),
      axis.title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 13),
      legend.position = "top", legend.title = element_blank(),
      legend.text = element_text(size = 13),
      strip.text = element_text(size = 13, face = "bold"),
      panel.grid.major = element_line(
        color = "grey80", linetype = "dotted", linewidth = 0.3
      ),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 1, color = "black")
    )
  
  # Density plot for combined pixels.
  long_df <- melt(
    df_plot %>% filter(PixelGroup == "Combined"),
    id.vars = c("Class", "PixelGroup"),
    measure.vars = c("Retrieved", "Predicted"),
    variable.name = "Type",
    value.name = "SIF"
  )
  
  long_df$Type <- factor(
    long_df$Type,
    levels = c("Retrieved", "Predicted"),
    labels = c("Retrieved", "Reconstructed")
  )
  
  mean_df <- aggregate(SIF ~ Class + Type + PixelGroup, data = long_df, mean)
  
  p_density <- ggplot(long_df, aes(x = SIF, fill = Type)) +
    geom_density(alpha = 0.5, color = NA) +
    facet_wrap(~Class, scales = "free") +
    geom_vline(
      data = mean_df,
      aes(xintercept = SIF, color = Type),
      linetype = "dashed", linewidth = 1
    ) +
    scale_fill_manual(
      values = c("Retrieved" = "#5D9BC6", "Reconstructed" = "#F29C11")
    ) +
    scale_color_manual(
      values = c("Retrieved" = "#5D9BC6", "Reconstructed" = "#F29C11")
    ) +
    labs(
      x = expression("SIF (mW" ~ m^{-2} ~ sr^{-1} ~ nm^{-1} ~ ")"),
      y = "Density (Combined)"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13),
      strip.text = element_text(size = 13, face = "bold"),
      legend.position = "top", legend.title = element_blank(),
      legend.text = element_text(size = 13),
      axis.line = element_line(linewidth = 1, color = "black"),
      panel.grid = element_blank()
    )
  
  # Combine scatter and density plots.
  combined_plot <- p_scatter / p_density +
    plot_layout(heights = c(2, 0.9))
  
  return(combined_plot)
}