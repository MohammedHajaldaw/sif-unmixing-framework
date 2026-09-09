# Function to plot class-wise SIF contribution comparisons across three resolutions
plot_classwise_spatial_validation <- function(
    predicted_30m, validation_30m, predicted_60m, validation_60m,
    predicted_300m, validation_300m, plot_title = "",
    metrics_pos = list(
      "Crops" = c(0, 1.1), "Mixed Vegetation" = c(0, 1.1),
      "Trees" = c(0, 1.1)
    )
) {
  
  # Overview:
  #   Creates a 3 × 3 grid of scatter plots comparing predicted and derived
  #   SIF contribution across three resolutions and three vegetation classes.
  
  # Requires:
  #   - predicted_*m (SpatRaster): Predicted SIF contribution rasters.
  #   - validation_*m (SpatRaster): Derived SIF contribution rasters.
  #   - plot_title (character): Overall plot title.
  #   - metrics_pos (list): Relative positions of metric annotations by class.
  
  # Effects:
  #   - Calculates validation metrics for each resolution and vegetation class.
  #   - Produces scatter plots with regression and 1:1 lines.
  
  # Returns:
  #   - A combined patchwork plot.
  
  class_names <- c("Crops", "Mixed Vegetation", "Trees")
  class_colors <- c(
    "Crops" = "#984EA3", "Mixed Vegetation" = "#a6d854", "Trees" = "#1C7255"
  )
  
  # Prepare data for one resolution.
  process_resolution <- function(predicted_raster, validation_raster, res_label) {
    predicted_raster[predicted_raster <= 0] <- NA
    validation_raster[validation_raster <= 0] <- NA
    
    df_all <- lapply(seq_len(3), function(i) {
      pred <- values(predicted_raster[[i]])
      retr <- values(validation_raster[[i]])
      valid <- !is.na(pred) & !is.na(retr)
      
      data.frame(
        Predicted = pred[valid], Retrieved = retr[valid],
        Class = class_names[i], Resolution = res_label
      )
    })
    
    bind_rows(df_all)
  }
  
  # Combine all resolutions.
  df_all <- bind_rows(
    process_resolution(predicted_30m, validation_30m, "30m"),
    process_resolution(predicted_60m, validation_60m, "60m"),
    process_resolution(predicted_300m, validation_300m, "300m")
  )
  
  plot_list <- list()
  
  for (res in c("30m", "60m", "300m")) {
    for (cls in class_names) {
      df_sub <- df_all %>% filter(Resolution == res, Class == cls)
      
      # Calculate validation metrics.
      n <- nrow(df_sub)
      error <- df_sub$Predicted - df_sub$Retrieved
      rmse <- sqrt(mean(error^2, na.rm = TRUE))
      bias <- mean(error, na.rm = TRUE)
      sst <- sum((df_sub$Retrieved - mean(df_sub$Retrieved, na.rm = TRUE))^2)
      r2 <- if (sst > 0) 1 - sum(error^2, na.rm = TRUE) / sst else NA
      nrmse <- rmse / diff(range(df_sub$Retrieved, na.rm = TRUE)) * 100
      
      metrics_label <- paste0(
        "R² = ", round(r2, 2), "\n",
        "RMSE = ", round(rmse, 3), "\n",
        "NRMSE = ", round(nrmse, 2), "%\n",
        "Bias = ", round(bias, 3)
      )
      
      # Determine metric and axis-label positions.
      pos <- metrics_pos[[cls]]
      x_pos <- min(df_sub$Retrieved, na.rm = TRUE) +
        diff(range(df_sub$Retrieved, na.rm = TRUE)) * pos[1]
      y_pos <- min(df_sub$Predicted, na.rm = TRUE) +
        diff(range(df_sub$Predicted, na.rm = TRUE)) * pos[2]
      
      show_x <- res == "60m" & cls == "Trees"
      show_y <- res == "30m" & cls == "Mixed Vegetation"
      
      p <- ggplot(df_sub, aes(x = Retrieved, y = Predicted)) +
        geom_point(color = class_colors[cls], size = 2, alpha = 0.6) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
        geom_abline(
          intercept = 0, slope = 1, linetype = "dashed",
          color = "gray60", linewidth = 0.8
        ) +
        annotate(
          "text", x = x_pos, y = y_pos, label = metrics_label,
          hjust = 0, vjust = 1, size = 4
        ) +
        annotate(
          "text", x = max(df_sub$Retrieved, na.rm = TRUE),
          y = min(df_sub$Predicted, na.rm = TRUE), label = paste0("n = ", n),
          hjust = 1, vjust = 0, size = 4
        ) +
        labs(
          title = NULL,
          x = if (show_x) expression(
            "Derived SIF Contribution (mW " * m^{-2} * " " * nm^{-1} * " " * sr^{-1} * ")"
          ) else NULL,
          y = if (show_y) expression(
            "Predicted SIF Contribution (mW " * m^{-2} * " " * nm^{-1} * " " * sr^{-1} * ")"
          ) else NULL
        ) +
        theme_minimal(base_size = 14) +
        theme(
          panel.grid.major = element_line(
            color = "gray80", linewidth = 0.2, linetype = "dotted"
          ),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black")
        )
      
      plot_list[[paste(res, cls, sep = "_")]] <- p
    }
  }
  
  # Arrange the 3 × 3 grid.
  final_plot <- wrap_plots(
    plot_list[c(
      "30m_Crops", "60m_Crops", "300m_Crops",
      "30m_Mixed Vegetation", "60m_Mixed Vegetation", "300m_Mixed Vegetation",
      "30m_Trees", "60m_Trees", "300m_Trees"
    )],
    ncol = 3
  ) + patchwork::plot_annotation(title = plot_title)
  
  return(final_plot)
}