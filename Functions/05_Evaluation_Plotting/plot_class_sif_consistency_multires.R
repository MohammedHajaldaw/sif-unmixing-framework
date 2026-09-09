# Function to plot class-wise SIF comparisons across three resolutions
plot_class_sif_consistency_multires <- function(
    predicted_30m, validation_30m, class_map_30m,
    predicted_60m, validation_60m, class_map_60m,
    predicted_300m, validation_300m, class_map_300m,
    plot_title = "",
    metrics_pos = list(
      "Crops" = c(0.05, 0.95),
      "Mixed Vegetation" = c(0.05, 0.95),
      "Trees" = c(0.05, 1.15)
    )
) {
  
  # Overview:
  #   Creates a 3 × 3 grid of scatter plots comparing reconstructed and
  #   retrieved SIF across three resolutions and three vegetation classes.
  
  # Requires:
  #   - predicted_*m (SpatRaster): Predicted class-specific SIF contributions.
  #   - validation_*m (SpatRaster): Retrieved SIF used for validation.
  #   - class_map_*m (SpatRaster): Vegetation classification raster.
  #   - plot_title (character): Overall plot title.
  #   - metrics_pos (list): Relative positions of the metric annotations by class.
  
  # Effects:
  #   - Calculates reconstructed SIF and validation metrics for each resolution
  #     and vegetation class.
  #   - Produces a 3 × 3 grid of scatter plots with regression and 1:1 lines.
  
  # Returns:
  #   - A combined patchwork plot.
  
  class_names <- c("Crops", "Mixed Vegetation", "Trees")
  class_colors <- c(
    "Crops" = "#984EA3",
    "Mixed Vegetation" = "#a6d854",
    "Trees" = "#1C7255"
  )
  
  # Prepare data for one resolution.
  process_resolution <- function(predicted_raster, validation_raster, class_map, res_label) {
    predicted_raster[predicted_raster <= 0 | is.na(class_map)] <- NA
    validation_raster[validation_raster <= 0 | is.na(class_map)] <- NA
    
    predicted_total <- if (nlyr(predicted_raster) > 1)
      app(predicted_raster, sum, na.rm = TRUE) else predicted_raster
    
    pred <- values(predicted_total)
    retr <- values(validation_raster)
    cls <- values(class_map)
    valid <- complete.cases(pred, retr, cls)
    
    data.frame(
      Predicted = pred[valid],
      Retrieved = retr[valid],
      Class = factor(cls[valid], labels = class_names),
      Resolution = res_label
    )
  }
  
  # Combine data from all resolutions.
  df_all <- bind_rows(
    process_resolution(predicted_30m, validation_30m, class_map_30m, "30m"),
    process_resolution(predicted_60m, validation_60m, class_map_60m, "60m"),
    process_resolution(predicted_300m, validation_300m, class_map_300m, "300m")
  )
  
  plot_list <- list()
  
  for (res in c("30m", "60m", "300m")) {
    for (cls in class_names) {
      df_sub <- df_all %>% filter(Resolution == res, Class == cls)
      
      # Calculate validation metrics.
      n <- nrow(df_sub)
      error <- df_sub$Predicted - df_sub$Retrieved
      rmse <- sqrt(mean(error^2))
      bias <- mean(error)
      sst <- sum((df_sub$Retrieved - mean(df_sub$Retrieved))^2)
      r2 <- if (sst > 0) 1 - sum(error^2) / sst else NA
      nrmse <- rmse / diff(range(df_sub$Retrieved)) * 100
      
      metrics_label <- paste0(
        "R² = ", round(r2, 2), "\n",
        "RMSE = ", round(rmse, 3), "\n",
        "NRMSE = ", round(nrmse, 2), "%\n",
        "Bias = ", round(bias, 3)
      )
      
      show_x <- res == "60m" && cls == "Trees"
      show_y <- res == "30m" && cls == "Mixed Vegetation"
      pos <- metrics_pos[[cls]]
      
      x_pos <- min(df_sub$Retrieved, na.rm = TRUE) +
        diff(range(df_sub$Retrieved, na.rm = TRUE)) * pos[1]
      
      y_pos <- min(df_sub$Predicted, na.rm = TRUE) +
        diff(range(df_sub$Predicted, na.rm = TRUE)) * pos[2]
      
      if (res == "300m" && cls == "Trees")
        y_pos <- y_pos + diff(range(df_sub$Predicted, na.rm = TRUE)) * 0.4
      
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
          "text",
          x = max(df_sub$Retrieved, na.rm = TRUE),
          y = min(df_sub$Predicted, na.rm = TRUE) -
            ifelse(
              res == "300m" && cls == "Trees",
              diff(range(df_sub$Predicted, na.rm = TRUE)) * 0.4, 0
            ),
          label = paste0("n = ", n), hjust = 1, vjust = 0, size = 4
        ) +
        labs(
          title = NULL,
          x = if (show_x)
            expression("Retrieved SIF (mW " * m^{-2} * " " * nm^{-1} * " " * sr^{-1} * ")")
          else NULL,
          y = if (show_y)
            expression("Reconstructed SIF (mW " * m^{-2} * " " * nm^{-1} * " " * sr^{-1} * ")")
          else NULL
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
  
  # Arrange rows by resolution and columns by vegetation class.
  final_plot <- wrap_plots(
    plot_list[c(
      "30m_Crops", "60m_Crops", "300m_Crops",
      "30m_Mixed Vegetation", "60m_Mixed Vegetation", "300m_Mixed Vegetation",
      "30m_Trees", "60m_Trees", "300m_Trees"
    )],
    ncol = 3
  ) + plot_annotation(title = plot_title)
  
  return(final_plot)
}
