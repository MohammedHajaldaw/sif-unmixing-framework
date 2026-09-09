# Function to plot 2020 class-specific SIF comparisons
evaluate_class_sif_cont <- function(sif_cont_pred, sif_cont_validation, plot_title = "") {
  
  # Overview:
  #   Creates a 1 × 3 scatter-plot panel comparing predicted and validation
  #   SIF contributions for the three vegetation classes.
  
  # Requires:
  #   - sif_cont_pred (SpatRaster): Predicted class-specific SIF contribution.
  #   - sif_cont_validation (SpatRaster): Validation class-specific SIF contribution.
  #   - plot_title (character): Overall plot title.
  
  # Effects:
  #   - Calculates R², RMSE, NRMSE, and bias for each vegetation class.
  #   - Produces class-specific scatter plots with regression and 1:1 lines.
  
  # Returns:
  #   - A combined patchwork plot containing the three class-specific panels.
  
  # Define vegetation classes and their plot colours.
  class_names <- c("Crops", "Mixed Vegetation", "Trees")
  class_colors <- c(
    "Crops" = "#984EA3",
    "Mixed Vegetation" = "#a6d854",
    "Trees" = "#1C7255"
  )
  
  # Extract valid predicted and validation values for each class.
  df_all <- lapply(seq_along(class_names), function(i) {
    pred <- values(sif_cont_pred[[i]])
    retr <- values(sif_cont_validation[[i]])
    valid <- !is.na(pred) & !is.na(retr) & pred > 0 & retr > 0
    
    data.frame(
      Predicted = pred[valid],
      Retrieved = retr[valid],
      Class = class_names[i]
    )
  }) |> bind_rows()
  
  plot_list <- list()
  
  # Create one scatter plot for each vegetation class.
  for (cls in class_names) {
    df_sub <- df_all %>% filter(Class == cls)
    n <- nrow(df_sub)
    
    # Calculate validation metrics.
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
    
    # Position the metric and sample-size annotations within each panel.
    x_range <- range(df_sub$Retrieved, na.rm = TRUE)
    y_range <- range(df_sub$Predicted, na.rm = TRUE)
    x_pos <- x_range[1] + diff(x_range) * 0.05
    y_pos <- y_range[2] - diff(y_range) * 0.05
    
    # Show axis labels only on the middle and left panels.
    show_x <- cls == "Mixed Vegetation"
    show_y <- cls == "Crops"
    
    p <- ggplot(df_sub, aes(x = Retrieved, y = Predicted)) +
      geom_point(color = class_colors[cls], size = 2, alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
      geom_abline(
        intercept = 0, slope = 1, linetype = "dashed",
        color = "gray60", linewidth = 0.8
      ) +
      annotate(
        "text", x = x_pos, y = y_pos, label = metrics_label,
        hjust = 0, vjust = 1, size = 4, color = "black"
      ) +
      annotate(
        "text",
        x = max(df_sub$Retrieved, na.rm = TRUE),
        y = min(df_sub$Predicted, na.rm = TRUE),
        label = paste0("n = ", n), hjust = 1, vjust = 0,
        size = 4, color = "black"
      ) +
      labs(
        title = cls,
        x = if (show_x)
          expression("Derived SIF Contribution (mW " * m^{-2} * " " * nm^{-1} * " " * sr^{-1} * ")")
        else NULL,
        y = if (show_y)
          expression("Predicted SIF Contribution (mW " * m^{-2} * " " * nm^{-1} * " " * sr^{-1} * ")")
        else NULL
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        panel.grid.major = element_line(
          color = "gray80", linewidth = 0.2, linetype = "dotted"
        ),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
      )
    
    plot_list[[cls]] <- p
  }
  
  # Arrange the class-specific plots from left to right.
  final_plot <- wrap_plots(
    plot_list[c("Crops", "Mixed Vegetation", "Trees")], ncol = 3
  ) + plot_annotation(title = plot_title)
  
  return(final_plot)
}
