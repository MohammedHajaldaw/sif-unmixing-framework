# Function to plot class-wise reconstructed versus retrieved SIF
plot_purity_sif_consistency_classwise <- function(predicted_raster, validation_raster, class_map,
                                                  pure_class_map,
                                                  plot_title = "Reconstructed vs Retrieved DESIS SIF",
                                                  panel_type = c("both", "scatter", "density"),
                                                  pixel_type_to_plot = c("Pure", "Mixed", "Combined"),
                                                  show_panel_titles = TRUE) {
  
  # Overview:
  #   Creates class-wise scatter and/or density plots comparing reconstructed
  #   and retrieved SIF for pure, mixed, or combined pixels.
  
  # Requires:
  #   - predicted_raster (SpatRaster): Predicted class-specific SIF contributions.
  #   - validation_raster (SpatRaster): Retrieved SIF used for validation.
  #   - class_map (SpatRaster): Vegetation classification raster.
  #   - pure_class_map (SpatRaster): Raster identifying pure vegetation pixels.
  #   - plot_title (character): Overall plot title.
  #   - panel_type (character): "both", "scatter", or "density".
  #   - pixel_type_to_plot (character): Pixel types to include.
  #   - show_panel_titles (logical): Whether to display individual panel titles.
  
  # Effects:
  #   - Identifies pure and mixed pixels.
  #   - Calculates validation metrics for each class and pixel type.
  #   - Produces scatter and density plots and combines them into a panel.
  
  # Returns:
  #   - A combined patchwork plot.
  
  panel_type <- match.arg(panel_type)
  pixel_type_to_plot <- match.arg(pixel_type_to_plot, several.ok = TRUE)
  
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
  
  # Remove invalid retrieved SIF and pixels outside the classification.
  validation_raster[validation_raster <= 0 | is.na(class_map)] <- NA
  predicted_total <- if (nlyr(predicted_raster) > 1)
    app(predicted_raster, sum, na.rm = TRUE) else predicted_raster
  
  pred <- values(predicted_total)
  retr <- values(validation_raster)
  cls <- values(class_map)
  pm <- values(pure_vs_mixed)
  valid <- complete.cases(pred, retr, cls, pm)
  
  df <- data.frame(
    Predicted = pred[valid],
    Retrieved = retr[valid],
    Class = factor(cls[valid], labels = class_names),
    Pixel = ifelse(pm[valid] == 1, "Pure", "Mixed")
  )
  
  # Add combined pure and mixed pixels.
  df <- bind_rows(df, df %>% mutate(Pixel = "Combined"))
  df <- df %>% filter(Pixel %in% pixel_type_to_plot)
  
  all_plots <- list()
  
  for (c in class_names) {
    df_class <- df %>% filter(Class == c)
    
    for (p in pixel_type_to_plot) {
      df_sub <- df_class %>% filter(Pixel == p)
      
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
      
      # Determine which panels receive axis labels.
      show_x <- c == "Trees" & p == pixel_type_to_plot[2]
      show_y <- c == "Mixed Vegetation" & p == pixel_type_to_plot[1]
      
      # Scatter plot.
      p_scatter <- ggplot(df_sub, aes(x = Retrieved, y = Predicted)) +
        geom_point(color = class_colors[c], size = 2, alpha = 0.6) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                    color = "gray60", linewidth = 0.8) +
        annotate(
          "text", x = min(df_sub$Retrieved), y = max(df_sub$Predicted),
          label = metrics_label, hjust = 0, vjust = 1, size = 4
        ) +
        annotate(
          "text", x = max(df_sub$Retrieved, na.rm = TRUE),
          y = min(df_sub$Predicted, na.rm = TRUE),
          label = paste0("n = ", n), hjust = 1, vjust = 0, size = 4
        ) +
        labs(
          title = if (show_panel_titles) paste(c, "-", p) else NULL,
          x = if (show_x) expression("Retrieved SIF (mW " * m^{-2} * " " * nm^{-1} * " " * sr^{-1} * ")") else NULL,
          y = if (show_y) expression("Reconstructed SIF (mW " * m^{-2} * " " * nm^{-1} * " " * sr^{-1} * ")") else NULL
        ) +
        theme_minimal(base_size = 14) +
        theme(
          panel.grid.major = element_line(color = "gray80", linewidth = 0.2, linetype = "dotted"),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black")
        )
      
      # Density plot.
      df_long <- df_sub %>%
        select(Predicted, Retrieved) %>%
        pivot_longer(
          cols = c(Predicted, Retrieved),
          names_to = "Type", values_to = "SIF"
        )
      
      p_density <- ggplot(df_long, aes(x = SIF, fill = Type)) +
        geom_density(alpha = 0.5, color = NA) +
        geom_vline(
          data = aggregate(SIF ~ Type, df_long, mean),
          aes(xintercept = SIF, color = Type),
          linetype = "dashed", linewidth = 1
        ) +
        scale_fill_manual(
          values = c("Retrieved" = "#5D9BC6", "Predicted" = "#F29C11")
        ) +
        scale_color_manual(
          values = c("Retrieved" = "#5D9BC6", "Predicted" = "#F29C11")
        ) +
        labs(
          title = if (show_panel_titles) paste(c, "-", p) else NULL,
          x = if (show_x) "SIF" else NULL,
          y = if (show_y) "Density" else NULL
        ) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "top", legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black")
        )
      
      all_plots[[paste(c, p, "scatter", sep = "_")]] <- p_scatter
      all_plots[[paste(c, p, "density", sep = "_")]] <- p_density
    }
  }
  
  # Combine plots.
  n_pixel_types <- length(pixel_type_to_plot)
  
  if (panel_type == "both") {
    final_plot <- wrap_plots(
      lapply(class_names, function(c) {
        lapply(pixel_type_to_plot, function(p) {
          list(
            all_plots[[paste(c, p, "scatter", sep = "_")]],
            all_plots[[paste(c, p, "density", sep = "_")]]
          )
        })
      }) %>% unlist(recursive = FALSE),
      ncol = 2 * n_pixel_types
    ) + plot_annotation(title = plot_title)
    
  } else if (panel_type == "scatter") {
    final_plot <- wrap_plots(
      lapply(class_names, function(c) {
        lapply(pixel_type_to_plot, function(p)
          all_plots[[paste(c, p, "scatter", sep = "_")]])
      }) %>% unlist(recursive = FALSE),
      ncol = n_pixel_types
    ) + plot_annotation(title = plot_title)
    
  } else {
    final_plot <- wrap_plots(
      lapply(class_names, function(c) {
        lapply(pixel_type_to_plot, function(p)
          all_plots[[paste(c, p, "density", sep = "_")]])
      }) %>% unlist(recursive = FALSE),
      ncol = n_pixel_types
    ) + plot_annotation(title = plot_title)
  }
  
  return(final_plot)
}