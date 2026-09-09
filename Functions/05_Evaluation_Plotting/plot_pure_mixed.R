# Function to plot pure and mixed vegetation pixels
plot_pure_mixed <- function(class_map, pure_class_map, title = NULL) {
  
  # Overview:
  #   Visualizes pure and mixed pixels as a map, donut chart, and
  #   class-specific stacked bar chart.
  
  # Requires:
  #   - class_map (SpatRaster): Vegetation classification raster.
  #   - pure_class_map (SpatRaster): Raster identifying pure vegetation pixels.
  #   - title (character, optional): Plot title for the charts.
  
  # Effects:
  #   - Identifies pure and mixed pixels.
  #   - Produces a map, donut chart, and stacked bar chart.
  
  # Returns:
  #   - No explicit return; produces all three plots.
  
  # Identify pure and mixed pixels.
  pure_vs_mixed <- crop(
    classify(
      ifel(is.na(pure_class_map), 0, pure_class_map),
      rcl = matrix(c(0, Inf, 1), ncol = 3, byrow = TRUE)
    ),
    class_map, touches = TRUE, mask = TRUE
  )
  
  # Map.
  pure_vs_mixed_bin <- pure_vs_mixed
  values(pure_vs_mixed_bin) <- ifelse(values(pure_vs_mixed_bin) == 0, 0, 1)
  pure_vs_mixed_bin <- as.factor(pure_vs_mixed_bin)
  levels(pure_vs_mixed_bin) <- data.frame(ID = c(1, 0), class = c("Pure", "Mixed"))
  
  print(
    tm_shape(pure_vs_mixed_bin) +
      tm_raster(
        palette = c("#4daf4a", "#e41a1c"), title = "Pixel Type",
        labels = c("Pure", "Mixed"), legend.show = TRUE
      ) +
      tm_layout(
        legend.outside = TRUE, legend.title.size = 1,
        legend.text.size = 0.9, frame = FALSE
      )
  )
  
  # Donut chart.
  df <- as.data.frame(freq(pure_vs_mixed)) %>%
    filter(!is.na(value)) %>%
    mutate(
      class = factor(value, levels = c(1, 0), labels = c("Pure", "Mixed")),
      percent = round(100 * count / sum(count), 1),
      label = paste0(count, " (", percent, "%)")
    )
  
  print(
    ggplot(df, aes(x = 2, y = count, fill = class)) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar(theta = "y", start = 0) +
      geom_text(
        aes(label = label), position = position_stack(vjust = 0.5),
        size = 5, color = "white", fontface = "bold"
      ) +
      scale_fill_manual(values = c("#4daf4a", "#e41a1c")) +
      xlim(0.5, 2.5) +
      theme_void() +
      theme(
        legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      ggtitle(ifelse(is.null(title), "Pure vs Mixed Pixels", title))
  )
  
  # Stacked bar chart.
  all_vals <- values(class_map)
  pure_vals <- values(pure_class_map)
  valid_idx <- !is.na(all_vals) & all_vals %in% c(1, 2, 4)
  all_vals <- all_vals[valid_idx]
  pure_vals <- pure_vals[valid_idx]
  
  class_labels <- c("1" = "Crops", "2" = "Mixed_Vegetation", "4" = "Trees")
  all_classes <- factor(
    all_vals, levels = c(1, 2, 4),
    labels = class_labels[as.character(c(1, 2, 4))]
  )
  
  df <- data.frame(
    class = all_classes,
    type = ifelse(!is.na(pure_vals) & pure_vals > 0, "Pure", "Mixed")
  )
  
  df_long <- df %>%
    group_by(class, type) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(class) %>%
    mutate(
      percent = round(100 * count / sum(count), 1),
      label = paste0(count, " (", percent, "%)")
    ) %>%
    ungroup()
  
  # Reverse factor levels so Crops appears at top.
  df_long$class <- factor(df_long$class, levels = c("Trees", "Mixed_Vegetation", "Crops"))
  
  print(
    ggplot(df_long, aes(x = class, y = count, fill = type)) +
      geom_bar(stat = "identity", width = 0.5, color = "white") +
      geom_text(
        aes(label = label), position = position_stack(vjust = 0.5),
        size = 5, color = "white", fontface = "bold"
      ) +
      scale_fill_manual(values = c("Pure" = "#4daf4a", "Mixed" = "#e41a1c")) +
      coord_flip() +
      labs(
        title = ifelse(is.null(title), "", title),
        x = "", y = "Number of Pixels", fill = "Pixel Type"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top", legend.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 14, color = "black", face = "bold"),
        axis.text.x = element_text(color = "black"),
        axis.title.x = element_text(color = "black", face = "bold", size = 14),
        axis.title.y = element_text(color = "black", face = "bold"),
        panel.grid.major = element_line(color = "gray80", linewidth = 0.3, linetype = "dashed"),
        panel.grid.minor = element_line(color = "gray90", linewidth = 0.2, linetype = "dotted"),
        axis.line = element_line(color = "black")
      )
  )
}