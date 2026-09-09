# Function to test SIF unmixing performance across purity thresholds
assess_purity_thresholds <- function(
    thresholds = c(.75, .80, .85, .90, .95, 1),
    class_map_l1_30m, class_map_l2_30m, class_frac_l1_30m, class_frac_l2_30m,
    desis_l1_idxs_30m, desis_l2_idxs_30m, desis_l1_sfmnn_30m, desis_l2_sfmnn_30m,
    hyplant_2020_idxs_30m, hyplant_2020_sfm_30m, class_frac_2020_30m,
    class_map_2020_30m, class_sif_mean_2020_sfm_30m, study_area_l1, study_area_l2
) {
  
  # Overview:
  #   Evaluates SIF unmixing performance across vegetation purity thresholds
  #   using independent HyPlant 2020 validation data at 30 m.
  
  # Requires:
  #   - thresholds (numeric): Purity thresholds to evaluate.
  #   - class_map_* (SpatRaster): Vegetation classification rasters.
  #   - class_frac_* (SpatRaster): Vegetation class fraction rasters.
  #   - *_idxs_* (SpatRaster): Vegetation index rasters.
  #   - *_sfm* (SpatRaster): Retrieved SIF rasters.
  #   - class_sif_mean_2020_sfm_30m (SpatRaster): Class-average HyPlant 2020 SIF.
  #   - study_area_l1, study_area_l2 (SpatVector): Flight-line study areas.
  
  # Effects:
  #   - Creates pure-pixel training data for each purity threshold.
  #   - Trains a Random Forest model and predicts HyPlant 2020 SIF.
  #   - Calculates validation metrics and produces sensitivity plots.
  
  # Returns:
  #   - A combined patchwork plot showing sensitivity to the purity threshold.
  
  classes <- c("Crops", "Mixed_Vegetation", "Trees")
  labels <- c("Crops", "Mixed Vegetation", "Trees")
  cols <- c(Crops = "#984EA3", `Mixed Vegetation` = "#a6d854", Trees = "#1C7255")
  
  # Count pixels by vegetation class.
  counts <- function(r) {
    f <- terra::freq(r)
    sapply(classes, \(x) sum(f$count[as.character(f$value) == x], na.rm = TRUE))
  }
  
  # Calculate validation metrics for one vegetation class.
  eval_class <- function(pred, ref, i, th, pure_n, total_n) {
    p <- terra::values(pred[[i]], mat = FALSE); r <- terra::values(ref[[i]], mat = FALSE)
    ok <- is.finite(p) & is.finite(r) & p > 0 & r > 0; p <- p[ok]; r <- r[ok]
    
    if (length(p) < 2)
      return(data.frame(
        Threshold = th, Class = labels[i], Pure_N = pure_n, Total_N = total_n,
        Pure_percent = 100 * pure_n / total_n, Validation_N = length(p),
        RMSE = NA, R2 = NA, Bias = NA
      ))
    
    e <- p - r; ss <- sum((r - mean(r))^2)
    data.frame(
      Threshold = th, Class = labels[i], Pure_N = pure_n, Total_N = total_n,
      Pure_percent = 100 * pure_n / total_n, Validation_N = length(p),
      RMSE = sqrt(mean(e^2)), R2 = if (ss > 0) 1 - sum(e^2) / ss else NA, Bias = mean(e)
    )
  }
  
  # Run the workflow for each purity threshold.
  results <- lapply(thresholds, \(th) {
    p1 <- get_pure_class(class_map_l1_30m, class_frac_l1_30m, th)
    p2 <- get_pure_class(class_map_l2_30m, class_frac_l2_30m, th)
    
    # Prepare 60 m training data from one flight line.
    process <- function(cls, idx, sif, area) {
      cls60 <- aggregate_crop_class(cls, 2, mode_resample, area)
      idx <- mask(idx, cls); sif <- mask(sif, cls); frac <- calc_class_frac(cls, cls60, area)
      list(
        idx = resample_crop(idx, cls60, area, method = "average"),
        sif = resample_crop(sif, cls60, area, method = "average"),
        frac = frac, class_sif = agg_sif_class(cls, sif, cls60, frac, area)
      )
    }
    
    a <- process(p1, desis_l1_idxs_30m, desis_l1_sfmnn_30m, study_area_l1)
    b <- process(p2, desis_l2_idxs_30m, desis_l2_sfmnn_30m, study_area_l2)
    
    # Merge training data and train the Random Forest model.
    class_sif <- merge(a$class_sif, b$class_sif)
    vi_raster <- merge(a$idx, b$idx); frac <- merge(a$frac, b$frac); sif_ret <- merge(a$sif, b$sif)
    
    model_fit <- train_rf_spatial_cv(
      sif_cont_train = class_sif, vi_raster = vi_raster, frac = frac, sif_ret = sif_ret,
      predictor_set = c("Indices", "Fractions", "Class_id", "SIFret"), run_sensitivity = TRUE
    )
    
    # Predict independent HyPlant 2020 SIF.
    pred <- predict_sif_cont(
      vi_raster = hyplant_2020_idxs_30m, frac = class_frac_2020_30m,
      sif_ret = hyplant_2020_sfm_30m, model_fit = model_fit, class_map = class_map_2020_30m
    )
    
    # Validate each vegetation class.
    pure_n <- counts(merge(p1, p2)); total_n <- counts(merge(class_map_l1_30m, class_map_l2_30m))
    do.call(rbind, lapply(1:3, \(i) eval_class(
      pred, class_sif_mean_2020_sfm_30m, i, th, pure_n[i], total_n[i]
    )))
  })
  
  # Prepare plot data.
  d <- do.call(rbind, results) |> mutate(
    Class = factor(Class, labels),
    Threshold = factor(Threshold, thresholds, sprintf("%.2f", thresholds)),
    count = sprintf("%d", Pure_N), pct = sprintf("(%.1f%%)", Pure_percent)
  )
  
  # Plot theme.
  theme0 <- theme_minimal(13) + theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = .4),
    axis.ticks = element_line(color = "black"), axis.text = element_text(color = "black", size = 12),
    legend.position = "top", legend.text = element_text(size = 12)
  )
  
  # Create a generic bar plot for validation metrics.
  barplot <- function(y, ylab, label, ylim = NULL) {
    ggplot(d, aes(Threshold, {{ y }}, fill = Class)) +
      geom_col(position = position_dodge(.75), width = .65) +
      geom_text(aes(label = label), position = position_dodge(.75), vjust = -.25, size = 3.8) +
      scale_fill_manual(values = cols) +
      scale_y_continuous(limits = ylim, expand = expansion(mult = c(0, .03))) +
      labs(x = NULL, y = ylab, fill = NULL) + theme0
  }
  
  # Panel A: pure pixels.
  p1 <- ggplot(d, aes(Threshold, Pure_percent, fill = Class)) +
    geom_col(position = position_dodge(.75), width = .65) +
    geom_text(aes(label = pct), position = position_dodge(.75), vjust = -.25, size = 3.8) +
    geom_text(aes(label = count), position = position_dodge(.75), vjust = -1.65, size = 3.8) +
    scale_y_continuous(
      limits = c(0, max(d$Pure_percent, na.rm = TRUE) * 1.3),
      breaks = seq(0, 100, 20), expand = expansion(mult = c(0, .02))
    ) +
    scale_fill_manual(values = cols) + labs(x = NULL, y = "Pure pixels", fill = NULL) + theme0
  
  # Panel B: R².
  p2 <- barplot(R2, expression(R^2), sprintf("%.2f", d$R2), c(0, 1.08))
  
  # Panel C: RMSE.
  p3 <- barplot(
    RMSE, expression(atop("RMSE", "(mW m"^{-2}*" nm"^{-1}*" sr"^{-1}*")")),
    sprintf("%.2f", d$RMSE), c(0, max(d$RMSE, na.rm = TRUE) * 1.18)
  )
  
  # Panel D: bias.
  lim <- max(abs(d$Bias), na.rm = TRUE) * 1.3
  p4 <- ggplot(d, aes(Threshold, Bias, fill = Class)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = .7) +
    geom_col(position = position_dodge(.75), width = .65) +
    geom_text(
      aes(label = sprintf("%+.2f", Bias), vjust = ifelse(Bias >= 0, -.3, 1.2)),
      position = position_dodge(.75), size = 3.8
    ) +
    scale_y_continuous(limits = c(-lim, lim), expand = expansion(mult = .05)) +
    scale_fill_manual(values = cols) +
    labs(
      x = "Purity threshold",
      y = expression(atop("Bias", "(mW m"^{-2}*" nm"^{-1}*" sr"^{-1}*")")), fill = NULL
    ) +
    theme0
  
  # Combine and display the sensitivity plots.
  (p1 / p2 / p3 / p4) + plot_layout(guides = "collect") +
    plot_annotation(title = "") &
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = .5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = .5, size = 12)
    )
}