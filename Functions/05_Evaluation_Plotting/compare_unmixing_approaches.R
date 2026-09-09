# Function to compare mixed-pixel SIF prediction performance
compare_unmixing_approaches <- function(
    sif_cont_validation, frac, sif_ret_validation, rf_pred, lm_pred,
    thresholds = c("10–90" = .10, "20–80" = .20, "30–70" = .30, "40–60" = .40)
) {
  
  # Overview:
  #   Compares Random Forest, fraction-baseline, and linear-regression
  #   predictions across different mixed-pixel composition ranges.
  
  # Requires:
  #   - sif_cont_validation (SpatRaster): Class-specific reference SIF contribution.
  #   - frac (SpatRaster): Vegetation class fractions.
  #   - sif_ret_validation (SpatRaster): Retrieved SIF used for the fraction baseline.
  #   - rf_pred (SpatRaster): Random Forest SIF contribution predictions.
  #   - lm_pred (SpatRaster): Linear regression SIF contribution predictions.
  #   - thresholds (numeric): Lower fraction thresholds defining mixed-pixel ranges.
  
  # Effects:
  #   - Calculates RMSE, R², and bias for each model, class, and composition range.
  #   - Produces line plots comparing model performance across mixed-pixel ranges.
  
  # Returns:
  #   - No explicit return; prints the combined performance plots.
  
  # Define vegetation classes, models, colours, and plotting symbols.
  classes <- c("Crops", "Mixed Vegetation", "Trees")
  models <- c("Random forest", "Fraction baseline", "Linear regression")
  cols <- c(
    "Random forest" = "#009E73",
    "Fraction baseline" = "#4D4D4D",
    "Linear regression" = "#D55E00"
  )
  shapes <- c("Random forest" = 16, "Fraction baseline" = 15, "Linear regression" = 17)
  
  # Calculate the fraction-weighted baseline SIF contribution.
  baseline <- sif_ret_validation * frac
  
  # Calculate validation metrics for one observed and predicted vector.
  metric <- function(o, p) {
    ok <- is.finite(o) & is.finite(p) & o > 0
    o <- o[ok]; p <- p[ok]
    
    if (length(o) < 2) return(c(RMSE = NA, R2 = NA, Bias = NA))
    
    e <- p - o; sst <- sum((o - mean(o))^2)
    c(
      RMSE = sqrt(mean(e^2)),
      R2 = if (sst > 0) 1 - sum(e^2) / sst else NA,
      Bias = mean(e)
    )
  }
  
  # Calculate metrics for each composition range, class, and model.
  d <- do.call(rbind, lapply(names(thresholds), function(th) {
    lo <- thresholds[[th]]; hi <- 1 - lo
    
    do.call(rbind, lapply(seq_along(classes), function(i) {
      f <- values(frac[[i]])[, 1]; o <- values(sif_cont_validation[[i]])[, 1]
      
      pred <- list(
        "Random forest" = values(rf_pred[[i]])[, 1],
        "Fraction baseline" = values(baseline[[i]])[, 1],
        "Linear regression" = values(lm_pred[[i]])[, 1]
      )
      
      ok <- is.finite(f) & f >= lo & f <= hi & is.finite(o) & o > 0
      
      do.call(rbind, lapply(names(pred), function(m) {
        p <- pred[[m]]; z <- ok & is.finite(p)
        data.frame(
          Threshold = th, Class = classes[i], Model = m,
          t(metric(o[z], p[z]))
        )
      }))
    }))
  }))
  
  # Set factor levels to control plotting order.
  d$Threshold <- factor(d$Threshold, levels = names(thresholds))
  d$Class <- factor(d$Class, levels = classes)
  d$Model <- factor(d$Model, levels = models)
  
  # Define the common plot theme.
  th <- theme_minimal(base_size = 14) + theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = .35),
    axis.line = element_line(color = "black", linewidth = .5),
    axis.ticks = element_line(color = "black", linewidth = .4),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    legend.position = "top", legend.text = element_text(size = 12, color = "black"),
    legend.key.width = unit(1.2, "cm"), legend.key.height = unit(.45, "cm"),
    plot.margin = margin(2, 8, 2, 5)
  )
  
  # Create a reusable line plot for each validation metric.
  lineplot <- function(y, ylab, ylim, strip = TRUE, xlab = NULL) {
    ggplot(d, aes(Threshold, .data[[y]], group = Model, color = Model, shape = Model)) +
      geom_line(linewidth = .8) + geom_point(size = 3.2) +
      facet_wrap(~Class, nrow = 1) +
      scale_color_manual(values = cols, limits = models) +
      scale_shape_manual(values = shapes, limits = models) +
      scale_y_continuous(limits = ylim, expand = expansion(mult = c(.02, .08))) +
      labs(x = xlab, y = ylab, color = NULL, shape = NULL) + th +
      theme(strip.text = if (strip) element_text(size = 14, face = "bold") else element_blank())
  }
  
  # Panel 1: R² across mixed-pixel composition ranges.
  p1 <- lineplot("R2", expression(R^2), c(0, 1.08), TRUE)
  
  # Panel 2: RMSE across mixed-pixel composition ranges.
  p2 <- lineplot(
    "RMSE",
    expression(atop("RMSE", "(mW m"^{-2}*" nm"^{-1}*" sr"^{-1}*")")),
    c(0, max(d$RMSE, na.rm = TRUE) * 1.18), FALSE
  )
  
  # Panel 3: Bias across mixed-pixel composition ranges.
  lim <- max(abs(d$Bias), na.rm = TRUE) * 1.25
  p3 <- lineplot(
    "Bias",
    expression(atop("Bias", "(mW m"^{-2}*" nm"^{-1}*" sr"^{-1}*")")),
    c(-lim, lim), FALSE, "Mixed-pixel composition (%)"
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = .6, color = "grey35")
  
  # Combine the three performance plots and collect the shared legend.
  print((p1 / p2 / p3) + plot_layout(guides = "collect") &
          theme(legend.position = "top"))
}
