# Function to test Random Forest predictor sets and SIF scaling
evaluate_predictor_sets <- function(
    sif_cont_train, vi_train, frac_train, sif_ret_train,
    vi_validation, frac_validation, sif_ret_validation,
    class_map, sif_cont_validation
) {
  
  # Overview:
  #   Tests different Random Forest predictor combinations and the effect of
  #   class-specific SIF scaling using independent validation data.
  
  # Requires:
  #   - sif_cont_train (SpatRaster): Class-specific training SIF contribution.
  #   - vi_train (SpatRaster): Training vegetation indices.
  #   - frac_train (SpatRaster): Training vegetation class fractions.
  #   - sif_ret_train (SpatRaster): Training retrieved SIF.
  #   - vi_validation (SpatRaster): Validation vegetation indices.
  #   - frac_validation (SpatRaster): Validation vegetation class fractions.
  #   - sif_ret_validation (SpatRaster): Validation retrieved SIF.
  #   - class_map (SpatRaster): Validation vegetation classification.
  #   - sif_cont_validation (SpatRaster): Validation class-specific SIF contribution.
  
  # Effects:
  #   - Trains one Random Forest model for each predictor set.
  #   - Tests predictions with and without class-specific SIF scaling.
  #   - Calculates RMSE, R², and bias for each vegetation class.
  #   - Produces comparative predictor-set plots.
  
  # Returns:
  #   - A list containing validation results, the combined plot, and fitted models.
  
  # Define predictor combinations to evaluate.
  ps <- list(
    c("Fractions", "Class_id", "Indices", "SIFret"),
    c("Fractions", "Class_id", "Indices"),
    c("Fractions", "Class_id", "SIFret"),
    c("Fractions", "Class_id")
  )
  model_names <- c(
    "All predictors", "Without retrieved SIF", "Without indices",
    "Without indices & retrieved SIF", "All predictors (no scaling)"
  )
  
  # Calculate validation metrics for one vegetation class.
  eval_class <- function(pred, ref, i, id, name, set, removed, scaling) {
    p <- values(pred[[i]], mat = FALSE); r <- values(ref[[i]], mat = FALSE)
    ok <- is.finite(p) & is.finite(r) & p > 0 & r > 0; p <- p[ok]; r <- r[ok]
    cls <- c("Crops", "Mixed Vegetation", "Trees")[i]
    
    if (length(p) < 2) return(data.frame(
      Model_ID = id, Model = name, Predictor_Set = paste(set, collapse = " & "),
      Removed = removed, Scaling = scaling, Class = cls, Validation_N = length(p),
      RMSE = NA_real_, R2 = NA_real_, Bias = NA_real_
    ))
    
    # Calculate prediction errors and validation metrics.
    e <- p - r; ss <- sum((r - mean(r))^2)
    data.frame(
      Model_ID = id, Model = name, Predictor_Set = paste(set, collapse = " & "),
      Removed = removed, Scaling = scaling, Class = cls, Validation_N = length(p),
      RMSE = sqrt(mean(e^2)), R2 = if (ss > 0) 1 - sum(e^2) / ss else NA_real_,
      Bias = mean(e)
    )
  }
  
  # Train one Random Forest model for each predictor combination.
  rf <- lapply(ps, \(set) train_rf_spatial_cv(
    sif_cont_train = sif_cont_train, vi_raster = vi_train, frac = frac_train,
    sif_ret = sif_ret_train, predictor_set = set
  ))
  
  # Predict the independent validation data and calculate class-wise metrics.
  res <- c(
    lapply(seq_along(ps), \(i) {
      set <- ps[[i]]
      rem <- setdiff(c("Fractions", "Class_id", "Indices", "SIFret"), set)
      rem <- if (length(rem)) paste(rem, collapse = " & ") else "None"
      
      # Apply class-specific SIF scaling during validation.
      pred <- predict_sif_cont(
        vi_raster = vi_validation, frac = frac_validation,
        sif_ret = sif_ret_validation, model_fit = rf[[i]],
        class_map = class_map, predictor_set = set, use_scaling = TRUE
      )
      
      do.call(rbind, lapply(1:3, \(j)
                            eval_class(
                              pred, sif_cont_validation, j, i, model_names[i],
                              set, rem, TRUE
                            )
      ))
    }),
    
    # Evaluate the full predictor model without SIF scaling.
    list({
      set <- ps[[1]]
      
      pred <- predict_sif_cont(
        vi_raster = vi_validation, frac = frac_validation,
        sif_ret = sif_ret_validation, model_fit = rf[[1]],
        class_map = class_map, predictor_set = set, use_scaling = FALSE
      )
      
      do.call(rbind, lapply(1:3, \(j)
                            eval_class(
                              pred, sif_cont_validation, j, 5, model_names[5],
                              set, "None", FALSE
                            )
      ))
    })
  )
  res <- do.call(rbind, res)
  
  # Rank models using mean RMSE and absolute bias across vegetation classes.
  ord <- res |>
    dplyr::summarise(
      .by = Model, rm = mean(RMSE, na.rm = TRUE),
      ab = mean(abs(Bias), na.rm = TRUE)
    ) |>
    dplyr::mutate(score = rank(rm) + rank(ab)) |>
    dplyr::arrange(score) |>
    dplyr::pull(Model)
  
  # Set factor order for consistent plotting.
  res$Model <- factor(res$Model, levels = ord)
  res$Class <- factor(res$Class, levels = c("Crops", "Mixed Vegetation", "Trees"))
  
  # Define colours for the predictor-set models.
  cols <- c(
    "All predictors" = "#222222",
    "Without retrieved SIF" = "#0072B2",
    "Without indices" = "#009E73",
    "Without indices & retrieved SIF" = "#D55E00",
    "All predictors (no scaling)" = "#CC79A7"
  )
  
  # Define the common plot theme.
  th <- theme_minimal(base_size = 15) + theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey85", linewidth = .35),
    axis.line = element_line(colour = "black", linewidth = .5),
    axis.ticks = element_line(colour = "black", linewidth = .5),
    axis.title = element_text(colour = "black", size = 15),
    axis.text = element_text(colour = "black", size = 13),
    axis.text.x = element_text(colour = "black", size = 13),
    axis.text.y = element_text(colour = "black", size = 12),
    legend.position = "top", legend.direction = "horizontal", legend.box = "horizontal",
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_blank(), plot.margin = margin(8, 12, 8, 12)
  )
  
  pd <- position_dodge(.86)
  
  # Panel 1: RMSE by vegetation class and predictor set.
  p1 <- ggplot(res, aes(Class, RMSE, fill = Model)) +
    geom_col(position = pd, width = .76) +
    geom_text(
      aes(label = sprintf("%.2f", RMSE)), position = pd,
      vjust = -.35, size = 3.8, colour = "black"
    ) +
    scale_fill_manual(values = cols, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, .12))) +
    labs(x = NULL, y = expression(atop("RMSE", "(mW m"^{-2}*" nm"^{-1}*" sr"^{-1}*")"))) + th
  
  # Panel 2: R² by vegetation class and predictor set.
  p2 <- ggplot(res, aes(Class, R2, fill = Model)) +
    geom_col(position = pd, width = .76) +
    geom_text(
      aes(label = sprintf("%.2f", R2)), position = pd,
      vjust = -.35, size = 3.8, colour = "black"
    ) +
    scale_fill_manual(values = cols, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, 1.08), breaks = seq(0, 1, .2),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = NULL, y = expression(R^2)) + th
  
  # Panel 3: Bias by vegetation class and predictor set.
  bl <- max(abs(res$Bias), na.rm = TRUE) * 1.25
  
  p3 <- ggplot(res, aes(Class, Bias, fill = Model)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = .7, colour = "grey35") +
    geom_col(position = pd, width = .76) +
    geom_text(
      aes(label = sprintf("%+.2f", Bias), vjust = ifelse(Bias >= 0, -.35, 1.2)),
      position = pd, size = 3.8, colour = "black"
    ) +
    scale_fill_manual(values = cols, drop = FALSE) +
    scale_y_continuous(limits = c(-bl, bl), expand = expansion(mult = .07)) +
    labs(x = NULL, y = expression(atop("Bias", "(mW m"^{-2}*" nm"^{-1}*" sr"^{-1}*")"))) + th
  
  # Combine the three metric panels and collect the shared legend.
  plt <- (p1 / p2 / p3) +
    patchwork::plot_layout(guides = "collect", heights = c(1, 1, 1)) &
    theme(
      legend.position = "top", legend.direction = "horizontal",
      legend.box = "horizontal", legend.text = element_text(size = 12, colour = "black")
    ) &
    guides(fill = guide_legend(
      nrow = 1, byrow = TRUE,
      keywidth = grid::unit(.9, "cm"), keyheight = grid::unit(.55, "cm")
    ))
  
  print(plt)
  invisible(list(results = res, plot = plt, models = rf))
}
