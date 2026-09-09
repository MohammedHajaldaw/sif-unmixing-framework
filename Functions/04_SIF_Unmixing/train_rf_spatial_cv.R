# Function to train a Random Forest model using spatial cross-validation
train_rf_spatial_cv <- function(sif_cont_train, vi_raster, frac, sif_ret,
                                predictor_set = c("Indices", "Fractions", "Class_id", "SIFret"),
                                run_sensitivity = FALSE) {
  
  # Overview:
  #   This function trains a single Random Forest model to predict SIF contribution.
  #   Training data are combined across vegetation classes, while spatial dependence
  #   is assessed using class-specific spherical semivariograms to define the CV block size.
  
  # Requires:
  #   - sif_cont_train (SpatRaster): Class-specific SIF contribution rasters used as response.
  #   - vi_raster (SpatRaster): Vegetation index raster stack used as predictors.
  #   - frac (SpatRaster): Raster stack containing vegetation class fractions.
  #   - sif_ret (SpatRaster): Retrieved SIF raster used as a predictor.
  #   - predictor_set (character): Predictor groups included in the Random Forest model.
  #   - run_sensitivity (logical): Whether to evaluate candidate spatial block sizes.
  
  # Effects:
  #   - Builds the combined training dataset.
  #   - Estimates class-specific spherical semivariogram ranges.
  #   - Defines the minimum semivariogram-based spatial block size.
  #   - Optionally evaluates candidate block sizes.
  #   - Assigns pixels to ten spatial CV folds.
  #   - Tunes and trains the final Random Forest model.
  
  # Returns:
  #   - A list containing the trained Random Forest model and prediction metadata.
  
  set.seed(500)
  names(sif_ret) <- "SIFret"
  combined_data <- list()
  num_layers <- nlyr(sif_cont_train)
  
  # Build the training dataset separately for each vegetation class.
  for (i in seq_len(num_layers)) {
    class_id <- sif_ret
    class_id[!is.na(class_id)] <- i
    names(class_id) <- "class_id"
    
    raster_list <- list()
    if ("Indices" %in% predictor_set) raster_list <- c(raster_list, list(vi_raster))
    raster_list <- c(raster_list, list(frac[[i]]))
    if ("Class_id" %in% predictor_set) raster_list <- c(raster_list, list(class_id))
    if ("SIFret" %in% predictor_set) raster_list <- c(raster_list, list(sif_ret))
    raster_list <- c(raster_list, list(sif_cont_train[[i]]))
    
    current_predictor <- do.call(c, raster_list)
    data <- as.data.frame(current_predictor, cells = TRUE, na.rm = TRUE)
    names(data)[c(which(names(data) == "cell"), ncol(data))] <- c("pixel_id", "sif_cont")
    
    xy <- terra::xyFromCell(frac[[i]], data$pixel_id)
    data$x <- xy[, 1]; data$y <- xy[, 2]
    
    fraction_column <- grep("frac", names(data), value = TRUE)
    if (length(fraction_column)) names(data)[names(data) %in% fraction_column] <- "frac"
    else stop("Unable to identify fraction columns in layer ", i)
    
    combined_data[[i]] <- data
  }
  
  # Combine all vegetation-class observations into the unified training dataset.
  combined_data <- do.call(rbind, combined_data)
  
  # Store class levels and maximum training SIF contribution for prediction scaling.
  class_levels <- seq_len(num_layers)
  max_train <- sapply(class_levels, function(i)
    max(terra::values(sif_cont_train[[i]]), na.rm = TRUE))
  names(max_train) <- as.character(class_levels)
  
  # Define semivariogram lag width and cutoff from the input spatial resolution.
  pixel_size <- max(terra::res(sif_ret))
  n_lags <- 10; pixels_per_lag <- 2
  lag_width <- pixel_size * pixels_per_lag
  cutoff <- lag_width * n_lags
  
  class_names <- sub(" .*", "", names(sif_cont_train))
  if (is.null(class_names) || any(class_names == ""))
    class_names <- paste0("Class ", seq_len(num_layers))
  
  class_ranges <- numeric(num_layers)
  
  # Estimate a spherical semivariogram range for each vegetation class.
  for (i in seq_len(num_layers)) {
    df <- as.data.frame(sif_cont_train[[i]], xy = TRUE, na.rm = TRUE)
    names(df) <- c("x", "y", "sif_cont")
    
    v <- variogram(sif_cont ~ 1, locations = ~x + y, data = df,
                   cutoff = cutoff, width = lag_width)
    fit_sph <- tryCatch(suppressWarnings(fit.variogram(v, vgm("Sph"))),
                        error = function(e) NULL)
    
    class_ranges[i] <- if (!is.null(fit_sph) && nrow(fit_sph) >= 2 &&
                           is.finite(fit_sph$range[2]) && fit_sph$range[2] > 0)
      fit_sph$range[2] else NA_real_
  }
  
  # Use the largest valid semivariogram range to determine the minimum spatial block size.
  valid_ranges <- class_ranges[is.finite(class_ranges) & class_ranges > 0]
  if (!length(valid_ranges)) stop("No valid spherical semivariogram ranges could be estimated.")
  
  max_variogram_range <- max(valid_ranges)
  minimum_block_size <- ceiling(max_variogram_range / pixel_size) * pixel_size
  candidate_sizes <- seq(minimum_block_size, minimum_block_size + 120, by = 60)
  
  # The minimum semivariogram-based block size is retained for the final spatial CV.
  best_block_size <- minimum_block_size
  tuning_results <- data.frame(block_size = numeric(), rmse = numeric())
  
  cat("----------------------------------------------------\n")
  cat("Starting Automated Block Size Search\n")
  cat("----------------------------------------------------\n")
  cat("Class-specific spherical semivariogram ranges:\n")
  print(data.frame(Class = class_names, Range_m = round(class_ranges)))
  cat("\nMaximum semivariogram range:", round(max_variogram_range), "m\n")
  cat("Minimum block size:", minimum_block_size, "m\n")
  cat("Candidate block sizes:", paste(candidate_sizes, collapse = ", "), "m\n")
  cat("----------------------------------------------------\n")
  
  # Evaluate candidate block sizes when sensitivity analysis is requested.
  if (run_sensitivity) {
    for (bs in candidate_sizes) {
      x_orig <- min(combined_data$x) - pixel_size / 2
      y_orig <- min(combined_data$y) - pixel_size / 2
      bx <- floor((combined_data$x - x_orig) / bs)
      by <- floor((combined_data$y - y_orig) / bs)
      s_block <- paste(bx, by, sep = "_")
      
      b_counts <- aggregate(pixel_id ~ s_block, combined_data, function(x) length(unique(x)))
      names(b_counts)[2] <- "n_px"
      if (nrow(b_counts) < 10) next
      
      b_counts <- b_counts[order(b_counts$n_px, decreasing = TRUE), ]
      f_counts <- rep(0, 10); b_folds <- numeric(nrow(b_counts))
      
      for (k in seq_len(nrow(b_counts))) {
        sf <- which.min(f_counts); b_folds[k] <- sf
        f_counts[sf] <- f_counts[sf] + b_counts$n_px[k]
      }
      
      names(b_folds) <- b_counts$s_block
      temp_fold <- b_folds[s_block]
      
      # Create ten-fold spatial CV indices for the candidate block size.
      temp_tc <- trainControl(method = "cv", number = 10,
                              index = lapply(1:10, function(f) which(temp_fold != f)),
                              indexOut = lapply(1:10, function(f) which(temp_fold == f)),
                              savePredictions = TRUE, verboseIter = FALSE)
      
      temp_data <- combined_data[, !(names(combined_data) %in% c("pixel_id", "x", "y"))]
      
      tryCatch({
        # Fit a temporary RF model and calculate its out-of-fold RMSE.
        t_model <- train(sif_cont ~ ., data = temp_data, method = "ranger",
                         trControl = temp_tc, tuneLength = 2, num.trees = 500, verbose = FALSE)
        
        bt <- t_model$bestTune; op <- t_model$pred
        bp <- op %>% inner_join(bt, by = intersect(names(op), names(bt)))
        t_rmse <- sqrt(mean((bp$pred - bp$obs)^2))
        
        tuning_results <- rbind(tuning_results, data.frame(block_size = bs, rmse = t_rmse))
        cat("Tested block size:", bs, "m | OOF RMSE:",
            round(t_rmse, 3), "mW m⁻² nm⁻¹ sr⁻¹\n")
      }, error = function(e) {})
    }
  } else {
    cat("Sensitivity analysis: OFF\n")
  }
  
  cat("----------------------------------------------------\n")
  cat(">>> Selected Minimum Semivariogram-Based Block Size:", best_block_size, "m <<<\n")
  cat("----------------------------------------------------\n\n")
  
  block_size <- best_block_size
  n_pixels <- length(unique(combined_data$pixel_id))
  x_origin <- min(combined_data$x) - pixel_size / 2
  y_origin <- min(combined_data$y) - pixel_size / 2
  
  # Assign observations to spatial blocks using the selected block size.
  combined_data$block_x <- floor((combined_data$x - x_origin) / block_size)
  combined_data$block_y <- floor((combined_data$y - y_origin) / block_size)
  combined_data$spatial_block <- paste(combined_data$block_x, combined_data$block_y, sep = "_")
  
  # Count pixels per block and assign blocks to ten folds while balancing pixel counts.
  block_pixel_counts <- aggregate(pixel_id ~ spatial_block, combined_data,
                                  function(x) length(unique(x)))
  names(block_pixel_counts)[2] <- "n_pixels"
  block_pixel_counts <- block_pixel_counts[order(block_pixel_counts$n_pixels, decreasing = TRUE), ]
  
  fold_pixel_counts <- rep(0, 10)
  block_folds <- numeric(nrow(block_pixel_counts))
  
  for (i in seq_len(nrow(block_pixel_counts))) {
    selected_fold <- which.min(fold_pixel_counts)
    block_folds[i] <- selected_fold
    fold_pixel_counts[selected_fold] <- fold_pixel_counts[selected_fold] + block_pixel_counts$n_pixels[i]
  }
  
  names(block_folds) <- block_pixel_counts$spatial_block
  combined_data$fold <- block_folds[combined_data$spatial_block]
  
  # Plot the class-specific semivariograms and fitted spherical models.
  par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3, 1))
  class_colors <- c("#984EA3", "#a6d854", "#1C7255")
  
  for (i in seq_len(num_layers)) {
    df <- as.data.frame(sif_cont_train[[i]], xy = TRUE, na.rm = TRUE)
    names(df) <- c("x", "y", "sif_cont")
    
    v <- variogram(sif_cont ~ 1, locations = ~x + y, data = df,
                   cutoff = cutoff, width = lag_width)
    fit_sph <- tryCatch(suppressWarnings(fit.variogram(v, vgm("Sph"))),
                        error = function(e) NULL)
    
    # Plot the empirical semivariogram if the spherical model cannot be fitted.
    if (is.null(fit_sph) || nrow(fit_sph) < 2 ||
        !is.finite(fit_sph$range[2]) || fit_sph$range[2] <= 0) {
      plot(v$dist, v$gamma, pch = 19, cex = 1, col = class_colors[i],
           xlim = c(0, cutoff), ylim = c(0, max(v$gamma, na.rm = TRUE) * 1.1),
           xlab = "Distance (m)", ylab = "Semivariance", main = class_names[i],
           bty = "l", cex.axis = 1.2, cex.lab = 1.25, cex.main = 1.3)
      next
    }
    
    # Add the fitted spherical model and estimated semivariogram range.
    range_m <- fit_sph$range[2]
    model_line <- variogramLine(fit_sph, dist_vector = seq(0, cutoff, length.out = 300))
    y_max <- max(c(v$gamma, model_line$gamma), na.rm = TRUE) * 1.1
    
    plot(v$dist, v$gamma, pch = 19, cex = 1, col = class_colors[i],
         xlim = c(0, cutoff), ylim = c(0, y_max),
         xlab = "Distance (m)", ylab = "Semivariance", main = class_names[i],
         bty = "l", cex.axis = 1.2, cex.lab = 1.25, cex.main = 1.3)
    lines(model_line$dist, model_line$gamma, lwd = 2, col = class_colors[i])
    abline(v = range_m, lty = 2, lwd = 2)
    text(cutoff * 0.50, y_max * 0.55,
         paste0("Range = ", round(range_m, 0), " m"), cex = 0.9)
  }
  
  # Display the block-size sensitivity results.
  plot.new()
  title(main = "Block-Size Sensitivity Analysis", line = 1)
  
  if (run_sensitivity && nrow(tuning_results) > 0) {
    text(0.32, 0.90, "Block size", cex = 1, font = 2)
    text(0.68, 0.90, "Out-Of-Fold RMSE", cex = 1, font = 2)
    text(0.68, 0.80, expression("(mW m"^{-2} * " nm"^{-1} * " sr"^{-1} * ")"), cex = 0.85)
    segments(0.15, 0.70, 0.85, 0.7, lwd = 1)
    
    for (j in seq_len(nrow(tuning_results))) {
      y <- 0.6 - (j - 1) * 0.12
      f <- ifelse(tuning_results$block_size[j] == best_block_size, 2, 1)
      text(0.32, y, paste0(tuning_results$block_size[j], " m"), cex = 1, font = f)
      text(0.68, y, sprintf("%.3f", tuning_results$rmse[j]), cex = 1, font = f)
    }
    
    segments(0.15, 0.25, 0.85, 0.25, lwd = 1)
    text(0.5, 0.09, paste0(
      "Maximum semivariogram range = ", round(max_variogram_range, 0), " m\n",
      "Minimum block >= range = ", minimum_block_size, " m\n",
      "Selected block = ", best_block_size, " m"), cex = 0.9)
  } else {
    text(0.5, 0.5, "Sensitivity analysis OFF", cex = 1)
  }
  
  par(mfrow = c(1, 1))
  
  # Plot the spatial distribution and number of pixels in each CV fold.
  unique_pixels <- combined_data[!duplicated(combined_data$pixel_id), ]
  fold_colors <- brewer.pal(10, "Paired")
  fold_pixel_counts <- aggregate(pixel_id ~ fold, combined_data,
                                 function(x) length(unique(x)))
  fold_pixel_counts <- fold_pixel_counts[order(fold_pixel_counts$fold), ]
  
  plot(unique_pixels$x, unique_pixels$y, col = fold_colors[unique_pixels$fold],
       pch = 15, asp = 1, axes = FALSE, xlab = "", ylab = "",
       main = "Spatial Cross-Validation Folds", bty = "n", font.main = 2)
  
  # Overlay the boundaries of the spatial blocks.
  block_ids <- unique(unique_pixels[, c("block_x", "block_y")])
  for (i in seq_len(nrow(block_ids))) {
    x0 <- x_origin + block_ids$block_x[i] * block_size
    y0 <- y_origin + block_ids$block_y[i] * block_size
    rect(x0, y0, x0 + block_size, y0 + block_size, border = "grey50", lwd = 0.8)
  }
  
  legend("topright", legend = c(
    paste0("Block size: ", block_size, " m"),
    paste0("Spatial blocks: ", length(unique(combined_data$spatial_block))),
    paste0("Valid pixels: ", n_pixels)), bty = "n", cex = 0.8)
  
  # Add the number of valid pixels assigned to each CV fold.
  par(new = TRUE, fig = c(0.08, 0.48, 0.08, 0.38), mar = c(3, 3.5, 1.5, 0.5))
  
  bp <- barplot(fold_pixel_counts[, 2], names.arg = paste0("Fold ", 1:10),
                col = fold_colors, border = NA,
                ylim = c(0, max(fold_pixel_counts[, 2]) * 1.30),
                yaxt = "n", las = 2, cex.axis = 0.7, cex.names = 0.65, cex.lab = 0.75)
  
  text(bp, fold_pixel_counts[, 2], fold_pixel_counts[, 2], pos = 3, cex = 0.65)
  title("Valid pixels per cross-validation fold", line = 0.9, cex.main = 0.8, font.main = 2)
  
  par(new = FALSE); par(fig = c(0, 1, 0, 1))
  
  # Define the final ten-fold spatial cross-validation.
  train_control <- trainControl(method = "cv", number = 10,
                                index = lapply(1:10, function(f) which(combined_data$fold != f)),
                                indexOut = lapply(1:10, function(f) which(combined_data$fold == f)),
                                savePredictions = TRUE, verboseIter = FALSE)
  
  # Remove spatial identifiers so they cannot be used as model predictors.
  combined_data[c("pixel_id", "x", "y", "block_x", "block_y", "spatial_block", "fold")] <- NULL
  
  # Train and tune the final Random Forest model.
  rf_model <- train(sif_cont ~ ., data = combined_data, method = "ranger",
                    trControl = train_control, tuneLength = 6, importance = "permutation")
  
  # Show the selected tuning parameters.
  best_tune <- rf_model$bestTune
  
  cat("\n====================================\n")
  cat("SELECTED TUNING PARAMETERS:\n")
  cat("====================================\n")
  print(best_tune)
  
  # Calculate out-of-fold RMSE and R-squared for the selected tuning parameters.
  oof_preds <- rf_model$pred
  best_preds <- oof_preds %>%
    inner_join(best_tune, by = intersect(names(oof_preds), names(best_tune)))
  
  spatial_rmse <- sqrt(mean((best_preds$pred - best_preds$obs)^2))
  ss_total <- sum((best_preds$obs - mean(best_preds$obs))^2)
  ss_res <- sum((best_preds$obs - best_preds$pred)^2)
  spatial_r2 <- 1 - ss_res / ss_total
  
  cat("\nSpatial CV Model Evaluation (Out-of-Fold) with Minimum Semivariogram-Based Block Size (",
      block_size, "m):\n",
      "OOF RMSE:", round(spatial_rmse, 3), "mW m-2 nm-1 sr-1\n",
      "OOF R-squared:", round(spatial_r2, 3), "\n")
  
  return(list(
    model = rf_model,
    max_train = max_train
  ))
}