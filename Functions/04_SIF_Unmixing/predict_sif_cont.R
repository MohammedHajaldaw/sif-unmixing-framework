# Function to predict SIF contribution using class-specific scaling factors
predict_sif_cont <- function(vi_raster, frac, sif_ret, model_fit, class_map,
                             predictor_set = c("Indices", "Fractions",
                                               "Class_id", "SIFret"),
                             use_scaling = TRUE) {
  
  # Overview:
  #   This function predicts the SIF contribution of each vegetation class.
  #   It calculates a scaling factor for each class based on the ratio
  #   between the maximum training and input SIF.
  
  # Requires:
  #   - vi_raster (SpatRaster): Vegetation index raster stack.
  #   - frac (SpatRaster): Raster stack containing class fractions.
  #   - sif_ret (SpatRaster): Retrieved SIF raster to be unmixed.
  #   - model_fit (list): Fitted model and training metadata.
  #   - class_map (SpatRaster): Vegetation class map.
  
  # Effects:
  #   - Calculates class-specific scaling factors.
  #   - Creates class-specific predictor rasters.
  #   - Generates model predictions.
  
  # Returns:
  #   - A SpatRaster containing the predicted SIF contribution per class.
  
  # Function to match a raster to the reference extent.
  match_extent <- function(r, study_area) {
    if (!terra::ext(r) == terra::ext(study_area)) r <- terra::crop(r, study_area)
    terra::mask(r, study_area)
  }
  
  # Ensure class_map is a SpatRaster and match input extents.
  if (!inherits(class_map, "SpatRaster")) stop("Error: class_map must be a SpatRaster.")
  vi_raster <- match_extent(vi_raster, frac[[1]])
  class_map <- match_extent(class_map, frac[[1]])
  sif_ret <- match_extent(sif_ret, frac[[1]])
  
  # Extract the fitted model and maximum training SIF contribution.
  model <- model_fit$model
  max_train <- model_fit$max_train
  
  # Get unique vegetation classes and initialize the final prediction.
  class_levels <- na.omit(terra::unique(class_map)[, 1])
  sif_cont_pred <- sif_ret * 0
  
  # Extract maximum retrieved SIF for each class.
  max_ret <- sapply(class_levels, function(x) {
    class_mask <- terra::ifel(class_map == x, 1, NA)
    max(terra::values(terra::mask(sif_ret, class_mask)), na.rm = TRUE)
  })
  
  # Scale retrieved SIF only when its maximum exceeds the training maximum.
  scale_fac <- ifelse(max_ret > max_train, max_train / max_ret, 1)
  message("Scaling factors:")
  print(round(scale_fac, 3))
  
  # Predict SIF contribution separately for each vegetation class.
  for (class_value in class_levels) {
    
    # Create a mask for the current class and skip classes without valid SIF.
    class_mask <- terra::ifel(class_map == class_value, 1, NA)
    sif_ret_masked <- sif_ret * class_mask
    if (sum(terra::values(sif_ret_masked), na.rm = TRUE) == 0) next
    
    # Get the scaling factor and apply it before prediction when selected.
    scale_i <- scale_fac[as.character(class_value)]
    sif_ret_scaled <- if (use_scaling) sif_ret_masked * scale_i else sif_ret_masked
    
    # Predict SIF contribution using the selected predictor groups.
    sif_cont_i <- NULL
    
    for (i in seq_len(nlyr(frac))) {
      class_id <- sif_ret_scaled
      class_id[!is.na(class_id)] <- i
      names(class_id) <- "class_id"
      
      # Build the predictor raster for the current class.
      pred_list <- list()
      if ("Indices" %in% predictor_set) pred_list <- c(pred_list, list(vi_raster))
      frac_i <- frac[[i]]
      names(frac_i) <- "frac"
      pred_list <- c(pred_list, list(frac_i))
      if ("Class_id" %in% predictor_set) pred_list <- c(pred_list, list(class_id))
      if ("SIFret" %in% predictor_set) {
        names(sif_ret_scaled) <- "SIFret"
        pred_list <- c(pred_list, list(sif_ret_scaled))
      }
      
      # Combine predictors and apply the fitted model.
      pred_i <- predict(do.call(c, pred_list), model, na.rm = TRUE)
      names(pred_i) <- paste0("SIFcont_", names(frac)[i])
      sif_cont_i <- if (is.null(sif_cont_i)) pred_i else c(sif_cont_i, pred_i)
    }
    
    # Reverse the scaling after prediction when scaling is used.
    sif_cont_i <- if (use_scaling) sif_cont_i / scale_i else sif_cont_i
    
    # Replace invalid predictions with zero and add them to the final output.
    values(sif_cont_i)[is.na(values(sif_cont_i)) | is.infinite(values(sif_cont_i))] <- 0
    sif_cont_pred <- sif_cont_pred + sif_cont_i
  }
  
  # Crop the final prediction to the study area.
  sif_cont_pred <- terra::crop(sif_cont_pred, sif_ret, touches = FALSE, mask = TRUE)
  
  return(sif_cont_pred)
}