# Function to aggregate SIF by vegetation class
agg_sif_class <- function(class_map, sif_raster, target_raster, fractions_raster, study_area) {
  
  # Overview:
  #   Calculates mean SIF for each vegetation class at the target resolution
  #   and weights the result by class fractions.
  
  # Requires:
  #   - class_map (SpatRaster): High-resolution vegetation classification.
  #   - sif_raster (SpatRaster): High-resolution SIF raster.
  #   - target_raster (SpatRaster): Raster defining the target resolution and grid.
  #   - fractions_raster (SpatRaster): Vegetation class fraction layers.
  #   - study_area (SpatVector): Study area boundary.
  
  # Effects:
  #   - Removes negative SIF values.
  #   - Calculates mean SIF for each vegetation class.
  #   - Resamples SIF to the target resolution and weights it by class fractions.
  #   - Crops and masks the result to the study area.
  
  # Returns:
  #   - A SpatRaster containing fraction-weighted mean SIF for each vegetation class.
  
  # Get unique vegetation classes.
  classes <- sort(unique(na.omit(values(class_map))))
  class_labels <- c("1" = "Crops", "2" = "Mixed_Vegetation", "4" = "Trees")
  
  # Remove negative SIF values.
  sif_raster[sif_raster < 0] <- NA
  
  # Calculate mean SIF for each vegetation class.
  sif_layers <- lapply(classes, function(class) {
    sif_class <- ifel(class_map == class, sif_raster, NA)
    sif_mean <- resample(sif_class, target_raster, method = "average")
    names(sif_mean) <- paste0(class_labels[as.character(class)], " Mean SIF")
    return(sif_mean)
  })
  
  # Combine and weight class-specific SIF by class fractions.
  result_raster <- rast(sif_layers) * fractions_raster
  values(result_raster)[is.na(values(result_raster))] <- 0
  
  # Crop and mask the result to the study area.
  result_raster <- crop(result_raster, study_area, touches = FALSE, mask = TRUE)
  
  return(result_raster)
}