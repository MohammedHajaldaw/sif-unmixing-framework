# Function to calculate vegetation class fractions at a target resolution
calc_class_frac <- function(class_map, target_raster, study_area) {
  
  # Overview:
  #   Calculates the fraction of each vegetation class within the target
  #   raster resolution and crops the result to the study area.
  
  # Requires:
  #   - class_map (SpatRaster): High-resolution vegetation classification raster.
  #   - target_raster (SpatRaster): Raster defining the target resolution and grid.
  #   - study_area (SpatVector): Study area boundary.
  
  # Effects:
  #   - Calculates the fraction of each vegetation class using average resampling.
  #   - Crops and masks the fraction layers to the study area.
  
  # Returns:
  #   - A SpatRaster containing one fractional cover layer per vegetation class.
  
  # Get unique vegetation classes.
  classes <- sort(unique(na.omit(values(class_map))))
  class_labels <- c("1" = "Crops", "2" = "Mixed_Vegetation", "4" = "Trees")
  
  # Calculate the fraction of each class at the target resolution.
  fraction_layers <- lapply(classes, function(class) {
    fraction <- resample(class_map == class, target_raster, method = "average") # Resample using average method
    names(fraction) <- paste0(class_labels[as.character(class)], " fractions")  # Assign meaningful names
    return(fraction)
  })
  
  # Combine fraction layers and crop to the study area.
  fraction_raster <- rast(fraction_layers)
  fraction_raster <- crop(fraction_raster, study_area, touches = FALSE, mask = TRUE)
  
  # Return the final fraction raster.
  return(fraction_raster)
}