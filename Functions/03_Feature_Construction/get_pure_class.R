# Function to identify pure vegetation classes
get_pure_class <- function(class_raster, fractions_raster, threshold = 0.85) {
  
  # Overview:
  #   Identifies pixels where the dominant vegetation class fraction meets
  #   the specified purity threshold and masks all other pixels.
  
  # Requires:
  #   - class_raster (SpatRaster): Vegetation classification raster.
  #   - fractions_raster (SpatRaster): Raster stack of vegetation class fractions.
  #   - threshold (numeric): Minimum fraction required for a pure pixel.
  
  # Effects:
  #   - Checks compatibility of the input rasters.
  #   - Calculates the maximum class fraction per pixel.
  #   - Masks pixels below the specified purity threshold.
  
  # Returns:
  #   - A SpatRaster containing the pure vegetation class pixels.
  
  # Check input compatibility.
  if (!inherits(class_raster, "SpatRaster") || !inherits(fractions_raster, "SpatRaster"))
    stop("Both inputs must be SpatRaster objects.")
  if (!all(dim(class_raster)[1:2] == dim(fractions_raster)[1:2]))
    stop("class_raster and fractions_raster must have the same spatial dimensions.")
  
  # Get the maximum class fraction per pixel.
  max_fraction <- max(fractions_raster)
  
  # Mask pixels below the purity threshold.
  pure_mask <- max_fraction >= threshold
  pure_class_map <- mask(class_raster, pure_mask, maskvalues = FALSE)
  
  return(pure_class_map)
}