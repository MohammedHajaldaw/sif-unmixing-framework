# Function to aggregate, crop, and assign vegetation classes
aggregate_crop_class <- function(input_raster, target_res, fun, crop_area, output_path = NULL) {
  
  # Overview:
  #   Aggregates and crops a vegetation classification raster, then assigns class labels.
  
  # Requires:
  #   - input_raster (SpatRaster): Vegetation classification raster.
  #   - target_res (numeric): Aggregation factor.
  #   - fun (function): Aggregation function, e.g. modal.
  #   - crop_area (SpatRaster or SpatVector): Area used for cropping and masking.
  #   - output_path (character, optional): Output file path.
  
  # Effects:
  #   - Aggregates the classification raster.
  #   - Crops and masks the result to the specified area.
  #   - Assigns labels to Crops, Mixed Vegetation, and Trees.
  #   - Optionally saves the processed raster.
  
  # Returns:
  #   - A SpatRaster containing the aggregated vegetation classification.
  
  # Aggregate the classification raster and crop it to the study area.
  raster <- aggregate(input_raster, fact = target_res, fun = fun)
  raster <- crop(raster, crop_area, touches = FALSE, mask = TRUE)
  
  # Assign labels to the three retained vegetation classes.
  levels(raster) <- list(data.frame(
    value = c(1, 2, 4),
    class = c("Crops", "Mixed_Vegetation", "Trees")
  ))
  
  # Optionally save the processed raster.
  if (!is.null(output_path)) writeRaster(raster, output_path, overwrite = TRUE)
  
  return(raster)
}