# Function to aggregate and crop a raster
aggregate_crop <- function(input_raster, target_res, fun, crop_area, output_path = NULL) {
  
  # Overview:
  #   Aggregates a raster to a specified resolution, then crops and masks it.
  
  # Requires:
  #   - input_raster (SpatRaster): Raster to aggregate.
  #   - target_res (numeric): Aggregation factor.
  #   - fun (function): Aggregation function, e.g. mean, max, or sum.
  #   - crop_area (SpatRaster or SpatVector): Area used for cropping and masking.
  #   - output_path (character, optional): Output file path.
  
  # Effects:
  #   - Aggregates the input raster.
  #   - Crops and masks the result to the specified area.
  #   - Optionally saves the processed raster.
  
  # Returns:
  #   - A SpatRaster containing the aggregated and cropped raster.
  
  # Aggregate the raster and crop it to the study area.
  raster <- aggregate(input_raster, fact = target_res, fun = fun)
  raster <- crop(raster, crop_area, touches = FALSE, mask = TRUE)
  
  # Optionally save the processed raster.
  if (!is.null(output_path)) writeRaster(raster, output_path, overwrite = TRUE)
  
  return(raster)
}