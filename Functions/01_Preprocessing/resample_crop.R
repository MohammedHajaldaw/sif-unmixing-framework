# Function to resample and crop a raster
resample_crop <- function(input_raster, target_raster, crop_area,
                          method = "near", output_path = NULL) {
  
  # Overview:
  #   Resamples an input raster to match a target raster, then crops and
  #   masks it to a specified study area.
  
  # Requires:
  #   - input_raster (SpatRaster): Raster to be resampled and cropped.
  #   - target_raster (SpatRaster): Raster defining the target resolution and extent.
  #   - crop_area (SpatRaster or SpatVector): Study area used for cropping and masking.
  #   - method (character): Resampling method, default is "near".
  #   - output_path (character, optional): Output file path.
  
  # Effects:
  #   - Resamples the input raster to the target raster.
  #   - Crops and masks the result to the specified study area.
  #   - Optionally saves the processed raster.
  
  # Returns:
  #   - A SpatRaster containing the resampled and cropped raster.
  
  # Resample the input raster to match the target raster.
  resampled_raster <- resample(input_raster, target_raster, method = method)
  
  # Crop and mask the resampled raster to the study area.
  cropped_raster <- crop(resampled_raster, crop_area, touches = FALSE, mask = TRUE
  )
  
  # Optionally save the processed raster.
  if (!is.null(output_path))
    writeRaster(cropped_raster, output_path, overwrite = TRUE)
  
  return(cropped_raster)
}
