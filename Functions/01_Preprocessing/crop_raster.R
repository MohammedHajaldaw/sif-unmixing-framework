# Function to crop and mask a raster to a specified area
crop_raster <- function(input_raster, crop_area, output_path = NULL) {
  
  # Overview:
  #   Crops and masks a raster to a specified area.
  
  # Requires:
  #   - input_raster (SpatRaster): Raster to be cropped.
  #   - crop_area (SpatRaster or SpatVector): Area used for cropping and masking.
  #   - output_path (character, optional): Output file path.
  
  # Effects:
  #   - Crops the raster to the extent of the crop area.
  #   - Masks values outside the crop area.
  #   - Optionally saves the cropped raster.
  
  # Returns:
  #   - A SpatRaster containing the cropped and masked raster.
  
  # Crop and mask the raster to the specified area.
  raster_crop <- crop(input_raster, crop_area,touches = FALSE,mask = TRUE)
  
  # Optionally save the cropped raster.
  if (!is.null(output_path))
    writeRaster(raster_crop, output_path, overwrite = TRUE)
  
  return(raster_crop)
}
