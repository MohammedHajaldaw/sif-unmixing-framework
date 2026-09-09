# Function to prepare spectral data for HyPlant or DESIS
prep_spectral_data <- function(input_raster, sensor = "HyPlant",
                               product = "reflectance", no_data_value = -999,
                               bands_metadata = NULL, output_path = NULL) {
  
  # Overview:
  #   Prepares HyPlant reflectance or DESIS reflectance/radiance data.
  
  # Requires:
  #   - input_raster (SpatRaster): Raw sensor raster.
  #   - sensor (character): "HyPlant" or "DESIS".
  #   - product (character): "reflectance" or "radiance".
  #   - no_data_value (numeric): HyPlant no-data value.
  #   - bands_metadata (list): DESIS metadata containing gain and offset values.
  #   - output_path (character, optional): Output file path.
  
  # Effects:
  #   - Replaces HyPlant no-data values and scales reflectance.
  #   - Applies DESIS gain and offset values.
  #   - Converts DESIS radiance to mW/m²/sr/nm.
  #   - Optionally saves the processed raster.
  
  # Returns:
  #   - A SpatRaster containing the prepared spectral data.
  
  # Validate sensor and product choices.
  if (!sensor %in% c("HyPlant", "DESIS"))
    stop("sensor must be 'HyPlant' or 'DESIS'.")
  
  if (!product %in% c("reflectance", "radiance"))
    stop("product must be 'reflectance' or 'radiance'.")
  
  # Radiance preparation is only available for DESIS.
  if (sensor == "HyPlant" && product == "radiance")
    stop("Radiance preparation is only available for DESIS.")
  
  # Prepare HyPlant reflectance.
  if (sensor == "HyPlant") {
    input_raster[input_raster == no_data_value] <- NA
    spectral_data <- input_raster / 10000
  }
  
  # Prepare DESIS reflectance or radiance using metadata calibration.
  if (sensor == "DESIS") {
    if (is.null(bands_metadata))
      stop("bands_metadata is required for DESIS.")
    
    # Extract band-specific gain and offset values from the metadata.
    band_gain <- bands_metadata[[2]]
    band_offset <- bands_metadata[[3]]
    
    # Apply gain and offset to each band.
    spectral_data <- rast(lapply(seq_len(nlyr(input_raster)), function(i)
      input_raster[[i]] * band_gain[i] + band_offset[i]))
    
    # Convert DESIS radiance to mW/m²/sr/nm.
    if (product == "radiance")
      spectral_data <- spectral_data * 10
  }
  
  # Optionally save the processed raster.
  if (!is.null(output_path))
    writeRaster(spectral_data, output_path, overwrite = TRUE)
  
  return(spectral_data)
}