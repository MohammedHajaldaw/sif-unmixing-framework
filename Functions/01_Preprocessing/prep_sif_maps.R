# Function to prepare SIF data for HyPlant or DESIS
prep_sif_maps <- function(sif_raster, target_raster, study_area = NULL,
                          sensor = "HyPlant", year = 2023, output_path = NULL) {
  
  # Overview:
  #   Prepares HyPlant SFM or DESIS SFMNN SIF data for the specified year.
  
  # Requires:
  #   - sif_raster (SpatRaster): Raw SIF raster.
  #   - target_raster (SpatRaster): Raster defining the target grid.
  #   - study_area (SpatVector, optional): Study area for cropping and masking.
  #   - sensor (character): "HyPlant" or "DESIS".
  #   - year (numeric): 2020 or 2023.
  #   - output_path (character, optional): Output file path.
  
  # Effects:
  #   - Extracts, scales, resamples, crops, and masks SIF data.
  #   - Applies sensor- and year-specific preprocessing.
  #   - Optionally saves the processed raster.
  
  # Returns:
  #   - A SpatRaster containing the prepared SIF data.
  
  if (!sensor %in% c("HyPlant", "DESIS")) stop("sensor must be 'HyPlant' or 'DESIS'.")
  if (!year %in% c(2020, 2023)) stop("year must be 2020 or 2023.")
  if (sensor == "DESIS" && year != 2023) stop("DESIS SIF data are only available for 2023.")
  if (!is.null(study_area) && !inherits(study_area, "SpatVector"))
    stop("study_area must be a SpatVector.")
  
  # Prepare HyPlant SFM SIF.
  if (sensor == "HyPlant") {
    sif <- sif_raster[[if (year == 2020) 11 else 2]]
    if (year == 2020) sif[sif == 32767] <- NA
    sif <- sif / 100
    sif[sif < -1 | sif > 5] <- NA
  }
  
  # Prepare DESIS SFMNN SIF and convert from 740 to 760 nm.
  if (sensor == "DESIS") {
    sif <- sif_raster[[1]] / 10000
    sif <- sif * 0.5093901
  }
  
  # Resample SIF to the target grid.
  sif <- resample(sif, target_raster, method = "near")
  
  # Crop to the study area or target raster extent.
  if (year == 2020 && sensor == "HyPlant")
    sif <- crop(sif, target_raster, touches = TRUE, mask = TRUE)
  else
    sif <- crop(sif, study_area, touches = FALSE, mask = TRUE)
  
  # Optionally save the processed raster.
  if (!is.null(output_path)) writeRaster(sif, output_path, overwrite = TRUE)
  
  return(sif)
}