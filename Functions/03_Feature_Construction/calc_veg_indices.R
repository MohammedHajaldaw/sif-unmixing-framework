# Function to calculate vegetation indices from reflectance data
calc_veg_indices <- function(reflectance_raster, wl, scale = TRUE, output_path = NULL) {
  
  # Overview:
  #   Calculates NDVI, REP, PRI, and WBI from a hyperspectral reflectance raster.
  
  # Requires:
  #   - reflectance_raster (SpatRaster): Hyperspectral reflectance raster.
  #   - wl (numeric): Wavelength values corresponding to the raster bands.
  #   - scale (logical): Retained for compatibility; not used in the calculation.
  #   - output_path (character, optional): File path for the output raster.
  
  # Effects:
  #   - Extracts reflectance bands at the required wavelengths.
  #   - Calculates four vegetation indices.
  #   - Replaces invalid values with NA.
  #   - Optionally saves the index raster.
  
  # Returns:
  #   - A SpatRaster containing NDVI, REP, PRI, and WBI.
  
  # Get the indices of the required reflectance bands.
  idx_530 <- which.min(abs(wl - 530))
  idx_531 <- which.min(abs(wl - 531))
  idx_570 <- which.min(abs(wl - 570))
  idx_670 <- which.min(abs(wl - 670))
  idx_700 <- which.min(abs(wl - 700))
  idx_740 <- which.min(abs(wl - 740))
  idx_800 <- which.min(abs(wl - 800))
  idx_900 <- which.min(abs(wl - 900))
  idx_970 <- which.min(abs(wl - 970))
  
  # Extract reflectance at the required wavelengths.
  R_530 <- reflectance_raster[, , idx_530, drop = FALSE]
  R_531 <- reflectance_raster[, , idx_531, drop = FALSE]
  R_570 <- reflectance_raster[, , idx_570, drop = FALSE]
  R_670 <- reflectance_raster[, , idx_670, drop = FALSE]
  R_700 <- reflectance_raster[, , idx_700, drop = FALSE]
  R_740 <- reflectance_raster[, , idx_740, drop = FALSE]
  R_800 <- reflectance_raster[, , idx_800, drop = FALSE]
  R_900 <- reflectance_raster[, , idx_900, drop = FALSE]
  R_970 <- reflectance_raster[, , idx_970, drop = FALSE]
  
  # Calculate vegetation indices.
  NDVI <- (R_800 - R_670) / (R_800 + R_670)
  Ri <- (R_670 + R_800) / 2
  REP <- 700 + 40 * ((Ri - R_700) / (R_740 - R_700))
  PRI <- (R_570 - R_531) / (R_570 + R_530)
  WBI <- R_970 / R_900
  
  # Replace invalid index values with NA.
  NDVI[is.infinite(NDVI) | is.nan(NDVI)] <- NA
  REP[is.infinite(REP) | is.nan(REP)] <- NA
  PRI[is.infinite(PRI) | is.nan(PRI)] <- NA
  WBI[is.infinite(WBI) | is.nan(WBI)] <- NA
  
  # Stack and name the vegetation indices.
  index_stack <- c(NDVI, REP, PRI, WBI)
  names(index_stack) <- c("NDVI", "REP", "PRI", "WBI")
  
  # Optionally save the vegetation index raster.
  if (!is.null(output_path))
    writeRaster(index_stack, output_path, overwrite = TRUE)
  
  return(index_stack)
}
