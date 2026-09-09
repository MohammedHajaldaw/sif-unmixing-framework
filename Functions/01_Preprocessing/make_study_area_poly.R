# Function to create a study area polygon from a reflectance raster
make_study_area_poly <- function(hyplant_refl, output_path) {
  
  # Overview:
  #   Creates a study area polygon from a HyPlant reflectance raster.
  
  # Requires:
  #   - hyplant_refl (SpatRaster): HyPlant reflectance raster.
  #   - output_path (character): Output file path for the study area shapefile.
  
  # Effects:
  #   - Converts the first reflectance band to a dissolved polygon.
  #   - Fills holes and ensures valid attribute names.
  #   - Saves the resulting polygon as an ESRI Shapefile.
  
  # Returns:
  #   - A SpatVector containing the study area polygon.
  
  # Convert the first reflectance band to a dissolved polygon.
  study_area_poly <- as.polygons(hyplant_refl[[1]], dissolve = TRUE)
  
  # Fill holes in the study area polygon.
  study_area_filled <- fillHoles(study_area_poly)
  
  # Ensure valid and unique attribute names.
  names(study_area_filled) <- make.names(names(study_area_filled), unique = TRUE)
  
  # Save the study area polygon.
  writeVector(study_area_filled, output_path, filetype = "ESRI Shapefile", overwrite = TRUE)
  
  return(study_area_filled)
}
