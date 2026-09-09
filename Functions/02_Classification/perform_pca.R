# Function to perform PCA on reflectance data
perform_pca <- function(reflectance_raster, wl, scale = TRUE, output_path = NULL) {
  
  # Overview:
  #   Calculates vegetation indices, combines them with 400–1000 nm reflectance,
  #   and performs PCA using a random pixel sample.
  
  # Requires:
  #   - reflectance_raster (SpatRaster): Hyperspectral reflectance raster.
  #   - wl (numeric): Wavelength values corresponding to raster bands.
  #   - scale (logical): Whether to standardize variables before PCA.
  #   - output_path (character, optional): Output file path.
  
  # Effects:
  #   - Extracts reflectance bands at selected wavelengths.
  #   - Calculates vegetation indices and removes invalid values.
  #   - Performs PCA and reports cumulative variance and loadings.
  #   - Applies the first nine PCA components to the full raster.
  #   - Optionally saves the PCA raster.
  
  # Returns:
  #   - A SpatRaster containing the first nine principal components.
  
  # Get the indices of the required reflectance bands.
  idx_480 <- which.min(abs(wl - 480))
  idx_530 <- which.min(abs(wl - 530))
  idx_531 <- which.min(abs(wl - 531))
  idx_550 <- which.min(abs(wl - 550))
  idx_570 <- which.min(abs(wl - 570))
  idx_670 <- which.min(abs(wl - 670))
  idx_681 <- which.min(abs(wl - 681))
  idx_700 <- which.min(abs(wl - 700))
  idx_709 <- which.min(abs(wl - 709))
  idx_710 <- which.min(abs(wl - 710))
  idx_740 <- which.min(abs(wl - 740))
  idx_750 <- which.min(abs(wl - 750))
  idx_754 <- which.min(abs(wl - 754))
  idx_800 <- which.min(abs(wl - 800))
  idx_900 <- which.min(abs(wl - 900))
  idx_970 <- which.min(abs(wl - 970))
  
  # Extract reflectance at the required wavelengths.
  R_480 <- reflectance_raster[, , idx_480, drop = FALSE]
  R_530 <- reflectance_raster[, , idx_530, drop = FALSE]
  R_531 <- reflectance_raster[, , idx_531, drop = FALSE]
  R_550 <- reflectance_raster[, , idx_550, drop = FALSE]
  R_570 <- reflectance_raster[, , idx_570, drop = FALSE]
  R_670 <- reflectance_raster[, , idx_670, drop = FALSE]
  R_681 <- reflectance_raster[, , idx_681, drop = FALSE]
  R_700 <- reflectance_raster[, , idx_700, drop = FALSE]
  R_709 <- reflectance_raster[, , idx_709, drop = FALSE]
  R_710 <- reflectance_raster[, , idx_710, drop = FALSE]
  R_740 <- reflectance_raster[, , idx_740, drop = FALSE]
  R_750 <- reflectance_raster[, , idx_750, drop = FALSE]
  R_754 <- reflectance_raster[, , idx_754, drop = FALSE]
  R_800 <- reflectance_raster[, , idx_800, drop = FALSE]
  R_900 <- reflectance_raster[, , idx_900, drop = FALSE]
  R_970 <- reflectance_raster[, , idx_970, drop = FALSE]
  
  # Calculate vegetation indices.
  SR <- R_800 / R_670
  NDVI <- (R_800 - R_670) / (R_800 + R_670)
  NDVIre <- (R_750 - R_710) / (R_750 + R_710)
  EVI <- 2.5 * ((R_800 - R_670) / (R_800 + 6 * R_670 - 7.5 * R_480 + 1))
  Ri <- (R_670 + R_800) / 2
  REP <- 700 + 40 * ((Ri - R_700) / (R_740 - R_700))
  MTCI <- (R_754 - R_709) / (R_709 - R_681)
  TCARI <- 3 * ((R_700 - R_670) - 0.2 * (R_700 - R_550) * (R_700 / R_670))
  PRI <- (R_570 - R_531) / (R_570 + R_530)
  cPRI <- PRI - 0.15 * (1 - exp(-0.5 * SR))
  WBI <- R_970 / R_900
  OSAVI <- (1 + 0.16) * ((R_800 - R_670) / (R_800 + R_670 + 0.16))
  
  # Replace invalid index values with NA.
  SR[is.infinite(SR) | is.nan(SR)] <- NA
  NDVI[is.infinite(NDVI) | is.nan(NDVI)] <- NA
  NDVIre[is.infinite(NDVIre) | is.nan(NDVIre)] <- NA
  EVI[is.infinite(EVI) | is.nan(EVI)] <- NA
  REP[is.infinite(REP) | is.nan(REP)] <- NA
  MTCI[is.infinite(MTCI) | is.nan(MTCI)] <- NA
  TCARI[is.infinite(TCARI) | is.nan(TCARI)] <- NA
  PRI[is.infinite(PRI) | is.nan(PRI)] <- NA
  cPRI[is.infinite(cPRI) | is.nan(cPRI)] <- NA
  WBI[is.infinite(WBI) | is.nan(WBI)] <- NA
  OSAVI[is.infinite(OSAVI) | is.nan(OSAVI)] <- NA
  
  # Stack and name the vegetation indices.
  index_stack <- c(SR, NDVI, NDVIre, EVI, REP, MTCI, TCARI, PRI, cPRI, WBI, OSAVI)
  names(index_stack) <- c(
    "SR", "NDVI", "NDVIre", "EVI", "REP", "MTCI",
    "TCARI", "PRI", "cPRI", "WBI", "OSAVI"
  )
  
  # Extract reflectance bands from 400 to 1000 nm.
  idx_range <- which.min(abs(wl - 400)):which.min(abs(wl - 1000))
  R_400_to_1000 <- reflectance_raster[, , idx_range, drop = FALSE]
  
  # Combine reflectance and vegetation indices for PCA.
  pca_stack <- c(R_400_to_1000, index_stack)
  
  # Sample pixels and perform PCA.
  sampled_points <- spatSample(
    pca_stack, 10000, na.rm = TRUE,
    as.df = TRUE, method = "random"
  )
  pca_result <- prcomp(sampled_points, scale. = scale)
  
  # Calculate cumulative variance explained.
  variance_explained <- summary(pca_result)$importance[2, ]
  cum_variance <- cumsum(variance_explained)
  n_components <- which(cum_variance >= 0.99)[1]
  
  # Print variance explained by the selected components.
  print(paste(
    "Cumulative variance of the selected components (", n_components, "): ",
    cum_variance[n_components], sep = ""
  ))
  
  # Assign input band names to the PCA loadings.
  rownames(pca_result$rotation) <- names(pca_stack)
  
  # Print the three strongest loadings for the first four components.
  cat("\nTop 3 contributing bands for the first 4 principal components:\n")
  for (i in 1:min(4, ncol(pca_result$rotation))) {
    cat(paste0("\nPC", i, ":\n"))
    sorted_loadings <- sort(abs(pca_result$rotation[, i]), decreasing = TRUE)
    print(round(head(sorted_loadings, 3), 4))
  }
  
  # Apply the first nine PCA components (use 1:n_components for threshold-based selection).
  principal_components <- predict(pca_stack, pca_result, index = 1:9)
  
  # Optionally save the PCA raster.
  if (!is.null(output_path))
    writeRaster(principal_components, output_path, overwrite = TRUE)
  
  return(principal_components)
}