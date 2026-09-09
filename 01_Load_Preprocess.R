# =============================================================================
# 01. DATA LOADING AND PREPROCESSING
# =============================================================================

# Load project setup
source("Setup.R")

# -----------------------------------------------------------------------------
# 01.1 Load source data and metadata
# -----------------------------------------------------------------------------

# Load HyPlant reflectance source data
hyplant_l1_refl_src <- rast("./Data/HyPlant/Reflectance/20230613-CKA-1359-0600-L5-S-DUAL_radiance_img_atm_pol-rect.bsq")
hyplant_l2_refl_src <- rast("./Data/HyPlant/Reflectance/20230613-CKA-1351-0600-L3-S-DUAL_radiance_img_atm_pol-rect.bsq")
hyplant_l3_refl_src <- rast("./Data/HyPlant/Reflectance/20200623-CKA-1329-0350-L5-W-DUAL_radiance_img_atm_pol-rect.bsq")

# Load HyPlant and DESIS SIF data
hyplant_l1_sfm_src <- rast("./Data/HyPlant/SIF/SFM/20230613-CKA-1359-0600-L5-S-FLUO_radiance_SFM_ALL-rect_georeferenced.tif")
hyplant_l2_sfm_src <- rast("./Data/HyPlant/SIF/SFM/20230613-CKA-1351-0600-L3-S-FLUO_radiance_SFM_ALL-rect.bil")
hyplant_l3_sfm_src <- rast("./Data/HyPlant/SIF/SFM/20200623-CKA-1329-0350-L5-W-FLUO_radiance_SIFO2A-rect_clip.tif")
desis_sfmnn_740_src <- rast("./Data/DESIS/SIF/SFMNN/DESIS-HSI-DT0867676616_004-20230613T123750-V0220-SPECTRAL_IMAGE__SFMNN_SIF.tif")

# Wavelength information
desis_wl <- read_excel("./Data/DESIS/L1C/L1C_Bands_Metadata.xlsx")[[4]]
hyplant_wl <- read.csv("./Data/HyPlant/Hyplant_wl.csv")[[1]]

# DESIS L2A product and bands metadata
desis_l2a <- rast("./Data/DESIS/L2A/DESIS-HSI-L2A-DT0867676616_004-20230613T123750-V0220-SPECTRAL_IMAGE.tif")
desis_l2a_bands_metadata <- read_excel("./Data/DESIS/L2A/L2A_Bands_Metadata.xlsx")


# -----------------------------------------------------------------------------
# 01.2  Preprocessing of HyPlant 2020 data
# -----------------------------------------------------------------------------

# Prepare HyPlant 2020 reflectance
hyplant_l3_refl_1m <- prep_spectral_data(
  hyplant_l3_refl_src, sensor = "HyPlant", no_data_value = 0
)

# Prepare the HyPlant 2020 SFM SIF map
hyplant_l3_sfm <- prep_sif_maps(
  hyplant_l3_sfm_src, hyplant_l3_refl_1m,
  sensor = "HyPlant", year = 2020
)

# Crop the reflectance map using the SFM map for flight line 3
hyplant_l3_refl_1m <- crop_raster(
  hyplant_l3_refl_1m, hyplant_l3_sfm
)

# Create study area polygon
study_area_l3 <- make_study_area_poly(
  hyplant_l3_refl_1m,
  output_path = "./Outputs/General/Study_Areas/study_area_l3.shp"
)

# Crop the flight line 3 SFM SIF map to the extent of the study area
hyplant_l3_sfm_1m <- crop_raster(
  hyplant_l3_sfm, study_area_l3,
  output_path = "./Outputs/HyPlant/SFM/HyPlant_l3_SFM_1m.tif"
)

# -----------------------------------------------------------------------------
# 01.3 Preprocessing of HyPlant 2023 data
# -----------------------------------------------------------------------------

# Prepare HyPlant 2023 reflectances
hyplant_l1_refl_1m <- prep_spectral_data(
  hyplant_l1_refl_src, sensor = "HyPlant", no_data_value = -999
)
hyplant_l2_refl_1m <- prep_spectral_data(
  hyplant_l2_refl_src, sensor = "HyPlant", no_data_value = -999
)

# Create study area polygons
study_area_l1 <- make_study_area_poly(
  hyplant_l1_refl_1m,
  output_path = "./Outputs/General/Study_Areas/Study_area_L1.shp"
)
study_area_l2 <- make_study_area_poly(
  hyplant_l2_refl_1m,
  output_path = "./Outputs/General/Study_Areas/Study_area_L2.shp"
)

# Combine L1 and L2 study areas into a multi-part polygon
study_area_l12 <- rbind(study_area_l1, study_area_l2)

# Crop reflectance maps to the extent of the study area
hyplant_l1_refl_1m <- crop_raster(
  hyplant_l1_refl_1m, study_area_l1,
  output_path = "./Outputs/HyPlant/Reflectance/HyPlant_L1_Reflectances_1m.tif"
)
hyplant_l2_refl_1m <- crop_raster(
  hyplant_l2_refl_1m, study_area_l2,
  output_path = "./Outputs/HyPlant/Reflectance/HyPlant_L2_Reflectances_1m.tif"
)

# -----------------------------------------------------------------------------
# 01.4 Classification
# -----------------------------------------------------------------------------

# Generate principal components from HyPlant reflectance
pc_l1 <- perform_pca(
  hyplant_l1_refl_1m, hyplant_wl, scale = TRUE,
  output_path = "./Outputs/General/Classification/HyPlant/PC_L1.tif"
)
pc_l2 <- perform_pca(
  hyplant_l2_refl_1m, hyplant_wl, scale = TRUE,
  output_path = "./Outputs/General/Classification/HyPlant/PC_L2.tif"
)
pc_l3 <- perform_pca(
  hyplant_l3_refl_1m, hyplant_wl, scale = TRUE,
  output_path = "./Outputs/General/Classification/HyPlant/PC_l3.tif"
)

# Read class labels shapefile
class_labels_l1 <- st_read("./Data/Classification/Four_classes_points_L1.shp")
class_labels_l2 <- st_read("./Data/Classification/Four_classes_points_L2.shp")
class_labels_l3 <- st_read("./Data/Classification/Four_classes_points_l3.shp")

# Perform supervised land cover classification using random forest
class_map_l1_1m <- rf_classification(
  pc_l1, class_labels_l1, study_area_l1,
  split_percentage = 0.75, num_trees = 100,
  output_path = "./Outputs/General/Classification/HyPlant/Class_map_L1_1m.tif"
)
class_map_l2_1m <- rf_classification(
  pc_l2, class_labels_l2, study_area_l2,
  split_percentage = 0.75, num_trees = 100,
  output_path = "./Outputs/General/Classification/HyPlant/Class_map_L2_1m.tif"
)
class_map_l3_1m <- rf_classification(
  pc_l3, class_labels_l3, study_area_l3,
  split_percentage = 0.75, num_trees = 100,
  output_path = "./Outputs/General/Classification/HyPlant/Class_map_l3_1m.tif"
)

# Aggregate classification maps to 30 m
class_map_l1_30m <- aggregate_crop(
  class_map_l1_1m, target_res = 30, fun = mode_resample,
  crop_area = study_area_l1,
  output_path = "./Outputs/General/Classification/HyPlant/Class_map_L1_30m.tif"
)
class_map_l2_30m <- aggregate_crop(
  class_map_l2_1m, target_res = 30, fun = mode_resample,
  crop_area = study_area_l2,
  output_path = "./Outputs/General/Classification/HyPlant/Class_map_L2_30m.tif"
)
class_map_l3_30m <- aggregate_crop(
  class_map_l3_1m, target_res = 30, fun = mode_resample,
  crop_area = study_area_l3,
  output_path = "./Outputs/General/Classification/HyPlant/Class_map_l3_30m.tif"
)

# Assign class names to classification values
class_levels <- data.frame(
  value = 1:4,
  class = c("Crops", "Mixed_Vegetation", "Non-Fluorescent", "Trees")
)

levels(class_map_l1_30m) <- class_levels
levels(class_map_l2_30m) <- class_levels
levels(class_map_l3_30m) <- class_levels


# -----------------------------------------------------------------------------
# 01.5 SIF and reflectance processing
# -----------------------------------------------------------------------------

# Prepare HyPlant and DESIS SIF maps
hyplant_l1_sfm_1m <- prep_sif_maps(
  hyplant_l1_sfm_src, class_map_l1_1m, study_area_l1,
  sensor = "HyPlant", year = 2023,
  output_path = "./Outputs/HyPlant/SFM/DESIS_L1_SFM_1m.tif"
)

hyplant_l2_sfm_1m <- prep_sif_maps(
  hyplant_l2_sfm_src, class_map_l2_1m, study_area_l2,
  sensor = "HyPlant", year = 2023,
  output_path = "./Outputs/HyPlant/SFM/HyPlant_L2_SFM_1m.tif"
)

desis_l1_sfmnn_30m <- prep_sif_maps(
  desis_sfmnn_740_src, class_map_l1_30m, study_area_l1,
  sensor = "DESIS", year = 2023,
  output_path = "./Outputs/DESIS/SIF/SFMNN/DESIS_L1_SFMNN_30m.tif"
)

desis_l2_sfmnn_30m <- prep_sif_maps(
  desis_sfmnn_740_src, class_map_l2_30m, study_area_l2,
  sensor = "DESIS", year = 2023,
  output_path = "./Outputs/DESIS/SIF/SFMNN/DESIS_L2_SFMNN_30m.tif"
)

# Resample HyPlant SIF maps to 30 m using the classification grids
hyplant_l1_sfm_30m <- resample_crop(
  hyplant_l1_sfm_1m, class_map_l1_30m, study_area_l1,
  method = "average",
  output_path = "./Outputs/HyPlant/SFM/HyPlant_L1_SFM_30m.tif"
)
hyplant_l2_sfm_30m <- resample_crop(
  hyplant_l2_sfm_1m, class_map_l2_30m, study_area_l2,
  method = "average",
  output_path = "./Outputs/HyPlant/SFM/HyPlant_L2_SFM_30m.tif"
)
hyplant_l3_sfm_30m <- resample_crop(
  hyplant_l3_sfm_1m, class_map_l3_30m, study_area_l3,
  method = "average",
  output_path = "./Outputs/HyPlant/SFM/HyPlant_l3_SFM_30m.tif"
)

# Resample reflectance rasters to 30 m
hyplant_l1_refl_30m <- resample_crop(
  hyplant_l1_refl_1m, class_map_l1_30m, study_area_l1,
  method = "average",
  output_path = "./Outputs/HyPlant/Reflectance/HyPlant_L1_Reflectances_30m.tif"
)
hyplant_l2_refl_30m <- resample_crop(
  hyplant_l2_refl_1m, class_map_l2_30m, study_area_l2,
  method = "average",
  output_path = "./Outputs/HyPlant/Reflectance/HyPlant_L2_Reflectances_30m.tif"
)
hyplant_l3_refl_30m <- resample_crop(
  hyplant_l3_refl_1m, class_map_l3_30m, study_area_l3,
  method = "average",
  output_path = "./Outputs/HyPlant/Reflectance/HyPlant_l3_Reflectances_30m.tif"
)

# -----------------------------------------------------------------------------
# 01.6 DESIS reflectance preprocessing
# -----------------------------------------------------------------------------

# Calculate the top-of-canopy reflectance for DESIS
desis_r_toc_30m <- prep_spectral_data(
  desis_l2a,
  sensor = "DESIS",
  product = "reflectance",
  bands_metadata = desis_l2a_bands_metadata,
  output_path = "./Outputs/DESIS/L2A/DESIS_TOC_Reflectance.tif"
)

# Resample and crop DESIS reflectance to the study areas
desis_l1_r_toc_30m <- resample_crop(
  desis_r_toc_30m, class_map_l1_30m, study_area_l1,
  method = "near",
  output_path = "./Outputs/DESIS/L2A/DESIS_L1_TOC_Reflectance.tif"
)

desis_l2_r_toc_30m <- resample_crop(
  desis_r_toc_30m, class_map_l2_30m, study_area_l2,
  method = "near",
  output_path = "./Outputs/DESIS/L2A/DESIS_L2_TOC_Reflectance.tif"
)