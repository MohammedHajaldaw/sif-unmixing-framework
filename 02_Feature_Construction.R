# =============================================================================
# 02. FEATURE CONSTRUCTION
# =============================================================================

# Run the required workflow stage for feature construction
source("Setup.R")
source("01_Load_Preprocess.R")

# Calculate DESIS vegetation indices
desis_l1_vi_30m <- calc_veg_indices(desis_l1_r_toc_30m, desis_wl, scale = FALSE,
  output_path = "./Outputs/General/Classification/DESIS/Unmixing_Indices_L1_1m.tif")
desis_l2_vi_30m <- calc_veg_indices(desis_l2_r_toc_30m, desis_wl, scale = FALSE,
  output_path = "./Outputs/General/Classification/DESIS/Unmixing_Indices_L2_1m.tif")

# Calculate HyPlant vegetation indices
hyplant_l1_vi_1m <- calc_veg_indices(hyplant_l1_refl_1m, hyplant_wl, scale = FALSE,
  output_path = "./Outputs/General/Classification/HyPlant/Unmixing_Indices_L1_1m.tif")
hyplant_l2_vi_1m <- calc_veg_indices(hyplant_l2_refl_1m, hyplant_wl, scale = FALSE,
  output_path = "./Outputs/General/Classification/HyPlant/Unmixing_Indices_L2_1m.tif")
hyplant_l3_vi_1m <- calc_veg_indices(hyplant_l3_refl_1m, hyplant_wl, scale = FALSE,
  output_path = "./Outputs/General/Classification/HyPlant/Unmixing_Indices_l3_1m.tif")

# Select the first 4 predictor layers
desis_l1_vi_30m <- desis_l1_vi_30m[, , 1:4, drop = FALSE]
desis_l2_vi_30m <- desis_l2_vi_30m[, , 1:4, drop = FALSE]
hyplant_l1_vi_1m <- hyplant_l1_vi_1m[, , 1:4, drop = FALSE]
hyplant_l2_vi_1m <- hyplant_l2_vi_1m[, , 1:4, drop = FALSE]
hyplant_l3_vi_1m <- hyplant_l3_vi_1m[, , 1:4, drop = FALSE]

# Remove the non-fluorescent class from the classification maps
reclass_matrix <- matrix(c(3, NA), ncol = 2, byrow = TRUE)
class_map_l1_1m <- classify(class_map_l1_1m, reclass_matrix)
class_map_l2_1m <- classify(class_map_l2_1m, reclass_matrix)
class_map_l3_1m <- classify(class_map_l3_1m, reclass_matrix)

# Update levels (drop "Non-Fluorescent")
levels(class_map_l1_1m) <- list(data.frame(value = c(1, 2, 4), class = c("Crops", "Mixed_Vegetation", "Trees")))
levels(class_map_l2_1m) <- list(data.frame(value = c(1, 2, 4), class = c("Crops", "Mixed_Vegetation", "Trees")))
levels(class_map_l3_1m) <- list(data.frame(value = c(1, 2, 4), class = c("Crops", "Mixed_Vegetation", "Trees")))

# Generate classification maps at 30m spatial resolution without the non-fluorescent class
class_map_l1_30m <- aggregate_crop_class(class_map_l1_1m, target_res = 30, fun = mode_resample, crop_area = study_area_l1)
class_map_l2_30m <- aggregate_crop_class(class_map_l2_1m, target_res = 30, fun = mode_resample, crop_area = study_area_l2)
class_map_l3_30m <- aggregate_crop_class(class_map_l3_1m, target_res = 30, fun = mode_resample, crop_area = study_area_l3)

# Aggregate classification maps to 60 m
class_map_l1_60m <- aggregate_crop_class(class_map_l1_30m, target_res = 2, fun = mode_resample, crop_area = study_area_l1,
                                         output_path = "./Outputs/General/Classification/HyPlant/Class_map_L1_60m.tif")
class_map_l2_60m <- aggregate_crop_class(class_map_l2_30m, target_res = 2, fun = mode_resample, crop_area = study_area_l2,
                                         output_path = "./Outputs/General/Classification/HyPlant/Class_map_L2_60m.tif")
class_map_l3_60m <- aggregate_crop_class(class_map_l3_30m, target_res = 2, fun = mode_resample, crop_area = study_area_l3,
                                         output_path = "./Outputs/General/Classification/HyPlant/Class_map_l3_60m.tif")

# Aggregate classification maps to 300 m
class_map_l1_300m <- aggregate_crop_class(class_map_l1_30m, target_res = 10, fun = mode_resample, crop_area = study_area_l1,
                                          output_path = "./Outputs/General/Classification/HyPlant/Class_map_L1_300m.tif")
class_map_l2_300m <- aggregate_crop_class(class_map_l2_30m, target_res = 10, fun = mode_resample, crop_area = study_area_l2,
                                          output_path = "./Outputs/General/Classification/HyPlant/Class_map_L2_300m.tif")
class_map_l3_300m <- aggregate_crop_class(class_map_l3_30m, target_res = 10, fun = mode_resample, crop_area = study_area_l3,
                                          output_path = "./Outputs/General/Classification/HyPlant/Class_map_l3_300m.tif")

# Resample HyPlant VIs to 30 m
hyplant_l1_vi_30m <- resample_crop(hyplant_l1_vi_1m, class_map_l1_30m, study_area_l1, method = "average",
                                   output_path = "./Outputs/General/Classification/HyPlant/Unmixing_Indices_L1_30m.tif")
hyplant_l2_vi_30m <- resample_crop(hyplant_l2_vi_1m, class_map_l2_30m, study_area_l2, method = "average",
                                   output_path = "./Outputs/General/Classification/HyPlant/Unmixing_Indices_L2_30m.tif")
hyplant_l3_vi_30m <- resample_crop(hyplant_l3_vi_1m, class_map_l3_30m, study_area_l3, method = "average",
                                   output_path = "./Outputs/General/Classification/HyPlant/Unmixing_Indices_l3_30m.tif")

# Resample DESIS and HyPlant VIs to 60 m
desis_l1_vi_60m <- resample_crop(desis_l1_vi_30m, class_map_l1_60m, study_area_l1, method = "average",
                                 output_path = "./Outputs/General/Classification/DESIS/Unmixing_Indices_L1_60m.tif")
desis_l2_vi_60m <- resample_crop(desis_l2_vi_30m, class_map_l2_60m, study_area_l2, method = "average",
                                 output_path = "./Outputs/General/Classification/DESIS/Unmixing_Indices_L2_60m.tif")
hyplant_l1_vi_60m <- resample_crop(hyplant_l1_vi_30m, class_map_l1_60m, study_area_l1, method = "average",
                                   output_path = "./Outputs/General/Classification/HyPlant/Unmixing_Indices_L1_60m.tif")
hyplant_l2_vi_60m <- resample_crop(hyplant_l2_vi_30m, class_map_l2_60m, study_area_l2, method = "average",
                                   output_path = "./Outputs/General/Classification/HyPlant/Unmixing_Indices_L2_60m.tif")
hyplant_l3_vi_60m <- resample_crop(hyplant_l3_vi_30m, class_map_l3_60m, study_area_l3, method = "average",
                                   output_path = "./Outputs/General/Classification/HyPlant/Unmixing_Indices_l3_60m.tif")

# Resample DESIS and HyPlant VIs to 300 m
desis_l1_vi_300m <- resample_crop(desis_l1_vi_30m, class_map_l1_300m, study_area_l1, method = "average",
                                  output_path = "./Outputs/General/Classification/DESIS/Unmixing_Indices_L1_300m.tif")
desis_l2_vi_300m <- resample_crop(desis_l2_vi_30m, class_map_l2_300m, study_area_l2, method = "average",
                                  output_path = "./Outputs/General/Classification/DESIS/Unmixing_Indices_L2_300m.tif")
hyplant_l1_vi_300m <- resample_crop(hyplant_l1_vi_30m, class_map_l1_300m, study_area_l1, method = "average",
                                    output_path = "./Outputs/General/Classification/HyPlant/Unmixing_Indices_L1_300m.tif")
hyplant_l2_vi_300m <- resample_crop(hyplant_l2_vi_30m, class_map_l2_300m, study_area_l2, method = "average",
                                    output_path = "./Outputs/General/Classification/HyPlant/Unmixing_Indices_L2_300m.tif")
hyplant_l3_vi_300m <- resample_crop(hyplant_l3_vi_30m, class_map_l3_300m, study_area_l3, method = "average",
                                    output_path = "./Outputs/General/Classification/HyPlant/Unmixing_Indices_l3_300m.tif")

# Resample DESIS and HyPlant SIF to 60 m
desis_l1_sfmnn_60m <- resample_crop(desis_l1_sfmnn_30m, class_map_l1_60m, study_area_l1, method = "average",
                                    output_path = "./Outputs/DESIS/SIF/SFMNN/DESIS_L1_SFMNN_60m.tif")
desis_l2_sfmnn_60m <- resample_crop(desis_l2_sfmnn_30m, class_map_l2_60m, study_area_l2, method = "average",
                                    output_path = "./Outputs/DESIS/SIF/SFMNN/DESIS_L2_SFMNN_60m.tif")
hyplant_l1_sfm_60m <- resample_crop(hyplant_l1_sfm_30m, class_map_l1_60m, study_area_l1, method = "average",
                                    output_path = "./Outputs/HyPlant/SFM/HyPlant_L1_SFM_60m.tif")
hyplant_l2_sfm_60m <- resample_crop(hyplant_l2_sfm_30m, class_map_l2_60m, study_area_l2, method = "average",
                                    output_path = "./Outputs/HyPlant/SFM/HyPlant_L2_SFM_60m.tif")
hyplant_l3_sfm_60m <- resample_crop(hyplant_l3_sfm_30m, class_map_l3_60m, study_area_l3, method = "average",
                                    output_path = "./Outputs/HyPlant/SFM/HyPlant_l3_SFM_60m.tif")

# Resample DESIS and HyPlant SIF to 300 m
desis_l1_sfmnn_300m <- resample_crop(desis_l1_sfmnn_30m, class_map_l1_300m, study_area_l1, method = "average",
                                     output_path = "./Outputs/DESIS/SIF/SFMNN/DESIS_L1_SFMNN_300m.tif")
desis_l2_sfmnn_300m <- resample_crop(desis_l2_sfmnn_30m, class_map_l2_300m, study_area_l2, method = "average",
                                     output_path = "./Outputs/DESIS/SIF/SFMNN/DESIS_L2_SFMNN_300m.tif")
hyplant_l1_sfm_300m <- resample_crop(hyplant_l1_sfm_30m, class_map_l1_300m, study_area_l1, method = "average",
                                     output_path = "./Outputs/HyPlant/SFM/HyPlant_L1_SFM_300m.tif")
hyplant_l2_sfm_300m <- resample_crop(hyplant_l2_sfm_30m, class_map_l2_300m, study_area_l2, method = "average",
                                     output_path = "./Outputs/HyPlant/SFM/HyPlant_L2_SFM_300m.tif")
hyplant_l3_sfm_300m <- resample_crop(hyplant_l3_sfm_30m, class_map_l3_300m, study_area_l3, method = "average",
                                     output_path = "./Outputs/HyPlant/SFM/HyPlant_l3_SFM_300m.tif")

# Calculate class fractions at 30 m
class_frac_l1_30m <- calc_class_frac(class_map_l1_1m, desis_l1_sfmnn_30m, study_area_l1)
class_frac_l2_30m <- calc_class_frac(class_map_l2_1m, desis_l2_sfmnn_30m, study_area_l2)
class_frac_l3_30m <- calc_class_frac(class_map_l3_1m, class_map_l3_30m, study_area_l3)

# Calculate class fractions at 60m spatial resolution for 2023
class_frac_l1_60m <- calc_class_frac(class_map_l1_30m, desis_l1_sfmnn_60m, study_area_l1)
class_frac_l2_60m <- calc_class_frac(class_map_l2_30m, desis_l2_sfmnn_60m, study_area_l2)

# Calculate class fractions at 300m spatial resolution in 2023
class_frac_l1_300m <- calc_class_frac(class_map_l1_30m, desis_l1_sfmnn_300m, study_area_l1)
class_frac_l2_300m <- calc_class_frac(class_map_l2_30m, desis_l2_sfmnn_300m, study_area_l2)

# Calculate class-specific SIF contributions at 30 m
class_sif_cont_l1_sfm_30m <- agg_sif_class(class_map_l1_1m, hyplant_l1_sfm_1m, class_map_l1_30m, class_frac_l1_30m, study_area_l1)
class_sif_cont_l2_sfm_30m <- agg_sif_class(class_map_l2_1m, hyplant_l2_sfm_1m, class_map_l2_30m, class_frac_l2_30m, study_area_l2)
class_sif_cont_l3_sfm_30m <- agg_sif_class(class_map_l3_1m, hyplant_l3_sfm_1m, class_map_l3_30m, class_frac_l3_30m, study_area_l3)

# Calculate class-specific SIF contributions at 60 m
class_sif_cont_l1_sfmnn_60m <- agg_sif_class(class_map_l1_30m, desis_l1_sfmnn_30m, class_map_l1_60m, class_frac_l1_60m, study_area_l1)
class_sif_cont_l2_sfmnn_60m <- agg_sif_class(class_map_l2_30m, desis_l2_sfmnn_30m, class_map_l2_60m, class_frac_l2_60m, study_area_l2)
class_sif_cont_l1_sfm_60m <- agg_sif_class(class_map_l1_30m, hyplant_l1_sfm_30m, class_map_l1_60m, class_frac_l1_60m, study_area_l1)
class_sif_cont_l2_sfm_60m <- agg_sif_class(class_map_l2_30m, hyplant_l2_sfm_30m, class_map_l2_60m, class_frac_l2_60m, study_area_l2)

# Calculate class-specific SIF contributions at 300 m
class_sif_cont_l1_sfmnn_300m <- agg_sif_class(class_map_l1_30m, desis_l1_sfmnn_30m, class_map_l1_300m, class_frac_l1_300m, study_area_l1)
class_sif_cont_l2_sfmnn_300m <- agg_sif_class(class_map_l2_30m, desis_l2_sfmnn_30m, class_map_l2_300m, class_frac_l2_300m, study_area_l2)
class_sif_cont_l1_sfm_300m <- agg_sif_class(class_map_l1_30m, hyplant_l1_sfm_30m, class_map_l1_300m, class_frac_l1_300m, study_area_l1)
class_sif_cont_l2_sfm_300m <- agg_sif_class(class_map_l2_30m, hyplant_l2_sfm_30m, class_map_l2_300m, class_frac_l2_300m, study_area_l2)

# Select pure pixels using the 85% purity threshold
pure_class_l1_30m <- get_pure_class(class_map_l1_30m, class_frac_l1_30m, threshold = 0.85)
pure_class_l2_30m <- get_pure_class(class_map_l2_30m, class_frac_l2_30m, threshold = 0.85)

# Aggregate pure classification maps to 60 m for model development
pure_class_l1_60m <- aggregate_crop_class(pure_class_l1_30m, target_res = 2, fun = mode_resample, crop_area = study_area_l1,
                                          output_path = "./Outputs/General/Classification/HyPlant/Class_map_L1_pure_60m.tif")
pure_class_l2_60m <- aggregate_crop_class(pure_class_l2_30m, target_res = 2, fun = mode_resample, crop_area = study_area_l2,
                                          output_path = "./Outputs/General/Classification/HyPlant/Class_map_L2_pure_60m.tif")

# Filter DESIS predictors to pure pixels
desis_l1_vi_pure_30m <- mask(desis_l1_vi_30m, pure_class_l1_30m)
desis_l2_vi_pure_30m <- mask(desis_l2_vi_30m, pure_class_l2_30m)

# Resample pure DESIS predictors to 60 m
desis_l1_vi_pure_60m <- resample_crop(desis_l1_vi_pure_30m, pure_class_l1_60m, study_area_l1, method = "average",
                                      output_path = "./Outputs/General/Classification/DESIS/Unmixing_Indices_L1_Pure_60m.tif")
desis_l2_vi_pure_60m <- resample_crop(desis_l2_vi_pure_30m, pure_class_l2_60m, study_area_l2, method = "average",
                                      output_path = "./Outputs/General/Classification/DESIS/Unmixing_Indices_L2_Pure_60m.tif")

class_frac_l1_pure_60m <- calc_class_frac(pure_class_l1_30m, pure_class_l1_60m, study_area_l1)
class_frac_l2_pure_60m <- calc_class_frac(pure_class_l2_30m, pure_class_l2_60m, study_area_l2)

# Filter DESIS SIF to pure pixels
desis_l1_sfmnn_pure_30m <- mask(desis_l1_sfmnn_30m, pure_class_l1_30m)
desis_l2_sfmnn_pure_30m <- mask(desis_l2_sfmnn_30m, pure_class_l2_30m)

# Resample pure DESIS SIF to 60 m
desis_l1_sfmnn_pure_60m <- resample_crop(desis_l1_sfmnn_pure_30m, pure_class_l1_60m, study_area_l1, method = "average",
                                         output_path = "./Outputs/DESIS/SIF/SFMNN/DESIS_L1_SFMNN_Pure_60m.tif")
desis_l2_sfmnn_pure_60m <- resample_crop(desis_l2_sfmnn_pure_30m, pure_class_l2_60m, study_area_l2, method = "average",
                                         output_path = "./Outputs/DESIS/SIF/SFMNN/DESIS_L2_SFMNN_Pure_60m.tif")

# Calculate pure-pixel class SIF contributions at 60 m
class_sif_cont_l1_sfmnn_pure_60m <- agg_sif_class(pure_class_l1_30m, desis_l1_sfmnn_pure_30m, pure_class_l1_60m, class_frac_l1_pure_60m, study_area_l1)
class_sif_cont_l2_sfmnn_pure_60m <- agg_sif_class(pure_class_l2_30m, desis_l2_sfmnn_pure_30m, pure_class_l2_60m, class_frac_l2_pure_60m, study_area_l2)

