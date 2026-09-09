# =============================================================================
# 04. PREDICTION AND PRODUCTS GENERATION
# =============================================================================

# Run the required workflow stages for prediction
source("Setup.R")
source("01_Load_Preprocess.R")
source("02_Feature_Construction.R")
source("03_Model_Development.R")

# Predict class-specific SIF contributions for DESIS 2023
desis_l12_sfmnn_pred_30m <- predict_sif_cont(desis_l12_vi_30m, class_frac_l12_30m, desis_l12_sfmnn_30m, rf_model, class_map_l12_30m)
desis_l12_sfmnn_pred_60m <- predict_sif_cont(desis_l12_vi_60m, class_frac_l12_60m, desis_l12_sfmnn_60m, rf_model, class_map_l12_60m)
desis_l12_sfmnn_pred_300m <- predict_sif_cont(desis_l12_vi_300m, class_frac_l12_300m, desis_l12_sfmnn_300m, rf_model, class_map_l12_300m)

# Reconstruct total SIF at 30 m resolution for DESIS
desis_l12_sfmnn_mix_30m <- sum(desis_l12_sfmnn_pred_30m)

# Prepare reconstructed and retrieved total SIF for vegetation-area evaluation
desis_l12_sfmnn_mix_30m_veg <- crop_raster(desis_l12_sfmnn_mix_30m, class_map_l12_30m)
desis_l12_sfmnn_30m_veg <- crop_raster(desis_l12_sfmnn_30m, class_map_l12_30m)

# Predict class-specific SIF contributions for HyPlant 2023
hyplant_l12_sfm_pred_30m <- predict_sif_cont(hyplant_l12_vi_30m, class_frac_l12_30m, hyplant_l12_sfm_30m, rf_model, class_map_l12_30m)
hyplant_l12_sfm_pred_60m <- predict_sif_cont(hyplant_l12_vi_60m, class_frac_l12_60m, hyplant_l12_sfm_60m, rf_model, class_map_l12_60m)
hyplant_l12_sfm_pred_300m <- predict_sif_cont(hyplant_l12_vi_300m, class_frac_l12_300m, hyplant_l12_sfm_300m, rf_model, class_map_l12_300m)

# Predict class-specific SIF contributions using linear regression
hyplant_l12_sfm_lm_pred_30m <- predict_sif_cont(hyplant_l12_vi_30m, class_frac_l12_30m, hyplant_l12_sfm_30m, lm_model, class_map_l12_30m)

# Predict class-specific SIF contributions for HyPlant 2020
hyplant_l3_sfm_pred_30m <- predict_sif_cont(hyplant_l3_vi_30m, class_frac_l3_30m, hyplant_l3_sfm_30m, rf_model, class_map_l3_30m)
hyplant_l3_sfm_lm_pred_30m <- predict_sif_cont(hyplant_l3_vi_30m, class_frac_l3_30m, hyplant_l3_sfm_30m, lm_model, class_map_l3_30m)

# Convert predicted SIF contributions to class-specific average SIF
desis_l12_sfmnn_prod_30m <- calculate_class_avg_sif(
  desis_l12_sfmnn_pred_30m, class_frac_l12_30m
)

# Save class-specific average SIF products
writeRaster(desis_l12_sfmnn_prod_30m, "./Outputs/DESIS_Unmixing_Products/SFMNN_l12_prod_30m.tif", overwrite = TRUE)

