# =============================================================================
# 03. MODEL DEVELOPMENT
# =============================================================================

# Run the required workflow stages for model development
source("Setup.R")
source("01_Load_Preprocess.R")
source("02_Feature_Construction.R")

# Merge L1 and L2 classification and fraction data
class_frac_l12_30m <- merge(class_frac_l1_30m, class_frac_l2_30m)
class_frac_l12_60m <- merge(class_frac_l1_60m, class_frac_l2_60m)
class_frac_l12_300m <- merge(class_frac_l1_300m, class_frac_l2_300m)

pure_class_map_l12_30m <- merge(pure_class_l1_30m, pure_class_l2_30m)

class_map_l12_30m <- merge(class_map_l1_30m, class_map_l2_30m)
class_map_l12_60m <- merge(class_map_l1_60m, class_map_l2_60m)
class_map_l12_300m <- merge(class_map_l1_300m, class_map_l2_300m)

# Prepare 60 m pure-pixel data for model development
class_sif_cont_l12_sfmnn_pure_60m <- merge( class_sif_cont_l1_sfmnn_pure_60m,  class_sif_cont_l2_sfmnn_pure_60m)
desis_l12_vi_pure_60m <- merge(desis_l1_vi_pure_60m, desis_l2_vi_pure_60m)
class_frac_l12_pure_60m <- merge(class_frac_l1_pure_60m, class_frac_l2_pure_60m)
desis_l12_sfmnn_pure_60m <- merge(desis_l1_sfmnn_pure_60m, desis_l2_sfmnn_pure_60m)

# Merge DESIS predictors and retrieved SIF
desis_l12_vi_30m <- merge(desis_l1_vi_30m, desis_l2_vi_30m)
desis_l12_vi_60m <- merge(desis_l1_vi_60m, desis_l2_vi_60m)
desis_l12_vi_300m <- merge(desis_l1_vi_300m, desis_l2_vi_300m)

desis_l12_sfmnn_30m <- merge(desis_l1_sfmnn_30m, desis_l2_sfmnn_30m)
desis_l12_sfmnn_60m <- merge(desis_l1_sfmnn_60m, desis_l2_sfmnn_60m)
desis_l12_sfmnn_300m <- merge(desis_l1_sfmnn_300m, desis_l2_sfmnn_300m)

# Merge HyPlant predictors and retrieved SIF for evaluation
class_sif_cont_l12_sfm_30m <- merge( class_sif_cont_l1_sfm_30m,  class_sif_cont_l2_sfm_30m)
class_sif_cont_l12_sfm_60m <- merge( class_sif_cont_l1_sfm_60m,  class_sif_cont_l2_sfm_60m)
class_sif_cont_l12_sfm_300m <- merge( class_sif_cont_l1_sfm_300m,  class_sif_cont_l2_sfm_300m)

hyplant_l12_vi_30m <- merge(hyplant_l1_vi_30m, hyplant_l2_vi_30m)
hyplant_l12_vi_60m <- merge(hyplant_l1_vi_60m, hyplant_l2_vi_60m)
hyplant_l12_vi_300m <- merge(hyplant_l1_vi_300m, hyplant_l2_vi_300m)

hyplant_l12_sfm_30m <- merge(hyplant_l1_sfm_30m, hyplant_l2_sfm_30m)
hyplant_l12_sfm_60m <- merge(hyplant_l1_sfm_60m, hyplant_l2_sfm_60m)
hyplant_l12_sfm_300m <- merge(hyplant_l1_sfm_300m, hyplant_l2_sfm_300m)

# Train unified RF model with spatial cross-validation
rf_model <- train_rf_spatial_cv(
  sif_cont_train =  class_sif_cont_l12_sfmnn_pure_60m,
  vi_raster = desis_l12_vi_pure_60m,
  frac = class_frac_l12_pure_60m,
  sif_ret = desis_l12_sfmnn_pure_60m,
  run_sensitivity = TRUE
)

# Train linear regression baseline with spatial cross-validation
lm_model <- train_lm_spatial_cv(
  sif_cont_train =  class_sif_cont_l12_sfmnn_pure_60m,
  vi_raster = desis_l12_vi_pure_60m,
  frac = class_frac_l12_pure_60m,
  sif_ret = desis_l12_sfmnn_pure_60m
)
