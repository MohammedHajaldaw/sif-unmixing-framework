# =============================================================================
# 06. SENSITIVITY ANALYSIS
# =============================================================================

# Run the required workflow stages for sensitivity analysis
source("Setup.R")
source("01_Load_Preprocess.R")
source("02_Feature_Construction.R")
source("03_Model_Development.R")

# Assess model performance across purity thresholds
assess_purity_thresholds(
  thresholds = c(.75, .80, .85, .90, .95),
  class_map_l1_30m, class_map_l2_30m,
  class_frac_l1_30m, class_frac_l2_30m,
  desis_l1_vi_30m, desis_l2_vi_30m,
  desis_l1_sfmnn_30m, desis_l2_sfmnn_30m,
  hyplant_l12_vi_60m, hyplant_l12_sfm_60m,
  class_frac_l12_60m, class_map_l12_60m,
  class_sif_cont_l12_sfm_60m,
  study_area_l1, study_area_l2
)

# Evaluate model performance across predictor sets
evaluate_predictor_sets(
  class_sif_cont_l12_sfmnn_pure_60m,
  desis_l12_vi_pure_60m,
  class_frac_l12_pure_60m,
  desis_l12_sfmnn_pure_60m,
  hyplant_l3_vi_30m,
  class_frac_l3_30m,
  hyplant_l3_sfm_30m,
  class_map_l3_30m,
  class_sif_cont_l3_sfm_30m
)

# Test predictor sensitivity to an extended 2020 HyPlant SIF range
evaluate_predictor_sets(
  class_sif_cont_l12_sfmnn_pure_60m,
  desis_l12_vi_pure_60m,
  class_frac_l12_pure_60m,
  desis_l12_sfmnn_pure_60m,
  hyplant_l3_vi_30m,
  class_frac_l3_30m,
  hyplant_l3_sfm_30m * 2,
  class_map_l3_30m,
  class_sif_cont_l3_sfm_30m * 2
)
