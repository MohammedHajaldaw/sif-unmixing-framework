# =============================================================================
# 05. EVALUATION
# =============================================================================

# Run the required workflow stages for evaluation
source("Setup.R")
source("01_Load_Preprocess.R")
source("02_Feature_Construction.R")
source("03_Model_Development.R")
source("04_Prediction_Products.R")

# Compare reconstructed and retrieved total SIF within the vegetation area
sif_diff_veg <- desis_l12_sfmnn_mix_30m_veg - desis_l12_sfmnn_30m_veg

# Compare reconstructed and retrieved total SIF maps
plot_spatial_comparison(desis_l12_sfmnn_mix_30m_veg, desis_l12_sfmnn_30m_veg,
                        title1 = "", title2 = "", interval = 0.3)

# Plot the spatial distribution of prediction differences
plot_sif_hist(sif_diff_veg, bins = 50,
              fill_palette = c("#0571b0", "lightgrey", "#f4a582", "#ca0020"),
              breaks = c(-0.38, -0.15, 0.15, 0.45, 0.62))

# Assess consistency between reconstructed and retrieved mixed-pixel SIF
plot_purity_sif_consistency(
  desis_l12_sfmnn_mix_30m, desis_l12_sfmnn_30m,
  class_map_l12_30m, pure_class_map_l12_30m,
  metrics_pos = c(0.05, 0.85), plot_title = ""
  # "Predicted vs Retrieved DESIS SIF at 30m Resolution"
)

# Visualize the relative SIF contribution of the three vegetation classes
plot_ternary(desis_l12_sfmnn_prod_30m)

# Plot the ternary colour scale
plot_ternary(type = "mesh", res = 10, blendPow = 300)

# Plot class-specific average SIF maps
plot_avg_sif_maps(
  desis_l12_sfmnn_prod_30m,
  class_colors = c("#984EA3", "#a6d854", "#1C7255"),
  gamma = 0.5, ncol = 3,
  titles = c("Crops", "Mixed Vegetation", "Trees"),
  units = NULL, axis_font_size = 1.2, title_font_size = 1.5
)

# Evaluate reconstructed and retrieved DESIS SIF by pixel purity
plot_purity_sif_consistency_classwise(
  predicted_raster = desis_l12_sfmnn_mix_30m,
  validation_raster = desis_l12_sfmnn_30m,
  class_map = class_map_l12_30m,
  pure_class_map = pure_class_map_l12_30m,
  plot_title = "", # "Predicted vs Retrieved DESIS SIF"
  panel_type = "scatter",
  pixel_type_to_plot = c("Pure", "Mixed", "Combined"),
  show_panel_titles = FALSE
)

# Evaluate class-specific SIF predictions across spatial resolutions
plot_class_sif_consistency_multires(
  predicted_30m = hyplant_l12_sfm_pred_30m,
  validation_30m = hyplant_l12_sfm_30m,
  class_map_30m = class_map_l12_30m,
  predicted_60m = hyplant_l12_sfm_pred_60m,
  validation_60m = hyplant_l12_sfm_60m,
  class_map_60m = class_map_l12_60m,
  predicted_300m = hyplant_l12_sfm_pred_300m,
  validation_300m = hyplant_l12_sfm_300m,
  class_map_300m = class_map_l12_300m,
  plot_title = "",
  metrics_pos = list(
    "Crops" = c(0, 1.1),
    "Mixed Vegetation" = c(0, 1.1),
    "Trees" = c(0, 1.23)
  )
)

# Evaluate predicted vs retrieved class-specific SIF contributions
plot_classwise_spatial_validation(
  predicted_30m = hyplant_l12_sfm_pred_30m,
  validation_30m =  class_sif_cont_l12_sfm_30m,
  predicted_60m = hyplant_l12_sfm_pred_60m,
  validation_60m =  class_sif_cont_l12_sfm_60m,
  predicted_300m = hyplant_l12_sfm_pred_300m,
  validation_300m =  class_sif_cont_l12_sfm_300m,
  plot_title = ""
)

# Compare RF and linear regression unmixing approaches
compare_unmixing_approaches(
  class_sif_cont_l12_sfm_30m,
  class_frac_l12_30m,
  hyplant_l12_sfm_30m,
  hyplant_l12_sfm_pred_30m,
  hyplant_l12_sfm_lm_pred_30m
)

# Evaluate 2020 HyPlant class-specific SIF predictions
evaluate_class_sif_cont(
  hyplant_l3_sfm_pred_30m,
  class_sif_cont_l3_sfm_30m,
  plot_title = ""
)
