# =============================================================================
# 07. INSPECTION
# =============================================================================

# Run the required workflow stages for inspection
source("Setup.R")
source("01_Load_Preprocess.R")
source("02_Feature_Construction.R")
source("03_Model_Development.R")

# Inspect pure and mixed pixels at the selected purity threshold (0.85)
plot_pure_mixed(class_map_l12_30m, pure_class_map_l12_30m)

# Inspect the DESIS predictor variables used for model training
plot_maps_viridis(
  raster_data = desis_l12_vi_pure_60m,
  # class_colors = c("#984EA3", "#a6d854", "#1C7255"),
  gamma = 0.5, ncol = 2,
  titles = c("Crops", "Mixed Vegetation", "Trees"),
  units = NULL, axis_font_size = 1.2, title_font_size = 1.5
)

# Inspect class fractions in pure pixels
plot_maps_viridis(
  raster_data = class_frac_l12_pure_60m,
  # class_colors = c("#984EA3", "#a6d854", "#1C7255"),
  gamma = 0.5, ncol = 3,
  titles = c("Crops", "Mixed Vegetation", "Trees"),
  units = NULL, axis_font_size = 1.2, title_font_size = 1.5
)

# Inspect retrieved pure-pixel DESIS SIF (SIFret)
plot(desis_l12_sfmnn_pure_60m, axes = FALSE, box = TRUE)
