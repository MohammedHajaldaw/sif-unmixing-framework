# =============================================================================
# SETUP
# =============================================================================

# Required packages
required_packages <- c(
  "readxl", "sp", "terra", "sf", "xml2", "zoo", "spatialEco", "tmap",
  "pracma", "ggplot2", "sensitivity", "dplyr", "lhs", "rhdf5", "ranger",
  "caret", "patchwork", "gridExtra", "gstat", "scales", "tidyr", "GGally",
  "rlang", "grid", "viridisLite", "reshape2", "ggtern", "RColorBrewer"
)

# Install and load packages
install_and_load <- function(package_name) {
  tryCatch({
    if (!package_name %in% installed.packages()) {
      install.packages(package_name)
    }
    library(package_name, character.only = TRUE)
    cat(package_name, "loaded successfully\n")
  }, error = function(e) {
    cat("Warning: Could not load", package_name, "-", e$message, "\n")
  })
}

invisible(lapply(required_packages, install_and_load))

# Set working directory
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  current_script_path <- rstudioapi::getActiveDocumentContext()$path
  setwd(dirname(current_script_path))
} else {
  stop("rstudioapi package is not available.")
}

# Source all project functions
function_files <- sort(list.files(
  "./Functions",
  pattern = "\\.R$",
  full.names = TRUE,
  recursive = TRUE
))

invisible(lapply(function_files, source))