# ==============================================================================
# Improving Satellite-Based Vegetation Monitoring:
# A Novel Machine Learning-Based Framework for Unmixing Sun-Induced Fluorescence
#
# Corresponding Author: Mohammed Hajaldaw
# Date: 09 September 2026
# License: GPLv3+
# ==============================================================================


# =============================================================================
# PROJECT SETUP
# =============================================================================

# Load required packages, set the working directory, and source all functions
source("Setup.R")


# =============================================================================
# ANALYSIS WORKFLOW
# =============================================================================

source("01_Load_Preprocess.R")
source("02_Feature_Construction.R")
source("03_Model_Development.R")
source("04_Prediction_Products.R")
source("05_Evaluation.R")
source("06_Sensitivity_Analysis.R")
source("07_Inspection.R")


# =============================================================================
# AI USE
# =============================================================================

# ChatGPT (OpenAI) was used to assist in coding, debugging, organization,
# review, and documentation. The corresponding author remains responsible 
# for the scientific methods, analysis, interpretation, and final code.
