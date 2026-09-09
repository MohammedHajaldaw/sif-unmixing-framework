# Improving Satellite-Based Vegetation Monitoring: A Novel Machine Learning-Based Framework for Unmixing Sun-Induced Fluorescence

This repository contains the R code developed for the machine learning-based framework presented in this study to estimate vegetation-specific Sun-Induced Fluorescence (SIF) from mixed satellite pixels.

The framework uses vegetation fraction maps, retrieved SIF, and vegetation indices to separate the SIF signal of different vegetation types within a mixed pixel. It generates vegetation-specific SIF emission maps for three vegetation classes: **Crops, Mixed Vegetation, and Trees**.

## Repository structure

```text
.
├── Main.R
├── Setup.R
├── 01_Load_Preprocess.R
├── 02_Feature_Construction.R
├── 03_Model_Development.R
├── 04_Prediction_Products.R
├── 05_Evaluation.R
├── 06_Sensitivity_Analysis.R
├── 07_Inspection.R
│
├── Functions/
│   ├── 01_Preprocessing/
│   ├── 02_Classification/
│   ├── 03_Feature_Construction/
│   ├── 04_SIF_Unmixing/
│   ├── 05_Evaluation_Plotting/
│   └── 06_Sensitivity_Analysis/
│
├── LICENSE
├── README.md
├── CITATION.cff
└── .gitignore
```

`Main.R` runs the complete workflow, while `Setup.R` loads the required packages and functions.

## Requirements

The workflow is implemented in **R**. Required packages are listed in `Setup.R`.

The analysis requires access to the DESIS and HyPlant datasets used in the study.

## Running the workflow

From the project directory:

```r
source("Main.R")
```

Run the workflow from a clean R session.

## Data availability

The datasets used in this study are not included because access is restricted by the data providers.

HyPlant and related SIF datasets can be requested from **Forschungszentrum Jülich (FZJ)**. DESIS imagery can be requested from the **German Aerospace Center (DLR)**.

## Citation

### Publication

> Hajaldaw, M. (2026). *Improving Satellite-Based Vegetation Monitoring: A Novel Machine Learning-Based Framework for Unmixing Sun-Induced Fluorescence*. Manuscript under review at the *International Journal of Applied Earth Observation and Geoinformation*.

### Source Code

> Hajaldaw, M. (2026). *Improving Satellite-Based Vegetation Monitoring: A Novel Machine Learning-Based Framework for Unmixing Sun-Induced Fluorescence*. Zenodo. DOI: To be added once the Zenodo DOI is available.

## License

This project is licensed under the **GNU General Public License version 3 or later (GPLv3+)**. See [`LICENSE`](LICENSE).

This license applies to the original code in this repository. Third-party software, datasets, and external materials remain subject to their respective licenses and terms.

## AI-assisted development

ChatGPT (OpenAI) was used to assist with coding, debugging, organization, review, and documentation. The corresponding author remains responsible for the scientific methods, analysis, interpretation, and final code.
