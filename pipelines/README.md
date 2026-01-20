# Pipelines

This directory contains data processing and model training pipelines.

These are research workflows for reproducing the EII data products,
distinct from the user-facing Python package in `src/eii/`.

## Structure

```
pipelines/
├── preprocessing/           # Data preparation scripts
│   ├── ecoregions/          # Ecoregion rasterization
│   ├── natural_lands/       # Natural/pristine areas identification
│   ├── npp/                 # NPP data processing
│   ├── protected_areas/     # Protected areas rasterization
│   ├── soilgrids/           # Soil data preparation
│   └── utils/               # Shared preprocessing utilities
├── modeling/                # Model training and evaluation
│   ├── functional_integrity/# NPP model training and inference
│   ├── structural_integrity/# Structural integrity modeling
│   └── eii_calculation/     # EII computation scripts
├── analysis/                # Model diagnostics and figures
└── export/                  # Asset export scripts
```

## Running Pipelines

Most pipelines are designed to run in Google Earth Engine (JavaScript)
or as Python scripts with Earth Engine Python API.

See individual subdirectory READMEs for specific instructions.
