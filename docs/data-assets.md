# Data & Assets

This page summarizes the EII data products and their Google Earth Engine assets.
For programmatic access, use `eii.client.ASSETS` and `eii.client.get_asset_path`.

## Core EII assets

| Asset key | Description | Bands | Resolution | Asset path |
|-----------|-------------|-------|------------|------------|
| `eii` | Combined Ecosystem Integrity Index | `eii` | 300 m | `projects/landler-open-data/assets/eii/global/eii_global_v1` |
| `components` | Multiband component stack | `functional_integrity`, `structural_integrity`, `compositional_integrity` | 300 m | `projects/landler-open-data/assets/eii/global/eii_global_v1` |
| `structural_core_area` | Structural integrity core area | `structural_integrity` | 300 m | `projects/landler-open-data/assets/eii/products/v1/structural_integrity/core_area` |
| `npp_predictions` | NPP model predictions | `potential_npp`, `actual_npp`, `relative_npp`, `npp_difference` | 300 m | `projects/landler-open-data/assets/eii/predictions/npp` |
| `npp_model` | Trained NPP model | n/a | n/a | `projects/landler-open-data/assets/eii/models/potential_npp_classifier` |

Use the asset registry for exact paths:

```python
from eii.client import ASSETS, get_asset_path

print(get_asset_path("eii"))
print(ASSETS["components"]["bands"])
```

## Source datasets (reference)

| Dataset | Purpose | Source |
|---------|---------|--------|
| CHELSA NPP | Climatic NPP proxy | `projects/landler-open-data/assets/datasets/chelsa/npp/chelsa_npp_1981_2010_v2-1` |
| CLMS NPP | Observed NPP | `projects/landler-open-data/assets/datasets/clms/npp/annual` |
| SoilGrids (sand) | Soil predictors | `projects/landler-open-data/assets/datasets/soilgrids/sand/sand_15-30cm_mean_gapfilled` |
| SoilGrids (clay) | Soil predictors | `projects/landler-open-data/assets/datasets/soilgrids/clay/clay_15-30cm_mean_gapfilled` |
| SoilGrids (pH) | Soil predictors | `projects/landler-open-data/assets/datasets/soilgrids/phh2o/phh2o_15-30cm_mean_gapfilled` |
| HMI masks | Natural area masks | `projects/landler-open-data/assets/datasets/natural_lands/hmi_masks/v1` |
| Ecoregions | Sampling stratification | `RESOLVE/ECOREGIONS/2017` |
| WDPA | Protected areas | `WCMC/WDPA/current/polygons` |

## Versioning and dates

The `ASSETS` registry includes a reference date for each product. Use these to
anchor analysis to a specific snapshot. See `eii.client.assets.REFERENCE_DATES`.

If you need full provenance and processing details, refer to the methodology pages.
