# SoilGrids Gap-Filling

Fills missing pixels in [SoilGrids](https://gee-community-catalog.org/projects/isric/) data using Gaussian-weighted spatial interpolation.

## Scripts

| Script | Description |
|--------|-------------|
| `interpolate_soilgrids.py` | GEE Python script to gap-fill selected variables and export as assets |
| `explore_gapfilling.js` | GEE Code Editor script to compare original vs gap-filled layers |

## Configuration

Edit `config.yaml` to set GEE project, output path, and variables:

```yaml
GEE_PROJECT="landler-open-data"
OUTPUT_ASSET_PREFIX="projects/landler-open-data/assets/datasets/soilgrids"
SEARCH_DISTANCE_M=5000

# Variables: variable_name,depth (one per line)
# Available: bdod, cec, cfvo, clay, nitrogen, phh2o, sand, silt, soc, ocd, ocs
# Depths: 0-5cm, 5-15cm, 15-30cm, 30-60cm, 60-100cm, 100-200cm
VARIABLES="
sand,15-30cm
clay,15-30cm
"
```

## Method

- Gaussian kernel convolution (default 5km search distance, σ = radius/3)
- Only NA pixels are filled; valid pixels unchanged
- Kernel normalized to account for missing neighbors

## Output

Assets saved to: `{OUTPUT_ASSET_PREFIX}/{variable}/{variable}_{depth}_mean_gapfilled`
