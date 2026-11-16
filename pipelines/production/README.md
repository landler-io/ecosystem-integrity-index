# Production

The global EII v2.0 run and its publication. Inputs come from `pipelines/preprocessing/`
(annual NPP, predictor stack, GHM, GLWD, ecoregions) and `pipelines/training/` (potential NPP,
natural CV, the stored constants in `src/eii/compute/`).

## 1. FI change term inputs

```bash
python pipelines/production/fi_change_term.py tiles fi_change
python pipelines/production/fi_change_term.py sensor fi_change --sensor T \
    --tiles @fi_change/tiles_T.txt --ecoregions ECOREGIONS.tif
python pipelines/production/fi_change_term.py sensor fi_change --sensor V \
    --tiles @fi_change/tiles_V.txt --ecoregions ECOREGIONS.tif
python pipelines/production/fi_change_term.py vrt fi_change
```

reads MODIS Terra MOD13Q1 061 (Planetary Computer, LP DAAC for missing composites) and VIIRS
Suomi NPP VNP13A1 002 (LP DAAC) per sinusoidal tile, and writes the drift-corrected EVI ratio r
per reference year on the 300 m EII grid (`r_T.vrt`, `r_V.vrt`) and the complete sensor years
(`completeness.json`). Readers and credentials (`PC_SDK_SUBSCRIPTION_KEY` or
`PC_SETTINGS_FILE`, `EDL_TOKEN_FILE` or `~/.netrc`) are in `evi_sources.py`. The MOD13Q1 HDF4
granules are converted with `gdal_translate`, which needs GDAL's HDF4 driver; set
`GDAL_TRANSLATE` to such a build if the one on `PATH` lacks it.

Reference years 2021-2025 are served by the TV series (Terra and Suomi NPP), so the T and V
sensors are all that v2.0 needs. The TNc and Nc series, with NOAA-20 VIIRS (VJ113A1) chained to
Suomi NPP, take over for reference years after the Suomi NPP and MODIS missions end. They are not
wired up yet: `--sensor VN` stores the NOAA-20 chain factor and baseline, but no stage writes
`r_Nc.vrt`, and `compute_eii_v2.py` stops with an error for a year served by TNc or Nc.

## 2. Global EII

```bash
python pipelines/production/compute_eii_v2.py run potential_npp_v2 eii_v2 \
    --clms-dir npp_annual_sums_v2 --natural-cv potential_npp_v2/natural_cv.vrt \
    --hmi eii_inputs/HMv20240801_2022s_AA_300.tif \
    --glwd-class eii_inputs/GLWD_v2_0_main_class.tif \
    --glwd-cover eii_inputs/GLWD_v2_0_area_pct.tif \
    --ecoregions predictor_stack/ecoregions/resolve_ecoregions_rasterized/resolve_ecoregions_rasterized_0000.tif \
    --fi-change fi_change --predictor-stack predictor_stack/stack.vrt
```

With these inputs and the default settings (`eii.compute.settings`) the run reproduces the
published v2.0 tiles: functional integrity by deviation tolerance with natural CV and per-biome
recentring, times the FI change term; structural integrity from the logistic habitat quality of
HMI (midpoint 0.4, steepness 10); compositional integrity from BII v1.1; EII by
`min_fuzzy_logic`. It first computes the recentring factors per reference year
(`<out>/recentring/`), then writes `<out>/<year>/<tile>.tif` for 2021-2025 with the bands
`eii`, `functional_integrity`, `structural_integrity`, `compositional_integrity`,
`limiting_component`, `wetland_flag`, `fi_state`, `fi_change` and `fi_change_sensors`, plus
`qa.csv` and `run_manifest.json`. GHM 2022 is downloaded to `--hmi` if missing; GLWD v2 is
figshare 10.6084/m9.figshare.28519994. BII is read from the requester-pays bucket through
GDAL: set `EII_BILLING_PROJECT` (or `--billing-project`) to a Google Cloud project you can bill;
credentials come from `GOOGLE_APPLICATION_CREDENTIALS`, the gcloud application-default login
or the Compute Engine metadata server. `compute_eii_v2.py factors` runs the recentring pass alone.
Existing tiles are skipped, so the run can be resumed.

## 3. Publication

```bash
python pipelines/production/publish_v2.py upload potential_npp potential_npp_v2
python pipelines/production/publish_v2.py upload eii eii_v2
python pipelines/production/publish_v2.py register potential_npp --public
python pipelines/production/publish_v2.py register eii --public
```

These steps write to the Landler bucket and Earth Engine project and record how the published
collections were made (locations in `pipelines/README.md`). `upload` writes ZSTD COGs to
`gs://landler-opendata-gee-assets/products/eii/v2_0/...` without overwriting and records
checksums and sample values in a publish manifest (`--manifest`); `register` creates the Earth Engine image collection (for `eii`, the one `eii.client` reads:
one image per tile and reference year) and checks it against the manifest.
