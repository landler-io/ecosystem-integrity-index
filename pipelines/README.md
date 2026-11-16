# Pipelines

The workflows behind the EII v2.0 data products, separate from the Python package in
`src/eii/`. Everything runs offline on local rasters with `eii.compute`; Earth Engine is used to
export a few source layers and to register the published COGs as image collections.

```
pipelines/
├── preprocessing/   # source ingests, annual NPP, predictor stack
├── training/        # everything that fits or calibrates a parameter of the product
├── production/      # global EII run and publication
└── validation/      # independent validation of the published product
```

Scripts take their inputs and outputs as arguments; relative paths in the examples resolve from
the repository root.

## Requirements

- Python packages: `pip install -e ".[pipelines]"`.
- Command line tools: the Google Cloud CLI (`gcloud`), the AWS CLI (CLMS download from the
  Copernicus Data Space), and the GDAL command line tools: `gdalbuildvrt` and `gdal_translate`
  with the HDF4 driver (MODIS granules; see `production/evi_sources.py`).
- The bucket `gs://landler-opendata-gee-assets` is requester pays: set `EII_BILLING_PROJECT` to
  a Google Cloud project you can bill (scripts with a command line also take
  `--billing-project`), and log in with `gcloud auth application-default login` or set
  `GOOGLE_APPLICATION_CREDENTIALS`.
- Credentials for the EVI sources of the FI change term (Planetary Computer, NASA Earthdata):
  `production/evi_sources.py`.

## End-to-end run order

Steps marked *provenance* write to the Landler bucket or Earth Engine project. They record how
the published inputs and products were made; their outputs are published (see below), and a
reproduction starts from the published copies.

1. **Preprocessing** ([preprocessing/README.md](preprocessing/README.md))
   1. `npp/00_get_cdse_clms_files.sh`, `npp/01_calculate_annual_npp.py`: CLMS DMP V2.0.2
      10-daily data to annual NPP sums 2019-2025 (2014-2018 optional).
   2. *Provenance*: `ingest_chelsa_v2_1.py`, `ingest_coastal_v2.py build|register`,
      `ingest_bii_v1_1.py run` (climate, coastal and biodiversity-intactness sources, published
      under `gs://landler-opendata-gee-assets/datasets/`) and `export_predictor_stack.py run`
      (terrain and upstream area, published under `.../products/eii/v2_0/predictor_stack/`).
   3. `sync_predictor_sources.py`, `build_predictor_stack.py`: the potential-NPP predictor
      stack on the CLMS NPP grid, starting from the published sources.
2. **Training** ([training/README.md](training/README.md))
   1. Potential-NPP training table: `build_candidate_pool_local.py`,
      `sample_candidate_pool_ee.py`, `sample_training_v2.py`, `extract_stack_at_points.py`,
      `training_table_localstack.py`, `training_table_dmp_v2.py`.
      The final table is published, so this step can be skipped.
   2. `fit_potential_npp_v2.py`: mean model, bias correction, prediction intervals, natural
      inter-annual CV model.
   3. `predict_potential_npp_v2.py` (and `--natural-cv`): global potential NPP tiles.
   4. `build_si_calibration_sample.py`, `calibrate_structural.py`: the calibration sample and
      the check of the logistic habitat-quality curve of structural integrity against
      WorldCover 2021.
3. **Production** ([production/README.md](production/README.md))
   1. `fi_change_term.py tiles|sensor|vrt`: per-sensor EVI ratios of the FI change term.
   2. `compute_eii_v2.py run`: EII per 10° tile and reference year 2021-2025, with the stored
      FI change constants in `src/eii/compute/fi_change_constants/`.
   3. *Provenance*: `publish_v2.py upload|register`: GCS COGs and Earth Engine collections.
4. **FI change constants** (optional check): `training/build_fi_change_controls.py` builds the
   control table and `training/derive_fi_change_constants.py` re-derives the stored constants (drift factors, tolerance models, k_c) and compares them with
   their checksums. Its input holds FI_state of a production run, which itself uses the
   constants, so production always uses the stored copies.
5. **Validation**: [validation/README.md](validation/README.md); `validation/report_figures.py`
   draws the figures of the methodology report.

## Reproducibility notes

Every step has its code in this repository. Some steps depend on Earth Engine, on data that are
not redistributed, or on records that are not published:

- **Candidate pool** (`training/sample_candidate_pool_ee.py`): interactive Earth Engine requests
  against the Earth Engine catalog (WorldClim 1.4, MERIT DEM), the community catalog (Global
  Aridity Index, Geomorpho90m) and public Landler assets (SoilGrids, CHELSA NPP, annual CLMS NPP
  v1.1). The candidates step needs a WDPA high-protection raster built from your own WDPA
  download (`build_candidate_pool_local.py --wdpa`; WDPA may not be redistributed).
- **Structural calibration sample** (`training/build_si_calibration_sample.py`): reproduces the
  published sample exactly. The land test and the record filter run in Earth Engine and need the
  v1 EII and natural-NPP assets; anthromes come from Anthromes v2 (Ellis et al. 2010).
- **FI change controls** (`training/build_fi_change_controls.py`): offline, from WorldCover,
  GHM 2022, RESOLVE, the MODIS and VIIRS EVI composites (Planetary Computer, NASA LP DAAC), the
  predictor stack and a production run. The exclusion list of earlier evaluation points is not
  published, and the published table's FI_state comes from a run before the potential-NPP gap
  fill (details in [training/README.md](training/README.md)). Production uses the stored
  constants either way.
- **CHELSA NPP and RESOLVE ecoregions** have no ingest script: `chelsa_npp` is CHELSA V2.1 NPP
  1981-2010 (`CHELSA_npp_1981-2010_V.2.1.tif`, CHELSA-BIOCLIM+, EnviDat
  doi:10.16904/envidat.332); the ecoregion raster is RESOLVE Ecoregions 2017 (Dinerstein et al.
  2017, Earth Engine `RESOLVE/ECOREGIONS/2017`) rasterised to `ECO_ID` at 5 km. Both are
  published in the bucket and synced by `sync_predictor_sources.py`.
- **Natural inter-annual CV** (`natural_cv`) is not distributed; `predict_potential_npp_v2.py
  --natural-cv` regenerates it from the fitted model.

## Published data

All `gs://` objects are in the requester-pays bucket `landler-opendata-gee-assets`.

| Data | Location |
| --- | --- |
| Predictor sources | `gs://landler-opendata-gee-assets/datasets/{chelsa_*,chelsa_npp,osm_coast_distance,eot20_tidal_range,woa23_sss,ecoregions}/v1_0/`, `.../products/eii/v2_0/predictor_stack/ee_export/` |
| BII v1.1 | `gs://landler-opendata-gee-assets/datasets/bii/v1_1/` |
| Candidate pool | `gs://landler-opendata-gee-assets/products/eii/v2_0/training/candidate_pool/npp_candidate_pool_v2_local_seed42_20261027.csv.gz` |
| Training table | `gs://landler-opendata-gee-assets/products/eii/v2_0/training/npp_training_samples_v2_wp11_localstack_dmpv2_olci_20261104.csv.gz` |
| FI change controls | `gs://landler-opendata-gee-assets/products/eii/v2_0/training/fi_change_controls.csv.gz` |
| Structural calibration sample | `gs://landler-opendata-gee-assets/products/eii/v2_0/training/si_calibration_sample.csv.gz` |
| Potential NPP | `projects/landler-open-data/assets/products/eii/v2_0/intermediate/functional/potential_npp_olci_gapfill` |
| EII v2 | `projects/landler-open-data/assets/products/eii/v2_0/release/global` |

The model card of the potential-NPP model is `src/eii/training/model_cards/v2.json`.
