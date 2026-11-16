# Preprocessing

Source data for EII v2: annual NPP, the potential-NPP predictor stack and the input datasets in
GCS. Outputs feed `pipelines/training/` and `pipelines/production/`.

## Run order

1. `npp/`: CLMS DMP 300 m V2.0.2 download and annual NPP sums 2019-2025, 2014-2018 optional
   (`00_get_cdse_clms_files.sh`, `01_calculate_annual_npp.py`; see [npp/README.md](npp/README.md)).
   Training reads the sums at the sample points; production reads them per tile.
2. Source datasets in `gs://landler-opendata-gee-assets/datasets/`. These ingests write to the
   Landler bucket and Earth Engine project and record how the published copies were made; a
   reproduction starts from the published datasets (step 3, `sync_predictor_sources.py`):
   - `ingest_chelsa_v2_1.py`: CHELSA V2.1 1981-2010 bio1, bio4, bio12, bio15 and PET as COGs;
   - `ingest_coastal_v2.py build|register`: coast distance (OSM), tidal range (EOT20) and sea
     surface salinity (WOA23) on the 300 m grid;
   - `ingest_bii_v1_1.py run|verify`: Biodiversity Intactness v1.1 2017-2025, the compositional
     component (read from the bucket by `compute_eii_v2.py`).

   Two bucket layers have no ingest script: `chelsa_npp` (CHELSA V2.1 NPP 1981-2010,
   `CHELSA_npp_1981-2010_V.2.1.tif` of CHELSA-BIOCLIM+, EnviDat doi:10.16904/envidat.332) and
   `ecoregions` (RESOLVE Ecoregions 2017, Earth Engine `RESOLVE/ECOREGIONS/2017`, rasterised to
   `ECO_ID` at 5 km in EPSG:4326).
3. Predictor stack:
   - `export_predictor_stack.py run`: Earth Engine export of the MERIT DEM and Geomorpho90m
     terrain bands and MERIT Hydro upstream area into the Landler bucket (published as
     `products/eii/v2_0/predictor_stack/ee_export/`);
   - `sync_predictor_sources.py <root>`: mirror the exports and the bucket layers to `<root>`
     (requester pays: `EII_BILLING_PROJECT` or `--billing-project`); a reproduction starts
     here;
   - `build_predictor_stack.py <root>`: one float32 GeoTIFF per 10° land tile on the CLMS NPP
     grid (1/336°) with the model-card predictors, plus `stack.vrt`. The stack also holds
     `slope`, which the model does not use; predictors are read by band name.

## Shared modules and tables

- `local_sources.py`: readers for the candidate-pool sources (GHM 1990-2022 from Zenodo,
  WorldCover 2021 from AWS, RESOLVE ecoregions from the bucket). The WDPA raster (IUCN Ia,
  Ib, II, 1 km inward buffer, on the GHM grid) is not distributed because the WDPA terms
  forbid redistribution; build it from your own WDPA download and pass it to
  `build_candidate_pool_local.py --wdpa` (selection rule in its docstring). Also used by
  `compute_eii_v2.py` to download GHM 2022, and holds `billing_project()`, the reader of
  `EII_BILLING_PROJECT`.
- `reference/eco_id_to_biome.json`, `reference/eco_id_to_realm.json`: RESOLVE 2017 `ECO_ID` to
  biome number and realm index.
- `utils/config_utils.py`: reader of `npp/config.cfg`.
