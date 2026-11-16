# CLMS NPP Pre-processing Pipeline

Converts CLMS Dry Matter Productivity (DMP) 10-daily dekadal data to annual Net Primary Production (NPP) sums.

## Configuration

Edit `config.cfg` to set input/output paths, year range, and processing parameters.

## Source product

The pipeline uses the consolidated CLMS DMP 300 m V2.0.2 reprocessing (RT6) from
`dmp_global_300m_10daily_v2`, one reprocessing across the PROBA-V (Nov 2013 to 2018)
and Sentinel-3 OLCI (2019 onward) periods; `PRODUCT_GLOB` selects it. The sensor
change at 2019 shifts annual sums by up to 25 % per biome, so EII v2 uses the OLCI
years only (`eii.compute.settings.OBSERVED_NPP_FIRST_YEAR`).

The default download and annual sums cover 2019–2025 (`DOWNLOAD_START_YEAR`,
`START_YEAR` to `END_YEAR`), the years EII v2 uses: training reads the OLCI years and
the first reference year, 2021, needs 2019. The PROBA-V years are optional: set the
start years to 2014 for annual sums from 2014, or `DOWNLOAD_START_YEAR=2013` for the
Nov–Dec 2013 dekads (seasonal analyses only, no annual sum). 2013–2025 is about
820 GB of raw data.

## Workflow

| Step | Script | Description |
|------|--------|-------------|
| 0 | `00_get_cdse_clms_files.sh` | Downloads raw DMP and QFLAG files from CDSE S3 |
| 1 | `01_calculate_annual_npp.py` | Applies quality mask, interpolates to daily, sums annual NPP |

## Prerequisites

- **Step 0**: AWS CLI and CDSE S3 keys ([setup guide](https://documentation.dataspace.copernicus.eu/APIs/S3.html))
  in `~/.aws/credentials` under the profile named by `AWS_PROFILE`:
  `aws_access_key_id` and `aws_secret_access_key`. The script generates its own AWS
  config with the CDSE endpoint. CDSE allows 4 concurrent S3 connections per user
  (HTTP 429 above that), so `S3_MAX_CONCURRENT_REQUESTS=4`; failed listings and
  syncs are retried with backoff. Reruns resume, and a lock prevents parallel runs.
- **Step 1**: Python with `numpy`, `rasterio`, `numba`

## Aggregation rules

- **Masking:** negative DMP values (-1 missing, -2 water) and QFLAG bits 1 (FAPAR
  invalid) and 2 (meteo invalid) are masked; QFLAG 255 (missing) includes both
  bits. The informational bits (8 EBF, 32 low-FAPAR correction, 64 and 128 FAPAR
  gap filling) are kept. Bit meanings are identical in OLCI V1.1.x and V2.0.2;
  PROBA-V V1.0.1 lacks bit 32.
- **Timing:** CLMS names each dekad by its last day. Values are mean daily rates
  for the dekad and are anchored at the dekad's midpoint before interpolation.
- **Completeness:** each year needs all 36 dekads of a single product version.
- **Grid:** outputs are written on the exact `clms_npp/v1_1` transform. V2.0.2
  differs from it only by float rounding (< 1e-9 degrees over the raster).

Run years in parallel with separate processes, for example
`python 01_calculate_annual_npp.py --start-year 2020 --end-year 2020`.

## Memory and interrupted runs

Aggregation reads one 512 × 512 spatial window across the year's observations,
computes its statistics, and immediately writes both output bands. It uses one
process, a 64 MiB GDAL cache, and no Dask workers or full-raster arrays.
`CHUNK_SIZE` must be a multiple of 16 between 16 and 512. Smaller windows reduce
memory use. Runtime depends on storage and CPU; processing is sequential.

During processing, output is named `annual_npp_YEAR.partial.tif`. It is renamed
to `annual_npp_YEAR.tif` only after a successful write and close. Failed runs exit
with an error; a rerun writes the partial file anew.

## Output

- Annual NPP as two-band INT32 GeoTIFFs: band 1 annual sum (gC m-2 yr-1), band 2
  intra-annual SD of daily values (gC m-2 day-1). Stored values are physical
  values times 10,000 (`scale_factor=0.0001`, nodata -1), as in `clms_npp/v1_1`.
