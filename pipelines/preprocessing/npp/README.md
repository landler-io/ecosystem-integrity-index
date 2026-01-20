# CLMS NPP Pre-processing Pipeline

Converts CLMS Dry Matter Productivity (DMP) 10-daily dekadal data to annual Net Primary Production (NPP) sums.

## Configuration

Edit `config.yaml` to set input/output paths, year range, and processing parameters.

## Workflow

| Step | Script | Description |
|------|--------|-------------|
| 0 | `00_get_cdse_clms_files.sh` | Downloads raw DMP and QFLAG files from CDSE S3 |
| 1 | `01_calculate_annual_npp.py` | Applies quality mask, interpolates to daily, sums annual NPP |
| 2 | `02_upload_to_gee.sh` | Uploads annual GeoTIFFs to GCS and ingests into Earth Engine |

## Prerequisites

- **Step 0**: CDSE account with S3 credentials ([setup guide](https://documentation.dataspace.copernicus.eu/APIs/S3.html))
- **Step 1**: Python with `xarray`, `rioxarray`, `dask`, `numba`
- **Step 2**: `gcloud` and `earthengine` CLI authenticated

## Recommended VM (Step 1)

Processing is I/O bound. Example GCP configuration:

| Component | Spec | Notes |
|-----------|------|-------|
| Machine | n2-highmem-16 | 16 vCPUs, 128GB RAM |
| Storage | 1TB SSD (pd-ssd) | ~900GB raw data + output |
| Config | 5 Dask workers, 30GB each | Adjust in `config.yaml` |

Runtime: ~30-45 min per year on this configuration.

## Output

- Annual NPP sums as INT32 GeoTIFFs (scaled by 100, units: g C / m² / year)
- Long-term average across configured year range
