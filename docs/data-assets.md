# Data & Assets

EII is computed per 10° tile and reference year
([`pipelines/production`](https://github.com/landler-io/ecosystem-integrity-index/tree/main/pipelines/production)).
The tiles are stored as ZSTD-compressed cloud-optimised GeoTIFFs (COGs) in Google Cloud Storage
and registered as Earth Engine image collections, which read the files in place.

All products share one grid: the CLMS 300 m grid, EPSG:4326, 1/336° pixels, 10° tiles, float32,
nodata -9999.

## Products

### EII

Earth Engine: `projects/landler-open-data/assets/products/eii/v2_0/release/global`

One image per tile and reference year, with the properties `reference_year` (2021-2025) and
`tile` (for example `e010n40`, named by its south-west corner).

| Band | Values |
|---|---|
| `eii` | Ecosystem Integrity Index, 0-1 |
| `functional_integrity` | FI = `fi_state` × `fi_change`, 0-1; equals `fi_state` where `fi_change_sensors` = 0 |
| `structural_integrity` | 0-1 |
| `compositional_integrity` | 0-1 |
| `limiting_component` | 1 functional, 2 structural, 3 compositional; nodata where a component is missing |
| `wetland_flag` | 0 none, 1 wetland, 2 peatland, 3 mangrove |
| `fi_state` | observed against potential NPP, 0-1 |
| `fi_change` | the pixel's EVI over the three-year window against its own baseline, 0-1; nodata where no sensor has a value |
| `fi_change_sensors` | number of satellite sensors behind `fi_change`, 0-2; 0 where FI = FI_state. Stored as float32 in the tile, returned as uint8 by the client |

`eii` is computed from `functional_integrity`. See
[Functional Integrity](methodology/functional_integrity.md#change-term) for the sensors that
serve each reference year.

Every tile carries the full run configuration in its GeoTIFF tags:
functional formula and multipliers, the recentring factors of its reference year, input versions,
the structural quality midpoint, the serving sensor series (`fi_change_series`), `fi_change_k`
and the SHA-256 of the stored change-term constants (`fi_change_constants_sha256`).

### Potential NPP

Earth Engine:
`projects/landler-open-data/assets/products/eii/v2_0/intermediate/functional/potential_npp_olci_gapfill`

One static image per tile, property `tile`.

| Band | Values |
|---|---|
| `potential_npp` | potential NPP, gC m⁻² yr⁻¹ (bias-corrected, before recentring) |
| `lower_68`, `upper_68` | 68 % prediction interval |
| `lower_90`, `upper_90` | 90 % prediction interval |

### Cloud Storage

The COGs behind both collections are in the bucket `landler-opendata-gee-assets`:

| Product | Path |
|---|---|
| EII | `gs://landler-opendata-gee-assets/products/eii/v2_0/release/global/<year>/<tile>.tif` |
| Potential NPP | `gs://landler-opendata-gee-assets/products/eii/v2_0/intermediate/functional/potential_npp_olci_gapfill/<tile>.tif` |

The bucket is Requester Pays: reading the files outside Earth Engine is billed to the reader's
Cloud project, for example `gcloud storage cp --billing-project=<your-project> gs://...` or
GDAL with `GS_USER_PROJECT=<your-project>`.

### Training table

`gs://landler-opendata-gee-assets/products/eii/v2_0/training/npp_training_samples_v2_wp11_localstack_dmpv2_olci_20261104.csv.gz`
(Requester Pays, about 42 MB)

The candidate sites of the [potential NPP model](methodology/potential_npp.md), one row per
site (218,721 rows): coordinates, biome and realm, ESA WorldCover class, HMI, the reference-site
flags (`reference`, `eligibility_source`), the annual CLMS NPP sums 2014-2025 (`sum_<year>`, in
stored units of gC m⁻² yr⁻¹ × 10,000), the evaluation-sample flag (`selected`), block and
cluster identifiers for cross-validation, and the predictor values at each site. The table
holds `slope`, which the model does not use; the model's predictors are listed in the
[model card](https://github.com/landler-io/ecosystem-integrity-index/tree/main/src/eii/training/model_cards).

### Access

The client resolves the collection and mosaics one reference year:

```python
from eii.client import get_asset_path, get_layers

get_asset_path("eii")  # the v2 collection
layers = get_layers("all", reference_year=2024)
parts = get_layers("fi_change", reference_year=2024)  # fi_state, fi_change, fi_change_sensors
```

Directly in Earth Engine:

```python
import ee

eii_2024 = (
    ee.ImageCollection("projects/landler-open-data/assets/products/eii/v2_0/release/global")
    .filter(ee.Filter.eq("reference_year", 2024))
    .mosaic()
)
potential = ee.ImageCollection(
    "projects/landler-open-data/assets/products/eii/v2_0/intermediate/functional/potential_npp_olci_gapfill"
).mosaic()
```

### Version 1

EII v1 remains available as a static, single-image product:

| Asset | Path |
|---|---|
| EII and components | `projects/landler-open-data/assets/eii/global/eii_global_v1` (bands `eii`, `functional_integrity`, `structural_integrity`, `compositional_integrity`) |
| Potential NPP | `projects/landler-open-data/assets/eii/intermediate/functional/predictions/natural_npp_v1` |

Select it per call with `version="v1"`, or per session with `eii.client.set_version("v1")` or
the `eii.client.use_version("v1")` context manager. v1 takes no `reference_year`, and
`exclude_wetlands` and the bands `limiting_component`, `wetland_flag` and the FI parts exist in
v2 only.

## Input datasets

| Input | Use | Dataset | Years |
|---|---|---|---|
| Observed NPP | Functional (state) | CLMS Dry Matter Productivity 300 m V2.0.2 (RT6), 10-daily, summed per year | 2019-2025 |
| Vegetation index | Functional (change) | MODIS Terra MOD13Q1 v061; VIIRS Suomi NPP VNP13A1 v002; VIIRS NOAA-20 VJ113A1 v002 (16-day EVI) | 2001-2025 |
| Climate | Potential NPP; change-term tolerance | CHELSA V2.1 bio1, bio4, bio12, bio15, Penman PET, NPP | 1981-2010 |
| Terrain | Potential NPP | MERIT DEM; Geomorpho90m CTI, TPI, northness, eastness | static |
| Coastal setting | Potential NPP | OpenStreetMap coastline distance; EOT20 tidal range; WOA23 salinity; MERIT Hydro upstream area | static |
| Reference sites | Potential NPP training | GHM 1990-2020 natural-land mask; WDPA strict protection; ESA WorldCover exclusions; GHM 2022 for gap-fill sites | |
| Land cover | Change-term drift groups and controls | ESA WorldCover 2021 | 2021 |
| Biomes | Intervals, recentring, drift groups | RESOLVE Ecoregions 2017 | static |
| Human modification | Structural; recentring reference | Global Human Modification 2022, 300 m (doi:10.5281/zenodo.14502573) | 2022 |
| Biodiversity intactness | Compositional; recentring reference | Biodiversity Intactness v1.1 (Impact Observatory, Vizzuality), 100 m | 2021-2025 |
| Wetlands | Wetland flag | GLWD v2 dominant class and cover (Lehner et al. 2025) | static |

Annual NPP sums are stored as integers equal to gC m⁻² yr⁻¹ × 10,000 (nodata -1). Multiply by
`eii.compute.settings.NPP_STORED_SCALE` (`1e-4`) for physical units.

## Licences

The EII code is licensed under Apache 2.0. The published data products (EII tiles, potential NPP
map, training tables) are licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/);
cite them as The Landbanking Group (2026), see the [README](https://github.com/landler-io/ecosystem-integrity-index#citation).
They are derived from the inputs below, whose own licence and attribution terms apply to those
inputs:

| Input | Licence or terms | Attribution |
|---|---|---|
| CLMS Dry Matter Productivity V2.0.2 | Copernicus data policy (free, full and open) | © European Union, Copernicus Land Monitoring Service |
| MODIS Terra MOD13Q1, VIIRS VNP13A1 and VJ113A1 | NASA Earth science data, no restrictions on use | NASA LP DAAC (Didan 2021; Didan & Barreto 2023) |
| CHELSA V2.1 | CC0 1.0 | Karger et al. (2017) |
| MERIT DEM | CC BY-NC 4.0 or ODbL 1.0 (dual licence) | Yamazaki et al. (2017) |
| MERIT Hydro | CC BY-NC 4.0 or ODbL 1.0 (dual licence) | Yamazaki et al. (2019) |
| Geomorpho90m | CC BY 4.0 | Amatulli et al. (2020) |
| OpenStreetMap coastline (`dist_coast`, coastal mask of `tidal_range` and `sss`) | ODbL 1.0 | © OpenStreetMap contributors |
| EOT20 | CC BY 4.0 | Hart-Davis et al. (2021) |
| World Ocean Atlas 2023 | public domain | NOAA National Centers for Environmental Information |
| Global Human Modification 1990-2020 and 2022 | CC BY 4.0 | Theobald (2024), doi:10.5281/zenodo.14449495 and doi:10.5281/zenodo.14502573 |
| WDPA | WDPA terms of use | UNEP-WCMC and IUCN, Protected Planet |
| ESA WorldCover 2021 | CC BY 4.0 | © ESA WorldCover project, Zanaga et al. (2022) |
| RESOLVE Ecoregions 2017 | CC BY 4.0 | Dinerstein et al. (2017) |
| Biodiversity Intactness v1.1 | CC BY 4.0 | Impact Observatory and Vizzuality |
| GLWD v2 | CC BY 4.0 | Lehner et al. (2025) |

The training table holds raw values of MERIT DEM elevation, MERIT Hydro upstream area,
Geomorpho90m terrain indices and the OpenStreetMap-derived coastline distance at each site.

Validation inputs keep their own licences and are not redistributed; see
[Validation](methodology/validation.md).
