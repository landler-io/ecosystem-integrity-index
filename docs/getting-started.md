# Getting Started

Install the package, authenticate Earth Engine and run the common workflows. The
[example notebooks](https://github.com/landler-io/ecosystem-integrity-index/tree/main/notebooks)
build on this page.

## 1. Prerequisites

- Python 3.10 or higher
- A Google Earth Engine account ([sign up](https://earthengine.google.com/signup/)) with a
  Cloud project

## 2. Install

The package is installed from GitHub:

```bash
pip install "ecosystem-integrity-index[client] @ git+https://github.com/landler-io/ecosystem-integrity-index"
```

Optional dependency sets:

- `[client]`: raster downloads and local data handling (GeoPandas, xarray, rasterio)
- `[compute]`: the EII calculation on your own numpy arrays
- `[training]`: sampling design and validation of the potential-NPP model
- `[docs]`: build this documentation locally

For development, clone the repository and install it in editable mode:

```bash
git clone https://github.com/landler-io/ecosystem-integrity-index.git
cd ecosystem-integrity-index
pip install -e ".[all]"
```

## 3. Authenticate

```python
import ee

ee.Authenticate()
ee.Initialize(project="your-gee-project-id")
```

The client reads EII v2 by default. v2 is annual: each map refers to a reference year between
2021 and 2025, and calls without `reference_year` use the latest year (2025). Pass the year
explicitly so that results stay the same when a new year is published.

## 4. Define an area of interest

The client accepts `ee.Geometry`, `ee.Feature`, `ee.FeatureCollection`, bounding-box tuples,
GeoPandas objects and shapely geometries.

```python
polygon = ee.Geometry.Rectangle([-60, -10, -55, -5])

bbox = (11.0, 47.5, 11.5, 48.0)  # (min_lon, min_lat, max_lon, max_lat)

countries = ee.FeatureCollection("FAO/GAUL/2015/level0")
brazil = countries.filter(ee.Filter.eq("ADM0_NAME", "Brazil")).geometry()
```

## 5. Summary statistics for an area

```python
from eii.client import get_stats

stats = get_stats(polygon, reference_year=2024, stats=["mean", "min", "max"])
```

The result has this shape:

```python
{
    "geometry_type": "Polygon",
    "values": {
        "eii": {"mean": ..., "min": ..., "max": ...},
        "functional_integrity": {"mean": ..., "min": ..., "max": ...},
        "structural_integrity": {"mean": ..., "min": ..., "max": ...},
        "compositional_integrity": {"mean": ..., "min": ..., "max": ...},
        "wetland_flag": {"share": ...},
    },
}
```

- `stats` takes any of `"mean"`, `"median"`, `"min"`, `"max"`, `"std"`; `percentiles` takes
  integers from 0 to 100, for example `[10, 90]`.
- `include_components=False` returns EII only.
- A point geometry returns the pixel values instead of statistics, and the wetland flag code.
- `wetland_flag.share` is the share of valid pixels where functional integrity is less reliable
  (see [Functional Integrity](methodology/functional_integrity.md#wetland-flag)).
  `exclude_wetlands=True` drops those pixels before the statistics are computed.

## 6. Change between years

```python
change = {
    year: get_stats(polygon, reference_year=year, include_components=False)["values"]["eii"]["mean"]
    for year in range(2021, 2026)
}
```

Structural integrity uses the same human-modification map for every year, so year-to-year change
comes from observed productivity, vegetation greenness and biodiversity intactness. Read
short-term changes with the [limitations](limitations.md#annual-change) in mind.

## 7. Download rasters

```python
from eii.client import get_raster

dataset = get_raster(polygon, reference_year=2024, include_components=True)
```

`get_raster` returns an `xarray.Dataset` (or writes a GeoTIFF with `output_format="geotiff"`) with
EII, the three components, `limiting_component` and `wetland_flag`. Large areas are split into
tiles automatically.

## 8. Earth Engine layers

```python
from eii.client import get_layers

layers = get_layers("all", reference_year=2024)
eii = layers["eii"]
limiting = layers["limiting_component"]  # 1 functional, 2 structural, 3 compositional

parts = get_layers("fi_change", reference_year=2024)
parts["fi_state"], parts["fi_change"], parts["fi_change_sensors"]
```

`layers="fi_change"` returns the two parts of functional integrity, FI_state and FI_change, and
the number of satellite sensors behind FI_change (0 where functional integrity equals FI_state).

## 9. Statistics for many regions

```python
from eii.client import get_zonal_stats

south_america = countries.filter(
    ee.Filter.inList("ADM0_NAME", ["Brazil", "Argentina", "Colombia", "Peru", "Chile"])
)
gdf = get_zonal_stats(
    south_america, reference_year=2024, stats=["mean", "median"], keep_columns=["ADM0_NAME"]
)
gdf[("eii", "mean")]
```

The result is a GeoDataFrame with one row per feature and MultiIndex columns such as
`("eii", "mean")`, `("functional_integrity", "median")` and `("wetland_flag", "share")`. The input
can be an `ee.FeatureCollection`, a GeoDataFrame or shapely geometries.

## 10. Version 1

EII v1, a single static map, stays available. Select it per call or per session:

```python
from eii.client import get_stats, set_version, use_version

get_stats(polygon, version="v1")

with use_version("v1"):
    v1_stats = get_stats(polygon)

set_version("v1")  # every later call in this session
```

v1 takes no `reference_year` and has no `limiting_component`, `wetland_flag` or FI parts. Its
methodology is documented in the
[v0.1.1 release](https://github.com/landler-io/ecosystem-integrity-index/tree/v0.1.1); the
differences are listed in [Changes from v1](v1-to-v2.md).

## Next steps

- [Methodology](methodology/index.md)
- [Data & Assets](data-assets.md)
- [API Reference](api/index.md)
