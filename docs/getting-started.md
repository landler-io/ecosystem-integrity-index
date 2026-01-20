# Getting Started

Install the package, authenticate Earth Engine, and run a few common workflows.
Consider this your consolidated primer, then start excuting our [example notebooks]("https://github.com/landler-io/ecosystem-integrity-index/tree/main/notebooks").

## 1) Prerequisites

- Python 3.10 or higher
- Google Earth Engine account ([sign up here](https://earthengine.google.com/signup/))

## 2) Install

```bash
git clone https://github.com/landler-io/ecosystem-integrity-index.git
cd ecosystem-integrity-index
pip install -e ".[client]"
```

Optional dependency sets:

- `.[client]` for raster downloads and local data handling
- `.[compute]` to compute EII components on the fly
- `.[training]` for model training utilities
- `.[docs]` to build the documentation locally

## 3) Authenticate Earth Engine

Authenticate once in a browser, then initialize with your project.

```python
import ee

ee.Authenticate()
ee.Initialize(project="your-gee-project-id")
```

## 4) Define an area of interest

```python
import ee

# Option 1: Rectangle from bounds
polygon = ee.Geometry.Rectangle([-60, -10, -55, -5])

# Option 2: Polygon from coordinates
polygon = ee.Geometry.Polygon([
    [[-60, -10], [-55, -10], [-55, -5], [-60, -5], [-60, -10]]
])

# Option 3: From a FeatureCollection
countries = ee.FeatureCollection("FAO/GAUL/2015/level0")
brazil = countries.filter(ee.Filter.eq("ADM0_NAME", "Brazil")).geometry()
```

## 5) Get summary statistics for an AOI

```python
from eii.client import get_stats

stats = get_stats(polygon, stats=["mean", "min", "max"])
print(stats)
```

Example output:

```python
{
    "geometry_type": "Polygon",
    "values": {
        "eii": {"mean": 0.72, "min": 0.1, "max": 1.0},
        "functional_integrity": {"mean": 0.8, "min": 0.2, "max": 1.0},
        "structural_integrity": {"mean": 0.6, "min": 0.1, "max": 1.0},
        "compositional_integrity": {"mean": 0.7, "min": 0.3, "max": 1.0}
    }
}
```

To return EII only:

```python
stats = get_stats(polygon, include_components=False)
```

## 6) Download rasters

```python
from eii.client import get_raster

dataset = get_raster(polygon, include_components=True, output_format="memory")
print(dataset)
```

Use `get_raster` to download EII and component bands, then analyze with xarray or GeoPandas.

## 7) Get component layers

```python
from eii.client import get_layers

layers = get_layers(layers="components")
```

## 8) Regional analysis (multi-region stats)

```python
import ee
from eii.client import get_zonal_stats

countries = ee.FeatureCollection("FAO/GAUL/2015/level0")
south_america = countries.filter(
    ee.Filter.inList("ADM0_NAME", [
        "Brazil", "Argentina", "Colombia", "Peru", "Chile"
    ])
)

stats = get_zonal_stats(south_america)
print(stats.first().getInfo())
```

### Export to CSV

```python
task = ee.batch.Export.table.toDrive(
    collection=stats,
    description="eii_south_america",
    fileFormat="CSV"
)
task.start()
```

### Visualization in Earth Engine

```python
import geemap

Map = geemap.Map()
Map.addLayer(stats, {"column": "eii"}, "EII")
Map
```

## Next Steps

- Read about the [Methodology](methodology/index.md)
- Check the [API Reference](api/index.md) for all available functions
