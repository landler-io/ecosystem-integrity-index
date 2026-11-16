# eii.client

Retrieval of the published EII layers through Google Earth Engine. The client reads EII v2 by
default; v2 is annual and `reference_year` (2021-2025) defaults to the latest year. Pass it
explicitly for reproducible results. Select v1 with `version="v1"` per call, or with
`set_version("v1")` or `use_version("v1")` per session.

## Data Retrieval

::: eii.client.retrieve
    options:
      members:
        - get_stats
        - get_raster
        - get_layers

## Wetland flag

The categorical `wetland_flag` band has the codes: 0 none, 1 wetland, 2 peatland,
3 mangrove (GLWD v2; see [Functional Integrity](../methodology/functional_integrity.md#wetland-flag)). It marks
pixels where the potential-NPP reference is known to run high, so functional
integrity (and EII) tends to be too low there. The flag never changes any score.

- `get_layers` returns it with `layers="all"` or `layers="wetland_flag"`.
- `get_raster(..., include_components=True)` appends it as the last band.
- `get_stats` and `get_zonal_stats` report flagged pixels and keep them in the
  statistics by default. `get_stats` adds `"wetland_flag"`: the code at a point, or
  `{"share": ...}` for an area (share of valid EII pixels that are flagged).
  `get_zonal_stats` adds the column `("wetland_flag", "share")`.
- `exclude_wetlands=True` masks flagged pixels before the reduction; the share is
  still computed before exclusion.

```python
from eii.client import get_stats

stats = get_stats(polygon, reference_year=2024, exclude_wetlands=True)
stats["values"]["wetland_flag"]  # {"share": ...}, computed before the exclusion
```

## Functional integrity parts

`get_layers("fi_change", reference_year=...)` returns `fi_state` (observed against potential
NPP), `fi_change` (EVI change against the pixel's 2009-2018 baseline, masked where no sensor
has a value) and `fi_change_sensors` (uint8, 0-2; 0 where functional integrity equals
`fi_state`). The `functional` layer is `fi_state` × `fi_change`. See
[Functional Integrity](../methodology/functional_integrity.md#change-term).

## Dataset version and asset paths

::: eii.client.assets
    options:
      show_source: false
