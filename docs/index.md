# Ecosystem Integrity Index (EII)

The Ecosystem Integrity Index is a global, annual map of terrestrial ecosystem condition at
300 m. It combines three dimensions of integrity into one score between 0 and 1:

- **Functional integrity**: is the land as productive as natural vegetation would be in the same
  place, and has its vegetation recently been lost? Observed net primary productivity (NPP) is
  compared with a modelled natural potential, and the pixel's vegetation greenness with its own
  2009-2018 baseline.
- **Structural integrity**: how little are the land and its surrounding landscape modified by
  people? Derived from the Human Modification Index.
- **Compositional integrity**: how much of the original species community remains? Taken from
  the Biodiversity Intactness Index.

A score of 1 means the pixel is intact in all three dimensions. The weakest dimension sets the
ceiling of the index, so a productive plantation with little native biodiversity cannot score
high.

## The product

| | |
|---|---|
| Version | EII v2.0 |
| Coverage | Global land, 60° S to 80° N |
| Resolution | CLMS 300 m grid (1/336°, EPSG:4326) |
| Reference years | 2021, 2022, 2023, 2024, 2025 |
| Bands | `eii`, `functional_integrity`, `structural_integrity`, `compositional_integrity`, `limiting_component`, `wetland_flag`, `fi_state`, `fi_change`, `fi_change_sensors` |
| Access | Python client (`eii.client`) or Earth Engine image collections |

Each reference year uses observed NPP and vegetation greenness of the three years ending in it
and biodiversity intactness of the same year. Results are computed per 10° tile and published as
cloud-optimised GeoTIFFs registered as Earth Engine collections.

```python
import ee
from eii.client import get_stats

ee.Initialize(project="your-gee-project-id")
get_stats((11.0, 47.5, 11.5, 48.0), reference_year=2024)
```

## Where to go next

- [Getting Started](getting-started.md): install, authenticate and run a first query.
- [Methodology](methodology/index.md): how each component and the index are computed.
- [Validation](methodology/validation.md): model skill and agreement with independent indicators.
- [Data & Assets](data-assets.md): collections, bands, files and source datasets.
- [Limitations](limitations.md): what the index can and cannot tell you.
- [Changes from v1](v1-to-v2.md): what differs from version 1.
- [API Reference](api/index.md): function documentation.
- [Methodology report (PDF)](report/eii_v2_methodology.pdf): the full method in one document.
