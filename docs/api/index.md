# API Reference

The `eii` package has three parts.

| Part | When to use | Key functions | Extra |
|---|---|---|---|
| **Client** | Retrieve published EII layers and statistics | `get_stats`, `get_zonal_stats`, `get_raster`, `get_layers` | `[client]` |
| **Compute** | Reproduce or adapt the calculation on your own arrays | `deviation_tolerance_integrity`, `fi_change`, `structural_integrity`, `combine_components` | `[compute]` |
| **Training** | Sampling design and validation of the potential NPP model | `draw_design`, `gap_fill`, `load_model_card` | `[training]` |

Install a part with its extra, for example
`pip install "ecosystem-integrity-index[client] @ git+https://github.com/landler-io/ecosystem-integrity-index"`.

The global layers are computed with the compute functions
([`pipelines/production`](https://github.com/landler-io/ecosystem-integrity-index/tree/main/pipelines/production))
and published as Earth Engine collections, which the client reads. The client needs an
initialised Earth Engine session; compute and training run on numpy arrays without Earth Engine.
Natural Capital modulation lives in `eii.modulation` and is exposed through the client
([Local Modulation](../methodology/local_modulation.md)).

## Modules

### [eii.client](client.md)

- `get_stats()`: values at a point or statistics for an area
- `get_raster()`: download EII and component rasters
- `get_layers()`: EII, components, limiting component, wetland flag and the FI parts as `ee.Image`
- `set_version()`, `use_version()`, `get_version()`, `get_asset_path()`: dataset version and
  asset paths
- `get_modulated_eii()` and related: [Local Modulation](../methodology/local_modulation.md)

### [eii.client.analysis](analysis.md)

- `get_zonal_stats()`: statistics for many features at once
- `compare_methods()`: recombine the components with each aggregation method

### [eii.compute](compute.md)

- `observed_npp_window()`, `observed_npp()`: observed NPP of a reference year
- `recentring_sums()`, `recentring_factors()`, `recentre_potential()`: per-biome recentring
- `deviation_tolerance_integrity()`: the state score of functional integrity
- `annual_evi()`, `correct_drift()`, `change_ratio()`, `combine_ratios()`, `tolerance_h()`,
  `fi_change()`: the change term
- `functional_with_change()`: functional integrity, FI_state × FI_change
- `wetland_flag()`: GLWD-based wetland flag
- `structural_integrity()`, `habitat_quality()`: landscape quality from HMI
- `compositional_integrity()`: clamped BII
- `combine_components()`, `limiting_component()`: EII and its limiting component

### [eii.training](training.md)

Sampling design, spatial cross-validation folds and the model card of the potential NPP model.

## Quick reference

```python
from eii.client import get_stats, get_raster, get_layers, get_zonal_stats
from eii.compute import deviation_tolerance_integrity, structural_integrity, combine_components

stats = get_stats((11.0, 47.5, 11.5, 48.0), reference_year=2024)
v1_stats = get_stats((11.0, 47.5, 11.5, 48.0), version="v1")
```
