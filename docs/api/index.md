# API Reference

This section contains the complete API documentation for the `eii` package.

## Client vs Compute

Use the client workflow when you want fast access to pre-computed global layers.
Use the compute workflow when you need alternative aggregation methods or custom inputs.

| Workflow | When to use | Key functions | Dependencies |
|----------|-------------|---------------|--------------|
| **Client (precomputed)** | Most users; quick stats and raster downloads | `get_stats`, `get_raster`, `get_layers` | `ecosystem-integrity-index[client]` |
| **Compute (on-the-fly)** | Custom aggregation or recalculation for an AOI | `calculate_eii`, `calculate_functional_integrity` | `ecosystem-integrity-index[compute]` |
| **Training** | Model development and validation | `train_npp_model`, `validate_model` | `ecosystem-integrity-index[training]` |

## Modules

### [eii.client](client.md)

Client interface for EII data access from Google Earth Engine.

- `get_stats()` - Extract EII statistics for a geometry
- `get_raster()` - Download EII raster for a geometry
- `get_layers()` - Get EII and/or integrity component layers as ee.Image
- `ASSETS` - Asset paths and metadata

### [eii.client.analysis](analysis.md)

Analysis utilities for EII data.

- `get_zonal_stats()` - Zonal statistics for feature collections
- `compare_methods()` - Compare different aggregation methods

### [eii.compute](compute.md)

Core computation utilities for on-the-fly EII.

- `calculate_eii()` - Compute EII from component layers
- `calculate_functional_integrity()` - NPP-based functional integrity
- `calculate_structural_integrity()` - Structural integrity (core area)
- `calculate_compositional_integrity()` - Compositional integrity (BII)

### [eii.training](training.md)

Model training utilities for advanced users.

## Quick Reference

```python
# Most common imports
from eii.client import get_stats, get_raster, get_layers, ASSETS
from eii.client import get_zonal_stats, compare_methods
from eii.compute import calculate_eii, calculate_functional_integrity
```
