# Changelog

## [2.0.0] - 2026-11-16

EII v2.0: an annual global index for the reference years 2021-2025. From this release the
package version matches the data version. The method is described in the
[documentation](https://landler-io.github.io/ecosystem-integrity-index/) and the
[methodology report](docs/report/eii_v2_methodology.pdf); the differences from v1 are listed in
[Changes from v1](docs/v1-to-v2.md).

### Added

- EII v2 data: one map per reference year 2021-2025 on the CLMS 300 m grid, published as
  cloud-optimised GeoTIFFs and registered as the Earth Engine collection
  `projects/landler-open-data/assets/products/eii/v2_0/release/global`.
- Functional integrity as the product of a state score (observed against potential NPP, with a
  tolerance from model uncertainty and natural variability) and a change term (EVI of MODIS
  Terra and VIIRS against the pixel's own 2009-2018 baseline).
- Structural integrity as landscape quality from the Human Modification Index with a logistic
  habitat-quality curve, calibrated against ESA WorldCover.
- New bands: `fi_state`, `fi_change`, `fi_change_sensors`, `limiting_component` and
  `wetland_flag`.
- Potential NPP with 68 % and 90 % prediction intervals as a separate Earth Engine collection,
  and the public training table of the potential NPP model.
- `reference_year` parameter in `get_stats`, `get_zonal_stats`, `get_raster`, `get_layers`,
  `compare_methods` and `get_modulated_eii`. It defaults to the latest year.
- `get_layers("limiting_component")`, `get_layers("wetland_flag")` and
  `get_layers("fi_change")` (returns `fi_state`, `fi_change` and `fi_change_sensors`).
- `exclude_wetlands` in `get_stats` and `get_zonal_stats`; the flagged share is reported by
  default.
- Dataset version selection: `version=` on every retrieval function, `set_version()`,
  `get_version()` and the `use_version()` context manager.
- `eii.compute`: numpy implementation of the v2 calculation (`deviation_tolerance_integrity`,
  recentring, the change-term functions, `structural_integrity`, `habitat_quality`,
  `compositional_integrity`, `combine_components`, `limiting_component`, `wetland_flag`).
- `eii.training`: sampling design, spatial cross-validation folds and the model card of the
  potential NPP model.
- `pipelines` extra with the dependencies of the reproduction pipelines in `pipelines/`.
- Methodology report (PDF) and rewritten documentation.

### Changed

- The default dataset is v2. v1 remains available with `version="v1"` or
  `set_version("v1")`.
- `compare_methods` recombines the published components; only `min_fuzzy_logic` is published
  as the `eii` band.

### Breaking

- Calls without `version` return v2 instead of v1, and v2 values differ from v1 (see
  [Changes from v1](docs/v1-to-v2.md)). Pin `version="v1"` to keep v1 results.
- The client is retrieval-only. Removed parameters: `compute_mode`, `aggregation_method`,
  `npp_year_range`, `include_seasonality`, `bii_year`, `absolute_diff_percentile` and the asset
  path overrides (`natural_npp_asset_path`, `observed_npp_asset_path`,
  `structural_asset_path`). Positional arguments after `geometry` shift accordingly; pass them
  by keyword.
- `eii.compute` is pure numpy. The Earth Engine compute functions are removed; the v1 Earth
  Engine code remains in the [v0.1.1 release](https://github.com/landler-io/ecosystem-integrity-index/tree/v0.1.1).
- Natural Capital modulation moved from `eii.compute.modulation` to `eii.modulation`.

### Removed

- `eii.utils` and `eii._utils` (Earth Engine helpers).
- `eii.compute.compositional`, `eii.compute.grid`, `eii.compute.inference`, `eii.compute.npp`,
  `eii.compute.stats` and `eii.training.model` (v1 Earth Engine computation and training).

## [0.1.1] - 2026-04-20

- fix data extraction for polygons smaller than the pixel resolution


## [0.1.0] - 2026-01-01

- Initial release
- Core EII data access functions
- Basic documentation
