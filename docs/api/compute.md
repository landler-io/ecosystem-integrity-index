# eii.compute

The EII calculation on local rasters: pure numpy functions, one per component and for their
combination. The global run over 10 degree tiles is in
[`pipelines/production`](https://github.com/landler-io/ecosystem-integrity-index/tree/main/pipelines/production);
its results are published as Earth Engine collections and read with [`eii.client`](client.md).

Install with `pip install "ecosystem-integrity-index[compute] @ git+https://github.com/landler-io/ecosystem-integrity-index"`.

## Functional integrity: state score

::: eii.compute.functional
    options:
      members:
        - observed_npp_window
        - observed_npp
        - deviation_tolerance_integrity
        - recentring_sums
        - recentring_factors
        - recentre_potential
        - wetland_flag

## Functional integrity: change term

::: eii.compute.change
    options:
      members:
        - encode_modis_evi
        - encode_viirs_evi
        - annual_evi
        - complete_years
        - natural_class
        - correct_drift
        - chain_factor
        - chained
        - baseline_median
        - change_ratio
        - combine_ratios
        - tolerance_h
        - fi_change
        - functional_with_change
        - series_for
        - constants_sha256

## Structural integrity

::: eii.compute.structural
    options:
      members:
        - structural_integrity
        - habitat_quality
        - landscape_quality
        - structural_halo_m
        - ground_disk
        - exponential_kernel

## Compositional integrity and aggregation

::: eii.compute.integrity
    options:
      members:
        - compositional_integrity
        - combine_components
        - limiting_component

## Settings

::: eii.compute.settings
