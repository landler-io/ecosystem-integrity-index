# eii.compute

Full EII computation utilities (on-the-fly).

## Aggregation

::: eii.compute.integrity
    options:
      members:
        - calculate_eii
        - combine_components

## Functional integrity (NPP)

::: eii.compute.npp
    options:
      members:
        - calculate_functional_integrity
        - load_natural_npp
        - load_natural_npp_tiles
        - load_npp_diff_percentiles
        - setup_predictor_stack
        - setup_response

## Structural integrity

::: eii.compute.structural
    options:
      members:
        - calculate_structural_integrity

## Compositional integrity

::: eii.compute.compositional
    options:
      members:
        - calculate_compositional_integrity
