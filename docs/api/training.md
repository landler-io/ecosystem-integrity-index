# eii.training

Sampling design, cross-validation folds and model card of the potential-NPP model. Sampling, the
training table, the fit and the global prediction are in
[`pipelines/training`](https://github.com/landler-io/ecosystem-integrity-index/tree/main/pipelines/training);
the method is described on the [potential NPP page](../methodology/potential_npp.md).

## Sampling design

::: eii.training.sampling
    options:
      members:
        - draw_design
        - allocate_quotas
        - thin_by_distance
        - select_samples
        - gap_fill

## Cross-validation folds

::: eii.training.validation
    options:
      members:
        - spatial_block_ids
        - group_kfold_folds
        - spatial_buffer_mask

## Model card

::: eii.training.model_card
    options:
      members:
        - ModelCard
        - load_model_card
