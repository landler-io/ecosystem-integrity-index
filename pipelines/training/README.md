# Training

Everything that fits or calibrates a parameter of EII v2.0: the potential-NPP model of
functional integrity (with its prediction intervals and natural inter-annual CV), the constants
of the FI change term and the habitat-quality curve of structural integrity.

## Potential NPP

Potential NPP is the NPP a pixel would have without human modification, predicted from climate,
terrain and coastal setting by a model trained at near-pristine reference sites.

### 1. Training table

The published table can be used directly:
`gs://landler-opendata-gee-assets/products/eii/v2_0/training/npp_training_samples_v2_wp11_localstack_dmpv2_olci_20261104.csv.gz`
(requester pays; 218,721 rows). It was built in these steps:

| Step | Script | Output |
| --- | --- | --- |
| Candidates and Poisson selection on the 300 m HMI grid | `build_candidate_pool_local.py candidates --wdpa WDPA.tif`, `select`, `points` | `points.csv.gz` |
| Predictors and annual CLMS NPP v1.1 sampled at the points in Earth Engine | `sample_candidate_pool_ee.py points.csv.gz pool.csv.gz` | `npp_candidate_pool_v2_local_seed42_20261027.csv.gz` (published, see `pipelines/README.md`) |
| Eligibility, thinning, quotas; `--keep-relaxed` also keeps eligible non-pristine rows (`reference == False`) | `sample_training_v2.py select --keep-relaxed` | `npp_training_samples_v2_wp11_blockcap_seed42_20261027.csv.gz` |
| Predictors from the local predictor stack | `extract_stack_at_points.py`, `training_table_localstack.py` | `npp_training_samples_v2_wp11_localstack_20261028.csv.gz` |
| Response from the CLMS DMP V2.0.2 annual sums, OLCI years 2019-2025 | `training_table_dmp_v2.py` | `npp_training_samples_v2_wp11_localstack_dmpv2_olci_20261104.csv.gz` |

The values sampled in Earth Engine only screen eligibility (valid NPP years) and define the
climate clusters of the quota draw; the last two steps replace every predictor and response
value. `sample_candidate_pool_ee.py` samples the 19 predictors of the v1 stack (WorldClim 1.4,
Global Aridity Index, MERIT DEM, Geomorpho90m, SoilGrids, CHELSA NPP) and the annual CLMS NPP
v1.1 sums 2014-2025 per chunk with interactive Earth Engine requests (resumable chunk cache; try
`--max-chunks 1` first); the asset ids are in its docstring. The annual CLMS NPP v1.1 collection is
no longer published, so a rebuild starts from the published candidate pool. The candidates step needs a WDPA high-protection raster on the
GHM grid that you build from your own WDPA download (WDPA may not be redistributed; the rule is
in the `build_candidate_pool_local.py` docstring); it reads the RESOLVE ecoregions from the
requester-pays bucket (`EII_BILLING_PROJECT`). `sample_training_v2.py` needs
`biome_land_shares.csv` (biome land shares from a global random sample of land points).

### 2. Fit

```bash
python pipelines/training/fit_potential_npp_v2.py TABLE.csv.gz models/potential_npp_v2
```

writes `lgbm_mean.txt` (LightGBM Tweedie mean model), `qrf.joblib` (quantile regression forest
of the intervals), `lgbm_interannual_cv.txt` (natural inter-annual CV) and `model_card.json`
(predictors, parameters, bias correction `c`, per-biome conformal factors of the 68 % and 90 %
intervals). Training rows are the reference rows plus gap-fill rows
(`eii.training.sampling.gap_fill`, weight 0.5); calibration uses buffered 4-degree block
cross-validation on the 30,000 `selected` rows. It takes a few minutes on 64 cores and
reproduces the published model to floating-point precision (identical forest; LightGBM
predictions within 1e-10 relative). Validation metrics are in `src/eii/training/model_cards/v2.json`.

### 3. Predict

```bash
python pipelines/training/predict_potential_npp_v2.py predictor_stack models/potential_npp_v2 potential_npp_v2
python pipelines/training/predict_potential_npp_v2.py predictor_stack models/potential_npp_v2 potential_npp_v2 --natural-cv
```

writes one GeoTIFF per 10° tile (`potential_npp`, `lower_68`, `upper_68`, `lower_90`,
`upper_90`, gC m-2 yr-1) and `potential_npp.vrt`, and the natural CV tiles and `natural_cv.vrt`.
The biome of the conformal factors comes from the RESOLVE raster synced with the predictor stack
and `pipelines/preprocessing/reference/eco_id_to_biome.json`.

## FI change term constants

```bash
python pipelines/training/derive_fi_change_constants.py fi_change_controls.csv.gz out/
```

derives `drift_factors.csv`, `h_model_{TV,TNc,Nc}.txt` and `k_c.json` from a table of
near-pristine control points (HMI < 0.05; WorldCover tree, shrub or grass) holding the
predictors, the production FI_state per reference year and annual mean EVI per sensor (columns in
the script docstring), and checks them against the checksums in `eii.compute.settings`. The
stored copies are in `src/eii/compute/fi_change_constants/`, and production uses them. The
control table of the stored constants,
`gs://landler-opendata-gee-assets/products/eii/v2_0/training/fi_change_controls.csv.gz`
(19,984 controls, requester pays), is published. It is built by

```bash
python pipelines/training/build_fi_change_controls.py select WORK --ecoregions ECOREGIONS.tif \
    --hmi HMv20240801_2022s_AA_300.tif --exclude exclusion.csv
python pipelines/training/build_fi_change_controls.py extract WORK --sensor T   # and V, N
python pipelines/training/build_fi_change_controls.py assemble WORK --run EII_RUN \
    --predictor-stack stack.vrt
```

`select` draws the controls (RESOLVE biome quotas, 5 km exclusion around earlier evaluation
points, WorldCover 2021 purity), `extract` reads the 16-day EVI composites of MODIS Terra
(MOD13Q1), VIIRS Suomi NPP (VNP13A1) and NOAA-20 (VJ113A1), and `assemble` adds annual EVI, the
predictors and FI_state. FI_state does not depend on the change-term constants, so a run with
the stored constants is a valid source. Sources and credentials are listed in the script
docstring. Three differences from the published table remain:

- The exclusion list of earlier evaluation points is not published; without it the selection
  differs.
- The published FI_state comes from a run on the potential NPP before the gap fill. Derived from
  the released FI_state, the drift factors, tolerance models and k_c are the same; the 10th
  percentiles recorded in `k_c.json` shift by up to 0.0013.
- Eleven pairs of controls have tied draw keys, so their row order is arbitrary, which changes
  the LightGBM tolerance models at bit level.

## Structural integrity: habitat-quality curve

```bash
python pipelines/training/calibrate_structural.py natural si_calibration/natural --hmi HMI.tif
python pipelines/training/calibrate_structural.py calibrate si_calibration_sample.csv.gz \
    si_calibration --hmi HMI.tif --natural si_calibration/natural/natural.vrt
```

compares structural integrity under the logistic habitat-quality curve (steepness 10, candidate
midpoints 0.3, 0.4, 0.5) with a WorldCover 2021 natural-cover target at sample points tagged with
RESOLVE biome and anthrome class, and applies the acceptance criteria C1 (no biome or anthrome
group above the land-cover bound) and C2 (wilderness keeps S near 1). Midpoint 0.4 is the gentlest
curve that passes both (`eii.compute.settings.QUALITY_MIDPOINT`). The sample points are the
29,585-point equal-area global land sample with RESOLVE biome and anthrome class, published as
`gs://landler-opendata-gee-assets/products/eii/v2_0/training/si_calibration_sample.csv.gz`
(requester pays). `build_si_calibration_sample.py` rebuilds it exactly (same points, order,
biome and anthrome):

```bash
python pipelines/training/build_si_calibration_sample.py si_calibration_sample.csv.gz \
    --ecoregions ECOREGIONS.tif --anthromes anthro2_a2000.tif
```

It draws 111,111 equal-area points (seed 42, 60° S-75° N), keeps the first 30,000 that are land
in ESA WorldCover 2021 at 1 km and complete in the record layers of the original sample (EII v1,
natural NPP v1, HMI 2022, BII 2020, the v1 predictors; 29,585 remain), and tags them with the
RESOLVE biome and the Anthromes v2 class of 2000 (Ellis et al. 2010, doi:10.7910/DVN/6KKAWS).
The land test and the record filter use interactive Earth Engine requests (about 35 minutes);
the asset ids are in the script docstring.
