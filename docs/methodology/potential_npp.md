# Potential NPP Model

Functional integrity needs, for every land pixel, the NPP that natural vegetation would reach
there. This page describes the model that provides it: where it learns, from which predictors,
which model class, how its uncertainty is calibrated and how well it performs.

## Design principles

- **Learn only from natural land.** The model is trained on the least modified land on Earth and
  applied everywhere, so its prediction for a cropland pixel is the NPP natural vegetation would
  reach under the same conditions.
- **Use only predictors that land use cannot change.** A predictor derived from present-day
  vegetation, soil management or land cover would let the model learn the current state of a
  pixel and hide its loss. Every predictor describes climate, terrain or coastal setting.
- **One global model.** Region or biome terms gain almost no skill and make predictions jump at
  region boundaries.
- **Calibrated uncertainty.** The prediction interval is the tolerance of functional integrity, so
  its coverage is calibrated per biome.

## Training data

**Candidate pool.** One candidate per 30-arc-second cell of global land (163.6 million). Water,
snow and ice, built-up land and cropland are excluded; bare and sparsely vegetated land stays, so
deserts are represented.

**Reference sites.** A candidate is a reference site if it lies inside strictly protected areas
(WDPA), or if its maximum Human Modification Index over 1990 to 2020 is below 0.1 and it is more
than about 3 km from land with HMI above 0.4. Sites are thinned to at least 2 km apart.

**Gap fill.** In some regions no natural land is left in the semi-arid zone. In the Sahel, on the
Loess Plateau and in semi-arid India the reference sites of a 4° block are bare desert (Lake Chad
block: 77 mm rain, NPP 5 gC m⁻² yr⁻¹), while the semi-arid land of the same block is all in use
(311 mm, NPP 278). Trained on such sites alone, the model would extrapolate the potential from
the desert, and observed NPP would be 2-3 times the potential over large areas. The training set
is therefore topped up with
*relaxed* sites: candidates that pass every filter except the natural-land mask and have HMI 2022
below 0.2. Strata are 4° blocks times aridity class (edges 0.05, 0.1, 0.2, 0.5 and 0.65 of the
aridity index). In every stratum with fewer than 40 reference sites, relaxed sites fill the
shortfall, lowest HMI first, and are fitted at weight 0.5. Half of all relaxed sites are set aside
beforehand and never fitted; they are the out-of-sample check below. The fill adds 8,896 sites
(6 % of training). In these regions "potential" therefore partly means lightly used land.

The model is fitted on all 142,480 reference sites plus the gap-fill sites.

**Response.** Mean annual CLMS DMP V2.0.2 NPP over 2019 to 2025 (the Sentinel-3 OLCI years), with
at least 5 of the 7 years valid.

**Evaluation sample.** Skill is measured on a stratified draw of 30,000 reference sites: a floor of
1,000 per biome, at most 40 sites per biome and 4° block, and at least 1 % of each biome's quota
per environmental cluster. Without the floor, small biomes such as mangroves and tropical
coniferous forest would barely enter the scores.

## Predictors

Fifteen predictors on a 300 m stack, in model order:

| Group | Predictor | Source |
|---|---|---|
| Climate | `mean_annual_temp` | CHELSA V2.1 1981-2010, bio1 |
| | `temp_seasonality` | CHELSA V2.1, bio4 |
| | `annual_precip` | CHELSA V2.1, bio12 |
| | `precip_seasonality` | CHELSA V2.1, bio15 |
| Water balance | `aridity` | CHELSA V2.1 annual precipitation / Penman PET |
| Climatic NPP | `chelsa_npp` | CHELSA V2.1 NPP 1981-2010 |
| Terrain | `elevation` | MERIT DEM |
| | `cti` | Geomorpho90m compound topographic index |
| | `tpi` | Geomorpho90m topographic position index |
| | `northness`, `eastness` | Geomorpho90m |
| Coastal setting | `dist_coast` | OpenStreetMap coastline, km, capped at 50 |
| | `tidal_range` | EOT20, 0 beyond 50 km of the coast, capped at 16 m |
| | `sss` | WOA23 sea surface salinity 1991-2020, 0 beyond 50 km |
| | `log_upa_10km` | MERIT Hydro upstream area, log10(1 + 10 km maximum), 0 beyond 50 km |

The coastal group lets the model learn mangrove and tidal-marsh productivity; it raises mangrove
skill (block-CV R² 0.22 to 0.30) without changing global skill.

### Predictors not used

Each candidate group was tested for skill (gain in within-region and leave-biome-out R²) and for
land-use leakage (whether it separates cropland from adjacent forest at the same elevation).

| Predictor | Reason |
|---|---|
| Soil pH (SoilGrids) | Mapped from covariates that include tree cover, MODIS EVI and NPP, and responds to liming and cultivation. It roughly doubled the land-cover pattern in the model residuals. Survey-based alternatives (GLiM lithology, HWSD) add no skill. |
| Sand, clay | No skill. |
| Latitude, longitude | No skill; worse transfer between regions. |
| Terrain ruggedness, regional TPI | No skill beyond the terrain set. |
| Slope | Gave steep ground a potential about 13 % higher in a Bavarian test region, where observed forest NPP is flat across terrain. |
| Biome indicator | Raises leave-biome-out bias. |
| Wetness (HAND, upstream area) | No skill; did not reduce the wetland bias. |

## Model

**Mean model.** One global LightGBM regressor with Tweedie loss:

| Parameter | Value |
|---|---|
| Objective | Tweedie, variance power 1.35 |
| Learning rate | 0.03 |
| Boosting rounds | 3,000 |
| Leaves | 31 |
| Minimum child samples | 20 |
| Column subsample | 0.8 |
| Row subsample | none |
| L2 penalty | none |

The Tweedie loss suits the skewed distribution of NPP and keeps the per-biome bias low: the
worst median biome deficit under block CV is 0.125 after bias correction (see [Skill](#skill)).

**Bias correction.** The Tweedie model under-predicts by about 2 % overall. The published
potential is

$$P = c \cdot \hat P, \qquad c = \frac{\sum A}{\sum \hat P} = 1.0203,$$

with $\hat P$ the block cross-validated out-of-fold predictions and $A$ the observed response,
both summed over the 30,000-site evaluation sample.

**Prediction intervals.** A quantile random forest (200 trees, minimum leaf 5, 5 features per
split, 70 % sample per tree) provides the interval shape. The bounds are asymmetric:

$$L = \max\left(0,\ P - k_b \max(q_{50} - q_{16},\ \varepsilon)\right), \qquad U = P + k_b \max(q_{84} - q_{50},\ \varepsilon)$$

$q$ are the forest's quantiles, $\varepsilon = 1$ gC m⁻² yr⁻¹ keeps the interval open where
quantiles coincide, and $k_b$ is a per-biome factor: the conformal-type quantile of the
pooled out-of-fold scores of the five block-CV folds, so that the interval covers 68 % of
out-of-fold observations in each biome. The 90 % interval uses the same construction with its own
factors. Biomes with fewer than 50 sites and pixels whose biome is unknown take the pooled factor
(but see the coastal labels below). The published models are refitted on all sites and the factors
applied to them. This is not split conformal prediction: the scores come from the fold models,
sites are spatially dependent and the lower bound is clipped at 0, so coverage is measured
(cross-fitted, each fold scored with factors from the other folds), not guaranteed.

**Natural inter-annual variability.** A companion LightGBM on the same 15 predictors and
reference sites (Tweedie 1.35, 15 leaves, minimum child samples 100, 1,000 rounds; no gap fill)
predicts the
coefficient of variation of the 2019 to 2025 annual sums (sites with at least 5 valid years and
mean NPP of at least 50 gC m⁻² yr⁻¹), times its own bias correction.

## Skill

Out-of-fold scores on the 30,000-site evaluation sample. **Block CV** holds out 4° blocks with a
1° buffer; **leave-realm-out** (LORO) holds out an entire biogeographic realm, which tests transfer
to land without nearby training data. Within-block R² removes the between-block signal and
measures local skill. Biome deficit is the median relative deficit $(P - A)/P$ per biome.

| | R² | Within-block R² | Pooled bias | Mean / max abs. biome deficit |
|---|---|---|---|---|
| Block CV, raw | 0.906 | 0.500 | -2.0 % | 0.050 / 0.141 |
| Block CV, bias-corrected | 0.906 | 0.497 | 0.0 % | 0.048 / 0.125 |
| LORO, raw | 0.867 | 0.423 | -5.0 % | 0.136 / 0.297 |
| LORO, bias-corrected | 0.861 | 0.412 | +0.4 % | 0.117 / 0.237 |

**Interval coverage** (target 68 % / 90 %):

| | 68 % interval | 90 % interval |
|---|---|---|
| Block CV, pooled (worst biome) | 0.68 (0.67) | 0.90 (0.89) |
| LORO, pooled (worst biome) | 0.70 (0.64) | 0.90 (0.86) |

**Gap-fill check.** On the held-out relaxed sites, observed NPP exceeds the potential by more
than 1.5 times in one 4° block, against seven blocks for a model fitted without the gap fill. The
Lake Chad block is at 1.49 times (2.57 without the fill); the West Loess Plateau block stays at
1.98. The gap fill does not change skill on reference sites (block-CV R² 0.906, 0.905 without).

The natural-CV model reaches block-CV R² 0.63 (within-block 0.25).

## Production

The model is applied on the 300 m predictor stack, per 10° tile
([`pipelines/training`](https://github.com/landler-io/ecosystem-integrity-index/tree/main/pipelines/training)).
The output has
five bands in gC m⁻² yr⁻¹: `potential_npp`, `lower_68`, `upper_68`, `lower_90`, `upper_90`, plus a
separate `natural_cv` map. Both are static across reference years. The five-band map is published
([Data & Assets](../data-assets.md#potential-npp)). The `natural_cv` map, needed to reproduce
FI_state, is not distributed: fit the model with `fit_potential_npp_v2.py` and predict it with
`predict_potential_npp_v2.py --natural-cv`. Pixels with any predictor
missing stay empty; this leaves twelve small-island tiles without MERIT terrain unpredicted.

The predictor list and model settings live in the model card
[`src/eii/training/model_cards/v2.json`](https://github.com/landler-io/ecosystem-integrity-index/tree/main/src/eii/training/model_cards),
the single source of `eii.training.settings.PREDICTOR_VARIABLES`. Sampling, the training table,
the fit and the prediction are in
[`pipelines/training`](https://github.com/landler-io/ecosystem-integrity-index/tree/main/pipelines/training).
The training table is public (see [Data & Assets](../data-assets.md#training-table)).

## Known weaknesses

- **Mediterranean and tropical coniferous forest** are predicted low (median deficit -0.12 in the
  Mediterranean biome under block CV); **desert** is predicted high (+0.12).
- **Unseen realms**: in a realm without training data the model under-predicts by about 5 %; the
  global bias correction covers only the 2 % seen under block CV.
- **Terrain**: in a Bavarian test region, forest observed/potential is still 4-6 % higher on flat
  lowland than on moderate terrain, mainly through the 1 km climate layers. Cropland losses on
  valley floors are understated accordingly.
- **Semi-arid regions without natural land** (Sahel, Loess Plateau, semi-arid India) take their
  potential partly from lightly used land through the gap fill. Irrigated cropland there has no
  natural check sites.
- **Europe** has few reference sites in its productive lowlands; leave-Europe-out R² is 0.61-0.71.
- **Wetlands and peatlands** are overestimated; see the
  [wetland flag](functional_integrity.md#wetland-flag).
- **Coastal biome labels**: the intended fallback for pixels without a biome is the pooled factor,
  but coastal land whose 5 km RESOLVE cell is ocean is labelled tundra by the raster mapping and
  takes the tundra interval factors (3.9 % of land, 0.8 % below 60° latitude). Before clipping at
  0, its intervals are about 2 % (68 %) and 11 % (90 %) narrower than with the pooled factor;
  potential NPP itself is unaffected.

## References

- Karger, D. N., et al. (2017). Climatologies at high resolution for the earth's land surface
  areas. *Scientific Data*, 4, 170122.
- Hart-Davis, M. G., et al. (2021). EOT20: a global ocean tide model from multi-mission
  satellite altimetry. *Earth System Science Data*, 13, 3869-3884.
- Amatulli, G., et al. (2020). Geomorpho90m, empirical evaluation and accuracy assessment of
  global high-resolution geomorphometric layers. *Scientific Data*, 7, 162.
- Yamazaki, D., et al. (2017). A high-accuracy map of global terrain elevations. *Geophysical
  Research Letters*, 44, 5844-5853.
- Meinshausen, N. (2006). Quantile regression forests. *Journal of Machine Learning Research*, 7,
  983-999.
