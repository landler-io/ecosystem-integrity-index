# Functional Integrity

Functional integrity (FI) asks whether land is as productive as natural vegetation would be in
the same place, and whether its vegetation has recently been lost. It is the product of two
scores:

$$\mathrm{FI} = \mathrm{FI}_{\text{state}} \times \mathrm{FI}_{\text{change}}$$

- **FI_state** compares observed **net primary productivity (NPP)** with the NPP that natural
  vegetation would reach under the local climate and terrain. NPP is the carbon that plants fix
  by photosynthesis minus what they respire, in gC m⁻² yr⁻¹. A clear deficit points to
  degradation, disturbance or conversion; a clear surplus points to irrigation, fertilisation or
  other inputs that natural vegetation would not receive.
- **FI_change** compares each pixel's vegetation greenness (EVI) over the last three years with
  its own 2009-2018 baseline, and scores a loss beyond what near-pristine land shows from noise
  alone.

Both scores are 1 inside a tolerance and fall linearly to 0 beyond it. The published band
`functional_integrity` is the product; `fi_state`, `fi_change` and `fi_change_sensors` hold its
parts.

The two scores see different things. FI_state measures the level of productivity against a
model, so its tolerance has to absorb the model's error, and only strong, lasting departures
lower it. FI_change measures change at the pixel against its own past. Model error cancels, so
its tolerance is much narrower (a median of about 5 % of baseline EVI on near-pristine land) and
it registers clearing, fire and other vegetation loss of the last few years that FI_state
misses. Land that was converted before the baseline and has stayed in the same use scores
FI_change = 1, and FI_state registers it only where its productivity departs from the potential
beyond the tolerance; long-standing conversion is carried by structural and compositional
integrity.

## Observed NPP

**Source.** Copernicus Land Monitoring Service (CLMS) Dry Matter Productivity, 300 m, 10-daily,
consolidated reprocessing V2.0.2 (RT6). The annual processing
([`pipelines/preprocessing`](https://github.com/landler-io/ecosystem-integrity-index/tree/main/pipelines/preprocessing)):

1. masks missing and water values and dekads whose quality flag marks FAPAR or meteorological
   input as invalid;
2. anchors each 10-day value at the middle of its dekad, interpolates to daily values and sums
   them per calendar year (all 36 dekads of one product version are required);
3. converts dry matter (kg DM ha⁻¹) to carbon (gC m⁻²).

Stored annual sums are integers equal to gC m⁻² yr⁻¹ × 10,000; `NPP_STORED_SCALE = 1e-4`
converts them.

**Sensor years.** DMP V2.0.2 comes from PROBA-V until the end of 2018 and from Sentinel-3 OLCI
from 2019. The two have no overlap year, and against MODIS the step between them is +26 % in
tundra, +18 % in boreal forest and +11 % in tropical moist forest, with a seasonal pattern that
one factor per biome cannot correct. EII therefore uses OLCI years only (2019 onward). The first
reference year is 2021, the first with a full three-year OLCI window.

**Window.** Observed NPP of a reference year is the mean of the three annual sums ending in it
(`observed_npp_window`; 2024 uses 2022 to 2024). Negative values are missing, and a pixel needs at
least two valid years (`observed_npp`). The three-year mean damps weather-driven fluctuation.

## Potential NPP and natural variability

The [potential NPP model](potential_npp.md) gives, per pixel:

- $P$, the expected NPP of natural vegetation under the local climate, terrain and coastal
  setting, as the mean over 2019 to 2025;
- $L$ and $U$, the lower and upper bounds of its 68 % prediction interval (calibrated per biome,
  asymmetric);
- $v$, the modelled natural inter-annual coefficient of variation of NPP.

$v$ is modelled from the same predictors at the same reference sites, never from the observed
variability of the pixel itself: on managed land observed variability reflects the land use
(harvests, rotations, irrigation), not natural variability.

## Recentring

$P$ is one static map, trained on the 2019 to 2025 mean, while each reference year is scored with
a three-year window. Regional weather and residual trends in the satellite record move the level
of observed NPP from window to window, most strongly at high latitudes (on intact boreal land,
CLMS NPP relative to MODIS rises from 1.07 in 2019 to 1.15 in 2024). Without a correction, intact
boreal forest and tundra would read as degraded in some windows, and the annual series would
swing with the record rather than with the land.

For each observed window and each RESOLVE biome $b$, $P$, $L$ and $U$ are therefore scaled by

$$f_b = \frac{\sum_{i \in R_b} A_i}{\sum_{i \in R_b} P_i},$$

the ratio of mean observed to mean potential NPP over the biome's reference pixels $R_b$.

- **Reference pixels** have HMI below 0.05 (GHM 2022) and Biodiversity Intactness of at least
  0.95 in the reference year. They are selected from the inputs, not from the calculated structural
  and compositional scores, so recentring does not depend on those scores. HMI and BII also
  determine structural and compositional integrity, so the three components still share this
  pressure information.
- **Ratio of means**, not median of ratios: about half of near-pristine desert pixels observe
  NPP 0.
- **Small biomes** with 200 reference pixels or fewer, and pixels without a biome, keep $f = 1$.
- **Two passes**: the global run first accumulates the per-biome sums over all tiles, then
  applies the factors. The factors are recorded in the tags of every output tile.

Recentring sets the level per biome and window; it does not touch the pixel-by-pixel comparison.
Local loss inside an intact biome still scores low. The price is that change uniform across a
whole biome within one window, such as a continent-wide drought, is absorbed: FI_state measures
function relative to intact land of the same biome in the same years.

## State score

With observed NPP $A$, recentred potential $P$ and bounds $L$, $U$, natural CV $v$ and the floor
$F_0$ = 50 gC m⁻² yr⁻¹:

$$D = \max(P, F_0), \qquad d = \frac{|A - P|}{D}$$

$$(u, k, z) = \begin{cases} \left(\dfrac{P - L}{D},\ k_d,\ z_d\right) & A \le P \ \text{(deficit)} \\[2ex]
\left(\dfrac{U - P}{D},\ k_s,\ z_s\right) & A > P \ \text{(surplus)} \end{cases}
\qquad t = \min\left(k\,\sqrt{u^2 + v^2},\ 0.999\,z\right)$$

$$\mathrm{FI}_{\text{state}} = \operatorname{clamp}\left(1 - \frac{d - t}{z - t},\ 0,\ 1\right)$$

- **Tolerance** $t$: the relative half-width of the prediction interval on the side of the
  deviation, combined with the natural CV. Inside it, FI_state = 1: a deviation that model error
  and weather can explain is not evidence of a change in function.
- **Deficit** multiplier $k_d$ = 1.0; **surplus** multiplier $k_s$ = 1.5. The surplus side is
  wider because surplus produces most false alarms on intact land.
- **Zero point** $z$: FI_state reaches 0 at $d = z$. A deficit reaches 0 at $z_d = 1$, total
  loss of NPP where $P \ge 50$. A surplus reaches 0 at $z_s = 2$, observed NPP three times the
  potential where $P \ge 50$ and $P + 100$ gC below it. Between
  the tolerance and the zero point the score falls linearly, so the surplus ramp is about half as
  steep as the deficit ramp.
- **Floor**: below 50 gC m⁻² yr⁻¹ the relative error of $P$ is about 100 %, so deviations are
  measured against 50 gC instead. Natural desert scores 1, while an irrigated field in the same
  desert scores 0 once it produces 100 gC m⁻² yr⁻¹ more than the potential. Below the floor
  FI_state depends on the absolute deviation, not the relative one, and a deficit cannot reach 0:
  with $P = 10$ and a lower bound of 0, even a total loss stays inside the tolerance.

### Why surplus counts, and why its ramp is longer

Most conversions that change NPP lower it, but irrigation, fertilisation and some plantations
raise it above anything natural vegetation reaches at the same site. Surplus carries most of the
land-use signal of FI_state: a deficit-only score separates wild from used land by 0.007-0.011 in
mean score, a two-sided score by 0.04-0.08.

A surplus reaches 0 only at three times the potential ($z_s = 2$), not at twice the potential as
symmetry with the deficit would suggest. Irrigation and local model error readily double NPP: a
zero point at twice the potential scores 17 % of the land in four South Asian tiles as 0, two
thirds of it irrigated cropland and much of the rest near-pristine land where the potential model
is locally too low (Thar desert, alpine meadows). With $z_s = 2$, mean FI_state on irrigated
cropland in Punjab, central India and the Central Valley is 0.32-0.75 (2024). Irrigated
agriculture in true deserts (Sindh, the Nile valley) stays near 0 for any linear ramp, because its
surplus is about ten times the floor.

### Worked example

$P$ = 500, $L$ = 400, $U$ = 620 gC m⁻² yr⁻¹, $v$ = 0.10.

| Observed $A$ | Side | $d$ | $t$ | $z$ | FI_state |
|---|---|---|---|---|---|
| 450 | deficit | 0.10 | 0.22 | 1 | 1.00 |
| 300 | deficit | 0.40 | 0.22 | 1 | 0.77 |
| 100 | deficit | 0.80 | 0.22 | 1 | 0.26 |
| 600 | surplus | 0.20 | 0.39 | 2 | 1.00 |
| 800 | surplus | 0.60 | 0.39 | 2 | 0.87 |
| 1000 | surplus | 1.00 | 0.39 | 2 | 0.62 |
| 1250 | surplus | 1.50 | 0.39 | 2 | 0.31 |
| 1500 | surplus | 2.00 | 0.39 | 2 | 0.00 |

## Change term

FI_change scores vegetation loss at a pixel against the pixel's own past. The input is the
Enhanced Vegetation Index (EVI) of the 16-day composites of three satellite sensors:

| Code | Sensor | Product | Native grid | Record |
|---|---|---|---|---|
| T | MODIS Terra | MOD13Q1 v061 | 250 m | from 2000 |
| V | VIIRS on Suomi NPP | VNP13A1 v002 | 500 m | 2012 to 2 November 2026 |
| N | VIIRS on NOAA-20 | VJ113A1 v002 | 500 m | from 2018 |

### Annual EVI

For each sensor and pixel, the annual value is the mean EVI of the 23 composites of the calendar
year:

- **Quality.** Terra composites with pixel reliability 0-1 are valid and 2 is snow. VIIRS
  composites are valid with MODLAND_QA ≤ 1 and reliability ≤ 7; reliability 8 is snow.
- **Snow** composites take the pixel's 5th percentile of valid EVI over the sensor's whole
  record, so winter snow cover does not read as vegetation loss.
- **Gaps** of up to two composites are interpolated linearly. A year needs at least 15 valid or
  snow composites.

### Drift correction

Sensor degradation, calibration changes and the global greening trend shift EVI over the years,
also on untouched land. Each annual value $x$ is corrected to $x / e^{\delta}$, with $\delta$ the stored
drift factor of its sensor, year and group: the median log anomaly, in that year, of the
near-pristine control
points of the group: 19,984 equal-area 300 m cells with HMI 2022 below 0.05 and tree, shrub or
grassland cover (ESA WorldCover 2021), with an equal quota per biome. Groups are biome ×
natural cover class; where a group has too few controls the factor falls back to the biome, then
to all controls. The factors for 2001-2025 are frozen. The factor of a new year is computed from
the controls without refitting earlier years, so stored baselines stay valid.

### Change ratio

$$r = \frac{\text{mean corrected EVI of the three-year window}}{\text{median corrected EVI of the baseline}}$$

- **Window**: the three years ending in the reference year, at least two valid.
- **Baseline**: 2009-2018 for Terra and 2012-2018 for VIIRS, at least six valid years. The
  baseline is fixed, so $r$ measures change since 2009-2018 and not against a moving reference.
- **Grid**: $r$ is computed on each sensor's native sinusoidal grid, resampled to the 300 m EII
  grid by area average, and only then combined. A 300 m cell rarely lies inside a single MODIS or
  VIIRS pixel, and nearest-neighbour resampling would change the score by more than the noise
  half-width on about 6 % of cells.
- **Combination**: the unweighted mean of the sensors' $r$ that have a value. The band
  `fi_change_sensors` records how many sensors contributed (0-2).

### Tolerance

The tolerance $h$ is the noise level of $r$ on near-pristine land: a LightGBM quantile model,
fitted on the controls, predicts the 68 % quantile of $|r - 1|$ from the 15 predictors of the
potential NPP model, plus the number of sensors where two sensors serve. It is floored at 0.02.
At the controls, the median $h$ is 0.052 with two sensors (90th percentile 0.088), and 0.065 for
NOAA-20 alone.

### Score

$$\mathrm{FI}_{\text{change}} = \operatorname{clamp}\left(1 - \frac{(1 - r) - k_c h}{1 - k_c h},\ 0,\ 1\right),
\qquad k_c = 1.0$$

FI_change is 1 while the EVI deficit $1 - r$ stays within $k_c h$, and falls linearly to 0 at
$r = 0$. A rise in EVI ($r > 1$) is not penalised. With $h$ = 0.05: $r$ = 0.97 scores 1.00,
$r$ = 0.80 scores 0.84, $r$ = 0.50 scores 0.53 and $r$ = 0.20 scores 0.21.

$k_c$ is the smallest multiplier on a fixed grid for which, at the near-pristine controls and in
every reference year 2021-2025, the 10th percentile of FI stays at or above 0.90, and the share of
controls whose FI crosses 0.9 between years is at most 1.5 times that of FI_state alone. Fitted
on all controls, $k_c$ = 1.0 for every sensor combination. Applied to held-out scores, the 10th
percentile is 0.908-0.916 for Terra with Suomi NPP and 0.903-0.913 for the NOAA-20 combinations.

### Product

$$\mathrm{FI} = \mathrm{FI}_{\text{state}} \times \mathrm{FI}_{\text{change}}, \qquad
\mathrm{FI} = \mathrm{FI}_{\text{state}} \ \text{where no sensor has a value}$$

FI_change is applied whether or not Hansen Global Forest Change registers loss at the pixel, so a
clearing can lower both FI and compositional integrity. Most dated forest disturbance the term
detects is also registered as forest loss, and compositional integrity does not respond to most
of these events, so the term is not masked there.

### Sensor schedule

A sensor serves a reference year only if it delivered all 23 composites in each year of the
three-year window; one missing composite removes the sensor from three reference years. The
VIIRS slot uses Suomi NPP while it serves. NOAA-20 replaces it chained to the Suomi NPP
record: its annual values are multiplied by the pixel's median ratio of Suomi NPP to NOAA-20
over 2018-2020, and divided by the Suomi NPP baseline. The two VIIRS sensors are never mixed
within one window. Each combination has its own stored tolerance model.

| Reference year | Sensors | Code |
|---|---|---|
| 2021-2025 | MODIS Terra and VIIRS Suomi NPP | TV |
| 2026 | MODIS Terra and VIIRS NOAA-20 chained to Suomi NPP (NOAA-20 alone if Terra's 2026 record is incomplete) | TNc |
| 2027 on | VIIRS NOAA-20 chained to Suomi NPP | Nc |
| No sensor with a complete window | none: FI = FI_state everywhere, `fi_change_sensors` = 0 | |

Suomi NPP stopped delivering data on 2 November 2026, and NASA plans to end MODIS Terra science
data collection around January 2027. NOAA-21 is the next VIIRS sensor; it is not part of the
schedule. The series that serves a year and the SHA-256 of the stored constants are written into
the tags of every tile (`fi_change_series`, `fi_change_k`, `fi_change_constants_sha256`).

### Global results, 2021-2025

All reference years 2021-2025 are served by Terra and Suomi NPP. Over land, weighted by area:

| | 2021 | 2022 | 2023 | 2024 | 2025 |
|---|---:|---:|---:|---:|---:|
| Mean FI_change | 0.993 | 0.992 | 0.992 | 0.991 | 0.991 |
| FI_change < 0.9 (share of cells with a value) | 1.6 % | 1.9 % | 2.1 % | 2.2 % | 2.3 % |
| Mean FI_state / mean FI | 0.956 / 0.951 | 0.956 / 0.950 | 0.956 / 0.950 | 0.957 / 0.951 | 0.957 / 0.950 |
| 10th percentile of FI at the controls | 0.922 | 0.924 | 0.925 | 0.923 | 0.917 |

Two sensors serve 85 % of land, one sensor 12 % (mostly VIIRS alone, in cloudy tropical forest,
tundra and mangroves), none 3 %, concentrated in deserts. The change term lowers mean EII by
0.0025-0.0028. Flagged cells concentrate in drylands (Morocco and Algeria, Iran and
Afghanistan, Pakistan and northern India), in fire regions of western North America, in
south-western Australia and in the Chaco. In drylands much of this is likely drought rather than
land use; FI does not separate the two ([Limitations](../limitations.md#climate-and-land-management)).

## Wetland flag

The potential model overestimates natural NPP on wetlands: at natural ICOS wetland sites observed
NPP is about two thirds of the potential. The predictors do not describe the local water supply
that sets wetland productivity, so a low FI on wetlands is partly model error. The `wetland_flag` band marks
these pixels; it never changes FI or EII.

| Code | Family | GLWD v2 classes | Minimum GLWD cover |
|---|---|---|---|
| 0 | none | dryland, open water (1-7), paddy rice (33) | |
| 1 | wetland | lacustrine, riverine, palustrine and ephemeral wetland (8-21), saltmarsh (29), delta (30), other coastal (31), salt pan (32) | 51 % |
| 2 | peatland | arctic/boreal, temperate and tropical peatland (22-27) | 33 % |
| 3 | mangrove | mangrove (28) | 51 % |

A pixel takes its family's code when the dominant GLWD v2 class (Lehner et al. 2025, ~460 m)
belongs to the family and the combined GLWD cover reaches the family minimum. The minimum cover
was chosen from a grid between 1 % and 51 %. Open water is already masked, and paddy rice is
cropland, where a low FI is real. The flag covers about 9 % of land, 16 of 24 ICOS wetland sites
and 2 of 51 ICOS forest sites.

The client reports the flagged share of an area by default; `exclude_wetlands=True` removes
flagged pixels from the statistics.

## Masks

FI is missing where observed NPP has fewer than two valid years in the window, or where potential
NPP is missing (sea, inland water, and pixels without complete predictors). Land is the extent of
the potential NPP map. FI_change is missing where no sensor has a value; FI then equals FI_state.

## What FI measures

- FI_state is a productivity-anomaly score and is near 1 on most land by design: in 2024 it is
  exactly 1 at 83 % of a global land sample.
- With the change term, FI also responds to vegetation loss since the baseline for as long as
  EVI stays below the pixel's baseline beyond the tolerance; the term has no expiry and was
  evaluated on events up to six years old.
- Long-standing conversion is carried by structural and compositional integrity; FI_change does
  not respond to it, and FI_state only where productivity departs strongly from the potential.
- The change term depends on satellite records that are ending: MODIS Terra with VIIRS on Suomi
  NPP for reference years 2021-2025, Terra with VIIRS on NOAA-20 chained to Suomi NPP for 2026,
  and NOAA-20 alone from 2027. For a reference year in which no sensor has complete data for the
  window, FI equals FI_state (band `fi_change_sensors` = 0).

## Parameters

| Parameter | Value | Setting in `eii.compute.settings` |
|---|---|---|
| Observed NPP window | 3 years ending in the reference year, from 2019 | `OBSERVED_NPP_WINDOW_YEARS`, `OBSERVED_NPP_FIRST_YEAR` |
| Minimum valid NPP years | 2 | `OBSERVED_NPP_MIN_VALID_YEARS` |
| Floor | 50 gC m⁻² yr⁻¹ | `POTENTIAL_NPP_FLOOR` |
| Deficit multiplier | 1.0 | `TOLERANCE_K` |
| Surplus multiplier | 1.5 | `TOLERANCE_K_SURPLUS` |
| Surplus zero point | relative surplus 2 (three times $P$) | `SURPLUS_ZERO` |
| Interval | 68 % | |
| Recentring reference | HMI < 0.05 and BII ≥ 0.95 | `RECENTRE_REFERENCE_HMI_MAX`, `RECENTRE_REFERENCE_BII_MIN` |
| Minimum reference pixels per biome | more than 200 | `RECENTRE_MIN_REFERENCE_PIXELS` |
| Change-term baseline | Terra 2009-2018, VIIRS 2012-2018, ≥ 6 valid years | `FI_CHANGE_BASELINES`, `FI_CHANGE_MIN_BASELINE_YEARS` |
| Change-term window | 3 years, ≥ 2 valid | `FI_CHANGE_WINDOW_YEARS`, `FI_CHANGE_MIN_WINDOW_YEARS` |
| Change-term multiplier $k_c$ | 1.0 | `FI_CHANGE_K` |
| Tolerance floor | 0.02 | `FI_CHANGE_H_FLOOR` |
| NOAA-20 chain years | 2018-2020, ≥ 2 valid | `FI_CHANGE_CHAIN_YEARS`, `FI_CHANGE_CHAIN_MIN_YEARS` |
| Composites per year / minimum valid / longest gap filled | 23 / 15 / 2 | `EVI_COMPOSITES_PER_YEAR`, `EVI_MIN_COMPOSITES`, `EVI_MAX_GAP` |
| Wetland flag cover | 51 % / 33 % / 51 % | `WETLAND_FLAG_MIN_COVER_PCT` |

The drift factors, the tolerance models and $k_c$ are stored in the package
(`eii/compute/fi_change_constants/`) and checked against their SHA-256 when loaded.

## Implementation

`eii.compute` implements each step on numpy arrays ([API](../api/compute.md)):

- FI_state: `observed_npp_window`, `observed_npp`, `recentring_sums`, `recentring_factors`,
  `recentre_potential`, `deviation_tolerance_integrity`;
- FI_change: `annual_evi`, `correct_drift`, `chain_factor`, `chained`, `baseline_median`,
  `change_ratio`, `combine_ratios`, `tolerance_h`, `fi_change`, `series_for`;
- the product: `functional_with_change`; the flag: `wetland_flag`.

The global run is in
[`pipelines/production`](https://github.com/landler-io/ecosystem-integrity-index/tree/main/pipelines/production).

```python
import numpy as np
from eii.compute import deviation_tolerance_integrity, fi_change, functional_with_change

fi_state = deviation_tolerance_integrity(
    actual=np.array([300.0, 800.0]),
    potential=np.array(500.0),
    lower=np.array(400.0),
    upper=np.array(620.0),
    natural_cv=np.array(0.1),
)
# array([0.773, 0.87 ], dtype=float32)
change = fi_change(r=np.array([0.8, np.nan]), h=np.array([0.05, 0.05]), k=1.0)
functional_with_change(fi_state, change)
# array([0.651, 0.87 ])
```

## References

- Copernicus Land Monitoring Service. Dry Matter Productivity 300 m, version 2.0.2.
- Didan, K. (2021). MODIS/Terra Vegetation Indices 16-Day L3 Global 250m SIN Grid V061. NASA
  LP DAAC.
- Didan, K., & Barreto, A. (2023). VIIRS/NPP and VIIRS/JPSS1 Vegetation Indices 16-Day L3
  Global 500m SIN Grid V002. NASA LP DAAC.
- Hansen, M. C., et al. (2013). High-resolution global maps of 21st-century forest cover change.
  *Science*, 342, 850-853.
- Lehner, B., et al. (2025). The Global Lakes and Wetlands Database (GLWD) version 2. *Earth
  System Science Data*, 17, 2277-2329.
