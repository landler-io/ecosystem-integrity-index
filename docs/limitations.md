# Limitations and Interpretation

EII condenses three dimensions of ecosystem condition into one globally consistent number. That
requires choices and data that shape how the index should be read. This page sets out what the
index measures, where it is weak, and how to use it for different purposes.

## Functional integrity

### A contemporary reference

Potential NPP is learned from today's least modified land, observed in 2019 to 2025. It describes
what natural vegetation produces under present conditions, including elevated CO₂, nitrogen
deposition and recent climate, not a pre-industrial baseline. An ecosystem that matches its
present-day natural potential scores 1 even if it once produced more.

### Extrapolation into transformed regions

Natural reference land is unevenly distributed. In regions with a long history of intensive use,
such as Western Europe, South Asia or eastern China, few reference sites remain in the productive
lowlands, and the model learns their potential partly from similar climates elsewhere. Held out as
a whole, Europe is predicted with R² of 0.61 to 0.71, against 0.9 globally. Read functional
integrity in such regions as an estimate informed by global relationships, not a locally
calibrated benchmark.

Where no natural semi-arid land is left at all (Sahel, Loess Plateau, semi-arid India), the
training set is topped up with lightly used land (HMI < 0.2) at half weight. There, potential NPP
partly describes lightly used rather than natural land.

### Known model biases

- Potential NPP is low in Mediterranean and tropical coniferous biomes and high in deserts (median
  deviations of up to about 12 % per biome at held-out reference sites under block CV, up to
  about 24 % when a whole realm is held out). The tolerance and recentring absorb part of this.
- Wetlands and peatlands are overestimated, so their functional integrity tends to be too low.
  The `wetland_flag` band marks them; use `exclude_wetlands=True` to leave them out of statistics.
- On flat lowland the potential is a few percent lower than on adjacent moderate terrain, which
  understates cropland losses on valley floors.
- See the [model page](methodology/potential_npp.md#known-weaknesses) for details.

### What functional integrity can detect

The state score (FI_state) is 1 within a tolerance that covers model uncertainty and natural
year-to-year variability; smaller deviations are not penalised. It responds to strong and lasting
changes in productivity: conversion with lower or higher NPP, irrigation in drylands, severe
degradation. A surplus is penalised less steeply than a deficit: where potential NPP is at least
50 gC m⁻² yr⁻¹, FI_state reaches 0 at total loss, but only at three times the potential for a
surplus. Below that floor deviations are measured in absolute terms (units of 50 gC), so a total
loss of very low potential NPP can still score 1. Irrigated agriculture in true deserts,
whose NPP is many times the 50 gC floor, still scores close to 0. Annual NPP does not separate
many conversions in which crops and natural vegetation produce similar amounts.

The change term (FI_change) responds to vegetation loss since the pixel's 2009-2018 baseline:
clearing, fire, severe dieback. Its tolerance is set so that near-pristine land keeps FI of at
least 0.9 in nine of ten cases, so small or partial losses within the noise of the EVI record
are not penalised, and a loss is no longer seen once the vegetation has regrown to its baseline
greenness. Land converted before the baseline and kept in the same use scores FI_change = 1.

### Near 1 on most land, by design

- Functional integrity is a productivity-anomaly measure. It is near 1 on most land by design:
  in 2024 the observed-against-potential score is exactly 1 at 83 % of a global land sample, and
  mean FI over land is 0.95.
- With the change term, FI = FI_state × FI_change, FI also responds to vegetation loss for as
  long as EVI stays below the pixel's baseline; it was evaluated on events up to six years old.
  FI_change is below 0.9 on about 2 % of land and lowers mean EII by less than 0.003.
- Long-standing conversion is carried by structural and compositional integrity. FI_change does
  not respond to it, FI_state only where productivity departs strongly from the potential, and a
  low FI on cropland should not be expected.
- FI near 1 almost everywhere is a property of the measure, not a gap in the data: annual NPP of
  long-used land is often similar to that of natural land, so a score that keeps intact land near
  1 cannot also separate the two by productivity.

### Satellite records behind the change term

The change term depends on satellite records that are ending: MODIS Terra with VIIRS on Suomi
NPP for reference years 2021-2025, Terra with VIIRS on NOAA-20 chained to Suomi NPP for 2026, and
NOAA-20 alone from 2027. Suomi NPP stopped delivering data on 2 November 2026 and MODIS Terra is
planned to stop around January 2027. NOAA-20 alone is noisier than two sensors, so FI_change
flags more cells at the same tolerance. For a reference year in which no sensor has complete data
for the window, FI equals FI_state (band `fi_change_sensors` = 0). The serving sensors are
recorded per tile; compare FI across years with that in mind.

On 3 % of land, mostly deserts, no sensor gives a change value, and FI equals FI_state. One
sensor (mostly VIIRS) serves 12 % of land, more in cloudy tropical forest, tundra and mangroves.

### Overlap with compositional integrity

FI_change is not masked where forest loss is mapped. Compositional integrity can respond to the
same clearing through the land-use maps behind Biodiversity Intactness, so a single recent event
can lower two components. The aggregation then lowers EII more than for a loss in one component.

### Special cases

- **Irrigated drylands score close to 0, by design.** Irrigated cropland in a desert (Sindh, the
  Nile valley and delta, canal-irrigated land in Egypt and Pakistan) produces many times the
  potential of the surrounding desert. Functional integrity measures departure from natural
  function, and the land has been replaced by a system of very different productivity, so F
  near 0 means "replaced", not "degraded". This agrees with HANPP accounting, which also treats
  NPP above the potential as human alteration.
- **River floodplains are under-predicted.** The potential model sees local climate and
  terrain but not water from upstream (upstream area is used only within 50 km of the coast).
  The natural state of large dryland floodplains (Nile, Indus) was wetland and riverine
  forest, far more productive than the climate predicts. There the surplus, and the low F, are
  partly a model limitation rather than a measured change.
- **Alpine meadows can read as surplus.** Productive meadows in the Pamir, Karakoram and on the
  Tibetan Plateau exceed a potential learned from the mostly sparse plateau, although they are
  near-pristine. The 1 km climate predictors cannot see the local water that makes them
  productive.
- **Bare rock and ice can read as deficit.** Where the potential model predicts some NPP on rock
  or glaciers (Karakoram, Tibet), observed NPP of 0 scores as a deficit.

### Biome-wide change is absorbed

Potential NPP is recentred for each three-year window so that intact land of each biome sits at
its potential. This removes drift in the satellite record and regional weather from the score. It
also means that change uniform across a whole biome within one window, such as a continent-wide
drought, does not lower functional integrity: the score measures function relative to intact land
of the same biome in the same years.

### Climate and land management

Functional integrity does not separate productivity changes caused by climate from those caused
by land management. Recentring absorbs biome-wide change in FI_state, but FI_change compares a
pixel with its own baseline and has no such correction beyond the global drift factors. Its
flags concentrate in drylands (North Africa, Iran and Afghanistan, Pakistan and northern India)
and in fire regions, where much of the loss is likely drought or fire rather than land use. Where
regional climate shifts persistently, check whether climate data support a climate-driven
explanation before reading a decline as local degradation.

## Structural integrity

- **Pressure, not habitat.** HMI measures human pressure. Low-pressure land that is not natural
  habitat, such as extensive pasture, counts by its low pressure, and small natural features
  inside high-HMI pixels are not seen. Unmapped degradation, such as selective logging and forest
  edges in the moist tropics, leaves structural integrity high.
- **Calibrated against a bound.** The quality curve was chosen so that structural integrity does
  not exceed a WorldCover land-cover bound, which itself overstates habitat (pasture and
  plantations count as natural). The midpoint (HMI 0.4) is a judgement within that bound.
- **Amount, not connectivity.** The landscape term measures how much good land lies near a pixel,
  not whether organisms can move between patches. Connectivity analyses need species-specific
  dispersal and matrix permeability.
- **Static input.** HMI 2022 is used for every reference year, so structural integrity does not
  register change after 2022.

## Compositional integrity

Biodiversity intactness is modelled from the PREDICTS database of biodiversity surveys. Surveys
concentrate in Western Europe, North America and a few tropical research areas; elsewhere the
model extrapolates. BII responds to mapped pressures (land conversion, intensity, roads,
population) and cannot see unmapped ones such as hunting, invasive species or species
reintroduction. It models the community expected under current pressure, so it can anticipate
losses or gains that have not yet happened on the ground (extinction debt, colonisation credit).

## Shared drivers of structure and composition

Structural integrity is computed from human modification, and biodiversity intactness is modelled
from similar pressure data. The two components are therefore correlated, and the effective
dimensionality of EII is closer to two than three: a pressure signal (structure and composition)
and a productivity signal (function). Functional integrity correlates only weakly with structure
and limits EII mostly on otherwise intact land.

## Annual change

Year-to-year change in EII comes from observed NPP, EVI and biodiversity intactness; structural
integrity is fixed at 2022.

- The three-year windows of NPP and EVI mean each year shares two years with its neighbours.
  Change between consecutive years is damped and lags real change by up to two years.
- Small year-to-year changes in functional integrity at a pixel can come from weather within the
  window and from the recentring factors. Trends over several years, and changes large enough to
  cross the tolerance, are more reliable.
- Compositional change follows the annual land-use inputs of the BII model.

## Resolution and boundaries

The grid spacing is 1/336°, so a pixel covers about 11 cos(latitude) ha (11 ha at the equator,
7.8 ha at 45°). The effective resolution is coarser: the climate predictors are about 1 km, the
change term uses 250-500 m MODIS and VIIRS EVI, and structural integrity averages over a
neighbourhood of several kilometres. Single pixels should not be read on their own. Roads, small
clearings and narrow riparian strips fall below the pixel size, and pixels at sharp boundaries (field edges, coastlines, urban fringes) mix both sides.
Structural integrity is computed at about 300 m and potential NPP depends on 1 km climate layers,
so its fine detail comes from terrain and observed NPP. Aggregating to a property, management
unit or region is more robust than reading single pixels.

## Aggregation

The minimum with fuzzy sum treats every dimension as necessary: a severe deficit in one cannot be
compensated by the others. This suits an overall condition assessment. For questions centred on
one function, such as productivity or carbon uptake, also read the component layers and the
`limiting_component` band.

## Guidance by application

- **National to continental reporting, portfolio screening.** EII gives consistent comparisons
  across diverse landscapes; pixel errors average out over large areas.
- **Regional prioritisation and planning.** Use EII as a relative ranking within a region. The
  limiting component points at the dimension to address.
- **Site assessment and project monitoring.** Use EII as one input among several, together with
  local data and expertise. [Local modulation](methodology/local_modulation.md) adds plot-level
  Natural Capital indicators.
- **Restoration monitoring.** Track components separately. Functional integrity can respond
  within a few years of recovery: FI_change returns to 1 once EVI regains its baseline level.
  Compositional integrity follows the mapped pressures that drive BII; structural integrity
  follows only once a new human-modification map is available.

## References

- Hudson, L. N., et al. (2017). The database of the PREDICTS (Projecting Responses of Ecological
  Diversity In Changing Terrestrial Systems) project. *Ecology and Evolution*, 7, 145-188.
- Scholes, R. J., & Biggs, R. (2005). A biodiversity intactness index. *Nature*, 434, 45-49.
