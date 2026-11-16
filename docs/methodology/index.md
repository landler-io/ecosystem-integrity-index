# Methodology Overview

The Ecosystem Integrity Index (EII) measures the condition of terrestrial ecosystems along three
dimensions: **function**, **structure** and **composition**. The framework follows the UNEP-WCMC
Ecosystem Integrity Index of Hill et al. (2026), with its own functional component, an
annual update cycle and a single, reproducible offline calculation.

## Pipeline at a glance

```
CHELSA climate, terrain, coast  ──►  potential NPP model  ──►  potential NPP + 68 % / 90 % intervals
                                     natural CV model     ──►  natural inter-annual CV
CLMS DMP V2.0.2 (3-year window) ─────────────────────────►  observed NPP
                                                               │
                         recentring per biome and window ◄─────┤
                                                               ▼
                                                            FI_state
MODIS Terra + VIIRS EVI (3-year window vs 2009-2018) ────►  FI_change
                                                               ▼
                                         functional integrity  F = FI_state × FI_change
GHM 2022 (q(HMI), 26.5 km kernel)  ─────────────────►  structural integrity  S
Biodiversity Intactness v1.1 (same year)  ───────────►  compositional integrity C
                                                               │
                                min × fuzzy sum  ◄─────────────┘
                                       ▼
                    EII, limiting component, wetland flag  (per reference year)
```

## The three components

### Functional integrity

*Does the land produce as much as natural vegetation would in the same place, and has its
vegetation recently been lost?*

[Functional integrity](functional_integrity.md) is the product of two scores. FI_state compares
observed net primary productivity (NPP) from Copernicus satellite data with the NPP that natural
vegetation would reach under the local climate and terrain. The
[potential NPP model](potential_npp.md) is trained on the world's least modified land.
Deviations in either direction count: a deficit points to degradation or disturbance, a surplus
to irrigation, fertilisation or other inputs. Deviations within the uncertainty of the model and
the natural year-to-year variability score 1. FI_change compares the pixel's vegetation
greenness (EVI from MODIS Terra and VIIRS) over the last three years with its own 2009-2018
baseline and lowers the score where vegetation has been lost.

### Structural integrity

*How intact is the land and the landscape around it?*

[Structural integrity](structural_integrity.md) scores each pixel by its own habitat quality and
by the quality of the surrounding landscape within about 25 km, following Hill et al. (2026) and
Beyer et al. (2020). Quality is a logistic curve in HMI, calibrated against WorldCover natural
land cover: lightly used land keeps nearly full quality, converted land little. Fully modified
land scores 0, unmodified land inside an unmodified landscape scores 1. For natural land, the score is
roughly the share of natural land around it ([in plain terms](structural_integrity.md#in-plain-terms)).

### Compositional integrity

*How much of the original species community remains?*

[Compositional integrity](compositional_integrity.md) is the Biodiversity Intactness Index (BII)
of the reference year: the modelled average abundance of originally present species relative to
an undisturbed baseline.

## Aggregation

The three components are combined with a [limiting-factor rule](ecosystem_integrity.md): the
lowest component sets the ceiling, and the other two lower the score further when they are also
degraded. Each pixel also records which component limits it.

## Annual cycle

EII is published as one map per reference year, 2021 to 2025. What changes between years:

| Input | Per reference year |
|---|---|
| Observed NPP | Mean of the three annual sums ending in the reference year |
| EVI (change term) | Mean of the three years ending in the reference year, against the fixed 2009-2018 baseline |
| Potential NPP, its intervals and natural CV | Static (one model, one map) |
| Recentring factors | Recomputed for each observed window |
| Human Modification Index | 2022 for every year |
| Biodiversity Intactness | Same year as the reference year |

## Further reading

- [Potential NPP model](potential_npp.md): training data, predictors, model and intervals.
- [Validation](validation.md): model skill and agreement with independent indicators.
- [Local Modulation](local_modulation.md): optional plot-level adjustment with Natural Capital KPIs.
- [Limitations](../limitations.md): interpretation guidance.

## References

- Beyer, H. L., Venter, O., Grantham, H. S., & Watson, J. E. M. (2020). Substantial losses in
  ecoregion intactness highlight urgency of globally coordinated action. *Conservation Letters*,
  13, e12692.
- Hill, S. L. L., et al. (2026). The ecosystem integrity index: a component indicator in the
  Kunming-Montreal global biodiversity framework reveals the integrity of terrestrial ecosystems
  globally. *Ecological Indicators*, 190, 115314.
