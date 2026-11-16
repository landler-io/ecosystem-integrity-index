# Compositional Integrity

Compositional integrity describes how much of an ecosystem's original species community remains.
A compositionally intact ecosystem keeps its native species at their natural abundances; the loss
of species, or a shift in their relative abundances, signals degradation.

## Data: Biodiversity Intactness Index

The Biodiversity Intactness Index (BII; Scholes & Biggs 2005) estimates the average abundance of
a region's originally present species relative to an undisturbed baseline. EII uses
**Biodiversity Intactness v1.1** (Gassert, Mazzarello & Hyde 2026; Impact Observatory and
Vizzuality; CC BY 4.0): annual maps from 2017 to 2025 at 100 m. The underlying models relate
species abundance in the PREDICTS database to land use, land-use intensity and other mapped
pressures.

Each reference year uses BII of the same year, so compositional integrity follows annual changes
in mapped land use.

## Method

1. BII of the reference year is averaged from 100 m onto the CLMS NPP grid. The source has no
   nodata value: sea, inland water and gaps are stored as 0. These zeros are treated as missing
   before averaging, so coastlines and lake shores are not pulled towards 0.
2. Values are clamped to [0, 1].
3. Pixels outside land (the extent of the potential NPP map) are masked.

$$C = \operatorname{clamp}(\mathrm{BII}, 0, 1)$$

Missing BII leaves compositional integrity, and therefore EII and the limiting component, masked
at that pixel.

`eii.compute.compositional_integrity` implements steps 2 and 3 ([API](../api/compute.md)).

## Interpretation

- BII is a modelled estimate, not a survey. It responds to mapped pressures (land conversion,
  intensity, roads, population) and cannot see unmapped ones such as hunting, invasive species or
  local restoration actions.
- PREDICTS surveys are concentrated in Europe, North America and a few tropical research areas;
  predictions elsewhere extrapolate from those relationships.
- Because BII is driven by pressure maps, compositional and structural integrity are correlated.
  See [Limitations](../limitations.md#shared-drivers-of-structure-and-composition).

## References

- Gassert, F., Mazzarello, J., & Hyde, S. (2026). Biodiversity Intactness v1.1. Impact Observatory
  and Vizzuality.
- Hudson, L. N., et al. (2017). The database of the PREDICTS (Projecting Responses of Ecological
  Diversity In Changing Terrestrial Systems) project. *Ecology and Evolution*, 7, 145-188.
- Scholes, R. J., & Biggs, R. (2005). A biodiversity intactness index. *Nature*, 434, 45-49.
