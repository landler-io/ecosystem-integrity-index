# Functional Integrity

Functional integrity assesses whether an ecosystem is operating at its full, natural potential. The primary engine of most ecosystems is photosynthesis, the process by which plants capture solar energy and convert it into organic matter. The rate at which this occurs is called **Net Primary Productivity (NPP)**.

A healthy, well-functioning ecosystem will have an NPP close to its natural potential, given the local climate, soil, and topography. A significant deviation from this potential suggests that the ecosystem's functions are impaired, for example due to land degradation, pollution, or unsustainable management practices.

## Methodology

To assess functional integrity, we compare an ecosystem’s **Actual NPP** to its **Potential Natural NPP**.

### 1. Actual NPP
We utilize data derived from state-of-the-art satellite constellations (i.e., Sentinel-3, Proba-V) that continuously observe the Earth's vegetation. This provides an estimate of the current NPP across the globe, every 10 days, at 300m spatial resolution.

### 2. Potential Natural NPP
Potential natural NPP is the NPP we would expect in an ecosystem under minimal human influence. We model this using a machine learning approach trained on the world's most pristine and protected natural areas.

**Predictor Variables:**
The model establishes the statistical relationship between NPP and key environmental factors:

*   **Climate:** WorldClim Bioclimatic variables (Temperature, Precipitation, Seasonality).
*   **Water Availability:** Global Aridity Index.
*   **Topography:** Elevation, Slope, TPI, TRI, CTI (derived from MERIT DEM and Geomorpho90m).
*   **Soil:** Sand and Clay content (SoilGrids).

### 3. Calculating the Score
Functional integrity combines a magnitude component (how far actual NPP deviates from potential) and a seasonality component (how unusual observed intra-annual variability is compared with natural variability).

**Magnitude integrity (0–1):**
1. Compute proportional deviation (symmetric around 1):
   - $$\mathrm{relative\_npp} = \frac{\mathrm{actual\_npp}}{\mathrm{potential\_npp}}$$
   - $$\mathrm{proportional\_score} = \frac{1}{1 + \left|\mathrm{relative\_npp} - 1\right|}$$
2. Compute absolute deviation:
   - $$\mathrm{npp\_difference} = \left|\mathrm{actual\_npp} - \mathrm{potential\_npp}\right|$$
   - Apply a **truncated linear scaling** using a global percentile ceiling (default p95):
     - $$\mathrm{absolute\_score} = 1 - \frac{\operatorname{clamp}(\mathrm{npp\_difference}, 0, p95)}{p95}$$
3. Smooth **derived layers** ($\mathrm{relative\_npp}$, $\mathrm{npp\_difference}$) with a 3×3 window to reduce hard-edge artifacts before scoring.
4. Combine: $$\mathrm{magnitude\_integrity} = \frac{\mathrm{proportional\_score} + \mathrm{absolute\_score}}{2}$$

**Seasonality integrity (0–1):**
1. Compare observed vs natural intra-annual std:
   - $$\mathrm{std\_ratio} = \frac{\mathrm{observed\_std}}{\mathrm{natural\_std}}$$
   - $$\mathrm{seasonality\_integrity} = \frac{1}{1 + \left|\mathrm{std\_ratio} - 1\right|}$$

**Final functional integrity:**
- $$\mathrm{functional\_integrity} = \frac{2}{3}\,\mathrm{magnitude\_integrity} + \frac{1}{3}\,\mathrm{seasonality\_integrity}$$

Scores are normalized to 0–1, where 1 represents functioning at natural potential.

See the [API Reference](../api/compute.md) for `calculate_functional_integrity`.
