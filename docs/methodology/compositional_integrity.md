# Compositional Integrity

Compositional integrity relates to the variety of life itself—the species that make up an ecosystem. A compositionally intact ecosystem is one that retains its native species in their natural abundances. The loss of species, or a significant change in their relative abundances, is a clear indicator of ecological degradation.

## Methodology

To measure compositional integrity, our index utilizes the **Biodiversity Intactness Index (BII)**.

### Biodiversity Intactness Index (BII)
The BII estimates the average abundance of a region’s originally-present species relative to an undisturbed baseline. It is a widely recognized metric in global conservation science (Scholes & Biggs, 2005).

**Data Source:**
We use high-resolution BII data produced by **Impact Observatory** (Gassert et al., 2022). This dataset provides a globally consistent and up-to-date measure of compositional integrity at 300m resolution.

### Implementation
The compositional component effectively ingests the BII data, ensuring it aligns temporally and spatially with the other integrity pillars.

*   **Metric:** Biodiversity Intactness Index (0-1 scale)
*   **Resolution:** 300m

See the [API Reference](../api/compute.md) for `calculate_compositional_integrity`.
