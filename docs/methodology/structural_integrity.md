# Structural Integrity

Structural integrity refers to the physical intactness of an ecosystem. This includes the size, shape, and connectivity of natural habitats. A landscape with high structural integrity is characterized by large, connected areas of natural vegetation, free from significant human modification.

When ecosystems become fragmented by roads, agriculture, or urban development, their structural integrity is compromised. This can isolate populations of species, disrupt ecological flows, and make the ecosystem more vulnerable to external pressures.

## Methodology

We assess structural integrity using a **quality-weighted core area** approach that captures both landscape configuration (fragmentation) and habitat quality.

### 1. Data Source: Human Modification Index

The analysis builds on the **Global Human Modification (GHM)** dataset (Theobald et al., 2020). This dataset integrates various human pressures:

- Built-up areas and urban centers
- Agricultural land use
- Infrastructure (roads, railways)
- Energy and mining activities

The HMI ranges from 0 (no modification) to 1 (complete modification).

### 2. Core Area Calculation

Unlike simple density metrics, our approach identifies **core habitat**—interior areas unaffected by edge effects. This is ecologically important because habitat edges experience altered microclimate, increased predation, and invasion by edge-adapted species.

**Process:**

1. **Habitat mask**: Pixels with HMI < 0.4 are classified as natural/semi-natural habitat
2. **Erosion**: The habitat mask is eroded by 300m (the edge depth) using a morphological minimum filter
3. **Core identification**: Pixels that survive erosion are "core" habitat; those removed are "edge" habitat

This approach is configuration-sensitive: a single large patch retains most of its area as core, while many small fragments of equal total area may have zero core (all edge).

### 3. Quality Weighting

Core pixels are weighted by their habitat quality class based on HMI:

| HMI Range | Quality Class | Weight | Description |
|-----------|--------------|--------|-------------|
| < 0.1 | Pristine | 4 | Minimal human modification |
| 0.1 - 0.2 | Low-impact | 3 | Light modification (e.g., extensive grazing) |
| 0.2 - 0.3 | Moderate | 2 | Moderate modification |
| 0.3 - 0.4 | Semi-natural | 1 | Modified but retains natural character |
| ≥ 0.4 | Modified | 0 | Not counted as habitat |

### 4. Neighborhood Aggregation

The weighted core values are averaged within a **5 km radius** neighborhood, then normalized by the maximum weight (4) to produce a 0-1 score.

**Score interpretation:**

| Score | Interpretation |
|-------|----------------|
| 1.0 | Landscape is entirely pristine core habitat |
| 0.75 | Landscape is entirely low-impact core habitat |
| 0.5 | Landscape is entirely moderate-quality core habitat |
| 0.25 | Landscape is entirely semi-natural core habitat |
| 0.0 | No core habitat (highly fragmented or modified) |

### 5. Parameters

| Parameter | Value | Justification |
|-----------|-------|---------------|
| **Neighborhood radius** | 5 km | Within 90th percentile of species' scale of effect (Miguet et al., 2016); balances farm-level actionability with landscape context |
| **HMI threshold** | 0.4 | Upper bound of Kennedy et al.'s (2019) "moderate modification" class; includes semi-natural habitat |
| **Edge depth** | 300 m | Supported by meta-analyses; 80% of edge effects occur within 300m (Ries et al., 2004) |

## Why Core Area?

Traditional approaches using smoothed habitat quality (e.g., Gaussian convolution of inverted HMI) capture habitat **extent** but not **configuration**. Two landscapes with identical habitat amount but different fragmentation patterns would receive the same score:

```
Landscape A: One large patch (40% habitat)
████████░░░░░░░░░░   → Smoothed score = 0.4
████████░░░░░░░░░░   → Core area score = HIGH (mostly interior)

Landscape B: Many small patches (40% habitat)
█░█░█░█░█░█░█░█░█░   → Smoothed score = 0.4 (same!)
░█░█░█░█░█░█░█░█░█   → Core area score = LOW (all edge)
```

Core area captures this distinction, making it a more meaningful measure of structural integrity.

## Limitations

Our structural metric measures habitat fragmentation and patch size distribution but does not directly measure **functional connectivity** (the ability of organisms to move between patches).

Core area captures configuration (one large vs. many small patches) but not connectivity (whether patches are within dispersal distance). For example, two large patches separated by 10 km would score high on core area but may be functionally isolated for low-mobility species.

True connectivity metrics require species-specific dispersal parameters and are computationally intensive at global scale. Core area provides a reasonable proxy by recognizing that larger, less fragmented patches generally support better population viability and ecological functioning.

## References

- Kennedy, C. M., et al. (2019). Managing the middle: A shift in conservation priorities based on the global human modification gradient. *Global Change Biology*, 25(3), 811-826.
- Miguet, P., et al. (2016). What determines the spatial extent of landscape effects on species? *Landscape Ecology*, 31(6), 1177-1194.
- Ries, L., et al. (2004). Ecological responses to habitat edges: mechanisms, models, and variability explained. *Annual Review of Ecology, Evolution, and Systematics*, 35, 491-522.
- Theobald, D. M., et al. (2020). Earth transformed: detailed mapping of global human modification from 1990 to 2017. *Earth System Science Data*, 12, 1953-1972.

See the [API Reference](../api/compute.md) for `calculate_structural_integrity`.
