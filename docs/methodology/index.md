# Methodology Overview

The Ecosystem Integrity Index (EII) is a holistic metric designed to quantify the health of terrestrial ecosystems globally. Unlike simpler metrics that rely on a single proxy, EII synthesizes three fundamental dimensions of ecosystem integrity: **Function**, **Structure**, and **Composition**.

This implementation builds and expands directly on the conceptual framework developed by the UNEP-WCMC: [Hill et al. (2022) The Ecosystem Integrity Index: a novel measure of terrestrial ecosystem integrity](https://www.biorxiv.org/content/10.1101/2022.08.21.504707v2.full).



## The Three Pillars of Integrity

Our framework defines ecosystem integrity through three distinct but interconnected pillars.

### 1. Functional Integrity
*Does the ecosystem function at its full potential?*

[**Functional Integrity**](functional_integrity.md) assesses the energy capture and productivity of the ecosystem. In our implementation, this is calculated by measuring the **deviation in Net Primary Productivity (NPP)**. Specifically, we compare the observed NPP against a modeled natural potential NPP for the given environmental conditions. A high score indicates that the ecosystem is functioning (photosynthesizing) at the level expected of a healthy, natural system.

### 2. Structural Integrity
*Is the physical habitat intact and connected?*

[**Structural Integrity**](structural_integrity.md) evaluates the physical landscape configuration. It measures the degree of human modification (infrastructure, agriculture, urbanization) and the connectivity of natural habitats. High structural integrity means large, unfragmented areas of natural vegetation where ecological flows can occur unimpeded.

### 3. Compositional Integrity
*Is the native biodiversity present?*

[**Compositional Integrity**](compositional_integrity.md) measures the variety of life within the ecosystem. We utilize the **Biodiversity Intactness Index (BII)** to estimate how much of the original biotic community (species abundance and diversity) remains compared to a pristine baseline.

## Aggregation Strategy

To produce the final EII score, these three components are combined using a **Limiting Factor** approach.

Instead of a simple average—which could mask a collapsed component (e.g., a highly productive monoculture has high function but low composition)—we use a logical framework where the final score is constrained by the lowest performing pillar. We further refine this with a fuzzy logic adjustment to account for the cumulative weight of degradation across multiple dimensions.

See [**Ecosystem Integrity Aggregation**](ecosystem_integrity.md) for the mathematical details.

## Local Modulation (Natural Capital)

An optional **Local Modulation** adjusts the base EII using plot-level Natural Capital (biodiversity, soil, water). It shifts EII by up to ±0.05 based on site-specific KPIs. See [**Local Modulation**](local_modulation.md) for the client API, KPI definitions, and data sources.
