# Ecosystem Integrity Index (EII)

The **Ecosystem Integrity Index (EII)** synthesizes the three pillars—Function, Structure, and Composition—into a single, scientifically robust score.

## The Aggregation Challenge
An ecosystem can be degraded in multiple ways, and a holistic assessment must account for all of them. Simple averaging often masks critical failures. For example, a forest plantation might have high productivity (Functional) and high connectivity (Structural), but near-zero biodiversity (Compositional). A simple average would suggest a "healthy" ecosystem, while in reality, it is ecologically poor.

## Methodology: Minimum with Fuzzy Logic

To address this, the EII employs a **Limiting Factor** approach combined with **Fuzzy Logic**.

### 1. Limiting Factor Principle
Following ecological theory (Liebig's Law of the Minimum), we recognize that an ecosystem's overall integrity is constrained by its most degraded component. Therefore, the EII score is primarily driven by the **minimum** of the three pillar scores.

### 2. Fuzzy Interaction
However, simply taking the minimum ignores the state of the other components. If two ecosystems both have a minimum score of 0.5, but one has high scores in the other pillars while the other has low scores, they should not receive the exact same final rating.

We use a fuzzy logic aggregation that starts with the minimum but modulates it based on the other factors:

$$ EII = M \times F $$

Where:
*   $M = \min(\text{Functional}, \text{Structural}, \text{Compositional})$
*   $F = \text{Fuzzy Sum of the remaining two components}$

The **Fuzzy Sum** of two values $A$ and $B$ is calculated as:
$$ \text{FuzzySum}(A, B) = A + B - (A \times B) $$

### Summary
This formula ensures that:
1.  **Strict Ceiling:** The EII can never exceed the score of the lowest performing pillar.
2.  **Penalty for Cumulative Degradation:** If multiple pillars are degraded, the score is lower than if only one pillar is degraded.

See the [API Reference](../api/compute.md) for `combine_components` and `calculate_eii`.
