# Ecosystem Integrity Index

The EII combines functional (F), structural (S) and compositional (C) integrity into one score.

## Why not an average

An ecosystem can be degraded in different ways, and an average hides a collapsed component. A
plantation may be highly productive (high F) and sit in a forested landscape (high S) while
holding little of the native species community (low C). Its average would suggest a healthy
ecosystem.

## Minimum with fuzzy sum

The EII follows a limiting-factor rule: the most degraded component sets the ceiling. The other
two components then lower the score further when they are degraded too. With the components
sorted as $x_{(1)} \le x_{(2)} \le x_{(3)}$:

$$\mathrm{EII} = x_{(1)} \cdot \left(x_{(2)} + x_{(3)} - x_{(2)}\,x_{(3)}\right)$$

The second factor is the fuzzy (probabilistic) sum of the other two components. It is 1 when
either of them is 1 and falls only when both are degraded.

| F | S | C | EII |
|---|---|---|---|
| 1.0 | 1.0 | 1.0 | 1.00 |
| 1.0 | 0.5 | 1.0 | 0.50 |
| 0.8 | 0.5 | 0.8 | 0.48 |
| 0.8 | 0.5 | 0.5 | 0.45 |
| 0.4 | 0.4 | 0.4 | 0.26 |

The rule guarantees that:

1. EII never exceeds the lowest component;
2. EII equals the lowest component when either of the other two is fully intact, and falls below
   it when both are degraded, so degradation in several components weighs more than in one.

## Missing components

EII requires all three components. Where any of them is missing, EII and the limiting component
are masked; a missing component is never treated as 0 or skipped. Each component band keeps its
own mask.

## Limiting component

The `limiting_component` band records the smallest component per pixel:

| Code | Component |
|---|---|
| 0 | none (a component is missing) |
| 1 | Functional |
| 2 | Structural |
| 3 | Compositional |

Ties prefer functional, then structural, then compositional. The band shows what drives a low
score: modified landscapes are usually limited by structure or composition, while functional
integrity mostly limits otherwise intact land where productivity departs from its potential.

## Alternative aggregations

`eii.compute.combine_components` also offers `minimum`, `product` and `geometric_mean`, and
`eii.client.compare_methods` recombines the published components with all four methods for an area.
The published EII uses `min_fuzzy_logic`.

```python
import numpy as np
from eii.compute import combine_components, limiting_component

f, s, c = np.array([1.0, 0.8]), np.array([0.5, 0.5]), np.array([1.0, 0.5])
combine_components(f, s, c)   # array([0.5 , 0.45], dtype=float32)
limiting_component(f, s, c)   # array([2, 2], dtype=uint8)
```

See the [API Reference](../api/compute.md) for `combine_components` and `limiting_component`.
