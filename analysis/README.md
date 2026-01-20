# Analysis

Reproducible research outputs and regional assessments.

This directory contains analysis scripts and results for specific studies,
separate from the reusable code in `src/eii/`.

## Structure

```
analysis/
├── figures/                 # Generated figures for publications
└── interactive/             # Ad-hoc exploratory notebooks
```

## Reproducibility

Each analysis subdirectory should contain:

- Analysis notebook or script
- `README.md` explaining the analysis
- `results/` folder with outputs (if not too large)

Large result files should be stored externally and linked in the README.
