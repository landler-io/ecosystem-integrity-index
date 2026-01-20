# Ecosystem Integrity Index (EII)

[![Documentation](https://img.shields.io/badge/docs-online-blue)](https://landler-io.github.io/ecosystem-integrity-index/)
[![License](https://img.shields.io/badge/license-Apache%202.0-green)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)


## Overview

The Ecosystem Integrity Index (EII) is a holistic metric designed to quantify the health of terrestrial ecosystems globally. Unlike simpler metrics that rely on a single proxy, EII synthesizes three fundamental dimensions of ecosystem integrity: **Function**, **Structure**, and **Composition**.

This implementation builds and expands on the conceptual framework developed by the UNEP-WCMC: [Hill et al. (2022) The Ecosystem Integrity Index: a novel measure of terrestrial ecosystem integrity](https://www.biorxiv.org/content/10.1101/2022.08.21.504707v2.full).


![EII Quicklooks](notebooks/figures/quicklook_eii_stats.png)


## Quick Start

```python
import ee
from eii.client import get_stats

# Initialize Earth Engine
ee.Authenticate()
ee.Initialize(project='your-project')

# Define area of interest
polygon = ee.Geometry.Rectangle([-60, -10, -55, -5])

# Get EII data
eii_data = get_stats(polygon)
print(eii_data)
```

## Installation

For development:

```bash
git clone https://github.com/landler-io/ecosystem-integrity-index.git
cd ecosystem-integrity-index
pip install -e ".[dev]"
```

## Documentation

Full documentation is available [here](https://landler-io.github.io/ecosystem-integrity-index)

- [Getting Started](docs/getting-started.md)
- [Methodology](docs/methodology/)
- [Limitations & Roadmap](docs/limitations.md)
- [API Reference](docs/api/)


## Repository Structure

```
ecosystem-integrity-index/
├── src/eii/              # Python package for data access and compute
├── pipelines/            # Data preprocessing, modeling, and exports
├── notebooks/            # Placeholder for example notebooks
├── docs/                 # Documentation
└── analysis/             # GEE web application (future)
```

## Citation

If you use EII in your research, please cite:

> Leutner, B. (2025) "The Ecosystem Integrity Index: A comprehensive, globally applicable, and unified metric for ecosystem health". Technical Whitepaper. The Landbanking Group, Munich, Germany. [PDF](https://media.thelandbankinggroup.com/ecosystem-integrity-index-eii.pdf)

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the Apache License 2.0 - see [LICENSE](LICENSE) for details.
