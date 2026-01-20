# Contributing to Ecosystem Integrity Index

Thank you for your interest in contributing to the Ecosystem Integrity Index project! This document provides guidelines and information for contributors.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [How to Contribute](#how-to-contribute)
- [Development Setup](#development-setup)
- [Pull Request Process](#pull-request-process)
- [Style Guidelines](#style-guidelines)
- [Testing](#testing)
- [Documentation](#documentation)

## Code of Conduct

This project adheres to a [Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to the maintainers.

## Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/your-username/ecosystem-integrity-index.git
   cd ecosystem-integrity-index
   ```
3. **Set up the development environment** (see below)
4. **Create a branch** for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## How to Contribute

### Reporting Bugs

- Check if the bug has already been reported in [Issues](https://github.com/landler-io/ecosystem-integrity-index/issues)
- If not, create a new issue using the bug report template
- Include as much detail as possible: steps to reproduce, expected behavior, actual behavior, environment details

### Suggesting Features

- Open an issue using the feature request template
- Describe the use case and why this feature would be valuable
- Be open to discussion about alternative approaches

### Contributing Code

- For small fixes, feel free to open a PR directly
- For larger changes, please open an issue first to discuss the approach
- Follow the style guidelines and ensure tests pass

### Contributing Documentation

- Documentation improvements are always welcome
- Follow the existing style and structure
- Test that documentation builds correctly

## Development Setup

### Prerequisites

- Python 3.10 or higher
- Google Earth Engine account (for running tests with GEE)

### Installation

```bash
# Clone the repository
git clone https://github.com/landler-io/ecosystem-integrity-index.git
cd ecosystem-integrity-index

# Create a virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install in development mode with all dependencies
pip install -e ".[all]"

# Install pre-commit hooks
pre-commit install
```

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=eii

# Run specific test file
pytest tests/test_gee/test_eii.py
```

### Building Documentation

```bash
# Serve documentation locally
mkdocs serve

# Build documentation
mkdocs build
```

## Pull Request Process

1. **Update your branch** with the latest main:
   ```bash
   git fetch origin
   git rebase origin/main
   ```

2. **Ensure all checks pass**:
   - Tests: `pytest`
   - Linting: `ruff check .`
   - Formatting: `ruff format .`
   - Type checking: `mypy src/`

3. **Write clear commit messages** following [Conventional Commits](https://www.conventionalcommits.org/):
   - `feat: add new EII extraction method`
   - `fix: correct polygon clipping edge case`
   - `docs: update installation instructions`

4. **Open a Pull Request**:
   - Use a descriptive title
   - Reference related issues
   - Describe what changes were made and why

5. **Respond to review feedback** in a timely manner

## Style Guidelines

### Python Code

- Follow [PEP 8](https://pep8.org/) style guidelines
- Use [ruff](https://github.com/astral-sh/ruff) for linting and formatting
- Maximum line length: 100 characters
- Use type hints for function signatures
- Write docstrings for public functions/classes (Google style)

```python
def get_eii_for_polygon(
    geometry: ee.Geometry,
    year: int,
    scale: int = 300,
) -> ee.Dictionary:
    """
    Extract EII statistics for a given polygon.

    Args:
        geometry: Earth Engine geometry defining the area of interest.
        year: Year for which to extract EII data.
        scale: Resolution in meters for the reduction. Defaults to 300.

    Returns:
        Dictionary containing EII statistics (mean, std, min, max, etc.).

    Raises:
        ValueError: If year is outside the available data range.

    Example:
        >>> polygon = ee.Geometry.Rectangle([-60, -10, -55, -5])
        >>> stats = get_eii_for_polygon(polygon, 2020)
    """
```

### JavaScript (Earth Engine)

- Use consistent indentation (2 spaces)
- Document functions with JSDoc comments
- Use meaningful variable names

### Documentation

- Use Markdown for all documentation
- Include code examples where appropriate
- Keep language clear and accessible

## Testing

### Test Requirements

- All new features should have corresponding tests
- Bug fixes should include a test that reproduces the bug
- Maintain or improve code coverage

### Test Structure

```
tests/
├── conftest.py          # Shared fixtures
├── test_gee/            # GEE-related tests
│   ├── test_eii.py
│   └── test_npp.py
└── test_analysis/       # Analysis module tests
    └── test_zonal_stats.py
```

### Mocking Earth Engine

For tests that don't require actual GEE access, use mocking:

```python
from unittest.mock import MagicMock, patch

@patch('ee.Image')
def test_something(mock_image):
    mock_image.return_value = MagicMock()
    # ... test code
```

## Documentation

### Structure

- `docs/`: Main documentation source
- `notebooks/`: Jupyter notebooks (also rendered in docs)
- Docstrings: API documentation (auto-generated)

### Building Locally

```bash
mkdocs serve
```

Then visit `http://localhost:8000` to preview.

## Questions?

If you have questions about contributing, feel free to:

- Open a [Discussion](https://github.com/landler-io/ecosystem-integrity-index/discussions)
- Ask in an Issue
- Contact the maintainers

Thank you for contributing! 🌍
