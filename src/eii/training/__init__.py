"""
NPP model training utilities.

This module provides functions for training and validating the NPP potential model,
including sampling from natural areas and model evaluation.

Requires additional dependencies: pip install eii[training]

Note:
    Most users will not need this module. Use eii.client for retrieving
    pre-computed EII data, or eii.compute for calculating EII with the
    existing trained model.

Example:
    >>> from eii.training import setup_training_grid, train_npp_model
    >>> grid_cells = setup_training_grid()
    >>> model = train_npp_model(training_data)
"""

from .model import (
    get_train_test_split,
    train_npp_model,
)
from .sampling import (
    setup_training_grid,
)
from .validation import validate_model

__all__ = [
    "setup_training_grid",
    "train_npp_model",
    "get_train_test_split",
    "validate_model",
]
