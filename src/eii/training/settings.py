"""
Settings for potential-NPP model training.
"""

from .model_card import load_model_card

# Quantile regression forest of the prediction intervals
RF_NUM_TREES = 200
RF_MIN_LEAF_POPULATION = 5
RF_VARIABLES_PER_SPLIT = 5
RF_BAG_FRACTION = 0.7

# Ordered predictors come from the v2 model card
PREDICTOR_VARIABLES = load_model_card("v2").predictors
