"""
NPP model training functions.
"""

from __future__ import annotations

import random

import ee

from .settings import (
    CV_BUFFER_DEG,
    CV_GRID_SIZE_DEG,
    MODEL_ASSETS_PATH,
    PREDICTOR_VARIABLES,
    RF_BAG_FRACTION,
    RF_MIN_LEAF_POPULATION,
    RF_NUM_TREES,
    RF_SEED,
    RF_VARIABLES_PER_SPLIT,
    TRAIN_TEST_SPLIT_RATIO,
)


def get_train_test_split(
    training_data: ee.FeatureCollection,
    split_ratio: float = TRAIN_TEST_SPLIT_RATIO,
    seed: int = RF_SEED,
    cv_grid_size: int = CV_GRID_SIZE_DEG,
    cv_buffer_size: float = CV_BUFFER_DEG,
) -> tuple[ee.FeatureCollection, ee.FeatureCollection]:
    """
    Perform spatially stratified train/test split.

    Uses a grid (e.g., 2 degrees) to create spatial blocks.
    Can optionally apply a negative buffer (margin) around each block
    to ensuring physical separation between training and validation sets.

    Args:
        training_data: FeatureCollection.
        split_ratio: Fraction for training (default 0.9).
        seed: Random seed for reproducibility.
        cv_grid_size: Grid size in degrees for cross-validation blocks.
        cv_buffer_size: Buffer size in degrees to exclude from block edges.
                        0.0 means no buffer. 0.5 means 0.5 deg excluded from all sides.

    Returns:
        Tuple of (training_set, validation_set).
    """

    print(CV_GRID_SIZE_DEG, CV_BUFFER_DEG)

    def add_block_info(feature):
        coords = feature.geometry().coordinates()
        lon = coords.get(0)
        lat = coords.get(1)

        # Shift to positive range for easier modulo/floor
        lon_shifted = ee.Number(lon).add(180)
        lat_shifted = ee.Number(lat).add(90)

        x = lon_shifted.divide(cv_grid_size).floor()
        y = lat_shifted.divide(cv_grid_size).floor()

        # relative position within the block [0, cv_grid_size)
        x_rel = lon_shifted.mod(cv_grid_size)
        y_rel = lat_shifted.mod(cv_grid_size)

        block_id = y.multiply(1000).add(x).toInt()

        # Keep if buffer <= pos <= size - buffer
        inner_cond = (
            x_rel.gte(cv_buffer_size)
            .And(x_rel.lte(ee.Number(cv_grid_size).subtract(cv_buffer_size)))
            .And(y_rel.gte(cv_buffer_size))
            .And(y_rel.lte(ee.Number(cv_grid_size).subtract(cv_buffer_size)))
        )

        return feature.set({"cv_block_id": block_id, "cv_keep": inner_cond})

    data_with_info = training_data.map(add_block_info)

    if cv_buffer_size > 0:
        data_to_split = data_with_info.filter(ee.Filter.eq("cv_keep", 1))
    else:
        data_to_split = data_with_info

    distinct_blocks = data_to_split.aggregate_array("cv_block_id").distinct().getInfo()

    random.seed(seed)
    random.shuffle(distinct_blocks)

    split_index = int(len(distinct_blocks) * split_ratio)
    training_blocks = distinct_blocks[:split_index]
    validation_blocks = distinct_blocks[split_index:]

    if len(distinct_blocks) > 5000:
        # Fallback to server-side hashing if too many blocks
        print(
            f"Warning: Large number of CV blocks ({len(distinct_blocks)}). Using hashed block split."
        )

        def add_hashed_split(f):
            bid = ee.Number(f.get("cv_block_id"))
            h = bid.multiply(12345).add(seed).mod(10000).divide(10000)
            return f.set("cv_random", h)

        data_hashed = data_to_split.map(add_hashed_split)

        training_set = data_hashed.filter(ee.Filter.lt("cv_random", split_ratio))
        validation_set = data_hashed.filter(ee.Filter.gte("cv_random", split_ratio))

    else:
        training_set = data_to_split.filter(ee.Filter.inList("cv_block_id", training_blocks))
        validation_set = data_to_split.filter(ee.Filter.inList("cv_block_id", validation_blocks))

    return training_set, validation_set


def get_train_test_split_server_side(
    training_data: ee.FeatureCollection,
    split_ratio: float = TRAIN_TEST_SPLIT_RATIO,
    seed: int = RF_SEED,
) -> tuple[ee.FeatureCollection, ee.FeatureCollection]:
    """
    Perform train/test split entirely server-side using random column.

    Args:
        training_data: FeatureCollection to split.
        split_ratio: Fraction for training (default 0.9).
        seed: Random seed for reproducibility.

    Returns:
        Tuple of (training_set, validation_set).
    """
    import ee

    data_with_random = training_data.randomColumn("split_random", seed)
    training_set = data_with_random.filter(ee.Filter.lt("split_random", split_ratio))
    validation_set = data_with_random.filter(ee.Filter.gte("split_random", split_ratio))

    return training_set, validation_set


def train_npp_model(
    training_data: ee.FeatureCollection,
    predictor_names: list[str] | None = None,
    response_property: str = "longterm_avg_npp_sum",
    output_asset_path: str | None = None,
    num_trees: int = RF_NUM_TREES,
    min_leaf_population: int = RF_MIN_LEAF_POPULATION,
    variables_per_split: int = RF_VARIABLES_PER_SPLIT,
    bag_fraction: float = RF_BAG_FRACTION,
    seed: int = RF_SEED,
    export: bool = True,
) -> tuple[ee.Classifier, ee.batch.Task | None]:
    """
    Train Random Forest model for NPP prediction.

    Args:
        training_data: FeatureCollection with predictor variables and response.
        predictor_names: List of predictor property names. If None, infers from
            first feature (requires getInfo).
        response_property: Name of the response variable property.
        output_asset_path: Asset path to export model. If None, generates default.
        num_trees: Number of trees in the forest.
        min_leaf_population: Minimum samples in a leaf.
        variables_per_split: Number of variables to consider per split.
        bag_fraction: Fraction of data to bag per tree.
        seed: Random seed.
        export: Whether to export the model to an asset.

    Returns:
        Tuple of (trained_model, export_task or None).
    """
    import ee

    if predictor_names is None:
        predictor_names = PREDICTOR_VARIABLES

    model = (
        ee.Classifier.smileRandomForest(
            numberOfTrees=num_trees,
            minLeafPopulation=min_leaf_population,
            variablesPerSplit=variables_per_split,
            bagFraction=bag_fraction,
            seed=seed,
        )
        .setOutputMode("REGRESSION")
        .train(
            features=training_data,
            classProperty=response_property,
            inputProperties=predictor_names,
        )
    )

    export_task = None
    if export:
        if output_asset_path is None:
            output_asset_path = f"{MODEL_ASSETS_PATH}/potential_npp_classifier"

        export_task = ee.batch.Export.classifier.toAsset(
            classifier=model,
            description="Export_NPP_Classifier",
            assetId=output_asset_path,
        )
        export_task.start()

    return model, export_task


def train_npp_models(
    training_data: ee.FeatureCollection,
    predictor_names: list[str] | None = None,
    output_mean_path: str | None = None,
    output_std_path: str | None = None,
    num_trees: int = RF_NUM_TREES,
    min_leaf_population: int = RF_MIN_LEAF_POPULATION,
    variables_per_split: int = RF_VARIABLES_PER_SPLIT,
    bag_fraction: float = RF_BAG_FRACTION,
    seed: int = RF_SEED,
    export: bool = True,
) -> dict:
    """
    Train Random Forest models for both NPP mean and NPP std prediction.

    Args:
        training_data: FeatureCollection with predictor variables, 'longterm_avg_npp_sum', and 'longterm_avg_npp_sd'.
        predictor_names: List of predictor property names.
        output_mean_path: Asset path for mean model. If None, generates default.
        output_std_path: Asset path for std model. If None, generates default.
        num_trees: Number of trees in the forest.
        min_leaf_population: Minimum samples in a leaf.
        variables_per_split: Number of variables to consider per split.
        bag_fraction: Fraction of data to bag per tree.
        seed: Random seed.
        export: Whether to export the models to assets.

    Returns:
        Dictionary with 'mean_model', 'std_model', 'mean_task', 'std_task'.
    """
    if output_mean_path is None:
        output_mean_path = f"{MODEL_ASSETS_PATH}/potential_npp_mean_classifier"
    if output_std_path is None:
        output_std_path = f"{MODEL_ASSETS_PATH}/potential_npp_std_classifier"

    # Train mean model
    mean_model, mean_task = train_npp_model(
        training_data=training_data,
        predictor_names=predictor_names,
        response_property="longterm_avg_npp_sum",
        output_asset_path=output_mean_path,
        num_trees=num_trees,
        min_leaf_population=min_leaf_population,
        variables_per_split=variables_per_split,
        bag_fraction=bag_fraction,
        seed=seed,
        export=export,
    )

    # Train std model
    std_model, std_task = train_npp_model(
        training_data=training_data,
        predictor_names=predictor_names,
        response_property="longterm_avg_npp_sd",
        output_asset_path=output_std_path,
        num_trees=num_trees,
        min_leaf_population=min_leaf_population,
        variables_per_split=variables_per_split,
        bag_fraction=bag_fraction,
        seed=seed,
        export=export,
    )

    return {
        "mean_model": mean_model,
        "std_model": std_model,
        "mean_task": mean_task,
        "std_task": std_task,
    }
