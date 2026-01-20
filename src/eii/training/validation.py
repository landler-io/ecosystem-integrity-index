import ee
import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


def validate_model(
    validation_set: ee.FeatureCollection,
    model_asset_path: str,
    response_vars: list[str] | None = None,
    prediction_names: list[str] | None = None,
) -> pd.DataFrame:
    """
    Validate a trained model using a FeatureCollection of validation points.

    This method applies the classifier directly to the features, which is
    extremely fast compared to image-based validation, as it utilizes the
    predictor values already stored in the table.

    Args:
        validation_set: FeatureCollection containing predictors and actual response.
        model_asset_path: Path to the trained GEE classifier.
        response_vars: List of property names for actual values (e.g. ['current_npp', 'npp_std'])
        prediction_names: List of property names for predicted values (matching model output)

    Returns:
        DataFrame containing metrics for each response variable.
    """
    if response_vars is None:
        response_vars = ["current_npp"]
    if prediction_names is None:
        prediction_names = ["classification"]
    print(f"Loading model from {model_asset_path}...")
    model = ee.Classifier.load(model_asset_path)

    predictions = validation_set.classify(model)

    return calculate_metrics(predictions, response_vars, prediction_names)


def export_validation_predictions(
    validation_set: ee.FeatureCollection,
    model_asset_path: str,
    output_asset_path: str,
    prediction_name: str = "classification",
    description: str = "export_validation_predictions",
):
    """
    Applies the model to the validation set and exports the result to a table asset.
    Useful for creating a persistent record of validation results that can be analyzed later.
    """
    print(f"Loading model from {model_asset_path}...")
    model = ee.Classifier.load(model_asset_path)

    predictions = validation_set.classify(model, outputName=prediction_name)

    print(f"Exporting predictions to {output_asset_path}...")
    task = ee.batch.Export.table.toAsset(
        collection=predictions, description=description, assetId=output_asset_path
    )
    task.start()
    return task


def calculate_metrics(
    fc: ee.FeatureCollection, response_vars: list[str], prediction_names: list[str]
) -> pd.DataFrame:
    """
    Computes validation metrics from a FeatureCollection containing actual and predicted values.
    Can be used with either on-the-fly predictions or a saved Table Asset.
    """
    cols_to_keep = response_vars + prediction_names

    print("Extracting predictions (client-side)...")
    data = fc.reduceColumns(
        reducer=ee.Reducer.toList(len(cols_to_keep)), selectors=cols_to_keep
    ).getInfo()["list"]

    df = pd.DataFrame(data, columns=cols_to_keep)

    metrics = []

    for actual_col, pred_col in zip(response_vars, prediction_names, strict=False):
        if actual_col not in df.columns or pred_col not in df.columns:
            print(f"Skipping {actual_col}: data missing")
            continue

        y_true = pd.to_numeric(df[actual_col], errors="coerce")
        y_pred = pd.to_numeric(df[pred_col], errors="coerce")

        # Drop NaNs
        valid_mask = y_true.notna() & y_pred.notna()
        y_true = y_true[valid_mask]
        y_pred = y_pred[valid_mask]

        if len(y_true) == 0:
            print(f"No valid data points for {actual_col}")
            continue

        mae = mean_absolute_error(y_true, y_pred)
        rmse = np.sqrt(mean_squared_error(y_true, y_pred))
        r2 = r2_score(y_true, y_pred)

        metrics.append({"Target": actual_col, "MAE": mae, "RMSE": rmse, "R2": r2, "N": len(y_true)})

    results_df = pd.DataFrame(metrics)
    return results_df
