"""
Utility helpers for estimating operating targets from historical results.

The main entry point, :func:`estimate_operating_targets_from_results`, filters a
results DataFrame to matching categorical parameters and then performs a
multi-dimensional interpolation over the remaining numerical inputs to estimate
``Pmax`` and ``dT_gc_out_min_offset`` for a new case. When linear interpolation
is not possible (e.g. outside the convex hull) the function falls back to a
nearest-neighbour estimate.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping, MutableMapping, Sequence

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.spatial import qhull

__all__ = ["estimate_operating_targets_from_results"]

# Columns that should be matched exactly before interpolation is attempted.
DEFAULT_CATEGORICAL_MATCH_COLUMNS: Sequence[str] = (
    "refrigerant",
    "n_comp",
    "comp_PR_uniform",
    "ihx_on",
    "comp_model_index",
    "sink",
    "source",
)

# Numerical columns that describe the operating point to interpolate over.
DEFAULT_NUMERIC_FEATURE_COLUMNS: Sequence[str] = (
    "T_source_in",
    "dT_GC_min_tar",
    "T_air_inlet",
    "T_air_tar",
    "Tref_evap_inlet",
    "dTmin_rhx",
    "motor_eff",
    "dT_rhx1_cold",
    "dT_rhx2_cold",
    "dT_rhx3_cold",
    "comp1_power_frac",
    "comp2_power_frac",
    "comp3_power_frac",
)

TARGET_COLUMNS: Sequence[str] = ("Pmax", "dT_gc_out_min_offset")


def _coerce_to_dataframe(results: pd.DataFrame | Path | str) -> pd.DataFrame:
    """Return a DataFrame regardless of whether ``results`` is already a DataFrame or a file path."""
    if isinstance(results, pd.DataFrame):
        return results.copy()
    if isinstance(results, (str, Path)):
        return pd.read_csv(Path(results))
    raise TypeError(
        "results must be a pandas.DataFrame or a path to a CSV file containing past results."
    )


def _ensure_numeric(series: pd.Series) -> pd.Series:
    """Convert a Series to float, preserving indices."""
    return pd.to_numeric(series, errors="coerce")


def _filter_by_categories(
    df: pd.DataFrame,
    inputs: Mapping[str, object],
    categorical_columns: Iterable[str],
) -> pd.DataFrame:
    """Filter the DataFrame so that categorical columns match the provided inputs exactly."""
    filtered = df
    for column in categorical_columns:
        if column in filtered.columns and column in inputs:
            filtered = filtered[filtered[column] == inputs[column]]
    return filtered


def _select_feature_columns(
    df: pd.DataFrame,
    inputs: Mapping[str, object],
    candidate_columns: Iterable[str],
) -> list[str]:
    """
    Determine which feature columns to use for interpolation.

    Columns that are not present in the DataFrame or do not have a corresponding
    value in ``inputs`` are skipped. Columns with no variance in the filtered
    dataset are ignored because they do not contribute to the interpolation.
    """
    selected: list[str] = []
    for column in candidate_columns:
        if column not in df.columns or column not in inputs:
            continue
        numeric_series = _ensure_numeric(df[column])
        if numeric_series.isnull().all():
            continue
        if numeric_series.nunique(dropna=True) <= 1:
            # Column provides no variation; ensure input is compatible but skip it.
            continue
        selected.append(column)
    return selected


def _prepare_training_data(
    df: pd.DataFrame,
    feature_columns: Sequence[str],
    target_columns: Sequence[str],
) -> tuple[np.ndarray, np.ndarray]:
    """
    Extract numerical arrays for interpolation.

    Rows containing NaNs in either the feature or target columns are removed.
    Duplicate feature rows are dropped to avoid issues with the interpolation
    backend; when duplicates exist the first occurrence is retained.
    """
    features = df.loc[:, feature_columns].apply(_ensure_numeric)
    targets = df.loc[:, target_columns].apply(_ensure_numeric)

    mask = ~(features.isnull().any(axis=1) | targets.isnull().any(axis=1))
    features = features.loc[mask]
    targets = targets.loc[mask]

    if features.empty:
        return np.empty((0, len(feature_columns))), np.empty((0, len(target_columns)))

    merged = pd.concat([features, targets], axis=1)
    merged = merged.drop_duplicates(subset=feature_columns)
    features = merged.loc[:, feature_columns]
    targets = merged.loc[:, target_columns]

    return features.to_numpy(dtype=float), targets.to_numpy(dtype=float)


def _as_scalar(value: object) -> float:
    """Convert interpolator outputs to a scalar float, handling array-like responses."""
    array = np.asarray(value).squeeze()
    if array.size == 0:
        return np.nan
    return float(array)


def estimate_operating_targets_from_results(
    results: pd.DataFrame | Path | str,
    inputs: Mapping[str, object],
    *,
    categorical_match_columns: Sequence[str] = DEFAULT_CATEGORICAL_MATCH_COLUMNS,
    numeric_feature_columns: Sequence[str] = DEFAULT_NUMERIC_FEATURE_COLUMNS,
    target_columns: Sequence[str] = TARGET_COLUMNS,
) -> MutableMapping[str, float]:
    """
    Estimate ``Pmax`` and ``dT_gc_out_min_offset`` for a given operating point.

    Parameters
    ----------
    results:
        Either a ``pandas.DataFrame`` of past optimisation results or a path to
        a CSV file containing those results.
    inputs:
        Mapping of column names to the desired operating point. The keys should
        include the fields listed in ``numeric_feature_columns`` and any
        categorical identifiers relevant to the interpolation.
    categorical_match_columns:
        Columns that must match exactly before interpolation is attempted.
        Defaults to a set of categorical parameters such as refrigerant and
        compressor configuration.
    numeric_feature_columns:
        Numerical columns used as interpolation features. Columns are ignored if
        they are absent from ``results`` or ``inputs`` or if they provide no
        variation within the filtered dataset.
    target_columns:
        Columns that will be estimated, defaulting to ``(\"Pmax\", \"dT_gc_out_min_offset\")``.

    Returns
    -------
    dict
        Mapping of each target column to its estimated float value.

    Raises
    ------
    ValueError
        If the filters remove every row or if the remaining dataset lacks the
        necessary information to perform interpolation.
    TypeError
        If ``results`` is neither a DataFrame nor a path-like object.
    """
    df = _coerce_to_dataframe(results)
    if df.empty:
        raise ValueError("The results DataFrame is empty; cannot estimate targets.")

    missing_targets = [col for col in target_columns if col not in df.columns]
    if missing_targets:
        raise ValueError(f"Required target columns {missing_targets} are missing from results.")

    df = df.dropna(subset=target_columns)

    filtered = _filter_by_categories(df, inputs, categorical_match_columns)
    if filtered.empty:
        raise ValueError(
            "No historical results match the provided categorical inputs; unable to perform interpolation."
        )

    feature_columns = _select_feature_columns(filtered, inputs, numeric_feature_columns)
    if not feature_columns:
        # No variation remained after filtering. Use mean of historical targets.
        target_means = filtered.loc[:, target_columns].apply(_ensure_numeric).mean()
        if target_means.isnull().any():
            raise ValueError("Not enough valid historical data to compute a fallback estimate.")
        return {col: float(target_means[col]) for col in target_columns}

    feature_array, target_array = _prepare_training_data(filtered, feature_columns, target_columns)

    if feature_array.shape[0] == 0:
        raise ValueError("No valid rows remain after removing NaNs; interpolation is impossible.")

    # Build the query point.
    try:
        query_point = np.array([float(inputs[col]) for col in feature_columns], dtype=float)
    except KeyError as exc:
        raise ValueError(f"Missing required input '{exc.args[0]}' for interpolation.") from exc

    if feature_array.shape[0] == 1:
        return {col: float(target_array[0, idx]) for idx, col in enumerate(target_columns)}

    ndim = feature_array.shape[1]
    estimates = np.full(len(target_columns), np.nan, dtype=float)

    for idx, target_col in enumerate(target_columns):
        values = target_array[:, idx]
        if feature_array.shape[0] > ndim:
            try:
                lin_interp = LinearNDInterpolator(feature_array, values, fill_value=np.nan)
                estimates[idx] = _as_scalar(lin_interp(query_point))
            except qhull.QhullError:
                estimates[idx] = np.nan
        if not np.isfinite(estimates[idx]):
            nearest_interp = NearestNDInterpolator(feature_array, values)
            estimates[idx] = _as_scalar(nearest_interp(query_point))

    return {col: float(estimates[i]) for i, col in enumerate(target_columns)}
