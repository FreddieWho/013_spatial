"""Cross-fitted readouts for shared and structure-specific spatial effects.

This module treats a biological structure as an externally supplied readout
target.  It never discovers structures by clustering and never assigns a
biological name to a latent factor.  Features are obtained from the centered
spot-by-gene effect implied by a complete frozen model fit.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np

from .field_exports import reconstructed_linear_effect

def _array(export: Mapping[str, object], key: str) -> np.ndarray:
    if key not in export:
        raise ValueError(f"field export is missing {key}")
    return np.asarray(export[key])


def _center_by_section(field: np.ndarray, section_ids: np.ndarray) -> np.ndarray:
    if field.ndim != 2 or section_ids.ndim != 1 or len(field) != len(section_ids):
        raise ValueError("field and section IDs are not spot-aligned")
    centered = np.asarray(field, dtype=float).copy()
    for section_id in np.unique(section_ids):
        indices = np.flatnonzero(section_ids == section_id)
        centered[indices] -= centered[indices].mean(axis=0, keepdims=True)
    return centered


def _check_compatible_exports(
    training_export: Mapping[str, object], heldout_export: Mapping[str, object]
) -> tuple[np.ndarray, np.ndarray]:
    train_field = _array(training_export, "field_mean")
    test_field = _array(heldout_export, "field_mean")
    train_genes = _array(training_export, "evaluation_gene_ids")
    test_genes = _array(heldout_export, "evaluation_gene_ids")
    if train_field.ndim != 2 or test_field.ndim != 2:
        raise ValueError("field means must be matrices")
    if train_field.shape[1] != test_field.shape[1]:
        raise ValueError("training and held-out factor dimensions differ")
    train_loading = _array(training_export, "evaluation_loading")
    test_loading = _array(heldout_export, "evaluation_loading")
    train_amplitude = _array(training_export, "factor_amplitude")
    test_amplitude = _array(heldout_export, "factor_amplitude")
    if train_loading.shape != test_loading.shape or not np.allclose(train_loading, test_loading):
        raise ValueError("training and held-out loadings differ")
    if train_amplitude.shape != test_amplitude.shape or not np.allclose(train_amplitude, test_amplitude):
        raise ValueError("training and held-out factor amplitudes differ")
    if train_genes.tolist() != test_genes.tolist():
        raise ValueError("training and held-out evaluation genes differ")
    train_ids = _array(training_export, "section_ids")
    test_ids = _array(heldout_export, "section_ids")
    if "spatial_rate_effect_centered" in training_export:
        train_effect = _array(training_export, "spatial_rate_effect_centered").astype(float)
    else:
        train_effect = reconstructed_linear_effect(dict(training_export))
    if "spatial_rate_effect_centered" in heldout_export:
        test_effect = _array(heldout_export, "spatial_rate_effect_centered").astype(float)
    else:
        test_effect = reconstructed_linear_effect(dict(heldout_export))
    if train_effect.ndim != 2 or test_effect.ndim != 2:
        raise ValueError("spatial effects must be matrices")
    if train_effect.shape[1] != test_effect.shape[1] or train_effect.shape[1] != len(train_genes):
        raise ValueError("training and held-out effect gene dimensions differ")
    if train_effect.shape[0] != len(train_ids) or test_effect.shape[0] != len(test_ids):
        raise ValueError("spatial effects are not spot-aligned")
    return train_effect, test_effect


def effect_coordinates(
    training_export: Mapping[str, object],
    heldout_export: Mapping[str, object],
    *,
    max_rank: int | None = None,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Build common low-rank coordinates from a rotation-stable effect.

    The right singular basis is learned from the training outer fold and then
    applied to the held-out fold.  The SVD is performed on the small factor
    cross-product, so the full spot-by-gene matrix is never materialised.
    """
    train_effect, test_effect = _check_compatible_exports(
        training_export, heldout_export
    )
    if max_rank is not None and max_rank < 1:
        raise ValueError("max_rank must be positive or None")
    rate_primary = (
        "spatial_rate_effect_centered" in training_export
        and "spatial_rate_effect_centered" in heldout_export
    )
    if rate_primary:
        # C = center(exp(F)) @ (W * amplitude).T is rank at most K.  Compute
        # its right singular basis through the K-dimensional factor side so a
        # 150k-by-4k export does not trigger a dense full SVD.
        train_field = _center_by_section(
            np.exp(np.clip(_array(training_export, "field_mean").astype(float), -30.0, 30.0)),
            _array(training_export, "section_ids"),
        )
        test_field = _center_by_section(
            np.exp(np.clip(_array(heldout_export, "field_mean").astype(float), -30.0, 30.0)),
            _array(heldout_export, "section_ids"),
        )
        loading = _array(training_export, "evaluation_loading").astype(float)
        amplitude = _array(training_export, "factor_amplitude").astype(float)
        if loading.shape[1] != train_field.shape[1] or amplitude.shape != (train_field.shape[1],):
            raise ValueError("rate effect factor dimensions differ")
        test_loading = _array(heldout_export, "evaluation_loading").astype(float)
        test_amplitude = _array(heldout_export, "factor_amplitude").astype(float)
        if not np.allclose(loading, test_loading) or not np.allclose(amplitude, test_amplitude):
            raise ValueError("training and held-out rate parameters differ")
        for field, effect in (
            (train_field, train_effect),
            (test_field, test_effect),
        ):
            stored = np.asarray(effect, dtype=float)
            indices = np.linspace(0, len(stored) - 1, num=min(8, len(stored)), dtype=int)
            expected = field[indices] @ (loading * amplitude[None, :]).T
            if not np.allclose(stored[indices], expected, rtol=2e-5, atol=2e-5):
                raise ValueError("stored rate effect does not match frozen factor parameters")
        field_u, field_s, field_vt = np.linalg.svd(train_field, full_matrices=False)
        tolerance = max(train_field.shape) * np.finfo(float).eps * max(float(field_s[0]), 1.0) if len(field_s) else 0.0
        active_field = field_s > tolerance
        if not np.any(active_field):
            raise ValueError("training held-out field has no spatial variation")
        low_rank_effect = (
            field_s[active_field, None]
            * field_vt[active_field]
        ) @ (loading * amplitude[None, :]).T
        _, effect_s, effect_vt = np.linalg.svd(low_rank_effect, full_matrices=False)
    else:
        _, effect_s, effect_vt = np.linalg.svd(train_effect, full_matrices=False)
    tolerance = max(train_effect.shape) * np.finfo(float).eps * max(float(effect_s[0]), 1.0) if len(effect_s) else 0.0
    active = effect_s > tolerance
    if not np.any(active):
        raise ValueError("training held-out field has no spatial variation")
    full_rank = int(np.sum(active))
    effect_s = effect_s[active]
    effect_vt = effect_vt[active]
    rank = int(len(effect_s))
    if max_rank is not None:
        rank = min(rank, max_rank)
    if rank < 1:
        raise ValueError("effect representation has no estimable rank")
    right_basis = effect_vt[:rank].T
    if rate_primary:
        projection = (loading * amplitude[None, :]).T @ right_basis
        train_z = train_field @ projection
        test_z = test_field @ projection
    else:
        train_z = train_effect @ right_basis
        test_z = test_effect @ right_basis
    return train_z, test_z, {
        "representation": (
            "centered_mnsf_spatial_rate_contribution"
            if "spatial_rate_effect_centered" in training_export
            else "rotation_stable_centered_linear_spot_by_gene_effect"
        ),
        "effect_definition": (
            "section_center(exp(field_mean)@"
            "(evaluation_loading*factor_amplitude).T)"
            if "spatial_rate_effect_centered" in training_export
            else "section_center(field_mean)@evaluation_loading.T"
        ),
        "rank": rank,
        "full_effect_rank": full_rank,
        "singular_values": effect_s[:rank].astype(float).tolist(),
        "max_rank": max_rank,
        "factor_naming": "forbidden",
    }


def _weights(weights: np.ndarray | None, n: int, *, name: str) -> np.ndarray:
    if weights is None:
        return np.ones(n, dtype=float)
    value = np.asarray(weights, dtype=float)
    if value.shape != (n,) or not np.isfinite(value).all() or np.any(value < 0):
        raise ValueError(f"{name} must be finite, non-negative and spot-aligned")
    if float(value.sum()) <= 0:
        raise ValueError(f"{name} must have positive total weight")
    return value


def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return np.average(values, axis=0, weights=weights)


def _standardize_features(
    train: np.ndarray, heldout: np.ndarray, weights: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    mean = _weighted_mean(train, weights)
    centered = train - mean
    scale = np.sqrt(_weighted_mean(centered * centered, weights))
    scale = np.where(scale > 1e-12, scale, 1.0)
    return (train - mean) / scale, (heldout - mean) / scale


def _ridge_fit_predict(
    train: np.ndarray,
    heldout: np.ndarray,
    target_train: np.ndarray,
    *,
    alpha: float,
    weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if alpha <= 0:
        raise ValueError("ridge alpha must be positive")
    mean_x = _weighted_mean(train, weights)
    mean_y = float(_weighted_mean(target_train, weights))
    centered_x = train - mean_x
    centered_y = target_train - mean_y
    sqrt_weights = np.sqrt(weights)
    weighted_x = centered_x * sqrt_weights[:, None]
    weighted_y = centered_y * sqrt_weights
    system = weighted_x.T @ weighted_x + alpha * np.eye(train.shape[1])
    coefficient = np.linalg.solve(system, weighted_x.T @ weighted_y)
    intercept = mean_y - mean_x @ coefficient
    return (
        heldout @ coefficient + intercept,
        coefficient,
        train @ coefficient + intercept,
    )


def _one_dimensional_readout(
    train_score: np.ndarray,
    heldout_score: np.ndarray,
    target_train: np.ndarray,
    weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    design = np.column_stack([np.ones(len(train_score)), train_score])
    sqrt_weights = np.sqrt(weights)
    coefficient, *_ = np.linalg.lstsq(
        design * sqrt_weights[:, None], target_train * sqrt_weights, rcond=None
    )
    heldout_design = np.column_stack([np.ones(len(heldout_score)), heldout_score])
    return heldout_design @ coefficient, design @ coefficient


def _auc(
    target: np.ndarray, prediction: np.ndarray, weights: np.ndarray
) -> float | None:
    classes = np.unique(target)
    if len(classes) != 2:
        return None
    positive = target == classes[1]
    n_positive = int(np.sum(positive))
    n_negative = len(target) - n_positive
    if n_positive == 0 or n_negative == 0:
        return None
    positive_weight = float(weights[positive].sum())
    negative_weight = float(weights[~positive].sum())
    if positive_weight <= 0 or negative_weight <= 0:
        return None
    positive_scores = prediction[positive]
    negative_scores = prediction[~positive]
    pairwise = (positive_scores[:, None] > negative_scores[None, :]).astype(float)
    pairwise += 0.5 * (positive_scores[:, None] == negative_scores[None, :])
    pair_weights = weights[positive, None] * weights[~positive][None, :]
    return float(np.sum(pairwise * pair_weights) / (positive_weight * negative_weight))


def _metrics(
    target: np.ndarray, prediction: np.ndarray, weights: np.ndarray
) -> dict[str, float | None]:
    target = np.asarray(target, dtype=float)
    prediction = np.asarray(prediction, dtype=float)
    mse = float(np.average((target - prediction) ** 2, weights=weights))
    target_mean = float(np.average(target, weights=weights))
    prediction_mean = float(np.average(prediction, weights=weights))
    covariance = float(np.average((target - target_mean) * (prediction - prediction_mean), weights=weights))
    target_sd = float(np.sqrt(np.average((target - target_mean) ** 2, weights=weights)))
    prediction_sd = float(np.sqrt(np.average((prediction - prediction_mean) ** 2, weights=weights)))
    correlation = (
        covariance / (target_sd * prediction_sd)
        if target_sd > 1e-12 and prediction_sd > 1e-12 else None
    )
    return {"mse": mse, "correlation": correlation, "auc": _auc(target, prediction, weights)}


def _task_names(labels: Mapping[str, object]) -> tuple[str, str]:
    names = tuple(str(name) for name in labels)
    if len(names) != 2:
        raise ValueError("shared structure readout requires exactly two structure targets")
    if len(set(names)) != 2:
        raise ValueError("structure target names must be unique")
    return names[0], names[1]


def evaluate_shared_specific_readout(
    training_z: np.ndarray,
    heldout_z: np.ndarray,
    training_labels: Mapping[str, object],
    heldout_labels: Mapping[str, object],
    *,
    alpha: float = 1.0,
    representation: Mapping[str, object] | None = None,
    training_weights: np.ndarray | None = None,
    heldout_weights: np.ndarray | None = None,
) -> dict[str, object]:
    """Compare a shared ecological score with shared plus specific residuals."""
    train_z = np.asarray(training_z, dtype=float)
    test_z = np.asarray(heldout_z, dtype=float)
    if train_z.ndim != 2 or test_z.ndim != 2 or train_z.shape[1] != test_z.shape[1]:
        raise ValueError("readout feature matrices have incompatible shapes")
    names = _task_names(training_labels)
    if tuple(str(name) for name in heldout_labels) != names:
        raise ValueError("training and held-out structure targets differ")
    train_y = []
    test_y = []
    for name in names:
        y_train = np.asarray(training_labels[name], dtype=float)
        y_test = np.asarray(heldout_labels[name], dtype=float)
        if y_train.shape != (len(train_z),) or y_test.shape != (len(test_z),):
            raise ValueError(f"label target is not aligned for {name}")
        if not np.isfinite(y_train).all() or not np.isfinite(y_test).all():
            raise ValueError(f"label target is non-finite for {name}")
        train_y.append(y_train)
        test_y.append(y_test)
    train_weights = _weights(training_weights, len(train_z), name="training_weights")
    test_weights = _weights(heldout_weights, len(test_z), name="heldout_weights")
    train_z, test_z = _standardize_features(train_z, test_z, train_weights)
    standardized_targets = []
    for y in train_y:
        mean = float(np.average(y, weights=train_weights))
        scale = float(np.sqrt(np.average((y - mean) ** 2, weights=train_weights)))
        standardized_targets.append((y - mean) / (scale if scale > 1e-12 else 1.0))
    shared_target = np.mean(np.vstack(standardized_targets), axis=0)
    shared_test_score, shared_coefficient, shared_train_score = _ridge_fit_predict(
        train_z, test_z, shared_target, alpha=alpha, weights=train_weights
    )
    shared_predictions: dict[str, np.ndarray] = {}
    shared_train_predictions: dict[str, np.ndarray] = {}
    for name, y in zip(names, train_y):
        prediction, train_prediction = _one_dimensional_readout(
            shared_train_score, shared_test_score, y, train_weights
        )
        shared_predictions[name] = prediction
        shared_train_predictions[name] = train_prediction

    shared_mean = float(np.average(shared_train_score, weights=train_weights))
    centered_shared_train = shared_train_score - shared_mean
    centered_shared_test = shared_test_score - shared_mean
    denominator = float(train_weights @ (centered_shared_train * centered_shared_train))
    if denominator > 1e-12:
        projection = (train_weights * centered_shared_train) @ train_z / denominator
        train_specific_z = train_z - centered_shared_train[:, None] * projection[None, :]
        test_specific_z = test_z - centered_shared_test[:, None] * projection[None, :]
    else:
        projection = np.zeros(train_z.shape[1], dtype=float)
        train_specific_z = train_z.copy()
        test_specific_z = test_z.copy()
    specific_predictions: dict[str, np.ndarray] = {}
    specific_train_predictions: dict[str, np.ndarray] = {}
    specific_norms: dict[str, float] = {}
    for name, y, baseline_train in zip(names, train_y, (shared_train_predictions[name] for name in names)):
        residual = y - baseline_train
        prediction_residual, coefficient, train_residual = _ridge_fit_predict(
            train_specific_z,
            test_specific_z,
            residual,
            alpha=alpha,
            weights=train_weights,
        )
        specific_predictions[name] = shared_predictions[name] + prediction_residual
        specific_train_predictions[name] = baseline_train + train_residual
        specific_norms[name] = float(np.linalg.norm(coefficient))

    shared_metrics = {
        name: _metrics(y, shared_predictions[name], test_weights)
        for name, y in zip(names, test_y)
    }
    specific_metrics = {
        name: _metrics(y, specific_predictions[name], test_weights)
        for name, y in zip(names, test_y)
    }
    shared_mse = float(np.mean([shared_metrics[name]["mse"] for name in names]))
    specific_mse = float(np.mean([specific_metrics[name]["mse"] for name in names]))
    auc_deltas = {
        name: (
            specific_metrics[name]["auc"] - shared_metrics[name]["auc"]
            if specific_metrics[name]["auc"] is not None and shared_metrics[name]["auc"] is not None
            else None
        )
        for name in names
    }
    return {
        "schema": "r04.structure_readout.v1",
        "target_count": 2,
        "target_names": list(names),
        "alpha": float(alpha),
        "representation": dict(representation or {}),
        "n_training_spots": int(len(train_z)),
        "n_heldout_spots": int(len(test_z)),
        "patient_or_section_weighting": "caller_supplied_sample_weights_or_equal_spots",
        "shared_only": {
            "metrics": shared_metrics,
            "mean_mse": shared_mse,
        },
        "shared_plus_specific_residual": {
            "metrics": specific_metrics,
            "mean_mse": specific_mse,
            "specific_residual_coefficient_norm": specific_norms,
        },
        "specific_increment": {
            "mean_mse_delta_specific_minus_shared": specific_mse - shared_mse,
            "auc_delta_by_target": auc_deltas,
            "interpretation": (
                "descriptive heldout readout; a positive AUC or negative MSE change "
                "supports reproducible target-specific residual information"
            ),
        },
        "factor_naming": "forbidden",
    }
