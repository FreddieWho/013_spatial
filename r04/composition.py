"""Composition adjustment for continuous fields.

The implementation is deliberately a regression residual, not a classifier.
The estimator accepts posterior draws from an upstream composition model when
available; the simple deterministic path is useful for smoke tests only and
records that uncertainty was not propagated.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CompositionAdjustment:
    adjusted: np.ndarray
    expected: np.ndarray
    uncertainty_available: bool
    fold_ids: np.ndarray
    adjusted_draws: np.ndarray | None = None


def ilr_transform(composition: np.ndarray, pseudocount: float = 1e-8) -> np.ndarray:
    """Isometric log-ratio coordinates for a positive simplex matrix."""
    values = np.asarray(composition, dtype=float)
    if values.ndim != 2 or values.shape[1] < 2:
        raise ValueError("composition must have at least two cell classes")
    values = np.maximum(values, 0.0) + pseudocount
    values /= values.sum(axis=1, keepdims=True)
    log_values = np.log(values)
    centered = log_values - log_values.mean(axis=1, keepdims=True)
    n = values.shape[1]
    basis = np.zeros((n, n - 1), dtype=float)
    for j in range(n - 1):
        basis[: j + 1, j] = 1.0 / np.sqrt((j + 1) * (j + 2))
        basis[j + 1, j] = -(j + 1) / np.sqrt((j + 1) * (j + 2))
    return centered @ basis


def _fit_linear(x: np.ndarray, y: np.ndarray, ridge: float) -> np.ndarray:
    design = np.column_stack([np.ones(len(x)), x])
    penalty = np.eye(design.shape[1]) * ridge
    penalty[0, 0] = 0.0
    return np.linalg.solve(design.T @ design + penalty, design.T @ y)


def crossfit_composition_adjustment(
    field: np.ndarray,
    composition: np.ndarray,
    technical: np.ndarray | None = None,
    groups: np.ndarray | None = None,
    n_folds: int = 5,
    ridge: float = 1e-4,
    composition_draws: np.ndarray | None = None,
) -> CompositionAdjustment:
    """Return ``field - E[field | ilr(composition), technical]``.

    ``groups`` are patient/block IDs.  They are the unit of the fold split;
    spots are never treated as independent biological replicates.
    """
    y = np.asarray(field, dtype=float)
    x = ilr_transform(composition)
    if y.ndim != 1 or len(y) != len(x):
        raise ValueError("field and composition have incompatible lengths")
    if technical is not None:
        technical = np.asarray(technical, dtype=float)
        if technical.ndim == 1:
            technical = technical[:, None]
        if len(technical) != len(y):
            raise ValueError("technical covariate length differs from field")
        x = np.column_stack([x, technical])
    if groups is None:
        groups = np.arange(len(y), dtype=str)
    groups = np.asarray(groups)
    if len(groups) != len(y):
        raise ValueError("group length differs from field")
    unique = np.unique(groups)
    if len(unique) < 2:
        raise ValueError("at least two groups are required for cross-fitting")
    folds = np.arange(len(unique)) % max(2, min(n_folds, len(unique)))
    fold_by_group = dict(zip(unique.tolist(), folds.tolist()))
    fold_ids = np.asarray([fold_by_group[g] for g in groups], dtype=int)
    expected = np.zeros_like(y)
    for fold in np.unique(fold_ids):
        train = fold_ids != fold
        test = ~train
        coef = _fit_linear(x[train], y[train], ridge)
        expected[test] = np.column_stack([np.ones(test.sum()), x[test]]) @ coef
    if composition_draws is None:
        return CompositionAdjustment(y - expected, expected, False, fold_ids)
    draws = np.asarray(composition_draws, dtype=float)
    if draws.ndim != 3 or draws.shape[1:] != composition.shape:
        raise ValueError("composition_draws must have shape (n_draws, n_spots, n_types)")
    adjusted_draws = []
    expected_draws = []
    for draw in draws:
        result = crossfit_composition_adjustment(
            y, draw, technical=technical, groups=groups, n_folds=n_folds, ridge=ridge,
        )
        adjusted_draws.append(result.adjusted)
        expected_draws.append(result.expected)
    adjusted_draws_array = np.stack(adjusted_draws, axis=0)
    expected_draws_array = np.stack(expected_draws, axis=0)
    return CompositionAdjustment(
        adjusted_draws_array.mean(axis=0), expected_draws_array.mean(axis=0), True,
        fold_ids, adjusted_draws_array,
    )
