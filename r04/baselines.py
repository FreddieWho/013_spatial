"""Non-spatial and null baselines for R-04."""

from __future__ import annotations

import numpy as np


def library_normalized_log_counts(counts: object) -> np.ndarray:
    values = counts.toarray() if hasattr(counts, "toarray") else np.asarray(counts)
    values = np.asarray(values, dtype=float)
    library = np.maximum(values.sum(axis=1, keepdims=True), 1.0)
    return np.log1p(values / library * 10_000.0)


def nonspatial_svd(counts: object, factors: int = 4) -> tuple[np.ndarray, np.ndarray]:
    """A matched-rank non-spatial factor baseline.

    This returns continuous scores and loadings.  It never turns scores into
    spot labels, so it can be compared to the spatial field likelihood without
    making clustering the scientific object.
    """
    values = library_normalized_log_counts(counts)
    u, singular, vt = np.linalg.svd(values - values.mean(axis=0, keepdims=True), full_matrices=False)
    k = min(factors, len(singular))
    return u[:, :k] * singular[:k], vt[:k].T


def coordinate_permutation(coords: np.ndarray, values: np.ndarray, seed: int) -> np.ndarray:
    """Coordinate null preserving values and destroying spatial assignment."""
    rng = np.random.default_rng(seed)
    permutation = rng.permutation(len(coords))
    return np.asarray(values)[permutation]


def explained_increment(full_score: np.ndarray, baseline_score: np.ndarray) -> float:
    """Simple held-out variance increment used by smoke tests and diagnostics."""
    full = np.asarray(full_score, dtype=float).ravel()
    baseline = np.asarray(baseline_score, dtype=float).ravel()
    if full.shape != baseline.shape:
        raise ValueError("score arrays must have equal shape")
    residual = full - baseline
    return float(np.var(full) - np.var(residual))
