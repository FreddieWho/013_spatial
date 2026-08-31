"""Patient/block-level metrics; spots are not biological replicates."""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
from scipy.special import gammaln


def negative_binomial_log_likelihood(
    counts: np.ndarray,
    log_mu: np.ndarray,
    dispersion: np.ndarray,
) -> float:
    """Return mean per-spot NB log likelihood for a fitted count mean."""
    y = np.asarray(counts, dtype=float)
    mean_log = np.asarray(log_mu, dtype=float)
    theta = np.asarray(dispersion, dtype=float)
    if y.shape != mean_log.shape or y.ndim != 2:
        raise ValueError("counts and log_mu must have equal two-dimensional shape")
    if theta.ndim != 1 or theta.shape[0] != y.shape[1] or np.any(theta <= 0):
        raise ValueError("dispersion must be positive and gene-aligned")
    mu = np.exp(np.clip(mean_log, -40.0, 40.0))
    denominator = theta[None, :] + mu
    log_probability = (
        gammaln(y + theta[None, :])
        - gammaln(theta[None, :])
        - gammaln(y + 1.0)
        + theta[None, :] * (np.log(theta)[None, :] - np.log(denominator))
        + y * (mean_log - np.log(denominator))
    )
    return float(np.mean(np.sum(log_probability, axis=1)))


def grouped_bootstrap_mean(values_by_group: Mapping[str, float], draws: int = 1000, seed: int = 20260807) -> tuple[float, float, float]:
    if len(values_by_group) < 2:
        raise ValueError("at least two groups are required for bootstrap")
    values = np.asarray(list(values_by_group.values()), dtype=float)
    rng = np.random.default_rng(seed)
    sampled = values[rng.integers(0, len(values), size=(draws, len(values)))]
    distribution = sampled.mean(axis=1)
    return float(values.mean()), float(np.quantile(distribution, 0.025)), float(np.quantile(distribution, 0.975))


def grouped_delta(full_by_group: Mapping[str, float], baseline_by_group: Mapping[str, float]) -> dict[str, float]:
    groups = set(full_by_group) & set(baseline_by_group)
    if len(groups) < 2:
        raise ValueError("full and baseline require two common biological groups")
    return {group: float(full_by_group[group] - baseline_by_group[group]) for group in sorted(groups)}


def posterior_spatial_variance_probability(draws: np.ndarray, tolerance: float = 1e-8) -> float:
    values = np.asarray(draws, dtype=float)
    if values.ndim != 2:
        raise ValueError("posterior field draws must have shape (draw, spot)")
    return float(np.mean(np.var(values, axis=1) > tolerance))


def empirical_null_fdr(observed: float, null_values: Sequence[float]) -> float:
    null = np.asarray(list(null_values), dtype=float)
    if len(null) == 0:
        raise ValueError("spatial null distribution cannot be empty")
    return float((1 + np.sum(null >= observed)) / (len(null) + 1))
