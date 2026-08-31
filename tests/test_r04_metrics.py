from __future__ import annotations

import numpy as np

from r04.metrics import (
    empirical_null_fdr,
    grouped_bootstrap_mean,
    grouped_delta,
    negative_binomial_log_likelihood,
    posterior_spatial_variance_probability,
)


def test_metrics_use_groups_and_not_spot_counts() -> None:
    mean, low, high = grouped_bootstrap_mean({"P1": 1.0, "P2": 2.0, "P3": 3.0}, draws=100, seed=1)
    assert mean == 2.0 and low <= mean <= high
    assert grouped_delta({"P1": 3, "P2": 4}, {"P1": 1, "P2": 2}) == {"P1": 2.0, "P2": 2.0}


def test_spatial_posterior_and_null_are_continuous_field_metrics() -> None:
    draws = np.asarray([[0, 0, 1], [0, 0, 2], [0, 0, 3]], float)
    assert posterior_spatial_variance_probability(draws) == 1.0
    assert empirical_null_fdr(2.0, [0.1, 0.2, 1.0]) == 0.25


def test_negative_binomial_metric_is_finite_and_gene_aligned() -> None:
    counts = np.asarray([[2, 4], [3, 1]], dtype=float)
    better = negative_binomial_log_likelihood(counts, np.log(np.asarray([[2.0, 4.0], [3.0, 1.0]])), np.asarray([3.0, 3.0]))
    worse = negative_binomial_log_likelihood(counts, np.log(np.ones((2, 2))), np.asarray([3.0, 3.0]))
    assert np.isfinite(better)
    assert better > worse
