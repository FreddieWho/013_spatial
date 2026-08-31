from __future__ import annotations

import numpy as np
import pytest

from r04.diagnostics import gene_batch_mean_scale, gene_batch_sum_scale
from r04.objectives import gene_batch_sum_objective
from r04.models import MNSFConfig, MNSFEstimator, SignedResidualGPConfig, SignedResidualGPEstimator
from r04.synthetic import make_overlapping_sections


def test_gene_batch_horvitz_thompson_estimator_targets_dense_mean() -> None:
    rng = np.random.default_rng(5)
    values = rng.normal(size=(37, 23))
    dense = float(values.mean())
    estimates = []
    for _ in range(2000):
        index = rng.choice(values.shape[1], size=7, replace=False)
        estimates.append(gene_batch_mean_scale(values[:, index].mean(), values.shape[1], len(index)))
    assert float(np.mean(estimates)) == pytest.approx(dense, abs=0.04)


def test_gene_batch_sum_scaling_targets_dense_per_spot_gene_sum() -> None:
    rng = np.random.default_rng(15)
    values = rng.uniform(0.1, 2.0, size=(37, 23))
    dense = float(values.sum(axis=1).mean())
    estimates = []
    for _ in range(2000):
        index = rng.choice(values.shape[1], size=7, replace=False)
        batch_sum = float(values[:, index].sum(axis=1).mean())
        estimates.append(gene_batch_sum_scale(batch_sum, values.shape[1], len(index)))
    assert float(np.mean(estimates)) == pytest.approx(dense, abs=0.05)


def test_tensorflow_dense_and_gene_batch_objectives_are_unbiased() -> None:
    tf = pytest.importorskip("tensorflow")
    values = tf.constant(np.arange(1, 1 + 12 * 9, dtype=np.float32).reshape(12, 9))
    dense = float(tf.reduce_mean(tf.reduce_sum(values, axis=1)).numpy())
    estimates = []
    for index in (np.array([0, 2, 5]), np.array([1, 4, 7]), np.array([3, 6, 8])):
        estimates.append(float(gene_batch_sum_objective(tf.gather(values, index, axis=1), 9, 3).numpy()))
    assert float(np.mean(estimates)) == pytest.approx(dense, abs=1e-6)


def test_tensorflow_full_batch_objective_has_matching_gradient() -> None:
    tf = pytest.importorskip("tensorflow")
    variable = tf.Variable(np.linspace(0.2, 1.1, 18, dtype=np.float32).reshape(6, 3))
    with tf.GradientTape() as tape:
        dense = gene_batch_sum_objective(variable, 3, 3)
    gradient = tape.gradient(dense, variable)
    expected = np.ones((6, 3), dtype=np.float32) / 6.0
    np.testing.assert_allclose(gradient.numpy(), expected, rtol=1e-6, atol=1e-6)


def test_gene_minibatch_smoke_records_unscaled_global_kl() -> None:
    sections, _ = make_overlapping_sections(sections=2, spots_per_section=16, genes=12, antagonistic=True)
    mnsf = MNSFEstimator(MNSFConfig(factors=2, inducing_points=5, steps=2, posterior_draws=2, gene_batch_size=5, seed=43)).fit(sections)
    signed = SignedResidualGPEstimator(SignedResidualGPConfig(factors=2, inducing_points=5, steps=2, posterior_draws=2, gene_batch_size=5, seed=47)).fit(sections)
    for estimator in (mnsf, signed):
        assert estimator.diagnostics_["gene_batch_size"] == 5
        assert estimator.diagnostics_["gene_batch_scaling"] == "uniform_without_replacement_ht"
        assert estimator.diagnostics_["global_spatial_kl_scaling"] == "unscaled_once_per_step"
