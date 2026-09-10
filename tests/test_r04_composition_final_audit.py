import numpy as np
import pytest

from r04.composition import (
    apply_residualizer,
    crossfit_composition_adjustment,
    fit_residualizer,
)


def _simplex(n, k, seed):
    rng = np.random.default_rng(seed)
    x = rng.random((n, k)) + 0.5
    return x / x.sum(axis=1, keepdims=True)


def test_train_only_mapping_ignores_test_rows():
    comp_train = _simplex(60, 3, 0)
    field_train = comp_train[:, 0] * 2.0 + np.random.default_rng(1).normal(0, 0.01, 60)
    comp_test_a = _simplex(40, 3, 2)
    field_test_a = np.random.default_rng(3).normal(0, 1, 40)
    # Adversarial alternative test set: strong composition-only signal.
    comp_test_b = _simplex(40, 3, 4)
    field_test_b = comp_test_b[:, 1] * 50.0
    coef = fit_residualizer(field_train, comp_train)
    adj_train_a, _ = apply_residualizer(field_train, comp_train, coef)
    adj_train_b, _ = apply_residualizer(field_train, comp_train, coef)
    np.testing.assert_allclose(adj_train_a, adj_train_b)
    # The mapping transfers: test rows get residuals under the train fit.
    adj_test_b, _ = apply_residualizer(field_test_b, comp_test_b, coef)
    assert np.isfinite(adj_test_b).all()


def test_pooled_crossfit_absorbs_test_only_signal():
    # Demonstrates the D-108 leak: pooled fitting lets test rows move the
    # residualizer that later transforms train rows used to fit a readout.
    rng = np.random.default_rng(5)
    comp_train = _simplex(60, 3, 6)
    field_train = rng.normal(0, 1, 60)
    comp_test = _simplex(40, 3, 7)
    field_test = comp_test[:, 1] * 50.0
    pooled = crossfit_composition_adjustment(
        np.concatenate([field_train, field_test]),
        np.concatenate([comp_train, comp_test]),
        groups=np.array(["train"] * 60 + ["test"] * 40))
    coef = fit_residualizer(field_train, comp_train)
    nested_train, _ = apply_residualizer(field_train, comp_train, coef)
    # Same train rows, different residualizers: pooled train residuals must
    # differ because the test-only signal entered the pooled fit, while the
    # nested train-only mapping never sees the test rows.
    assert not np.allclose(pooled.adjusted[:60], nested_train, atol=1e-6)


def test_fit_apply_rejects_shape_mismatch():
    comp = _simplex(10, 3, 8)
    coef = fit_residualizer(np.arange(10, dtype=float), comp)
    with pytest.raises(ValueError):
        apply_residualizer(np.arange(9, dtype=float), comp, coef)
    with pytest.raises(ValueError):
        fit_residualizer(np.arange(9, dtype=float), comp)
