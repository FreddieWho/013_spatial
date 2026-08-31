import numpy as np
import pytest

from r04.checkpoint_audit import nuisance_gauge_summary
from scripts.r04_checkpoint_audit import _audit_payload


def test_nuisance_gauge_shift_preserves_additive_rate() -> None:
    baseline = np.asarray([2.0, 3.0, 5.0])
    loading = np.asarray([[0.2], [0.3], [0.5]])
    nuisance = np.asarray([[1.5], [2.0], [3.0]])
    shift = 1.0
    before = baseline[None, :] + nuisance @ loading.T
    after = (
        (baseline + shift * loading[:, 0])[None, :]
        + (nuisance - shift) @ loading.T
    )
    np.testing.assert_allclose(before, after)

    summary = nuisance_gauge_summary(
        np.log(np.expm1(baseline)),
        np.log(loading),
        np.log(np.expm1(nuisance)),
    )
    assert summary["rank"] == 1
    assert summary["h_min"][0] == pytest.approx(1.50001, abs=1e-6)


def test_checkpoint_audit_payload_records_ridge() -> None:
    payload = _audit_payload(
        manifest_hash="manifest",
        gene_count=4000,
        factors=3,
        fold=4,
        folds=5,
        group_seed=20260807,
        inducing_points=16,
        lengthscale=3.0,
        ridge=1e-5,
        mc_draws=4,
        checkpoints=[],
    )

    assert payload["schema"] == "r04.checkpoint_audit.v1"
    assert payload["ridge"] == 1e-5
