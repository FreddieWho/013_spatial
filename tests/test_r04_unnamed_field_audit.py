import numpy as np

from scripts.r04_unnamed_field_audit import (
    _canonical_corrs,
    _subspace_basis,
    _verdict,
)


def _random_loading(seed, genes=60):
    rng = np.random.default_rng(seed)
    return rng.normal(size=(genes, 3))


def test_rotation_sign_permutation_invariance():
    rng = np.random.default_rng(0)
    base = _random_loading(1)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    rotated = (base @ q)[:, [2, 0, 1]] * np.array([-1.0, 1.0, -1.0])
    for rank in (1, 2, 3):
        cc = _canonical_corrs(_subspace_basis(base, rank), _subspace_basis(rotated, rank))
        assert all(v > 0.999 for v in cc)


def test_verdict_thresholds():
    assert _verdict(0.01, 0.05) == "DOES_NOT_SURVIVE"
    assert _verdict(0.60, 0.05) == "SURVIVES"
    assert _verdict(0.30, 0.05) == "NOT_IDENTIFIABLE"
    assert _verdict(0.04, 0.05) == "DOES_NOT_SURVIVE"


def test_shared_direction_detected_and_noise_rejected(tmp_path):
    rng = np.random.default_rng(2)
    shared = rng.normal(size=(200, 1)) * 5.0
    load_a = np.column_stack([shared.ravel(), rng.normal(size=(200, 2))])
    load_b = np.column_stack([shared.ravel(), rng.normal(size=(200, 2))])
    cc1 = _canonical_corrs(_subspace_basis(load_a, 1), _subspace_basis(load_b, 1))
    assert cc1[0] > 0.9
    assert _verdict(cc1[0], 0.2) == "SURVIVES"
    noise_b = rng.normal(size=(200, 3))
    cc_noise = _canonical_corrs(_subspace_basis(load_a, 1), _subspace_basis(noise_b, 1))
    assert _verdict(cc_noise[0], 0.2) == "DOES_NOT_SURVIVE"
