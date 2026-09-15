import numpy as np
import pytest
from scipy import sparse

from r16 import axis_factory as F
from r16 import census as C


def test_spatial_mean_weights_excludes_self_and_rowsum_one():
    rng = np.random.default_rng(0)
    coords = rng.normal(size=(50, 2))
    w = F.spatial_mean_weights(coords, k=5)
    assert w.shape == (50, 50)
    assert w.diagonal().sum() == 0.0
    assert np.allclose(np.asarray(w.sum(axis=1)).ravel(), 1.0)
    assert np.allclose(w.data, 0.2)


def test_augment_endpoints():
    rng = np.random.default_rng(1)
    x = rng.normal(size=(20, 5))
    m = rng.normal(size=(20, 5))
    assert np.allclose(F.augment(x, m, 0.0), x)
    assert np.allclose(F.augment(x, m, 1.0), m)
    mid = F.augment(x, m, 0.5)
    assert mid.shape == x.shape and np.isfinite(mid).all()


def test_candidate_pcs_shapes_and_determinism():
    rng = np.random.default_rng(2)
    x = rng.normal(size=(60, 40))
    s1, l1, v1 = F.candidate_pcs(x, n_pc=10)
    s2, l2, v2 = F.candidate_pcs(x, n_pc=10)
    assert s1.shape == (60, 10) and l1.shape == (40, 10)
    assert np.allclose(np.abs(s1), np.abs(s2))  # sign ambiguity allowed
    assert float(v1.sum()) <= 1.0 + 1e-9


def _grid(n_side=16):
    xs, ys = np.meshgrid(np.arange(n_side), np.arange(n_side))
    return np.column_stack([xs.ravel(), ys.ravel()]).astype(float)


def test_lisa_finds_blob_rejects_noise():
    rng = np.random.default_rng(3)
    coords = _grid()
    w = C.spatial_weights(coords, C.KNN_SPATIAL)
    blob = ((coords[:, 0] - 8) ** 2 + (coords[:, 1] - 8) ** 2 < 16).astype(float)
    out = F.lisa_tile_stats(blob, w, 499, rng)
    assert out["tile_ok"] and out["max_hotspot"] >= 20
    noise = rng.normal(size=len(coords))
    out_n = F.lisa_tile_stats(noise, w, 499, rng)
    assert not out_n["tile_ok"]


def test_bh_fdr_basic():
    p = np.array([0.001, 0.002, 0.5, 0.9])
    rej = F.bh_fdr(p, alpha=0.05)
    assert rej.tolist() == [True, True, False, False]
    assert not F.bh_fdr(np.array([0.5, 0.9]), alpha=0.05).any()


def test_abs_cosine_sign_invariance():
    rng = np.random.default_rng(4)
    a = rng.normal(size=50)
    assert F.abs_cosine(a, -a) == pytest.approx(1.0)
    assert F.abs_cosine(a, a) == pytest.approx(1.0)
    assert F.abs_cosine(np.zeros(50), a) == 0.0
