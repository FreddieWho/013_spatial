"""Unit tests for r16/joint.py (pure functions, synthetic data)."""

import numpy as np
from scipy import sparse

from r16 import joint as J


def test_pool_log1p_shapes():
    rng = np.random.default_rng(0)
    mats = [sparse.csr_matrix(rng.poisson(2, size=(30, 50))),
            sparse.csr_matrix(rng.poisson(2, size=(40, 50)))]
    x = J.pool_log1p(mats)
    assert x.shape == (70, 50) and x.dtype == np.float32
    assert np.allclose(x, np.log1p(np.vstack([m.toarray() for m in mats])))


def test_neighbor_mean_excludes_self_and_matches_manual():
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [10.0, 10.0]])
    x = np.arange(8, dtype=float).reshape(4, 2)
    m = J.neighbor_mean_log1p(x, coords, k=2)
    # spot 0 neighbors are spots 1,2 (self excluded)
    assert np.allclose(m[0], (x[1] + x[2]) / 2)
    assert m.shape == x.shape


def test_augment_blocks_weighting():
    rng = np.random.default_rng(1)
    a = rng.normal(size=(10, 6))
    b = rng.normal(size=(10, 6))
    out = J.augment_blocks(a, b, 0.5)
    assert out.shape == (10, 12)
    w = float(np.sqrt(0.5))
    assert np.allclose(out, np.hstack([w * a, w * b]))
    out0 = J.augment_blocks(a, b, 0.0)
    assert np.allclose(out0[:, :6], a) and (out0[:, 6:] == 0).all()


def test_patient_halves_deterministic_and_disjoint():
    pats = [f"P{i:02d}" for i in range(30)]
    a1, b1 = J.patient_halves(pats)
    a2, b2 = J.patient_halves(pats)
    assert a1 == a2 and b1 == b2
    assert len(a1) == 15 and len(b1) == 15 and not (a1 & b1)
    assert (a1 | b1) == set(pats)


def test_best_match_table_identity_and_orthogonal():
    rng = np.random.default_rng(2)
    A = rng.normal(size=(5, 8))
    best = J.best_match_table(A, A.copy())
    assert np.allclose(best, 1.0)
    B = np.zeros((3, 8))
    best0 = J.best_match_table(A, B)
    assert (best0 == 0.0).all()


def test_combat_removes_batch_shift():
    rng = np.random.default_rng(3)
    base = rng.normal(size=(60, 20))
    x = np.vstack([base[:30] + 5.0, base[30:] - 5.0]).astype(np.float32)
    batches = np.array(["s1"] * 30 + ["s2"] * 30)
    xc = J.combat_correct(x, batches)
    assert xc.shape == x.shape
    assert np.all(np.isfinite(xc))
    # batch means should be much closer after correction, per gene
    gap_before = np.abs(x[:30].mean(axis=0) - x[30:].mean(axis=0)).mean()
    gap_after = np.abs(xc[:30].mean(axis=0) - xc[30:].mean(axis=0)).mean()
    assert gap_after < 0.2 * gap_before
    # deterministic
    xc2 = J.combat_correct(x, batches)
    assert np.array_equal(xc, xc2)
    # within-batch structure preserved (rank correlation of a gene)
    from scipy.stats import spearmanr
    r = spearmanr(x[:30, 0], xc[:30, 0]).statistic
    assert r > 0.9
