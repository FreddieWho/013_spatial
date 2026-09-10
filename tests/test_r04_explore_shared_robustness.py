import numpy as np
import pytest

from scripts.r04_explore_shared_robustness import (
    cosine,
    recover_gene_projection,
    shared_tls_auc,
    top_overlap,
)


def _synthetic(n_per_patient=80, seed=0):
    rng = np.random.default_rng(seed)
    patients = np.asarray(["p1"] * n_per_patient + ["p2"] * n_per_patient)
    signal = np.concatenate([np.zeros(n_per_patient), np.ones(n_per_patient)])
    z = np.column_stack([signal + 0.3 * rng.standard_normal(2 * n_per_patient),
                         rng.standard_normal(2 * n_per_patient)])
    order = rng.permutation(2 * n_per_patient)
    base = np.concatenate([np.zeros(60), np.ones(20), np.zeros(20), np.ones(60)])
    labels = {
        "TLS": base[order],
        "TUMOR_STROMA_BOUNDARY": base[rng.permutation(2 * n_per_patient)],
    }
    return z, labels, patients


def test_recover_gene_projection_is_exact():
    rng = np.random.default_rng(0)
    xc = rng.standard_normal((200, 20))
    true_p, _ = np.linalg.qr(rng.standard_normal((20, 2)))
    z = xc @ true_p
    p_hat = recover_gene_projection(xc, z)
    assert np.allclose(p_hat, true_p, rtol=1e-9, atol=1e-9)


def test_recover_gene_projection_rejects_rank_deficiency():
    xc = np.ones((50, 8))
    with pytest.raises(ValueError):
        recover_gene_projection(xc, np.ones((50, 3)))


def test_cosine_and_overlap_edges():
    assert cosine(np.array([1.0, 0.0]), np.array([0.0, 1.0])) == pytest.approx(0.0)
    assert cosine(np.array([1.0, 1.0]), np.array([2.0, 2.0])) == pytest.approx(1.0)
    assert cosine(np.zeros(3), np.ones(3)) is None
    assert top_overlap(np.array([3.0, 1.0, 2.0]), np.array([3.0, 1.0, 2.0]), 2) == 1.0
    assert 0.0 <= top_overlap(np.array([3.0, 1.0, 2.0]), np.array([1.0, 2.0, 3.0]), 1) <= 1.0


def test_shared_tls_auc_runs_and_bounds():
    z, labels, patients = _synthetic()
    train = patients == "p1"
    test = patients == "p2"
    w_train = np.full(int(train.sum()), 1.0 / train.sum())
    w_test = np.full(int(test.sum()), 1.0 / test.sum())
    auc, beta, scale = shared_tls_auc(
        z[train], z[test],
        {n: labels[n][train] for n in labels},
        {n: labels[n][test] for n in labels},
        w_train, w_test,
    )
    assert 0.0 <= auc <= 1.0
    assert beta.shape == (2,)
    assert scale.shape == (2,)
    assert np.all(scale > 0)
