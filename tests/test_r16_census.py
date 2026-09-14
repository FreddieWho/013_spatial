import numpy as np
import pytest
from scipy import sparse

from r16 import census as C


def grid_coords(n_side=20):
    xs, ys = np.meshgrid(np.arange(n_side), np.arange(n_side))
    return np.column_stack([xs.ravel(), ys.ravel()])


def test_moran_permutation_detects_clustered_vs_random():
    rng = np.random.default_rng(0)
    coords = grid_coords()
    w = C.spatial_weights(coords, k=4)
    clustered = (coords[:, 0] < 10).astype(float)
    _, p_cl = C.moran_permutation_p(clustered, w, 99, rng)
    assert p_cl <= 0.01
    random_mask = rng.permutation(clustered)
    _, p_rd = C.moran_permutation_p(random_mask, w, 99, rng)
    assert p_rd > 0.01


def test_filter_min_size_and_jaccard():
    labels = np.array([0] * 25 + [1] * 5 + [2] * 30)
    out = C.filter_min_size(labels, 20)
    assert (out[25:30] == -1).all()
    assert (out[:25] == 0).all() and (out[30:] == 2).all()
    assert C.jaccard(np.array([1, 1, 0]), np.array([1, 1, 0])) == 1.0
    assert C.jaccard(np.array([1, 1, 0]), np.array([0, 0, 1])) == 0.0


def test_signatures_and_cosine():
    x = np.vstack([np.tile([1.0, -1.0], (10, 1)), np.tile([-1.0, 1.0], (10, 1))])
    labels = np.array([0] * 10 + [1] * 10)
    sigs = C.cluster_signatures(x, labels)
    assert set(sigs) == {0, 1}
    assert C.cosine(sigs[0], sigs[0]) == pytest.approx(1.0)
    assert C.cosine(sigs[0], sigs[1]) == pytest.approx(-1.0)


def test_leiden_deterministic_on_two_blobs():
    rng = np.random.default_rng(1)
    a = rng.normal(0, 0.1, size=(60, 5))
    b = rng.normal(5, 0.1, size=(60, 5))
    scores = np.vstack([a, b])
    g = C.knn_igraph(scores, 10)
    l1 = C.leiden_labels(g, 0.5)
    l2 = C.leiden_labels(g, 0.5)
    assert (l1 == l2).all()
    # two blobs should be separable: adjusted overlap of label sets with truth
    truth = np.array([0] * 60 + [1] * 60)
    best = max(C.jaccard(truth == 0, l1 == c) for c in np.unique(l1))
    assert best > 0.9
