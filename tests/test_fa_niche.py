import numpy as np
from scipy import sparse

from field_analyzer import niche_graph as ng


def toy():
    # 3x3 grid of cells, 3 niches in blocks
    xy = np.array([[x, y] for y in range(3) for x in range(3)], dtype=float)
    labels = np.array([0, 0, 0, 0, 1, 1, 1, 2, 2])
    graph = ng.knn_adjacency(xy, k=2)
    return xy, labels, graph


def test_contact_counts_are_symmetric():
    xy, labels, graph = toy()
    C, row, total = ng.contact_counts(labels, graph, 3)
    assert np.max(np.abs(C - C.T)) == 0.0
    assert row.sum() == C.sum()


def test_enrichment_reports_symmetry_error_and_diagonal_enrichment():
    xy, labels, graph = toy()
    res = ng.contact_enrichment(labels, graph, 3)
    assert res["counts_symmetry_error"] == 0.0
    # same-niche contact should be the most enriched
    diag = np.diag(res["enrichment"])
    off = res["enrichment"][~np.eye(3, dtype=bool)]
    assert np.nanmax(diag) >= np.nanmax(off)


def test_permutation_null_is_much_smaller_than_blocked_labels():
    xy, labels, graph = toy()
    obs = ng.contact_enrichment(labels, graph, 3)["statistic"]
    perm = [ng.contact_enrichment(ng.permutation_null_labels(labels, seed=i), graph, 3)["statistic"]
            for i in range(10)]
    assert obs > float(np.median(perm))


def test_rotation_null_keeps_every_label_count():
    xy, labels, _ = toy()
    rotated = ng.rotation_null_labels(xy, labels, angle_deg=90.0)
    assert np.array_equal(np.sort(rotated), np.sort(labels))


def test_kmeans_returns_requested_number_of_labels():
    rng = np.random.default_rng(0)
    x = np.vstack([rng.normal(0, 0.1, (30, 3)), rng.normal(5, 0.1, (30, 3))])
    labels = ng.kmeans_niches(x, k=2, seed=0)
    assert len(set(labels.tolist())) == 2


def test_argmax_niches_is_label_consistent_across_sections():
    from field_analyzer.niche_graph import argmax_niches

    names = ["A", "B"]
    s1 = np.array([[2.0, 0.0], [0.0, 2.0], [0.1, 0.1]])
    s2 = np.array([[1.0, 0.2], [0.2, 1.0], [0.0, 0.0]])
    l1, n1 = argmax_niches(s1, names, mix_threshold=0.5)
    l2, n2 = argmax_niches(s2, names, mix_threshold=0.5)
    assert n1 == n2
    assert l1[0] == l2[0] == 0      # both call the first row "A"
    assert l1[1] == l2[1] == 1
    assert l1[2] == n1.index("MIXED") and l2[2] == n2.index("MIXED")
