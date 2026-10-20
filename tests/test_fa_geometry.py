import numpy as np
import pytest

from field_analyzer import geometry


def grid(n=20):
    yy, xx = np.mgrid[0:n, 0:n]
    xy = np.column_stack([xx.ravel(), yy.ravel()]).astype(float)
    return xy, xx.ravel(), yy.ravel()


def test_isolated_point_gets_nan_not_its_own_value():
    xy, _, _ = grid(10)
    support = np.ones(len(xy), dtype=bool)
    adj = geometry.build_adjacency(xy, support, k=4)
    values = np.arange(len(xy), dtype=float)
    # a lone point 50 pitches away from everyone
    xy2 = np.vstack([xy, [[50.0, 50.0]]])
    support2 = np.r_[support, True]
    adj2 = geometry.build_adjacency(xy2, support2, k=4)
    out = geometry.neighbor_sum(adj2.graph, np.r_[values, 999.0], valid=support2)
    assert np.isnan(out[-1]), "isolated point must stay undefined, never copy its own value"
    assert adj2.n_isolated == 1
    assert adj2.isolated.tolist() == [len(xy2) - 1]


def test_no_edges_into_non_tissue():
    xy, xx, yy = grid(20)
    hole = ((xx - 10) ** 2 + (yy - 10) ** 2) <= 9
    support = ~hole
    adj = geometry.build_adjacency(xy, support, k=6)
    region = geometry.region_classes(support, support.copy())
    assert geometry.assert_edges_avoid_non_tissue(adj.graph, region) == 0
    coo = adj.graph.tocoo()
    assert not (region[coo.row] == geometry.NON_TISSUE).any()


def test_direction_undefined_outside_valid_mask():
    xy, xx, yy = grid(30)
    support = np.ones(len(xy), dtype=bool)
    adj = geometry.build_adjacency(xy, support, k=6)
    values = np.exp(-((xx - 15.0) ** 2 + (yy - 15.0) ** 2) / 50.0)
    valid = np.ones(len(xy), dtype=bool)
    valid[:100] = False
    res = geometry.local_direction(values, valid, xy, adj)
    assert np.isnan(res["angle_mod_pi"][~valid]).all()
    assert not res["measurable"][~valid].any()
    assert res["measurable"][valid].mean() > 0.9


def test_distance_pairs_respect_support():
    xy, xx, yy = grid(15)
    support = np.ones(len(xy), dtype=bool)
    support[0] = False
    pairs = geometry.distance_pairs(xy, support, max_distance=2.5, n_bins=5)
    assert pairs["counts"].sum() == len(pairs["distances"])
    assert 0 not in set(pairs["pairs"].ravel().tolist())


def test_pitch_estimate_is_empirical():
    xy, _, _ = grid(10)
    pitch = geometry.estimate_pitch(xy)
    assert pitch == pytest.approx(1.0)


def test_duplicate_coordinates_refused():
    xy = np.array([[0.0, 0.0], [0.0, 0.0], [1.0, 0.0]])
    with pytest.raises(Exception):
        geometry.estimate_pitch(xy)
