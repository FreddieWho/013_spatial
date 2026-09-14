import numpy as np
import pytest
from scipy import sparse

from r16 import axes as A


def line_graph(n=10):
    rows = np.arange(n - 1)
    w = sparse.csr_matrix((np.ones(2 * (n - 1)),
                           (np.r_[rows, rows + 1], np.r_[rows + 1, rows])),
                          shape=(n, n))
    return w


def test_signed_hop_distance_on_line():
    w = line_graph(10)
    mask = np.array([True] * 3 + [False] * 7)
    signed = A.signed_hop_distance(w, mask)
    assert signed[2] == -1 and signed[0] == -3      # inside, distance to outside
    assert signed[3] == 1 and signed[9] == 7        # outside, distance to mask


def test_signed_hop_distance_degenerate_masks():
    w = line_graph(6)
    assert (A.signed_hop_distance(w, np.ones(6, bool)) == 0).all()
    assert (A.signed_hop_distance(w, np.zeros(6, bool)) == 0).all()


def test_axis_scores_missing_genes():
    x = np.ones((5, 4))
    s, n = A.axis_scores(x, {"A": 0, "B": 1}, ["A", "MISSING"])
    assert n == 1 and np.allclose(s, 1.0)
    s0, n0 = A.axis_scores(x, {"A": 0}, ["NOPE"])
    assert n0 == 0 and (s0 == 0).all()


def test_contour_stability_bounds_and_monotonic():
    rng = np.random.default_rng(0)
    s = rng.normal(size=500)
    stab, js = A.contour_stability(s)
    assert 0 < stab <= 1 and len(js) == 2
    m6 = A.contour_mask(s, 0.6)
    m8 = A.contour_mask(s, 0.8)
    assert m8.sum() < m6.sum() and (m8 <= m6).all()


def test_binned_profile():
    vals = np.array([1.0, 2.0, 3.0, 4.0])
    signed = np.array([-1, -1, 1, 2])
    prof = A.binned_profile(vals, signed, np.array([-1, 0, 1, 2]))
    assert prof[0] == (pytest.approx(1.5), 2)
    assert np.isnan(prof[1][0]) and prof[1][1] == 0
    assert prof[3] == (pytest.approx(4.0), 1)
