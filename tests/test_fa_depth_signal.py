import numpy as np
import pytest
from scipy import sparse

from field_analyzer import depth_conditioning as dc
from field_analyzer import route_b, signal_power, synthetic as syn


def test_residualise_removes_quadratic_depth_trend():
    rng = np.random.default_rng(0)
    d = np.linspace(-2, 2, 300)
    y = 1.5 + 3.0 * d - 0.8 * d ** 2 + rng.normal(scale=1e-6, size=d.shape)
    res, r2 = dc.residualise(y, d, degree=2)
    assert r2 > 0.999
    assert np.max(np.abs(res)) < 1e-4


def test_residualise_with_constant_depth_is_a_no_op():
    y = np.arange(20.0)
    res, r2 = dc.residualise(y, np.zeros(20), degree=2)
    assert np.allclose(res, y - y.mean())
    assert r2 == pytest.approx(0.0, abs=1e-9)


def test_cross_correlation_peak_finds_known_shift():
    n = syn.GRID
    xy = syn.grid_xy(n)
    base = syn.gaussian(xy, (24.0, 32.0), 4.0)
    shifted = syn.gaussian(xy, (24.0 + 6.0, 32.0), 4.0)
    peak = route_b.cross_correlation_peak(base, shifted, n, max_lag=12)
    assert abs(peak["magnitude"] - 6.0) <= 1.0
    assert abs(abs(peak["dx"]) - 6.0) <= 1.0


def test_cross_correlation_peak_is_zero_for_identical_fields():
    n = syn.GRID
    f = syn.gaussian(syn.grid_xy(n), (30.0, 30.0), 5.0)
    peak = route_b.cross_correlation_peak(f, f, n, max_lag=8)
    assert peak["dx"] == 0 and peak["dy"] == 0


def test_markers_per_field_changes_the_gene_axis():
    s = syn.build_s2(0)
    default = len(s.gene_names())
    s.markers_per_field = 40
    wider = len(s.gene_names())
    assert wider - default == 3 * (40 - syn.MARKERS_PER_FIELD)
    counts = syn.counts_from_fields(s, seed=0)
    assert counts.shape[1] == wider


def test_score_fields_returns_one_score_per_field():
    s = syn.build_s2(0)
    s.markers_per_field = 24
    scores = signal_power.score_fields(s, depth_mu=2.0, seed=0)
    assert set(scores) == set(s.field_gene_sets())
    for values in scores.values():
        assert values.shape == (syn.GRID * syn.GRID,)
        assert np.isfinite(values).all()


def test_condition_section_needs_no_hidden_values():
    """Smoke: a tiny fake table runs end-to-end and reports depth R^2 values."""
    from types import SimpleNamespace

    rng = np.random.default_rng(0)
    n_obs, n_genes = 200, 220
    counts = sparse.csr_matrix(rng.poisson(4.0, size=(n_obs, n_genes)).astype(float))
    genes = np.array([f"g{i}" for i in range(n_genes)])
    xy = np.column_stack([np.repeat(np.arange(20), 10), np.tile(np.arange(10), 20)]).astype(float)
    table = SimpleNamespace(
        counts=counts, genes=genes, xy=xy, aux={"in_tissue": np.ones(n_obs, dtype=int)}, n_obs=n_obs
    )
    fields = {"f1": [f"g{i}" for i in range(10)], "f2": [f"g{i}" for i in range(10, 20)]}
    from field_analyzer.roles import build_field_splits

    splits = build_field_splits(fields, genes, seed=0)
    res = dc.condition_section(table, fields, splits, unit_id="u", section_id="s")
    assert set(res.map_before) == {"f1", "f2"}
    assert all(np.isfinite(v) for v in res.depth_r2_a.values())
