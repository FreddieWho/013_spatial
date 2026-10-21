import numpy as np
import pytest
from scipy import sparse

from field_analyzer import factor_model as fm


def make_section(n_spots=200, n_genes=120, seed=0, signal=1.0, scale=20.0):
    rng = np.random.default_rng(seed)
    xy = np.column_stack([rng.uniform(0, 40, n_spots), rng.uniform(0, 40, n_spots)])
    depth = scale * np.exp(rng.normal(0, 0.2, n_spots))
    # one smooth latent factor drives a block of genes
    f = np.exp(-((xy[:, 0] - 20) ** 2 + (xy[:, 1] - 20) ** 2) / 200.0)
    mean = depth[:, None] * np.exp(signal * np.outer(f, rng.normal(0, 1, n_genes)))
    counts = sparse.csr_matrix(rng.poisson(mean).astype(float))
    return fm.SectionData(
        section_id="s", unit_id="u", counts=counts,
        genes=np.array([f"g{i}" for i in range(n_genes)]),
        depth=depth, xy=xy, visible=np.ones(n_spots, dtype=bool),
    )


def test_depth_offset_uses_only_non_modelled_genes():
    rng = np.random.default_rng(0)
    all_counts = sparse.csr_matrix(rng.poisson(3.0, size=(50, 200)).astype(float))
    modelled = np.zeros(200, dtype=bool)
    modelled[:40] = True
    depth = fm.estimate_depth(all_counts, modelled)
    lib = np.asarray(all_counts[:, ~modelled].sum(axis=1)).ravel()
    assert np.allclose(depth, lib / np.median(lib[lib > 0]))


def test_fit_section_shapes_and_finiteness():
    d = make_section()
    fit = fm.fit_section(d, None, n_factors=3, n_smooth=2, n_iter=6)
    assert fit.factors.shape == (3, d.counts.shape[0])
    assert fit.loadings.shape == (d.counts.shape[1], 3)
    assert np.isfinite(fit.poisson_deviance)
    assert np.all(fit.loadings >= 0)


def test_spatial_smoothing_reduces_graph_roughness():
    """Roughness must be measured on the spatial graph, not on the storage order."""
    d = make_section()
    raw = fm.fit_section(d, None, n_factors=2, n_smooth=0, n_iter=6)
    smooth = fm.fit_section(d, None, n_factors=2, n_smooth=6, n_iter=6)
    smoother = fm._graph_smoother(d.xy, k=6)
    adj = (smoother > 0).astype(float)
    def roughness(f):
        diff = f[:, :, None] - f[:, None, :]
        d2 = np.asarray(diff ** 2).mean(axis=0)
        n_edges = adj.sum()
        return float((d2 * np.asarray(adj.todense())).sum() / max(n_edges, 1))
    assert roughness(smooth.factors) <= roughness(raw.factors) + 1e-9


def test_local_gene_knn_beats_depth_only_on_a_smooth_field():
    d = make_section(seed=1, signal=2.0, scale=40.0)
    rng = np.random.default_rng(1)
    holdout = rng.random(d.counts.shape[0]) < 0.2
    base = fm.depth_only_loglik(d, holdout=holdout)
    knn = fm.local_gene_knn_loglik(d, holdout=holdout, k=6)
    assert knn["deviance"] < base["deviance"]


def test_heldout_loglik_only_scores_heldout_spots():
    d = make_section(seed=2)
    rng = np.random.default_rng(2)
    holdout = rng.random(d.counts.shape[0]) < 0.1
    lam, _ = fm.fit_shared_loadings([d], n_factors=2, n_outer=1, n_iter=5, n_smooth=1)
    res = fm.heldout_spot_loglik(d, lam, holdout=holdout, n_factors=2, n_smooth=1, n_iter=5)
    assert res["n_holdout"] == int(holdout.sum())
    assert np.isfinite(res["deviance"])
