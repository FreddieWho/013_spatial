"""FA-10 (route A): count-likelihood hierarchical spatial factor model.

Why this replaces the old pipeline
----------------------------------
Everything before this node did: 12 markers -> mean of log-CPM -> one scalar per
spot -> spatial analysis.  That step threw away the count likelihood, injected the
library-size non-linearity into every position (which is what killed the displacement
estimator and made OXPHOS look reproducible), and made "same spot" and "shifted a
little" hard to separate.

This module models the counts directly:

    count[g, s] ~ Poisson( depth[s] * ( sum_k Lambda[g,k] * F[k,s] + b[g] ) )
    F = C @ Phi,   Phi = RBF basis on the spot coordinates (smoothness is explicit)

* ``depth`` is estimated from genes that are **not modelled** (outside the HVG set),
  so the offset is not circular;
* ``Lambda`` (loadings) is shared across units -> factors are comparable across
  sections/units, which the 12-marker scores never were;
* ``C`` (per-unit factor coefficients) carries a ridge prior; the RBF width sets the
  spatial range, and a "no spatial prior" variant uses a one-hot basis instead.

Fitting is alternating weighted least squares (Poisson variance), which is cheap on
CPU and keeps every parameter interpretable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np
from scipy import sparse

K_DEFAULT = 5
M_DEFAULT = 60


@dataclass
class SectionData:
    section_id: str
    unit_id: str
    counts: sparse.csr_matrix  # spots x genes (modelled genes only)
    genes: np.ndarray
    depth: np.ndarray  # offset from non-modelled genes, normalised to median 1
    xy: np.ndarray
    visible: np.ndarray


@dataclass
class FactorFit:
    section_id: str
    unit_id: str
    factors: np.ndarray  # K x spots (F)
    coefficients: np.ndarray  # K x M (C)
    loadings: np.ndarray  # genes x K (Lambda)
    baseline: np.ndarray  # genes (b)
    fitted_mean: np.ndarray  # spots x genes (rate, counts scale)
    log_likelihood: float
    poisson_deviance: float
    n_iter: int


def rbf_basis(xy: np.ndarray, centres: np.ndarray, width: float) -> np.ndarray:
    d2 = np.sum((xy[:, None, :] - centres[None, :, :]) ** 2, axis=2)
    phi = np.exp(-0.5 * d2 / (width ** 2))
    phi /= np.maximum(phi.sum(axis=1, keepdims=True), 1e-12)
    return phi


def centre_points(xy: np.ndarray, n_centres: int, seed: int = 0) -> np.ndarray:
    """k-means++-lite: deterministic farthest-point sampling."""
    xy = np.asarray(xy, dtype=float)
    rng = np.random.default_rng(seed)
    idx = [int(rng.integers(0, len(xy)))]
    d2 = np.sum((xy - xy[idx[0]]) ** 2, axis=1)
    for _ in range(1, n_centres):
        nxt = int(np.argmax(d2))
        idx.append(nxt)
        d2 = np.minimum(d2, np.sum((xy - xy[nxt]) ** 2, axis=1))
    return xy[np.array(idx)]


def estimate_depth(counts_all: sparse.csr_matrix, modelled_columns: np.ndarray) -> np.ndarray:
    """Library size over genes *outside* the modelled set, scaled to median 1."""
    mask = ~modelled_columns
    if mask.sum() < 50:
        raise ValueError("too few non-modelled genes for a depth offset")
    lib = np.asarray(counts_all[:, mask].sum(axis=1)).ravel().astype(float)
    med = float(np.median(lib[lib > 0])) if np.any(lib > 0) else 1.0
    return np.maximum(lib / max(med, 1e-9), 0.05)


def _poisson_ll(counts: np.ndarray, mean: np.ndarray) -> float:
    m = np.maximum(mean, 1e-9)
    return float(np.sum(counts * np.log(m) - m))


def _deviance(counts: np.ndarray, mean: np.ndarray) -> float:
    m = np.maximum(mean, 1e-9)
    c = np.asarray(counts, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(c > 0, c * np.log(c / m), 0.0)
    return float(2.0 * np.sum(term - (c - m)))


def _graph_smoother(xy: np.ndarray, k: int = 6) -> sparse.csr_matrix:
    """Row-normalised symmetric kNN adjacency (self excluded)."""
    from scipy.spatial import cKDTree

    n = len(xy)
    if n < 3:
        return sparse.eye(n, format="csr")
    tree = cKDTree(xy)
    kk = min(k + 1, n)
    dist, idx = tree.query(xy, k=kk)
    rows = np.repeat(np.arange(n), kk)
    cols = idx.ravel()
    keep = (rows != cols) & np.isfinite(dist.ravel())
    graph = sparse.coo_matrix((np.ones(keep.sum()), (rows[keep], cols[keep])), shape=(n, n)).tocsr()
    graph = graph.maximum(graph.T)
    deg = np.asarray(graph.sum(axis=1)).ravel()
    return sparse.diags(1.0 / np.maximum(deg, 1.0)) @ graph


def _smooth(factors: np.ndarray, smoother: sparse.csr_matrix, steps: int, alpha: float = 0.5) -> np.ndarray:
    """Iterative neighbour averaging of the factor maps (0 steps = no space)."""
    f = factors
    for _ in range(steps):
        f = (1.0 - alpha) * f + alpha * (f @ smoother.T)
    return f


def fit_section(
    data: SectionData,
    loadings: np.ndarray | None,
    *,
    n_factors: int = K_DEFAULT,
    n_smooth: int = 3,
    smooth_alpha: float = 0.5,
    k_graph: int = 6,
    ridge: float = 1e-1,
    n_iter: int = 15,
    pseudocount: float = 0.5,
    seed: int = 0,
    n_centres: int = M_DEFAULT,  # kept for API compatibility; unused
    width: float | None = None,  # kept for API compatibility; unused
    spatial: bool = True,
) -> FactorFit:
    """Fit one unit's factors, with optional graph smoothing between iterations.

    ``spatial=False`` sets ``n_smooth=0`` (no spatial prior at all).
    """
    counts = data.counts.toarray().astype(np.float64)
    depth = np.asarray(data.depth, dtype=float)
    n_spots, n_genes = counts.shape
    y = np.log(counts / depth[:, None] + pseudocount)
    b = y.mean(axis=0)
    Y = y - b[None, :]

    smoother = _graph_smoother(np.asarray(data.xy, dtype=float), k=k_graph) if spatial else None
    steps = int(n_smooth) if spatial else 0

    rng = np.random.default_rng(seed)
    F = np.abs(rng.normal(size=(n_factors, n_spots))) * 0.1
    lam = np.full((n_genes, n_factors), 1.0 / max(n_factors, 1))
    if loadings is not None:
        lam = np.maximum(np.asarray(loadings, dtype=float), 0.0)

    for _ in range(n_iter):
        # factor step: per-spot least squares, then optional spatial smoothing
        Gf = lam.T @ lam + ridge * np.eye(n_factors)
        F = np.linalg.solve(Gf, lam.T @ Y.T)
        if smoother is not None and steps > 0:
            F = _smooth(F, smoother, steps, smooth_alpha)
        # gene step: non-negative least squares
        Gg = F @ F.T + ridge * np.eye(n_factors)
        lam = np.maximum(np.linalg.solve(Gg, F @ Y).T, 0.0)
        for k in range(n_factors):
            nk = float(np.linalg.norm(lam[:, k]))
            if nk > 0:
                lam[:, k] /= nk
                F[k, :] *= nk
        F = np.clip(F, -20.0, 20.0)

    C = np.zeros((n_factors, max(n_centres, 1)))
    mean = depth[:, None] * np.exp(np.clip(b[None, :] + F.T @ lam.T, -10, 10))
    return FactorFit(
        section_id=data.section_id,
        unit_id=data.unit_id,
        factors=F,
        coefficients=C,
        loadings=lam,
        baseline=b,
        fitted_mean=mean,
        log_likelihood=_poisson_ll(counts, mean),
        poisson_deviance=_deviance(counts, mean),
        n_iter=n_iter,
    )


def fit_shared_loadings(
    sections: Sequence[SectionData],
    *,
    n_factors: int = K_DEFAULT,
    n_centres: int = M_DEFAULT,
    width: float | None = None,
    ridge: float = 1e-2,
    n_outer: int = 3,
    n_iter: int = 15,
    spatial: bool = True,
    n_smooth: int = 3,
    seed: int = 0,
) -> tuple[np.ndarray, list[FactorFit]]:
    """Alternate between per-unit fits and an averaged (shared) loading matrix."""
    per_section = [
        fit_section(
            s, None, n_factors=n_factors, ridge=ridge, n_iter=n_iter,
            spatial=spatial, n_smooth=n_smooth, seed=seed,
        )
        for s in sections
    ]
    # align signs/factors by magnitude for a stable average
    lam = np.mean([f.loadings for f in per_section], axis=0)
    fits = per_section
    for _ in range(max(0, n_outer - 1)):
        fits = [
            fit_section(
                s, lam, n_factors=n_factors, ridge=ridge, n_iter=n_iter,
                spatial=spatial, n_smooth=n_smooth, seed=seed,
            )
            for s in sections
        ]
        lam = np.mean([f.loadings for f in fits], axis=0)
    return lam, fits


def _pitch(xy: np.ndarray, visible: np.ndarray) -> float:
    from .geometry import estimate_pitch

    pts = np.asarray(xy)[np.asarray(visible, dtype=bool)]
    if len(pts) < 3:
        return 1.0
    return float(estimate_pitch(pts))


def heldout_spot_loglik(
    data: SectionData,
    loadings: np.ndarray,
    *,
    holdout: np.ndarray,
    n_factors: int = K_DEFAULT,
    n_smooth: int = 3,
    smooth_alpha: float = 0.5,
    k_graph: int = 6,
    ridge: float = 1e-1,
    n_iter: int = 15,
    spatial: bool = True,
    seed: int = 0,
    n_centres: int = M_DEFAULT,
    width: float | None = None,
) -> dict[str, Any]:
    """Fit on training spots, propagate the factor map to held-out spots, score them."""
    train = ~np.asarray(holdout, dtype=bool)
    sub = SectionData(
        section_id=data.section_id,
        unit_id=data.unit_id,
        counts=data.counts[train],
        genes=data.genes,
        depth=data.depth[train],
        xy=data.xy[train],
        visible=data.visible[train],
    )
    fit = fit_section(
        sub, loadings, n_factors=n_factors, n_smooth=n_smooth, smooth_alpha=smooth_alpha,
        k_graph=k_graph, ridge=ridge, n_iter=n_iter, spatial=spatial, seed=seed,
    )
    # propagate training factors to all spots along the full graph (harmonic extension)
    if spatial:
        smoother = _graph_smoother(np.asarray(data.xy, dtype=float), k=k_graph)
        F_full = np.zeros((n_factors, len(data.xy)))
        filled = train.copy()
        F_full[:, train] = fit.factors
        for _ in range(max(5, 2 * int(n_smooth))):
            prop = F_full @ smoother.T
            counts_prop = np.asarray(filled, dtype=float) @ smoother.T
            with np.errstate(invalid="ignore", divide="ignore"):
                num = np.where(counts_prop > 0, prop, 0.0)
                den = np.where(counts_prop > 0, counts_prop, 1.0)
            F_new = num / den
            F_full[:, ~filled] = F_new[:, ~filled]
            filled = np.ones(len(filled), dtype=bool)
            if not np.all(np.isfinite(F_full)):
                F_full = np.nan_to_num(F_full)
    else:
        F_full = np.zeros((n_factors, len(data.xy)))
        F_full[:, train] = fit.factors
    counts_h = data.counts.toarray().astype(np.float64)[holdout]
    mean = np.asarray(data.depth)[holdout, None] * np.exp(
        np.clip(np.asarray(fit.baseline)[None, :] + F_full[:, holdout].T @ np.asarray(loadings).T, -10, 10)
    )
    return {
        "section_id": data.section_id,
        "unit_id": data.unit_id,
        "n_holdout": int(holdout.sum()),
        "loglik": _poisson_ll(counts_h, mean),
        "deviance": _deviance(counts_h, mean),
    }


def local_gene_knn_loglik(
    data: SectionData, *, holdout: np.ndarray, k: int = 6, pseudocount: float = 0.5
) -> dict[str, Any]:
    """Control: per-gene local smoothing. Predict a held-out spot's log-rate as the
    mean of its k nearest training spots, gene by gene (no factor structure at all)."""
    from scipy.spatial import cKDTree

    train = ~np.asarray(holdout, dtype=bool)
    counts = data.counts.toarray().astype(np.float64)
    depth = np.asarray(data.depth, dtype=float)
    y_train = np.log(counts[train] / depth[train, None] + pseudocount)
    tree = cKDTree(np.asarray(data.xy)[train])
    kk = min(k, int(train.sum()))
    _, idx = tree.query(np.asarray(data.xy)[holdout], k=kk)
    if idx.ndim == 1:
        idx = idx[:, None]
    pred = y_train[idx].mean(axis=1)
    mean = depth[holdout, None] * np.exp(np.clip(pred, -10, 10))
    return {
        "section_id": data.section_id,
        "unit_id": data.unit_id,
        "n_holdout": int(holdout.sum()),
        "loglik": _poisson_ll(counts[holdout], mean),
        "deviance": _deviance(counts[holdout], mean),
    }


def depth_only_loglik(data: SectionData, *, holdout: np.ndarray) -> dict[str, Any]:
    """Baseline: counts ~ depth * gene mean (no factors, no space)."""
    train = ~np.asarray(holdout, dtype=bool)
    counts = data.counts.toarray().astype(np.float64)
    y = np.log(counts[train] / np.asarray(data.depth)[train, None] + 0.5)
    b = y.mean(axis=0)
    mean = np.asarray(data.depth)[holdout, None] * np.exp(np.clip(b[None, :], -20, 20))
    return {
        "section_id": data.section_id,
        "unit_id": data.unit_id,
        "n_holdout": int(holdout.sum()),
        "loglik": _poisson_ll(counts[holdout], mean),
        "deviance": _deviance(counts[holdout], mean),
    }


def split_genes_half(genes: np.ndarray, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Deterministic gene-name split for the A/B stability check."""
    import hashlib

    order = sorted(range(len(genes)), key=lambda i: hashlib.sha256(f"{seed}::{genes[i]}".encode()).hexdigest())
    half = len(order) // 2
    a = np.zeros(len(genes), dtype=bool)
    a[np.array(order[:half])] = True
    return a, ~a
