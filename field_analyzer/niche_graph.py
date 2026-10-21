"""FA-11 (route C): niche transition matrix on single-cell data.

Continuous "fields" turned out to be too smooth and too noisy to support structural
claims.  This module replaces them with a *discrete* object: niches, and the
transition matrix between adjacent niches.

    labels  : K niches from a fixed marker-program score matrix (k-means, seeded)
    graph   : symmetric kNN graph on cell centroids
    T[i,j]  : P(neighbour is j | cell is i)      (row-normalised adjacency counts)
    asymmetry statistic S = ||T - T^T||_F over off-diagonal entries

"Connection" becomes a testable statement: **is the transition matrix asymmetric
beyond what the same labels would give under a null that keeps their spatial
autocorrelation but destroys the specific adjacency?**

Two cheap nulls are used:
  * rotation  : rotate the label field about the centroid (keeps autocorrelation,
                breaks which cells are adjacent)
  * permutation: shuffle labels among cells (keeps composition, destroys space)
"""

from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import h5py
import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

PROGRAMS: dict[str, list[str]] = {
    # the colon panel has no stroma/vessel canonical markers; these are what exists
    "T_NK": ["CD3D", "CD3E", "CD3G", "CD2", "CD5", "CD6", "CD7", "CD8A", "CD8B"],
    "B_plasma": ["MS4A1", "CD79A", "CD79B", "BANK1"],
    "myeloid": ["CD14", "CD163", "CD177", "MS4A7"],
    "epithelial": ["CEACAM5", "CEACAM6", "CEACAM7", "MUC12", "TFF1", "OLFM4", "AQP8",
                   "ANXA13", "CDHR5", "MS4A12", "EPHB3", "SOX9", "ASCL2"],
    "proliferative": ["MKI67"],
    "Treg": ["FOXP3"],
}


@dataclass
class XeniumSlice:
    counts: sparse.csr_matrix
    genes: np.ndarray
    xy: np.ndarray
    cell_ids: np.ndarray

    @property
    def n_cells(self) -> int:
        return self.counts.shape[0]


def read_xenium(directory) -> XeniumSlice:
    """Read the 10x Xenium ``cell_feature_matrix.h5`` + ``cells.csv.gz`` pair."""
    from pathlib import Path

    directory = Path(directory)
    with h5py.File(directory / "cell_feature_matrix.h5", "r") as f:
        m = f["matrix"]
        shape = tuple(int(v) for v in m["shape"][:])  # features x cells
        data = m["data"][:]
        indices = m["indices"][:]
        indptr = m["indptr"][:]
        genes = np.array([x.decode() for x in m["features"]["name"][:]])
        barcodes = np.array([x.decode() for x in m["barcodes"][:]])
        mat = sparse.csc_matrix((data, indices, indptr), shape=shape)
    counts = mat.T.tocsr()  # cells x genes

    xs, ys, ids = [], [], []
    with gzip.open(directory / "cells.csv.gz", "rt") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            ids.append(row["cell_id"])
            xs.append(float(row["x_centroid"]))
            ys.append(float(row["y_centroid"]))
    xy = np.column_stack([xs, ys])
    id_index = {cid: i for i, cid in enumerate(ids)}
    order = np.array([id_index[b] for b in barcodes], dtype=np.int64)
    return XeniumSlice(counts=counts, genes=genes, xy=xy[order], cell_ids=barcodes)


def program_scores(slice_: XeniumSlice, programs: Mapping[str, Sequence[str]] = PROGRAMS) -> tuple[np.ndarray, list[str]]:
    """Per-cell z-scored program scores from log-normalised counts."""
    genes = np.asarray(slice_.genes).astype(str)
    lib = np.maximum(np.asarray(slice_.counts.sum(axis=1)).ravel(), 1.0)
    names, cols = [], []
    for name, members in programs.items():
        idx = np.flatnonzero(np.isin(genes, list(members)))
        if len(idx) == 0:
            continue
        sub = slice_.counts[:, idx].toarray().astype(np.float64)
        norm = np.log1p(sub * 1e4 / lib[:, None])
        score = norm.mean(axis=1)
        sd = float(np.std(score))
        names.append(name)
        cols.append((score - float(np.mean(score))) / (sd if sd > 0 else 1.0))
    return np.column_stack(cols) if cols else np.zeros((slice_.n_cells, 0)), names


def kmeans_niches(scores: np.ndarray, k: int = 5, *, seed: int = 0, iters: int = 30) -> np.ndarray:
    """Plain k-means (sklearn is unusable under numpy 2 here); deterministic seed."""
    x = np.asarray(scores, dtype=np.float64)
    rng = np.random.default_rng(seed)
    n = len(x)
    if n == 0 or x.shape[1] == 0:
        return np.zeros(n, dtype=int)
    centres = x[rng.choice(n, size=k, replace=False)].copy()
    labels = np.zeros(n, dtype=int)
    for _ in range(iters):
        d2 = ((x[:, None, :] - centres[None, :, :]) ** 2).sum(axis=2)
        new = d2.argmin(axis=1)
        if np.array_equal(new, labels):
            break
        labels = new
        for j in range(k):
            sel = labels == j
            if sel.any():
                centres[j] = x[sel].mean(axis=0)
    return labels


def knn_adjacency(xy: np.ndarray, k: int = 8) -> sparse.csr_matrix:
    xy = np.asarray(xy, dtype=np.float64)
    n = len(xy)
    tree = cKDTree(xy)
    kk = min(k + 1, n)
    dist, idx = tree.query(xy, k=kk)
    rows = np.repeat(np.arange(n), kk)
    cols = idx.ravel()
    keep = (rows != cols) & np.isfinite(dist.ravel())
    graph = sparse.coo_matrix((np.ones(int(keep.sum())), (rows[keep], cols[keep])), shape=(n, n)).tocsr()
    graph = graph.maximum(graph.T)
    graph.setdiag(0)
    graph.eliminate_zeros()
    return graph


def transition_matrix(labels: np.ndarray, graph: sparse.csr_matrix, k: int) -> np.ndarray:
    """Row-normalised niche transition counts."""
    labels = np.asarray(labels, dtype=int)
    coo = graph.tocoo()
    counts = np.zeros((k, k), dtype=np.float64)
    np.add.at(counts, (labels[coo.row], labels[coo.col]), 1.0)
    row = counts.sum(axis=1, keepdims=True)
    return np.where(row > 0, counts / np.maximum(row, 1e-12), 0.0)


def contact_counts(labels: np.ndarray, graph: sparse.csr_matrix, k: int) -> tuple[np.ndarray, np.ndarray, float]:
    """Undirected niche contact counts, per-niche degree totals and the edge total."""
    labels = np.asarray(labels, dtype=int)
    coo = graph.tocoo()
    C = np.zeros((k, k), dtype=np.float64)
    np.add.at(C, (labels[coo.row], labels[coo.col]), 1.0)
    row = C.sum(axis=1)
    total = float(graph.nnz)  # one entry per direction on a symmetric graph
    return C, row, total


def contact_enrichment(labels: np.ndarray, graph: sparse.csr_matrix, k: int) -> dict[str, Any]:
    """Observed / expected contacts under label independence, given the graph degrees.

    .. warning::
       A **symmetric** graph makes the directed count matrix exactly symmetric
       (``max|C - C^T| == 0``), so any "i -> j asymmetry" computed on it is only the
       row-normalisation artefact of unequal niche sizes.  This function therefore
       reports *co-localisation enrichment* (a symmetric quantity), not direction.
       Direction needs an external axis (e.g. distance to a tumour boundary).
    """
    C, row, total = contact_counts(labels, graph, k)
    expected = np.outer(row, row) / max(total, 1e-12)
    with np.errstate(divide="ignore", invalid="ignore"):
        E = np.where(expected > 0, C / expected, np.nan)
    off = ~np.eye(k, dtype=bool)
    logE = np.log(np.clip(E, 1e-6, None))
    return {
        "statistic": float(np.linalg.norm(logE[off])),
        "max_enrichment": float(np.nanmax(E[off])),
        "min_enrichment": float(np.nanmin(E[off])),
        "enrichment": E,
        "counts": C,
        "counts_symmetry_error": float(np.max(np.abs(C - C.T))),
    }


def asymmetry(T: np.ndarray) -> dict[str, Any]:
    """Kept for backwards compatibility; see the warning in :func:`contact_enrichment`."""
    A = T - T.T
    return {"frobenius": float(np.linalg.norm(A)), "matrix": A}


def rotation_null_labels(xy: np.ndarray, labels: np.ndarray, angle_deg: float, *, seed: int = 0) -> np.ndarray:
    """Rotate the label field about the centroid, then re-assign by nearest cell."""
    xy = np.asarray(xy, dtype=np.float64)
    centre = xy.mean(axis=0)
    ang = np.deg2rad(angle_deg)
    rot = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    rotated = (xy - centre) @ rot.T + centre
    tree = cKDTree(xy)
    _, idx = tree.query(rotated, k=1)
    return np.asarray(labels)[idx]


def permutation_null_labels(labels: np.ndarray, *, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return np.asarray(labels)[rng.permutation(len(labels))]


def null_distribution(
    xy: np.ndarray,
    labels: np.ndarray,
    graph: sparse.csr_matrix,
    k: int,
    *,
    n_rot: int = 60,
    n_perm: int = 60,
    seed: int = 0,
) -> dict[str, Any]:
    obs = contact_enrichment(labels, graph, k)
    rot_stats, perm_stats = [], []
    for i in range(n_rot):
        lab = rotation_null_labels(xy, labels, angle_deg=360.0 * (i + 1) / (n_rot + 1), seed=seed + i)
        rot_stats.append(contact_enrichment(lab, graph, k)["statistic"])
    for i in range(n_perm):
        lab = permutation_null_labels(labels, seed=seed + 1000 + i)
        perm_stats.append(contact_enrichment(lab, graph, k)["statistic"])
    rot_stats = np.asarray(rot_stats)
    perm_stats = np.asarray(perm_stats)
    return {
        "observed_statistic": obs["statistic"],
        "observed_max_enrichment": obs["max_enrichment"],
        "observed_min_enrichment": obs["min_enrichment"],
        "counts_symmetry_error": obs["counts_symmetry_error"],
        "enrichment_matrix": obs["enrichment"],
        "rotation_null_mean": float(rot_stats.mean()),
        "rotation_null_p95": float(np.quantile(rot_stats, 0.95)),
        "rotation_z": float((obs["statistic"] - rot_stats.mean()) / max(rot_stats.std(), 1e-12)),
        "permutation_null_mean": float(perm_stats.mean()),
        "permutation_null_p95": float(np.quantile(perm_stats, 0.95)),
        "permutation_z": float((obs["statistic"] - perm_stats.mean()) / max(perm_stats.std(), 1e-12)),
        "n_rotation": int(n_rot),
        "n_permutation": int(n_perm),

    }


def half_consistency(xy: np.ndarray, labels: np.ndarray, graph: sparse.csr_matrix, k: int) -> dict[str, Any]:
    """Split by the median x coordinate; do the two halves agree on the enrichment pattern?"""
    x = np.asarray(xy)[:, 0]
    left = x <= np.median(x)
    out = {}
    for name, sel in (("left", left), ("right", ~left)):
        sub_graph = graph[sel][:, sel]
        out[name] = contact_enrichment(np.asarray(labels)[sel], sub_graph, k)["enrichment"]
    off = ~np.eye(k, dtype=bool)
    a = np.log(np.clip(out["left"][off], 1e-6, None))
    b = np.log(np.clip(out["right"][off], 1e-6, None))
    return {
        "corr_of_log_enrichment": float(np.corrcoef(a, b)[0, 1]) if np.std(a) > 0 and np.std(b) > 0 else float("nan"),
        "sign_agreement": float(np.mean(np.sign(a) == np.sign(b))),
        "n_pairs": int(off.sum()),
    }
