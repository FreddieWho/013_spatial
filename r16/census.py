"""R-16 composition pattern census core (design doc v1.1 §3).

Per-section Leiden clustering on HVG-10000 expression, cluster signatures,
cross-section signature matching, resolution stability, Moran's I coherence
with value-permutation null (valid for autocorrelation tests, I-021 refined).

Deterministic: all stochastic steps take explicit seeds.
"""

from __future__ import annotations

import json
from pathlib import Path

import igraph as ig
import leidenalg
import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

SEED = 20260914
RESOLUTIONS = (0.25, 0.5, 1.0)
MIN_CLUSTER_SIZE = 20
PCA_DIMS = 30
KNN_EXPR = 15
KNN_SPATIAL = 6
MATCH_COSINE = 0.75
STABILITY_JACCARD = 0.5
NULL_DRAWS = 200


# ---------------------------------------------------------------- cache I/O


def load_section(cache_dir: Path, stem: str):
    """Return (counts csr, barcodes, genes, coords) from a census cache entry."""
    js = json.loads((cache_dir / f"{stem}.json").read_text())
    z = np.load(cache_dir / f"{stem}.npz")
    mat = sparse.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    return mat, js["barcodes"], js["genes"], z["coords"].astype(float), js


def list_stems(cache_dir: Path) -> list[str]:
    return sorted(p.name[:-5] for p in cache_dir.glob("*.json") if p.name != "cache_manifest.json")


# ---------------------------------------------------------------- clustering


def log1p_zscore(mat: sparse.csr_matrix) -> np.ndarray:
    x = mat.toarray().astype(np.float32)
    np.log1p(x, out=x)
    mu = x.mean(axis=0)
    sd = x.std(axis=0)
    sd[sd < 1e-8] = 1.0
    return (x - mu) / sd


def pca_scores(x: np.ndarray, n_comps: int = PCA_DIMS, seed: int = SEED) -> np.ndarray:
    n_comps = min(n_comps, min(x.shape) - 1)
    return PCA(n_components=n_comps, random_state=seed).fit_transform(x)


def knn_igraph(scores: np.ndarray, k: int) -> ig.Graph:
    nn = NearestNeighbors(n_neighbors=k + 1).fit(scores)
    _, idx = nn.kneighbors(scores)
    edges = [(i, j) for i in range(len(scores)) for j in idx[i, 1:]]
    g = ig.Graph(n=len(scores), edges=edges)
    return g.simplify()


def leiden_labels(graph: ig.Graph, resolution: float, seed: int = SEED) -> np.ndarray:
    part = leidenalg.find_partition(
        graph, leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=resolution, seed=seed, n_iterations=-1,
    )
    return np.asarray(part.membership, dtype=int)


def filter_min_size(labels: np.ndarray, min_size: int = MIN_CLUSTER_SIZE) -> np.ndarray:
    out = labels.copy()
    for c in np.unique(labels):
        if (labels == c).sum() < min_size:
            out[labels == c] = -1
    return out


# ---------------------------------------------------------------- signatures


def cluster_signatures(x_z: np.ndarray, labels: np.ndarray) -> dict[int, np.ndarray]:
    return {c: x_z[labels == c].mean(axis=0) for c in np.unique(labels) if c >= 0}


def top_markers(x_log: np.ndarray, labels: np.ndarray, cluster: int, genes: list[str],
                topn: int = 10) -> list[str]:
    inn = x_log[labels == cluster]
    out = x_log[labels != cluster]
    if len(inn) == 0 or len(out) == 0:
        return []
    diff = inn.mean(axis=0) - out.mean(axis=0)
    se = np.sqrt(inn.var(axis=0) / len(inn) + out.var(axis=0) / len(out)) + 1e-8
    t = diff / se
    order = np.argsort(-t)[:topn]
    return [genes[i] for i in order if t[i] > 0]


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    return float(a @ b / (na * nb))


def jaccard(a: np.ndarray, b: np.ndarray) -> float:
    u = np.logical_or(a, b).sum()
    return float(np.logical_and(a, b).sum() / u) if u else 0.0


def match_units_to_groups(signatures: list) -> tuple[list, dict]:
    """Average-linkage hierarchical matching of per-(section, cluster) units.

    Used when cluster ids are section-local (e.g. per-section GraphST
    training): raw ids must NEVER be pooled across sections. Mirrors the
    Tier-2 v2.0 convention: cosine distance, average linkage, cut at
    1 - MATCH_COSINE; purity = min pairwise cosine within group (singleton
    groups get 1.0).
    Returns (group_ids aligned to input order, {group_id: purity}).
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    n = len(signatures)
    if n == 0:
        return [], {}
    if n == 1:
        return [0], {0: 1.0}
    mat = np.vstack([np.asarray(s, dtype=float) for s in signatures])
    nrm = np.linalg.norm(mat, axis=1, keepdims=True)
    nrm[nrm == 0] = 1.0
    sim = (mat / nrm) @ (mat / nrm).T
    np.fill_diagonal(sim, 1.0)
    dist = 1 - np.clip(sim, -1, 1)
    Z = linkage(squareform(dist, checks=False), method="average")
    lab = fcluster(Z, t=1 - MATCH_COSINE, criterion="distance")
    groups: dict = {}
    for i, g in enumerate(lab):
        groups.setdefault(int(g), []).append(i)
    gids = [0] * n
    purities = {}
    for gi, (g, members) in enumerate(sorted(groups.items())):
        for i in members:
            gids[i] = gi
        if len(members) > 1:
            sub = sim[np.ix_(members, members)].copy()
            np.fill_diagonal(sub, 1.0)
            purities[gi] = float(sub.min())
        else:
            purities[gi] = 1.0
    return gids, purities


# ---------------------------------------------------------------- Moran's I


def spatial_weights(coords: np.ndarray, k: int = KNN_SPATIAL) -> sparse.csr_matrix:
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, idx = nn.kneighbors(coords)
    rows = np.repeat(np.arange(len(coords)), k)
    cols = idx[:, 1:].ravel()
    w = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(coords),) * 2)
    return w.maximum(w.T)


def moran_i(values: np.ndarray, w: sparse.csr_matrix) -> float:
    z = values - values.mean()
    denom = float(z @ z)
    if denom == 0:
        return 0.0
    s0 = w.sum()
    return float(len(values) / s0 * (z @ (w @ z)) / denom)


def moran_permutation_p(values: np.ndarray, w: sparse.csr_matrix, draws: int,
                        rng: np.random.Generator) -> tuple[float, float]:
    """Value-permutation null: valid for autocorrelation tests (I-021 refined)."""
    obs = moran_i(values, w)
    null = np.empty(draws)
    for i in range(draws):
        null[i] = moran_i(rng.permutation(values), w)
    p = float((np.sum(null >= obs) + 1) / (draws + 1))
    return obs, p
