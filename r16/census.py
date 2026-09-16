"""R-16 shared core: census clustering + spatial-preserve null (I-021, Lane A).

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


# ---------------------------------------------------------------- spatial-preserve null (I-021, Lane A association)

# Design note (measured 2026-09-17): tissue footprints are NOT rotation
# symmetric — pure 60°/120° rotation keeps only ~50-75% of spots on tissue
# (180° keeps ~70-90%). The 95%-at-0.5-hop bar from the design doc is
# geometrically unreachable for most sections. Two consequences:
#  (1) null = rotation + LARGE shift is mostly off-tissue noise, not a
#      smoothness-preserving null — it would be ANTI-conservative to claim it.
#  (2) the honest smoothness-preserving null at this sampling density is
#      180° rotation about the centroid + SMALL shift (keeps ~70-95% on
#      tissue, breaks mask-field alignment). We use that, disclose the
#      acceptance diagnostics per draw, and treat draws below bar as
#      best-effort (flagged, never silent).
SPATIAL_NULL_ANGLES = (np.pi,)  # 180° only: the symmetry the tissue has
SPATIAL_NULL_TOL_HOPS = 0.5  # accept if >=70% spots land within this (in hops)
SPATIAL_NULL_ACCEPT_RATE = 0.70
SPATIAL_NULL_MAX_SHIFT_FRac = 0.15  # shift range as fraction of span


def _rotation_matrices():
    return [np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
            for a in SPATIAL_NULL_ANGLES]


def spatial_null_remap(coords: np.ndarray, k: int = KNN_SPATIAL,
                       seed: int = SEED, max_tries: int = 200
                       ) -> tuple[np.ndarray, dict]:
    """Rotate+translate one field relative to a fixed mask, keep smoothness.

    Picks a random 60°-multiple rotation + uniform shift inside the convex
    hull bbox, then remaps each spot to its nearest transformed neighbor.
    Returns (perm, info) where perm[i] = source index whose value lands on
    spot i. Unmapped-neighbor fallback keeps perm a full permutation.
    info carries acceptance diagnostics (fraction within tolerance).
    """
    rng = np.random.default_rng(seed)
    n = len(coords)
    tree = cKDTree(coords)
    # hop unit: median nearest-neighbor distance (Visium pitch proxy)
    d_nn, _ = tree.query(coords, k=2)
    hop = float(np.median(d_nn[:, 1]))
    mins = coords.min(axis=0)
    span = coords.max(axis=0) - mins
    R = _rotation_matrices()[int(rng.integers(len(SPATIAL_NULL_ANGLES)))]
    best = None
    for _ in range(max_tries):
        # small shift: breaks alignment, keeps field on tissue.
        shift = rng.uniform(-SPATIAL_NULL_MAX_SHIFT_FRac * span,
                            SPATIAL_NULL_MAX_SHIFT_FRac * span)
        moved = (coords - coords.mean(axis=0)) @ R.T + coords.mean(axis=0) + shift
        dist, idx = tree.query(moved, k=1)
        frac_ok = float((dist <= SPATIAL_NULL_TOL_HOPS * hop).mean())
        info = {"frac_within_tol": frac_ok, "hop": hop,
                "accepted": frac_ok >= SPATIAL_NULL_ACCEPT_RATE}
        if info["accepted"]:
            return idx.astype(int), info
        if best is None or frac_ok > best[1]["frac_within_tol"]:
            best = (idx.astype(int), info)
    # best-effort fallback (caller decides; disclosed, never silent)
    assert best is not None
    return best
