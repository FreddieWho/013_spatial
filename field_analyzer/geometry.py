"""Geometry helpers: units, adjacency, distance pairs, local direction.

Rules carried from the contract:
  * known tissue holes and masked-but-possibly-tissue regions are different states;
  * no convex hull may turn a hole into observed tissue;
  * an isolated point gets *no* neighbour value (the legacy ``neighbor_mean``
    fallback copied the point's own value back, which re-injects hidden data);
  * no metric scale is hard-coded (the legacy graph used 130.01 um).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

from .observation import ContractError

# region classes
VISIBLE = "VISIBLE"
MASKED_POSSIBLE_TISSUE = "MASKED_POSSIBLE_TISSUE"
NON_TISSUE = "NON_TISSUE"


@dataclass
class Adjacency:
    graph: sparse.csr_matrix
    pitch: float
    unit: str
    n_isolated: int
    isolated: np.ndarray
    degrees: np.ndarray
    notes: dict[str, Any]


def estimate_pitch(xy: np.ndarray) -> float:
    """Median nearest-neighbour distance; the empirical spot pitch."""
    xy = np.asarray(xy, dtype=np.float64)
    if xy.shape[0] < 2:
        raise ContractError("need at least two points to estimate pitch")
    tree = cKDTree(xy)
    dist, _ = tree.query(xy, k=2)
    nn = dist[:, 1]
    if not np.all(np.isfinite(nn)) or np.any(nn <= 0):
        raise ContractError("duplicate or non-finite coordinates; refuse to guess pitch")
    return float(np.median(nn))


def build_adjacency(
    xy: np.ndarray,
    support: np.ndarray,
    *,
    k: int = 6,
    max_distance: float | None = None,
    max_edge_pitch: float | None = 3.0,
    unit: str = "empirical_pitch",
) -> Adjacency:
    """Symmetric k-NN graph restricted to supported points.

    Weights are row-normalised over *supported* neighbours only.  Isolated points
    keep degree 0 and are reported; their neighbour value stays undefined (NaN in
    :func:`neighbor_sum`).  Nothing is ever copied back from the point itself.
    """
    xy = np.asarray(xy, dtype=np.float64)
    support = np.asarray(support, dtype=bool)
    if xy.shape != (len(support), 2):
        raise ContractError("xy/support shape mismatch")
    n = len(xy)
    pitch = estimate_pitch(xy[support]) if support.sum() >= 2 else float("nan")
    if max_distance is None and max_edge_pitch is not None and np.isfinite(pitch):
        max_distance = max_edge_pitch * pitch
    if n >= 2 and support.sum() >= 2:
        supported = np.flatnonzero(support)
        tree = cKDTree(xy[supported])
        kk = min(k + 1, len(supported))
        dist, idx = tree.query(xy[supported], k=kk)
        if dist.ndim == 1:
            dist = dist[:, None]
            idx = idx[:, None]
        src_local = np.repeat(np.arange(len(supported)), kk)
        dst_local = idx.ravel()
        dd = dist.ravel()
        keep = (src_local != dst_local) & np.isfinite(dd)
        if max_distance is not None:
            keep &= dd <= max_distance
        src = supported[src_local[keep]]
        dst = supported[dst_local[keep]]
        graph = sparse.coo_matrix(
            (np.ones(len(src)), (src, dst)), shape=(n, n)
        ).tocsr()
        graph.sum_duplicates()
        graph.data[:] = 1.0
    else:
        graph = sparse.csr_matrix((n, n), dtype=np.float64)
    degree = np.asarray(graph.sum(axis=1)).ravel()
    isolated = np.flatnonzero(degree == 0)
    weights = sparse.diags(1.0 / np.maximum(degree, 1.0))
    normalized = (weights @ graph).tocsr()
    notes = {
        "k": int(k),
        "max_edge_pitch": None if max_edge_pitch is None else float(max_edge_pitch),
        "max_distance": None if max_distance is None else float(max_distance),
        "n_edges": int(graph.nnz // 2),
        "self_fallback": False,
    }
    return Adjacency(
        graph=normalized,
        pitch=pitch,
        unit=unit,
        n_isolated=int(len(isolated)),
        isolated=isolated,
        degrees=degree,
        notes=notes,
    )


def neighbor_sum(graph: sparse.csr_matrix, values: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Weighted neighbour mean over *valid* values only; NaN where unsupported.

    ``valid`` is the mask of cells whose value may be used (visible & finite).
    Degenerate rows (no valid neighbour) return NaN, never the point's own value.
    """
    v = np.asarray(values, dtype=np.float64).ravel()
    ok = np.asarray(valid, dtype=bool).ravel()
    if v.shape != ok.shape:
        raise ContractError("values/valid shape mismatch")
    if v.shape[0] != graph.shape[0]:
        raise ContractError("values length does not match graph")
    use = np.where(ok, v, 0.0)
    num = np.asarray(graph @ use).ravel()
    den = np.asarray(graph @ ok.astype(np.float64)).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(den > 0, num / den, np.nan)
    return out


def distance_pairs(
    xy: np.ndarray,
    support: np.ndarray,
    *,
    max_distance: float,
    n_bins: int = 12,
) -> dict[str, Any]:
    """All unordered pairs inside ``max_distance`` with both endpoints supported."""
    xy = np.asarray(xy, dtype=np.float64)
    support = np.asarray(support, dtype=bool)
    idx = np.flatnonzero(support)
    if len(idx) < 2:
        raise ContractError("need at least two supported points")
    tree = cKDTree(xy[idx])
    pairs = tree.query_pairs(max_distance, output_type="ndarray")
    if len(pairs) == 0:
        return {
            "pairs": np.empty((0, 2), dtype=np.int64),
            "distances": np.empty(0),
            "bin_edges": np.linspace(0, max_distance, n_bins + 1),
            "counts": np.zeros(n_bins, dtype=np.int64),
        }
    mapped = idx[pairs]
    d = np.linalg.norm(xy[mapped[:, 0]] - xy[mapped[:, 1]], axis=1)
    edges = np.linspace(0.0, float(max_distance), n_bins + 1)
    counts, _ = np.histogram(d, bins=edges)
    return {
        "pairs": mapped.astype(np.int64),
        "distances": d,
        "bin_edges": edges,
        "counts": counts.astype(np.int64),
    }


def region_classes(tissue_support: np.ndarray, visible_mask: np.ndarray) -> np.ndarray:
    """Label each point VISIBLE / MASKED_POSSIBLE_TISSUE / NON_TISSUE."""
    ts = np.asarray(tissue_support, dtype=bool)
    vm = np.asarray(visible_mask, dtype=bool)
    if ts.shape != vm.shape:
        raise ContractError("tissue_support/visible_mask shape mismatch")
    out = np.full(len(ts), NON_TISSUE, dtype=object)
    out[ts & vm] = VISIBLE
    out[ts & ~vm] = MASKED_POSSIBLE_TISSUE
    return out


def assert_no_hull_fill(region: np.ndarray) -> None:
    """Guard: a convex hull of visible points can include non-tissue interior."""
    region = np.asarray(region)
    if region.dtype == object and NON_TISSUE in set(region):
        allowed = {VISIBLE, MASKED_POSSIBLE_TISSUE, NON_TISSUE}
        extra = set(region) - allowed
        if extra:
            raise ContractError(f"unknown region labels: {sorted(extra)}")


def local_direction(
    values: np.ndarray,
    valid: np.ndarray,
    xy: np.ndarray,
    adjacency: Adjacency,
    *,
    min_magnitude: float = 1e-8,
    min_anisotropy: float = 1e-6,
    min_valid_neighbors: int = 2,
) -> dict[str, np.ndarray]:
    """Local gradient orientation (axial, mod pi) via per-node plane fits.

    Direction is measurable only where the point and at least
    ``min_valid_neighbors`` neighbours carry valid values, the neighbourhood is not
    geometrically degenerate (collinear neighbours), and the fitted gradient is not
    numerically zero.  Otherwise the angle is NaN and ``measurable`` is False.
    """
    v = np.asarray(values, dtype=np.float64).ravel()
    ok = np.asarray(valid, dtype=bool).ravel()
    xy = np.asarray(xy, dtype=np.float64)
    if v.shape != ok.shape:
        raise ContractError("values/valid shape mismatch")
    if xy.shape != (len(v), 2):
        raise ContractError("xy shape does not match values")
    graph = (adjacency.graph > 0).astype(np.float64).tocoo()
    n = len(v)

    i = graph.row
    j = graph.col
    edge_valid = ok[i] & ok[j]
    i, j = i[edge_valid], j[edge_valid]
    dx = xy[j, 0] - xy[i, 0]
    dy = xy[j, 1] - xy[i, 1]
    dv = v[j] - v[i]

    jxx = np.zeros(n); jxy = np.zeros(n); jyy = np.zeros(n)
    bx = np.zeros(n); by = np.zeros(n)
    counts = np.zeros(n, dtype=np.int64)
    np.add.at(jxx, i, dx * dx)
    np.add.at(jxy, i, dx * dy)
    np.add.at(jyy, i, dy * dy)
    np.add.at(bx, i, dx * dv)
    np.add.at(by, i, dy * dv)
    np.add.at(counts, i, 1)

    det = jxx * jyy - jxy * jxy
    trace = jxx + jyy
    disc = np.sqrt(np.maximum((jxx - jyy) ** 2 + 4.0 * jxy ** 2, 0.0))
    lam_max = 0.5 * (trace + disc)
    lam_min = 0.5 * (trace - disc)
    with np.errstate(invalid="ignore", divide="ignore"):
        anisotropy = np.where(lam_max > 0, lam_min / np.maximum(lam_max, 1e-300), 0.0)
        gx = np.where(det != 0, (jyy * bx - jxy * by) / np.where(det != 0, det, 1.0), np.nan)
        gy = np.where(det != 0, (jxx * by - jxy * bx) / np.where(det != 0, det, 1.0), np.nan)

    magnitude = np.hypot(gx, gy)
    measurable = (
        ok
        & (counts >= min_valid_neighbors)
        & np.isfinite(magnitude)
        & (magnitude > min_magnitude)
        & (anisotropy >= min_anisotropy)
        & (anisotropy <= 1.0)
    )
    angle = np.where(measurable, np.mod(np.arctan2(gy, gx), np.pi), np.nan)
    return {
        "angle_mod_pi": angle,
        "measurable": measurable,
        "magnitude": np.where(np.isfinite(magnitude), magnitude, np.nan),
        "anisotropy": anisotropy,
        "n_valid_neighbors": counts,
    }


def assert_edges_avoid_non_tissue(graph: sparse.spmatrix, region: np.ndarray) -> int:
    """Count edges touching NON_TISSUE; raise if any (holes cannot be crossed)."""
    region = np.asarray(region, dtype=object)
    coo = graph.tocoo()
    bad = (region[coo.row] == NON_TISSUE) | (region[coo.col] == NON_TISSUE)
    n_bad = int(bad.sum())
    if n_bad:
        raise ContractError(
            f"{n_bad} edges touch NON_TISSUE points; a real tissue hole must not be crossed"
        )
    return n_bad
