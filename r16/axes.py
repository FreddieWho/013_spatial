"""R-16 Tier-1 composition axes (design v2.0, D-115).

Continuous marker-module score rulers: axis score per spot = mean z-scored
log1p over the axis's marker genes (z within section). Cross-patient
comparability is constructive: same genes, same formula — no matching step.
Discrete geometry is derived by contouring at DISCLOSED quantiles (0.6/0.7/0.8
swept, stability reported). Validated end-to-end by the CXCL13→B-contour
positive control (25/30 patients peak inside/edge).
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path

from r16.census import jaccard

QUANTILES = (0.6, 0.7, 0.8)
PLASMA_GENES = ["IGHM", "IGKC", "JCHAIN", "MZB1", "SDC1", "IGLC3", "IGHA2"]
# ILC retired (D-117): rare, poorly defined on Visium spots; empirically weakest axis.
RETIRED_AXES = frozenset({"ILC"})
WORKING_AXES = ("T", "B", "Mye", "Epi", "Stromal", "Plasma")


def axis_scores(x_z: np.ndarray, gidx: dict[str, int], axis_genes: list[str]):
    """Return (score per spot, n genes actually used)."""
    cols = [gidx[g] for g in axis_genes if g in gidx]
    if not cols:
        return np.zeros(x_z.shape[0], dtype=np.float32), 0
    return x_z[:, cols].mean(axis=1), len(cols)


def contour_mask(scores: np.ndarray, q: float) -> np.ndarray:
    return scores >= np.quantile(scores, q)


def contour_stability(scores: np.ndarray, quantiles=QUANTILES):
    """Mean and per-step Jaccard between adjacent-quantile masks."""
    masks = [contour_mask(scores, q) for q in quantiles]
    js = [jaccard(masks[i], masks[i + 1]) for i in range(len(masks) - 1)]
    return float(np.mean(js)), js


def signed_hop_distance(w: sparse.csr_matrix, mask: np.ndarray) -> np.ndarray:
    """Signed BFS hops to mask boundary: outside positive, inside negative."""
    W = w.astype(bool).astype(np.float32)
    n = W.shape[0]
    idx_in = np.flatnonzero(mask)
    idx_out = np.flatnonzero(~mask)
    if len(idx_in) == 0 or len(idx_out) == 0:
        return np.zeros(n)
    d_out = shortest_path(W, directed=False, indices=idx_in).min(axis=0)
    d_in = shortest_path(W, directed=False, indices=idx_out).min(axis=0)
    return np.where(mask, -d_in, d_out)


def binned_profile(values: np.ndarray, signed: np.ndarray, bins: np.ndarray):
    """Per-hop (value mean, n spots); empty bins -> (nan, 0)."""
    out = []
    for b in bins:
        m = signed == b
        out.append((float(values[m].mean()), int(m.sum())) if m.any() else (np.nan, 0))
    return out
