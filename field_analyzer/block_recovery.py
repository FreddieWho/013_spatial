"""FA-05: empirical recoverable range under a blocked region (BLOCK_RECOVER).

Task semantics (05_COMMON_PROTOCOL §1): a continuous disc H is unobserved; the model
may use everything in O = Ω \ H.  Query coordinates are known.  We score on a *fixed*
core Q (r <= 2 pitch) so that changing the mask radius does not change the target set.

Only classical references are implemented here: they bound what "recovery" means
before any learned relation model is claimed to add something.

  background     : visible mean of the target field (no spatial information)
  IDW(k=8)       : inverse-distance weighted nearest visible spots of the target field
  kernel         : gaussian kernel smoothing of the target field, bandwidth = 2 pitch
  multifield     : ridge on kernel-smoothed *other* fields (same visible budget)

A learned cross-field operator is NOT included; that is the FA-03B/05 follow-up.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.spatial import cKDTree

from . import geometry


@dataclass
class Predictions:
    name: str
    values: np.ndarray  # full length, NaN where not computable


def predict_background(visible_values: np.ndarray, n: int) -> Predictions:
    v = np.asarray(visible_values, dtype=float)
    finite = np.isfinite(v)
    mean = float(np.nanmean(v[finite])) if finite.any() else float("nan")
    return Predictions("background_mean", np.full(n, mean))


def predict_idw(
    xy_visible: np.ndarray, visible_values: np.ndarray, query_xy: np.ndarray, *, k: int = 8
) -> Predictions:
    v = np.asarray(visible_values, dtype=float)
    ok = np.isfinite(v)
    if ok.sum() < k:
        return Predictions("idw_k8", np.full(len(query_xy), np.nan))
    tree = cKDTree(xy_visible[ok])
    kk = min(k, int(ok.sum()))
    dist, idx = tree.query(query_xy, k=kk)
    if dist.ndim == 1:
        dist = dist[:, None]
        idx = idx[:, None]
    weights = 1.0 / np.maximum(dist, 1e-6)
    vals = v[ok][idx]
    out = np.sum(weights * vals, axis=1) / np.sum(weights, axis=1)
    return Predictions(f"idw_k{k}", out)


def predict_kernel(
    xy_visible: np.ndarray,
    visible_values: np.ndarray,
    query_xy: np.ndarray,
    *,
    k: int = 16,
    bandwidth: float | None = None,
) -> Predictions:
    """Gaussian kernel smoothing with an *adaptive* bandwidth.

    With a fixed bandwidth the estimate is undefined once the hole is wider than the
    kernel reach (that is a real "cannot estimate", not an error), but then the
    comparison across mask radii is impossible.  The default therefore takes the
    k-th neighbour distance as the bandwidth; pass ``bandwidth`` for the fixed
    version, and read ``n_core`` to see whether anything was computable.
    """
    v = np.asarray(visible_values, dtype=float)
    ok = np.isfinite(v)
    if ok.sum() < k:
        return Predictions("kernel", np.full(len(query_xy), np.nan))
    tree = cKDTree(xy_visible[ok])
    kk = min(k, int(ok.sum()))
    dist, idx = tree.query(query_xy, k=kk)
    if dist.ndim == 1:
        dist = dist[:, None]
        idx = idx[:, None]
    vals = v[ok][idx]
    bw = np.full(len(query_xy), float(bandwidth)) if bandwidth else np.maximum(dist[:, -1], 1e-6)
    w = np.exp(-0.5 * (dist / bw[:, None]) ** 2)
    out = np.sum(w * vals, axis=1) / np.sum(w, axis=1)
    return Predictions(f"kernel_k{kk}", out)


def predict_multifield(
    query_xy: np.ndarray,
    feature_predictions: Mapping[str, np.ndarray],
    target_visible_values: np.ndarray,
    visible_mask: np.ndarray,
    *,
    alpha: float = 1.0,
) -> Predictions:
    """Ridge on smoothed *other* fields, fitted on visible rows only."""
    names = sorted(feature_predictions)
    X = np.column_stack([np.asarray(feature_predictions[k], dtype=float) for k in names])
    y = np.asarray(target_visible_values, dtype=float)
    ok = visible_mask & np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    if ok.sum() < max(10, len(names) + 2):
        return Predictions("multifield_ridge", np.full(len(query_xy), np.nan))
    Xt = np.column_stack([np.ones(ok.sum()), X[ok]])
    A = Xt.T @ Xt + alpha * np.eye(Xt.shape[1])
    A[0, 0] -= alpha
    beta = np.linalg.solve(A, Xt.T @ y[ok])
    Xall = np.column_stack([np.ones(len(query_xy)), X])
    return Predictions("multifield_ridge", Xall @ beta)


def geometric_centers(xy: np.ndarray, support: np.ndarray, n_centers: int, seed: int = 0) -> np.ndarray:
    """Deterministic geometry-only centre pick: hash of rounded coordinates."""
    idx = np.flatnonzero(support)
    key = np.array([hash((int(xy[i, 0] // 8), int(xy[i, 1] // 8), seed)) for i in idx], dtype=np.int64)
    order = idx[np.argsort(key)]
    picked: list[np.ndarray] = []
    min_sep = 4.0 * geometry.estimate_pitch(xy[support])
    for i in order:
        if all(np.linalg.norm(xy[i] - xy[j]) >= min_sep for j in picked):
            picked.append(np.array([i]))
        if len(picked) >= n_centers:
            break
    return xy[np.array([p[0] for p in picked])]


def evaluate_section(
    xy: np.ndarray,
    support: np.ndarray,
    field_scores: Mapping[str, np.ndarray],
    *,
    centers: np.ndarray,
    radii_pitch: Sequence[float],
    core_radius_pitch: float = 2.0,
    bandwidth_pitch: float = 2.0,
) -> list[dict[str, Any]]:
    pitch = geometry.estimate_pitch(xy[support])
    rows: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    names = sorted(field_scores)
    for ci, center in enumerate(centers):
        dist_center = np.linalg.norm(xy - center, axis=1)
        core = dist_center <= core_radius_pitch * pitch
        if core.sum() < 5:
            skipped.append({"center_index": ci, "reason": "core has fewer than 5 points", "n_core": int(core.sum())})
            continue
        for radius_pitch in radii_pitch:
            hidden = (dist_center <= radius_pitch * pitch) & support
            visible = support & ~hidden
            if visible.sum() < 50 or hidden.sum() < 10:
                skipped.append(
                    {
                        "center_index": ci,
                        "radius_pitch": float(radius_pitch),
                        "reason": "fewer than 50 visible points or fewer than 10 hidden points",
                        "n_visible": int(visible.sum()),
                        "n_hidden": int(hidden.sum()),
                    }
                )
                continue
            vis_xy = xy[visible]
            q_xy = xy[core]
            smoothed = {
                name: predict_kernel(vis_xy, field_scores[name][visible], xy, k=16).values
                for name in names
            }
            for target in names:
                truth = field_scores[target][core]
                preds = {
                    "background_mean": predict_background(field_scores[target][visible], len(xy)),
                    "idw_k8": predict_idw(vis_xy, field_scores[target][visible], xy),
                    "kernel_k16": Predictions("kernel_k16", smoothed[target]),
                    "kernel_local_bw2": predict_kernel(
                        vis_xy, field_scores[target][visible], xy, k=16, bandwidth=bandwidth_pitch * pitch
                    ),
                }
                others = {k: smoothed[k] for k in names if k != target}
                if others:
                    preds["multifield_ridge"] = predict_multifield(
                        xy, others, field_scores[target], visible
                    )
                for name, pred in preds.items():
                    values = pred.values[core]
                    ok = np.isfinite(values) & np.isfinite(truth)
                    if ok.sum() < 5:
                        continue
                    rows.append(
                        {
                            "center_index": ci,
                            "center_x": float(center[0]),
                            "center_y": float(center[1]),
                            "radius_pitch": float(radius_pitch),
                            "core_radius_pitch": float(core_radius_pitch),
                            "target_field": target,
                            "predictor": name,
                            "n_core": int(ok.sum()),
                            "rmse": float(np.sqrt(np.mean((values[ok] - truth[ok]) ** 2))),
                            "std_truth": float(np.std(truth[ok])),
                            "mae": float(np.mean(np.abs(values[ok] - truth[ok]))),
                            "pearson": float(np.corrcoef(values[ok], truth[ok])[0, 1])
                            if np.std(values[ok]) > 0 and np.std(truth[ok]) > 0
                            else float("nan"),
                            "distance_to_nearest_visible": float(
                                np.min(np.linalg.norm(vis_xy - center, axis=1))
                            ),
                        }
                    )
    return rows, skipped
