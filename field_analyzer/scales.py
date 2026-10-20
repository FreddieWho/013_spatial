"""Scale and cross-field relation measurability (FA-01).

Deliberately separate quantities (05_COMMON_PROTOCOL §7, plan §4):
  * raw field variogram                      -> the field's own correlation range
  * residual variogram after a smooth trend   -> a *residual* scale, never the field's
  * directional variogram                     -> anisotropy, only where pairs exist
  * cross-field variogram / correlation       -> relation strength vs lag
  * empirical recoverable range               -> (FA-05)

States instead of guesses:
  OK                         plateau reached inside the field of view
  NOT_IDENTIFIED_WITHIN_FOV  no plateau; the fitted range exceeds the observed span
  UNIDENTIFIABLE             spatial variance ~0, so a length parameter is meaningless
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np
from scipy.optimize import curve_fit
from scipy.spatial import cKDTree

OK = "OK"
NOT_IDENTIFIED = "NOT_IDENTIFIED_WITHIN_FOV"
UNIDENTIFIABLE = "UNIDENTIFIABLE"
NOT_MEASURABLE = "NOT_MEASURABLE"

DIRECTION_LABELS = ("0-45", "45-90", "90-135", "135-180")
MIN_PAIRS_FIT = 30
MIN_PAIRS_DIRECTION = 30
MIN_USABLE_BINS = 5


@dataclass
class PairSet:
    i: np.ndarray
    j: np.ndarray
    distance: np.ndarray
    angle_mod_pi: np.ndarray


def pairs_within(xy: np.ndarray, valid: np.ndarray, max_distance: float) -> PairSet:
    """All unordered pairs with both endpoints valid and distance <= max_distance."""
    xy = np.asarray(xy, dtype=np.float64)
    valid = np.asarray(valid, dtype=bool)
    idx = np.flatnonzero(valid)
    if len(idx) < 3:
        raise ValueError("need at least three valid points")
    tree = cKDTree(xy[idx])
    raw = tree.query_pairs(max_distance, output_type="ndarray")
    if len(raw) == 0:
        return PairSet(
            i=np.empty(0, dtype=np.int64),
            j=np.empty(0, dtype=np.int64),
            distance=np.empty(0),
            angle_mod_pi=np.empty(0),
        )
    mapped = idx[raw]
    delta = xy[mapped[:, 1]] - xy[mapped[:, 0]]
    d = np.linalg.norm(delta, axis=1)
    ang = np.mod(np.arctan2(delta[:, 1], delta[:, 0]), np.pi)
    return PairSet(i=mapped[:, 0], j=mapped[:, 1], distance=d, angle_mod_pi=ang)


def _bin_index(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    idx = np.digitize(values, edges) - 1
    idx = np.clip(idx, 0, len(edges) - 2)
    return idx


def empirical_variogram(
    values: np.ndarray, pairs: PairSet, edges: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """gamma(h) = 0.5 * mean((v_i - v_j)^2) per distance bin, with pair counts."""
    v = np.asarray(values, dtype=np.float64)
    n_bins = len(edges) - 1
    if len(pairs.distance) == 0:
        return np.full(n_bins, np.nan), np.zeros(n_bins), np.zeros(n_bins)
    diff = v[pairs.i] - v[pairs.j]
    squared = diff ** 2
    finite = np.isfinite(squared)
    b = _bin_index(pairs.distance, edges)
    counts = np.bincount(b[finite], minlength=n_bins).astype(np.float64)
    sums = np.bincount(b[finite], weights=squared[finite], minlength=n_bins)
    gamma = np.where(counts > 0, 0.5 * sums / np.maximum(counts, 1), np.nan)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return gamma, counts, centers


def cross_field_curves(
    v1: np.ndarray, v2: np.ndarray, pairs: PairSet, edges: np.ndarray
) -> dict[str, np.ndarray]:
    """Co-variogram plus a *normalised cross-covariance* of two fields.

    ``cross_corr(h)`` is the symmetrised, standardised cross-covariance
    ``0.5 * (z1_i z2_j + z1_j z2_i)`` averaged per lag bin.  At h -> 0 it equals the
    same-point correlation; a peak at h > 0 is evidence of a spatial offset between
    the fields.  The co-variogram gamma12(h) is reported in field units as a
    secondary quantity.  The correlation of *increments* is deliberately not
    reported: it is easy to misread as a field-field correlation.
    """
    a = np.asarray(v1, dtype=np.float64)
    b = np.asarray(v2, dtype=np.float64)
    n_bins = len(edges) - 1
    if len(pairs.distance) == 0:
        nan = np.full(n_bins, np.nan)
        return {"cross_corr": nan, "cross_gamma": nan, "counts": np.zeros(n_bins), "centers": nan}
    finite_values = np.isfinite(a) & np.isfinite(b)
    if finite_values.sum() < 10:
        nan = np.full(n_bins, np.nan)
        return {"cross_corr": nan, "cross_gamma": nan, "counts": np.zeros(n_bins), "centers": nan}
    mu_a, mu_b = np.nanmean(a[finite_values]), np.nanmean(b[finite_values])
    sd_a, sd_b = np.nanstd(a[finite_values]), np.nanstd(b[finite_values])
    z1 = (a - mu_a) / sd_a if sd_a > 0 else np.zeros_like(a)
    z2 = (b - mu_b) / sd_b if sd_b > 0 else np.zeros_like(b)
    centred = (
        np.isfinite(z1[pairs.i]) & np.isfinite(z2[pairs.j])
        & np.isfinite(z1[pairs.j]) & np.isfinite(z2[pairs.i])
    )
    sym = 0.5 * (z1[pairs.i] * z2[pairs.j] + z1[pairs.j] * z2[pairs.i])
    raw = 0.5 * (a[pairs.i] - a[pairs.j]) * (b[pairs.i] - b[pairs.j])
    b_idx = _bin_index(pairs.distance, edges)
    counts = np.bincount(b_idx[centred], minlength=n_bins).astype(float)
    corr_sum = np.bincount(b_idx[centred], weights=sym[centred], minlength=n_bins)
    gamma_ok = np.isfinite(raw)
    gamma_sum = np.bincount(b_idx[gamma_ok], weights=raw[gamma_ok], minlength=n_bins)
    gamma_cnt = np.bincount(b_idx[gamma_ok], minlength=n_bins).astype(float)
    cross_corr = np.where(counts > 0, corr_sum / np.maximum(counts, 1), np.nan)
    cross_gamma = np.where(gamma_cnt > 0, gamma_sum / np.maximum(gamma_cnt, 1), np.nan)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return {"cross_corr": cross_corr, "cross_gamma": cross_gamma, "counts": counts, "centers": centers}


def directional_variogram(
    values: np.ndarray, pairs: PairSet, edges: np.ndarray, n_directions: int = 4
) -> dict[str, np.ndarray]:
    v = np.asarray(values, dtype=np.float64)
    n_bins = len(edges) - 1
    gamma = np.full((n_directions, n_bins), np.nan)
    counts = np.zeros((n_directions, n_bins))
    if len(pairs.distance) == 0:
        return {"gamma": gamma, "counts": counts}
    diff = (v[pairs.i] - v[pairs.j]) ** 2
    finite = np.isfinite(diff)
    b = _bin_index(pairs.distance, edges)
    sector = np.minimum((pairs.angle_mod_pi / (np.pi / n_directions)).astype(int), n_directions - 1)
    for s in range(n_directions):
        sel = finite & (sector == s)
        cnt = np.bincount(b[sel], minlength=n_bins).astype(float)
        sums = np.bincount(b[sel], weights=diff[sel], minlength=n_bins)
        gamma[s] = np.where(cnt > 0, 0.5 * sums / np.maximum(cnt, 1), np.nan)
        counts[s] = cnt
    return {"gamma": gamma, "counts": counts}


def _exp_model(h: np.ndarray, nugget: float, sill: float, rng: float) -> np.ndarray:
    return nugget + sill * (1.0 - np.exp(-h / rng))


def fit_variogram(centers: np.ndarray, gamma: np.ndarray, counts: np.ndarray) -> dict[str, Any]:
    """Fit nugget+sill+exponential range, or say why no range can be reported."""
    finite = np.isfinite(gamma) & (counts >= MIN_PAIRS_FIT)
    if finite.sum() < MIN_USABLE_BINS:
        return {
            "status": NOT_MEASURABLE,
            "reason": f"only {int(finite.sum())} usable lag bins (need {MIN_USABLE_BINS})",
            "nugget": None,
            "sill": None,
            "range": None,
            "n_bins_used": int(finite.sum()),
        }
    h = centers[finite]
    g = gamma[finite]
    w = np.sqrt(counts[finite])
    span = float(h.max())
    total_variance = float(np.nanvar(g)) + float(np.nanmean(g))
    if not np.isfinite(total_variance) or total_variance <= 1e-12:
        return {
            "status": UNIDENTIFIABLE,
            "reason": "spatial variance is numerically zero; a length parameter has no meaning",
            "nugget": None,
            "sill": None,
            "range": None,
            "n_bins_used": int(finite.sum()),
        }
    init = [max(float(g.min()), 0.0), max(float(g.max() - g.min()), 1e-9), span / 2.0]
    bounds = ([0.0, 0.0, span / 50.0], [np.inf, np.inf, 50.0 * span])
    try:
        popt, _ = curve_fit(
            _exp_model,
            h,
            g,
            p0=init,
            sigma=1.0 / np.maximum(w, 1e-9),
            bounds=bounds,
            maxfev=20000,
        )
    except Exception as exc:  # pragma: no cover - numerical failure path
        return {
            "status": NOT_MEASURABLE,
            "reason": f"fit failed: {type(exc).__name__}",
            "nugget": None,
            "sill": None,
            "range": None,
            "n_bins_used": int(finite.sum()),
        }
    nugget, sill, rng = (float(x) for x in popt)
    plateau = nugget + 0.8 * sill
    reached = float(np.nanmax(g)) >= plateau
    if sill <= 1e-9 or (nugget > 0 and sill <= 0.05 * nugget):
        # a sill that is negligible next to the nugget means there is no spatial
        # structure above noise; the length parameter is then meaningless and this
        # must be reported as unidentifiable rather than as "no plateau".
        status = UNIDENTIFIABLE
        reason = (
            f"fitted sill {sill:.3g} is negligible relative to nugget {nugget:.3g}; "
            "no spatial structure above noise"
        )
    elif not reached or rng >= 0.8 * span:
        status = NOT_IDENTIFIED
        reason = (
            f"no plateau inside FOV: max gamma {float(np.nanmax(g)):.3g} < 80% of sill {plateau:.3g}"
            if not reached
            else f"fitted range {rng:.3g} >= 80% of observed span {span:.3g}"
        )
    else:
        status = OK
        reason = "plateau reached inside the field of view"
    return {
        "status": status,
        "reason": reason,
        "nugget": nugget,
        "sill": sill,
        "range": rng,
        "n_bins_used": int(finite.sum()),
        "observed_span": span,
        "max_gamma": float(np.nanmax(g)),
    }


def quadratic_trend(xy: np.ndarray, values: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Smooth global trend (quadratic surface) fitted on valid points only."""
    xy = np.asarray(xy, dtype=np.float64)
    v = np.asarray(values, dtype=np.float64)
    valid = np.asarray(valid, dtype=bool)
    sel = valid & np.isfinite(v)
    if sel.sum() < 6:
        return np.full_like(v, np.nan, dtype=np.float64)
    x = xy[sel, 0] - xy[sel, 0].mean()
    y = xy[sel, 1] - xy[sel, 1].mean()
    X = np.column_stack([np.ones(sel.sum()), x, y, x * x, x * y, y * y])
    beta, *_ = np.linalg.lstsq(X, v[sel], rcond=None)
    xt = xy[:, 0] - xy[sel, 0].mean()
    yt = xy[:, 1] - xy[sel, 1].mean()
    Xall = np.column_stack([np.ones(len(xy)), xt, yt, xt * xt, xt * yt, yt * yt])
    return Xall @ beta


def lag_edges(pitch: float, max_lag_pitch: float, n_bins: int) -> np.ndarray:
    if not np.isfinite(pitch) or pitch <= 0:
        raise ValueError("pitch must be positive and finite")
    return np.linspace(0.0, max_lag_pitch * pitch, n_bins + 1)


def summarise_curve(gamma: np.ndarray, counts: np.ndarray, centers: np.ndarray, pitch: float) -> list[dict[str, Any]]:
    rows = []
    for k in range(len(gamma)):
        rows.append(
            {
                "lag_pitch": float(centers[k] / pitch) if np.isfinite(pitch) else float("nan"),
                "gamma": float(gamma[k]) if np.isfinite(gamma[k]) else None,
                "n_pairs": int(counts[k]),
                "measurable": bool(np.isfinite(gamma[k]) and counts[k] >= MIN_PAIRS_FIT),
            }
        )
    return rows
