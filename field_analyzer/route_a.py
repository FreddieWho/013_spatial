"""Route A: finite deformable geometry (A0 radial reference vs A1 deformable).

A0 is a *legitimate strong reference*, not a straw man:
  * ring: circular ring (centre, radius, thickness, amplitude)
  * blob: one Gaussian (centre, sigma, amplitude)
A1 adds a bounded amount of shape freedom:
  * ring: radius/thickness modulated by angular orders 1 and 2 (Fourier)
  * blob: up to three Gaussian peaks, count chosen by BIC

Both are fitted on visible data only.  Truth is used for scoring, never for
choosing the anchor, the number of peaks, or the initialisation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

import numpy as np
from scipy.optimize import least_squares

MAX_ATOMS = 3


@dataclass
class FitResult:
    model: str
    params: dict[str, Any]
    n_params: int
    visible_sse: float
    bic: float
    prediction: np.ndarray = field(repr=False, default=None)


# --------------------------------------------------------------------------- #
# models on the 2-D grid
# --------------------------------------------------------------------------- #
def ring_model(xy: np.ndarray, p: np.ndarray, deformable: bool) -> np.ndarray:
    cx, cy, r0, t0, amp = p[:5]
    if deformable:
        c1c, c1s, c2c, c2s, d1c, d1s, d2c, d2s = p[5:13]
    else:
        c1c = c1s = c2c = c2s = d1c = d1s = d2c = d2s = 0.0
    dx = xy[:, 0] - cx
    dy = xy[:, 1] - cy
    d = np.hypot(dx, dy)
    theta = np.arctan2(dy, dx)
    radial_mod = 1.0 + c1c * np.cos(theta) + c1s * np.sin(theta) + c2c * np.cos(2 * theta) + c2s * np.sin(2 * theta)
    thick_mod = 1.0 + d1c * np.cos(theta) + d1s * np.sin(theta) + d2c * np.cos(2 * theta) + d2s * np.sin(2 * theta)
    radial_mod = np.clip(radial_mod, 0.2, 3.0)
    thick_mod = np.clip(thick_mod, 0.2, 3.0)
    r = np.maximum(r0 * radial_mod, 1e-6)
    t = np.maximum(t0 * thick_mod, 1e-6)
    return amp * np.exp(-0.5 * ((d - r) / t) ** 2)


def gaussian_sum_model(xy: np.ndarray, p: np.ndarray, k: int) -> np.ndarray:
    out = np.zeros(xy.shape[0], dtype=np.float64)
    for i in range(k):
        cx, cy, sigma, amp = p[4 * i:4 * i + 4]
        d2 = (xy[:, 0] - cx) ** 2 + (xy[:, 1] - cy) ** 2
        out += amp * np.exp(-d2 / (2.0 * max(sigma, 1e-3) ** 2))
    return out


# --------------------------------------------------------------------------- #
# fitting
# --------------------------------------------------------------------------- #
def _initial_centers(field: np.ndarray, xy: np.ndarray, k: int) -> list[np.ndarray]:
    """Coarse peak picks from the *visible* field (no truth involved)."""
    order = np.argsort(field)[::-1]
    picked: list[np.ndarray] = []
    for idx in order:
        c = xy[idx]
        if all(np.linalg.norm(c - p) > 4.0 for p in picked):
            picked.append(c)
        if len(picked) >= k:
            break
    return picked


def _bounds_and_start_ring(field: np.ndarray, xy: np.ndarray, deformable: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    peak = xy[int(np.argmax(field))]
    center0 = np.array([peak[0], peak[1]])
    d = np.linalg.norm(xy - center0, axis=1)
    r0 = float(np.percentile(d, 25) + 4.0)
    t0 = 3.0
    p0 = np.array([center0[0], center0[1], r0, t0, max(float(np.max(field)), 1e-6)])
    lower = np.array([xy[:, 0].min(), xy[:, 1].min(), 2.0, 0.5, 0.0])
    upper = np.array([xy[:, 0].max(), xy[:, 1].max(), d.max(), 20.0, 10.0 * max(float(np.max(field)), 1.0)])
    if deformable:
        p0 = np.r_[p0, np.zeros(8)]
        lower = np.r_[lower, -0.9 * np.ones(8)]
        upper = np.r_[upper, 0.9 * np.ones(8)]
    return p0, lower, upper


def fit_ring(field: np.ndarray, xy: np.ndarray, *, deformable: bool, valid: np.ndarray | None = None) -> FitResult:
    valid = np.ones(len(field), dtype=bool) if valid is None else valid
    y = field[valid]
    x = xy[valid]
    p0, lb, ub = _bounds_and_start_ring(y if len(y) else field, x, deformable)
    p0 = np.clip(p0, lb + 1e-6, ub - 1e-6)

    def residual(p: np.ndarray) -> np.ndarray:
        return ring_model(x, p, deformable) - y

    sol = least_squares(residual, p0, bounds=(lb, ub), max_nfev=4000, xtol=1e-10, ftol=1e-10)
    pred = ring_model(xy, sol.x, deformable)
    n_params = 13 if deformable else 5
    sse = float(np.sum(sol.fun ** 2))
    n = max(int(valid.sum()), 1)
    bic = n * np.log(max(sse / n, 1e-12)) + n_params * np.log(n)
    params = {
        "cx": float(sol.x[0]),
        "cy": float(sol.x[1]),
        "radius": float(sol.x[2]),
        "thickness": float(sol.x[3]),
        "amplitude": float(sol.x[4]),
    }
    if deformable:
        params["radial_fourier"] = [float(v) for v in sol.x[5:9]]
        params["thickness_fourier"] = [float(v) for v in sol.x[9:13]]
    return FitResult("A1_ring" if deformable else "A0_ring", params, n_params, sse, float(bic), pred)


def fit_gaussians(field: np.ndarray, xy: np.ndarray, *, k: int, valid: np.ndarray | None = None) -> FitResult:
    valid = np.ones(len(field), dtype=bool) if valid is None else valid
    y = field[valid]
    x = xy[valid]
    centers = _initial_centers(np.where(valid, field, -np.inf), xy, k)
    while len(centers) < k:
        centers.append(xy[np.random.default_rng(0).integers(0, len(xy))])
    p0: list[float] = []
    lb: list[float] = []
    ub: list[float] = []
    amp0 = max(float(np.percentile(y, 95)), 1e-6)
    for c in centers:
        p0.extend([c[0], c[1], 5.0, amp0 / k])
        lb.extend([xy[:, 0].min(), xy[:, 1].min(), 0.5, 0.0])
        ub.extend([xy[:, 0].max(), xy[:, 1].max(), 30.0, 10.0 * amp0])

    def residual(p: np.ndarray) -> np.ndarray:
        return gaussian_sum_model(x, p, k) - y

    sol = least_squares(residual, np.array(p0), bounds=(np.array(lb), np.array(ub)), max_nfev=6000)
    pred = gaussian_sum_model(xy, sol.x, k)
    n_params = 4 * k
    sse = float(np.sum(sol.fun ** 2))
    n = max(int(valid.sum()), 1)
    bic = n * np.log(max(sse / n, 1e-12)) + n_params * np.log(n)
    peaks = [
        {"x": float(sol.x[4 * i]), "y": float(sol.x[4 * i + 1]), "sigma": float(sol.x[4 * i + 2]),
         "amp": float(sol.x[4 * i + 3])}
        for i in range(k)
    ]
    peaks.sort(key=lambda d: -d["amp"])
    return FitResult(f"A{1 if k > 1 else 0}_gauss{k}", {"peaks": peaks}, n_params, sse, float(bic), pred)


def spatial_blocks(xy: np.ndarray, n_blocks: int = 4) -> np.ndarray:
    """Quadrant blocks (median split on x and y) for honest spatial cross-validation."""
    mx = np.median(xy[:, 0])
    my = np.median(xy[:, 1])
    bx = (xy[:, 0] > mx).astype(int)
    by = (xy[:, 1] > my).astype(int)
    return bx * 2 + by


def _cv_sse(field: np.ndarray, xy: np.ndarray, fit_fn, n_blocks: int = 4) -> float:
    blocks = spatial_blocks(xy, n_blocks)
    total = 0.0
    for b in range(int(blocks.max()) + 1):
        test = blocks == b
        if test.sum() < 10:
            continue
        train = ~test
        pred = fit_fn(train)
        total += float(np.sum((pred[test] - field[test]) ** 2))
    return total


@dataclass
class Candidate:
    name: str
    kind: str  # "ring" | "blob"
    k: int
    params: dict[str, Any]
    prediction: np.ndarray
    in_sse: float
    cv_sse: float
    n_params: int
    bic: float


def _ring_candidate(field, xy, deformable: bool) -> Candidate:
    fit = fit_ring(field, xy, deformable=deformable)

    def train_predict(train_mask: np.ndarray) -> np.ndarray:
        sub = fit_ring(field, xy, deformable=deformable, valid=train_mask)
        return sub.prediction

    cv = _cv_sse(field, xy, train_predict)
    return Candidate(fit.model, "ring", 1, fit.params, fit.prediction, fit.visible_sse, cv, fit.n_params, fit.bic)


def _blob_candidate(field, xy, k: int) -> Candidate:
    fit = fit_gaussians(field, xy, k=k)

    def train_predict(train_mask: np.ndarray) -> np.ndarray:
        sub = fit_gaussians(field, xy, k=k, valid=train_mask)
        return sub.prediction

    cv = _cv_sse(field, xy, train_predict)
    return Candidate(fit.model, "blob", k, fit.params, fit.prediction, fit.visible_sse, cv, fit.n_params, fit.bic)


def candidates_for(kind: str, field: np.ndarray, xy: np.ndarray, *, k_max: int = MAX_ATOMS) -> list[Candidate]:
    if kind == "ring":
        return [_ring_candidate(field, xy, False), _ring_candidate(field, xy, True)]
    out = [_blob_candidate(field, xy, k) for k in range(1, k_max + 1)]
    return out


def select_by_cv(candidates: Sequence[Candidate]) -> Candidate:
    """Prefer the simplest model whose blocked-CV error is within 5% of the best.

    BIC is reported but not used for selection: residuals are spatially correlated,
    so an iid-penalty over-selects deformable / multi-peak models.
    """
    best = min(candidates, key=lambda c: c.cv_sse)
    tol = 1.05 * best.cv_sse
    for cand in sorted(candidates, key=lambda c: (c.n_params, c.cv_sse)):
        if cand.cv_sse <= tol:
            return cand
    return best


# --------------------------------------------------------------------------- #
# scoring (truth used here only)
# --------------------------------------------------------------------------- #
def scaled_rmse(pred: np.ndarray, truth: np.ndarray) -> float:
    denom = float(pred @ pred)
    scale = float(pred @ truth) / denom if denom > 0 else 0.0
    return float(np.sqrt(np.mean((scale * pred - truth) ** 2)))


def pearson(pred: np.ndarray, truth: np.ndarray) -> float:
    if np.std(pred) == 0 or np.std(truth) == 0:
        return float("nan")
    return float(np.corrcoef(pred, truth)[0, 1])


def angular_amplitude(params: dict[str, Any]) -> float:
    fourier = params.get("radial_fourier")
    if not fourier:
        return 0.0
    return float(np.hypot(fourier[2], fourier[3])) if len(fourier) >= 4 else 0.0


def centered(field: np.ndarray, xy: np.ndarray) -> np.ndarray:
    total = float(np.sum(field))
    if total <= 0:
        return xy.mean(axis=0)
    return np.sum(xy * field[:, None], axis=0) / total


def peak_coords(field: np.ndarray, xy: np.ndarray, *, frac: float = 0.5, min_sep: float = 4.0) -> np.ndarray:
    """Genuine local maxima (3x3 neighbourhood) above ``frac`` of the maximum.

    The grid is a square lattice stored row-major, so the field reshapes to
    (n, n) and local maxima are found on that lattice.
    """
    n = int(round(np.sqrt(len(field))))
    if n * n != len(field):
        raise ValueError("peak_coords expects a square lattice")
    img = field.reshape(n, n)
    peak = img.max()
    if peak <= 0:
        return np.empty((0, 2))
    padded = np.pad(img, 1, mode="constant", constant_values=-np.inf)
    is_local = np.ones_like(img, dtype=bool)
    for di in (0, 1, 2):
        for dj in (0, 1, 2):
            if di == 1 and dj == 1:
                continue
            is_local &= img >= padded[di:di + n, dj:dj + n] - 1e-12
    is_local &= img > frac * peak
    order = np.argsort(img[is_local])[::-1]
    coords = np.argwhere(is_local)[order]
    picked: list[np.ndarray] = []
    for (i, j) in coords:
        c = xy[i * n + j]
        if all(np.linalg.norm(c - p) >= min_sep for p in picked):
            picked.append(c)
    return np.asarray(picked) if picked else np.empty((0, 2))


def peak_distances(pred: np.ndarray, truth: np.ndarray, xy: np.ndarray) -> dict[str, float]:
    p = peak_coords(pred, xy)
    t = peak_coords(truth, xy)
    if len(p) == 0 or len(t) == 0:
        return {"pred_to_truth": float("nan"), "truth_to_pred": float("nan"), "n_pred": len(p), "n_truth": len(t)}
    d_pt = [float(np.min(np.linalg.norm(t - c, axis=1))) for c in p]
    d_tp = [float(np.min(np.linalg.norm(p - c, axis=1))) for c in t]
    return {
        "pred_to_truth": float(np.mean(d_pt)),
        "truth_to_pred": float(np.mean(d_tp)),
        "n_pred": len(p),
        "n_truth": len(t),
    }
