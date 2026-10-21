"""Route B: cross-field relation operators (B0 stationary reference vs B1 lag-resolved).

The capability under test is narrow and falsifiable: **can a model tell "same-point
coupling" from "the same coupling displaced by a few pitches"?**  A stationary LMC /
multi-output GP cannot represent a displacement (its cross-covariance depends on
distance only), so it is a legitimate strong reference for prediction but is
expected to return displacement = 0.

B1 uses an explicit 2-D lag basis: partner(s) ~ Σ_l w_l · anchor(s − l) + b.
Weights are fitted by ridge regression on visible data only; the displacement is
read off the weight map. No neural kernel, no per-position free parameters.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


def grid_shift(field: np.ndarray, n: int, dy: int, dx: int) -> np.ndarray:
    """Shift a row-major square field by (dy, dx) with zero fill (no wraparound)."""
    img = field.reshape(n, n)
    out = np.zeros_like(img)
    ys_src = slice(max(0, -dy), n - max(0, dy))
    ys_dst = slice(max(0, dy), n - max(0, -dy))
    xs_src = slice(max(0, -dx), n - max(0, dx))
    xs_dst = slice(max(0, dx), n - max(0, -dx))
    out[ys_dst, xs_dst] = img[ys_src, xs_src]
    return out.ravel()


@dataclass
class LagFit:
    lags: list[tuple[int, int]]
    weights: np.ndarray
    intercept: float
    alpha: float
    prediction: np.ndarray
    displacement: tuple[float, float]
    weight_peak: tuple[int, int]
    weight_map: dict[tuple[int, int], float]

    def to_dict(self) -> dict[str, Any]:
        return {
            "displacement": [float(self.displacement[0]), float(self.displacement[1])],
            "displacement_magnitude": float(np.hypot(*self.displacement)),
            "weight_peak": list(self.weight_peak),
            "alpha": float(self.alpha),
            "intercept": float(self.intercept),
            "top_weights": sorted(
                ({"lag": list(k), "weight": float(v)} for k, v in self.weight_map.items()),
                key=lambda d: -abs(d["weight"]),
            )[:5],
        }


def _design(anchor: np.ndarray, n: int, lags: list[tuple[int, int]]) -> np.ndarray:
    cols = [grid_shift(anchor, n, dy, dx) for dy, dx in lags]
    return np.column_stack(cols)


def fit_lag_operator(
    anchor: np.ndarray,
    partner: np.ndarray,
    n: int,
    *,
    max_lag: int = 12,
    step: int = 4,
    alpha: float = 1.0,
    train_mask: np.ndarray | None = None,
) -> LagFit:
    lags = [(dy, dx) for dy in range(-max_lag, max_lag + 1, step) for dx in range(-max_lag, max_lag + 1, step)]
    X = _design(anchor, n, lags)
    y = partner
    train = np.ones(len(y), dtype=bool) if train_mask is None else np.asarray(train_mask, dtype=bool)
    Xt = X[train]
    yt = y[train]
    A = Xt.T @ Xt + alpha * np.eye(Xt.shape[1])
    b = Xt.T @ yt
    coef = np.linalg.solve(A, b)
    intercept = float(np.mean(yt - Xt @ coef))
    pred = X @ coef + intercept
    weight_map = {lag: float(w) for lag, w in zip(lags, coef)}
    peak = max(weight_map, key=lambda k: abs(weight_map[k]))
    # displacement = positive-weight centroid (more stable than the single argmax)
    pos = [(lag, w) for lag, w in weight_map.items() if w > 0]
    if pos:
        total = sum(w for _, w in pos)
        dy = sum(lag[0] * w for lag, w in pos) / total
        dx = sum(lag[1] * w for lag, w in pos) / total
    else:
        dy = dx = 0.0
    return LagFit(
        lags=lags,
        weights=coef,
        intercept=intercept,
        alpha=alpha,
        prediction=pred,
        displacement=(float(dy), float(dx)),
        weight_peak=peak,
        weight_map=weight_map,
    )


def fit_isotropic_operator(
    anchor: np.ndarray,
    partner: np.ndarray,
    n: int,
    *,
    max_lag: int = 12,
    n_rings: int = 4,
    alpha: float = 1.0,
    train_mask: np.ndarray | None = None,
) -> LagFit:
    """B0 reference: weights depend on distance only (no direction, no displacement).

    Implemented as a lag basis restricted to ring-averaged shifts, so the model can
    smooth/mix but cannot represent a directional offset.
    """
    lags = [(dy, dx) for dy in range(-max_lag, max_lag + 1, 2) for dx in range(-max_lag, max_lag + 1, 2)]
    radii = np.array([np.hypot(dy, dx) for dy, dx in lags])
    edges = np.linspace(0, max_lag * np.sqrt(2), n_rings + 1)
    ring_id = np.clip(np.digitize(radii, edges) - 1, 0, n_rings - 1)
    X = _design(anchor, n, lags)
    R = np.zeros((len(lags), n_rings))
    for k, r in enumerate(ring_id):
        R[k, r] = 1.0
    Xr = X @ R
    train = np.ones(len(partner), dtype=bool) if train_mask is None else np.asarray(train_mask, dtype=bool)
    Xt, yt = Xr[train], partner[train]
    A = Xt.T @ Xt + alpha * np.eye(n_rings)
    coef = np.linalg.solve(A, Xt.T @ yt)
    intercept = float(np.mean(yt - Xt @ coef))
    pred = Xr @ coef + intercept
    weight_map = {lags[k]: float(coef[ring_id[k]]) for k in range(len(lags))}
    return LagFit(
        lags=lags,
        weights=coef,
        intercept=intercept,
        alpha=alpha,
        prediction=pred,
        displacement=(0.0, 0.0),  # stationary by construction
        weight_peak=(0, 0),
        weight_map=weight_map,
    )


def classify_displacement(fit: LagFit, threshold: float = 3.0) -> str:
    return "offset" if np.hypot(*fit.displacement) >= threshold else "same_point"


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def blocked_cv(
    anchor: np.ndarray,
    partner: np.ndarray,
    n: int,
    fit_fn,
    *,
    n_blocks: int = 4,
) -> float:
    """Held-out SSE with quadrant blocks (spatially honest)."""
    yy, xx = np.mgrid[0:n, 0:n]
    blocks = ((xx.ravel() > np.median(xx)) * 2 + (yy.ravel() > np.median(yy))).astype(int)
    total = 0.0
    for b in range(int(blocks.max()) + 1):
        test = blocks == b
        fit = fit_fn(train_mask=~test)
        total += float(np.sum((fit.prediction[test] - partner[test]) ** 2))
    return total


def capability_sweep(
    n: int = 64,
    *,
    sigmas=(3.0, 6.0, 9.0),
    offsets=(5.0, 10.0, 15.0),
    noise_levels=(0.0, 0.15),
    max_lag: int = 20,
    step: int = 4,
    seed: int = 0,
) -> list[dict[str, Any]]:
    """When *can* a lag operator identify a displacement? Sweep field width, shift, noise.

    Recovery criterion: the weight peak lies within 2 pitch of the true shift.
    """
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:n, 0:n]
    xy = np.column_stack([xx.ravel(), yy.ravel()]).astype(float)
    center = np.array([n / 3, n / 2])
    rows: list[dict[str, Any]] = []
    for sigma in sigmas:
        for off in offsets:
            for noise in noise_levels:
                anchor = np.exp(-np.sum((xy - center) ** 2, axis=1) / (2 * sigma ** 2))
                shifted = np.exp(-np.sum((xy - (center + np.array([off, 0.0]))) ** 2, axis=1) / (2 * sigma ** 2))
                partner = shifted + noise * rng.normal(size=len(xy))
                fit = fit_lag_operator(anchor, partner, n, max_lag=max_lag, step=step)
                b0 = fit_isotropic_operator(anchor, partner, n, max_lag=max_lag)
                true_peak = (0, int(round(off / step) * step))
                err = float(np.hypot(fit.weight_peak[0] - true_peak[0], fit.weight_peak[1] - true_peak[1]))
                cv_b0 = blocked_cv(anchor, partner, n, lambda train_mask: fit_isotropic_operator(anchor, partner, n, max_lag=max_lag, train_mask=train_mask))
                cv_b1 = blocked_cv(anchor, partner, n, lambda train_mask: fit_lag_operator(anchor, partner, n, max_lag=max_lag, step=step, train_mask=train_mask))
                rows.append(
                    {
                        "sigma_pitch": sigma,
                        "offset_pitch": off,
                        "noise_sd": noise,
                        "peak_lag": list(fit.weight_peak),
                        "peak_error_pitch": err,
                        "recovered": bool(err <= 2.0),
                        "B0_cv_sse": cv_b0,
                        "B1_cv_sse": cv_b1,
                        "cv_ratio_B1_over_B0": cv_b1 / cv_b0 if cv_b0 > 0 else float("nan"),
                    }
                )
    return rows


def cross_correlation_peak(
    anchor: np.ndarray, partner: np.ndarray, n: int, *, max_lag: int = 14, min_overlap: int = 500
) -> dict[str, Any]:
    """Direct 2-D cross-correlation peak: the displacement, without any lag regression.

    This is the simplest estimator that can see a shift when the fields are broad and
    heavily overlapping.  It is reported alongside the ridge/NNLS operators because
    those can put their largest weight at lag 0 even when this peak is non-zero.
    """
    a = np.asarray(anchor, dtype=float).reshape(n, n)
    b = np.asarray(partner, dtype=float).reshape(n, n)
    best = {"corr": -2.0, "dy": 0, "dx": 0}
    for dy in range(-max_lag, max_lag + 1):
        for dx in range(-max_lag, max_lag + 1):
            a_sub = a[max(0, dy):n + min(0, dy), max(0, dx):n + min(0, dx)]
            b_sub = b[max(0, -dy):n + min(0, -dy), max(0, -dx):n + min(0, -dx)]
            if a_sub.size < min_overlap:
                continue
            if np.std(a_sub) == 0 or np.std(b_sub) == 0:
                continue
            c = float(np.corrcoef(a_sub.ravel(), b_sub.ravel())[0, 1])
            if c > best["corr"]:
                best = {"corr": c, "dy": dy, "dx": dx}
    best["magnitude"] = float(np.hypot(best["dy"], best["dx"]))
    return best
