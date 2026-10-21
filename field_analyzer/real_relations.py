"""FA-08A/08B: real-data displacement and cross-unit holdout.

Two questions that follow directly from FA-04 (A/B molecular validation) and
FA-07 (depth conditioning / estimator correction):

FA-08A  displacement on real sections
    For each of the six field pairs compute the 2-D cross-correlation peak between
    irregularly spaced spots (kd-tree, no rasterisation), and ask whether a non-zero
    displacement is identifiable at all.  Reported for raw scores and for
    depth-residualised scores, and checked for A/B-half agreement.

FA-08B  cross-unit holdout of the surviving relations
    Three relations survived depth conditioning (EMT~IFN, B_AXIS~EMT, B_AXIS~IFN).
    Fit a spot-level linear map on 7 units, test on the held-out unit, and compare
    against the trivial predictor (section mean = 0 after within-section z-scoring).
    Leave-one-unit-out is the strongest generalisation test the 8 units allow.

No human labels; nothing here creates patient-level claims.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.spatial import cKDTree

from . import geometry, real_validation as rv
from .depth_conditioning import residualise
from .pipeline import build_inputs


def zscore(values: np.ndarray) -> np.ndarray:
    v = np.asarray(values, dtype=float)
    ok = np.isfinite(v)
    if ok.sum() < 10:
        return np.full_like(v, np.nan)
    mu, sd = float(np.mean(v[ok])), float(np.std(v[ok]))
    if sd <= 0:
        return np.full_like(v, np.nan)
    out = np.full_like(v, np.nan)
    out[ok] = (v[ok] - mu) / sd
    return out


@dataclass
class SectionScores:
    section_id: str
    unit_id: str
    xy: np.ndarray
    visible: np.ndarray
    a: dict[str, np.ndarray]
    b: dict[str, np.ndarray]
    log_lib: np.ndarray


def section_scores(
    table: Any, fields: Mapping[str, Sequence[str]], splits: Mapping[str, Any], *, unit_id: str, section_id: str
) -> SectionScores:
    xy = np.asarray(table.xy, dtype=float)
    support = np.ones(table.n_obs, dtype=bool)
    if "in_tissue" in table.aux:
        support = np.asarray(table.aux["in_tissue"]).astype(bool)
    spec_a, spec_b = rv.half_specs(fields, splits)
    inp_a = build_inputs(table.counts, table.genes, xy, support, support, spec_a)
    inp_b = build_inputs(table.counts, table.genes, xy, support, support, spec_b)
    a = {fid: inp_a.arrays[f"field::{fid}"] for fid in spec_a.fields}
    b = {fid: inp_b.arrays[f"field::{fid}"] for fid in spec_b.fields}
    return SectionScores(
        section_id=section_id,
        unit_id=unit_id,
        xy=xy,
        visible=support,
        a=a,
        b=b,
        log_lib=np.log1p(inp_a.arrays["library"]),
    )


# --------------------------------------------------------------------------- #
# FA-08A: displacement on irregular spots
# --------------------------------------------------------------------------- #
def displacement_curve(
    xy: np.ndarray,
    source: np.ndarray,
    target: np.ndarray,
    *,
    pitch: float,
    max_lag_pitch: float = 12.0,
    step_pitch: float = 2.0,
    tolerance_pitch: float = 0.6,
    min_pairs: int = 200,
    null_rotate_deg: float | None = None,
) -> dict[str, Any]:
    """2-D cross-correlation over displacement vectors, on irregular spots.

    For each candidate displacement d, every source spot is matched to the nearest
    target spot at s + d (within ``tolerance_pitch``); the correlation over matched
    pairs is recorded.  Nearest-neighbour matching avoids rasterising the tissue.
    """
    src_xy = np.asarray(xy, dtype=float)
    s = np.asarray(source, dtype=float)
    t = np.asarray(target, dtype=float)
    ok = np.isfinite(s) & np.isfinite(t)
    if ok.sum() < min_pairs:
        return {"status": "NOT_MEASURABLE", "reason": "too few finite spots"}
    src_xy, s, t = src_xy[ok], s[ok], t[ok]
    if null_rotate_deg is not None:
        # null control: keep the target field's own smoothness but destroy its
        # alignment with the source by rotating the target's spatial layout.
        centre = src_xy.mean(axis=0)
        ang = np.deg2rad(null_rotate_deg)
        rot = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
        target_xy = (src_xy - centre) @ rot.T + centre
    else:
        target_xy = src_xy
    tree = cKDTree(target_xy)
    max_lag = max_lag_pitch * pitch
    step = step_pitch * pitch
    tol = tolerance_pitch * pitch
    lags = np.arange(-max_lag, max_lag + 1e-9, step)
    corr0 = float(np.corrcoef(s, t)[0, 1]) if np.std(s) > 0 and np.std(t) > 0 else float("nan")

    best = {"corr": -2.0, "dy": 0.0, "dx": 0.0, "n": 0}
    grid = np.full((len(lags), len(lags)), np.nan)
    for iy, dy in enumerate(lags):
        for ix, dx in enumerate(lags):
            shifted = src_xy + np.array([dx, dy])
            dist, idx = tree.query(shifted, k=1, distance_upper_bound=tol)
            matched = np.isfinite(dist)
            if matched.sum() < min_pairs:
                continue
            a_vals = s[matched]
            b_vals = t[idx[matched]]
            if np.std(a_vals) == 0 or np.std(b_vals) == 0:
                continue
            c = float(np.corrcoef(a_vals, b_vals)[0, 1])
            grid[iy, ix] = c
            if c > best["corr"]:
                best = {"corr": c, "dy": float(dy), "dx": float(dx), "n": int(matched.sum())}
    magnitude = float(np.hypot(best["dx"], best["dy"]))
    # Convention: the search shifts the SOURCE by (dx, dy) and matches the target at
    # the same coordinates, so the reported vector is the displacement of the TARGET
    # relative to the SOURCE (negated).  Positive dx therefore means "the partner
    # field sits at +dx relative to the reference field".
    return {
        "status": "OK",
        "corr_at_zero": corr0,
        "peak_corr": best["corr"],
        "peak_dx_pitch": -best["dx"] / pitch,
        "peak_dy_pitch": -best["dy"] / pitch,
        "peak_magnitude_pitch": magnitude / pitch,
        "peak_n_pairs": best["n"],
        "contrast_vs_zero": best["corr"] - corr0,
        "grid": grid,
    }


def displacement_rows(scores: SectionScores, *, conditioned: bool, **kwargs: Any) -> list[dict[str, Any]]:
    pitch = geometry.estimate_pitch(scores.xy[scores.visible])
    fields = sorted(scores.a)
    rows: list[dict[str, Any]] = []
    a_scores = dict(scores.a)
    b_scores = dict(scores.b)
    if conditioned:
        a_scores = {f: residualise(a_scores[f], scores.log_lib)[0] for f in a_scores}
        b_scores = {f: residualise(b_scores[f], scores.log_lib)[0] for f in b_scores}
    for i in range(len(fields)):
        for j in range(i + 1, len(fields)):
            f1, f2 = fields[i], fields[j]
            for half, table_scores in (("A", a_scores), ("B", b_scores)):
                res = displacement_curve(
                    scores.xy, table_scores[f1], table_scores[f2], pitch=pitch, **kwargs
                )
                row = {
                    "unit_id": scores.unit_id,
                    "section_id": scores.section_id,
                    "field_pair": f"{f1}~{f2}",
                    "half": half,
                    "conditioned": conditioned,
                    "pitch": pitch,
                }
                row.update({k: v for k, v in res.items() if k != "grid"})
                rows.append(row)
    return rows


def summarise_displacement(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    import collections

    pairs = sorted({r["field_pair"] for r in rows})
    out: dict[str, Any] = {"by_pair": {}}
    for pair in pairs:
        entry: dict[str, Any] = {}
        for half in ("A", "B"):
            sel = [r for r in rows if r["field_pair"] == pair and r["half"] == half and r.get("status") == "OK"]
            if not sel:
                entry[half] = {"n": 0}
                continue
            # unit-equal: per-unit median first
            by_unit: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
            for r in sel:
                by_unit[r["unit_id"]].append(r)
            def unit_median(key: str) -> float:
                vals = [float(np.median([r[key] for r in rs])) for rs in by_unit.values()]
                return float(np.median(vals)) if vals else float("nan")
            entry[half] = {
                "n_sections": len(sel),
                "peak_magnitude_pitch": unit_median("peak_magnitude_pitch"),
                "contrast_vs_zero": unit_median("contrast_vs_zero"),
                "corr_at_zero": unit_median("corr_at_zero"),
                "contrast_positive_sections": int(
                    sum(1 for r in sel if r.get("contrast_vs_zero", 0) > 0.02)
                ),
            }
        out["by_pair"][pair] = entry
    return out


# --------------------------------------------------------------------------- #
# FA-08B: leave-one-unit-out relation holdout
# --------------------------------------------------------------------------- #
@dataclass
class LouoResult:
    target: str
    predictor: str
    held_unit: str
    n_train_sections: int
    r2_test: float
    r2_within: float
    r_half_agreement: float
    slope: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "target": self.target,
            "predictor": self.predictor,
            "held_unit": self.held_unit,
            "n_train_sections": self.n_train_sections,
            "r2_test": self.r2_test,
            "r2_within": self.r2_within,
            "r_half_agreement": self.r_half_agreement,
            "slope": self.slope,
        }


def _fit_linear(x: np.ndarray, y: np.ndarray, alpha: float = 1.0) -> tuple[float, float]:
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 10:
        return float("nan"), float("nan")
    X = np.column_stack([np.ones(ok.sum()), x[ok]])
    A = X.T @ X + alpha * np.eye(2)
    A[0, 0] -= alpha
    coef = np.linalg.solve(A, X.T @ y[ok])
    return float(coef[0]), float(coef[1])


def _r2(pred: np.ndarray, truth: np.ndarray) -> float:
    ok = np.isfinite(pred) & np.isfinite(truth)
    if ok.sum() < 10:
        return float("nan")
    ss_res = float(np.sum((pred[ok] - truth[ok]) ** 2))
    ss_tot = float(np.sum((truth[ok] - np.mean(truth[ok])) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def louo_holdout(
    sections: Sequence[SectionScores],
    *,
    pairs: Sequence[tuple[str, str]],
    conditioned: bool = True,
    alpha: float = 1.0,
) -> list[LouoResult]:
    """For each pair and each held-out unit: fit on the other units, test on it.

    Scores are z-scored **within each section** before pooling, so the held-out test
    asks whether the relation transfers, not whether the section is bright.
    """
    units = sorted({s.unit_id for s in sections})
    z: dict[str, dict[str, dict[str, np.ndarray]]] = {}
    for sec in sections:
        a = dict(sec.a)
        b = dict(sec.b)
        if conditioned:
            a = {f: residualise(a[f], sec.log_lib)[0] for f in a}
            b = {f: residualise(b[f], sec.log_lib)[0] for f in b}
        z[sec.section_id] = {
            "A": {f: zscore(a[f]) for f in a},
            "B": {f: zscore(b[f]) for f in b},
            "unit": sec.unit_id,
        }

    results: list[LouoResult] = []
    for target, predictor in pairs:
        for held in units:
            train = [sid for sid, d in z.items() if d["unit"] != held]
            test = [sid for sid, d in z.items() if d["unit"] == held]
            if len(train) < 3 or not test:
                continue
            xs = np.concatenate([z[s]["A"][predictor] for s in train])
            ys = np.concatenate([z[s]["A"][target] for s in train])
            intercept, slope = _fit_linear(xs, ys, alpha=alpha)
            if not np.isfinite(slope):
                continue
            pred = np.concatenate([intercept + slope * z[s]["B"][predictor] for s in test])
            truth = np.concatenate([z[s]["B"][target] for s in test])
            # oracle upper bound: slope fitted *inside* the held-out unit
            x_test = np.concatenate([z[s]["A"][predictor] for s in test])
            inter_w, slope_w = _fit_linear(x_test, truth, alpha=alpha)
            within = inter_w + slope_w * x_test
            # how well the two gene halves agree on the target itself
            half_a_target = np.concatenate([z[s]["A"][target] for s in test])
            results.append(
                LouoResult(
                    target=target,
                    predictor=predictor,
                    held_unit=held,
                    n_train_sections=len(train),
                    r2_test=_r2(pred, truth),
                    r2_within=_r2(within, truth),
                    r_half_agreement=float(np.corrcoef(half_a_target[ok_half], truth[ok_half])[0, 1])
                    if (ok_half := np.isfinite(half_a_target) & np.isfinite(truth)).sum() > 10
                    else float("nan"),
                    slope=slope,
                )
            )
    return results
