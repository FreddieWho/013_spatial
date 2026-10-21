"""FA-07A: depth conditioning on the real sections.

Question: after removing the part of each field score that is explained by
sequencing depth, do the cross-field relations and the A/B agreement survive?

Method (all on visible spots only):
  * depth proxy d = log1p(library), where the library excludes every analysed field
    gene (same target-independent denominator used for the scores);
  * residualise each score on [1, d, d^2] by least squares (quadratic keeps mild
    non-linearity without letting the fit chase structure);
  * recompute map agreement (A vs B halves) and cross-field same-point correlations
    on the residuals, and report before/after plus the depth R^2.

This is a descriptive conditioning step, not a causal correction: removing a depth
trend can also remove real biology that happens to track depth.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np

from . import real_validation as rv


def depth_design(log_lib: np.ndarray, degree: int = 2) -> np.ndarray:
    x = np.asarray(log_lib, dtype=float)
    ok = np.isfinite(x)
    mu, sd = float(np.mean(x[ok])), float(np.std(x[ok]) or 1.0)
    z = (x - mu) / sd
    cols = [np.ones_like(z)]
    for k in range(1, degree + 1):
        cols.append(z ** k)
    return np.column_stack(cols)


def residualise(values: np.ndarray, log_lib: np.ndarray, degree: int = 2) -> tuple[np.ndarray, float]:
    """Least-squares residual of ``values`` on a polynomial in log library."""
    y = np.asarray(values, dtype=float)
    ok = np.isfinite(y) & np.isfinite(log_lib)
    out = np.full_like(y, np.nan)
    if ok.sum() < 10:
        return out, float("nan")
    X = depth_design(np.asarray(log_lib, dtype=float), degree)[ok]
    beta, *_ = np.linalg.lstsq(X, y[ok], rcond=None)
    pred = X @ beta
    ss_res = float(np.sum((y[ok] - pred) ** 2))
    ss_tot = float(np.sum((y[ok] - np.mean(y[ok])) ** 2))
    out[ok] = y[ok] - pred
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return out, r2


@dataclass
class ConditionedSection:
    section_id: str
    unit_id: str
    map_before: dict[str, float]
    map_after: dict[str, float]
    relation_before_a: dict[tuple[str, str], float]
    relation_after_a: dict[tuple[str, str], float]
    relation_before_b: dict[tuple[str, str], float]
    relation_after_b: dict[tuple[str, str], float]
    depth_r2_a: dict[str, float]
    depth_r2_b: dict[str, float]


def condition_section(
    table: Any,
    fields: Mapping[str, Sequence[str]],
    splits: Mapping[str, Any],
    *,
    unit_id: str,
    section_id: str,
    degree: int = 2,
) -> ConditionedSection:
    xy = np.asarray(table.xy, dtype=float)
    support = np.ones(table.n_obs, dtype=bool)
    if "in_tissue" in table.aux:
        support = np.asarray(table.aux["in_tissue"]).astype(bool)
    spec_a, spec_b = rv.half_specs(fields, splits)
    from .pipeline import build_inputs

    inp_a = build_inputs(table.counts, table.genes, xy, support, support, spec_a)
    inp_b = build_inputs(table.counts, table.genes, xy, support, support, spec_b)
    log_lib = np.log1p(inp_a.arrays["library"])
    log_lib_b = np.log1p(inp_b.arrays["library"])

    scores_a = {fid: inp_a.arrays[f"field::{fid}"] for fid in spec_a.fields}
    scores_b = {fid: inp_b.arrays[f"field::{fid}"] for fid in spec_b.fields}

    map_before = {fid: rv.pearson(scores_a[fid], scores_b[fid]) for fid in scores_a}

    res_a: dict[str, np.ndarray] = {}
    res_b: dict[str, np.ndarray] = {}
    r2_a: dict[str, float] = {}
    r2_b: dict[str, float] = {}
    for fid in scores_a:
        res_a[fid], r2_a[fid] = residualise(scores_a[fid], log_lib, degree)
        res_b[fid], r2_b[fid] = residualise(scores_b[fid], log_lib_b, degree)
    map_after = {fid: rv.pearson(res_a[fid], res_b[fid]) for fid in res_a}

    return ConditionedSection(
        section_id=section_id,
        unit_id=unit_id,
        map_before=map_before,
        map_after=map_after,
        relation_before_a=rv.same_point_relations(scores_a),
        relation_after_a=rv.same_point_relations(res_a),
        relation_before_b=rv.same_point_relations(scores_b),
        relation_after_b=rv.same_point_relations(res_b),
        depth_r2_a=r2_a,
        depth_r2_b=r2_b,
    )


def summarise(results: Sequence[ConditionedSection]) -> dict[str, Any]:
    field_ids = sorted(results[0].map_before) if results else []
    pairs = sorted(results[0].relation_before_a) if results else []

    def unit_median(values_by_section: Mapping[str, float]) -> float:
        buckets: dict[str, list[float]] = {}
        for sec, val in values_by_section.items():
            if np.isfinite(val):
                buckets.setdefault(next(r.unit_id for r in results if r.section_id == sec), []).append(float(val))
        per = [float(np.median(v)) for v in buckets.values()]
        return float(np.median(per)) if per else float("nan")

    def pair_stat(pair: tuple[str, str], kind: str) -> float:
        vals: dict[str, list[float]] = {}
        for r in results:
            source = {
                "before_a": r.relation_before_a,
                "after_a": r.relation_after_a,
                "before_b": r.relation_before_b,
                "after_b": r.relation_after_b,
            }[kind]
            v = source[pair]
            if np.isfinite(v):
                vals.setdefault(r.unit_id, []).append(float(v))
        per = [float(np.median(v)) for v in vals.values()]
        return float(np.median(per)) if per else float("nan")

    def sign_consistent(pair: tuple[str, str]) -> int:
        """Units where the conditioned half-A and half-B signs agree."""
        units: dict[str, list[tuple[float, float]]] = {}
        for r in results:
            a, b = r.relation_after_a[pair], r.relation_after_b[pair]
            if np.isfinite(a) and np.isfinite(b):
                units.setdefault(r.unit_id, []).append((a, b))
        n = 0
        for vals in units.values():
            a = float(np.median([v[0] for v in vals]))
            b = float(np.median([v[1] for v in vals]))
            if np.sign(a) == np.sign(b):
                n += 1
        return n

    return {
        "n_sections": len(results),
        "n_units": len({r.unit_id for r in results}),
        "depth_r2_median_A": {
            fid: unit_median({r.section_id: r.depth_r2_a[fid] for r in results}) for fid in field_ids
        },
        "map_agreement_before": {
            fid: unit_median({r.section_id: r.map_before[fid] for r in results}) for fid in field_ids
        },
        "map_agreement_after": {
            fid: unit_median({r.section_id: r.map_after[fid] for r in results}) for fid in field_ids
        },
        "relation_by_pair": {
            "~".join(p): {
                "before_A": pair_stat(p, "before_a"),
                "after_A": pair_stat(p, "after_a"),
                "before_B": pair_stat(p, "before_b"),
                "after_B": pair_stat(p, "after_b"),
                "sign_consistent_units_after": sign_consistent(p),
            }
            for p in pairs
        },
        "notes": [
            "depth proxy is log1p(library) with the target-independent denominator; residual = score - quadratic(depth)",
            "conditioning is descriptive: it can also remove biology that tracks depth",
            "unit medians: per-unit median first, then median over units",
        ],
    }
