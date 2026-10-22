"""FA-12: cross-field neighbour versus same-field neighbour.

The question is whether a neighbour's *other* field predicts a hidden spot better
than that spot's own field at its neighbours.  Same-spot coupling is reported, but
it is not spatial.

Decision, fixed before the real run
------------------------------------
* Mask 10% of spots in every section (seed 0). Neighbour features use only the
  other 90%. A linear map is fit on the other units and tested on the held-out
  unit's masked spots.
* A pair supports a spatial relation only if at least 6 of 8 units have
  cross-neighbour R² above same-neighbour R² and above both rotation nulls
  (137° and 90°).
* The all-spot association is a sensitivity check. One protocol passing is not
  an upgrade.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import numpy as np
from scipy.spatial import cKDTree

from . import geometry
from .depth_conditioning import residualise
from .real_relations import SectionScores, _fit_linear, _r2, zscore

DECISION_UNITS = 6
NULL_ANGLES = (137.0, 90.0)
PAIRS = (
    ("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
    ("B_AXIS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
    ("B_AXIS", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
    ("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
)


def misalign_values(xy: np.ndarray, values: np.ndarray, angle_deg: float) -> np.ndarray:
    """Move the field by a rotation about its centroid, then sample back at ``xy``.

    Finite values are resampled from the rotated positions. Missing spots stay
    missing. This keeps the value multiset's smoothness but breaks alignment
    with another field.
    """
    xy = np.asarray(xy, dtype=float)
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan)
    ok = np.isfinite(values)
    if ok.sum() < 3:
        return out
    centre = xy[ok].mean(axis=0)
    ang = np.deg2rad(angle_deg)
    rot = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
    moved = (xy[ok] - centre) @ rot.T + centre
    _, idx = cKDTree(moved).query(xy[ok], k=1)
    out[ok] = values[ok][idx]
    return out


def holdout_mask(n: int, fraction: float, seed: int, section_id: str) -> np.ndarray:
    """Deterministic per-section mask. Fraction 0 hides nothing."""
    if fraction <= 0:
        return np.zeros(n, dtype=bool)
    digest = 0
    for ch in section_id:
        digest = (digest * 131 + ord(ch)) % 1_000_003
    rng = np.random.default_rng(seed + digest)
    mask = rng.random(n) < fraction
    # keep both sides non-empty on any section large enough to score
    if n >= 20 and (mask.sum() < 5 or (~mask).sum() < 5):
        cut = max(5, int(round(n * fraction)))
        mask = np.zeros(n, dtype=bool)
        mask[:cut] = True
        rng.shuffle(mask)
    return mask


def _prepare(section: SectionScores, conditioned: bool) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    a = {f: np.asarray(v, dtype=float) for f, v in section.a.items()}
    b = {f: np.asarray(v, dtype=float) for f, v in section.b.items()}
    if conditioned:
        a = {f: residualise(v, section.log_lib)[0] for f, v in a.items()}
        b = {f: residualise(v, section.log_lib)[0] for f, v in b.items()}
    visible = np.asarray(section.visible, dtype=bool)
    for store in (a, b):
        for f, v in store.items():
            v = zscore(v)
            v[~visible] = np.nan
            store[f] = v
    return a, b


def _collect(sections: Sequence[Mapping[str, Any]], ids: Sequence[str], key: str, selector: str) -> tuple[np.ndarray, np.ndarray]:
    xs, ys = [], []
    for sec in sections:
        if sec["section_id"] not in ids:
            continue
        sel = sec["holdout"] if selector == "holdout" else ~sec["holdout"]
        sel = sel & np.isfinite(sec["features"][key]) & np.isfinite(sec["y"])
        if sel.any():
            xs.append(sec["features"][key][sel])
            ys.append(sec["y"][sel])
    if not xs:
        return np.array([]), np.array([])
    return np.concatenate(xs), np.concatenate(ys)


def _r2_transferred(train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray, test_y: np.ndarray, alpha: float) -> float:
    intercept, slope = _fit_linear(train_x, train_y, alpha=alpha)
    if not np.isfinite(slope):
        return float("nan")
    return _r2(intercept + slope * test_x, test_y)


def compare_pairs(
    sections: Sequence[SectionScores],
    pairs: Sequence[tuple[str, str]] = PAIRS,
    *,
    conditioned: bool = True,
    holdout_fraction: float = 0.1,
    seed: int = 0,
    k: int = 6,
    angles: Sequence[float] = NULL_ANGLES,
    alpha: float = 1.0,
    protocol: str = "masked",
) -> list[dict[str, Any]]:
    """Leave-one-unit-out comparison of three predictors and the rotation nulls."""
    units = sorted({s.unit_id for s in sections})
    prepared: list[dict[str, Any]] = []
    for sec in sections:
        a, b = _prepare(sec, conditioned)
        graph = geometry.build_adjacency(sec.xy, sec.visible, k=k).graph
        holdout = holdout_mask(len(sec.xy), holdout_fraction, seed, sec.section_id)
        holdout &= np.asarray(sec.visible, dtype=bool)
        prepared.append({
            "section_id": sec.section_id,
            "unit_id": sec.unit_id,
            "xy": np.asarray(sec.xy, dtype=float),
            "visible": np.asarray(sec.visible, dtype=bool),
            "graph": graph,
            "holdout": holdout,
            "a": a,
            "b": b,
        })

    rows: list[dict[str, Any]] = []
    for target, predictor in pairs:
        built = []
        for sec in prepared:
            observed = sec["visible"] & ~sec["holdout"]
            y = sec["b"][target]
            cross_src = sec["a"][predictor]
            features = {
                "same_neighbour": geometry.neighbor_sum(sec["graph"], y, observed & np.isfinite(y)),
                "cross_neighbour": geometry.neighbor_sum(sec["graph"], cross_src, observed & np.isfinite(cross_src)),
                "same_spot": cross_src,
            }
            for angle in angles:
                mis = misalign_values(sec["xy"], cross_src, angle)
                features[f"rot_{int(angle)}"] = geometry.neighbor_sum(
                    sec["graph"], mis, observed & np.isfinite(mis)
                )
            built.append({**sec, "y": y, "features": features})

        by_unit: dict[str, list[str]] = {}
        for sec in built:
            by_unit.setdefault(sec["unit_id"], []).append(sec["section_id"])
        selector = "holdout" if holdout_fraction > 0 else "train"
        # fraction 0: holdout mask is empty, so both train (~holdout) and a fake
        # holdout selector would be wrong. Use all finite visible spots for both.
        if holdout_fraction <= 0:
            for sec in built:
                sec["holdout"] = sec["visible"] & np.isfinite(sec["y"])
            selector_test, selector_train = "holdout", "holdout"
        else:
            selector_test, selector_train = "holdout", "train"

        for held in units:
            test_ids = by_unit.get(held, [])
            train_ids = [sid for u, sids in by_unit.items() if u != held for sid in sids]
            if len({sec["unit_id"] for sec in built if sec["section_id"] in train_ids}) < 3:
                continue
            row: dict[str, Any] = {
                "protocol": protocol,
                "target": target,
                "predictor": predictor,
                "held_unit": held,
                "n_test": int(_collect(built, test_ids, "same_spot", selector_test)[1].size),
            }
            for key in ("same_neighbour", "cross_neighbour", "same_spot", *(f"rot_{int(a)}" for a in angles)):
                tr_x, tr_y = _collect(built, train_ids, key, selector_train)
                te_x, te_y = _collect(built, test_ids, key, selector_test)
                row[f"r2_{key}"] = _r2_transferred(tr_x, tr_y, te_x, te_y, alpha) if te_y.size else float("nan")
                row["n_test"] = int(te_y.size)
            cross = row["r2_cross_neighbour"]
            nulls = [row[f"r2_rot_{int(a)}"] for a in angles]
            row["beats_same"] = bool(np.isfinite(cross) and cross > row["r2_same_neighbour"])
            row["beats_null"] = bool(np.isfinite(cross) and all(np.isfinite(v) and cross > v for v in nulls))
            row["beats_same_spot"] = bool(np.isfinite(cross) and cross > row["r2_same_spot"])
            row["passes"] = bool(row["beats_same"] and row["beats_null"])
            rows.append(row)
    return rows


def summarise(rows: Sequence[Mapping[str, Any]], *, min_units: int = DECISION_UNITS) -> dict[str, Any]:
    """Unit-level decision. A protocol upgrades only at ``min_units`` passes."""
    protocols = sorted({r["protocol"] for r in rows})
    pairs = sorted({(r["target"], r["predictor"]) for r in rows})
    out: dict[str, Any] = {"by_protocol": {}, "upgrade_af01": False, "passing_pairs": []}
    passing = []
    for protocol in protocols:
        entry: dict[str, Any] = {}
        for target, predictor in pairs:
            sel = [r for r in rows if r["protocol"] == protocol and r["target"] == target and r["predictor"] == predictor]
            n_pass = int(sum(1 for r in sel if r["passes"]))
            n = len(sel)
            upgraded = n_pass >= min_units and n >= min_units
            entry[f"{target}<-{predictor}"] = {
                "n_units": n,
                "n_pass": n_pass,
                "n_beat_same": int(sum(1 for r in sel if r["beats_same"])),
                "n_beat_null": int(sum(1 for r in sel if r["beats_null"])),
                "n_beat_same_spot": int(sum(1 for r in sel if r["beats_same_spot"])),
                "median_r2_cross": _median(sel, "r2_cross_neighbour"),
                "median_r2_same": _median(sel, "r2_same_neighbour"),
                "median_r2_same_spot": _median(sel, "r2_same_spot"),
                "median_r2_rot_137": _median(sel, "r2_rot_137"),
                "median_r2_rot_90": _median(sel, "r2_rot_90"),
                "upgrades": upgraded,
            }
            if upgraded:
                passing.append({"protocol": protocol, "target": target, "predictor": predictor})
        out["by_protocol"][protocol] = entry
    # both protocols must pass the same pair; a mask-only win is not an upgrade
    masked = { (p["target"], p["predictor"]) for p in passing if p["protocol"] == "masked" }
    full = { (p["target"], p["predictor"]) for p in passing if p["protocol"] == "full" }
    agreed = masked & full
    out["passing_pairs"] = [{"target": t, "predictor": p} for t, p in sorted(agreed)]
    out["upgrade_af01"] = bool(agreed)
    out["decision_rule"] = (
        f"upgrade only if the same pair passes on >={min_units} units in both the "
        "masked protocol and the full-spot sensitivity protocol"
    )
    return out


def _median(rows: Sequence[Mapping[str, Any]], key: str) -> float:
    vals = [float(r[key]) for r in rows if key in r and np.isfinite(r[key])]
    return float(np.median(vals)) if vals else float("nan")
