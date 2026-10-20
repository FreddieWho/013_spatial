"""FA-04: molecular internal validation on the confirmed real sections.

Two halves of every field definition are scored *independently* and with the same
target-independent library:

    spec_A: fields = A-half genes, witness = B-half union
    spec_B: fields = B-half genes, witness = A-half union

Both therefore share one external denominator that excludes every analysed field
gene, so neither half can enter the other through the library.

Three questions are answered per section, and one is deliberately NOT answered:

  1. map agreement   -- do the two halves give the same spatial picture of a field?
  2. relation agreement -- do cross-field same-point correlations agree between halves?
  3. scale agreement -- do the two halves imply the same correlation range?
  4. (NOT tested) displacement -- FA-03B showed it is not identifiable from this
     observation, so no displacement claim is made on real data.

No human TLS labels, no patient-level claims: the 8 known units are within-study
units, not independent patients.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np

from . import scales
from .pipeline import InputSpec, build_inputs


def half_specs(fields: Mapping[str, Sequence[str]], splits: Mapping[str, Any]) -> tuple[InputSpec, InputSpec]:
    a_defs = {fid: list(splits[fid].a_genes) for fid in fields}
    b_defs = {fid: list(splits[fid].b_genes) for fid in fields}
    a_union = frozenset(g for v in a_defs.values() for g in v)
    b_union = frozenset(g for v in b_defs.values() for g in v)
    spec_a = InputSpec(fields=a_defs, witness_genes=b_union, k_neighbors=6, min_library_genes=100)
    spec_b = InputSpec(fields=b_defs, witness_genes=a_union, k_neighbors=6, min_library_genes=100)
    return spec_a, spec_b


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 10:
        return float("nan")
    if np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return float("nan")
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def same_point_relations(scores: Mapping[str, np.ndarray]) -> dict[tuple[str, str], float]:
    ids = sorted(scores)
    out: dict[tuple[str, str], float] = {}
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            out[(ids[i], ids[j])] = pearson(scores[ids[i]], scores[ids[j]])
    return out


def field_range(
    score: np.ndarray,
    xy: np.ndarray,
    valid: np.ndarray,
    *,
    max_lag_pitch: float = 12.0,
    n_bins: int = 12,
) -> dict[str, Any]:
    from . import geometry

    support = np.asarray(valid, dtype=bool) & np.isfinite(score)
    if support.sum() < 50:
        return {"status": scales.NOT_MEASURABLE, "range": None, "n_usable_bins": 0}
    pitch = geometry.estimate_pitch(xy[support])
    pairs = scales.pairs_within(xy, support, max_lag_pitch * pitch)
    edges = scales.lag_edges(pitch, max_lag_pitch, n_bins)
    gamma, counts, centers = scales.empirical_variogram(score, pairs, edges)
    fit = scales.fit_variogram(centers, gamma, counts)
    fit["range_pitch"] = (fit["range"] / pitch) if fit.get("range") else None
    fit["pitch"] = pitch
    return fit


@dataclass
class SectionValidation:
    section_id: str
    unit_id: str
    n_obs: int
    n_visible: int
    map_agreement: dict[str, float]
    relation_a: dict[tuple[str, str], float]
    relation_b: dict[tuple[str, str], float]
    range_a: dict[str, dict[str, Any]]
    range_b: dict[str, dict[str, Any]]
    pitch: float
    depth_corr: dict[str, float]


def validate_section(
    table: Any,
    fields: Mapping[str, Sequence[str]],
    splits: Mapping[str, Any],
    *,
    unit_id: str,
    section_id: str,
    max_lag_pitch: float = 12.0,
) -> SectionValidation:
    xy = np.asarray(table.xy, dtype=float)
    support = np.ones(table.n_obs, dtype=bool)
    if "in_tissue" in table.aux:
        support = np.asarray(table.aux["in_tissue"]).astype(bool)
    spec_a, spec_b = half_specs(fields, splits)
    inputs_a = build_inputs(table.counts, table.genes, xy, support, support, spec_a)
    inputs_b = build_inputs(table.counts, table.genes, xy, support, support, spec_b)

    scores_a = {fid: inputs_a.arrays[f"field::{fid}"] for fid in spec_a.fields}
    scores_b = {fid: inputs_b.arrays[f"field::{fid}"] for fid in spec_b.fields}
    agreement = {fid: pearson(scores_a[fid], scores_b[fid]) for fid in scores_a}
    rel_a = same_point_relations(scores_a)
    rel_b = same_point_relations(scores_b)
    range_a = {fid: field_range(scores_a[fid], xy, support, max_lag_pitch=max_lag_pitch) for fid in scores_a}
    range_b = {fid: field_range(scores_b[fid], xy, support, max_lag_pitch=max_lag_pitch) for fid in scores_b}
    log_lib = np.log1p(inputs_a.arrays["library"])
    depth_corr = {fid: pearson(scores_a[fid], log_lib) for fid in scores_a}
    return SectionValidation(
        section_id=section_id,
        unit_id=unit_id,
        n_obs=int(table.n_obs),
        n_visible=int(support.sum()),
        map_agreement=agreement,
        relation_a=rel_a,
        relation_b=rel_b,
        range_a=range_a,
        range_b=range_b,
        pitch=float(inputs_a.meta["pitch"]),
        depth_corr=depth_corr,
    )


def summarise(results: Sequence[SectionValidation]) -> dict[str, Any]:
    """Unit-level summary: one number per unit, not per spot."""
    field_ids = sorted(results[0].map_agreement) if results else []
    pair_ids = sorted(results[0].relation_a) if results else []

    def per_unit(values_by_section: dict[str, float], unit_of: dict[str, str]) -> dict[str, float]:
        buckets: dict[str, list[float]] = {}
        for section, value in values_by_section.items():
            if np.isfinite(value):
                buckets.setdefault(unit_of[section], []).append(value)
        return {unit: float(np.mean(vals)) for unit, vals in buckets.items()}

    unit_of = {r.section_id: r.unit_id for r in results}
    map_agreement = {
        fid: per_unit({r.section_id: r.map_agreement[fid] for r in results}, unit_of) for fid in field_ids
    }
    relation_gap = {}
    relation_sign = {}
    for pair in pair_ids:
        gap = {}
        sign = {}
        for r in results:
            a, b = r.relation_a[pair], r.relation_b[pair]
            if np.isfinite(a) and np.isfinite(b):
                gap[r.section_id] = abs(a - b)
                sign[r.section_id] = float(np.sign(a) == np.sign(b))
        relation_gap["~".join(pair)] = per_unit(gap, unit_of)
        relation_sign["~".join(pair)] = per_unit(sign, unit_of)

    def range_gap(which: str) -> dict[str, dict[str, float]]:
        out: dict[str, dict[str, float]] = {}
        for fid in field_ids:
            gaps: dict[str, float] = {}
            for r in results:
                ra = r.range_a[fid]
                rb = r.range_b[fid]
                if ra.get("range_pitch") and rb.get("range_pitch"):
                    gaps[r.section_id] = abs(ra["range_pitch"] - rb["range_pitch"])
            out[fid] = per_unit(gaps, unit_of)
        return out

    def median_of(nested: dict[str, dict[str, float]]) -> dict[str, float]:
        return {
            key: float(np.median(list(vals.values()))) if vals else float("nan")
            for key, vals in nested.items()
        }

    depth_corr = {
        fid: per_unit({r.section_id: r.depth_corr[fid] for r in results}, unit_of) for fid in field_ids
    }
    ok_a = {fid: sum(1 for r in results if r.range_a[fid]["status"] == scales.OK) for fid in field_ids}
    ok_b = {fid: sum(1 for r in results if r.range_b[fid]["status"] == scales.OK) for fid in field_ids}
    return {
        "n_sections": len(results),
        "n_units": len({r.unit_id for r in results}),
        "map_agreement_median_by_unit": median_of(map_agreement),
        "map_agreement_median_overall": float(
            np.median([v for fid in field_ids for v in map_agreement[fid].values()])
        ),
        "relation_gap_median": median_of(relation_gap),
        "relation_sign_agreement_median": median_of(relation_sign),
        "range_gap_pitch_median": median_of(range_gap("both")),
        "depth_corr_median": median_of(depth_corr),
        "n_sections_with_ok_range_A": ok_a,
        "n_sections_with_ok_range_B": ok_b,
        "notes": [
            "unit-level numbers aggregate sections within a within-study unit, then take the median over units",
            "the two halves share one target-independent library; neither half enters the other's denominator",
            "depth_corr is the Pearson correlation of each field score with log library size; a high value means the field partly tracks sequencing depth",
            "no displacement claim: FA-03B showed it is not identifiable from this observation",
        ],
    }
