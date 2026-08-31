"""Candidate matching and fail-closed R-04 gates."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, replace
from typing import Iterable

import numpy as np
from scipy.optimize import linear_sum_assignment

from .spatial import spatial_lengthscale_is_interior
from .types import CandidateField, FieldFit, GateResult


def cosine_similarity(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float).ravel()
    right = np.asarray(right, dtype=float).ravel()
    denominator = np.linalg.norm(left) * np.linalg.norm(right)
    if denominator == 0:
        return 0.0
    return float(np.dot(left, right) / denominator)


def signed_correlation(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float).ravel()
    right = np.asarray(right, dtype=float).ravel()
    if len(left) != len(right) or np.std(left) == 0 or np.std(right) == 0:
        return 0.0
    return float(np.corrcoef(left, right)[0, 1])


def field_fit_hash(fit: FieldFit) -> str:
    payload = {
        "model_id": fit.model_id,
        "factor_id": fit.factor_id,
        "section_id": fit.section_id,
        "input_hash": fit.input_hash,
        "loading": np.asarray(fit.loading, dtype=np.float64).tolist(),
        "field_mean": np.asarray(fit.field_mean, dtype=np.float64).tolist(),
        "field_sd": np.asarray(fit.field_sd, dtype=np.float64).tolist(),
        "lengthscale": fit.lengthscale,
        "spatial_variance": fit.spatial_variance,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _group_fits(fits: Iterable[FieldFit]) -> dict[str, list[FieldFit]]:
    groups: dict[str, list[FieldFit]] = {}
    for fit in fits:
        groups.setdefault(fit.factor_id, []).append(fit)
    for values in groups.values():
        values.sort(key=lambda fit: fit.section_id)
    return groups


def _group_input_hash(group: list[FieldFit]) -> str:
    hashes = {fit.input_hash for fit in group}
    return next(iter(hashes)) if len(hashes) == 1 else ""


def _group_stats(left: list[FieldFit], right: list[FieldFit]) -> tuple[float, float, tuple[str, ...]] | None:
    if not left or not right:
        return None
    left_hashes = {fit.input_hash for fit in left}
    right_hashes = {fit.input_hash for fit in right}
    if len(left_hashes) != 1 or left_hashes != right_hashes:
        return None
    left_by_section = {fit.section_id: fit for fit in left}
    right_by_section = {fit.section_id: fit for fit in right}
    common = tuple(sorted(set(left_by_section) & set(right_by_section)))
    if not common:
        return None
    loading = abs(cosine_similarity(left[0].loading, right[0].loading))
    field_correlations = [
        abs(signed_correlation(left_by_section[section].field_mean, right_by_section[section].field_mean))
        for section in common
    ]
    return loading, float(np.mean(field_correlations)), common


def match_model_factors(
    mnsf: Iterable[FieldFit],
    signed: Iterable[FieldFit],
    *,
    loading_threshold: float = 0.60,
    field_threshold: float = 0.50,
) -> list[CandidateField]:
    """Build a tiered factor union across sections; never assigns a spot to a class."""
    left_groups = _group_fits(mnsf)
    right_groups = _group_fits(signed)
    used: set[str] = set()
    candidates: list[CandidateField] = []
    for factor_id in sorted(left_groups):
        item_group = left_groups[factor_id]
        best: tuple[float, str, float, float, tuple[str, ...]] | None = None
        for other_id in sorted(right_groups):
            if other_id in used:
                continue
            stats = _group_stats(item_group, right_groups[other_id])
            if stats is None:
                continue
            loading, field, common = stats
            score = min(loading, field)
            if best is None or score > best[0]:
                best = (score, other_id, loading, field, common)
        input_hash = _group_input_hash(item_group)
        if best and best[2] >= loading_threshold and best[3] >= field_threshold:
            other_id = best[1]
            other_group = right_groups[other_id]
            used.add(other_id)
            fit_hashes = [field_fit_hash(fit) for fit in [*item_group, *other_group]]
            candidates.append(CandidateField(
                candidate_id=f"CAND_{len(candidates)+1:04d}", tier="BOTH",
                member_factor_ids=(factor_id, other_id),
                loading_cosine=best[2], field_correlation=best[3],
                restart_fraction=0.0, patient_fraction=0.0,
                raw_status="PROVISIONAL_UNAGGREGATED", input_hash=input_hash,
                field_hash=hashlib.sha256("".join(sorted(fit_hashes)).encode("ascii")).hexdigest(),
                diagnostics={"matched_sections": list(best[4]), "matching_scope": "factor_across_sections"},
            ))
        else:
            fit_hashes = sorted(field_fit_hash(fit) for fit in item_group)
            candidates.append(CandidateField(
                candidate_id=f"CAND_{len(candidates)+1:04d}", tier="MNSF_ONLY",
                member_factor_ids=(factor_id,), loading_cosine=0.0,
                field_correlation=0.0, restart_fraction=0.0, patient_fraction=0.0,
                raw_status="PROVISIONAL_UNAGGREGATED", input_hash=input_hash,
                field_hash=hashlib.sha256("".join(fit_hashes).encode("ascii")).hexdigest(),
                diagnostics={"n_sections": len(item_group), "matching_scope": "factor_across_sections"},
            ))
    for factor_id in sorted(right_groups):
        if factor_id not in used:
            item_group = right_groups[factor_id]
            fit_hashes = sorted(field_fit_hash(fit) for fit in item_group)
            candidates.append(CandidateField(
                candidate_id=f"CAND_{len(candidates)+1:04d}", tier="SIGNED_ONLY",
                member_factor_ids=(factor_id,), loading_cosine=0.0,
                field_correlation=0.0, restart_fraction=0.0, patient_fraction=0.0,
                raw_status="PROVISIONAL_UNAGGREGATED", input_hash=_group_input_hash(item_group),
                field_hash=hashlib.sha256("".join(fit_hashes).encode("ascii")).hexdigest(),
                diagnostics={"n_sections": len(item_group), "matching_scope": "factor_across_sections"},
            ))
    return candidates


def _patient_field_similarity(left: list[FieldFit], right: list[FieldFit]) -> tuple[float, float, tuple[str, ...]] | None:
    left_hashes = {fit.input_hash for fit in left}
    right_hashes = {fit.input_hash for fit in right}
    if len(left_hashes) != 1 or len(right_hashes) != 1 or left_hashes != right_hashes:
        return None
    left_by_section = {fit.section_id: fit for fit in left}
    right_by_section = {fit.section_id: fit for fit in right}
    common = tuple(sorted(set(left_by_section) & set(right_by_section)))
    if not common:
        return None
    by_patient: dict[str, list[float]] = {}
    signed_by_patient: dict[str, list[float]] = {}
    for section_id in common:
        left_fit = left_by_section[section_id]
        right_fit = right_by_section[section_id]
        patient = str(left_fit.diagnostics.get("patient_id", section_id))
        signed = signed_correlation(left_fit.field_mean, right_fit.field_mean)
        by_patient.setdefault(patient, []).append(abs(signed))
        signed_by_patient.setdefault(patient, []).append(signed)
    patient_values = [float(np.mean(values)) for values in by_patient.values()]
    signed_values = [float(np.mean(values)) for values in signed_by_patient.values()]
    return float(np.median(patient_values)), float(np.median(signed_values)), common


def _one_to_one_match_model_factors(
    mnsf: Iterable[FieldFit],
    signed: Iterable[FieldFit],
    *,
    loading_threshold: float,
    field_threshold: float,
) -> list[CandidateField]:
    """Match factor-level continuous signatures with a conjunctive Hungarian gate."""
    left_groups = _group_fits(mnsf)
    right_groups = _group_fits(signed)
    left_ids = sorted(left_groups)
    right_ids = sorted(right_groups)
    stats: dict[tuple[str, str], dict[str, object]] = {}
    score = np.full((len(left_ids), len(right_ids)), -1.0, dtype=float)
    for i, left_id in enumerate(left_ids):
        for j, right_id in enumerate(right_ids):
            left_group, right_group = left_groups[left_id], right_groups[right_id]
            field_stats = _patient_field_similarity(left_group, right_group)
            if field_stats is None:
                continue
            gene_signed = cosine_similarity(left_group[0].loading, right_group[0].loading)
            gene = abs(gene_signed)
            field, signed_field, common = field_stats
            pair_score = min(gene, field)
            stats[(left_id, right_id)] = {
                "gene": gene,
                "gene_signed": gene_signed,
                "field": field,
                "signed_field": signed_field,
                "common": common,
                "score": pair_score,
                "sign_flip": -1 if gene_signed < 0 else 1,
            }
            if gene >= loading_threshold and field >= field_threshold:
                score[i, j] = pair_score
    assignments: dict[str, tuple[str, dict[str, object], float]] = {}
    if score.size:
        rows, columns = linear_sum_assignment(-score)
        for row, column in zip(rows, columns):
            left_id, right_id = left_ids[row], right_ids[column]
            pair = stats.get((left_id, right_id))
            if pair is not None and score[row, column] >= 0:
                alternatives = sorted(
                    (value["score"] for key, value in stats.items() if key[0] == left_id and key[1] != right_id),
                    reverse=True,
                )
                margin = float(pair["score"] - alternatives[0]) if alternatives else float(pair["score"])
                assignments[left_id] = (right_id, pair, margin)
    candidates: list[CandidateField] = []
    used_right: set[str] = set()

    def append_single(group_id: str, group: list[FieldFit], tier: str) -> None:
        hashes = sorted(field_fit_hash(fit) for fit in group)
        candidates.append(CandidateField(
            candidate_id=f"CAND_{len(candidates) + 1:04d}",
            tier=tier,
            member_factor_ids=(group_id,),
            loading_cosine=0.0,
            field_correlation=0.0,
            restart_fraction=0.0,
            patient_fraction=0.0,
            raw_status="PROVISIONAL_UNAGGREGATED",
            input_hash=_group_input_hash(group),
            field_hash=hashlib.sha256("".join(hashes).encode("ascii")).hexdigest(),
            diagnostics={"n_sections": len(group), "matching_scope": "factor_one_to_one_hungarian"},
        ))

    for left_id in left_ids:
        assignment = assignments.get(left_id)
        if assignment is None:
            append_single(left_id, left_groups[left_id], "MNSF_ONLY")
            continue
        right_id, pair, margin = assignment
        used_right.add(right_id)
        members = [*left_groups[left_id], *right_groups[right_id]]
        hashes = sorted(field_fit_hash(fit) for fit in members)
        candidates.append(CandidateField(
            candidate_id=f"CAND_{len(candidates) + 1:04d}",
            tier="BOTH",
            member_factor_ids=(left_id, right_id),
            loading_cosine=float(pair["gene"]),
            field_correlation=float(pair["field"]),
            restart_fraction=0.0,
            patient_fraction=0.0,
            raw_status="PROVISIONAL_UNAGGREGATED",
            input_hash=_group_input_hash(left_groups[left_id]),
            field_hash=hashlib.sha256("".join(hashes).encode("ascii")).hexdigest(),
            diagnostics={
                "matched_sections": list(pair["common"]),
                "matching_scope": "factor_one_to_one_hungarian",
                "patient_level_field_correlation": pair["field"],
                "assignment_margin": margin,
                "sign_flip": pair["sign_flip"],
                "signed_loading_cosine": pair["gene_signed"],
                "signed_field_correlation_after_loading_flip": float(pair["signed_field"] * pair["sign_flip"]),
                "loading_threshold": loading_threshold,
                "field_threshold": field_threshold,
            },
        ))
    for right_id in right_ids:
        if right_id not in used_right:
            append_single(right_id, right_groups[right_id], "SIGNED_ONLY")
    return candidates


def match_model_factors(
    mnsf: Iterable[FieldFit],
    signed: Iterable[FieldFit],
    *,
    loading_threshold: float = 0.60,
    field_threshold: float = 0.50,
) -> list[CandidateField]:
    return _one_to_one_match_model_factors(
        mnsf,
        signed,
        loading_threshold=loading_threshold,
        field_threshold=field_threshold,
    )


def candidate_registry_hash(candidates: Iterable[CandidateField]) -> str:
    payload = [asdict(candidate) for candidate in candidates]
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def aggregate_candidate_stability(
    candidate: CandidateField,
    *,
    restart_total: int,
    restart_hits: int,
    patient_total: int,
    patient_hits: int,
) -> CandidateField:
    """Attach observed stability; never infer 100% from one provisional fit."""
    if restart_total <= 0 or patient_total <= 0:
        raise ValueError("stability denominators must be positive")
    if not (0 <= restart_hits <= restart_total and 0 <= patient_hits <= patient_total):
        raise ValueError("stability counts are outside their denominators")
    return replace(
        candidate,
        restart_fraction=restart_hits / restart_total,
        patient_fraction=patient_hits / patient_total,
        raw_status="AGGREGATED_OBSERVED",
        diagnostics={
            **candidate.diagnostics,
            "restart_hits": restart_hits,
            "restart_total": restart_total,
            "patient_hits": patient_hits,
            "patient_total": patient_total,
        },
    )


def evaluate_candidate_gate(
    candidate: CandidateField,
    *,
    heldout_delta: float,
    heldout_ci_low: float,
    spatial_variance_probability: float,
    null_fdr: float,
    lengthscale: float,
    max_patient_weight: float,
    independent_lineages: int,
    groups_per_lineage: int,
    single_model: bool = False,
    uncertainty_complete: bool = False,
    lengthscale_supported: bool = False,
) -> GateResult:
    reasons: list[str] = []
    if candidate.restart_fraction < 0.8:
        reasons.append("restart_stability_below_4_of_5")
    if candidate.patient_fraction < 0.7:
        reasons.append("patient_stability_below_70_percent")
    if heldout_delta <= 0 or heldout_ci_low <= 0:
        reasons.append("no_positive_grouped_heldout_increment")
    if spatial_variance_probability <= 0.95:
        reasons.append("posterior_spatial_variance_probability_not_above_0.95")
    if null_fdr > (0.05 if single_model else 0.10):
        reasons.append("spatial_null_fdr_above_threshold")
    if not spatial_lengthscale_is_interior(lengthscale):
        reasons.append("lengthscale_at_boundary")
    if max_patient_weight > 0.50:
        reasons.append("single_patient_dominates")
    if not uncertainty_complete:
        reasons.append("uncertainty_does_not_cover_loading_dispersion_and_lengthscale")
    if not candidate.input_hash:
        reasons.append("candidate_input_hash_missing")
    if not lengthscale_supported:
        reasons.append("lengthscale_not_data_supported")
    required_lineages = 3 if single_model else 2
    if independent_lineages < required_lineages or groups_per_lineage < 3:
        reasons.append("insufficient_independent_lineages")
    status = "PASS_R04_STABLE_FIELD" if not reasons else "SCIENTIFIC_FAIL_CANDIDATE_GATE"
    return GateResult(status=status, reasons=tuple(reasons), metrics={
        "heldout_delta": heldout_delta,
        "heldout_ci_low": heldout_ci_low,
        "spatial_variance_probability": spatial_variance_probability,
        "null_fdr": null_fdr,
        "lengthscale": lengthscale,
        "max_patient_weight": max_patient_weight,
        "uncertainty_complete": float(uncertainty_complete),
        "lengthscale_supported": float(lengthscale_supported),
    })
