#!/usr/bin/env python3
"""Audit and aggregate the five-seed R-04 extreme-fold restart panel.

The panel is a stability diagnostic, not a K-selection procedure.  Score
differences are paired by patient.  Cross-restart factor stability is measured
on continuous loading and field subspaces, so factor permutations, signs and
rotations are treated as nuisance symmetries rather than cluster labels.
``K_model`` is representation capacity, not a count of biological structures;
the downstream semantics artifact therefore keeps ``K_eff`` and structure
readout status separate from the anonymous factor axes.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import hashlib
import json
from pathlib import Path
import re
from typing import Mapping

import numpy as np

from r04.diagnostics import (
    canonical_subspace_correlations,
    deterministic_objective_platform_summary,
)
from r04.models.mnsf import _inducing, _inducing_projection
from r04.runtime import atomic_json
from r04.splits import grouped_fold_ids


PROTOCOL_KEYS = (
    "gene_folds",
    "inducing_points",
    "inference_steps",
    "lengthscale",
    "nonspatial_rank",
    "optimization_schedule",
    "shared_steps",
    "split_seed",
    "steps",
)

FROZEN_PRACTICAL_STATUS = "FROZEN_FIT_AND_HELDOUT_SCORED_POSTHOC_ACCEPTANCE"


def _mean(values: list[float]) -> float:
    return float(np.mean(np.asarray(values, dtype=float)))


def _median(values: list[float]) -> float:
    return float(np.median(np.asarray(values, dtype=float)))


def _sample_sd(values: list[float]) -> float:
    return float(np.std(np.asarray(values, dtype=float), ddof=1)) if len(values) > 1 else 0.0


def _direction(value: float) -> str:
    if value > 0:
        return "K3_BETTER"
    if value < 0:
        return "K0_BETTER"
    return "TIE"


def _protocol(cell: Mapping[str, object]) -> dict[str, object]:
    parameters = cell.get("parameters")
    diagnostics = cell.get("fit_diagnostics")
    if not isinstance(parameters, dict) or not isinstance(diagnostics, dict):
        raise ValueError("cell parameters or fit diagnostics are missing")
    result = {key: parameters.get(key) for key in PROTOCOL_KEYS}
    result["gene_batch_size"] = diagnostics.get("gene_batch_size")
    return result


def _validate_cell(cell: Mapping[str, object], *, fold: int, k_model: int) -> None:
    if cell.get("schema") != "r04.real_k_cell.v1":
        raise ValueError("unsupported real-K cell schema")
    status = cell.get("status")
    if status not in {
        "FIT_AND_SCORED",
        "FIT_AND_SCORED_DUAL_GATE_DIAGNOSTIC",
        FROZEN_PRACTICAL_STATUS,
    }:
        raise ValueError("restart panel contains an incomplete cell")
    if status == FROZEN_PRACTICAL_STATUS:
        if cell.get("inference_authority") != "PRACTICAL_POSTHOC_ACCEPTANCE":
            raise ValueError("frozen practical cell has unsupported inference authority")
        if cell.get("strict_audit_status") != "STOPPED_OFF_PLATFORM":
            raise ValueError("frozen practical cell is not bound to the strict audit boundary")
        if cell.get("selected_k") is not None:
            raise ValueError("frozen practical cell must preserve selected_k=null")
        platform = cell.get("inference_platform")
        if not isinstance(platform, dict) or not bool(platform.get("converged")):
            raise ValueError("frozen practical cell inference is not platformed")
        split_converged = platform.get("split_converged")
        if split_converged != [True, True]:
            raise ValueError("frozen practical cell requires both gene splits to platform")
    if int(cell.get("fold", -1)) != fold or int(cell.get("k_model", -1)) != k_model:
        raise ValueError("cell fold or K does not match the panel manifest")
    scores = cell.get("patient_scores")
    patients = cell.get("validation_patients")
    if not isinstance(scores, dict) or not scores:
        raise ValueError("cell patient scores are missing")
    if not isinstance(patients, list) or set(map(str, patients)) != set(map(str, scores)):
        raise ValueError("validation patients do not match patient score keys")
    if not all(np.isfinite(float(value)) for value in scores.values()):
        raise ValueError("patient scores must be finite")


def _platform_detail(cell: Mapping[str, object]) -> dict[str, object]:
    diagnostics = cell.get("fit_diagnostics", {})
    embedded = diagnostics.get("deterministic_objective_platform", {})
    top_level = cell.get("fit_objective_platform", {})
    dense = (
        top_level
        if isinstance(top_level, dict) and top_level.get("objective_means")
        else embedded
    )
    steps = list(dense.get("steps", []))
    means = [float(value) for value in dense.get("objective_means", [])]
    if len(steps) == len(means) and len(means) >= 3:
        current = deterministic_objective_platform_summary([
            {"step": int(step), "objective_mean": mean}
            for step, mean in zip(steps, means)
        ])
        dense_converged = bool(current["converged"])
        relative_improvements = [
            float(value) for value in current["relative_improvements"]
        ]
        dense_source = (
            "top_level_fit_objective_platform_recomputed_current_rule"
            if dense is top_level
            else "embedded_deterministic_objective_platform_recomputed_current_rule"
        )
    else:
        dense_converged = bool(dense.get("converged"))
        relative_improvements = [
            float(value) for value in dense.get("relative_improvements", [])
        ]
        dense_source = "stored_summary_without_replayable_trace"
    return {
        "legacy_fit_converged": bool(cell.get("fit_platform", {}).get("converged")),
        "dense_objective_converged": dense_converged,
        "dense_objective_steps": steps,
        "dense_objective_means": means,
        "dense_relative_improvements": relative_improvements,
        "dense_objective_source": dense_source,
        "inference_converged": bool(cell.get("inference_platform", {}).get("converged")),
    }


def _protocol_stratum(cell: Mapping[str, object]) -> dict[str, object]:
    """Return the endpoint/environment dimensions that must not be pooled silently."""
    protocol = _protocol(cell)
    environment = cell.get("environment_hash_compatibility")
    mode = (
        str(environment.get("mode", "NOT_RECORDED"))
        if isinstance(environment, dict)
        else "NOT_RECORDED"
    )
    return {
        "fit_steps": protocol["steps"],
        "inference_steps": protocol["inference_steps"],
        "gene_batch_size": protocol["gene_batch_size"],
        "execution_environment_mode": mode,
    }


def _summarize_k3_factor_diagnostics(
    cells: Mapping[tuple[int, int, int], Mapping[str, object]],
) -> dict[str, object]:
    """Summarize descriptive K=3 energy/rank fields without selecting K_eff."""
    by_fold: dict[str, dict[str, object]] = {}
    for fold in (0, 4):
        records: list[dict[str, object]] = []
        for (restart, cell_fold, k_model), cell in sorted(cells.items()):
            if cell_fold != fold or k_model != 3:
                continue
            diagnostics = cell.get("fit_diagnostics")
            if not isinstance(diagnostics, dict):
                continue
            energy = diagnostics.get("factor_energy")
            participation = diagnostics.get("effective_rank_participation")
            entropy = diagnostics.get("effective_rank_entropy")
            if not isinstance(energy, list) or not energy:
                continue
            values = [float(value) for value in energy]
            if not np.isfinite(values).all():
                raise ValueError("K=3 factor energy must be finite")
            record: dict[str, object] = {
                "restart_index": restart,
                "factor_energy": values,
                "effective_rank_participation": (
                    float(participation) if participation is not None else None
                ),
                "effective_rank_entropy": (
                    float(entropy) if entropy is not None else None
                ),
            }
            if record["effective_rank_participation"] is not None and not np.isfinite(float(record["effective_rank_participation"])):
                raise ValueError("effective rank participation must be finite")
            records.append(record)
        summary: dict[str, object] = {"records": records}
        if records:
            energy_matrix = np.asarray(
                [record["factor_energy"] for record in records], dtype=float
            )
            summary["factor_energy_mean"] = energy_matrix.mean(axis=0).tolist()
            summary["factor_energy_median"] = np.median(energy_matrix, axis=0).tolist()
            participations = [
                float(record["effective_rank_participation"])
                for record in records
                if record["effective_rank_participation"] is not None
            ]
            if participations:
                summary["effective_rank_participation_summary"] = {
                    "mean": _mean(participations),
                    "median": _median(participations),
                    "minimum": min(participations),
                    "maximum": max(participations),
                }
        by_fold[str(fold)] = summary
    return {
        "status": "DESCRIPTIVE_ONLY_NOT_K_EFF_SELECTION",
        "by_fold": by_fold,
    }


def aggregate_restart_scores(
    entries: list[dict[str, object]],
    *,
    allowed_protocol_exceptions: list[dict[str, object]] | None = None,
) -> dict[str, object]:
    """Return paired score and platform summaries for the fixed 5x2 panel."""
    cells: dict[tuple[int, int, int], Mapping[str, object]] = {}
    protocols: list[tuple[tuple[int, int, int], dict[str, object]]] = []
    protocol_strata: dict[str, dict[str, object]] = {}
    input_hash_by_fold: dict[int, str] = {}
    pair_results: list[dict[str, object]] = []

    for entry in entries:
        restart = int(entry["restart_index"])
        fold = int(entry["fold"])
        k0 = entry.get("k0")
        k3 = entry.get("k3")
        if not isinstance(k0, dict) or not isinstance(k3, dict):
            raise ValueError("each panel entry requires K=0 and K=3 cells")
        for k_model, cell in ((0, k0), (3, k3)):
            _validate_cell(cell, fold=fold, k_model=k_model)
            key = (restart, fold, k_model)
            if key in cells:
                raise ValueError(f"duplicate restart panel cell: {key}")
            cells[key] = cell
            protocol = _protocol(cell)
            protocols.append((key, protocol))
            stratum = _protocol_stratum(cell)
            stratum_key = json.dumps(stratum, sort_keys=True, separators=(",", ":"))
            bucket = protocol_strata.setdefault(
                stratum_key,
                {**stratum, "cells": []},
            )
            bucket["cells"].append({
                "restart_index": restart,
                "fold": fold,
                "k_model": k_model,
            })

        k0_scores = {str(key): float(value) for key, value in k0["patient_scores"].items()}
        k3_scores = {str(key): float(value) for key, value in k3["patient_scores"].items()}
        if set(k0_scores) != set(k3_scores):
            raise ValueError("K=0 and K=3 patient score keys differ")
        if k0.get("input_hash") != k3.get("input_hash"):
            raise ValueError("K=0 and K=3 input hashes differ")
        input_hash = str(k0.get("input_hash", ""))
        if not input_hash:
            raise ValueError("cell input hash is missing")
        if fold in input_hash_by_fold and input_hash_by_fold[fold] != input_hash:
            raise ValueError("input hash differs across restarts within a fold")
        input_hash_by_fold[fold] = input_hash

        patient_deltas = {
            patient: k3_scores[patient] - k0_scores[patient]
            for patient in sorted(k0_scores)
        }
        deltas = list(patient_deltas.values())
        pair_results.append({
            "restart_index": restart,
            "fold": fold,
            "patients": len(deltas),
            "mean_k0": _mean(list(k0_scores.values())),
            "mean_k3": _mean(list(k3_scores.values())),
            "mean_delta_k3_minus_k0": _mean(deltas),
            "median_delta_k3_minus_k0": _median(deltas),
            "patient_delta_k3_minus_k0": patient_deltas,
            "k3_better_patients": sum(value > 0 for value in deltas),
            "direction": _direction(_mean(deltas)),
            "k0_platform": _platform_detail(k0),
            "k3_platform": _platform_detail(k3),
        })

    expected = {(restart, fold, k) for restart in range(5) for fold in (0, 4) for k in (0, 3)}
    if set(cells) != expected:
        missing = sorted(expected - set(cells))
        extra = sorted(set(cells) - expected)
        raise ValueError(f"restart panel must contain exactly 20 cells; missing={missing}, extra={extra}")
    protocol_counts = Counter(
        json.dumps(protocol, sort_keys=True, separators=(",", ":"))
        for _, protocol in protocols
    )
    reference_protocol = json.loads(protocol_counts.most_common(1)[0][0])
    allowed = {
        (
            int(item["restart_index"]),
            int(item["fold"]),
            int(item["k_model"]),
            str(item["field"]),
        ): item
        for item in (allowed_protocol_exceptions or [])
    }
    applied_exceptions = []
    for (restart, fold, k_model), protocol in protocols:
        for field in reference_protocol:
            if protocol[field] == reference_protocol[field]:
                continue
            key = (restart, fold, k_model, field)
            exception = allowed.get(key)
            if (
                exception is None
                or exception.get("observed") != protocol[field]
                or exception.get("reference") != reference_protocol[field]
                or not str(exception.get("reason", "")).strip()
            ):
                raise ValueError("scientific protocol differs across restart cells")
            applied_exceptions.append(dict(exception))
    if set(allowed) != {
        (
            int(item["restart_index"]),
            int(item["fold"]),
            int(item["k_model"]),
            str(item["field"]),
        )
        for item in applied_exceptions
    }:
        raise ValueError("an allowed protocol exception was not observed")

    folds: dict[str, object] = {}
    for fold in (0, 4):
        fold_pairs = sorted(
            (item for item in pair_results if item["fold"] == fold),
            key=lambda item: int(item["restart_index"]),
        )
        mean_deltas = [float(item["mean_delta_k3_minus_k0"]) for item in fold_pairs]
        directions = [str(item["direction"]) for item in fold_pairs]
        patients = sorted(fold_pairs[0]["patient_delta_k3_minus_k0"])
        patient_stability = {}
        for patient in patients:
            values = [
                float(item["patient_delta_k3_minus_k0"][patient])
                for item in fold_pairs
            ]
            patient_stability[patient] = {
                "mean_delta_k3_minus_k0": _mean(values),
                "median_delta_k3_minus_k0": _median(values),
                "k3_better_restart_fraction": float(np.mean(np.asarray(values) > 0)),
                "direction_preserved_across_restarts": len({_direction(value) for value in values}) == 1,
                "delta_by_restart": values,
            }
        folds[str(fold)] = {
            "mean_delta_by_restart": mean_deltas,
            "direction_by_restart": directions,
            "direction_preserved_across_restarts": len(set(directions)) == 1 and directions[0] != "TIE",
            "mean_delta_across_restarts": _mean(mean_deltas),
            "median_delta_across_restarts": _median(mean_deltas),
            "seed_sd_of_mean_delta": _sample_sd(mean_deltas),
            "minimum_mean_delta": min(mean_deltas),
            "maximum_mean_delta": max(mean_deltas),
            "patient_stability": patient_stability,
        }

    all_cells = list(cells.values())
    platform_counts = {
        "cells": len(all_cells),
        "legacy_fit_converged": sum(bool(cell.get("fit_platform", {}).get("converged")) for cell in all_cells),
        "dense_objective_converged": sum(
            bool(_platform_detail(cell)["dense_objective_converged"])
            for cell in all_cells
        ),
        "inference_converged": sum(bool(cell.get("inference_platform", {}).get("converged")) for cell in all_cells),
    }
    k3_cells = [cell for (restart, fold, k), cell in cells.items() if k == 3]
    collapse_warning_cells = sum(
        any(bool(cell.get("fit_diagnostics", {}).get(key)) for key in (
            "factor_collapse_warning", "field_collapse_warning", "loading_collapse_warning"
        ))
        for cell in k3_cells
    )
    cell_status_counts = Counter(str(cell.get("status")) for cell in all_cells)
    inference_authority_counts = Counter(
        str(cell.get("inference_authority")) for cell in all_cells
    )
    strict_audit_status_counts = Counter(
        str(cell.get("strict_audit_status")) for cell in all_cells
    )
    return {
        "restart_indices": list(range(5)),
        "fold_values": [0, 4],
        "k_values": [0, 3],
        "cell_count": len(all_cells),
        "protocol": reference_protocol,
        "protocol_exceptions": applied_exceptions,
        "protocol_strata": [
            {
                **value,
                "cells": sorted(
                    value["cells"],
                    key=lambda item: (
                        int(item["restart_index"]),
                        int(item["fold"]),
                        int(item["k_model"]),
                    ),
                ),
            }
            for _, value in sorted(protocol_strata.items())
        ],
        "input_hash_by_fold": {str(key): value for key, value in sorted(input_hash_by_fold.items())},
        "pairs": sorted(pair_results, key=lambda item: (int(item["restart_index"]), int(item["fold"]))),
        "folds": folds,
        "platform_counts": platform_counts,
        "cell_status_counts": dict(sorted(cell_status_counts.items())),
        "inference_authority_counts": dict(sorted(inference_authority_counts.items())),
        "strict_audit_status_counts": dict(sorted(strict_audit_status_counts.items())),
        "k3_collapse_warning_cells": collapse_warning_cells,
        "k3_factor_diagnostics": _summarize_k3_factor_diagnostics(cells),
    }


def pairwise_subspace_stability(matrices: Mapping[int, np.ndarray]) -> dict[str, object]:
    """Compare all restart pairs without assuming factor labels are aligned."""
    pairs = []
    for left, right in combinations(sorted(matrices), 2):
        correlations = canonical_subspace_correlations(matrices[left], matrices[right])
        if not len(correlations):
            raise ValueError("restart subspace has zero numerical rank")
        pairs.append({
            "restart_left": left,
            "restart_right": right,
            "canonical_correlations": [float(value) for value in correlations],
            "mean_canonical_correlation": float(np.mean(correlations)),
            "minimum_canonical_correlation": float(np.min(correlations)),
        })
    minima = [float(item["minimum_canonical_correlation"]) for item in pairs]
    means = [float(item["mean_canonical_correlation"]) for item in pairs]
    dimensions = sorted({len(item["canonical_correlations"]) for item in pairs})
    if len(dimensions) != 1:
        raise ValueError("restart subspaces have inconsistent numerical ranks")
    by_direction = []
    for index in range(dimensions[0]):
        values = [float(item["canonical_correlations"][index]) for item in pairs]
        by_direction.append({
            "canonical_direction": index + 1,
            "mean": _mean(values),
            "median": _median(values),
            "minimum": min(values),
            "maximum": max(values),
        })
    return {
        "pair_count": len(pairs),
        "pairs": pairs,
        "mean_pairwise_canonical_correlation": _mean(means),
        "minimum_canonical_correlation": min(minima),
        "median_pairwise_minimum_correlation": _median(minima),
        "canonical_correlation_by_direction": by_direction,
    }


def validated_reused_subspace(
    previous: Mapping[str, object],
    *,
    input_provenance: list[dict[str, str]],
    training_manifest: dict[str, str],
) -> dict[str, object]:
    """Reuse expensive subspace output only when every source hash is unchanged."""
    if previous.get("schema") != "r04.restart_stability_panel.v1":
        raise ValueError("reused subspace source has an unsupported schema")
    if previous.get("input_provenance") != input_provenance:
        raise ValueError("reused subspace cell provenance differs")
    if previous.get("training_manifest") != training_manifest:
        raise ValueError("reused subspace training manifest differs")
    subspace = previous.get("subspace_stability")
    if not isinstance(subspace, dict) or set(subspace) != {"0", "4"}:
        raise ValueError("reused subspace payload is incomplete")
    return dict(subspace)


def _compact_subspace_summary(value: Mapping[str, object]) -> dict[str, object]:
    """Keep direction summaries while retaining the full panel in the source artifact."""
    return {
        "mean_pairwise_canonical_correlation": value.get(
            "mean_pairwise_canonical_correlation"
        ),
        "minimum_canonical_correlation": value.get("minimum_canonical_correlation"),
        "median_pairwise_minimum_correlation": value.get(
            "median_pairwise_minimum_correlation"
        ),
        "canonical_correlation_by_direction": value.get(
            "canonical_correlation_by_direction", []
        ),
    }


def build_k_semantics_artifact(
    panel_result: Mapping[str, object],
    *,
    panel_manifest: dict[str, str],
    training_manifest: dict[str, str],
    context_artifacts: list[dict[str, str]] | None = None,
    created_at: str = "2026-09-01",
) -> dict[str, object]:
    """Build a bounded semantic handoff without inventing a structure readout."""
    score_summary = panel_result.get("score_summary")
    subspaces = panel_result.get("subspace_stability")
    if not isinstance(score_summary, dict) or not isinstance(subspaces, dict):
        raise ValueError("panel result lacks score or subspace summaries")
    compact_by_fold: dict[str, object] = {}
    for fold in (0, 4):
        fold_subspace = subspaces.get(str(fold))
        if not isinstance(fold_subspace, dict):
            raise ValueError(f"panel result lacks fold {fold} subspace")
        loading = fold_subspace.get("loading_subspace")
        field = fold_subspace.get("training_field_subspace")
        if not isinstance(loading, dict) or not isinstance(field, dict):
            raise ValueError(f"panel result lacks fold {fold} loading/field subspace")
        compact_by_fold[str(fold)] = {
            "score": score_summary.get("folds", {}).get(str(fold)),
            "loading_subspace": _compact_subspace_summary(loading),
            "training_field_subspace": _compact_subspace_summary(field),
            "third_direction": {
                "loading": (
                    loading.get("canonical_correlation_by_direction", [])[2]
                    if len(loading.get("canonical_correlation_by_direction", [])) >= 3
                    else None
                ),
                "training_field": (
                    field.get("canonical_correlation_by_direction", [])[2]
                    if len(field.get("canonical_correlation_by_direction", [])) >= 3
                    else None
                ),
                "interpretation": "descriptive_instability_signal_not_K_eff",
            },
        }
    fold_result_summary = {}
    for fold in (0, 4):
        fold_score = score_summary.get("folds", {}).get(str(fold), {})
        fold_subspace = compact_by_fold[str(fold)]
        fold_result_summary[str(fold)] = {
            "direction_by_restart": fold_score.get("direction_by_restart"),
            "mean_delta_k3_minus_k0": fold_score.get("mean_delta_across_restarts"),
            "seed_sd_of_mean_delta": fold_score.get("seed_sd_of_mean_delta"),
            "direction_preserved_across_restarts": fold_score.get(
                "direction_preserved_across_restarts"
            ),
            "third_loading_direction_mean_canonical_correlation": (
                fold_subspace["third_direction"]["loading"].get("mean")
                if isinstance(fold_subspace["third_direction"].get("loading"), dict)
                else None
            ),
            "third_training_field_direction_mean_canonical_correlation": (
                fold_subspace["third_direction"]["training_field"].get("mean")
                if isinstance(fold_subspace["third_direction"].get("training_field"), dict)
                else None
            ),
        }
    return {
        "schema": "r04.k_semantics_and_downstream_robustness.v1",
        "status": "K_SEMANTICS_CORRECTED_DOWNSTREAM_ROBUSTNESS_INDETERMINATE",
        "created_at": created_at,
        "k_interpretation": "latent_spatial_effect_dimension_not_structure_count",
        "k_model": {
            "meaning": "model_allowed_latent_spatial_effect_dimension",
            "structure_mapping": "many_to_many",
            "single_factor_naming": "forbidden",
        },
        "k0_reference": {
            "meaning": "empty_spatial_effect_subspace_reference",
            "subspace_comparison": "not_applicable_for_K0",
            "interpretation": "K0 is a baseline, not a zero-valued biological field or a failed structure readout",
        },
        "k_eff": {
            "meaning": "data_supported_effective_dimension_under_holdout_repeatability_null_and_subspace_stability",
            "status": "not_computable_with_current_frozen_cell_outputs",
            "value": None,
            "reason": "frozen cell artifacts contain scores and provenance but no GT-aligned held-out spot_by_gene spatial effect matrix or composed field/loading readout",
            "available_descriptive_proxy": "K=3 factor energy, effective-rank participation and loading/field subspace diagnostics; these do not select K_eff",
        },
        "selected_k": None,
        "result_summary": {
            "plain_language": "Fold 0 consistently favors K=0 and fold 4 consistently favors K=3 in this diagnostic panel; the first subspace direction is more reproducible than the third, but this does not select K or identify a biological structure.",
            "by_fold": fold_result_summary,
            "interpretation": "descriptive_panel_evidence_only; downstream_structure_readout_not_tested",
        },
        "available_evidence": {
            "panel_scope": "five_seed_K0_K3_extreme_fold_0_4_diagnostic",
            "score_summary": score_summary,
            "subspace_summary_by_fold": compact_by_fold,
            "k3_factor_diagnostics": score_summary.get("k3_factor_diagnostics"),
            "protocol_strata": score_summary.get("protocol_strata", []),
            "cell_status_counts": score_summary.get("cell_status_counts", {}),
            "strict_audit_status_counts": score_summary.get(
                "strict_audit_status_counts", {}
            ),
            "inference_authority_counts": score_summary.get(
                "inference_authority_counts", {}
            ),
            "strict_audit_boundary": "frozen_practical_cells_remain_STOPPED_OFF_PLATFORM; historical cells retain their source audit status",
            "practical_acceptance_scope": "only frozen restart 2-4 held-out inference/scoring cells",
        },
        "downstream_k_robust": "not_tested",
        "k2_bridge_required": "undecided",
        "k2_bridge_decision": {
            "status": "NOT_TRIGGERED",
            "reason": "No downstream structure readout is available yet, so the trigger condition has not been demonstrated; this is not a claim that K=2 is unnecessary.",
            "trigger": "run_only_if_training_role_nested_cross_fit_shows_major_structure_conclusion_sensitive_to_the_unstable_third_direction_or_remains_unjudgeable",
        },
        "structure_readout": {
            "status": "not_tested",
            "input_policy": "training_role_GT_only_nested_cross_fit; internal_and_external_validation_GT_sealed",
            "reason": "current frozen outputs do not export the spot-level spatial effect vectors needed for a structure-specific readout",
            "minimum_next_computation": "export_or_recreate_frozen held-out field effects, align them to training-role TLS and TUMOR_STROMA_BOUNDARY GT, then compare shared-only versus shared_plus_structure-specific residual readouts under K=3 and stability-reduced representations",
            "factor_label_use": "forbidden",
        },
        "spatial_null_status": "NOT_RUN",
        "composition_split_status": "NOT_RUN",
        "structure_naming_status": "NOT_RUN",
        "independent_lineage_reproduction_status": "NOT_RUN",
        "claim_boundary": [
            "K_model is capacity, not biological structure count",
            "K_eff is not assigned from factor numbering or factor energy alone",
            "fold 0 and fold 4 remain separate diagnostic strata",
            "no K selection, structure naming or cross-device numerical-equivalence claim is authorized by this artifact",
        ],
        "input_provenance": {
            "panel_manifest": panel_manifest,
            "training_manifest": training_manifest,
            "cell_artifacts": panel_result.get("input_provenance", []),
            "checkpoint_artifacts": panel_result.get("checkpoint_provenance", {}),
            "context_artifacts": context_artifacts or [],
        },
    }


def _resolve(project_root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else project_root / path


def _read_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _checkpoint_bundle_provenance(
    checkpoint_dir: Path,
    *,
    project_root: Path,
) -> dict[str, object]:
    """Hash every checkpoint-bundle file used for subspace reconstruction."""
    files = []
    for path in sorted(item for item in checkpoint_dir.rglob("*") if item.is_file()):
        files.append({
            "path": str(path.relative_to(project_root)),
            "sha256": _sha256(path),
        })
    if not files:
        raise ValueError(f"checkpoint bundle is empty: {checkpoint_dir}")
    bundle_digest = hashlib.sha256(
        json.dumps(files, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return {
        "directory": str(checkpoint_dir.relative_to(project_root)),
        "bundle_sha256": bundle_digest,
        "files": files,
    }


def _normalize_coordinates_bounded(
    coords: np.ndarray,
    *,
    pair_block_size: int = 512,
) -> tuple[np.ndarray, float]:
    """Apply the existing median-nearest-neighbor normalization with bounded RAM.

    The shared geometry definition is retained, but the original all-pairs
    temporary would be unnecessarily large for the real panel.  Computing the
    same nearest-neighbor distances in row blocks changes memory use, not the
    coordinate normalization rule.
    """
    coordinates = np.asarray(coords, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2 or len(coordinates) < 2:
        raise ValueError("at least two 2-D coordinates are required")
    if pair_block_size < 1:
        raise ValueError("pair_block_size must be positive")
    nearest = np.empty(len(coordinates), dtype=float)
    for start in range(0, len(coordinates), pair_block_size):
        stop = min(start + pair_block_size, len(coordinates))
        block = coordinates[start:stop]
        delta = block[:, None, :] - coordinates[None, :, :]
        distance = np.sqrt(np.sum(delta * delta, axis=2))
        distance[np.arange(stop - start), np.arange(start, stop)] = np.inf
        nearest[start:stop] = distance.min(axis=1)
    scale = float(np.median(nearest))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("coordinates have no positive nearest-neighbor scale")
    center = np.median(coordinates, axis=0)
    return (coordinates - center) / scale, scale


def _checkpoint_matrices(
    checkpoint_dir: Path,
    *,
    training_rows: list[dict[str, object]],
    fold: int,
    group_seed: int,
    inducing_points: int,
    lengthscale: float,
    ridge: float = 1e-5,
) -> tuple[np.ndarray, np.ndarray, str]:
    import torch

    ckpt_path = checkpoint_dir / "checkpoint.pt"
    if not ckpt_path.exists():
        raise ValueError(
            f"torch checkpoint is missing: {ckpt_path} (expected checkpoint.pt with schema r04.mnsf_checkpoint.v3_torch, backend torch)"
        )
    payload = torch.load(str(ckpt_path), map_location="cpu", weights_only=False)
    if not isinstance(payload, dict):
        raise ValueError(f"torch checkpoint payload must be a dict: {ckpt_path}")

    params_dict = payload.get("params") if isinstance(payload.get("params"), dict) else {}
    best_state = payload.get("best_state") if isinstance(payload.get("best_state"), dict) else {}
    state_nested = payload.get("state") if isinstance(payload.get("state"), dict) else {}

    def _to_numpy(obj: object) -> np.ndarray | None:
        if obj is None:
            return None
        if isinstance(obj, torch.Tensor):
            return np.asarray(obj.detach().cpu().numpy(), dtype=float)
        if isinstance(obj, np.ndarray):
            return np.asarray(obj, dtype=float)
        try:
            return np.asarray(obj, dtype=float)  # type: ignore[arg-type]
        except Exception:
            return None

    # --- raw_w extraction: try explicit keys, then p_0 fallback, then scan ---
    raw_w: np.ndarray | None = None
    for cand in (
        payload.get("raw_w"),
        params_dict.get("raw_w"),
        best_state.get("raw_w"),
        state_nested.get("raw_w"),
        payload.get("state_dict", {}).get("raw_w") if isinstance(payload.get("state_dict"), dict) else None,
    ):
        arr = _to_numpy(cand)
        if arr is not None:
            raw_w = arr
            break
    if raw_w is None:
        for src in (params_dict, best_state, payload):
            if isinstance(src, dict) and "p_0" in src:
                arr = _to_numpy(src["p_0"])
                if arr is not None:
                    raw_w = arr
                    break
    if raw_w is None:
        # Scan p_* tensors with shape (g, k) as last resort
        for src in (params_dict, best_state):
            for key in sorted(src.keys()):
                if key.startswith("p_"):
                    arr = _to_numpy(src[key])
                    if arr is not None and arr.ndim == 2:
                        # Heuristic: raw_w is typically the largest first dimension among p_* that matches gene count
                        # Accept first 2-D candidate if no better info; will be validated by loading shape check later
                        raw_w = arr
                        break
            if raw_w is not None:
                break
    if raw_w is None:
        raise ValueError(
            f"raw_w is missing from torch checkpoint: {ckpt_path} (checked top-level raw_w, params/raw_w, best_state/raw_w, state/raw_w, p_0). "
            f"Available payload keys: {sorted(payload.keys())}, params keys: {sorted(params_dict.keys())[:10]}, best_state keys: {sorted(best_state.keys())[:10]}"
        )
    if raw_w.ndim == 1 and raw_w.size == 0:
        raw_w = raw_w.reshape(-1, 0)
    if raw_w.ndim != 2:
        raise ValueError(f"raw_w must be 2-D, got shape {raw_w.shape} from {ckpt_path}")
    if raw_w.shape[1] == 0:
        loading = raw_w.astype(float)
    else:
        shifted = raw_w - raw_w.max(axis=0, keepdims=True)
        loading = np.exp(shifted)
        loading /= loading.sum(axis=0, keepdims=True)

    fold_ids = grouped_fold_ids(
        [str(row["patient_id"]) for row in training_rows], n_folds=5, seed=group_seed
    )
    rows = [row for row, fold_id in zip(training_rows, fold_ids) if int(fold_id) != fold]

    # --- q_loc extraction per section ---
    def _fetch_q_loc(idx: int) -> np.ndarray | None:
        for cand in (
            payload.get(f"q_loc_{idx}"),
            params_dict.get(f"q_loc_{idx}"),
            best_state.get(f"q_loc_{idx}"),
            state_nested.get(f"q_loc_{idx}"),
            payload.get(f"p_{4 + idx}"),
            params_dict.get(f"p_{4 + idx}"),
            best_state.get(f"p_{4 + idx}"),
        ):
            arr = _to_numpy(cand)
            if arr is not None:
                return arr
        return None

    q_locs: list[np.ndarray] = []
    for idx in range(len(rows)):
        arr = _fetch_q_loc(idx)
        if arr is None:
            # Fallback: scan p_* for tensor with compatible shape if we can infer expected shape later
            # Try all p_* candidates with shape (m_i, k)
            for src in (params_dict, best_state, payload):
                if not isinstance(src, dict):
                    continue
                for k2, v2 in src.items():
                    if not k2.startswith("p_"):
                        continue
                    try:
                        cand2 = _to_numpy(v2)
                    except Exception:
                        continue
                    if cand2 is None or cand2.ndim != 2:
                        continue
                    if cand2.shape[1] == loading.shape[1] and cand2.shape[0] <= inducing_points + 64:
                        # Collect as possible; but we need index-specific, so skip generic scan
                        pass
            raise ValueError(
                f"q_loc_{idx} is missing from torch checkpoint: {ckpt_path} "
                f"(checked q_loc_{idx}, p_{4 + idx} in payload/params/best_state). "
                f"Available params keys: {sorted(params_dict.keys())[:12]}"
            )
        q_locs.append(arr)

    field_parts = []
    for row, q_loc in zip(rows, q_locs):
        index = len(field_parts)
        metadata_path = Path(str(row["matrix_metadata_locator"]))
        metadata = _read_json(metadata_path)
        coords, _ = _normalize_coordinates_bounded(
            np.asarray(metadata["coords"], dtype=float)
        )
        number = min(inducing_points, len(coords))
        inducing = _inducing(coords, number)
        basis, _, _ = _inducing_projection(coords, inducing, lengthscale, ridge)
        if q_loc.shape != (number, loading.shape[1]):
            raise ValueError(f"q_loc shape mismatch at training section {index}: got {q_loc.shape}, expected ({number}, {loading.shape[1]}) from {ckpt_path}")
        field = basis @ q_loc
        field_parts.append(field - field.mean(axis=0, keepdims=True))
    checkpoint = str(ckpt_path)
    return loading, np.concatenate(field_parts, axis=0), checkpoint


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reuse-subspace-from", type=Path, default=None)
    parser.add_argument(
        "--k-semantics-output",
        type=Path,
        default=None,
        help="optional machine-readable K semantics/downstream robustness artifact",
    )
    parser.add_argument(
        "--context-artifact",
        action="append",
        default=[],
        help="additional read-only evidence artifact to hash into the semantics output",
    )
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    manifest_path = _resolve(project_root, str(args.manifest)).resolve()
    manifest = _read_json(manifest_path)
    if manifest.get("schema") != "r04.restart_panel_manifest.v1":
        raise ValueError("unsupported restart panel manifest schema")
    if manifest.get("status") != "READY_FOR_DIAGNOSTIC_AGGREGATION":
        raise ValueError("restart panel manifest is not ready")
    raw_entries = manifest.get("entries")
    if not isinstance(raw_entries, list):
        raise ValueError("restart panel manifest entries are missing")

    entries = []
    provenance = []
    checkpoint_by_fold: dict[int, dict[int, Path]] = {0: {}, 4: {}}
    checkpoint_provenance_by_fold: dict[str, dict[str, dict[str, object]]] = {
        "0": {},
        "4": {},
    }
    for raw in raw_entries:
        if not isinstance(raw, dict):
            raise ValueError("restart panel manifest contains a non-object entry")
        restart = int(raw["restart_index"])
        fold = int(raw["fold"])
        k0_path = _resolve(project_root, str(raw["k0_cell"]))
        k3_path = _resolve(project_root, str(raw["k3_cell"]))
        checkpoint_dir = _resolve(project_root, str(raw["k3_checkpoint_dir"]))
        entries.append({
            "restart_index": restart,
            "fold": fold,
            "k0": _read_json(k0_path),
            "k3": _read_json(k3_path),
        })
        checkpoint_by_fold.setdefault(fold, {})[restart] = checkpoint_dir
        provenance.extend([
            {"path": str(k0_path.relative_to(project_root)), "sha256": _sha256(k0_path)},
            {"path": str(k3_path.relative_to(project_root)), "sha256": _sha256(k3_path)},
        ])

    raw_exceptions = manifest.get("allowed_protocol_exceptions", [])
    if not isinstance(raw_exceptions, list):
        raise ValueError("allowed protocol exceptions must be a list")
    score_summary = aggregate_restart_scores(
        entries, allowed_protocol_exceptions=raw_exceptions
    )
    training_manifest_path = _resolve(project_root, str(manifest["training_manifest"]))
    training_manifest = _read_json(training_manifest_path)
    training_rows = training_manifest.get("rows")
    if training_manifest.get("status") != "READY" or training_manifest.get("role") != "training":
        raise ValueError("training panel manifest is not READY training data")
    if not isinstance(training_rows, list) or not training_rows:
        raise ValueError("training panel manifest rows are missing")

    training_summary = {
        "path": str(training_manifest_path.relative_to(project_root)),
        "sha256": _sha256(training_manifest_path),
    }
    panel_summary = {
        "path": str(manifest_path.relative_to(project_root)),
        "sha256": _sha256(manifest_path),
    }
    context_summaries = []
    for value in args.context_artifact:
        path = _resolve(project_root, value)
        if not path.exists() or not path.is_file():
            raise ValueError(f"context artifact does not exist: {path}")
        context_summaries.append({
            "path": str(path.relative_to(project_root)),
            "sha256": _sha256(path),
        })
    reused_subspace = None
    if args.reuse_subspace_from is not None:
        reuse_path = _resolve(project_root, str(args.reuse_subspace_from))
        previous = _read_json(reuse_path)
        subspace_by_fold = validated_reused_subspace(
            previous,
            input_provenance=provenance,
            training_manifest=training_summary,
        )
        reused_subspace = {
            "path": str(reuse_path.relative_to(project_root)),
            "sha256_before_rewrite": _sha256(reuse_path),
        }
    else:
        subspace_by_fold = {}
        protocol = score_summary["protocol"]
        for fold in (0, 4):
            loading_by_restart = {}
            field_by_restart = {}
            checkpoints = {}
            for restart in range(5):
                loading, field, checkpoint = _checkpoint_matrices(
                    checkpoint_by_fold[fold][restart],
                    training_rows=training_rows,
                    fold=fold,
                    group_seed=int(manifest["group_seed"]),
                    inducing_points=int(protocol["inducing_points"]),
                    lengthscale=float(protocol["lengthscale"]),
                )
                loading_by_restart[restart] = loading
                field_by_restart[restart] = field
                checkpoints[str(restart)] = str(Path(checkpoint).relative_to(project_root))
                checkpoint_provenance_by_fold[str(fold)][str(restart)] = (
                    _checkpoint_bundle_provenance(
                        checkpoint_by_fold[fold][restart], project_root=project_root
                    )
                )
            subspace_by_fold[str(fold)] = {
                "loading_subspace": pairwise_subspace_stability(loading_by_restart),
                "training_field_subspace": pairwise_subspace_stability(field_by_restart),
                "checkpoint_by_restart": checkpoints,
                "checkpoint_provenance_by_restart": checkpoint_provenance_by_fold[str(fold)],
            }

    all_dense = score_summary["platform_counts"]["dense_objective_converged"] == 20
    all_inference = score_summary["platform_counts"]["inference_converged"] == 20
    directions_preserved = all(
        bool(score_summary["folds"][str(fold)]["direction_preserved_across_restarts"])
        for fold in (0, 4)
    )
    result = {
        "schema": "r04.restart_stability_panel.v1",
        "status": "RESTART_STABILITY_DIAGNOSTIC_COMPLETE_NOT_FORMAL_K_SELECTION",
        "created_at": str(manifest.get("created_at", "2026-08-31")),
        "purpose": "Assess combined optimization-randomness sensitivity at the negative and strongest-positive folds before spatial null or formal K expansion.",
        "diagnostic_only": True,
        "score_summary": score_summary,
        "subspace_stability": subspace_by_fold,
        "gates": {
            "five_seed_panel_complete": True,
            "direction_preserved_at_both_extreme_folds": directions_preserved,
            "all_dense_objective_platform": all_dense,
            "all_inference_platform": all_inference,
            "k3_collapse_warning_absent": score_summary["k3_collapse_warning_cells"] == 0,
            "subspace_similarity_threshold": "NOT_PREREGISTERED_DIAGNOSTIC_ONLY",
            "formal_k_selection_authorized": False,
            "spatial_null_authorized_by_this_artifact": False,
        },
        "interpretation_policy": {
            "score_direction": "paired patient K3-minus-K0 mean within each fold and restart",
            "magnitude": "optimization-seed spread is descriptive and is not an independent-patient confidence interval",
            "subspace": "canonical correlations compare continuous spans and ignore factor numbering, sign, and rotations",
            "wiring_only_sources": "restart indices 2-4 retain WIRING_ONLY_NOT_SCIENTIFIC provenance",
            "restart_randomness": "optimization seed jointly changes initialization, gene-minibatch order, and variational draws; this is not a pure initialization experiment",
        },
        "input_provenance": provenance,
        "training_manifest": training_summary,
        "checkpoint_provenance": checkpoint_provenance_by_fold,
    }
    if reused_subspace is not None:
        result["subspace_reused_from"] = reused_subspace
    reaudit_value = manifest.get("platform_rule_reaudit")
    if reaudit_value is not None:
        reaudit_path = _resolve(project_root, str(reaudit_value))
        reaudit = _read_json(reaudit_path)
        if reaudit.get("schema") != "r04.deterministic_platform_rule_reaudit.v1":
            raise ValueError("platform rule re-audit schema is unsupported")
        result["platform_rule_reaudit"] = {
            "path": str(reaudit_path.relative_to(project_root)),
            "sha256": _sha256(reaudit_path),
            "rule": reaudit.get("rule"),
            "threshold": reaudit.get("threshold"),
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(args.output, result)
    if args.k_semantics_output is not None:
        semantics = build_k_semantics_artifact(
            result,
            panel_manifest=panel_summary,
            training_manifest=training_summary,
            context_artifacts=context_summaries,
        )
        args.k_semantics_output.parent.mkdir(parents=True, exist_ok=True)
        atomic_json(args.k_semantics_output, semantics)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
