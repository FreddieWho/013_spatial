#!/usr/bin/env python3
"""Assemble the current R-04 K/readout/resource evidence handoff.

This does not recompute scores or alter historical artifacts.  It validates
the upstream schemas, copies their machine-readable summaries, and records
the hashes and the current conditional K=2 bridge status in one small handoff.
Validation GT remains represented only by the training-role readout status;
internal/external validation GT is never loaded here.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Mapping

from r04.runtime import atomic_json


def _read_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source(path: Path, *, project_root: Path, value: Mapping[str, object]) -> dict[str, object]:
    return {
        "path": str(path.resolve().relative_to(project_root.resolve())),
        "sha256": _sha256(path),
        "schema": value.get("schema"),
        "status": value.get("status"),
    }


def _inner_summary(value: object) -> dict[str, object]:
    if not isinstance(value, dict):
        return {"status": "NOT_AVAILABLE"}
    return {
        "status": value.get("status"),
        "valid_inner_folds": value.get("valid_inner_folds"),
        "total_inner_folds": value.get("total_inner_folds"),
        "specific_minus_shared_mean_mse_delta": value.get(
            "specific_minus_shared_mean_mse_delta"
        ),
        "mean_auc_delta_by_target": value.get("mean_auc_delta_by_target"),
    }


def build_final_artifact(
    aggregate: Mapping[str, object],
    structure_readout: Mapping[str, object],
    k2_stage: Mapping[str, object],
    k2_manifest: Mapping[str, object],
    *,
    source_artifacts: list[dict[str, object]],
    created_at: str = "2026-09-01",
) -> dict[str, object]:
    """Build the current evidence handoff without inventing a K result."""
    if aggregate.get("schema") != "r04.k_semantics_and_downstream_robustness.v1":
        raise ValueError("unsupported aggregate semantics schema")
    if structure_readout.get("schema") != "r04.structure_readout_panel.v1":
        raise ValueError("unsupported structure readout schema")
    if k2_stage.get("schema") != "r04.k2_bridge_stage_a.v1":
        raise ValueError("unsupported K=2 stage schema")
    if k2_manifest.get("schema") != "r04.k2_bridge_stage_a_manifest.v1":
        raise ValueError("unsupported K=2 manifest schema")
    if k2_manifest.get("status") != "FROZEN_PRECOMPUTE_MANIFEST":
        raise ValueError("K=2 manifest is not frozen")
    if k2_stage.get("status") != "BLOCKED_CPU_RUNTIME":
        raise ValueError("this handoff expects the current CPU runtime block")

    evidence = aggregate.get("available_evidence")
    if not isinstance(evidence, dict):
        raise ValueError("aggregate lacks available evidence")
    score_summary = evidence.get("score_summary")
    subspace_summary = evidence.get("subspace_summary_by_fold")
    entries = structure_readout.get("entries")
    if not isinstance(score_summary, dict) or not isinstance(subspace_summary, dict):
        raise ValueError("aggregate lacks score/subspace summaries")
    if not isinstance(entries, list) or len(entries) != 10:
        raise ValueError("structure readout must contain ten K=3 training entries")

    readout_summary = []
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("structure readout entry is not an object")
        coverage = entry.get("coverage_training")
        paired = coverage.get("paired") if isinstance(coverage, dict) else None
        paired_patients = paired.get("patients", []) if isinstance(paired, dict) else []
        readout_summary.append({
            "restart_index": entry.get("restart_index"),
            "fold": entry.get("fold"),
            "paired_patient_count": len(paired_patients),
            "paired_patients": list(paired_patients),
            "full_k3": _inner_summary(
                entry.get("full_k3", {}).get("inner_crossfit")
                if isinstance(entry.get("full_k3"), dict)
                else None
            ),
            "stable_rank2_sensitivity": _inner_summary(
                entry.get("stable_rank2_sensitivity", {}).get("inner_crossfit")
                if isinstance(entry.get("stable_rank2_sensitivity"), dict)
                else None
            ),
        })

    source = copy.deepcopy(source_artifacts)
    return {
        "schema": "r04.k_semantics_and_downstream_robustness.v2",
        "status": "K_SEMANTICS_CORRECTED_K2_BRIDGE_BLOCKED_CPU_RUNTIME",
        "created_at": created_at,
        "source_artifacts": source,
        "k_interpretation": "latent_spatial_effect_dimension_not_structure_count",
        "k_model": copy.deepcopy(aggregate.get("k_model")),
        "k_eff": copy.deepcopy(aggregate.get("k_eff")),
        "selected_k": None,
        "score_summary": copy.deepcopy(score_summary),
        "subspace_stability_by_fold": copy.deepcopy(subspace_summary),
        "third_direction_interpretation": (
            "descriptive instability signal; not a biological structure count and "
            "not sufficient to assign K_eff"
        ),
        "downstream_k_robust": structure_readout.get("downstream_k_robust"),
        "k2_bridge_required": structure_readout.get("k2_bridge_required"),
        "structure_readout": {
            "status": structure_readout.get("status"),
            "scope": structure_readout.get("scope"),
            "downstream_k_robust": structure_readout.get("downstream_k_robust"),
            "k2_bridge_required": structure_readout.get("k2_bridge_required"),
            "reason": structure_readout.get("k2_bridge_reason"),
            "training_role_entry_count": len(entries),
            "per_entry_summary": readout_summary,
            "internal_external_validation_gt": structure_readout.get(
                "internal_external_validation_gt"
            ),
            "factor_naming": "forbidden",
        },
        "k2_bridge": {
            "required": structure_readout.get("k2_bridge_required"),
            "manifest_status": k2_manifest.get("status"),
            "manifest_entry_count": len(k2_manifest.get("entries", [])),
            "stage_status": k2_stage.get("status"),
            "blocked_reason": k2_stage.get("blocked_reason"),
            "attempted_cell": k2_stage.get("attempted_cell"),
            "scientific_result": "NOT_AVAILABLE_NO_K2_CELL_COMPLETED",
            "gpu_required_for_bounded_completion": k2_stage.get(
                "gpu_required_for_bounded_completion"
            ),
        },
        "spatial_null_status": "NOT_RUN",
        "composition_split_status": "NOT_RUN",
        "structure_naming_status": "FORBIDDEN_NOT_RUN",
        "independent_lineage_reproduction_status": "NOT_RUN",
        "claim_boundary": [
            "K_model is representation capacity, not a count of biological structures",
            "K_eff and selected_k remain unresolved",
            "the readout is training-role descriptive evidence only and does not establish structure identity",
            "the K=2 bridge has no scientific result because its first CPU cell stopped before a checkpoint",
            "spatial nulls, composition controls and independent lineage replication remain pending",
            "internal and external validation GT remains sealed",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aggregate", type=Path, required=True)
    parser.add_argument("--structure-readout", type=Path, required=True)
    parser.add_argument("--k2-stage", type=Path, required=True)
    parser.add_argument("--k2-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    project_root = Path(".").resolve()
    paths = [
        args.aggregate.resolve(),
        args.structure_readout.resolve(),
        args.k2_stage.resolve(),
        args.k2_manifest.resolve(),
    ]
    values = [_read_json(path) for path in paths]
    result = build_final_artifact(
        values[0], values[1], values[2], values[3],
        source_artifacts=[
            _source(path, project_root=project_root, value=value)
            for path, value in zip(paths, values)
        ],
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(args.output, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
