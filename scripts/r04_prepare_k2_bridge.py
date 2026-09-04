#!/usr/bin/env python3
"""Freeze the minimal K=2 bridge manifest before any bridge training.

The bridge changes only ``K_model``.  It reuses the K=3 panel's exact
optimization seeds, outer folds, fit lengths and source fit diagnostics so
the comparison is a capacity bridge rather than a new tuning campaign.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from r04.runtime import atomic_json


RESTARTS = (0, 1, 2)
FOLDS = (0, 4)


def _read_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _resolve(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--panel-manifest", type=Path, required=True)
    parser.add_argument("--effect-panel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    panel_path = _resolve(root, str(args.panel_manifest))
    effect_panel_path = _resolve(root, str(args.effect_panel))
    panel = _read_json(panel_path)
    effect_panel = _read_json(effect_panel_path)
    if panel.get("schema") != "r04.restart_panel_manifest.v1":
        raise ValueError("unexpected frozen K=3 panel manifest schema")
    if effect_panel.get("status") != "COMPLETE_TRAINING_ROLE_ONLY":
        raise ValueError("training effect panel must be complete before K=2 bridge")
    panel_entries = {
        (int(item["restart_index"]), int(item["fold"])): item
        for item in effect_panel.get("entries", [])
        if isinstance(item, dict)
    }
    source_entries = {
        (int(item["restart_index"]), int(item["fold"])): item
        for item in panel.get("entries", [])
        if isinstance(item, dict)
    }
    expected_keys = {(restart, fold) for restart in RESTARTS for fold in FOLDS}
    if not expected_keys.issubset(panel_entries) or not expected_keys.issubset(source_entries):
        raise ValueError("K=3 source panels are incomplete for stage A")

    output_root = _resolve(root, str(args.output_root))
    entries: list[dict[str, object]] = []
    for restart, fold in sorted(expected_keys):
        source = source_entries[(restart, fold)]
        source_cell_path = _resolve(root, str(source["k3_cell"]))
        source_cell = _read_json(source_cell_path)
        if source_cell.get("k_model") != 3:
            raise ValueError("bridge source is not K=3")
        source_effect = panel_entries[(restart, fold)]
        fit_cell_path = _resolve(root, str(source_effect["source_fit_cell"]["path"]))
        fit_cell = _read_json(fit_cell_path)
        fit_diagnostics = fit_cell.get("fit_diagnostics")
        parameters = source_cell.get("parameters")
        if not isinstance(fit_diagnostics, dict) or not isinstance(parameters, dict):
            raise ValueError("K=3 source fit protocol is incomplete")
        required = {
            "gene_folds": 2,
            "inducing_points": 16,
            "inference_steps": 2400,
            "lengthscale": 3.0,
            "nonspatial_rank": 1,
            "optimization_schedule": "joint",
            "shared_steps": 0,
            "split_seed": 20260817,
        }
        for key, expected in required.items():
            if parameters.get(key) != expected:
                raise ValueError(f"K=3 source protocol differs at {key}")
        if fit_diagnostics.get("gene_batch_size") != 512:
            raise ValueError("K=3 source gene minibatch differs from frozen protocol")
        fit_steps = fit_diagnostics.get("steps")
        if not isinstance(fit_steps, int) or fit_steps < 1:
            raise ValueError("K=3 source fit step is missing")
        source_scores = source_cell.get("patient_scores")
        if not isinstance(source_scores, dict) or not source_scores:
            raise ValueError("K=3 source patient scores are missing")
        entries.append({
            "restart_index": restart,
            "fold": fold,
            "optimization_seed": int(parameters["optimization_seed"]),
            "fit_steps": fit_steps,
            "k2_output_dir": str(output_root / f"r{restart}_fold{fold}"),
            "k3_source_cell": {"path": str(source_cell_path), "sha256": _sha256(source_cell_path)},
            "k3_source_fit_cell": {"path": str(fit_cell_path), "sha256": _sha256(fit_cell_path)},
            "k3_source_input_hash": source_cell.get("input_hash"),
            "k3_source_config_hash": source_cell.get("config_hash"),
            "k3_source_patient_scores": {str(key): float(value) for key, value in source_scores.items()},
            "fit_protocol": {
                "k_model": 2,
                "comparison_k_model": 3,
                "steps": fit_steps,
                "inference_steps": 2400,
                "inducing_points": 16,
                "lengthscale": 3.0,
                "learning_rate": 0.01,
                "gene_batch_size": 512,
                "gene_folds": 2,
                "split_seed": 20260817,
                "optimization_seed": int(parameters["optimization_seed"]),
                "optimization_schedule": "joint",
                "shared_steps": 0,
                "nonspatial_rank": 1,
                "diagnostic_interval": int(fit_diagnostics.get("diagnostic_interval", 0)),
                "evaluation_interval": int(fit_diagnostics.get("evaluation_interval", 0)),
                "evaluation_mc_draws": int(fit_diagnostics.get("evaluation_mc_draws", 4)),
                "environment_hash": str(fit_diagnostics.get("environment_hash", "")),
                "execution_mode": "eager",
            },
            "source_effect_export": {
                "path": str(source_effect["structure_exports"][0]["npz_path"]),
                "sha256": source_effect["structure_exports"][0].get("npz_sha256"),
            },
        })

    output_path = _resolve(root, str(args.output))
    payload = {
        "schema": "r04.k2_bridge_stage_a_manifest.v1",
        "status": "FROZEN_PRECOMPUTE_MANIFEST",
        "created_at": "2026-09-01",
        "trigger": "training_role_readout_cannot_judge_third_direction_with_two_paired_patients_per_extreme_fold",
        "scope": "stage_a_extreme_folds_three_seeds",
        "selected_k": None,
        "k_interpretation": "latent_spatial_effect_dimension_not_structure_count",
        "panel_manifest": {"path": str(panel_path), "sha256": _sha256(panel_path)},
        "training_effect_panel": {"path": str(effect_panel_path), "sha256": _sha256(effect_panel_path)},
        "training_manifest": effect_panel.get("training_manifest"),
        "gene_list": {"path": str(root / "infra/r04/gene_universe.txt"), "sha256": _sha256(root / "infra/r04/gene_universe.txt")},
        "output_root": str(output_root),
        "gpu_visible": False,
        "validation_gt": "NOT_READ",
        "internal_external_validation_gt": "SEALED_NOT_READ",
        "entries": entries,
        "stop_conditions": {
            "stop_after_stage_a": "If K=2 and K=3 give the same main structure-readout direction and K=3 third direction has no stable incremental evidence, do not expand to all folds/seeds.",
            "expand_to_stage_b": "Only if K=3 has a clear seed-consistent patient-level held-out increment and a concordant structure-specific readout increment beyond uncertainty.",
            "global_rank_unresolved": "If fold 0 and fold 4 require opposite effective ranks, report fold/patient heterogeneity or a conditional overcomplete representation; do not force one global K.",
        },
        "run_command": (
            "env CUDA_VISIBLE_DEVICES='' PYTHONPATH=. "
            "python scripts/r04_k2_bridge_stage_a.py --manifest "
            + str(_resolve(root, str(args.output)))
        ),
        "claim_boundary": "This manifest authorizes a conditional K=2 capacity bridge only; it cannot establish a biological structure count or a named factor.",
    }
    atomic_json(output_path, payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
