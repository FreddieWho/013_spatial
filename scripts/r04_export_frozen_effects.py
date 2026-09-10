#!/usr/bin/env python3
"""Recreate and persist held-out R-04 fields from frozen K=3 checkpoints.

This command restores an audited fit checkpoint and runs only the already
specified held-out inference.  It never updates the fitted parameters, reads
ground truth, or chooses K.  Outputs are intended for the later training-role
structure-readout diagnostic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np

from r04.continuation_platform_audit import checkpoint_bundle_fingerprint
from r04.runtime import atomic_json, resource_status
from scripts.r04_real_k_search import (
    _grouped_folds,
    _load_training,
    _score_fold,
)


EXPECTED_INFERENCE_STEPS = 2400
EXPECTED_GENE_FOLDS = 2
EXPECTED_GROUP_SEED = 20260807
EXPECTED_SPLIT_SEED = 20260817
EXPECTED_INDUCING_POINTS = 16
EXPECTED_LENGTHSCALE = 3.0
EXPECTED_LEARNING_RATE = 0.01
EXPECTED_GENE_BATCH_SIZE = 512


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


def _parse_ints(value: str, *, name: str) -> list[int]:
    try:
        values = sorted({int(token.strip()) for token in value.split(",") if token.strip()})
    except ValueError as exc:
        raise SystemExit(f"{name} must contain comma-separated integers") from exc
    if not values:
        raise SystemExit(f"{name} must not be empty")
    return values


def _checkpoint_step(cell: dict[str, object], checkpoint_dir: Path) -> int:
    diagnostics = cell.get("fit_diagnostics")
    if not isinstance(diagnostics, dict) or not isinstance(diagnostics.get("steps"), int):
        raise ValueError("source cell has no integer fit step")
    step = int(diagnostics["steps"])
    pointer = (checkpoint_dir / "checkpoint").read_text(encoding="utf-8")
    marker = 'model_checkpoint_path: "ckpt-'
    if marker not in pointer:
        raise ValueError("checkpoint pointer has no ckpt step")
    pointer_step = pointer.split(marker, 1)[1].split('"', 1)[0]
    if not pointer_step.isdigit() or int(pointer_step) != step:
        raise ValueError("source cell and checkpoint pointer step differ")
    return step


def _validate_source(
    *,
    root: Path,
    entry: dict[str, object],
    panel: dict[str, object],
) -> tuple[Path, dict[str, object], Path, dict[str, object], Path, int]:
    cell_path = _resolve(root, str(entry["k3_cell"]))
    checkpoint_dir = _resolve(root, str(entry["k3_checkpoint_dir"]))
    if not cell_path.is_file() or not checkpoint_dir.is_dir():
        raise ValueError("source cell or checkpoint directory is missing")
    cell = _read_json(cell_path)
    if cell.get("schema") != "r04.real_k_cell.v1" or cell.get("k_model") != 3:
        raise ValueError("source cell is not a K=3 real-cell artifact")
    if cell.get("fold") != entry.get("fold"):
        raise ValueError("panel entry and source cell fold differ")
    parameters = cell.get("parameters")
    if not isinstance(parameters, dict):
        raise ValueError("source cell parameters are missing")
    expected = {
        "gene_folds": EXPECTED_GENE_FOLDS,
        "inducing_points": EXPECTED_INDUCING_POINTS,
        "lengthscale": EXPECTED_LENGTHSCALE,
        "inference_steps": EXPECTED_INFERENCE_STEPS,
        "split_seed": EXPECTED_SPLIT_SEED,
        "optimization_schedule": "joint",
        "shared_steps": 0,
    }
    for key, value in expected.items():
        if parameters.get(key) != value:
            raise ValueError(f"source cell protocol differs at {key}")
    diagnostics = cell.get("fit_diagnostics")
    if not isinstance(diagnostics, dict) or diagnostics.get("gene_batch_size") != EXPECTED_GENE_BATCH_SIZE:
        raise ValueError("source cell gene minibatch differs from frozen protocol")
    if cell.get("inference_platform", {}).get("converged") is not True:
        raise ValueError("source cell inference platform is not converged")
    metadata = _read_json(checkpoint_dir / "checkpoint_metadata.json")
    torch_source = metadata.get("backend") == "torch"
    if torch_source:
        # Torch-era single-file checkpoint (D-100): no TF pointer/bundle.
        if metadata.get("checkpoint_schema") != "r04.mnsf_checkpoint.v3_torch":
            raise ValueError("source Torch checkpoint schema is not supported")
        if metadata.get("factors") != cell.get("k_model"):
            raise ValueError("source Torch checkpoint factor count differs from cell")
        checkpoint_file = checkpoint_dir / "checkpoint.pt"
        if not checkpoint_file.is_file():
            raise ValueError("source Torch checkpoint file is missing")
        fit_diag = cell.get("fit_diagnostics")
        if not isinstance(fit_diag, dict):
            raise ValueError("source cell has no fit diagnostics")
        if metadata.get("objective_input_hash") != fit_diag.get("objective_input_hash"):
            raise ValueError("source Torch checkpoint objective input hash differs")
        if not metadata.get("environment_hash"):
            raise ValueError("source Torch checkpoint metadata has no environment hash")
    elif metadata.get("checkpoint_schema") != "r04.mnsf_checkpoint.v2":
        raise ValueError("source checkpoint schema is not supported")
    fit_cell_path = checkpoint_dir.parents[2] / "cell.json"
    if not fit_cell_path.is_file():
        # The first restart diagnostic stores its fit cell below ``cells/``;
        # the panel entry is the authoritative fit/score record in that case.
        fit_cell_path = cell_path
    fit_cell = _read_json(fit_cell_path)
    fit_diagnostics = fit_cell.get("fit_diagnostics")
    if not isinstance(fit_diagnostics, dict):
        raise ValueError("source fit cell has no diagnostics")
    for key in ("input_hash", "config_hash", "environment_hash"):
        if fit_cell.get(key) and metadata.get(key) != fit_cell.get(key):
            raise ValueError(f"source fit cell and checkpoint {key} differ")
    for key in ("input_hash", "config_hash"):
        if cell.get(key) and fit_cell.get(key) and cell.get(key) != fit_cell.get(key):
            if key == "input_hash":
                raise ValueError(f"panel cell and source fit cell {key} differ")
    if torch_source:
        diagnostics = cell.get("fit_diagnostics")
        if not isinstance(diagnostics, dict) or not isinstance(diagnostics.get("steps"), int):
            raise ValueError("source cell has no integer fit step")
        step = int(diagnostics["steps"])
        bundle = {"sha256": _sha256(checkpoint_dir / "checkpoint.pt"), "kind": "torch_single_file_v3"}
    else:
        step = _checkpoint_step(cell, checkpoint_dir)
        bundle = checkpoint_bundle_fingerprint(checkpoint_dir, step=step)
        if not bundle.get("sha256"):
            raise ValueError("source checkpoint bundle has no hash")
    return cell_path, cell, fit_cell_path, fit_cell, checkpoint_dir, step, str(metadata.get("environment_hash")), bundle


def _write_panel_status(path: Path, payload: dict[str, object]) -> None:
    atomic_json(path, payload)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--panel-manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--restart-indices", default="0,1,2,3,4")
    parser.add_argument("--fold-values", default="0,4")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--training-only",
        action="store_true",
        help="export only frozen training-role effects; do not run held-out inference or score replay",
    )
    args = parser.parse_args()

    # R-04 currently forbids GPU use.  This is set before the estimator is
    # imported and is recorded in the panel provenance.
    os.environ["CUDA_VISIBLE_DEVICES"] = ""

    root = args.project_root.resolve()
    panel_path = _resolve(root, str(args.panel_manifest))
    output_root = _resolve(root, str(args.output_root))
    panel = _read_json(panel_path)
    entries = panel.get("entries")
    if not isinstance(entries, list) or not entries:
        raise SystemExit("panel manifest has no entries")
    restart_indices = set(_parse_ints(args.restart_indices, name="restart-indices"))
    fold_values = set(_parse_ints(args.fold_values, name="fold-values"))
    selected = [
        entry for entry in entries
        if isinstance(entry, dict)
        and entry.get("restart_index") in restart_indices
        and entry.get("fold") in fold_values
    ]
    if len(selected) != len(restart_indices) * len(fold_values):
        raise SystemExit("requested restart/fold entries are incomplete")
    output_root.mkdir(parents=True, exist_ok=True)
    if resource_status(output_root) == "BLOCKED_STORAGE":
        _write_panel_status(output_root / "panel.json", {"schema": "r04.frozen_effect_export_panel.v1", "status": "BLOCKED_STORAGE"})
        return 2

    manifest_path = _resolve(root, str(panel.get("training_manifest", "")))
    gene_path = root / "infra/r04/gene_universe.txt"
    sections, genes, manifest_hash = _load_training(manifest_path, gene_path)
    fold_ids = _grouped_folds(sections, folds=5, seed=EXPECTED_GROUP_SEED)
    panel_output = {
        "schema": "r04.frozen_effect_export_panel.v1",
        "status": "PREFLIGHT" if args.dry_run else "RUNNING",
        "purpose": (
            "training_role_structure_readout_effect_export"
            if args.training_only
            else "heldout_field_export_for_training_role_structure_readout"
        ),
        "scope": "training_role_only" if args.training_only else "training_and_outer_heldout",
        "k_model": 3,
        "selected_k": None,
        "panel_manifest": {"path": str(panel_path), "sha256": _sha256(panel_path)},
        "training_manifest": {"path": str(manifest_path), "sha256": _sha256(manifest_path)},
        "gene_list": {"path": str(gene_path), "sha256": _sha256(gene_path), "count": len(genes)},
        "group_seed": EXPECTED_GROUP_SEED,
        "split_seed": EXPECTED_SPLIT_SEED,
        "entries": [],
        "validation_gt": "NOT_READ",
        "internal_external_validation_gt": "SEALED_NOT_READ",
        "gpu_visible": False,
        "fit_updates": 0,
    }
    _write_panel_status(output_root / "panel.json", panel_output)

    for entry in sorted(selected, key=lambda item: (int(item["restart_index"]), int(item["fold"]))):
        restart = int(entry["restart_index"])
        fold = int(entry["fold"])
        (
            cell_path,
            source_cell,
            fit_cell_path,
            fit_cell,
            checkpoint_dir,
            fit_steps,
            source_environment_hash,
            checkpoint_bundle,
        ) = _validate_source(
            root=root, entry=entry, panel=panel
        )
        source_cell_hash = _sha256(cell_path)
        record = {
            "restart_index": restart,
            "fold": fold,
            "source_cell": {"path": str(cell_path), "sha256": source_cell_hash},
            "source_fit_cell": {"path": str(fit_cell_path), "sha256": _sha256(fit_cell_path)},
            "source_checkpoint": {"path": str(checkpoint_dir), **checkpoint_bundle},
            "source_environment_hash": source_environment_hash,
            "source_status": source_cell.get("status"),
            "source_inference_authority": source_cell.get("inference_authority"),
            "strict_audit_status": source_cell.get("strict_audit_status"),
        }
        cell_root = output_root / f"r{restart}_fold{fold}"
        export_root = cell_root / "exports"
        structure_export_root = cell_root / "structure_exports"
        replay_checkpoint_root = cell_root / "replay_checkpoints"
        cell_root.mkdir(parents=True, exist_ok=True)
        if args.dry_run:
            record["status"] = "PREFLIGHT_OK"
            panel_output["entries"].append(record)
            _write_panel_status(output_root / "panel.json", panel_output)
            continue
        training = [
            section for index, section in enumerate(sections) if int(fold_ids[index]) != fold
        ]
        validation = [
            section for index, section in enumerate(sections) if int(fold_ids[index]) == fold
        ]
        parameters = source_cell["parameters"]
        fit_diagnostics = fit_cell["fit_diagnostics"]
        seed = int(parameters["optimization_seed"])
        try:
            result = _score_fold(
                training,
                validation,
                genes=genes,
                k_model=3,
                fold=fold,
                output_dir=cell_root,
                steps=fit_steps,
                inference_steps=EXPECTED_INFERENCE_STEPS,
                inducing_points=EXPECTED_INDUCING_POINTS,
                lengthscale=EXPECTED_LENGTHSCALE,
                learning_rate=EXPECTED_LEARNING_RATE,
                learning_rate_final=None,
                learning_rate_decay_steps=None,
                gene_batch_size=EXPECTED_GENE_BATCH_SIZE,
                diagnostic_interval=int(fit_diagnostics["diagnostic_interval"]),
                evaluation_interval=int(fit_diagnostics["evaluation_interval"]),
                evaluation_mc_draws=int(fit_diagnostics["evaluation_mc_draws"]),
                optimization_schedule="joint",
                shared_steps=0,
                gene_folds=EXPECTED_GENE_FOLDS,
                split_seed=EXPECTED_SPLIT_SEED,
                optimization_seed=seed,
                checkpoint_steps=fit_steps,
                resume_checkpoint_dir=checkpoint_dir,
                environment_hash=source_environment_hash,
                inference_export_root=export_root,
                structure_export_root=structure_export_root,
                checkpoint_output_dir=replay_checkpoint_root / "k3" / f"fold{fold}",
                export_provenance={
                    "source_cell_path": str(cell_path),
                    "source_cell_sha256": source_cell_hash,
                    "source_fit_cell_path": str(fit_cell_path),
                    "source_fit_cell_sha256": _sha256(fit_cell_path),
                    "source_checkpoint_dir": str(checkpoint_dir),
                    "source_checkpoint_bundle_sha256": checkpoint_bundle["sha256"],
                    "source_checkpoint_step": fit_steps,
                    "source_status": source_cell.get("status"),
                    "source_inference_authority": source_cell.get("inference_authority"),
                    "strict_audit_status": source_cell.get("strict_audit_status"),
                    "fit_updates": 0,
                    "gpu_visible": False,
                    "validation_gt_read": "NOT_READ",
                    "internal_external_validation_gt": "SEALED_NOT_READ",
                },
                # Use the only currently implemented Torch execution path.
                execution_mode="eager",
                training_only=args.training_only,
            )
        except Exception as exc:
            record.update({"status": "FAILED", "error": type(exc).__name__, "message": str(exc)})
            panel_output["entries"].append(record)
            panel_output["status"] = "FAILED"
            _write_panel_status(output_root / "panel.json", panel_output)
            raise
        replay_fit_diagnostics = result.get("fit_diagnostics")
        if not isinstance(replay_fit_diagnostics, dict):
            raise RuntimeError("frozen replay has no fit diagnostics")
        if replay_fit_diagnostics.get("optimizer_steps_this_call") != 0:
            raise RuntimeError("frozen replay performed fit updates")
        if args.training_only:
            if result.get("patient_scores") != {}:
                raise RuntimeError("training-only export unexpectedly produced held-out scores")
            if result.get("input_hash") != source_cell.get("input_hash"):
                raise RuntimeError("training-only replay input hash differs from source cell")
            if result.get("config_hash") != source_cell.get("config_hash"):
                raise RuntimeError("training-only replay config hash differs from source cell")
            score_replay_status = "NOT_RUN_TRAINING_ONLY_SOURCE_SCORES_RETAINED"
            record_status = "TRAINING_EFFECT_EXPORTED"
            score_match = None
        else:
            replay_scores = result.get("patient_scores")
            source_scores = source_cell.get("patient_scores")
            if not isinstance(replay_scores, dict) or not isinstance(source_scores, dict):
                raise RuntimeError("source or replay patient scores are missing")
            if set(replay_scores) != set(source_scores) or not all(
                np.isclose(float(replay_scores[key]), float(source_scores[key]), rtol=1e-6, atol=1e-5)
                for key in replay_scores
            ):
                diffs = {
                    key: {
                        "replay": float(replay_scores[key]),
                        "source": float(source_scores[key]),
                        "abs_diff": abs(float(replay_scores[key]) - float(source_scores[key])),
                    }
                    for key in replay_scores if key in source_scores
                }
                worst = max(diffs, key=lambda k: diffs[k]["abs_diff"]) if diffs else None
                record.update({
                    "status": "FAILED_SCORE_MISMATCH_FIELDS_RETAINED",
                    "replay_patient_scores": {k: float(v) for k, v in replay_scores.items()},
                    "score_diffs": diffs,
                    "worst_key": worst,
                })
                panel_output["entries"].append(record)
                _write_panel_status(output_root / "panel.json", panel_output)
                raise RuntimeError(
                    "frozen replay patient scores differ from source cell: "
                    f"worst={worst} diff={diffs[worst]['abs_diff']:.6f}" if worst else
                    "frozen replay patient scores differ from source cell"
                )
            score_replay_status = "MATCHED_SOURCE_CELL"
            record_status = "EXPORTED"
            score_match = True
        record.update({
            "status": record_status,
            "source_cell_status": source_cell.get("status"),
            "exports": result.get("inference_field_exports", []),
            "structure_exports": result.get("structure_field_exports", []),
            "replay_input_hash": result.get("input_hash"),
            "replay_config_hash": result.get("config_hash"),
            "replay_fit_diagnostics": result.get("fit_diagnostics"),
            "replay_inference_platform": result.get("inference_platform"),
            "fit_updates": replay_fit_diagnostics.get("optimizer_steps_this_call"),
            "replay_score_match": score_match,
            "score_replay_status": score_replay_status,
        })
        atomic_json(cell_root / "export_cell.json", record)
        panel_output["entries"].append(record)
        _write_panel_status(output_root / "panel.json", panel_output)
    panel_output["status"] = (
        "PREFLIGHT_COMPLETE"
        if args.dry_run
        else "COMPLETE_TRAINING_ROLE_ONLY"
        if args.training_only
        else "COMPLETE"
    )
    panel_output["manifest_hash_recomputed"] = manifest_hash
    _write_panel_status(output_root / "panel.json", panel_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
