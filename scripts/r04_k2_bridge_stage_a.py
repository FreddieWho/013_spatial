#!/usr/bin/env python3
"""Run the frozen, minimal K=2 bridge on extreme folds.

This is conditional R-04 work, not a general K search.  The only model change
from the K=3 source protocol is ``factors=2``; no GT is read by the fit or
held-out score path.  Training-role GT is read only after all six cells have
finished, for the descriptive inner readout comparison.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np

os.environ["CUDA_VISIBLE_DEVICES"] = ""

from r04.field_exports import read_heldout_field_export
from r04.runtime import atomic_json, resource_status
from r04.structure_gt import labels_for_spots, load_training_gt_catalog
from scripts.r04_real_k_search import _grouped_folds, _load_training, _score_fold
from scripts.r04_run_structure_readout import _readout_one_rank


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


def _score_summary(k2: dict[str, float], k3: dict[str, float]) -> dict[str, object]:
    common = sorted(set(k2) & set(k3))
    if not common:
        raise ValueError("K=2 and K=3 have no common patient scores")
    deltas = np.asarray([float(k2[p]) - float(k3[p]) for p in common], dtype=float)
    return {
        "common_patients": common,
        "k2_minus_k3_by_patient": {p: float(k2[p] - k3[p]) for p in common},
        "mean_k2_minus_k3": float(np.mean(deltas)),
        "median_k2_minus_k3": float(np.median(deltas)),
        "k2_better_patient_count": int(np.sum(deltas > 0)),
        "k3_better_patient_count": int(np.sum(deltas < 0)),
        "tie_patient_count": int(np.sum(deltas == 0)),
        "interpretation": "descriptive paired held-out NB score difference; not a K selection rule",
    }


def _target_labels(catalog: dict[str, object], export: dict[str, object]) -> dict[str, np.ndarray]:
    return {
        name: labels_for_spots(
            catalog,
            name,
            export["section_ids"],
            export["barcodes"],
        )
        for name in ("TLS", "TUMOR_STROMA_BOUNDARY")
    }


def _readout_pair(
    k2_export: dict[str, object],
    k3_export: dict[str, object],
    catalog: dict[str, object],
) -> dict[str, object]:
    k2_labels = _target_labels(catalog, k2_export)
    k3_labels = _target_labels(catalog, k3_export)
    k2 = _readout_one_rank(k2_export, None, k2_labels, max_rank=None)
    k3 = _readout_one_rank(k3_export, None, k3_labels, max_rank=None)
    return {
        "k2_full": k2,
        "k3_full": k3,
        "k3_stable_rank2": _readout_one_rank(
            k3_export, None, k3_labels, max_rank=2
        ),
        "paired_training_patient_count": len(
            set(k2.get("training_patients", [])) & set(k3.get("training_patients", []))
        ),
        "interpretation": (
            "training-role inner cross-fit comparison only; structure names are targets, "
            "not factor labels, and outer validation GT remains sealed"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    root = args.project_root.resolve()
    manifest_path = _resolve(root, str(args.manifest))
    manifest = _read_json(manifest_path)
    if manifest.get("schema") != "r04.k2_bridge_stage_a_manifest.v1":
        raise SystemExit("unexpected K=2 bridge manifest schema")
    if manifest.get("status") != "FROZEN_PRECOMPUTE_MANIFEST":
        raise SystemExit("K=2 bridge manifest is not frozen")
    output_root = _resolve(root, str(manifest["output_root"]))
    output_path = _resolve(
        root,
        str(args.output) if args.output is not None else str(output_root / "stage_a.json"),
    )
    if resource_status(output_root) == "BLOCKED_STORAGE":
        atomic_json(output_path, {"schema": "r04.k2_bridge_stage_a.v1", "status": "BLOCKED_STORAGE"})
        return 2

    training_manifest = _resolve(root, str(manifest["training_manifest"]["path"]))
    gene_path = _resolve(root, str(manifest["gene_list"]["path"]))
    sections, genes, manifest_hash = _load_training(training_manifest, gene_path)
    fold_ids = _grouped_folds(sections, folds=5, seed=20260807)
    k3_effect_panel = _read_json(_resolve(root, str(manifest["training_effect_panel"]["path"])))
    k3_effect_entries = {
        (int(item["restart_index"]), int(item["fold"])): item
        for item in k3_effect_panel["entries"]
        if isinstance(item, dict)
    }
    payload: dict[str, object] = {
        "schema": "r04.k2_bridge_stage_a.v1",
        "status": "RUNNING",
        "manifest": {"path": str(manifest_path), "sha256": _sha256(manifest_path)},
        "selected_k": None,
        "k_interpretation": "latent_spatial_effect_dimension_not_structure_count",
        "scope": manifest["scope"],
        "gpu_visible": False,
        "validation_gt": "NOT_READ",
        "internal_external_validation_gt": "SEALED_NOT_READ",
        "training_manifest_hash_recomputed": manifest_hash,
        "entries": [],
        "stop_conditions": manifest["stop_conditions"],
        "claim_boundary": manifest["claim_boundary"],
    }
    atomic_json(output_path, payload)

    for item in manifest["entries"]:
        restart = int(item["restart_index"])
        fold = int(item["fold"])
        protocol = item["fit_protocol"]
        output_dir = _resolve(root, str(item["k2_output_dir"]))
        checkpoint_dir = output_dir / "checkpoints" / "k2" / f"fold{fold}"
        structure_root = output_dir / "structure_exports"
        output_dir.mkdir(parents=True, exist_ok=True)
        training = [section for index, section in enumerate(sections) if int(fold_ids[index]) != fold]
        validation = [section for index, section in enumerate(sections) if int(fold_ids[index]) == fold]
        try:
            result = _score_fold(
                training,
                validation,
                genes=genes,
                k_model=2,
                fold=fold,
                output_dir=output_dir,
                steps=int(protocol["steps"]),
                inference_steps=int(protocol["inference_steps"]),
                inducing_points=int(protocol["inducing_points"]),
                lengthscale=float(protocol["lengthscale"]),
                learning_rate=float(protocol["learning_rate"]),
                learning_rate_final=None,
                learning_rate_decay_steps=None,
                gene_batch_size=int(protocol["gene_batch_size"]),
                diagnostic_interval=int(protocol["diagnostic_interval"]),
                evaluation_interval=int(protocol["evaluation_interval"]),
                evaluation_mc_draws=int(protocol["evaluation_mc_draws"]),
                optimization_schedule=str(protocol["optimization_schedule"]),
                shared_steps=int(protocol["shared_steps"]),
                gene_folds=int(protocol["gene_folds"]),
                split_seed=int(protocol["split_seed"]),
                optimization_seed=int(protocol["optimization_seed"]),
                checkpoint_steps=int(protocol["steps"]),
                checkpoint_output_dir=checkpoint_dir,
                environment_hash=str(protocol["environment_hash"]),
                structure_export_root=structure_root,
                export_heldout_effect=False,
                execution_mode="eager",
                export_provenance={
                    "bridge_stage": "A",
                    "comparison_k_model": 3,
                    "source_k3_cell_path": item["k3_source_cell"]["path"],
                    "source_k3_cell_sha256": item["k3_source_cell"]["sha256"],
                    "source_k3_fit_cell_path": item["k3_source_fit_cell"]["path"],
                    "source_k3_fit_cell_sha256": item["k3_source_fit_cell"]["sha256"],
                    "validation_gt_read": "NOT_READ",
                    "internal_external_validation_gt": "SEALED_NOT_READ",
                    "gpu_visible": False,
                },
            )
        except Exception as exc:
            failed = {
                "restart_index": restart,
                "fold": fold,
                "status": "FAILED",
                "error": type(exc).__name__,
                "message": str(exc),
            }
            payload["entries"].append(failed)
            payload["status"] = "FAILED"
            atomic_json(output_path, payload)
            raise
        diagnostics = result.get("fit_diagnostics")
        if not isinstance(diagnostics, dict) or diagnostics.get("optimizer_steps_this_call") != int(protocol["steps"]):
            raise RuntimeError("K=2 stage A fit did not complete the frozen step count")
        inference_platform = result.get("inference_platform")
        if not isinstance(inference_platform, dict) or inference_platform.get("converged") is not True:
            record = {
                "restart_index": restart,
                "fold": fold,
                "status": "INFERENCE_NOT_PLATFORMED",
                "fit_updates": diagnostics.get("optimizer_steps_this_call"),
                "inference_platform": inference_platform,
            }
            payload["entries"].append(record)
            payload["status"] = "INCOMPLETE_INFERENCE_NOT_PLATFORMED"
            atomic_json(output_path, payload)
            raise RuntimeError("K=2 stage A inference platform did not pass")
        record = {
            "restart_index": restart,
            "fold": fold,
            "status": "K2_FIT_AND_SCORED",
            "fit_updates": diagnostics.get("optimizer_steps_this_call"),
            "fit_platform": result.get("fit_platform"),
            "inference_platform": inference_platform,
            "inference_platform_required": True,
            "patient_scores": result.get("patient_scores"),
            "score_comparison": _score_summary(
                {str(k): float(v) for k, v in result["patient_scores"].items()},
                {str(k): float(v) for k, v in item["k3_source_patient_scores"].items()},
            ),
            "cell": {"path": str(output_dir / "cell.json")},
            "structure_exports": result.get("structure_field_exports", []),
        }
        cell_payload = {
            **result,
            "stage": "K2_BRIDGE_STAGE_A",
            "manifest_path": str(manifest_path),
            "manifest_sha256": _sha256(manifest_path),
            "validation_gt": "NOT_READ",
            "internal_external_validation_gt": "SEALED_NOT_READ",
            "gpu_visible": False,
            "fit_updates": diagnostics.get("optimizer_steps_this_call"),
        }
        atomic_json(output_dir / "cell.json", cell_payload)
        record["cell"]["sha256"] = _sha256(output_dir / "cell.json")
        payload["entries"].append(record)
        atomic_json(output_path, payload)

    if len(payload["entries"]) != len(manifest["entries"]):
        raise RuntimeError("K=2 stage A entries are incomplete")

    # GT is intentionally loaded only after all model fits and held-out scores.
    catalog = load_training_gt_catalog(
        _resolve(root, "infra/structure-registry/structure_instances.tsv"),
        training_manifest,
    )
    readouts = []
    for item, record in zip(manifest["entries"], payload["entries"]):
        restart = int(item["restart_index"])
        fold = int(item["fold"])
        k2_metadata = record["structure_exports"][0]
        k3_metadata = k3_effect_entries[(restart, fold)]["structure_exports"][0]
        k2_export = read_heldout_field_export(_resolve(root, str(k2_metadata["npz_path"])))
        k3_export = read_heldout_field_export(_resolve(root, str(k3_metadata["npz_path"])))
        readouts.append({
            "restart_index": restart,
            "fold": fold,
            "readout": _readout_pair(k2_export, k3_export, catalog),
        })
    payload["structure_readout"] = {
        "status": "COMPLETE_DESCRIPTIVE_TRAINING_ROLE_INNER_CROSSFIT",
        "entries": readouts,
        "internal_external_validation_gt": "SEALED_NOT_READ",
        "factor_naming": "forbidden",
    }
    payload["stage_a_decision"] = "NOT_CLOSED_INSUFFICIENT_PAIRED_GT"
    payload["downstream_k_robust"] = "not_tested"
    payload["k2_bridge_required"] = True
    payload["next_action"] = (
        "Do not expand to stage B from this diagnostic alone; the two-patient paired "
        "training GT support cannot establish a general K=2/K=3 structure conclusion."
    )
    payload["status"] = "COMPLETE_STAGE_A_NOT_CLOSED"
    atomic_json(output_path, payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
