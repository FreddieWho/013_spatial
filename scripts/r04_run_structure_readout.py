#!/usr/bin/env python3
"""Run the minimal training-role structure readout on frozen effect exports.

The effect panel must be complete before this command reads any structure GT.
Only the training-role confirmatory annotations selected by the registry are
loaded; internal/external validation GT remains sealed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from r04.field_exports import read_heldout_field_export
from r04.runtime import atomic_json
from r04.structure_gt import labels_for_spots, load_training_gt_catalog
from r04.structure_readout import effect_coordinates, evaluate_shared_specific_readout


STRUCTURES = ("TLS", "TUMOR_STROMA_BOUNDARY")


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


def _balanced_patient_weights(patient_ids: np.ndarray) -> np.ndarray:
    values = np.asarray(patient_ids).astype(str)
    unique = np.unique(values)
    if len(unique) == 0:
        raise ValueError("no patients are available for readout")
    weights = np.zeros(len(values), dtype=float)
    for patient in unique:
        indices = np.flatnonzero(values == patient)
        weights[indices] = 1.0 / (len(unique) * len(indices))
    return weights


def _coverage(labels: dict[str, np.ndarray], patient_ids: np.ndarray) -> dict[str, object]:
    result: dict[str, object] = {}
    for structure, values in labels.items():
        finite = np.isfinite(values)
        result[structure] = {
            "annotated_spots": int(np.sum(finite)),
            "patients": sorted(set(patient_ids[finite].astype(str))),
            "positive_spots": int(np.sum(values[finite] == 1)),
            "negative_spots": int(np.sum(values[finite] == 0)),
        }
    paired = np.isfinite(np.vstack([labels[name] for name in STRUCTURES])).all(axis=0)
    result["paired"] = {
        "spots": int(np.sum(paired)),
        "patients": sorted(set(patient_ids[paired].astype(str))),
    }
    return result


def _inner_patient_folds(patient_ids: np.ndarray, *, max_folds: int = 3) -> list[tuple[np.ndarray, np.ndarray, list[str]]]:
    """Return deterministic patient-level folds for training-role cross-fitting."""
    values = np.asarray(patient_ids).astype(str)
    patients = sorted(set(values.tolist()))
    n_folds = min(max_folds, len(patients))
    if n_folds < 2:
        return []
    fold_by_patient = {
        patient: index % n_folds for index, patient in enumerate(patients)
    }
    result = []
    for fold in range(n_folds):
        test = np.asarray([fold_by_patient[patient] == fold for patient in values])
        train = ~test
        result.append((train, test, [patient for patient in patients if fold_by_patient[patient] == fold]))
    return result


def _aggregate_inner_readout(folds: list[dict[str, object]]) -> dict[str, object]:
    valid = [item for item in folds if item.get("status") == "COMPUTED"]
    if not valid:
        return {
            "status": "NOT_TESTABLE_NO_VALID_INNER_FOLD",
            "valid_inner_folds": 0,
            "total_inner_folds": len(folds),
        }

    def mean_path(path: tuple[str, ...]) -> float | None:
        values: list[float] = []
        for item in valid:
            current: object = item
            for key in path:
                if not isinstance(current, dict) or key not in current:
                    return None
                current = current[key]
            if current is not None:
                values.append(float(current))
        return float(np.mean(values)) if values else None

    auc_delta = {}
    for structure in STRUCTURES:
        values = [
            float(item["readout"]["specific_increment"]["auc_delta_by_target"][structure])
            for item in valid
            if item["readout"]["specific_increment"]["auc_delta_by_target"][structure] is not None
        ]
        auc_delta[structure] = float(np.mean(values)) if values else None
    return {
        "status": "COMPUTED_DESCRIPTIVE_TRAINING_ROLE_INNER_CROSSFIT",
        "valid_inner_folds": len(valid),
        "total_inner_folds": len(folds),
        "shared_only_mean_mse": mean_path(("readout", "shared_only", "mean_mse")),
        "shared_plus_specific_mean_mse": mean_path(("readout", "shared_plus_specific_residual", "mean_mse")),
        "specific_minus_shared_mean_mse_delta": mean_path(("readout", "specific_increment", "mean_mse_delta_specific_minus_shared")),
        "mean_auc_delta_by_target": auc_delta,
        "fold_results": folds,
        "interpretation": (
            "descriptive training-role cross-fit only; no outer validation GT was read, "
            "and no threshold was introduced after observing these values"
        ),
    }


def _readout_one_rank(
    training_export: dict[str, object],
    heldout_export: dict[str, object] | None,
    training_labels: dict[str, np.ndarray],
    *,
    max_rank: int | None,
) -> dict[str, object]:
    has_outer_heldout = heldout_export is not None
    effect_target = heldout_export if heldout_export is not None else training_export
    train_z, test_z, representation = effect_coordinates(
        training_export, effect_target, max_rank=max_rank
    )
    train_valid = np.isfinite(np.vstack([training_labels[name] for name in STRUCTURES])).all(axis=0)
    if not np.any(train_valid):
        return {
            "status": "NOT_TESTABLE_INSUFFICIENT_PAIRED_TRAINING_GT",
            "representation": representation,
            "training_coverage": int(np.sum(train_valid)),
            "outer_heldout_gt": "SEALED_NOT_READ" if has_outer_heldout else "NOT_IN_SCOPE",
        }
    paired_z = train_z[train_valid]
    paired_labels = {
        name: np.asarray(training_labels[name], dtype=float)[train_valid]
        for name in STRUCTURES
    }
    paired_patients = np.asarray(training_export["patient_ids"])[train_valid]
    folds = _inner_patient_folds(paired_patients)
    if not folds:
        return {
            "status": "NOT_TESTABLE_INSUFFICIENT_PAIRED_TRAINING_PATIENTS",
            "representation": representation,
            "training_coverage": int(np.sum(train_valid)),
            "training_patients": sorted(set(paired_patients.astype(str))),
            "outer_heldout_gt": "SEALED_NOT_READ" if has_outer_heldout else "NOT_IN_SCOPE",
        }
    fold_results: list[dict[str, object]] = []
    for inner_index, (inner_train, inner_test, test_patients) in enumerate(folds):
        fold_output: dict[str, object] = {
            "inner_fold": inner_index,
            "test_patients": test_patients,
            "training_spots": int(np.sum(inner_train)),
            "heldout_spots": int(np.sum(inner_test)),
        }
        if any(
            len(np.unique(paired_labels[name][inner_train])) < 2
            or len(np.unique(paired_labels[name][inner_test])) < 2
            for name in STRUCTURES
        ):
            fold_output["status"] = "NOT_TESTABLE_NO_TARGET_VARIATION"
            fold_results.append(fold_output)
            continue
        readout = evaluate_shared_specific_readout(
            paired_z[inner_train],
            paired_z[inner_test],
            {name: paired_labels[name][inner_train] for name in STRUCTURES},
            {name: paired_labels[name][inner_test] for name in STRUCTURES},
            representation=representation,
            training_weights=_balanced_patient_weights(paired_patients[inner_train]),
            heldout_weights=_balanced_patient_weights(paired_patients[inner_test]),
        )
        fold_output.update({"status": "COMPUTED", "readout": readout})
        fold_results.append(fold_output)
    return {
        "status": "COMPUTED_DESCRIPTIVE_TRAINING_ROLE_INNER_CROSSFIT",
        "representation": representation,
        "training_coverage": int(np.sum(train_valid)),
        "training_patients": sorted(set(paired_patients.astype(str))),
        "outer_heldout_spots": int(len(test_z)) if has_outer_heldout else 0,
        "outer_heldout_patients": (
            sorted(set(np.asarray(heldout_export["patient_ids"]).astype(str)))
            if has_outer_heldout else []
        ),
        "outer_heldout_gt": "SEALED_NOT_READ" if has_outer_heldout else "NOT_IN_SCOPE",
        "inner_crossfit": _aggregate_inner_readout(fold_results),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--effect-panel", type=Path, required=True)
    parser.add_argument("--registry", type=Path, default=Path("infra/structure-registry/structure_instances.tsv"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    panel_path = _resolve(root, str(args.effect_panel))
    panel = _read_json(panel_path)
    if panel.get("schema") != "r04.frozen_effect_export_panel.v1" or panel.get("status") not in {
        "COMPLETE", "COMPLETE_TRAINING_ROLE_ONLY"
    }:
        raise SystemExit("effect panel must be complete before reading structure GT")
    training_only = panel.get("status") == "COMPLETE_TRAINING_ROLE_ONLY"
    entries = panel.get("entries")
    if not isinstance(entries, list) or not entries:
        raise SystemExit("effect panel has no entries")

    loaded: list[tuple[dict[str, object], dict[str, object], dict[str, object]]] = []
    for record in entries:
        allowed_entry_statuses = {"EXPORTED", "TRAINING_EFFECT_EXPORTED"}
        if not isinstance(record, dict) or record.get("status") not in allowed_entry_statuses:
            raise SystemExit("effect panel contains a non-exported entry")
        exports = record.get("structure_exports")
        expected_export_count = 1 if training_only else 2
        if not isinstance(exports, list) or len(exports) != expected_export_count:
            raise SystemExit("structure export count does not match panel scope")
        by_role: dict[str, dict[str, object]] = {}
        for metadata in exports:
            if not isinstance(metadata, dict):
                raise SystemExit("invalid structure export metadata")
            provenance = metadata.get("provenance")
            if not isinstance(provenance, dict) or provenance.get("validation_gt_read") != "NOT_READ":
                raise SystemExit("structure export provenance is not GT-sealed")
            role = str(provenance.get("representation_role", ""))
            key = "training" if role == "training_outer_fold_readout_fit" else "heldout" if role == "heldout_outer_fold_readout_apply" else ""
            if not key or key in by_role:
                raise SystemExit("structure export roles are incomplete or duplicated")
            by_role[key] = read_heldout_field_export(_resolve(root, str(metadata["npz_path"])))
        if "training" not in by_role or (not training_only and "heldout" not in by_role):
            raise SystemExit("structure export roles are incomplete")
        loaded.append((record, by_role["training"], by_role.get("heldout")))

    # No GT has been read before this point.  The effect panel is now fully
    # hash-checked and the readout representation is fixed by code/schema.
    registry_path = _resolve(root, str(args.registry))
    training_manifest = _resolve(root, str(panel["training_manifest"]["path"]))
    catalog = load_training_gt_catalog(registry_path, training_manifest)
    output = {
        "schema": "r04.structure_readout_panel.v1",
        "status": "RUNNING",
        "selected_k": None,
        "k_interpretation": "latent_spatial_effect_dimension_not_structure_count",
        "downstream_k_robust": "not_tested",
        "k2_bridge_required": True,
        "k2_bridge_reason": (
            "training-role paired GT has only two patients per extreme outer fold; "
            "the descriptive inner cross-fit cannot establish third-direction robustness"
        ),
        "spatial_null_status": "NOT_RUN",
        "composition_split_status": "NOT_RUN",
        "structure_naming_status": "FORBIDDEN_NOT_RUN",
        "effect_panel": {"path": str(panel_path), "sha256": _sha256(panel_path)},
        "training_gt_catalog": {
            "path": str(registry_path),
            "sha256": _sha256(registry_path),
            "schema": catalog["schema"],
            "training_manifest": str(training_manifest),
            "training_manifest_sha256": _sha256(training_manifest),
        },
        "structures": list(STRUCTURES),
        "scope": "training_role_inner_crossfit_only" if training_only else "training_role_and_outer_apply",
        "internal_external_validation_gt": "SEALED_NOT_READ",
        "factor_naming": "forbidden",
        "weighting": "patient_equal_within_each_outer_fold",
        "entries": [],
        "claim_boundary": (
            "descriptive readout-level evidence only; does not prove causal structure "
            "identity, independence from composition/markers, or a unique K"
        ),
    }
    output_path = _resolve(root, str(args.output))
    atomic_json(output_path, output)

    for record, training_export, heldout_export in loaded:
        train_labels = {
            name: labels_for_spots(
                catalog,
                name,
                training_export["section_ids"],
                training_export["barcodes"],
            )
            for name in STRUCTURES
        }
        entry_output = {
            "restart_index": record["restart_index"],
            "fold": record["fold"],
            "coverage_training": _coverage(train_labels, np.asarray(training_export["patient_ids"])),
            "full_k3": _readout_one_rank(
                training_export, heldout_export, train_labels, max_rank=None
            ),
            "stable_rank2_sensitivity": _readout_one_rank(
                training_export, heldout_export, train_labels, max_rank=2
            ),
            "outer_heldout_gt": "SEALED_NOT_READ" if heldout_export is not None else "NOT_IN_SCOPE",
        }
        output["entries"].append(entry_output)
        atomic_json(output_path, output)
    output["status"] = "COMPLETE_DESCRIPTIVE_NOT_FORMAL_STRUCTURE_CLAIM"
    atomic_json(output_path, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
