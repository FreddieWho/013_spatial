#!/usr/bin/env python3
"""Exploratory deep dive of the training-role structure readout.

For each outer fold and representation rank (1/2/3), records absolute
shared-only and shared-plus-specific metrics per structure, residual
coefficient norms, singular values, and per-inner-fold tables.  Reads frozen
Torch training-role field exports and training-role GT only.

Exploratory grade: descriptive decomposition of already-computed readouts;
no confirmatory claim may be based on these values.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from r04.field_exports import read_heldout_field_export
from r04.runtime import atomic_json
from r04.structure_gt import labels_for_spots, load_training_gt_catalog
from r04.structure_readout import effect_coordinates, evaluate_shared_specific_readout
from scripts.r04_explore_readout_null import (
    STRUCTURES,
    _balanced_patient_weights,
    _has_target_variation,
    _inner_patient_folds,
    _paired_subset,
    _read_json,
)


def rank_readout(
    paired_z: np.ndarray,
    paired_labels: dict[str, np.ndarray],
    paired_patients: np.ndarray,
    representation: dict[str, object] | None = None,
) -> dict[str, object]:
    """Per-inner-fold readout records plus aggregate absolute metrics."""
    folds = _inner_patient_folds(paired_patients)
    fold_records = []
    for inner_index, (inner_train, inner_test, test_patients) in enumerate(folds):
        record: dict[str, object] = {
            "inner_fold": inner_index,
            "test_patients": test_patients,
            "training_spots": int(np.sum(inner_train)),
            "heldout_spots": int(np.sum(inner_test)),
        }
        if not _has_target_variation(paired_labels, inner_train, inner_test):
            record["status"] = "NOT_TESTABLE_NO_TARGET_VARIATION"
            fold_records.append(record)
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
        record.update({"status": "COMPUTED", "readout": readout})
        fold_records.append(record)
    valid = [item for item in fold_records if item.get("status") == "COMPUTED"]

    def mean_metric(section: str, metric: str, structure: str) -> float | None:
        values = [
            item["readout"][section]["metrics"][structure][metric]
            for item in valid
            if item["readout"][section]["metrics"][structure][metric] is not None
        ]
        return float(np.mean(values)) if values else None

    def mean_norm(structure: str) -> float | None:
        values = [
            float(item["readout"]["shared_plus_specific_residual"]
                  ["specific_residual_coefficient_norm"][structure])
            for item in valid
        ]
        return float(np.mean(values)) if values else None

    return {
        "valid_inner_folds": len(valid),
        "total_inner_folds": len(folds),
        "absolute_metrics": {
            structure: {
                "shared_auc": mean_metric("shared_only", "auc", structure),
                "specific_auc": mean_metric("shared_plus_specific_residual", "auc", structure),
                "shared_mse": mean_metric("shared_only", "mse", structure),
                "specific_mse": mean_metric("shared_plus_specific_residual", "mse", structure),
                "shared_correlation": mean_metric("shared_only", "correlation", structure),
                "specific_correlation": mean_metric(
                    "shared_plus_specific_residual", "correlation", structure),
                "specific_residual_coefficient_norm": mean_norm(structure),
            }
            for structure in STRUCTURES
        },
        "fold_records": fold_records,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--effect-panel", type=Path, required=True)
    parser.add_argument("--registry", type=Path,
                        default=Path("infra/structure-registry/structure_instances.tsv"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ranks", default="1,2,3")
    args = parser.parse_args()
    root = args.project_root.resolve()
    panel = _read_json(root / args.effect_panel)
    if panel.get("schema") != "r04.frozen_effect_export_panel.v1":
        raise SystemExit("effect panel schema mismatch")
    catalog = load_training_gt_catalog(
        root / args.registry, root / str(panel["training_manifest"]["path"])
        if isinstance(panel.get("training_manifest"), dict)
        else root / str(panel["training_manifest"]),
    )
    ranks = []
    for token in str(args.ranks).split(","):
        token = token.strip()
        if not token.isdigit() or int(token) < 1:
            raise SystemExit(f"invalid rank: {token}")
        ranks.append(int(token))
    output: dict[str, object] = {
        "schema": "r04.explore_readout_deepdive.v1",
        "status": "EXPLORATORY_RANK_CURVE",
        "evidence_grade": "exploratory_post_hoc_no_confirmatory_claim",
        "effect_panel": {"path": str(args.effect_panel)},
        "structures": list(STRUCTURES),
        "entries": [],
    }
    output_path = root / args.output
    atomic_json(output_path, output)
    for record in panel.get("entries", []):
        exports = record.get("structure_exports", [])
        if record.get("status") != "TRAINING_EFFECT_EXPORTED" or len(exports) != 1:
            raise SystemExit("effect panel entry is not a single training export")
        export = read_heldout_field_export(root / str(exports[0]["npz_path"]))
        train_labels = {
            name: labels_for_spots(
                catalog, name, export["section_ids"], export["barcodes"]
            )
            for name in STRUCTURES
        }
        entry: dict[str, object] = {
            "restart_index": record["restart_index"],
            "fold": record["fold"],
            "ranks": {},
        }
        for rank in ranks:
            train_z, _, representation = effect_coordinates(
                export, export, max_rank=rank
            )
            paired_z, paired_labels, paired_patients = _paired_subset(
                train_z, train_labels, np.asarray(export["patient_ids"])
            )
            readout = rank_readout(
                paired_z, paired_labels, paired_patients,
                representation=representation,
            )
            readout["representation"] = representation
            readout["n_paired_spots"] = int(len(paired_patients))
            readout["n_paired_patients"] = sorted(set(paired_patients.tolist()))
            entry["ranks"][f"rank{rank}"] = readout
        output["entries"].append(entry)
        atomic_json(output_path, output)
    output["status"] = "EXPLORATORY_COMPLETE"
    atomic_json(output_path, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
