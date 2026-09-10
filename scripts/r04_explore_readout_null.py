#!/usr/bin/env python3
"""Exploratory label-permutation null for the training-role structure readout.

Reads frozen Torch training-role field exports and training-role GT only.
For each outer fold and representation rank, the observed inner-crossfit
specific increment is compared against a null distribution obtained by
permuting structure labels jointly within patient (preserves per-patient
prevalence and inter-structure correlation, breaks spot-label association).

Exploratory grade: the statistic is evaluated post hoc; no confirmatory
claim may be based on these p-values.  Observed values are cross-checked
against the frozen readout panel and fail closed on mismatch.
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
from scripts.r04_run_structure_readout import (
    _balanced_patient_weights,
    _inner_patient_folds,
)

STRUCTURES = ("TLS", "TUMOR_STROMA_BOUNDARY")


def _read_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _paired_subset(
    train_z: np.ndarray,
    labels: dict[str, np.ndarray],
    patient_ids: np.ndarray,
) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]:
    mask = np.isfinite(np.vstack([labels[name] for name in STRUCTURES])).all(axis=0)
    if not np.any(mask):
        raise ValueError("no paired annotated spots")
    return (
        np.asarray(train_z, dtype=float)[mask],
        {name: np.asarray(labels[name], dtype=float)[mask] for name in STRUCTURES},
        np.asarray(patient_ids).astype(str)[mask],
    )


def _has_target_variation(
    paired_labels: dict[str, np.ndarray], train: np.ndarray, test: np.ndarray
) -> bool:
    return all(
        len(np.unique(paired_labels[name][train])) >= 2
        and len(np.unique(paired_labels[name][test])) >= 2
        for name in STRUCTURES
    )


def aggregate_specific_increment(
    paired_z: np.ndarray,
    paired_labels: dict[str, np.ndarray],
    paired_patients: np.ndarray,
    representation: dict[str, object] | None = None,
) -> dict[str, object]:
    """Mean inner-crossfit specific increment, same aggregation as the panel."""
    folds = _inner_patient_folds(paired_patients)
    auc_values: dict[str, list[float]] = {name: [] for name in STRUCTURES}
    mse_values: list[float] = []
    valid = 0
    for inner_train, inner_test, _ in folds:
        if not _has_target_variation(paired_labels, inner_train, inner_test):
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
        for name in STRUCTURES:
            value = readout["specific_increment"]["auc_delta_by_target"][name]
            if value is not None:
                auc_values[name].append(float(value))
        mse_values.append(
            float(readout["specific_increment"]["mean_mse_delta_specific_minus_shared"])
        )
        valid += 1
    return {
        "valid_inner_folds": valid,
        "total_inner_folds": len(folds),
        "mean_auc_delta_by_target": {
            name: (float(np.mean(values)) if values else None)
            for name, values in auc_values.items()
        },
        "mean_mse_delta": float(np.mean(mse_values)) if mse_values else None,
    }


def permutation_null_specific_increment(
    paired_z: np.ndarray,
    paired_labels: dict[str, np.ndarray],
    paired_patients: np.ndarray,
    *,
    seed: int,
    draws: int,
    representation: dict[str, object] | None = None,
) -> dict[str, object]:
    """Null distribution of the aggregate increment under within-patient label permutation."""
    if draws < 1:
        raise ValueError("draws must be positive")
    rng = np.random.default_rng(seed)
    patients = np.asarray(paired_patients).astype(str)
    groups = [np.flatnonzero(patients == value) for value in np.unique(patients)]
    null_auc: dict[str, list[float]] = {name: [] for name in STRUCTURES}
    null_mse: list[float] = []
    valid_draws = 0
    for _ in range(draws):
        permutation = np.arange(len(patients))
        for indices in groups:
            permutation[indices] = rng.permutation(indices)
        permuted = {
            name: np.asarray(values, dtype=float)[permutation]
            for name, values in paired_labels.items()
        }
        aggregate = aggregate_specific_increment(
            paired_z, permuted, patients, representation=representation
        )
        if aggregate["valid_inner_folds"] == 0:
            continue
        valid_draws += 1
        for name in STRUCTURES:
            value = aggregate["mean_auc_delta_by_target"][name]
            if value is not None:
                null_auc[name].append(float(value))
        if aggregate["mean_mse_delta"] is not None:
            null_mse.append(float(aggregate["mean_mse_delta"]))
    return {
        "draws": draws,
        "valid_draws": valid_draws,
        "null_auc_delta_by_target": {name: null_auc[name] for name in STRUCTURES},
        "null_mse_delta": null_mse,
    }


def empirical_upper_p(null_values: list[float], observed: float | None) -> float | None:
    if observed is None or not null_values:
        return None
    return float(np.mean(np.asarray(null_values) >= observed))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--effect-panel", type=Path, required=True)
    parser.add_argument("--readout-panel", type=Path, required=True,
                        help="frozen readout panel for observed-value cross-check")
    parser.add_argument("--registry", type=Path,
                        default=Path("infra/structure-registry/structure_instances.tsv"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260905)
    parser.add_argument("--draws", type=int, default=200)
    parser.add_argument("--ranks", default="full,2")
    args = parser.parse_args()
    root = args.project_root.resolve()
    panel = _read_json(root / args.effect_panel)
    if panel.get("schema") != "r04.frozen_effect_export_panel.v1":
        raise SystemExit("effect panel schema mismatch")
    readout_panel = _read_json(root / args.readout_panel)
    panel_by_fold = {(e["restart_index"], e["fold"]): e for e in readout_panel.get("entries", [])}
    catalog = load_training_gt_catalog(
        root / args.registry, root / str(panel["training_manifest"]["path"])
        if isinstance(panel.get("training_manifest"), dict)
        else root / str(panel["training_manifest"]),
    )
    rank_specs: list[tuple[str, int | None]] = []
    for token in str(args.ranks).split(","):
        token = token.strip()
        if token == "full":
            rank_specs.append(("full", None))
        elif token == "2":
            rank_specs.append(("rank2", 2))
        else:
            raise SystemExit(f"invalid rank spec: {token}")
    output: dict[str, object] = {
        "schema": "r04.explore_readout_null.v1",
        "status": "EXPLORATORY_LABEL_PERMUTATION_NULL",
        "evidence_grade": "exploratory_post_hoc_no_confirmatory_claim",
        "seed": args.seed,
        "draws": args.draws,
        "permutation_unit": "spots_within_patient_labels_permuted_jointly",
        "effect_panel": {"path": str(args.effect_panel)},
        "readout_panel": {"path": str(args.readout_panel)},
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
        key = (record["restart_index"], record["fold"])
        panel_entry = panel_by_fold.get(key)
        if panel_entry is None:
            raise SystemExit(f"readout panel has no entry for {key}")
        train_labels = {
            name: labels_for_spots(
                catalog, name, export["section_ids"], export["barcodes"]
            )
            for name in STRUCTURES
        }
        for rank_name, max_rank in rank_specs:
            train_z, _, representation = effect_coordinates(
                export, export, max_rank=max_rank
            )
            paired_z, paired_labels, paired_patients = _paired_subset(
                train_z, train_labels, np.asarray(export["patient_ids"])
            )
            observed = aggregate_specific_increment(
                paired_z, paired_labels, paired_patients,
                representation=representation,
            )
            expected = panel_entry["full_k3" if rank_name == "full" else "stable_rank2_sensitivity"]["inner_crossfit"]
            for name in STRUCTURES:
                got = observed["mean_auc_delta_by_target"][name]
                want = (expected.get("mean_auc_delta_by_target") or {}).get(name)
                if got is None or want is None or not np.isclose(got, want, rtol=1e-9, atol=1e-12):
                    raise SystemExit(
                        f"observed cross-check failed for fold {record['fold']} {rank_name} {name}: "
                        f"{got} != panel {want}"
                    )
            null = permutation_null_specific_increment(
                paired_z, paired_labels, paired_patients,
                seed=args.seed + 1000 * int(record["fold"]) + (0 if rank_name == "full" else 77),
                draws=args.draws, representation=representation,
            )
            if null["valid_draws"] < args.draws // 2:
                raise SystemExit("too few valid null draws")
            summary = {
                "restart_index": record["restart_index"],
                "fold": record["fold"],
                "rank": rank_name,
                "n_paired_spots": int(len(paired_patients)),
                "n_paired_patients": sorted(set(paired_patients.tolist())),
                "observed": observed,
                "null_quantiles_auc": {
                    name: {
                        q: float(np.quantile(values, float(q)))
                        for q in ("0.025", "0.5", "0.975")
                    } if values else None
                    for name, values in null["null_auc_delta_by_target"].items()
                },
                "null_quantiles_mse": {
                    q: float(np.quantile(null["null_mse_delta"], float(q)))
                    for q in ("0.025", "0.5", "0.975")
                } if null["null_mse_delta"] else None,
                "empirical_upper_p_auc": {
                    name: empirical_upper_p(
                        null["null_auc_delta_by_target"][name],
                        observed["mean_auc_delta_by_target"][name],
                    )
                    for name in STRUCTURES
                },
                "empirical_lower_p_auc": {
                    name: (
                        float(np.mean(np.asarray(null["null_auc_delta_by_target"][name])
                                      <= observed["mean_auc_delta_by_target"][name]))
                        if null["null_auc_delta_by_target"][name]
                        and observed["mean_auc_delta_by_target"][name] is not None
                        else None
                    )
                    for name in STRUCTURES
                },
                "valid_null_draws": null["valid_draws"],
                "null_draws_requested": null["draws"],
            }
            output["entries"].append(summary)
            atomic_json(output_path, output)
    output["status"] = "EXPLORATORY_COMPLETE"
    atomic_json(output_path, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
