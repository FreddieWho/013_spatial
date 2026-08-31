#!/usr/bin/env python3
"""Audit and aggregate the five-seed R-04 extreme-fold restart panel.

The panel is a stability diagnostic, not a K-selection procedure.  Score
differences are paired by patient.  Cross-restart factor stability is measured
on continuous loading and field subspaces, so factor permutations, signs and
rotations are treated as nuisance symmetries rather than cluster labels.
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
from r04.spatial import normalize_coordinates
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
    if cell.get("status") not in {
        "FIT_AND_SCORED",
        "FIT_AND_SCORED_DUAL_GATE_DIAGNOSTIC",
    }:
        raise ValueError("restart panel contains an incomplete cell")
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


def aggregate_restart_scores(
    entries: list[dict[str, object]],
    *,
    allowed_protocol_exceptions: list[dict[str, object]] | None = None,
) -> dict[str, object]:
    """Return paired score and platform summaries for the fixed 5x2 panel."""
    cells: dict[tuple[int, int, int], Mapping[str, object]] = {}
    protocols: list[tuple[tuple[int, int, int], dict[str, object]]] = []
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
            protocols.append((key, _protocol(cell)))

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
    return {
        "restart_indices": list(range(5)),
        "fold_values": [0, 4],
        "k_values": [0, 3],
        "cell_count": len(all_cells),
        "protocol": reference_protocol,
        "protocol_exceptions": applied_exceptions,
        "input_hash_by_fold": {str(key): value for key, value in sorted(input_hash_by_fold.items())},
        "pairs": sorted(pair_results, key=lambda item: (int(item["restart_index"]), int(item["fold"]))),
        "folds": folds,
        "platform_counts": platform_counts,
        "k3_collapse_warning_cells": collapse_warning_cells,
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
    import tensorflow as tf

    checkpoint = tf.train.latest_checkpoint(str(checkpoint_dir))
    if checkpoint is None:
        raise ValueError(f"checkpoint is missing: {checkpoint_dir}")
    variables = dict(tf.train.list_variables(checkpoint))
    raw_w_name = "raw_w/.ATTRIBUTES/VARIABLE_VALUE"
    if raw_w_name not in variables:
        raise ValueError(f"raw_w is missing from checkpoint: {checkpoint}")
    reader = tf.train.load_checkpoint(checkpoint)
    raw_w = np.asarray(reader.get_tensor(raw_w_name), dtype=float)
    shifted = raw_w - raw_w.max(axis=0, keepdims=True)
    loading = np.exp(shifted)
    loading /= loading.sum(axis=0, keepdims=True)

    fold_ids = grouped_fold_ids(
        [str(row["patient_id"]) for row in training_rows], n_folds=5, seed=group_seed
    )
    rows = [row for row, fold_id in zip(training_rows, fold_ids) if int(fold_id) != fold]
    q_loc_names = sorted(
        (
            (int(match.group(1)), name)
            for name in variables
            if (match := re.fullmatch(r"q_loc_(\d+)/\.ATTRIBUTES/VARIABLE_VALUE", name))
        ),
        key=lambda item: item[0],
    )
    if len(q_loc_names) != len(rows):
        raise ValueError(
            f"checkpoint q_loc count {len(q_loc_names)} does not match training sections {len(rows)}"
        )
    field_parts = []
    for row, (index, name) in zip(rows, q_loc_names):
        if index != len(field_parts):
            raise ValueError("checkpoint q_loc indices are not contiguous")
        metadata_path = Path(str(row["matrix_metadata_locator"]))
        metadata = _read_json(metadata_path)
        coords, _ = normalize_coordinates(np.asarray(metadata["coords"], dtype=float))
        number = min(inducing_points, len(coords))
        inducing = _inducing(coords, number)
        basis, _, _ = _inducing_projection(coords, inducing, lengthscale, ridge)
        q_loc = np.asarray(reader.get_tensor(name), dtype=float)
        if q_loc.shape != (number, loading.shape[1]):
            raise ValueError(f"q_loc shape mismatch at training section {index}")
        field = basis @ q_loc
        field_parts.append(field - field.mean(axis=0, keepdims=True))
    return loading, np.concatenate(field_parts, axis=0), checkpoint


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reuse-subspace-from", type=Path, default=None)
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    manifest = _read_json(args.manifest)
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
            subspace_by_fold[str(fold)] = {
                "loading_subspace": pairwise_subspace_stability(loading_by_restart),
                "training_field_subspace": pairwise_subspace_stability(field_by_restart),
                "checkpoint_by_restart": checkpoints,
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
        "created_at": "2026-08-31",
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
