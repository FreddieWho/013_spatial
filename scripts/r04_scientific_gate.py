#!/usr/bin/env python3
"""Summarise objective and synthetic calibration gates without real-data claims."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from r04.runtime import atomic_json
from scripts.r04_k0_gate import evaluate_no_field_shards
from scripts.r04_k_aggregate import aggregate_k_shards
from scripts.r04_k_contract import EXECUTION_CONTRACT, GENERATOR_CONTRACT, LEAKAGE_CONTROL


def _read(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _provenance(path: Path, value: dict[str, object]) -> dict[str, object]:
    payload = path.read_bytes()
    return {
        "path": str(path),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "bytes": len(payload),
        "schema": value.get("schema"),
    }


def _input_validation_reasons(
    objective: dict[str, object],
    calibrations: list[dict[str, object]],
    k_calibration: dict[str, object] | None,
    k0_gate: dict[str, object] | None,
) -> list[str]:
    reasons: list[str] = []
    if objective.get("schema") != "r04.objective_oracle.v1":
        reasons.append("objective_schema_invalid")
    if objective.get("fit_infer_dense_scaling") != "per_spot_gene_sum":
        reasons.append("objective_fit_infer_binding_invalid")
    required_objective_fields = {
        "dense_objective",
        "full_batch_objective",
        "mean_gene_batch_objective",
        "objective_abs_error",
        "mean_batch_gradient_max_abs_error",
    }
    if not required_objective_fields.issubset(objective):
        reasons.append("objective_gradient_oracle_missing")
    else:
        try:
            objective_values = [float(objective[key]) for key in required_objective_fields]
            if (
                not all(math.isfinite(value) for value in objective_values)
                or objective.get("gradient_finite") is not True
                or int(objective.get("batch_count", 0)) != 20
                or float(objective["objective_abs_error"]) > 1e-5
                or float(objective["mean_batch_gradient_max_abs_error"]) > 1e-5
                or abs(float(objective["dense_objective"]) - float(objective["full_batch_objective"])) > 1e-6
            ):
                reasons.append("objective_or_gradient_oracle_tolerance_failed")
        except (TypeError, ValueError, OverflowError):
            reasons.append("objective_or_gradient_oracle_tolerance_failed")
    if any(item.get("schema") != "r04.identifiability_calibration.v2" for item in calibrations):
        reasons.append("calibration_schema_invalid")
    expected_seeds = {20260807, 20260817, 20260827}
    expected_identities = {
        ("two_independent", domain, seed)
        for domain in ("disconnected", "crescent", "branch")
        for seed in expected_seeds
    } | {
        (mode, "disconnected", seed)
        for mode in ("no_field", "one_disconnected", "two_collinear")
        for seed in expected_seeds
    }
    observed_identities: set[tuple[str, str, int]] = set()
    for item in calibrations:
        truth = item.get("truth", {})
        if not isinstance(truth, dict):
            reasons.append("calibration_truth_invalid")
            continue
        try:
            identity = (str(truth["mode"]), str(truth["domain"]), int(item["seed"]))
        except (KeyError, TypeError, ValueError):
            reasons.append("calibration_identity_missing")
            continue
        observed_identities.add(identity)
        if item.get("model") != "mnsf" or truth.get("generator") != "mnsf_correct":
            reasons.append("calibration_model_generator_binding_invalid")
        expected_config = {
            "inducing_points": 16,
            "nonspatial_rank": 0,
            "lengthscale": 3.0,
            "learning_rate": 0.05,
            "null_draws": 200,
            "platform_window": 50,
        }
        try:
            calibration_steps = int(item.get("steps", 0))
        except (TypeError, ValueError, OverflowError):
            calibration_steps = -1
        if (
            item.get("status") != "CALIBRATION_COMPLETE_NOT_VALIDATED"
            or item.get("n_factors_model") != 2
            or calibration_steps < 300
            or item.get("config") != expected_config
        ):
            reasons.append("calibration_execution_contract_invalid")
        if truth.get("generator_version") not in {
            "independent_gene_effects_v2", "crescent_supported_fields_v3"
        }:
            reasons.append("calibration_generator_version_invalid")
    if observed_identities != expected_identities or len(calibrations) != len(expected_identities):
        reasons.append("calibration_matrix_incomplete_or_duplicated")
    if k_calibration is None or k_calibration.get("schema") != "r04.k_calibration_aggregate.v1":
        reasons.append("k_calibration_schema_invalid")
    elif (
        k_calibration.get("k_values") != [1, 2, 3, 4]
        or k_calibration.get("k_eff_calibration_truth") != 2
        or len(k_calibration.get("data_seeds", [])) != 3
        or k_calibration.get("calibration_cell_count") != 24
        or k_calibration.get("independent_resampling_unit") != "data_seed"
        or len(k_calibration.get("input_data_hash_by_seed", {})) != 3
        or k_calibration.get("generator_replay_verified") is not True
        or k_calibration.get("execution_contract") != EXECUTION_CONTRACT
        or k_calibration.get("generator_contract") != GENERATOR_CONTRACT
        or k_calibration.get("leakage_control") != LEAKAGE_CONTROL
    ):
        reasons.append("k_calibration_binding_invalid")
    elif (
        len(k_calibration.get("input_provenance", [])) < 12
        or len({item.get("sha256") for item in k_calibration.get("input_provenance", [])})
        != len(k_calibration.get("input_provenance", []))
    ):
        reasons.append("k_calibration_provenance_invalid")
    if k0_gate is None or k0_gate.get("schema") != "r04.k0_gate.v1":
        reasons.append("k0_gate_schema_invalid")
    elif (
        k0_gate.get("candidate_k") != [0, 1, 2]
        or len(k0_gate.get("seeds", [])) != 3
        or k0_gate.get("calibration_cell_count") != 12
        or k0_gate.get("independent_resampling_unit") != "data_seed"
        or len(k0_gate.get("input_data_hash_by_seed", {})) != 3
        or k0_gate.get("generator_replay_verified") is not True
        or k0_gate.get("execution_contract") != EXECUTION_CONTRACT
        or k0_gate.get("generator_contract") != GENERATOR_CONTRACT
        or k0_gate.get("leakage_control") != LEAKAGE_CONTROL
    ):
        reasons.append("k0_gate_binding_invalid")
    elif (
        len(k0_gate.get("input_provenance", [])) < 6
        or len({item.get("sha256") for item in k0_gate.get("input_provenance", [])})
        != len(k0_gate.get("input_provenance", []))
    ):
        reasons.append("k0_gate_provenance_invalid")
    return sorted(set(reasons))


def _nested_provenance_reasons(
    aggregate: dict[str, object] | None,
    *,
    label: str,
) -> list[str]:
    if aggregate is None:
        return [f"{label}_provenance_missing"]
    items = aggregate.get("input_provenance", [])
    if not isinstance(items, list):
        return [f"{label}_provenance_invalid"]
    reasons: list[str] = []
    for item in items:
        if not isinstance(item, dict):
            reasons.append(f"{label}_provenance_invalid")
            continue
        path = Path(str(item.get("path", "")))
        try:
            payload = path.read_bytes()
        except OSError:
            reasons.append(f"{label}_source_missing")
            continue
        if (
            hashlib.sha256(payload).hexdigest() != item.get("sha256")
            or len(payload) != item.get("bytes")
        ):
            reasons.append(f"{label}_source_hash_mismatch")
    return sorted(set(reasons))


def _recomputed_aggregate_reasons(
    aggregate: dict[str, object] | None,
    *,
    label: str,
) -> list[str]:
    """Re-run aggregate validation over the exact nested source artifacts."""
    if aggregate is None:
        return [f"{label}_recomputation_missing"]
    items = aggregate.get("input_provenance")
    if not isinstance(items, list) or not items:
        return [f"{label}_recomputation_missing"]
    try:
        shards = [
            json.loads(Path(str(item["path"])).read_text(encoding="utf-8"))
            for item in items
        ]
        if not all(isinstance(shard, dict) for shard in shards):
            raise ValueError("nested calibration source is not a JSON object")
        recomputed = (
            aggregate_k_shards(shards)
            if label == "k_calibration"
            else evaluate_no_field_shards(shards)
        )
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return [f"{label}_source_semantics_invalid"]
    if any(aggregate.get(key) != value for key, value in recomputed.items()):
        return [f"{label}_recomputed_aggregate_mismatch"]
    return []


def _calibration_result_reasons(calibration: dict[str, object]) -> list[str]:
    """Validate one planted-control result, including finite subspace values."""
    reasons: list[str] = []
    platform = calibration.get("platform", {})
    subspace = calibration.get("subspace", {})
    truth = calibration.get("truth", {})
    if not isinstance(platform, dict) or platform.get("converged") is not True:
        reasons.append("platform_not_reached")
    try:
        raw_required = truth.get("k_eff_true", truth.get("k_true", 0))
        required_float = float(raw_required)
        required = int(required_float)
        if (
            isinstance(raw_required, bool)
            or not math.isfinite(required_float)
            or required_float != required
            or required < 0
        ):
            raise ValueError
    except (AttributeError, TypeError, ValueError, OverflowError):
        return [*reasons, "calibration_truth_dimension_invalid"]
    if required < 1:
        return reasons
    correlations = (
        subspace.get("canonical_correlations", [])
        if isinstance(subspace, dict)
        else []
    )
    cutoffs = (
        subspace.get("direction_matched_null_95")
        if isinstance(subspace, dict)
        else None
    )
    if cutoffs is None and isinstance(subspace, dict):
        legacy_cutoff = subspace.get("null_95_cutoff")
        cutoffs = [legacy_cutoff] * len(correlations) if legacy_cutoff is not None else None
    try:
        correlation_values = [float(value) for value in correlations]
        cutoff_values = [float(value) for value in cutoffs]
    except (TypeError, ValueError, OverflowError):
        correlation_values = []
        cutoff_values = []
    if (
        len(correlation_values) < required
        or len(cutoff_values) < required
        or not all(
            math.isfinite(correlation_values[index])
            and math.isfinite(cutoff_values[index])
            and correlation_values[index] > cutoff_values[index]
            for index in range(required)
        )
    ):
        reasons.append("planted_effect_subspace_not_recovered_above_null")
    return reasons


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--objective", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, action="append", required=True)
    parser.add_argument("--k-calibration", type=Path)
    parser.add_argument("--k0-gate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    objective = _read(args.objective)
    calibrations = [_read(path) for path in args.calibration]
    k_calibration = _read(args.k_calibration) if args.k_calibration else None
    k0_gate = _read(args.k0_gate) if args.k0_gate else None
    reasons = _input_validation_reasons(objective, calibrations, k_calibration, k0_gate)
    reasons.extend(_nested_provenance_reasons(k_calibration, label="k_calibration"))
    reasons.extend(_nested_provenance_reasons(k0_gate, label="k0_gate"))
    reasons.extend(_recomputed_aggregate_reasons(k_calibration, label="k_calibration"))
    reasons.extend(_recomputed_aggregate_reasons(k0_gate, label="k0_gate"))
    if objective.get("status") != "PASS_OBJECTIVE_ORACLE":
        reasons.append("objective_oracle_failed")
    if k_calibration is not None and k_calibration.get("final_k_status") != "K_SELECTED":
        reasons.append("k_not_identifiable")
    if k0_gate is None or k0_gate.get("status") != "K0_SELECTED":
        reasons.append("no_field_false_positive_gate_not_passed")
    if not calibrations:
        reasons.append("no_calibration_artifacts")
    summaries: list[dict[str, object]] = []
    for calibration in calibrations:
        platform = calibration.get("platform", {})
        subspace = calibration.get("subspace", {})
        truth = calibration.get("truth", {})
        reasons.extend(_calibration_result_reasons(calibration))
        summaries.append({
            "model": calibration.get("model"),
            "seed": calibration.get("seed"),
            "domain": truth.get("domain") if isinstance(truth, dict) else None,
            "mode": truth.get("mode") if isinstance(truth, dict) else None,
            "k_true": truth.get("k_true") if isinstance(truth, dict) else None,
            "k_eff_true": truth.get("k_eff_true") if isinstance(truth, dict) else None,
            "canonical_correlations": subspace.get("canonical_correlations") if isinstance(subspace, dict) else None,
            "null_95_cutoff": subspace.get("null_95_cutoff") if isinstance(subspace, dict) else None,
            "direction_matched_null_95": subspace.get("direction_matched_null_95") if isinstance(subspace, dict) else None,
            "platform_converged": platform.get("converged") if isinstance(platform, dict) else False,
        })
    status = "PASS_SCIENTIFIC_CALIBRATION" if not reasons else "BLOCKED_SCIENTIFIC_CALIBRATION"
    result = {
        "schema": "r04.scientific_gate.v2",
        "status": status,
        "reasons": sorted(set(reasons)),
        "objective_oracle": objective.get("status"),
        "input_provenance": {
            "objective": _provenance(args.objective, objective),
            "calibrations": [
                _provenance(path, value)
                for path, value in zip(args.calibration, calibrations)
            ],
            "k_calibration": _provenance(args.k_calibration, k_calibration)
            if args.k_calibration is not None and k_calibration is not None else None,
            "k0_gate": _provenance(args.k0_gate, k0_gate)
            if args.k0_gate is not None and k0_gate is not None else None,
        },
        "calibrations": summaries,
        "k_status": "K_NOT_IDENTIFIABLE" if k_calibration is None or k_calibration.get("final_k_status") != "K_SELECTED" else "SYNTHETIC_K_RULE_CALIBRATED",
        "k_calibration": {
            "status": k_calibration.get("status") if k_calibration else "NOT_RUN",
            "final_k_status": k_calibration.get("final_k_status") if k_calibration else "K_NOT_IDENTIFIABLE",
            "one_se_candidate": k_calibration.get("one_se_candidate") if k_calibration else None,
        },
        "k0_gate": {
            "status": k0_gate.get("status") if k0_gate else "NOT_RUN",
            "one_se_candidate": k0_gate.get("one_se_candidate") if k0_gate else None,
        },
        "formal_restart_authorized": False,
        "next_gate": "WIRE_AND_TEST_REAL_DATA_K0_UPWARD_SEARCH_WITH_DURABLE_CHECKPOINTS",
        "gpu_used": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(args.output, result)
    return 0 if not reasons else 2


if __name__ == "__main__":
    raise SystemExit(main())
