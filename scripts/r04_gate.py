#!/usr/bin/env python3
"""Evaluate one frozen candidate from precomputed grouped metrics."""

from __future__ import annotations

import argparse
from pathlib import Path

from r04.candidates import evaluate_candidate_gate
from r04.runtime import atomic_json
from r04.serialization import read_candidate_registry, read_json
from r04.types import GateResult


def _write_result(path: Path, result: GateResult) -> int:
    atomic_json(path, {
        "schema": "r04.gate.v1", "status": result.status,
        "reasons": list(result.reasons), "metrics": dict(result.metrics),
    })
    return 0 if result.status == "PASS_R04_STABLE_FIELD" else 2


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--candidate-id", required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        registry_hash, candidates = read_candidate_registry(args.registry)
    except (OSError, ValueError, TypeError) as exc:
        result = GateResult("BLOCKED_CANDIDATE_REGISTRY_INVALID", (str(exc),))
        return _write_result(args.output, result)
    candidate = next((value for value in candidates if value.candidate_id == args.candidate_id), None)
    if candidate is None:
        result = GateResult("BLOCKED_CANDIDATE_NOT_IN_REGISTRY", ("candidate_id_not_in_frozen_registry",))
        return _write_result(args.output, result)
    try:
        metrics = read_json(args.metrics)
    except (OSError, ValueError, TypeError) as exc:
        result = GateResult("BLOCKED_GATE_METRICS_INVALID", (str(exc),))
        return _write_result(args.output, result)
    if not isinstance(metrics, dict):
        result = GateResult("BLOCKED_GATE_METRICS_INVALID", ("metrics must be a JSON object",))
        return _write_result(args.output, result)
    required_provenance = {"candidate_registry_hash", "input_hash", "field_hash", "uncertainty_complete", "lengthscale_supported"}
    required_metrics = {"heldout_delta", "heldout_ci_low", "spatial_variance_probability", "null_fdr", "lengthscale", "max_patient_weight", "independent_lineages", "groups_per_lineage"}
    missing = sorted((required_provenance | required_metrics) - set(metrics))
    if missing:
        result = GateResult("BLOCKED_GATE_PROVENANCE", ("required_gate_fields_missing:" + ",".join(missing),))
    elif not candidate.input_hash or not candidate.field_hash:
        result = GateResult("BLOCKED_FIELD_INPUT_HASH_MISMATCH", ("candidate_field_or_input_hash_missing",))
    elif metrics["candidate_registry_hash"] != registry_hash:
        result = GateResult("BLOCKED_CANDIDATE_HASH_MISMATCH", ("candidate_registry_hash_mismatch",))
    elif metrics["input_hash"] != candidate.input_hash or metrics["field_hash"] != candidate.field_hash:
        result = GateResult("BLOCKED_FIELD_INPUT_HASH_MISMATCH", ("field_or_input_hash_mismatch",))
    else:
        gate_keys = required_metrics | {"single_model", "uncertainty_complete", "lengthscale_supported"}
        try:
            result = evaluate_candidate_gate(candidate, **{key: metrics[key] for key in gate_keys if key in metrics})
        except (TypeError, ValueError) as exc:
            result = GateResult("BLOCKED_GATE_METRICS_INVALID", (str(exc),))
    return _write_result(args.output, result)


if __name__ == "__main__":
    raise SystemExit(main())
