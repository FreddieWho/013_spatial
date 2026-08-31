#!/usr/bin/env python3
"""Run held-out inference from an already fitted R-04 checkpoint.

This is the first dual-gate diagnostic: the dense objective platform is
reported separately from the legacy minibatch fit gate.  It never selects K or
aggregates a formal campaign.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from r04.diagnostics import deterministic_objective_platform_summary
from r04.practical_acceptance import validate_practical_inference_acceptance
from r04.runtime import atomic_json, resource_status
from r04.continuation_platform_audit import (
    checkpoint_bundle_fingerprint,
    replace_final_live_evaluation,
)
from scripts.r04_real_k_search import (
    _grouped_folds,
    _load_training,
    _optimization_seed,
    _parse_gene_batch_size,
    _score_fold,
)

EXPECTED_FROZEN_INFERENCE_PROTOCOL = {
    "folds": 5,
    "group_seed": 20260807,
    "gene_folds": 2,
    "split_seed": 20260817,
    "inference_steps": 2400,
    "inducing_points": 16,
    "lengthscale": 3.0,
    "learning_rate": 0.01,
}


def _read_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _checkpoint_environment_hash(checkpoint_dir: Path) -> str:
    metadata = _read_json(checkpoint_dir / "checkpoint_metadata.json")
    value = metadata.get("environment_hash")
    if not isinstance(value, str) or not value:
        raise RuntimeError("checkpoint metadata has no environment hash")
    return value


def _objective_trace(audit: dict[str, object]) -> list[dict[str, object]]:
    diagnostics = audit.get("fit_diagnostics")
    if isinstance(diagnostics, dict):
        trace = diagnostics.get("evaluation_trace")
        if isinstance(trace, list):
            return trace
    checkpoints = audit.get("checkpoints")
    if isinstance(checkpoints, list):
        return [
            {
                "step": item["global_step"],
                "objective_mean": item["objective_mean"],
            }
            for item in checkpoints
            if isinstance(item, dict)
            and "global_step" in item
            and "objective_mean" in item
        ]
    raise ValueError("fit audit has no deterministic objective trace")


def _require_continuation_audit(
    source_audit: dict[str, object],
    audit_path: Path | None,
    fit_audit_path: Path | None = None,
    source_checkpoint_dir: Path | None = None,
    checkpoint_repeat_audit_path: Path | None = None,
    practical_acceptance_path: Path | None = None,
) -> dict[str, object] | None:
    diagnostics = source_audit.get("fit_diagnostics")
    legacy_current_panel = (
        source_audit.get("schema") == "r04.late_lr_cell.v1"
        and source_audit.get("restart_index") in {2, 3, 4}
        and isinstance(diagnostics, dict)
        and diagnostics.get("start_step") == 7200
        and diagnostics.get("steps") == 8400
    )
    declared = bool(source_audit.get("continuation_panel")) or legacy_current_panel
    if declared and (audit_path is None or fit_audit_path is None):
        raise RuntimeError("continuation-panel source cells require a panel audit")
    if audit_path is None and fit_audit_path is None:
        return None
    if audit_path is None or fit_audit_path is None:
        raise RuntimeError("continuation panel audit requires the fit audit path")
    audit = _read_json(audit_path)
    strict_all_pass = (
        audit.get("schema") == "r04.continuation_platform_audit.v2"
        and audit.get("status") == "ALL_PASS"
        and audit.get("paired_gate") is True
        and audit.get("inference_authorized") is True
    )
    practical_scope: dict[str, object] | None = None
    if strict_all_pass:
        if practical_acceptance_path is not None:
            raise RuntimeError("practical acceptance is not valid for an ALL_PASS audit")
    else:
        if practical_acceptance_path is None:
            raise RuntimeError("continuation panel audit is not ALL_PASS")
        acceptance = _read_json(practical_acceptance_path)
        try:
            practical_scope = validate_practical_inference_acceptance(
                acceptance, audit, audit_path
            )
        except ValueError as exc:
            raise RuntimeError(f"invalid practical inference acceptance: {exc}") from exc
    digest = hashlib.sha256(fit_audit_path.read_bytes()).hexdigest()
    identity = {
        "restart_index": source_audit.get("restart_index"),
        "k_model": source_audit.get("k_model"),
        "fold": source_audit.get("fold"),
    }
    matches = [
        item for item in audit.get("cells", [])
        if isinstance(item, dict)
        and item.get("identity") == identity
        and item.get("source_sha256") == digest
    ]
    if len(matches) != 1:
        raise RuntimeError("continuation panel audit has no matching source hash")
    dense = matches[0].get("dense")
    accepted_exception = (
        practical_scope is not None
        and matches[0].get("cell_id") == practical_scope.get("accepted_cell_id")
        and matches[0].get("identity") == practical_scope.get("accepted_identity")
    )
    if (
        matches[0].get("exact_repeat") is not True
        or not isinstance(dense, dict)
        or (dense.get("converged") is not True and not accepted_exception)
    ):
        raise RuntimeError("continuation panel source cell is not on platform")
    if checkpoint_repeat_audit_path is None:
        recorded_repeat_path = matches[0].get("checkpoint_repeat_audit_path")
        if not isinstance(recorded_repeat_path, str) or not recorded_repeat_path:
            raise RuntimeError("continuation panel audit has no repeat artifact path")
        checkpoint_repeat_audit_path = Path(recorded_repeat_path)
    if not checkpoint_repeat_audit_path.is_file():
        raise RuntimeError("continuation checkpoint repeat audit is missing")
    repeat_digest = hashlib.sha256(checkpoint_repeat_audit_path.read_bytes()).hexdigest()
    if repeat_digest != matches[0].get("checkpoint_repeat_audit_sha256"):
        raise RuntimeError("checkpoint repeat audit hash does not match the panel audit")
    if source_checkpoint_dir is None:
        raise RuntimeError("continuation panel audit requires the source checkpoint directory")
    expected_bundle = matches[0].get("checkpoint_bundle")
    if not isinstance(expected_bundle, dict) or not isinstance(expected_bundle.get("step"), int):
        raise RuntimeError("continuation panel audit checkpoint metadata is invalid")
    checkpoint_step = expected_bundle["step"]
    if not isinstance(diagnostics, dict) or diagnostics.get("steps") != checkpoint_step:
        raise RuntimeError("continuation panel audit checkpoint step mismatch")
    actual_bundle = checkpoint_bundle_fingerprint(source_checkpoint_dir, step=checkpoint_step)
    if not isinstance(expected_bundle, dict) or expected_bundle.get("sha256") != actual_bundle["sha256"]:
        raise RuntimeError("continuation panel audit checkpoint hash mismatch")
    result = dict(matches[0])
    result["inference_authority"] = (
        practical_scope["inference_authority"]
        if practical_scope is not None
        else "STRICT_ALL_PASS_AUDIT"
    )
    result["strict_audit_status"] = audit.get("status")
    result["strict_audit_sha256"] = hashlib.sha256(audit_path.read_bytes()).hexdigest()
    if practical_acceptance_path is not None:
        result["practical_acceptance_path"] = str(practical_acceptance_path)
        result["practical_acceptance_sha256"] = hashlib.sha256(
            practical_acceptance_path.read_bytes()
        ).hexdigest()
    return result


def _require_frozen_inference_endpoint(
    continuation_row: dict[str, object] | None,
    *,
    fit_steps: int,
    expected_start_step: int,
    inference_protocol: dict[str, object] | None = None,
) -> None:
    """Prevent a panel-authorized frozen endpoint from receiving more fit updates."""
    if continuation_row is None:
        return
    bundle = continuation_row.get("checkpoint_bundle")
    if not isinstance(bundle, dict) or not isinstance(bundle.get("step"), int):
        raise RuntimeError("continuation panel audit checkpoint metadata is invalid")
    checkpoint_step = bundle["step"]
    if fit_steps != checkpoint_step or expected_start_step != checkpoint_step:
        raise RuntimeError("continuation inference must use the audited frozen checkpoint without more fit steps")
    if inference_protocol is not None and inference_protocol != EXPECTED_FROZEN_INFERENCE_PROTOCOL:
        raise RuntimeError("continuation inference protocol does not match the frozen panel protocol")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument("--gene-list", type=Path, required=True)
    parser.add_argument("--source-checkpoint-dir", type=Path, required=True)
    parser.add_argument("--fit-audit-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--factors", type=int, required=True)
    parser.add_argument("--fold", type=int, default=0)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--group-seed", type=int, default=20260807)
    parser.add_argument("--restart-index", type=int, default=0)
    parser.add_argument("--gene-folds", type=int, default=2)
    parser.add_argument("--split-seed", type=int, default=20260817)
    parser.add_argument("--fit-steps", type=int, default=7200)
    parser.add_argument("--expected-start-step", type=int, default=7200)
    parser.add_argument("--inference-steps", type=int, default=1200)
    parser.add_argument("--checkpoint-steps", type=int, default=600)
    parser.add_argument("--inducing-points", type=int, default=16)
    parser.add_argument("--lengthscale", type=float, default=3.0)
    parser.add_argument("--learning-rate", type=float, default=0.01)
    parser.add_argument("--gene-batch-size", default="512")
    parser.add_argument("--diagnostic-interval", type=int, default=100)
    parser.add_argument("--continuation-panel-audit", type=Path, default=None)
    parser.add_argument("--practical-acceptance", type=Path, default=None)
    parser.add_argument("--checkpoint-repeat-audit-json", type=Path, default=None)
    parser.add_argument(
        "--allow-source-environment-hash",
        action="store_true",
        help="explicitly reuse the source checkpoint environment hash for a validated remote runtime",
    )
    args = parser.parse_args()
    if args.factors < 0 or args.restart_index < 0 or args.folds < 2 or args.gene_folds < 2:
        raise SystemExit("factors must be non-negative; folds and gene-folds must be at least 2")
    if args.fit_steps < args.expected_start_step or args.inference_steps < 1:
        raise SystemExit("invalid fit or inference step range")
    if not args.source_checkpoint_dir.exists():
        raise SystemExit(f"source checkpoint directory does not exist: {args.source_checkpoint_dir}")
    if not args.fit_audit_json.exists():
        raise SystemExit(f"fit audit does not exist: {args.fit_audit_json}")
    try:
        gene_batch_size = _parse_gene_batch_size(args.gene_batch_size)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    args.output_dir.mkdir(parents=True, exist_ok=True)
    if resource_status(args.output_dir) == "BLOCKED_STORAGE":
        atomic_json(args.output_dir / "r04_resume_infer.json", {"status": "BLOCKED_STORAGE"})
        return 2
    resume_environment_hash = (
        _checkpoint_environment_hash(args.source_checkpoint_dir)
        if args.allow_source_environment_hash else None
    )
    sections, genes, manifest_hash = _load_training(args.manifest_json, args.gene_list)
    fold_ids = _grouped_folds(sections, folds=args.folds, seed=args.group_seed)
    if args.fold not in set(int(value) for value in fold_ids):
        raise SystemExit(f"fold not present: {args.fold}")
    training = [
        section for index, section in enumerate(sections)
        if int(fold_ids[index]) != args.fold
    ]
    validation = [
        section for index, section in enumerate(sections)
        if int(fold_ids[index]) == args.fold
    ]
    source_audit = _read_json(args.fit_audit_json)
    continuation_row = _require_continuation_audit(
        source_audit,
        args.continuation_panel_audit,
        args.fit_audit_json,
        args.source_checkpoint_dir,
        args.checkpoint_repeat_audit_json,
        args.practical_acceptance,
    )
    _require_frozen_inference_endpoint(
        continuation_row,
        fit_steps=args.fit_steps,
        expected_start_step=args.expected_start_step,
        inference_protocol={
            "folds": args.folds,
            "group_seed": args.group_seed,
            "gene_folds": args.gene_folds,
            "split_seed": args.split_seed,
            "inference_steps": args.inference_steps,
            "inducing_points": args.inducing_points,
            "lengthscale": args.lengthscale,
            "learning_rate": args.learning_rate,
        },
    )
    seed = _optimization_seed(
        args.group_seed, args.factors, args.fold, args.restart_index
    )
    if continuation_row is not None:
        expected_identity = {
            "restart_index": args.restart_index,
            "k_model": args.factors,
            "fold": args.fold,
        }
        if continuation_row.get("identity") != expected_identity:
            raise RuntimeError("continuation inference arguments do not match the audited cell")
        if continuation_row.get("optimization_seed") != seed:
            raise RuntimeError("continuation inference seed does not match the audited cell")
        fit_objective_platform = dict(continuation_row["dense"])
    else:
        objective_trace = _objective_trace(source_audit)
        repeat_audit = None
        if args.checkpoint_repeat_audit_json is not None:
            repeat_audit = _read_json(args.checkpoint_repeat_audit_json)
        elif isinstance(source_audit.get("checkpoint_repeat_audit"), dict):
            repeat_audit = source_audit["checkpoint_repeat_audit"]
        if repeat_audit is not None:
            objective_trace = replace_final_live_evaluation(objective_trace, repeat_audit)
        fit_objective_platform = deterministic_objective_platform_summary(
            objective_trace
        )
    cell = _score_fold(
        training,
        validation,
        genes=genes,
        k_model=args.factors,
        fold=args.fold,
        output_dir=args.output_dir,
        steps=args.fit_steps,
        inference_steps=args.inference_steps,
        inducing_points=args.inducing_points,
        lengthscale=args.lengthscale,
        learning_rate=args.learning_rate,
        learning_rate_final=None,
        learning_rate_decay_steps=None,
        gene_batch_size=gene_batch_size,
        diagnostic_interval=args.diagnostic_interval,
        evaluation_interval=0,
        evaluation_mc_draws=4,
        optimization_schedule="joint",
        shared_steps=0,
        gene_folds=args.gene_folds,
        split_seed=args.split_seed,
        optimization_seed=seed,
        checkpoint_steps=args.checkpoint_steps,
        resume_checkpoint_dir=args.source_checkpoint_dir,
        environment_hash=resume_environment_hash,
    )
    diagnostics = cell.get("fit_diagnostics")
    if not isinstance(diagnostics, dict) or diagnostics.get("start_step") != args.expected_start_step:
        raise RuntimeError("resumed fit did not start at the expected checkpoint step")
    if continuation_row is not None and diagnostics.get("steps") != diagnostics.get("start_step"):
        raise RuntimeError("frozen inference unexpectedly performed fit updates")
    practical = (
        continuation_row is not None
        and continuation_row.get("inference_authority") == "PRACTICAL_POSTHOC_ACCEPTANCE"
    )
    frozen_inference_provenance = None
    if continuation_row is not None:
        checkpoint_bundle = continuation_row.get("checkpoint_bundle")
        if not isinstance(checkpoint_bundle, dict):
            raise RuntimeError("continuation checkpoint provenance is missing")
        frozen_inference_provenance = {
            "identity": {
                "restart_index": args.restart_index,
                "k_model": args.factors,
                "fold": args.fold,
            },
            "source_sha256": continuation_row.get("source_sha256"),
            "checkpoint_bundle_sha256": checkpoint_bundle.get("sha256"),
            "strict_audit_sha256": continuation_row.get("strict_audit_sha256"),
            "practical_acceptance_sha256": continuation_row.get(
                "practical_acceptance_sha256"
            ),
            "protocol": dict(EXPECTED_FROZEN_INFERENCE_PROTOCOL),
        }
    cell["status"] = (
        "FROZEN_FIT_AND_HELDOUT_SCORED_POSTHOC_ACCEPTANCE"
        if practical else "FIT_AND_SCORED_DUAL_GATE_DIAGNOSTIC"
    )
    cell["fit_objective_platform"] = fit_objective_platform
    cell["legacy_minibatch_fit_gate"] = cell.get("fit_platform")
    cell["fit_converged"] = bool(fit_objective_platform["converged"])
    cell["manifest_hash"] = manifest_hash
    cell["source_checkpoint_dir"] = str(args.source_checkpoint_dir)
    cell["fit_audit_json"] = str(args.fit_audit_json)
    cell["inference_steps_inherited_fit"] = args.inference_steps
    cell["restart_index"] = args.restart_index
    cell["optimization_seed"] = seed
    cell["environment_hash_compatibility"] = {
        "mode": (
            "SOURCE_CHECKPOINT_HASH_EXPLICIT_OVERRIDE"
            if resume_environment_hash is not None else "RUNTIME_HASH"
        ),
        "source_checkpoint_environment_hash": resume_environment_hash,
    }
    if continuation_row is not None:
        cell["inference_authority"] = continuation_row.get("inference_authority")
        cell["strict_audit_status"] = continuation_row.get("strict_audit_status")
        cell["selected_k"] = None
        cell["frozen_inference_provenance"] = frozen_inference_provenance
    atomic_json(args.output_dir / "cell.json", cell)
    atomic_json(args.output_dir / "r04_resume_infer.json", {
        "schema": "r04.resume_infer.v1",
        "status": (
            "FROZEN_HELDOUT_DIAGNOSTIC_NOT_K_SELECTION"
            if practical else "DUAL_GATE_DIAGNOSTIC_NOT_SCIENTIFIC"
        ),
        "k_model": args.factors,
        "fold": args.fold,
        "fit_objective_platform": fit_objective_platform,
        "legacy_minibatch_fit_gate": cell.get("fit_platform"),
        "fit_converged": bool(fit_objective_platform["converged"]),
        "inference_platform": cell.get("inference_platform"),
        "source_checkpoint_dir": str(args.source_checkpoint_dir),
        "fit_audit_json": str(args.fit_audit_json),
        "continuation_panel_audit": str(args.continuation_panel_audit) if continuation_row else None,
        "practical_acceptance": str(args.practical_acceptance) if practical else None,
        "inference_authority": (
            continuation_row.get("inference_authority") if continuation_row else None
        ),
        "strict_audit_status": (
            continuation_row.get("strict_audit_status") if continuation_row else None
        ),
        "frozen_inference_provenance": frozen_inference_provenance,
        "restart_index": args.restart_index,
        "optimization_seed": seed,
        "environment_hash_compatibility": {
            "mode": (
                "SOURCE_CHECKPOINT_HASH_EXPLICIT_OVERRIDE"
                if resume_environment_hash is not None else "RUNTIME_HASH"
            ),
            "source_checkpoint_environment_hash": resume_environment_hash,
            "runtime_environment_hash": (
                cell.get("fit_diagnostics", {}).get("runtime_environment_hash")
                if isinstance(cell.get("fit_diagnostics"), dict) else None
            ),
        },
        "wiring_only": not practical,
        "selected_k": None,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
