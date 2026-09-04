#!/usr/bin/env python3
"""Torch-era checkpoint-only held-out inference/scoring entry for the R-04 representative chain.

Runs zero-fit-update scoring from an audited frozen Torch v3 checkpoint:

1. verifies the endpoint through a ``r04.checkpoint_audit.v1`` artifact produced by
   ``scripts/r04_checkpoint_audit.py`` (objective-input hash must be ``VERIFIED``
   at the requested frozen step);
2. resumes the checkpoint with zero remaining fit steps and executes the fixed
   2400-step, two-gene-split held-out inference/scoring protocol via
   ``scripts.r04_real_k_search._score_fold``;
3. records environment-hash compatibility (explicit source-hash override is only
   for a validated remote runtime), the endpoint audit sha256, and
   ``selected_k=null``.

This entry exists because the strict continuation-panel gate in
``scripts/r04_resume_infer.py`` is bound to the TensorFlow-era panel audit
(TF checkpoint bundle fingerprints and the restart 2-4 panel cells) and cannot
audit a Torch representative checkpoint.  It changes no model math, protocol
parameter, fold, gene split, or threshold, and its output is a representative
diagnostic that never authorizes K selection.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from r04.io_contract import sections_content_hash, sections_objective_input_hash
from r04.runtime import atomic_json, resource_status
from scripts.r04_real_k_search import (
    _grouped_folds,
    _load_training,
    _optimization_seed,
    _parse_gene_batch_size,
    _score_fold,
)

EXPECTED_INFERENCE_PROTOCOL = {
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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument("--gene-list", type=Path, required=True)
    parser.add_argument("--source-checkpoint-dir", type=Path, required=True)
    parser.add_argument("--endpoint-audit-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--factors", type=int, required=True)
    parser.add_argument("--fold", type=int, default=0)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--group-seed", type=int, default=20260807)
    parser.add_argument("--restart-index", type=int, default=0)
    parser.add_argument("--gene-folds", type=int, default=2)
    parser.add_argument("--split-seed", type=int, default=20260817)
    parser.add_argument("--fit-steps", type=int, default=8400)
    parser.add_argument("--expected-start-step", type=int, default=8400)
    parser.add_argument("--inference-steps", type=int, default=2400)
    parser.add_argument("--checkpoint-steps", type=int, default=600)
    parser.add_argument("--inducing-points", type=int, default=16)
    parser.add_argument("--lengthscale", type=float, default=3.0)
    parser.add_argument("--learning-rate", type=float, default=0.01)
    parser.add_argument("--gene-batch-size", default="512")
    parser.add_argument("--diagnostic-interval", type=int, default=100)
    parser.add_argument(
        "--allow-source-environment-hash",
        action="store_true",
        help="explicitly reuse the source checkpoint environment hash for a validated remote runtime",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="run all audit/hash/protocol checks without the expensive inference",
    )
    args = parser.parse_args()

    if args.factors < 0 or args.restart_index < 0 or args.folds < 2 or args.gene_folds < 2:
        raise SystemExit("factors must be non-negative; folds and gene-folds must be at least 2")
    if args.fit_steps != args.expected_start_step:
        raise SystemExit("checkpoint-only entry requires --fit-steps == --expected-start-step")
    if args.inference_steps < 1:
        raise SystemExit("inference-steps must be positive")
    protocol = {
        "folds": args.folds,
        "group_seed": args.group_seed,
        "gene_folds": args.gene_folds,
        "split_seed": args.split_seed,
        "inference_steps": args.inference_steps,
        "inducing_points": args.inducing_points,
        "lengthscale": args.lengthscale,
        "learning_rate": args.learning_rate,
    }
    if protocol != EXPECTED_INFERENCE_PROTOCOL:
        raise SystemExit(
            f"inference protocol drift is not allowed for this entry: {protocol} != {EXPECTED_INFERENCE_PROTOCOL}"
        )
    if not args.source_checkpoint_dir.exists():
        raise SystemExit(f"source checkpoint directory does not exist: {args.source_checkpoint_dir}")
    if not args.endpoint_audit_json.exists():
        raise SystemExit(f"endpoint audit does not exist: {args.endpoint_audit_json}")
    try:
        gene_batch_size = _parse_gene_batch_size(args.gene_batch_size)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    # Endpoint audit: Torch-native checkpoint audit must verify the frozen step.
    audit = _read_json(args.endpoint_audit_json)
    if audit.get("schema") != "r04.checkpoint_audit.v1":
        raise SystemExit("endpoint audit schema must be r04.checkpoint_audit.v1")
    audit_entries = audit.get("checkpoints")
    if not isinstance(audit_entries, list) or not audit_entries:
        raise SystemExit("endpoint audit has no checkpoint entries")
    endpoint_matches = [
        item for item in audit_entries
        if isinstance(item, dict) and item.get("global_step") == args.fit_steps
    ]
    if len(endpoint_matches) != 1:
        raise SystemExit("endpoint audit does not contain exactly one frozen-step entry")
    endpoint = endpoint_matches[0]
    if endpoint.get("objective_input_hash_status") != "VERIFIED":
        raise SystemExit("endpoint audit objective input hash is not VERIFIED")
    audit_protocol = {
        key: audit.get(key)
        for key in ("factors", "fold", "folds", "group_seed", "inducing_points", "lengthscale")
    }
    expected_audit_protocol = {
        "factors": args.factors,
        "fold": args.fold,
        "folds": args.folds,
        "group_seed": args.group_seed,
        "inducing_points": args.inducing_points,
        "lengthscale": args.lengthscale,
    }
    if audit_protocol != expected_audit_protocol:
        raise SystemExit("endpoint audit protocol does not match this run")
    audit_sha256 = _sha256(args.endpoint_audit_json)

    # Source checkpoint metadata checks (fail closed).
    metadata = _read_json(args.source_checkpoint_dir / "checkpoint_metadata.json")
    if (
        metadata.get("checkpoint_schema") != "r04.mnsf_checkpoint.v3_torch"
        or metadata.get("backend") != "torch"
    ):
        raise SystemExit("source checkpoint schema/backend mismatch")
    if int(metadata.get("factors", -1)) != args.factors:
        raise SystemExit("source checkpoint factor count mismatch")
    source_environment_hash = metadata.get("environment_hash")
    if not isinstance(source_environment_hash, str) or not source_environment_hash:
        raise SystemExit("source checkpoint metadata has no environment hash")
    resume_environment_hash = (
        source_environment_hash if args.allow_source_environment_hash else None
    )

    sections, genes, manifest_hash = _load_training(args.manifest_json, args.gene_list)
    if audit.get("manifest_hash") != manifest_hash:
        raise SystemExit("endpoint audit manifest hash does not match the supplied manifest")
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

    # Replicate the fail-closed hash checks the resume path performs in fit().
    input_hash = sections_content_hash(training)
    objective_input_hash = sections_objective_input_hash(training)
    if metadata.get("input_hash") != input_hash:
        raise SystemExit("source checkpoint input hash does not match the training sections")
    if metadata.get("objective_input_hash") != objective_input_hash:
        raise SystemExit("source checkpoint objective input hash does not match the training sections")
    for entry in audit_entries:
        if isinstance(entry, dict) and entry.get("objective_input_hash") not in (None, objective_input_hash):
            raise SystemExit("endpoint audit objective input hash mismatch")

    seed = _optimization_seed(args.group_seed, args.factors, args.fold, args.restart_index)
    if args.validate_only:
        atomic_json(Path("/tmp/r04_torch_frozen_infer_validate.json"), {
            "schema": "r04.torch_frozen_infer.v1",
            "status": "VALIDATE_ONLY_OK",
            "k_model": args.factors,
            "fold": args.fold,
            "fit_steps": args.fit_steps,
            "input_hash": input_hash,
            "objective_input_hash": objective_input_hash,
            "manifest_hash": manifest_hash,
            "endpoint_audit_sha256": audit_sha256,
            "environment_hash_compatibility": {
                "mode": (
                    "SOURCE_CHECKPOINT_HASH_EXPLICIT_OVERRIDE"
                    if resume_environment_hash is not None else "RUNTIME_HASH"
                ),
                "source_checkpoint_environment_hash": source_environment_hash,
            },
            "selected_k": None,
        })
        return 0

    args.output_dir.mkdir(parents=True, exist_ok=True)
    if resource_status(args.output_dir) == "BLOCKED_STORAGE":
        atomic_json(args.output_dir / "r04_torch_frozen_infer.json", {"status": "BLOCKED_STORAGE"})
        return 2

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
    if diagnostics.get("optimizer_steps_this_call") != 0:
        raise RuntimeError("checkpoint-only entry unexpectedly performed fit updates")

    cell["status"] = "FIT_AND_SCORED_TORCH_REPRESENTATIVE_DIAGNOSTIC"
    cell["selected_k"] = None
    cell["wiring_only"] = True
    cell["restart_index"] = args.restart_index
    cell["optimization_seed"] = seed
    cell["manifest_hash"] = manifest_hash
    cell["source_checkpoint_dir"] = str(args.source_checkpoint_dir)
    cell["endpoint_audit_json"] = str(args.endpoint_audit_json)
    cell["endpoint_audit_sha256"] = audit_sha256
    cell["endpoint_dense_objectives"] = [
        {
            "global_step": item.get("global_step"),
            "objective_mean": item.get("objective_mean"),
            "objective_sd": item.get("objective_sd"),
            "objective_input_hash_status": item.get("objective_input_hash_status"),
        }
        for item in audit_entries
        if isinstance(item, dict)
    ]
    cell["environment_hash_compatibility"] = {
        "mode": (
            "SOURCE_CHECKPOINT_HASH_EXPLICIT_OVERRIDE"
            if resume_environment_hash is not None else "RUNTIME_HASH"
        ),
        "source_checkpoint_environment_hash": source_environment_hash,
        "runtime_environment_hash": diagnostics.get("runtime_environment_hash"),
    }
    atomic_json(args.output_dir / "cell.json", cell)
    atomic_json(args.output_dir / "r04_torch_frozen_infer.json", {
        "schema": "r04.torch_frozen_infer.v1",
        "status": "TORCH_REPRESENTATIVE_DIAGNOSTIC_NOT_K_SELECTION",
        "k_model": args.factors,
        "fold": args.fold,
        "fit_steps": args.fit_steps,
        "inference_steps": args.inference_steps,
        "restart_index": args.restart_index,
        "optimization_seed": seed,
        "fit_platform": cell.get("fit_platform"),
        "inference_platform": cell.get("inference_platform"),
        "source_checkpoint_dir": str(args.source_checkpoint_dir),
        "endpoint_audit_json": str(args.endpoint_audit_json),
        "endpoint_audit_sha256": audit_sha256,
        "environment_hash_compatibility": cell["environment_hash_compatibility"],
        "wiring_only": True,
        "selected_k": None,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
