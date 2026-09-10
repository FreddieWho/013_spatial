#!/usr/bin/env python3
"""Held-out full-gene field export only, from a frozen Torch checkpoint.

Resumes the source checkpoint with zero fit steps and runs a single
full-gene held-out inference (no gene cross-fit splits).  Complements
scripts/r04_export_frozen_effects.py for cases where the split loop already
completed or timed out: same estimator configuration, same protocol values
read from the source cell (fail closed on missing keys), fresh replay
checkpoint directory, GT never read.

Acceptance follows D-103 minus the loss-equivalence item (a no-split
inference has no source endpoint to compare against; convergence platform,
zero fit updates and hash identity still apply).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from r04.diagnostics import platform_convergence_summary
from r04.field_exports import write_heldout_field_export
from r04.models import MNSFConfig, MNSFEstimator
from r04.runtime import atomic_json, resource_status
from scripts.r04_real_k_search import _grouped_folds, _load_training

# Frozen protocol value; inert under zero-step resume (asserted below).
# Kept as a literal so a silent default change elsewhere cannot leak in.
FROZEN_LEARNING_RATE = 0.01


def _read_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--cell", type=Path, required=True)
    parser.add_argument("--checkpoint-dir", type=Path, required=True)
    parser.add_argument("--training-manifest", type=Path, required=True)
    parser.add_argument("--gene-list", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True,
                        help="heldout_full.npz destination")
    parser.add_argument("--replay-checkpoint-dir", type=Path, required=True,
                        help="fresh directory for replay checkpoints (never the source)")
    parser.add_argument("--fold", type=int, required=True)
    parser.add_argument("--expect-input-hash", type=str, required=True)
    parser.add_argument("--expect-config-hash", type=str, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    cell = _read_json(root / args.cell)
    if cell.get("schema") != "r04.real_k_cell.v1":
        raise SystemExit("source cell schema mismatch")
    if cell.get("fold") != args.fold:
        raise SystemExit("cell fold differs from requested fold")
    params = cell.get("parameters")
    if not isinstance(params, dict):
        raise SystemExit("source cell parameters are missing")

    def need(container: dict, key: str):
        if key not in container or container[key] is None:
            raise SystemExit(f"source cell is missing required field: {key}")
        return container[key]

    k_model = cell.get("k_model")
    steps = need(params, "steps")
    inference_steps = need(params, "inference_steps")
    seed = need(params, "optimization_seed")
    fit_diag = cell.get("fit_diagnostics", {})
    if not isinstance(fit_diag, dict):
        raise SystemExit("source cell fit diagnostics are missing")
    source_env = fit_diag.get("environment_hash")
    if not source_env:
        raise SystemExit("source cell has no environment hash for explicit override")
    if cell.get("input_hash") != args.expect_input_hash:
        raise SystemExit("source cell input hash differs from expected")
    if cell.get("config_hash") != args.expect_config_hash:
        raise SystemExit("source cell config hash differs from expected")

    metadata = _read_json(root / args.checkpoint_dir / "checkpoint_metadata.json")
    if metadata.get("checkpoint_schema") != "r04.mnsf_checkpoint.v3_torch":
        raise SystemExit("source checkpoint is not a Torch v3 checkpoint")
    if not (root / args.checkpoint_dir / "checkpoint.pt").is_file():
        raise SystemExit("source checkpoint.pt is missing")
    replay_dir = root / args.replay_checkpoint_dir
    if replay_dir.exists():
        raise SystemExit(f"replay checkpoint dir already exists: {replay_dir}")

    sections, genes, manifest_hash = _load_training(
        root / args.training_manifest, root / args.gene_list)
    fold_ids = _grouped_folds(sections, folds=5, seed=20260807)
    training = [s for i, s in enumerate(sections) if int(fold_ids[i]) != args.fold]
    validation = [s for i, s in enumerate(sections) if int(fold_ids[i]) == args.fold]
    if not training or not validation:
        raise SystemExit("empty training or validation split")

    out_path = root / args.output
    if resource_status(out_path.parent) == "BLOCKED_STORAGE":
        raise SystemExit("BLOCKED_STORAGE")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    estimator = MNSFEstimator(MNSFConfig(
        factors=k_model,
        inducing_points=need(params, "inducing_points"),
        nonspatial_rank=1,
        lengthscale=need(params, "lengthscale"),
        learning_rate=FROZEN_LEARNING_RATE,
        learning_rate_final=None,
        learning_rate_decay_steps=None,
        steps=steps,
        posterior_draws=8,
        gene_batch_size=fit_diag.get("gene_batch_size"),
        diagnostic_interval=fit_diag.get("diagnostic_interval"),
        evaluation_interval=fit_diag.get("evaluation_interval"),
        evaluation_mc_draws=fit_diag.get("evaluation_mc_draws"),
        optimization_schedule=need(params, "optimization_schedule"),
        shared_steps=need(params, "shared_steps"),
        checkpoint_dir=str(replay_dir),
        resume_checkpoint_dir=str(root / args.checkpoint_dir),
        checkpoint_steps=steps,
        seed=seed,
        environment_hash=source_env,
        execution_mode="eager",
    )).fit(training)
    if estimator.diagnostics_.get("optimizer_steps_this_call") != 0:
        raise SystemExit("heldout-only replay performed fit updates")
    if estimator.fit_input_hash_ != args.expect_input_hash:
        raise SystemExit("replay input hash differs from expected source")
    if estimator.fit_config_hash_ != args.expect_config_hash:
        raise SystemExit("replay config hash differs from expected source")

    heldout_fields = estimator.infer(validation, steps=inference_steps, posterior_draws=8)
    platform = platform_convergence_summary(
        estimator.inference_diagnostics_.get("loss_trace", []),
        window=50, stable_windows=2)
    if platform.get("converged") is not True:
        raise SystemExit("heldout-only inference platform not converged")
    record = write_heldout_field_export(
        out_path,
        sections=validation,
        field_groups=heldout_fields,
        evaluation_loading=estimator.loading_,
        factor_amplitude=estimator.factor_amplitude_,
        frozen_gene_ids=genes,
        adaptation_gene_ids=(),
        evaluation_gene_ids=genes,
        adaptation_indices=np.asarray([], dtype=int),
        evaluation_indices=np.arange(len(genes), dtype=int),
        provenance={
            "k_model": k_model,
            "fold": args.fold,
            "gene_split": None,
            "optimization_seed": seed,
            "inference_steps": inference_steps,
            "fit_input_hash": estimator.fit_input_hash_,
            "fit_config_hash": estimator.fit_config_hash_,
            "fit_environment_hash": estimator.fit_environment_hash_,
            "environment_hash_override": True,
            "source_environment_hash": source_env,
            "inference_platform": platform,
            "representation_role": "heldout_outer_fold_readout_apply",
            "source_cell_path": str(args.cell),
            "source_cell_sha256": _sha256(root / args.cell),
            "source_checkpoint_dir": str(args.checkpoint_dir),
            "fit_updates": 0,
            "gpu_visible": False,
            "validation_gt_read": "NOT_READ",
            "internal_external_validation_gt": "SEALED_NOT_READ",
            "driver": "r04_export_heldout_only.py",
            "acceptance": "D-103 minus loss-equivalence (no-split inference has no source endpoint)",
        },
    )
    atomic_json(out_path.with_suffix(".record.json"), {
        "schema": "r04.heldout_only_export_record.v1",
        "status": "EXPORTED",
        "npz_sha256": record.get("npz_sha256"),
        "inference_platform": platform,
        "optimizer_steps_this_call": estimator.diagnostics_.get("optimizer_steps_this_call"),
    })
    print(f"HELDOUT_EXPORT_OK steps={inference_steps} converged={platform.get('converged')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
