#!/usr/bin/env python3
"""Fit-only late-learning-rate continuation from an explicit R-04 checkpoint.

This is a wiring-only optimizer diagnostic.  It does not perform inference or
scoring, and requires an explicit source checkpoint so that Adam state, input
hashes, and the step sequence are preserved across the branch.  The historical
default is K=0; ``--factors`` also permits a matched continuation for K>0.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from r04.diagnostics import (
    deterministic_objective_platform_summary,
    platform_convergence_summary,
)
from r04.models import MNSFConfig, MNSFEstimator
from r04.runtime import atomic_json, resource_status
from scripts.r04_real_k_search import (
    _grouped_folds,
    _load_training,
    _optimization_seed,
    _parse_gene_batch_size,
)


def _diagnostic_payload(
    *,
    factors: int,
    fold: int,
    steps: int,
    expected_start_step: int,
    learning_rate: float,
    learning_rate_final: float | None,
    learning_rate_decay_steps: int | None,
    gene_batch_size: int | None,
    evaluation_interval: int,
    evaluation_mc_draws: int,
    fit_converged: bool,
    source_checkpoint_dir: str,
    restart_index: int,
    optimization_seed: int,
    fit_objective_platform: dict[str, object] | None = None,
    legacy_minibatch_fit_gate: dict[str, object] | None = None,
) -> dict[str, object]:
    """Build the top-level wiring-only summary without losing the model K."""
    return {
        "schema": "r04.late_lr_diagnostic.v2",
        "status": "WIRING_ONLY_NOT_SCIENTIFIC",
        "k_model": factors,
        "fold": fold,
        "steps": steps,
        "expected_start_step": expected_start_step,
        "learning_rate": learning_rate,
        "learning_rate_final": learning_rate_final,
        "learning_rate_decay_steps": learning_rate_decay_steps,
        "gene_batch_size": gene_batch_size,
        "evaluation_interval": evaluation_interval,
        "evaluation_mc_draws": evaluation_mc_draws,
        "fit_converged": fit_converged,
        "fit_objective_platform": fit_objective_platform or {"converged": fit_converged},
        "legacy_minibatch_fit_gate": legacy_minibatch_fit_gate or {},
        "cell_count": 1,
        "selected_k": None,
        "wiring_only": True,
        "optimizer_state_inherited": True,
        "source_checkpoint_dir": source_checkpoint_dir,
        "restart_index": restart_index,
        "optimization_seed": optimization_seed,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument("--gene-list", type=Path, required=True)
    parser.add_argument("--source-checkpoint-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fold", type=int, default=0)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--group-seed", type=int, default=20260807)
    parser.add_argument("--restart-index", type=int, default=0)
    parser.add_argument("--factors", type=int, default=0)
    parser.add_argument("--expected-start-step", type=int, default=4800)
    parser.add_argument("--steps", type=int, default=7200)
    parser.add_argument("--checkpoint-steps", type=int, default=600)
    parser.add_argument("--inducing-points", type=int, default=16)
    parser.add_argument("--lengthscale", type=float, default=3.0)
    parser.add_argument("--learning-rate", type=float, default=0.01)
    parser.add_argument("--learning-rate-final", type=float, default=None)
    parser.add_argument("--learning-rate-decay-steps", type=int, default=None)
    parser.add_argument("--gene-batch-size", default="512")
    parser.add_argument("--diagnostic-interval", type=int, default=100)
    parser.add_argument("--evaluation-interval", type=int, default=0)
    parser.add_argument("--evaluation-mc-draws", type=int, default=4)
    args = parser.parse_args()

    if (
        args.folds < 2
        or args.factors < 0
        or args.restart_index < 0
        or args.expected_start_step < 0
        or args.steps <= args.expected_start_step
    ):
        raise SystemExit("invalid folds or continuation step range")
    if args.learning_rate <= 0:
        raise SystemExit("learning-rate must be positive")
    if args.learning_rate_final is not None and args.learning_rate_final <= 0:
        raise SystemExit("learning-rate-final must be positive")
    if args.learning_rate_final is not None and (
        args.learning_rate_decay_steps is None or args.learning_rate_decay_steps < 1
    ):
        raise SystemExit("learning-rate-decay-steps is required with learning-rate-final")
    if args.learning_rate_final is None and args.learning_rate_decay_steps is not None:
        raise SystemExit("learning-rate-final is required with learning-rate-decay-steps")
    if args.evaluation_interval < 0 or args.evaluation_mc_draws < 1:
        raise SystemExit("evaluation interval must be non-negative and draws must be positive")
    try:
        gene_batch_size = _parse_gene_batch_size(args.gene_batch_size)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    args.output_dir.mkdir(parents=True, exist_ok=True)
    if resource_status(args.output_dir) == "BLOCKED_STORAGE":
        atomic_json(args.output_dir / "r04_late_lr_diagnostic.json", {"status": "BLOCKED_STORAGE"})
        return 2
    if not args.source_checkpoint_dir.exists():
        raise SystemExit(f"source checkpoint directory does not exist: {args.source_checkpoint_dir}")

    atomic_json(args.output_dir / "r04_late_lr_diagnostic.json", {
        "schema": "r04.late_lr_diagnostic.v1",
        "status": "RUNNING",
        "k_model": args.factors,
        "fold": args.fold,
        "steps": args.steps,
        "expected_start_step": args.expected_start_step,
        "learning_rate": args.learning_rate,
        "learning_rate_final": args.learning_rate_final,
        "learning_rate_decay_steps": args.learning_rate_decay_steps,
        "gene_batch_size": gene_batch_size,
        "evaluation_interval": args.evaluation_interval,
        "evaluation_mc_draws": args.evaluation_mc_draws,
        "source_checkpoint_dir": str(args.source_checkpoint_dir),
        "restart_index": args.restart_index,
        "wiring_only": True,
    })

    sections, genes, manifest_hash = _load_training(args.manifest_json, args.gene_list)
    fold_ids = _grouped_folds(sections, folds=args.folds, seed=args.group_seed)
    available_folds = sorted(set(int(value) for value in fold_ids))
    if args.fold not in available_folds:
        raise SystemExit(f"fold {args.fold} is not present: {available_folds}")
    training = [
        section for index, section in enumerate(sections) if int(fold_ids[index]) != args.fold
    ]
    seed = _optimization_seed(
        args.group_seed, args.factors, args.fold, args.restart_index
    )
    checkpoint_dir = args.output_dir / "checkpoints" / f"k{args.factors}" / f"fold{args.fold}"
    estimator = MNSFEstimator(MNSFConfig(
        factors=args.factors,
        inducing_points=args.inducing_points,
        nonspatial_rank=1,
        lengthscale=args.lengthscale,
        learning_rate=args.learning_rate,
        learning_rate_final=args.learning_rate_final,
        learning_rate_decay_steps=args.learning_rate_decay_steps,
        steps=args.steps,
        posterior_draws=8,
        gene_batch_size=gene_batch_size,
        diagnostic_interval=args.diagnostic_interval,
        evaluation_interval=args.evaluation_interval,
        evaluation_mc_draws=args.evaluation_mc_draws,
        checkpoint_dir=str(checkpoint_dir),
        resume_checkpoint_dir=str(args.source_checkpoint_dir),
        checkpoint_steps=args.checkpoint_steps,
        seed=seed,
    )).fit(training)
    diagnostics = estimator.diagnostics_
    if diagnostics.get("start_step") != args.expected_start_step:
        raise RuntimeError(
            f"unexpected continuation start step: {diagnostics.get('start_step')}"
        )
    fit_platform = platform_convergence_summary(
        diagnostics.get("loss_trace", []), window=50, stable_windows=2
    )
    evaluation_trace = diagnostics.get("evaluation_trace", [])
    fit_objective_platform = deterministic_objective_platform_summary(
        evaluation_trace if isinstance(evaluation_trace, list) else []
    )
    cell = {
        "schema": "r04.late_lr_cell.v2",
        "status": "FIT_ONLY_DIAGNOSTIC",
        "continuation_panel": True,
        "k_model": args.factors,
        "fold": args.fold,
        "n_training_sections": len(training),
        "training_patients": sorted({str(section.patient_id) for section in training}),
        "fit_platform": fit_platform,
        "fit_objective_platform": fit_objective_platform,
        "legacy_minibatch_fit_gate": fit_platform,
        "fit_converged": bool(fit_objective_platform["converged"]),
        "fit_diagnostics": diagnostics,
        "inference_platform": {"status": "NOT_RUN", "converged": False},
        "inference_diagnostics": {"status": "NOT_RUN"},
        "input_hash": estimator.fit_input_hash_,
        "config_hash": estimator.fit_config_hash_,
        "checkpoint_dir": str(checkpoint_dir),
        "source_checkpoint_dir": str(args.source_checkpoint_dir),
        "optimizer_state_inherited": True,
        "restart_index": args.restart_index,
        "optimization_seed": seed,
        "manifest_hash": manifest_hash,
    }
    atomic_json(args.output_dir / "cell.json", cell)
    atomic_json(args.output_dir / "r04_late_lr_diagnostic.json", _diagnostic_payload(
        factors=args.factors,
        fold=args.fold,
        steps=args.steps,
        expected_start_step=args.expected_start_step,
        learning_rate=args.learning_rate,
        learning_rate_final=args.learning_rate_final,
        learning_rate_decay_steps=args.learning_rate_decay_steps,
        gene_batch_size=gene_batch_size,
        evaluation_interval=args.evaluation_interval,
        evaluation_mc_draws=args.evaluation_mc_draws,
        fit_converged=bool(fit_objective_platform["converged"]),
        source_checkpoint_dir=str(args.source_checkpoint_dir),
        restart_index=args.restart_index,
        optimization_seed=seed,
        fit_objective_platform=fit_objective_platform,
        legacy_minibatch_fit_gate=fit_platform,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
