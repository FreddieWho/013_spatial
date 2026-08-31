#!/usr/bin/env python3
"""Read-only deterministic objective audit for existing R-04 checkpoints."""

from __future__ import annotations

import argparse
from pathlib import Path

from r04.checkpoint_audit import evaluate_checkpoint
from r04.runtime import atomic_json
from scripts.r04_real_k_search import _grouped_folds, _load_training


def _audit_payload(
    *,
    manifest_hash: str,
    gene_count: int,
    factors: int,
    fold: int,
    folds: int,
    group_seed: int,
    inducing_points: int,
    lengthscale: float,
    ridge: float,
    mc_draws: int,
    checkpoints: list[dict[str, object]],
) -> dict[str, object]:
    return {
        "schema": "r04.checkpoint_audit.v1",
        "status": "DIAGNOSTIC_ONLY",
        "manifest_hash": manifest_hash,
        "gene_count": gene_count,
        "factors": factors,
        "fold": fold,
        "folds": folds,
        "group_seed": group_seed,
        "inducing_points": inducing_points,
        "lengthscale": lengthscale,
        "ridge": ridge,
        "mc_draws": mc_draws,
        "checkpoints": checkpoints,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument("--gene-list", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, action="append", required=True)
    parser.add_argument("--factors", type=int, required=True)
    parser.add_argument("--fold", type=int, default=0)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--group-seed", type=int, default=20260807)
    parser.add_argument("--inducing-points", type=int, default=16)
    parser.add_argument("--lengthscale", type=float, default=3.0)
    parser.add_argument("--ridge", type=float, default=1e-5)
    parser.add_argument("--mc-draws", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260820)
    args = parser.parse_args()
    sections, genes, manifest_hash = _load_training(args.manifest_json, args.gene_list)
    fold_ids = _grouped_folds(sections, folds=args.folds, seed=args.group_seed)
    if args.fold not in set(int(value) for value in fold_ids):
        raise SystemExit(f"fold not present: {args.fold}")
    training = [
        section for index, section in enumerate(sections)
        if int(fold_ids[index]) != args.fold
    ]
    checkpoints = []
    for checkpoint in args.checkpoint:
        if not checkpoint.exists() and not Path(f"{checkpoint}.index").exists():
            raise SystemExit(f"checkpoint does not exist: {checkpoint}")
        checkpoints.append(evaluate_checkpoint(
            training,
            checkpoint,
            factors=args.factors,
            inducing_points=args.inducing_points,
            lengthscale=args.lengthscale,
            ridge=args.ridge,
            mc_draws=args.mc_draws,
            seed=args.seed,
        ))
    atomic_json(args.output_json, _audit_payload(
        manifest_hash=manifest_hash,
        gene_count=len(genes),
        factors=args.factors,
        fold=args.fold,
        folds=args.folds,
        group_seed=args.group_seed,
        inducing_points=args.inducing_points,
        lengthscale=args.lengthscale,
        ridge=args.ridge,
        mc_draws=args.mc_draws,
        checkpoints=checkpoints,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
