#!/usr/bin/env python3
"""Run bounded synthetic calibration for continuous-field identifiability."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from r04.diagnostics import (
    canonical_permutation_cutoffs,
    canonical_subspace_correlations,
    mnsf_spatial_effect_matrix,
    platform_convergence_summary,
)
from r04.models import (
    MNSFConfig,
    MNSFEstimator,
    SignedResidualGPConfig,
    SignedResidualGPEstimator,
)
from r04.runtime import atomic_json
from r04.synthetic import make_identifiability_calibration_sections


def _field_matrix(estimator: object, factors: int) -> np.ndarray:
    groups = getattr(estimator, "fields_")
    return np.column_stack([
        np.concatenate([group[factor].field_mean for group in groups])
        for factor in range(factors)
    ])


def _fit_one(
    sections: list,
    *,
    model: str,
    factors: int,
    inducing_points: int,
    nonspatial_rank: int,
    lengthscale: float,
    learning_rate: float,
    steps: int,
    seed: int,
) -> object:
    if model == "mnsf":
        return MNSFEstimator(MNSFConfig(
            factors=factors,
            inducing_points=inducing_points,
            nonspatial_rank=nonspatial_rank,
            lengthscale=lengthscale,
            learning_rate=learning_rate,
            steps=steps,
            posterior_draws=8,
            seed=seed,
        )).fit(sections)
    if model == "signed":
        return SignedResidualGPEstimator(SignedResidualGPConfig(
            factors=factors,
            inducing_points=inducing_points,
            lengthscale=lengthscale,
            learning_rate=learning_rate,
            steps=steps,
            posterior_draws=8,
            seed=seed,
        )).fit(sections)
    raise ValueError(f"unsupported model: {model}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--model", choices=("mnsf", "signed"), default="mnsf")
    parser.add_argument("--domain", choices=("crescent", "branch", "disconnected"), default="disconnected")
    parser.add_argument("--mode", choices=("no_field", "one_disconnected", "two_independent", "two_collinear"), default="two_independent")
    parser.add_argument("--generator", choices=("mnsf_correct", "signed_correct", "misspecified"), default="mnsf_correct")
    parser.add_argument("--factors", type=int, default=2)
    parser.add_argument("--inducing-points", type=int, default=16)
    parser.add_argument("--nonspatial-rank", type=int, default=0)
    parser.add_argument("--lengthscale", type=float, default=3.0)
    parser.add_argument("--learning-rate", type=float, default=0.05)
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20260807)
    parser.add_argument("--null-draws", type=int, default=200)
    parser.add_argument("--platform-window", type=int, default=50)
    args = parser.parse_args()
    if (
        args.factors < 1
        or args.inducing_points < 2
        or args.nonspatial_rank < 0
        or args.lengthscale <= 0
        or args.learning_rate <= 0
        or args.steps < 1
        or args.null_draws < 1
    ):
        raise SystemExit("invalid calibration configuration")

    sections, truth, metadata = make_identifiability_calibration_sections(
        mode=args.mode,
        domain=args.domain,
        generator=args.generator,
        sections=3,
        spots_per_section=60,
        genes=24,
        seed=args.seed,
    )
    estimator = _fit_one(
        sections,
        model=args.model,
        factors=args.factors,
        inducing_points=args.inducing_points,
        nonspatial_rank=args.nonspatial_rank,
        lengthscale=args.lengthscale,
        learning_rate=args.learning_rate,
        steps=args.steps,
        seed=args.seed,
    )
    estimate = _field_matrix(estimator, args.factors)
    groups = np.concatenate([
        np.repeat(section.section_id, len(section.barcode)) for section in sections
    ])
    diagnostics = dict(getattr(estimator, "diagnostics_", {}))
    field_correlations = (
        canonical_subspace_correlations(truth, estimate)
        if truth.shape[1] else np.empty(0)
    )
    field_cutoffs = (
        canonical_permutation_cutoffs(
            truth,
            estimate,
            groups=groups,
            seed=args.seed + 400,
            draws=args.null_draws,
        )
        if truth.shape[1] else np.empty(0)
    )

    effect_correlations = np.empty(0)
    effect_cutoffs = np.empty(0)
    if args.model == "mnsf" and args.generator == "mnsf_correct" and truth.shape[1]:
        true_effect = mnsf_spatial_effect_matrix(
            truth,
            np.asarray(metadata["true_loading"], dtype=float),
            np.asarray(metadata["true_amplitude"], dtype=float),
        )
        estimated_effect = mnsf_spatial_effect_matrix(
            estimate,
            np.asarray(estimator.loading_, dtype=float),
            np.asarray(estimator.factor_amplitude_, dtype=float),
        )
        effect_correlations = canonical_subspace_correlations(true_effect, estimated_effect)
        effect_cutoffs = canonical_permutation_cutoffs(
            true_effect,
            estimated_effect,
            groups=groups,
            seed=args.seed + 800,
            draws=args.null_draws,
        )

    main_correlations = effect_correlations if len(effect_correlations) else field_correlations
    main_cutoffs = effect_cutoffs if len(effect_cutoffs) else field_cutoffs
    platform = platform_convergence_summary(
        diagnostics.get("loss_trace", []),
        window=args.platform_window,
        stable_windows=2,
    )
    result = {
        "schema": "r04.identifiability_calibration.v2",
        "status": "CALIBRATION_COMPLETE_NOT_VALIDATED",
        "model": args.model,
        "steps": args.steps,
        "seed": args.seed,
        "n_factors_model": args.factors,
        "truth": metadata,
        "config": {
            "inducing_points": args.inducing_points,
            "nonspatial_rank": args.nonspatial_rank,
            "lengthscale": args.lengthscale,
            "learning_rate": args.learning_rate,
            "null_draws": args.null_draws,
            "platform_window": args.platform_window,
        },
        "subspace": {
            "target": "mnsf_spot_by_gene_spatial_rate" if len(effect_correlations) else "continuous_field",
            "canonical_correlations": main_correlations.tolist(),
            "direction_matched_null_95": main_cutoffs.tolist(),
            "all_directions_above_null": bool(
                len(main_correlations) and np.all(main_correlations > main_cutoffs)
            ),
            "field_canonical_correlations": field_correlations.tolist(),
            "field_direction_matched_null_95": field_cutoffs.tolist(),
        },
        "platform": platform,
        "model_diagnostics": diagnostics,
        "formal_restart_authorized": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(args.output, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
