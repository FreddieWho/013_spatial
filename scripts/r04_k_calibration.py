#!/usr/bin/env python3
"""Leakage-free held-out K calibration on continuous-field controls."""

from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln

from r04.diagnostics import (
    canonical_permutation_cutoffs,
    canonical_subspace_correlations,
    mnsf_spatial_effect_matrix,
    platform_convergence_summary,
)
from r04.metrics import negative_binomial_log_likelihood
from r04.models import MNSFConfig, MNSFEstimator
from r04.models.mnsf import predict_log_mean_from_fields
from r04.io_contract import sections_content_hash
from r04.runtime import atomic_json
from r04.synthetic import make_identifiability_calibration_sections


def gene_crossfit_splits(
    genes: int,
    *,
    folds: int,
    seed: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Split genes so evaluation counts are never used to adapt held-out fields."""
    if genes < 2 or folds < 2 or folds > genes:
        raise ValueError("gene cross-fit needs 2 <= folds <= genes")
    permutation = np.random.default_rng(seed).permutation(genes)
    evaluation_folds = [np.sort(values) for values in np.array_split(permutation, folds)]
    universe = np.arange(genes, dtype=int)
    return [
        (np.setdiff1d(universe, evaluation, assume_unique=True), evaluation)
        for evaluation in evaluation_folds
    ]


def required_recovered_directions(k_model: int, k_eff_true: int) -> int:
    if k_model < 1 or k_eff_true < 0:
        raise ValueError("K dimensions must be non-negative with positive K_model")
    return min(k_model, k_eff_true)


def _subset_genes(section: object, indices: np.ndarray) -> object:
    counts = section.counts[:, indices]
    return replace(
        section,
        counts=counts,
        gene_id=tuple(section.gene_id[index] for index in indices),
    )


def _fit_intercept_nb(sections: list) -> tuple[np.ndarray, np.ndarray]:
    counts = np.concatenate([
        section.counts.toarray() if hasattr(section.counts, "toarray") else np.asarray(section.counts)
        for section in sections
    ], axis=0).astype(float)
    library = np.concatenate([
        np.asarray(section.library_size, dtype=float)
        if section.library_size is not None
        else np.asarray(section.counts.sum(axis=1)).ravel().astype(float)
        for section in sections
    ])
    rate = np.maximum(counts.sum(axis=0) / np.maximum(library.sum(), 1e-12), 1e-8)
    mean = np.maximum(library[:, None] * rate[None, :], 1e-8)
    dispersion = np.empty(counts.shape[1], dtype=float)
    for gene in range(counts.shape[1]):
        y = counts[:, gene]
        mu = mean[:, gene]

        def objective(log_theta: float) -> float:
            theta = float(np.exp(log_theta))
            denominator = theta + mu
            value = (
                gammaln(y + theta)
                - gammaln(theta)
                - gammaln(y + 1.0)
                + theta * (np.log(theta) - np.log(denominator))
                + y * (np.log(mu) - np.log(denominator))
            )
            return -float(np.sum(value))

        fitted = minimize_scalar(objective, bounds=(-5.0, 14.0), method="bounded")
        dispersion[gene] = float(np.exp(fitted.x))
    return rate, dispersion


def _score_fold(
    sections: list,
    truth_by_section: list[np.ndarray],
    truth_metadata: dict[str, object],
    *,
    holdout: int,
    factors: int,
    steps: int,
    data_seed: int,
    optimization_seed: int,
    split_seed: int,
    gene_folds: int,
    null_draws: int,
    inference_steps: int,
) -> dict[str, object]:
    training = [section for index, section in enumerate(sections) if index != holdout]
    validation = sections[holdout]
    baseline_rate, baseline_dispersion = _fit_intercept_nb(training)
    estimator = MNSFEstimator(MNSFConfig(
        factors=factors,
        inducing_points=16,
        nonspatial_rank=0,
        lengthscale=3.0,
        learning_rate=0.05,
        steps=steps,
        posterior_draws=4,
        seed=optimization_seed,
    )).fit(training)
    fit_platform = platform_convergence_summary(
        estimator.diagnostics_.get("loss_trace", []), window=50, stable_windows=2
    )
    scores: list[float] = []
    baseline_scores: list[float] = []
    recovery: list[dict[str, object]] = []
    splits = gene_crossfit_splits(
        len(validation.gene_id), folds=gene_folds, seed=split_seed
    )
    true_loading = np.asarray(truth_metadata["true_loading"], dtype=float)
    true_amplitude = np.asarray(truth_metadata["true_amplitude"], dtype=float)
    for split_index, (adaptation, evaluation) in enumerate(splits):
        projected = estimator.subset_to_genes(
            tuple(validation.gene_id[index] for index in adaptation)
        )
        inferred = projected.infer(
            [_subset_genes(validation, adaptation)],
            steps=inference_steps,
            posterior_draws=4,
        )
        log_mu = predict_log_mean_from_fields(
            [validation],
            inferred,
            estimator.loading_,
            estimator.factor_amplitude_,
            estimator.gene_baseline_,
        )
        counts = (
            validation.counts.toarray()
            if hasattr(validation.counts, "toarray")
            else np.asarray(validation.counts)
        )
        scores.append(negative_binomial_log_likelihood(
            counts[:, evaluation], log_mu[:, evaluation], estimator.dispersion_[evaluation]
        ))
        library = (
            np.asarray(validation.library_size, dtype=float)
            if validation.library_size is not None else counts.sum(axis=1).astype(float)
        )
        baseline_log_mu = (
            np.log(np.maximum(library, 1e-8))[:, None]
            + np.log(baseline_rate[evaluation])[None, :]
        )
        baseline_scores.append(negative_binomial_log_likelihood(
            counts[:, evaluation], baseline_log_mu, baseline_dispersion[evaluation]
        ))
        estimated_fields = np.column_stack([field.field_mean for field in inferred[0]])
        estimated_effect = mnsf_spatial_effect_matrix(
            estimated_fields,
            estimator.loading_[evaluation],
            estimator.factor_amplitude_,
        )
        true_effect = mnsf_spatial_effect_matrix(
            truth_by_section[holdout],
            true_loading[evaluation],
            true_amplitude,
        )
        correlations = canonical_subspace_correlations(true_effect, estimated_effect)
        cutoffs = canonical_permutation_cutoffs(
            true_effect,
            estimated_effect,
            seed=split_seed + split_index * 101 + 7000,
            draws=null_draws,
        )
        infer_trace = inferred[0][0].diagnostics.get("loss_trace", [])
        recovery.append({
            "adaptation_genes": adaptation.tolist(),
            "evaluation_genes": evaluation.tolist(),
            "canonical_correlations": correlations.tolist(),
            "direction_matched_null_95": cutoffs.tolist(),
            "inference_platform": platform_convergence_summary(
                infer_trace, window=50, stable_windows=2
            ),
        })
    required = required_recovered_directions(
        factors, int(truth_metadata["k_eff_true"])
    )
    recovery_pass = all(
        len(item["canonical_correlations"]) >= required
        and all(
            item["canonical_correlations"][index]
            > item["direction_matched_null_95"][index]
            for index in range(required)
        )
        for item in recovery
    )
    return {
        "holdout": holdout,
        "data_seed": data_seed,
        "optimization_seed": optimization_seed,
        "split_seed": split_seed,
        "k_model": factors,
        "scores": scores,
        "baseline_scores": baseline_scores,
        "fit_platform": fit_platform,
        "recovery": recovery,
        "recovery_pass": recovery_pass,
        "required_recovered_directions": required,
    }


def _parse_ints(value: str) -> list[int]:
    result = sorted(set(int(item) for item in value.split(",") if item.strip()))
    if not result or any(item < 1 for item in result):
        raise argparse.ArgumentTypeError("expected comma-separated positive integers")
    return result


def _parse_nonnegative_ints(value: str) -> list[int]:
    result = sorted(set(int(item) for item in value.split(",") if item.strip()))
    if not result or any(item < 0 for item in result):
        raise argparse.ArgumentTypeError("expected comma-separated non-negative integers")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--seeds", type=_parse_ints, default=[20260813])
    parser.add_argument("--k-values", type=_parse_ints, default=[1, 2, 3, 4])
    parser.add_argument("--section-folds", type=int, default=2)
    parser.add_argument("--holdouts", type=_parse_nonnegative_ints)
    parser.add_argument("--gene-folds", type=int, default=2)
    parser.add_argument("--null-draws", type=int, default=50)
    parser.add_argument("--inference-steps", type=int, default=250)
    parser.add_argument("--mode", choices=("no_field", "one_disconnected", "two_independent", "two_collinear"), default="two_independent")
    args = parser.parse_args()
    if (
        args.steps < 150
        or args.inference_steps < 150
        or (args.holdouts is None and args.section_folds < 2)
        or args.gene_folds < 2
    ):
        raise SystemExit("K calibration needs at least 150 steps and two folds")
    if len(args.seeds) != 1:
        raise SystemExit(
            "each K calibration artifact must use exactly one data seed; "
            "write separate shards before aggregation"
        )
    generator_seed = args.seeds[0]
    sections, truth, metadata = make_identifiability_calibration_sections(
        mode=args.mode,
        domain="disconnected",
        generator="mnsf_correct",
        sections=3,
        spots_per_section=50,
        genes=24,
        seed=generator_seed,
    )
    truth_by_section: list[np.ndarray] = []
    offset = 0
    for section in sections:
        truth_by_section.append(truth[offset:offset + len(section.barcode)])
        offset += len(section.barcode)
    holdouts = (
        args.holdouts
        if args.holdouts is not None
        else list(range(min(args.section_folds, len(sections))))
    )
    if any(holdout >= len(sections) for holdout in holdouts):
        raise SystemExit("holdout index is outside the generated sections")
    runs: list[dict[str, object]] = []
    for factors in args.k_values:
        for seed in args.seeds:
            for holdout in holdouts:
                runs.append(_score_fold(
                    sections,
                    truth_by_section,
                    metadata,
                    holdout=holdout,
                    factors=factors,
                    steps=args.steps,
                    data_seed=seed,
                    optimization_seed=seed + factors * 100 + holdout,
                    split_seed=seed + holdout * 1009,
                    gene_folds=args.gene_folds,
                    null_draws=args.null_draws,
                    inference_steps=args.inference_steps,
                ))
    scores_by_k = {
        factors: [score for run in runs if run["k_model"] == factors for score in run["scores"]]
        for factors in args.k_values
    }
    baseline_scores_by_k = {
        factors: [score for run in runs if run["k_model"] == factors for score in run["baseline_scores"]]
        for factors in args.k_values
    }
    means = {k: float(np.mean(values)) for k, values in scores_by_k.items()}
    ses = {
        k: float(np.std(values, ddof=1) / np.sqrt(len(values))) if len(values) > 1 else 0.0
        for k, values in scores_by_k.items()
    }
    best_k = max(means, key=means.get)
    predictive_eligible = [
        k for k in args.k_values if means[k] >= means[best_k] - ses[best_k]
    ]
    support = {
        k: bool(all(
            run["fit_platform"]["converged"] and run["recovery_pass"]
            for run in runs if run["k_model"] == k
        ))
        for k in args.k_values
    }
    eligible = [k for k in predictive_eligible if support[k]]
    selected = min(eligible) if eligible else None
    boundary = selected == max(args.k_values) if selected is not None else best_k == max(args.k_values)
    final_status = "K_NOT_IDENTIFIABLE"
    if selected is not None and not boundary and len(args.seeds) >= 3:
        final_status = "K_SELECTED"
    elif boundary:
        final_status = "K_SEARCH_BOUNDARY"
    result = {
        "schema": "r04.k_calibration.v2",
        "status": "K_CALIBRATION_COMPLETE",
        "truth": metadata,
        "input_data_hash": sections_content_hash(sections),
        "generator_contract": {
            "generator": "mnsf_correct",
            "domain": "disconnected",
            "sections": 3,
            "spots_per_section": 50,
            "genes": 24,
        },
        "steps": args.steps,
        "inference_steps": args.inference_steps,
        "section_folds": args.section_folds,
        "gene_folds": args.gene_folds,
        "null_draws": args.null_draws,
        "seeds": args.seeds,
        "k_values": args.k_values,
        "holdouts": holdouts,
        "leakage_control": "heldout_section_fields_adapted_on_disjoint_gene_folds",
        "runs": runs,
        "scores_by_k": {str(k): value for k, value in scores_by_k.items()},
        "baseline_scores_by_k": {str(k): value for k, value in baseline_scores_by_k.items()},
        "heldout_score_mean_by_k": {str(k): value for k, value in means.items()},
        "heldout_score_se_by_k": {str(k): value for k, value in ses.items()},
        "support_by_k": {str(k): value for k, value in support.items()},
        "one_se_candidate": {
            "best_predictive_k": best_k,
            "predictive_eligible_k": predictive_eligible,
            "eligible_k": eligible,
            "selected_k": selected,
            "search_boundary": boundary,
        },
        "final_k_status": final_status,
        "formal_restart_authorized": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(args.output, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
