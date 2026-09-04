#!/usr/bin/env python3
"""Patient-grouped, molecule-only K screening for the real R-04 training set.

This is a discovery-stage screen.  It never reads GT or validation roles.  K=0
uses the same mNSF count model, offset, dispersion and non-spatial nuisance
rank as positive K; it only removes the coordinate-dependent GP component.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import replace
from pathlib import Path

import numpy as np

from r04.diagnostics import platform_convergence_summary
from r04.field_exports import write_heldout_field_export
from r04.loaders import load_section_from_row
from r04.metrics import negative_binomial_log_likelihood
from r04.models import MNSFConfig, MNSFEstimator
from r04.models.mnsf import predict_log_mean_from_fields
from r04.runtime import atomic_json, derived_seed, resource_status
from r04.splits import grouped_fold_ids
from scripts.r04_fit import _align_sections
from scripts.r04_k_calibration import gene_crossfit_splits


def _read_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _subset_genes(section: object, indices: np.ndarray) -> object:
    return replace(
        section,
        counts=section.counts[:, indices],
        gene_id=tuple(section.gene_id[index] for index in indices),
    )


def _counts(section: object) -> np.ndarray:
    return (
        section.counts.toarray()
        if hasattr(section.counts, "toarray")
        else np.asarray(section.counts)
    )


def _library(section: object) -> np.ndarray:
    if section.library_size is not None:
        return np.asarray(section.library_size, dtype=float)
    return _counts(section).sum(axis=1).astype(float)


def _load_training(
    manifest_path: Path,
    gene_path: Path,
    *,
    max_genes: int | None = None,
    max_spots_per_section: int | None = None,
) -> tuple[list[object], tuple[str, ...], str]:
    manifest = _read_json(manifest_path)
    if manifest.get("status") != "READY" or manifest.get("role") != "training":
        raise ValueError("real K search requires a READY pure training manifest")
    genes = tuple(line.strip() for line in gene_path.read_text(encoding="utf-8").splitlines() if line.strip())
    if not genes:
        raise ValueError("frozen gene panel is empty")
    rows = manifest.get("rows")
    if not isinstance(rows, list) or not rows:
        raise ValueError("training manifest has no rows")
    sections = [load_section_from_row(row) for row in rows]
    sections, _ = _align_sections(sections, genes)
    if max_genes is not None:
        if max_genes < 2 or max_genes > len(genes):
            raise ValueError("max-genes must be between 2 and the frozen panel size")
        genes = genes[:max_genes]
        sections, _ = _align_sections(sections, genes)
    if max_spots_per_section is not None:
        if max_spots_per_section < 2:
            raise ValueError("max-spots-per-section must be at least 2")
        sections = [
            replace(
                section,
                barcode=section.barcode[:max_spots_per_section],
                coords=section.coords[:max_spots_per_section],
                counts=section.counts[:max_spots_per_section],
                library_size=(
                    section.library_size[:max_spots_per_section]
                    if section.library_size is not None else None
                ),
            )
            for section in sections
        ]
    patients = {section.patient_id for section in sections}
    if len(patients) < 5:
        raise ValueError("at least five patients are required for grouped K screening")
    return sections, genes, str(manifest.get("input_manifest_hash", ""))


def _grouped_folds(sections: list[object], *, folds: int, seed: int) -> np.ndarray:
    return grouped_fold_ids(
        [section.patient_id for section in sections], n_folds=folds, seed=seed
    )


def _optimization_seed(
    master_seed: int, k_model: int, fold: int, restart_index: int
) -> int:
    if restart_index < 0:
        raise ValueError("restart index must be non-negative")
    if restart_index == 0:
        return derived_seed(master_seed, "real-k", k_model, fold)
    return derived_seed(master_seed, "real-k-restart", restart_index, k_model, fold)


def _parse_fold_values(value: str, available: list[int]) -> list[int]:
    tokens = [token.strip() for token in value.split(",") if token.strip()]
    if not tokens:
        raise ValueError("fold-values must contain at least one integer")
    try:
        selected = sorted({int(token) for token in tokens})
    except ValueError as exc:
        raise ValueError("fold-values must contain only integers") from exc
    missing = sorted(set(selected) - set(available))
    if missing:
        raise ValueError(f"requested fold values not present: {missing}")
    return selected


def _parse_gene_batch_size(value: str) -> int | None:
    normalized = value.strip().lower()
    if normalized in {"full", "all", "none"}:
        return None
    try:
        batch_size = int(normalized)
    except ValueError as exc:
        raise ValueError("gene-batch-size must be a positive integer or full") from exc
    if batch_size < 1:
        raise ValueError("gene-batch-size must be positive or full")
    return batch_size


def _is_wiring_only(
    *,
    fold_values: str | None,
    max_folds: int | None,
    max_genes: int | None,
    max_spots_per_section: int | None,
) -> bool:
    return any(value is not None for value in (
        fold_values,
        max_folds,
        max_genes,
        max_spots_per_section,
    ))


def _mark_wiring_only(result: dict[str, object], wiring_only: bool) -> dict[str, object]:
    if wiring_only:
        result["screen_selected_k"] = result.get("selected_k")
        result["selected_k"] = None
        result["status"] = "WIRING_ONLY_NOT_SCIENTIFIC"
    result["wiring_only"] = wiring_only
    return result


def _combine_inference_platforms(
    summaries: list[dict[str, object]],
) -> dict[str, object]:
    if not summaries:
        raise ValueError("gene cross-fit produced no inference platform summaries")
    split_converged = [item.get("converged") is True for item in summaries]
    return {
        "converged": all(split_converged),
        "split_count": len(summaries),
        "split_converged": split_converged,
    }


def _score_fold(
    training: list[object],
    validation: list[object],
    *,
    genes: tuple[str, ...],
    k_model: int,
    fold: int,
    output_dir: Path,
    steps: int,
    inference_steps: int,
    inducing_points: int,
    lengthscale: float,
    learning_rate: float,
    learning_rate_final: float | None,
    learning_rate_decay_steps: int | None,
    gene_batch_size: int | None,
    diagnostic_interval: int,
    evaluation_interval: int,
    evaluation_mc_draws: int,
    optimization_schedule: str,
    shared_steps: int,
    gene_folds: int,
    split_seed: int,
    optimization_seed: int,
    checkpoint_steps: int,
    resume_checkpoint_dir: Path | None = None,
    environment_hash: str | None = None,
    checkpoint_output_dir: Path | None = None,
    inference_export_root: Path | None = None,
    structure_export_root: Path | None = None,
    export_provenance: dict[str, object] | None = None,
    execution_mode: str = "auto",
    training_only: bool = False,
    export_heldout_effect: bool = True,
) -> dict[str, object]:
    if resource_status(output_dir) == "BLOCKED_STORAGE":
        raise RuntimeError("BLOCKED_STORAGE")
    checkpoint_dir = checkpoint_output_dir or (
        output_dir / "checkpoints" / f"k{k_model}" / f"fold{fold}"
    )
    estimator = MNSFEstimator(MNSFConfig(
        factors=k_model,
        inducing_points=inducing_points,
        nonspatial_rank=1,
        lengthscale=lengthscale,
        learning_rate=learning_rate,
        learning_rate_final=learning_rate_final,
        learning_rate_decay_steps=learning_rate_decay_steps,
        steps=steps,
        posterior_draws=8,
        gene_batch_size=gene_batch_size,
        diagnostic_interval=diagnostic_interval,
        evaluation_interval=evaluation_interval,
        evaluation_mc_draws=evaluation_mc_draws,
        optimization_schedule=optimization_schedule,
        shared_steps=shared_steps,
        checkpoint_dir=str(checkpoint_dir),
        resume_checkpoint_dir=(
            str(resume_checkpoint_dir) if resume_checkpoint_dir is not None else None
        ),
        checkpoint_steps=checkpoint_steps,
        seed=optimization_seed,
        environment_hash=environment_hash or "",
        execution_mode=execution_mode,
    )).fit(training)
    fit_platform = platform_convergence_summary(
        estimator.diagnostics_.get("loss_trace", []), window=50, stable_windows=2
    )
    execution_mode_requested = str(
        estimator.diagnostics_.get("execution_mode_requested", execution_mode)
    )
    execution_mode_effective = str(
        estimator.diagnostics_.get(
            "execution_mode_effective",
            estimator.diagnostics_.get("execution_mode", "eager"),
        )
    )
    score_by_patient: dict[str, list[float]] = {}
    section_scores: list[dict[str, object]] = []
    inference_platforms: list[dict[str, object]] = []
    inference_diagnostics_by_split: list[dict[str, object]] = []
    inference_field_exports: list[dict[str, object]] = []
    structure_field_exports: list[dict[str, object]] = []
    if structure_export_root is not None:
        structure_field_exports.append(write_heldout_field_export(
            structure_export_root / "training_full.npz",
            sections=training,
            field_groups=estimator.fields_,
            evaluation_loading=estimator.loading_,
            factor_amplitude=estimator.factor_amplitude_,
            frozen_gene_ids=genes,
            adaptation_gene_ids=(),
            evaluation_gene_ids=genes,
            adaptation_indices=np.asarray([], dtype=int),
            evaluation_indices=np.arange(len(genes), dtype=int),
            provenance={
                "k_model": k_model,
                "fold": fold,
                "gene_split": None,
                "optimization_seed": optimization_seed,
                "inference_steps": 0,
                "fit_input_hash": estimator.fit_input_hash_,
                "fit_config_hash": estimator.fit_config_hash_,
                "fit_environment_hash": estimator.fit_environment_hash_,
                "fit_updates": estimator.diagnostics_.get("optimizer_steps_this_call"),
                "representation_role": "training_outer_fold_readout_fit",
                **(export_provenance or {}),
            },
        ))
    if training_only:
        return {
            "schema": "r04.real_k_cell.v1",
            "status": "FIT_ONLY_TRAINING_EFFECT_EXPORTED",
            "k_model": k_model,
            "fold": fold,
            "n_training_sections": len(training),
            "n_validation_sections": len(validation),
            "training_patients": sorted({str(section.patient_id) for section in training}),
            "validation_patients": sorted({str(section.patient_id) for section in validation}),
            "fit_platform": fit_platform,
            "execution_mode": execution_mode_effective,
            "execution_mode_requested": execution_mode_requested,
            "execution_mode_effective": execution_mode_effective,
            "inference_platform": {"converged": False, "status": "NOT_RUN_TRAINING_ONLY"},
            "fit_diagnostics": estimator.diagnostics_,
            "inference_diagnostics": {},
            "inference_diagnostics_by_split": [],
            "inference_field_exports": [],
            "structure_field_exports": structure_field_exports,
            "patient_scores": {},
            "section_scores": [],
            "input_hash": estimator.fit_input_hash_,
            "config_hash": estimator.fit_config_hash_,
            "checkpoint_dir": str(checkpoint_dir),
            "parameters": {
                "steps": steps,
                "inference_steps": inference_steps,
                "inducing_points": inducing_points,
                "lengthscale": lengthscale,
                "nonspatial_rank": 1,
                "gene_folds": gene_folds,
                "split_seed": split_seed,
                "optimization_seed": optimization_seed,
                "optimization_schedule": optimization_schedule,
                "shared_steps": shared_steps,
                "training_only": True,
                "execution_mode": execution_mode_effective,
                "execution_mode_requested": execution_mode_requested,
                "execution_mode_effective": execution_mode_effective,
            },
        }
    splits = gene_crossfit_splits(len(genes), folds=gene_folds, seed=split_seed)
    for split_index, (adaptation, evaluation) in enumerate(splits):
        adaptation_genes = tuple(genes[index] for index in adaptation)
        adaptation_sections = [_subset_genes(section, adaptation) for section in validation]
        projected = estimator.subset_to_genes(adaptation_genes)
        inferred = projected.infer(
            adaptation_sections,
            steps=inference_steps,
            posterior_draws=8,
        )
        loading = estimator.loading_[evaluation]
        baseline = estimator.gene_baseline_[evaluation]
        nuisance_loading = (
            estimator.nonspatial_loading_[evaluation]
            if estimator.nonspatial_loading_ is not None else None
        )
        log_mu = predict_log_mean_from_fields(
            adaptation_sections,
            inferred,
            loading,
            estimator.factor_amplitude_,
            baseline,
            nonspatial_component=projected.inference_nonspatial_component_,
            nonspatial_loading=nuisance_loading,
        )
        inference_platform = platform_convergence_summary(
            projected.inference_diagnostics_.get("loss_trace", []),
            window=50,
            stable_windows=2,
        )
        inference_platforms.append(inference_platform)
        inference_diagnostics_by_split.append({
            "gene_split": split_index,
            "adaptation_gene_count": int(len(adaptation)),
            "evaluation_gene_count": int(len(evaluation)),
            "platform": inference_platform,
            "diagnostics": projected.inference_diagnostics_,
        })
        if inference_export_root is not None:
            export_path = inference_export_root / f"gene_split_{split_index}.npz"
            inference_field_exports.append(write_heldout_field_export(
                export_path,
                sections=validation,
                field_groups=inferred,
                evaluation_loading=estimator.loading_[evaluation],
                factor_amplitude=projected.factor_amplitude_,
                frozen_gene_ids=genes,
                adaptation_gene_ids=adaptation_genes,
                evaluation_gene_ids=tuple(genes[index] for index in evaluation),
                adaptation_indices=adaptation,
                evaluation_indices=evaluation,
                provenance={
                    "k_model": k_model,
                    "fold": fold,
                    "gene_split": split_index,
                    "gene_split_seed": split_seed,
                    "optimization_seed": optimization_seed,
                    "inference_steps": inference_steps,
                    "fit_input_hash": estimator.fit_input_hash_,
                    "fit_config_hash": estimator.fit_config_hash_,
                    "fit_environment_hash": estimator.fit_environment_hash_,
                    "inference_platform": inference_platform,
                    "representation_role": "heldout_training_outer_fold",
                    **(export_provenance or {}),
                },
            ))
        offset = 0
        for section in validation:
            n_spots = len(section.coords)
            score = negative_binomial_log_likelihood(
                _counts(section)[:, evaluation],
                log_mu[offset : offset + n_spots],
                estimator.dispersion_[evaluation],
            )
            patient = str(section.patient_id)
            score_by_patient.setdefault(patient, []).append(score)
            section_scores.append({
                "section_id": section.section_id,
                "patient_id": patient,
                "gene_split": split_index,
                "score": float(score),
                "evaluation_gene_count": int(len(evaluation)),
            })
            offset += n_spots
        if offset != len(log_mu):
            raise RuntimeError("validation prediction spot count does not match sections")
    if structure_export_root is not None and export_heldout_effect:
        heldout_fields = estimator.infer(
            validation,
            steps=inference_steps,
            posterior_draws=8,
        )
        structure_inference_platform = platform_convergence_summary(
            estimator.inference_diagnostics_.get("loss_trace", []),
            window=50,
            stable_windows=2,
        )
        structure_field_exports.append(write_heldout_field_export(
            structure_export_root / "heldout_full.npz",
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
                "fold": fold,
                "gene_split": None,
                "optimization_seed": optimization_seed,
                "inference_steps": inference_steps,
                "fit_input_hash": estimator.fit_input_hash_,
                "fit_config_hash": estimator.fit_config_hash_,
                "fit_environment_hash": estimator.fit_environment_hash_,
                "inference_platform": structure_inference_platform,
                "representation_role": "heldout_outer_fold_readout_apply",
                **(export_provenance or {}),
            },
        ))
    patient_scores = {
        patient: float(np.mean(values))
        for patient, values in score_by_patient.items()
    }
    return {
        "schema": "r04.real_k_cell.v1",
        "status": "FIT_AND_SCORED",
        "k_model": k_model,
        "fold": fold,
        "n_training_sections": len(training),
        "n_validation_sections": len(validation),
        "training_patients": sorted({str(section.patient_id) for section in training}),
        "validation_patients": sorted({str(section.patient_id) for section in validation}),
        "fit_platform": fit_platform,
        "execution_mode": execution_mode_effective,
        "execution_mode_requested": execution_mode_requested,
        "execution_mode_effective": execution_mode_effective,
        "inference_platform": _combine_inference_platforms(inference_platforms),
        "fit_diagnostics": estimator.diagnostics_,
        "inference_diagnostics": projected.inference_diagnostics_,
        "inference_diagnostics_by_split": inference_diagnostics_by_split,
        "inference_field_exports": inference_field_exports,
        "structure_field_exports": structure_field_exports,
        "patient_scores": patient_scores,
        "section_scores": section_scores,
        "input_hash": estimator.fit_input_hash_,
        "config_hash": estimator.fit_config_hash_,
        "checkpoint_dir": str(checkpoint_dir),
        "parameters": {
            "steps": steps,
            "inference_steps": inference_steps,
            "inducing_points": inducing_points,
            "lengthscale": lengthscale,
            "nonspatial_rank": 1,
            "gene_folds": gene_folds,
            "split_seed": split_seed,
            "optimization_seed": optimization_seed,
            "optimization_schedule": optimization_schedule,
            "shared_steps": shared_steps,
            "execution_mode": execution_mode_effective,
            "execution_mode_requested": execution_mode_requested,
            "execution_mode_effective": execution_mode_effective,
        },
    }


def _bootstrap_ci(values: np.ndarray, *, seed: int, draws: int = 2000) -> tuple[float, float]:
    if len(values) < 2:
        raise ValueError("grouped bootstrap requires at least two patients")
    rng = np.random.default_rng(seed)
    sampled = values[rng.integers(0, len(values), size=(draws, len(values)))]
    means = sampled.mean(axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def _aggregate(cells: list[dict[str, object]], k_values: list[int], *, seed: int) -> dict[str, object]:
    by_k: dict[int, dict[str, float]] = {}
    for cell in cells:
        if cell.get("status") != "FIT_AND_SCORED":
            raise ValueError("cannot aggregate a failed or missing K cell")
        k = int(cell["k_model"])
        patients = cell["patient_scores"]
        if not isinstance(patients, dict):
            raise ValueError("cell patient scores are missing")
        for patient, score in patients.items():
            if patient in by_k.setdefault(k, {}):
                raise ValueError(f"duplicate patient score for K={k}: {patient}")
            by_k[k][str(patient)] = float(score)
    if set(by_k) != set(k_values):
        raise ValueError("K cells are incomplete")
    support_by_k: dict[int, bool] = {}
    for k in k_values:
        k_cells = [cell for cell in cells if int(cell["k_model"]) == k]
        has_platform = all(
            "fit_platform" in cell and "inference_platform" in cell
            for cell in k_cells
        )
        support_by_k[k] = (
            all(bool(cell["fit_platform"].get("converged")) for cell in k_cells)
            and all(bool(cell["inference_platform"].get("converged")) for cell in k_cells)
            if has_platform else True
        )
    means = {k: float(np.mean(list(by_k[k].values()))) for k in k_values}
    ses = {
        k: float(np.std(list(by_k[k].values()), ddof=1) / np.sqrt(len(by_k[k])))
        for k in k_values
    }
    k0 = by_k[0]
    delta_ci: dict[int, dict[str, float]] = {}
    for k in k_values:
        common = sorted(set(k0) & set(by_k[k]))
        deltas = np.asarray([by_k[k][patient] - k0[patient] for patient in common])
        low, high = _bootstrap_ci(deltas, seed=seed + k)
        delta_ci[k] = {
            "mean": float(np.mean(deltas)),
            "ci_low": low,
            "ci_high": high,
            "patients": len(common),
        }
    best = max(k_values, key=means.get)
    predictive_eligible = [k for k in k_values if means[k] >= means[best] - ses[best]]
    selected = min(predictive_eligible)
    if selected > 0 and delta_ci[selected]["ci_low"] <= 0:
        selected = 0
    if not all(support_by_k.values()):
        selected = None
        status = "K_NOT_IDENTIFIABLE_NOT_CONVERGED"
        boundary = False
    else:
        boundary = best == max(k_values) or selected == max(k_values)
        if boundary:
            status = "K_SEARCH_BOUNDARY"
        elif selected == 0:
            status = "K0_SELECTED_REAL"
        else:
            status = "K_SELECTED_REAL_ROUGH_SCREEN"
    next_k_values = (
        [max(k_values) + 1, max(k_values) + 2]
        if boundary else []
    )
    return {
        "schema": "r04.real_k_aggregate.v1",
        "status": status,
        "k_values": k_values,
        "patient_count": len(by_k[0]),
        "mean_score_by_k": {str(k): means[k] for k in k_values},
        "score_se_by_k": {str(k): ses[k] for k in k_values},
        "support_by_k": {str(k): support_by_k[k] for k in k_values},
        "delta_to_k0_by_k": {str(k): value for k, value in delta_ci.items()},
        "best_predictive_k": best,
        "predictive_eligible_k": predictive_eligible,
        "selected_k": selected,
        "search_boundary": boundary,
        "next_k_values": next_k_values,
        "null_status": "NOT_RUN",
        "formal_candidate_status": "NOT_AUTHORIZED_UNTIL_SPATIAL_NULL_AND_RESTART",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument("--gene-list", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--k-values", default="0,1,2,3,4")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--gene-folds", type=int, default=2)
    parser.add_argument("--group-seed", type=int, default=20260807)
    parser.add_argument("--split-seed", type=int, default=20260817)
    parser.add_argument("--restart-index", type=int, default=0)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--inference-steps", type=int, default=100)
    parser.add_argument("--checkpoint-steps", type=int, default=100)
    parser.add_argument("--inducing-points", type=int, default=16)
    parser.add_argument("--lengthscale", type=float, default=3.0)
    parser.add_argument("--learning-rate", type=float, default=0.05)
    parser.add_argument("--learning-rate-final", type=float, default=None)
    parser.add_argument("--learning-rate-decay-steps", type=int, default=None)
    parser.add_argument(
        "--gene-batch-size",
        default="256",
        help="positive integer, or full/all for the complete gene panel",
    )
    parser.add_argument("--diagnostic-interval", type=int, default=100)
    parser.add_argument("--evaluation-interval", type=int, default=0)
    parser.add_argument("--evaluation-mc-draws", type=int, default=4)
    parser.add_argument(
        "--optimization-schedule",
        choices=("joint", "staged_shared_first"),
        default="joint",
        help="joint baseline or shared-parameter warmup followed by joint fitting",
    )
    parser.add_argument("--shared-steps", type=int, default=0)
    parser.add_argument("--fold-values", default=None, help="comma-separated explicit fold values; wiring-only")
    parser.add_argument("--max-folds", type=int, default=None, help="wiring-only limit; output is not scientific")
    parser.add_argument("--max-genes", type=int, default=None, help="wiring-only panel trim")
    parser.add_argument("--max-spots-per-section", type=int, default=None, help="wiring-only spot trim")
    args = parser.parse_args()
    k_values = sorted({int(value) for value in args.k_values.split(",") if value.strip()})
    if not k_values or k_values[0] != 0 or any(k < 0 for k in k_values):
        raise SystemExit("k-values must include non-negative K=0")
    if args.folds < 2 or args.gene_folds < 2:
        raise SystemExit("folds and gene-folds must be at least 2")
    if args.fold_values is not None and args.max_folds is not None:
        raise SystemExit("fold-values and max-folds are mutually exclusive")
    if args.learning_rate <= 0:
        raise SystemExit("learning-rate must be positive")
    if args.restart_index < 0:
        raise SystemExit("restart-index must be non-negative")
    if args.learning_rate_final is not None and args.learning_rate_final <= 0:
        raise SystemExit("learning-rate-final must be positive")
    if args.learning_rate_final is not None and (
        args.learning_rate_decay_steps is None or args.learning_rate_decay_steps < 1
    ):
        raise SystemExit(
            "learning-rate-decay-steps is required when learning-rate-final is set"
        )
    if args.learning_rate_final is None and args.learning_rate_decay_steps is not None:
        raise SystemExit(
            "learning-rate-final is required when learning-rate-decay-steps is set"
        )
    if args.diagnostic_interval < 0:
        raise SystemExit("diagnostic-interval must be non-negative")
    if args.shared_steps < 0:
        raise SystemExit("shared-steps must be non-negative")
    if args.evaluation_interval < 0:
        raise SystemExit("evaluation-interval must be non-negative")
    if args.evaluation_mc_draws < 1:
        raise SystemExit("evaluation-mc-draws must be positive")
    try:
        gene_batch_size = _parse_gene_batch_size(args.gene_batch_size)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    wiring_only = _is_wiring_only(
        fold_values=args.fold_values,
        max_folds=args.max_folds,
        max_genes=args.max_genes,
        max_spots_per_section=args.max_spots_per_section,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if resource_status(args.output_dir) == "BLOCKED_STORAGE":
        atomic_json(args.output_dir / "r04_real_k_search.json", {"status": "BLOCKED_STORAGE"})
        return 2
    atomic_json(args.output_dir / "r04_real_k_search.json", {
        "schema": "r04.real_k_search.v1",
        "status": "RUNNING",
        "pid": os.getpid(),
        "k_values": k_values,
        "folds": args.folds,
        "restart_index": args.restart_index,
        "steps": args.steps,
        "inference_steps": args.inference_steps,
        "learning_rate": args.learning_rate,
        "learning_rate_final": args.learning_rate_final,
        "learning_rate_decay_steps": args.learning_rate_decay_steps,
        "gene_batch_size": gene_batch_size,
        "diagnostic_interval": args.diagnostic_interval,
        "evaluation_interval": args.evaluation_interval,
        "evaluation_mc_draws": args.evaluation_mc_draws,
        "optimization_schedule": args.optimization_schedule,
        "shared_steps": args.shared_steps,
        "fold_values": args.fold_values,
        "wiring_only": wiring_only,
    })
    sections, genes, manifest_hash = _load_training(
        args.manifest_json,
        args.gene_list,
        max_genes=args.max_genes,
        max_spots_per_section=args.max_spots_per_section,
    )
    fold_ids = _grouped_folds(sections, folds=args.folds, seed=args.group_seed)
    available_fold_values = sorted(set(int(value) for value in fold_ids))
    if args.fold_values is not None:
        try:
            fold_values = _parse_fold_values(args.fold_values, available_fold_values)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
    else:
        fold_values = available_fold_values
    if args.max_folds is not None:
        fold_values = fold_values[: args.max_folds]
    cells: list[dict[str, object]] = []
    for k_model in k_values:
        for fold in fold_values:
            training = [section for index, section in enumerate(sections) if fold_ids[index] != fold]
            validation = [section for index, section in enumerate(sections) if fold_ids[index] == fold]
            cell_dir = args.output_dir / "cells" / f"k{k_model}" / f"fold{fold}"
            cell_dir.mkdir(parents=True, exist_ok=True)
            seed = _optimization_seed(
                args.group_seed, k_model, fold, args.restart_index
            )
            try:
                cell = _score_fold(
                    training,
                    validation,
                    genes=genes,
                    k_model=k_model,
                    fold=fold,
                    output_dir=args.output_dir,
                    steps=args.steps,
                    inference_steps=args.inference_steps,
                    inducing_points=args.inducing_points,
                    lengthscale=args.lengthscale,
                    learning_rate=args.learning_rate,
                    learning_rate_final=args.learning_rate_final,
                    learning_rate_decay_steps=args.learning_rate_decay_steps,
                    gene_batch_size=gene_batch_size,
                    diagnostic_interval=args.diagnostic_interval,
                    evaluation_interval=args.evaluation_interval,
                    evaluation_mc_draws=args.evaluation_mc_draws,
                    optimization_schedule=args.optimization_schedule,
                    shared_steps=args.shared_steps,
                    gene_folds=args.gene_folds,
                    split_seed=args.split_seed,
                    optimization_seed=seed,
                    checkpoint_steps=args.checkpoint_steps,
                )
            except Exception as exc:
                cell = {
                    "schema": "r04.real_k_cell.v1",
                    "status": "FAILED",
                    "k_model": k_model,
                    "fold": fold,
                    "error": type(exc).__name__,
                    "message": str(exc),
                }
                atomic_json(cell_dir / "cell.json", cell)
                atomic_json(args.output_dir / "r04_real_k_search.json", {
                    "schema": "r04.real_k_search.v1",
                    "status": "FAILED",
                    "pid": os.getpid(),
                    "failed_cell": {"k_model": k_model, "fold": fold},
                    "error": type(exc).__name__,
                    "message": str(exc),
                })
                raise
            atomic_json(cell_dir / "cell.json", cell)
            cells.append(cell)
    aggregate = _aggregate(cells, k_values, seed=args.group_seed)
    result = {
        **aggregate,
        "schema": "r04.real_k_search.v1",
        "manifest_hash": manifest_hash,
        "n_sections": len(sections),
        "n_patients": len({str(section.patient_id) for section in sections}),
        "folds": args.folds,
        "observed_folds": fold_values,
        "group_seed": args.group_seed,
        "restart_index": args.restart_index,
        "split_seed": args.split_seed,
        "steps": args.steps,
        "inference_steps": args.inference_steps,
        "gene_folds": args.gene_folds,
        "inducing_points": args.inducing_points,
        "lengthscale": args.lengthscale,
        "learning_rate": args.learning_rate,
        "learning_rate_final": args.learning_rate_final,
        "learning_rate_decay_steps": args.learning_rate_decay_steps,
        "optimization_schedule": args.optimization_schedule,
        "shared_steps": args.shared_steps,
        "gene_batch_size": gene_batch_size,
        "diagnostic_interval": args.diagnostic_interval,
        "evaluation_interval": args.evaluation_interval,
        "evaluation_mc_draws": args.evaluation_mc_draws,
        "fold_values": args.fold_values,
        "nonspatial_rank": 1,
        "cell_count": len(cells),
    }
    result = _mark_wiring_only(result, wiring_only)
    atomic_json(args.output_dir / "r04_real_k_search.json", result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
