#!/usr/bin/env python3
"""Run bounded planted-factor diagnostics before any real-data restart."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from r04.models import MNSFConfig, MNSFEstimator
from r04.runtime import atomic_json
from r04.synthetic import make_overlapping_sections


def _matrix(groups: list, factors: int) -> np.ndarray:
    return np.column_stack([
        np.concatenate([group[factor].field_mean for group in groups])
        for factor in range(factors)
    ])


def _match(truth: np.ndarray, estimate: np.ndarray) -> list[float]:
    truth = truth - truth.mean(axis=0, keepdims=True)
    estimate = estimate - estimate.mean(axis=0, keepdims=True)
    corr = np.corrcoef(truth.T, estimate.T)[: truth.shape[1], truth.shape[1] :]
    rows, columns = linear_sum_assignment(-np.abs(corr))
    return [float(abs(corr[row, column])) for row, column in zip(rows, columns)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--steps", type=int, default=20)
    args = parser.parse_args()
    seeds = (20260807, 20260817, 20260827)
    results = []
    for seed in seeds:
        sections, truth = make_overlapping_sections(
            sections=3, spots_per_section=60, genes=24, antagonistic=True, seed=seed
        )
        estimator = MNSFEstimator(MNSFConfig(
            factors=2, inducing_points=12, lengthscale=1.0,
            steps=args.steps, posterior_draws=4, gene_batch_size=None, seed=seed,
        )).fit(sections)
        results.append({
            "seed": seed,
            "matched_field_absolute_correlation": _match(truth, _matrix(estimator.fields_, 2)),
            "diagnostics": estimator.diagnostics_,
        })
    atomic_json(args.output, {
        "schema": "r04.planted_calibration.v1",
        "status": "CALIBRATION_COMPLETE_NOT_VALIDATED",
        "steps": args.steps,
        "seeds": list(seeds),
        "n_factors": 2,
        "results": results,
        "formal_restart_authorized": False,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
