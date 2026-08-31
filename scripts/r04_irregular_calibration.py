#!/usr/bin/env python3
"""Calibrate continuous-field recovery on irregular synthetic supports."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from r04.models import MNSFConfig, MNSFEstimator, SignedResidualGPConfig, SignedResidualGPEstimator
from r04.runtime import atomic_json
from r04.serialization import write_field_fits
from r04.synthetic import make_irregular_continuous_field_sections


def _field_matrix(groups, factors: int) -> np.ndarray:
    return np.column_stack([
        np.concatenate([group[factor].field_mean for group in groups])
        for factor in range(factors)
    ])


def _matching(truth: np.ndarray, estimate: np.ndarray) -> dict[str, object]:
    truth = truth - truth.mean(axis=0, keepdims=True)
    estimate = estimate - estimate.mean(axis=0, keepdims=True)
    correlation = np.corrcoef(truth.T, estimate.T)[: truth.shape[1], truth.shape[1] :]
    rows, columns = linear_sum_assignment(-np.abs(correlation))
    matched = [float(abs(correlation[row, column])) for row, column in zip(rows, columns)]
    return {
        "absolute_correlation": correlation.tolist(),
        "matched_factor_correlations": matched,
        "minimum_matched_absolute_correlation": min(matched) if matched else 0.0,
        "mean_matched_absolute_correlation": float(np.mean(matched)) if matched else 0.0,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--domain", choices=("crescent", "branch", "disconnected"), default="crescent")
    parser.add_argument("--steps", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260807)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sections, truth = make_irregular_continuous_field_sections(
        domain=args.domain, sections=3, spots_per_section=80, genes=24, seed=args.seed
    )
    mnsf = MNSFEstimator(MNSFConfig(
        factors=2, inducing_points=16, lengthscale=0.5, steps=args.steps,
        posterior_draws=12, seed=args.seed,
    )).fit(sections)
    signed = SignedResidualGPEstimator(SignedResidualGPConfig(
        factors=2, inducing_points=16, lengthscale=0.5, steps=args.steps,
        posterior_draws=12, seed=args.seed + 1000,
    )).fit(sections)
    mnsf_fields = [fit for group in mnsf.fields_ for fit in group]
    signed_fields = [fit for group in signed.fields_ for fit in group]
    write_field_fits(args.output_dir / "mnsf_fields.json", mnsf_fields)
    write_field_fits(args.output_dir / "signed_fields.json", signed_fields)
    atomic_json(args.output_dir / "irregular_calibration.json", {
        "schema": "r04.irregular_calibration.v1",
        "status": "CALIBRATION_COMPLETE_NOT_VALIDATED",
        "domain": args.domain,
        "n_sections": len(sections),
        "n_spots": int(len(truth)),
        "n_genes": len(sections[0].gene_id),
        "steps": args.steps,
        "seed": args.seed,
        "field_definition": "explicit smooth coordinate functions on an irregular support",
        "mnsf": {"matching": _matching(truth, _field_matrix(mnsf.fields_, 2)), "diagnostics": mnsf.diagnostics_},
        "signed": {"matching": _matching(truth, _field_matrix(signed.fields_, 2)), "diagnostics": signed.diagnostics_},
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
