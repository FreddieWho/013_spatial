#!/usr/bin/env python3
"""Run the small deterministic R-04 dual-model smoke test."""

from __future__ import annotations

import argparse
from pathlib import Path

from r04.candidates import match_model_factors
from r04.models import MNSFConfig, MNSFEstimator, SignedResidualGPConfig, SignedResidualGPEstimator
from r04.runtime import atomic_json
from r04.serialization import write_candidate_registry, write_field_fits
from r04.synthetic import make_overlapping_sections


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--steps", type=int, default=30)
    parser.add_argument("--factors", type=int, default=2)
    parser.add_argument("--inducing-points", type=int, default=16)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sections, truth = make_overlapping_sections(sections=4, spots_per_section=36, genes=32, antagonistic=True)
    try:
        mnsf = MNSFEstimator(MNSFConfig(factors=args.factors, inducing_points=args.inducing_points, steps=args.steps, posterior_draws=20)).fit(sections)
        signed = SignedResidualGPEstimator(SignedResidualGPConfig(factors=args.factors, inducing_points=args.inducing_points, steps=args.steps, posterior_draws=20)).fit(sections)
    except Exception as exc:  # status is explicit; no fallback to clustering
        atomic_json(args.output_dir / "r04_smoke.json", {"status": "BLOCKED_COMPUTE", "error": type(exc).__name__, "message": str(exc)})
        return 2
    mnsf_fits = [fit for section in mnsf.fields_ for fit in section]
    signed_fits = [fit for section in signed.fields_ for fit in section]
    candidates = match_model_factors(mnsf_fits, signed_fits)
    write_field_fits(args.output_dir / "mnsf_fields.json", mnsf_fits)
    write_field_fits(args.output_dir / "signed_fields.json", signed_fits)
    registry_hash = write_candidate_registry(args.output_dir / "candidate_registry.json", candidates)
    atomic_json(args.output_dir / "r04_smoke.json", {
        "status": "PASS_SYNTHETIC_SMOKE",
        "candidate_registry_hash": registry_hash,
        "n_candidates": len(candidates),
        "truth_shape": list(truth.shape),
        "mnsf_diagnostics": mnsf.diagnostics_,
        "signed_diagnostics": signed.diagnostics_,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
