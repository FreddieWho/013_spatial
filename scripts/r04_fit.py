#!/usr/bin/env python3
"""Fit the frozen R-04 dual model on one pure training-role manifest."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy import sparse

from r04.candidates import match_model_factors
from r04.loaders import load_section_from_row
from r04.models import (
    MNSFConfig,
    MNSFEstimator,
    SignedResidualGPConfig,
    SignedResidualGPEstimator,
)
from r04.runtime import atomic_json, resource_status
from r04.serialization import (
    read_json,
    write_candidate_registry,
    write_field_fits,
    write_frozen_model,
)


ROLE_NAMES = (
    "training",
    "internal_validation",
    "external_validation",
    "serial_section_validation",
)


def _model_checkpoint_dir(base: Path | None, model: str) -> str | None:
    if base is None:
        return None
    if model not in {"mnsf", "signed"}:
        raise ValueError(f"unsupported checkpoint model: {model}")
    return str(base / model)


def _manifest_roles(manifest: dict[str, object]) -> tuple[str, ...]:
    rows = manifest.get("rows")
    if not isinstance(rows, list) or not rows:
        raise ValueError("manifest rows are missing or empty")
    roles = tuple(sorted({str(row.get("primary_role", "")) for row in rows if isinstance(row, dict)}))
    if not roles or "" in roles or any(role not in ROLE_NAMES for role in roles):
        raise ValueError(f"manifest has missing or unsupported roles: {roles}")
    return roles


def _align_sections(sections, genes: tuple[str, ...], *, allow_missing: bool = False):
    """Align a frozen panel without zero-filling absent genes.

    Fitting requires complete coverage.  Validation may opt into a reduced
    observed panel; callers must then record the returned panel metadata.
    """
    aligned = []
    coverage = []
    for section in sections:
        positions = {gene: index for index, gene in enumerate(section.gene_id)}
        observed_genes = tuple(gene for gene in genes if gene in positions)
        missing = tuple(gene for gene in genes if gene not in positions)
        if missing and not allow_missing:
            raise ValueError(f"frozen gene list is not fully observed in {section.section_id}")
        if not observed_genes:
            raise ValueError(f"no frozen genes are observed in {section.section_id}")
        columns = [positions[gene] for gene in observed_genes]
        counts = section.counts[:, columns] if sparse.issparse(section.counts) else np.asarray(section.counts)[:, columns]
        aligned.append(type(section)(
            section.section_id,
            section.patient_id,
            section.block_id,
            section.lineage,
            section.barcode,
            section.coords,
            counts,
            observed_genes,
            section.library_size,
            section.metadata,
        ))
        coverage.append({
            "section_id": section.section_id,
            "observed_genes": len(observed_genes),
            "frozen_genes": len(genes),
            "coverage": len(observed_genes) / len(genes),
            "missing_genes": list(missing),
        })
    return aligned, coverage


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument("--gene-list", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--model", choices=("both", "mnsf", "signed"), default="both")
    parser.add_argument("--expected-role", choices=ROLE_NAMES, default="training")
    parser.add_argument("--factors", type=int, default=4)
    parser.add_argument("--inducing-points", type=int, default=64)
    parser.add_argument("--lengthscale", type=float, default=2.0)
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--posterior-draws", type=int, default=200)
    parser.add_argument("--gene-batch-size", type=int, default=None)
    parser.add_argument("--checkpoint-dir", type=Path)
    parser.add_argument("--checkpoint-steps", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260807)
    args = parser.parse_args()
    if args.checkpoint_steps < 1:
        raise SystemExit("checkpoint-steps must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if resource_status(args.output_dir) == "BLOCKED_STORAGE":
        atomic_json(args.output_dir / "r04_run.json", {"status": "BLOCKED_STORAGE"})
        return 2

    try:
        manifest = read_json(args.manifest_json)
        if manifest.get("status") != "READY":
            raise ValueError("manifest_not_ready")
        roles = _manifest_roles(manifest)
        if roles != (args.expected_role,):
            raise ValueError(f"fit requires a pure {args.expected_role} manifest, observed {roles}")
        genes = tuple(line.strip() for line in args.gene_list.read_text(encoding="utf-8").splitlines() if line.strip())
        if not genes:
            raise ValueError("empty_frozen_gene_list")
        sections, coverage = _align_sections(
            [load_section_from_row(row) for row in manifest["rows"]], genes
        )
        mnsf = None
        signed = None
        if args.model in {"both", "mnsf"}:
            mnsf = MNSFEstimator(MNSFConfig(
                factors=args.factors,
                inducing_points=args.inducing_points,
                lengthscale=args.lengthscale,
                steps=args.steps,
                posterior_draws=args.posterior_draws,
                gene_batch_size=args.gene_batch_size,
                checkpoint_dir=_model_checkpoint_dir(args.checkpoint_dir, "mnsf"),
                checkpoint_steps=args.checkpoint_steps,
                seed=args.seed,
            )).fit(sections)
        if args.model in {"both", "signed"}:
            signed = SignedResidualGPEstimator(SignedResidualGPConfig(
                factors=args.factors,
                inducing_points=args.inducing_points,
                lengthscale=args.lengthscale,
                steps=args.steps,
                posterior_draws=args.posterior_draws,
                gene_batch_size=args.gene_batch_size,
                checkpoint_dir=_model_checkpoint_dir(args.checkpoint_dir, "signed"),
                checkpoint_steps=args.checkpoint_steps,
                seed=args.seed + 1000,
            )).fit(sections)
    except Exception as exc:
        atomic_json(args.output_dir / "r04_run.json", {
            "schema": "r04.run.v1",
            "status": "BLOCKED_COMPUTE" if type(exc).__name__ != "ValueError" else "PREFLIGHT_FAILED",
            "error": type(exc).__name__,
            "message": str(exc),
        })
        return 2

    mnsf_fits = [fit for group in (mnsf.fields_ if mnsf else []) for fit in group]
    signed_fits = [fit for group in (signed.fields_ if signed else []) for fit in group]
    candidates = match_model_factors(mnsf_fits, signed_fits)
    if mnsf is not None:
        write_frozen_model(args.output_dir / "mnsf_model.json", mnsf)
    if signed is not None:
        write_frozen_model(args.output_dir / "signed_model.json", signed)
    write_field_fits(args.output_dir / "mnsf_fields.json", mnsf_fits)
    write_field_fits(args.output_dir / "signed_fields.json", signed_fits)
    registry_hash = write_candidate_registry(args.output_dir / "candidate_registry.json", candidates)
    atomic_json(args.output_dir / "r04_run.json", {
        "schema": "r04.run.v1",
        "status": "FIT_COMPLETE_NOT_VALIDATED",
        "candidate_registry_hash": registry_hash,
        "n_sections": len(sections),
        "n_genes": len(genes),
        "model": args.model,
        "expected_role": args.expected_role,
        "observed_roles": list(roles),
        "input_manifest_hash": manifest.get("input_manifest_hash", manifest.get("parent_input_manifest_hash", "")),
        "coverage": coverage,
        "factors": args.factors,
        "inducing_points": args.inducing_points,
        "lengthscale": args.lengthscale,
        "steps": args.steps,
        "posterior_draws": args.posterior_draws,
        "gene_batch_size": args.gene_batch_size,
        "checkpoint_dir": str(args.checkpoint_dir) if args.checkpoint_dir else None,
        "checkpoint_steps": args.checkpoint_steps,
        "seed": args.seed,
        "model_diagnostics": {
            "mnsf": mnsf.diagnostics_ if mnsf else None,
            "signed": signed.diagnostics_ if signed else None,
        },
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
