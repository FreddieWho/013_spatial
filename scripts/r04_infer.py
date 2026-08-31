#!/usr/bin/env python3
"""Infer continuous fields for validation sections from a frozen R-04 model."""

from __future__ import annotations

import argparse
from pathlib import Path

from r04.loaders import load_section_from_row
from r04.io_contract import stable_json_hash
from r04.runtime import atomic_json, resource_status
from r04.serialization import read_frozen_model, read_json, write_field_fits
from scripts.r04_fit import ROLE_NAMES, _align_sections, _manifest_roles


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-json", type=Path, required=True)
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument("--gene-list", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--expected-role",
        choices=ROLE_NAMES,
        action="append",
        default=None,
        help="allowed validation role; repeat for multiple pure validation partitions",
    )
    parser.add_argument("--steps", type=int, default=None)
    parser.add_argument("--posterior-draws", type=int, default=None)
    parser.add_argument(
        "--allow-missing-genes",
        action="store_true",
        help="run conditional inference on the common observed panel; never zero-fills",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if resource_status(args.output_dir) == "BLOCKED_STORAGE":
        atomic_json(args.output_dir / "r04_infer.json", {"status": "BLOCKED_STORAGE"})
        return 2

    expected_roles = tuple(args.expected_role or (
        "internal_validation",
        "external_validation",
        "serial_section_validation",
    ))
    try:
        manifest = read_json(args.manifest_json)
        if manifest.get("status") != "READY":
            raise ValueError("manifest_not_ready")
        roles = _manifest_roles(manifest)
        if any(role not in expected_roles for role in roles):
            raise ValueError(f"inference role {roles} is outside expected roles {expected_roles}")
        genes = tuple(
            line.strip()
            for line in args.gene_list.read_text(encoding="utf-8").splitlines()
            if line.strip()
        )
        if not genes:
            raise ValueError("empty_frozen_gene_list")
        model = read_frozen_model(args.model_json)
        model_genes = tuple(getattr(model, "gene_id_", ()))
        if model_genes != genes:
            raise ValueError("gene list does not match the frozen model artifact")
        raw_sections = [load_section_from_row(row) for row in manifest["rows"]]
        observed_sets = [set(section.gene_id) for section in raw_sections]
        common_genes = tuple(gene for gene in genes if all(gene in observed for observed in observed_sets))
        missing_any = len(common_genes) != len(genes)
        if missing_any and not args.allow_missing_genes:
            raise ValueError("validation has missing frozen genes; pass --allow-missing-genes for conditional inference")
        if not common_genes:
            raise ValueError("validation has no common observed frozen genes")
        coverage = [
            {
                "section_id": section.section_id,
                "observed_genes": len(set(section.gene_id) & set(genes)),
                "frozen_genes": len(genes),
                "coverage": len(set(section.gene_id) & set(genes)) / len(genes),
                "missing_genes": [gene for gene in genes if gene not in set(section.gene_id)],
            }
            for section in raw_sections
        ]
        sections, _ = _align_sections(raw_sections, common_genes)
        panel_diagnostics = {}
        if missing_any:
            model = model.subset_to_genes(common_genes)
            panel_diagnostics = dict(getattr(model, "diagnostics_", {}))
            retained = panel_diagnostics.get("retained_loading_energy", [])
            if any(float(value) < 0.90 for value in retained) or panel_diagnostics.get("conditional_rank_loss", False):
                raise ValueError("conditional_panel_not_testable: loading energy or rank criterion failed")
        inferred = model.infer(
            sections,
            steps=args.steps,
            posterior_draws=args.posterior_draws,
        )
    except Exception as exc:
        atomic_json(args.output_dir / "r04_infer.json", {
            "schema": "r04.infer.v1",
            "status": "BLOCKED_INFERENCE" if type(exc).__name__ != "ValueError" else "PREFLIGHT_FAILED",
            "error": type(exc).__name__,
            "message": str(exc),
        })
        return 2

    fits = [fit for group in inferred for fit in group]
    write_field_fits(args.output_dir / "fields.json", fits)
    atomic_json(args.output_dir / "r04_infer.json", {
        "schema": "r04.infer.v1",
        "status": "INFERENCE_COMPLETE_NOT_VALIDATED",
        "model_id": model.model_id,
        "fit_input_hash": getattr(model, "fit_input_hash_", ""),
        "input_manifest_hash": manifest.get("input_manifest_hash", manifest.get("parent_input_manifest_hash", "")),
        "expected_roles": list(expected_roles),
        "observed_roles": list(roles),
        "n_sections": len(sections),
        "n_genes": len(genes),
        "n_fields": len(fits),
        "coverage": coverage,
        "conditional_panel": bool(len(set(genes) - set(sections[0].gene_id))),
        "conditional_panel_hash": stable_json_hash(list(sections[0].gene_id)),
        "conditional_panel_diagnostics": panel_diagnostics,
        "offset_source": "full_library_size_metadata" if all(section.library_size is not None for section in sections) else "observed_gene_sum",
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
