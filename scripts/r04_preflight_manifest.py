#!/usr/bin/env python3
"""Preflight all R-04 molecule-only sections and freeze training genes."""

from __future__ import annotations

import argparse
import json
import tempfile
from collections import Counter
from dataclasses import asdict
from pathlib import Path

import numpy as np

from r04.gene_selection import select_training_genes
from r04.io_contract import section_content_hash
from r04.loaders import load_section_from_row
from r04.runtime import atomic_json, resource_status


def _write_gene_list(path: Path, genes: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        handle.write("\n".join(genes) + "\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument("--gene-list", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--n-genes", type=int, default=4000)
    args = parser.parse_args()

    if resource_status(args.project_root) == "BLOCKED_STORAGE":
        atomic_json(args.output, {"status": "BLOCKED_STORAGE"})
        return 2
    manifest = json.loads(args.manifest_json.read_text(encoding="utf-8"))
    if manifest.get("status") != "READY":
        atomic_json(args.output, {"status": "BLOCKED_INPUT_CONTRACT", "reason": "manifest_not_ready"})
        return 2

    sections = []
    section_reports = []
    errors = []
    for row in manifest["rows"]:
        try:
            section = load_section_from_row(row)
            libraries = np.asarray(section.counts.sum(axis=1)).ravel()
            if len(set(section.barcode)) != len(section.barcode):
                raise ValueError("duplicate barcodes")
            if len(set(section.gene_id)) != len(section.gene_id):
                raise ValueError("duplicate gene IDs")
            if len(np.unique(section.coords, axis=0)) != len(section.coords):
                raise ValueError("duplicate spatial coordinates")
            if np.any(libraries <= 0):
                raise ValueError("non-positive library size")
            sections.append(section)
            section_reports.append({
                "section_id": section.section_id,
                "patient_id": section.patient_id,
                "lineage": section.lineage,
                "primary_role": row.get("primary_role", ""),
                "n_spots": int(section.counts.shape[0]),
                "zero_library_spots_excluded": int(section.metadata.get("zero_library_spots_excluded", "0")),
                "n_genes": int(section.counts.shape[1]),
                "library_min": float(libraries.min()),
                "library_median": float(np.median(libraries)),
                "library_max": float(libraries.max()),
                "gene_id_sha256": section_content_hash(section),
            })
        except Exception as exc:
            errors.append({"section_id": row.get("section_id", ""), "error": type(exc).__name__, "message": str(exc)})

    if errors:
        atomic_json(args.output, {
            "schema": "r04.preflight.v1",
            "status": "BLOCKED_INPUT_CONTRACT",
            "errors": errors,
            "section_reports": section_reports,
        })
        return 2

    training = [section for section in sections if next(row for row in manifest["rows"] if row["section_id"] == section.section_id).get("primary_role") == "training"]
    report = select_training_genes(training, n_genes=args.n_genes)
    _write_gene_list(args.gene_list, report.gene_id)
    selected = set(report.gene_id)
    lineage_coverage = {}
    for section in sections:
        coverage = len(selected.intersection(section.gene_id)) / len(selected)
        lineage_coverage.setdefault(section.lineage, []).append(coverage)
    coverage_summary = {lineage: {"min": min(values), "mean": float(np.mean(values)), "n_sections": len(values)} for lineage, values in lineage_coverage.items()}
    low_coverage = {lineage: value for lineage, value in coverage_summary.items() if value["min"] < 0.8}
    status = "BLOCKED_GENE_COVERAGE" if low_coverage else "READY"
    atomic_json(args.output, {
        "schema": "r04.preflight.v1",
        "status": status,
        "manifest_hash": manifest.get("input_manifest_hash"),
        "n_sections": len(sections),
        "role_counts": dict(Counter(row.get("primary_role", "") for row in manifest["rows"])),
        "lineage_counts": dict(Counter(section.lineage for section in sections)),
        "gene_selection": asdict(report),
        "lineage_gene_coverage": coverage_summary,
        "low_coverage_lineages": low_coverage,
        "section_reports": section_reports,
        "zero_library_spots_excluded_total": sum(report["zero_library_spots_excluded"] for report in section_reports),
    })
    return 0 if status == "READY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
