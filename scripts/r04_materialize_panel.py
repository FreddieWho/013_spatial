#!/usr/bin/env python3
"""Materialize a frozen gene panel once for repeated R-04 CPU runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from scipy import sparse

from r04.io_contract import stable_json_hash
from r04.loaders import load_section_from_row
from r04.runtime import atomic_json, resource_status


ROLE_NAMES = (
    "training",
    "internal_validation",
    "external_validation",
    "serial_section_validation",
)


def materialize_panel(
    manifest: dict[str, object],
    genes: tuple[str, ...],
    output_dir: Path,
    *,
    expected_role: str,
    start_row: int = 0,
    end_row: int | None = None,
) -> dict[str, object]:
    if manifest.get("status") != "READY":
        raise ValueError("panel materialization requires a READY manifest")
    rows = manifest.get("rows")
    if not isinstance(rows, list) or not rows:
        raise ValueError("manifest rows are missing or empty")
    roles = {str(row.get("primary_role", "")) for row in rows if isinstance(row, dict)}
    if roles != {expected_role}:
        raise ValueError(f"panel materialization requires a pure {expected_role} manifest, observed {sorted(roles)}")
    if not genes or len(set(genes)) != len(genes):
        raise ValueError("frozen gene panel must be non-empty and unique")
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = output_dir / "sections"
    cache_dir.mkdir(parents=True, exist_ok=True)
    if start_row < 0 or start_row >= len(rows):
        raise ValueError("start_row is outside the manifest")
    selected_rows = rows[start_row:end_row]
    cached_rows = []
    for index, row in enumerate(selected_rows, start=start_row):
        section = load_section_from_row(row)
        positions = {gene: position for position, gene in enumerate(section.gene_id)}
        missing = [gene for gene in genes if gene not in positions]
        if missing:
            raise ValueError(f"frozen panel is not fully observed in {section.section_id}")
        columns = [positions[gene] for gene in genes]
        counts = section.counts[:, columns].tocsr()
        stem = f"{index:04d}_{section.section_id}"
        matrix_path = cache_dir / f"{stem}.npz"
        metadata_path = cache_dir / f"{stem}.json"
        sparse.save_npz(matrix_path, counts)
        atomic_json(metadata_path, {
            "schema": "r04.panel_cache_metadata.v1",
            "section_id": section.section_id,
            "barcodes": list(section.barcode),
            "coords": section.coords.tolist(),
            "library_size": section.library_size.tolist() if section.library_size is not None else None,
            "gene_id": list(genes),
            "counts_shape": list(counts.shape),
            "source_matrix_locator": row["matrix_locator"],
        })
        cached = dict(row)
        cached.update({
            "matrix_locator": str(matrix_path.resolve()),
            "matrix_kind": "r04_panel_npz",
            "matrix_metadata_locator": str(metadata_path.resolve()),
            "panel_gene_hash": stable_json_hash(list(genes)),
        })
        cached_rows.append(cached)
    output = {
        "schema": "r04.role_manifest.v1",
        "status": "READY",
        "role": expected_role,
        "parent_manifest_hash": manifest.get("input_manifest_hash", ""),
        "input_manifest_hash": stable_json_hash(cached_rows),
        "panel_gene_hash": stable_json_hash(list(genes)),
        "n_rows": len(cached_rows),
        "n_patients": len({row.get("patient_id") for row in cached_rows}),
        "row_start": start_row,
        "row_end": start_row + len(cached_rows),
        "rows": cached_rows,
    }
    atomic_json(output_dir / f"{expected_role}_panel_manifest.json", output)
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument("--gene-list", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-role", choices=ROLE_NAMES, default="training")
    parser.add_argument("--start-row", type=int, default=0)
    parser.add_argument("--end-row", type=int, default=None)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if resource_status(args.output_dir) == "BLOCKED_STORAGE":
        atomic_json(args.output_dir / "panel_materialization.json", {"status": "BLOCKED_STORAGE"})
        return 2
    try:
        manifest = json.loads(args.manifest_json.read_text(encoding="utf-8"))
        genes = tuple(line.strip() for line in args.gene_list.read_text(encoding="utf-8").splitlines() if line.strip())
        output = materialize_panel(
            manifest,
            genes,
            args.output_dir,
            expected_role=args.expected_role,
            start_row=args.start_row,
            end_row=args.end_row,
        )
        atomic_json(args.output_dir / "panel_materialization.json", {
            "schema": "r04.panel_materialization.v1",
            "status": "READY",
            "role": args.expected_role,
            "manifest": f"{args.expected_role}_panel_manifest.json",
            "input_manifest_hash": output["input_manifest_hash"],
            "panel_gene_hash": output["panel_gene_hash"],
            "row_start": output["row_start"],
            "row_end": output["row_end"],
        })
        return 0
    except Exception as exc:
        atomic_json(args.output_dir / "panel_materialization.json", {
            "schema": "r04.panel_materialization.v1",
            "status": "BLOCKED_MATERIALIZATION",
            "error": type(exc).__name__,
            "message": str(exc),
        })
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
