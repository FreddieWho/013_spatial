#!/usr/bin/env python3
"""Combine independently materialized, contiguous R-04 panel partitions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from r04.io_contract import stable_json_hash
from r04.runtime import atomic_json


def combine(paths: list[Path], output: Path, *, expected_rows: int) -> dict[str, object]:
    if not paths:
        raise ValueError("at least one panel manifest is required")
    manifests = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
    if any(value.get("status") != "READY" for value in manifests):
        raise ValueError("all panel manifests must be READY")
    roles = {value.get("role") for value in manifests}
    parents = {value.get("parent_manifest_hash") for value in manifests}
    panels = {value.get("panel_gene_hash") for value in manifests}
    if len(roles) != 1 or len(parents) != 1 or len(panels) != 1:
        raise ValueError("panel manifests disagree on role, parent, or gene panel")
    ordered = sorted(manifests, key=lambda value: int(value.get("row_start", -1)))
    rows = []
    expected_start = 0
    for value in ordered:
        start = int(value.get("row_start", -1))
        end = int(value.get("row_end", -1))
        if start != expected_start or end < start or end - start != len(value.get("rows", [])):
            raise ValueError("panel manifest row ranges are not contiguous")
        rows.extend(value["rows"])
        expected_start = end
    if expected_start != expected_rows or len(rows) != expected_rows:
        raise ValueError(f"expected {expected_rows} rows, got {len(rows)}")
    section_ids = [row.get("section_id") for row in rows]
    if len(section_ids) != len(set(section_ids)):
        raise ValueError("combined panel manifest contains duplicate sections")
    result = {
        "schema": "r04.role_manifest.v1",
        "status": "READY",
        "role": next(iter(roles)),
        "parent_manifest_hash": next(iter(parents)),
        "input_manifest_hash": stable_json_hash(rows),
        "panel_gene_hash": next(iter(panels)),
        "n_rows": len(rows),
        "n_patients": len({row.get("patient_id") for row in rows}),
        "rows": rows,
    }
    atomic_json(output, result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--part", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-rows", type=int, required=True)
    args = parser.parse_args()
    try:
        combine(args.part, args.output, expected_rows=args.expected_rows)
        return 0
    except Exception as exc:
        atomic_json(args.output.with_suffix(args.output.suffix + ".error.json"), {
            "schema": "r04.role_manifest.v1",
            "status": "BLOCKED_COMBINATION",
            "error": type(exc).__name__,
            "message": str(exc),
        })
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
