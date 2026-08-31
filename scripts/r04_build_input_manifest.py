#!/usr/bin/env python3
"""Join frozen R-02 splits with explicit matrix/coordinate locators.

The locator table is intentionally an input supplied by the operator. This
script never derives a sample identity or guesses a path from a filename.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from r04.io_contract import stable_json_hash
from r04.runtime import atomic_json


ROLE_ALLOWLIST = {"discovery", "training", "internal_validation", "external_validation", "serial_section_validation"}


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_input_rows(
    split_path: Path,
    locator_path: Path,
    project_root: Path,
    *,
    require_existing: bool = False,
) -> list[dict[str, str]]:
    splits = _read(split_path)
    locators = _read(locator_path)
    required_split = {"physical_unit_id", "logical_unit_id", "patient_id", "block_id", "leakage_group_id", "primary_role", "outer_fold", "block_level_eligible", "split_status"}
    required_locator = {"physical_unit_id", "matrix_locator", "coordinate_locator", "matrix_kind", "read_authorization"}
    if not splits or not required_split <= set(splits[0]):
        raise ValueError("outer split table lacks required R04 identity columns")
    if not locators or not required_locator <= set(locators[0]):
        raise ValueError("R04 locator table lacks explicit matrix/coordinate columns")
    locator_by_id = {row["physical_unit_id"]: row for row in locators}
    if len(locator_by_id) != len(locators):
        raise ValueError("R04 locator table contains duplicate physical_unit_id")
    rows: list[dict[str, str]] = []
    for split in splits:
        if split["split_status"] != "FROZEN_R02" or split["block_level_eligible"].lower() != "yes":
            continue
        if split["primary_role"] not in ROLE_ALLOWLIST:
            continue
        if not split["patient_id"]:
            raise ValueError(f"missing patient ID in eligible split row {split['physical_unit_id']}")
        locator = locator_by_id.get(split["physical_unit_id"])
        if locator is None:
            raise ValueError(f"no explicit R04 locator for eligible physical unit {split['physical_unit_id']}")
        matrix = Path(locator["matrix_locator"])
        if not matrix.is_absolute():
            matrix = project_root / matrix
        coordinate = locator["coordinate_locator"]
        if coordinate and not coordinate.startswith("IN_H5AD:"):
            coordinate_path = Path(coordinate)
            if not coordinate_path.is_absolute():
                coordinate_path = project_root / coordinate_path
            coordinate = str(coordinate_path)
        if require_existing:
            if not matrix.exists():
                raise FileNotFoundError(f"matrix locator does not exist: {split['physical_unit_id']}: {matrix}")
            if coordinate and not coordinate.startswith("IN_H5AD:") and not Path(coordinate).exists():
                raise FileNotFoundError(f"coordinate locator does not exist: {split['physical_unit_id']}: {coordinate}")
        rows.append({
            "section_id": split["physical_unit_id"],
            "patient_id": split["patient_id"],
            "block_id": split["block_id"],
            "lineage": split["logical_unit_id"],
            "matrix_locator": str(matrix),
            "coordinate_locator": coordinate,
            "matrix_kind": locator["matrix_kind"],
            "read_authorization": locator["read_authorization"],
            "count_key": locator.get("count_key", "X"),
            "outer_fold": split["outer_fold"],
            "leakage_group_id": split["leakage_group_id"],
            "primary_role": split["primary_role"],
        })
    if not rows:
        raise ValueError("no eligible R04 rows were produced")
    return sorted(rows, key=lambda row: (row["lineage"], row["patient_id"], row["section_id"]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outer-splits", type=Path, required=True)
    parser.add_argument("--locator-table", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = build_input_rows(args.outer_splits, args.locator_table, args.project_root.resolve(), require_existing=True)
    atomic_json(args.output, {
        "schema": "r04.locator_manifest.v1",
        "status": "READY",
        "input_manifest_hash": stable_json_hash(rows),
        "rows": rows,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
