#!/usr/bin/env python3
"""Build the explicit R-04 locator table from frozen registry assets.

HTAN locators come only from the pinned replay index.  Every non-HTAN asset
must be listed in the operator-maintained override table; no path is derived
from a sample filename or a target label.
"""

from __future__ import annotations

import argparse
import csv
import tempfile
from pathlib import Path


REQUIRED_SPLIT = {
    "physical_unit_id", "logical_unit_id", "patient_id", "block_id",
    "leakage_group_id", "primary_role", "outer_fold",
    "block_level_eligible", "split_status",
}
REQUIRED_LOCATOR = {
    "physical_unit_id", "matrix_locator", "coordinate_locator",
    "matrix_kind", "read_authorization",
}
ROLE_ALLOWLIST = {
    "discovery", "training", "internal_validation", "external_validation",
    "serial_section_validation",
}


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", newline="", dir=path.parent, delete=False) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    temporary.replace(path)


def build_locator_rows(
    *,
    outer_splits: Path,
    replay_index: Path,
    overrides: Path,
    project_root: Path,
    require_existing: bool = True,
) -> list[dict[str, str]]:
    splits = _read(outer_splits)
    replay = _read(replay_index)
    override_rows = _read(overrides)
    if not splits or not REQUIRED_SPLIT <= set(splits[0]):
        raise ValueError("outer split table lacks required R04 columns")
    if not replay or not {"physical_unit_id", "path"} <= set(replay[0]):
        raise ValueError("replay index lacks explicit physical_unit_id/path columns")
    if not override_rows or not REQUIRED_LOCATOR <= set(override_rows[0]):
        raise ValueError("locator override table lacks explicit locator columns")

    replay_by_id = {row["physical_unit_id"]: row for row in replay}
    override_by_id = {row["physical_unit_id"]: row for row in override_rows}
    if len(replay_by_id) != len(replay) or len(override_by_id) != len(override_rows):
        raise ValueError("locator source contains duplicate physical_unit_id")

    rows: list[dict[str, str]] = []
    for split in splits:
        if split["split_status"] != "FROZEN_R02" or split["block_level_eligible"].lower() != "yes":
            continue
        if split["primary_role"] not in ROLE_ALLOWLIST:
            continue
        physical_id = split["physical_unit_id"]
        locator = override_by_id.get(physical_id)
        if locator is None:
            replay_row = replay_by_id.get(physical_id)
            if replay_row is not None:
                locator = {
                    "physical_unit_id": physical_id,
                    "matrix_locator": replay_row["path"],
                    "coordinate_locator": "IN_H5AD:obsm/spatial",
                    "matrix_kind": "h5ad_counts",
                    "read_authorization": "HTAN_COUNTS_AND_COORDINATES",
                    "count_key": "X",
                }
        if locator is None:
            raise ValueError(f"no explicit locator for eligible physical unit {physical_id}")

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
                raise FileNotFoundError(f"matrix locator does not exist: {physical_id}: {matrix}")
            if coordinate and not coordinate.startswith("IN_H5AD:") and not Path(coordinate).exists():
                raise FileNotFoundError(f"coordinate locator does not exist: {physical_id}: {coordinate}")

        rows.append({
            "physical_unit_id": physical_id,
            "logical_unit_id": split["logical_unit_id"],
            "patient_id": split["patient_id"],
            "block_id": split["block_id"],
            "leakage_group_id": split["leakage_group_id"],
            "primary_role": split["primary_role"],
            "outer_fold": split["outer_fold"],
            "matrix_locator": str(matrix),
            "coordinate_locator": coordinate,
            "matrix_kind": locator["matrix_kind"],
            "read_authorization": locator["read_authorization"],
            "count_key": locator.get("count_key", "X"),
        })
    if not rows:
        raise ValueError("no eligible R04 locator rows were produced")
    return sorted(rows, key=lambda row: (row["logical_unit_id"], row["patient_id"], row["physical_unit_id"]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outer-splits", type=Path, required=True)
    parser.add_argument("--replay-index", type=Path, required=True)
    parser.add_argument("--overrides", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-missing", action="store_true")
    args = parser.parse_args()
    rows = build_locator_rows(
        outer_splits=args.outer_splits,
        replay_index=args.replay_index,
        overrides=args.overrides,
        project_root=args.project_root.resolve(),
        require_existing=not args.allow_missing,
    )
    _write(args.output, rows, list(rows[0]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
