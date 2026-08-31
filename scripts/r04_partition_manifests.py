#!/usr/bin/env python3
"""Create role-pure R-04 manifests from the frozen molecule-only manifest.

Training and validation are separate statistical operations.  This utility
keeps the parent manifest hash in every partition and fails closed when a
partition contains more than one role.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from r04.io_contract import stable_json_hash
from r04.runtime import atomic_json


ROLE_NAMES = (
    "training",
    "internal_validation",
    "external_validation",
    "serial_section_validation",
)


def partition_manifest(manifest: dict) -> dict[str, dict]:
    if manifest.get("status") != "READY":
        raise ValueError("role partition requires a READY parent manifest")
    rows = manifest.get("rows")
    if not isinstance(rows, list) or not rows:
        raise ValueError("parent manifest has no rows")
    parent_hash = str(manifest.get("input_manifest_hash", ""))
    if not parent_hash:
        raise ValueError("parent manifest hash is missing")
    result: dict[str, dict] = {}
    for role in ROLE_NAMES:
        selected = [row for row in rows if row.get("primary_role") == role]
        if not selected:
            continue
        observed = {str(row.get("primary_role", "")) for row in selected}
        if observed != {role}:
            raise ValueError(f"role partition is contaminated: {role}: {sorted(observed)}")
        partition_hash = stable_json_hash(selected)
        result[role] = {
            "schema": "r04.role_manifest.v1",
            "status": "READY",
            "role": role,
            "parent_manifest_hash": parent_hash,
            "partition_hash": partition_hash,
            "n_rows": len(selected),
            "n_patients": len({row.get("patient_id") for row in selected}),
            "rows": selected,
        }
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest_json.read_text(encoding="utf-8"))
    partitions = partition_manifest(manifest)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    index = {
        "schema": "r04.role_manifest_index.v1",
        "status": "READY",
        "parent_manifest_hash": manifest["input_manifest_hash"],
        "role_counts": dict(Counter(row.get("primary_role", "") for row in manifest["rows"])),
        "partitions": {},
    }
    for role, payload in partitions.items():
        filename = f"{role}_manifest.json"
        atomic_json(args.output_dir / filename, payload)
        index["partitions"][role] = {
            "path": filename,
            "partition_hash": payload["partition_hash"],
            "n_rows": payload["n_rows"],
            "n_patients": payload["n_patients"],
        }
    atomic_json(args.output_dir / "index.json", index)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
