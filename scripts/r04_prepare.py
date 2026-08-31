#!/usr/bin/env python3
"""Build a molecule-only R-04 input manifest."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path

from r04.io_contract import read_manifest, sha256_file, stable_json_hash
from r04.runtime import atomic_json, resource_status


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    args = parser.parse_args()
    if resource_status(args.project_root) == "BLOCKED_STORAGE":
        atomic_json(args.output, {"status": "BLOCKED_STORAGE"})
        return 2
    rows = read_manifest(args.manifest)
    payload = []
    for row in rows:
        value = asdict(row)
        for key in ("matrix_locator", "coordinate_locator"):
            path = Path(value[key])
            if path.exists() and path.is_file():
                value[f"{key}_sha256"] = sha256_file(path)
            else:
                value[f"{key}_sha256"] = "UNAVAILABLE"
        payload.append(value)
    atomic_json(args.output, {
        "schema": "r04.input_manifest.v1",
        "status": "READY",
        "input_manifest_hash": stable_json_hash(payload),
        "rows": payload,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
