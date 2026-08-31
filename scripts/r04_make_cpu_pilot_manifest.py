#!/usr/bin/env python3
"""Select the deterministic six-patient CPU pilot from frozen training rows."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from r04.runtime import atomic_json


def select_rows(manifest: dict, preflight: dict) -> list[dict[str, str]]:
    reports = {row["section_id"]: row for row in preflight["section_reports"]}
    training = [row for row in manifest["rows"] if row.get("primary_role") == "training"]
    by_patient: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in training:
        by_patient[row["patient_id"]].append(row)
    representatives = [max(rows, key=lambda row: reports[row["section_id"]]["n_spots"]) for rows in by_patient.values()]
    representatives.sort(key=lambda row: (reports[row["section_id"]]["n_spots"], row["patient_id"], row["section_id"]))
    if len(representatives) < 6:
        raise ValueError("at least six training patients are required for the CPU pilot")
    indices = [0, 1, len(representatives) // 2 - 1, len(representatives) // 2, len(representatives) - 2, len(representatives) - 1]
    selected = [representatives[index] for index in indices]
    if len({row["patient_id"] for row in selected}) != 6:
        raise ValueError("CPU pilot selection contains repeated patients")
    return selected


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument("--preflight-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest_json.read_text(encoding="utf-8"))
    preflight = json.loads(args.preflight_json.read_text(encoding="utf-8"))
    if manifest.get("status") != "READY" or preflight.get("status") != "READY":
        raise ValueError("CPU pilot requires READY input manifest and preflight gate")
    rows = select_rows(manifest, preflight)
    atomic_json(args.output, {
        "schema": "r04.locator_manifest.v1",
        "status": "READY",
        "input_manifest_hash": f"{manifest['input_manifest_hash']}::CPU_PILOT6",
        "pilot_policy": "two sections from each of the low, middle and high spot-count patient strata",
        "rows": rows,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
