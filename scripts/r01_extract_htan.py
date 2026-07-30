#!/usr/bin/env python3
"""Extract bundled HTAN CRC identity metadata without opening H5AD assets."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any


FIELDS = [
    "source_record_id",
    "study_id",
    "asset_name",
    "local_present",
    "raw_sample_key",
    "raw_patient_id",
    "raw_block_id",
    "raw_project",
    "duplicate_flag",
    "identity_conflict",
    "evidence_grade",
    "record_status",
    "notes",
]


def _metadata_asset_name(row: dict[str, str]) -> str:
    raw_path = (row.get("trimmed_adata") or "").strip()
    return Path(raw_path).name if raw_path else ""


def extract_htan_rows(
    metadata_csv: Path,
    asset_dir: Path,
) -> list[dict[str, Any]]:
    with metadata_csv.open(newline="", encoding="utf-8-sig") as handle:
        metadata_rows = list(csv.DictReader(handle))

    metadata_by_asset: dict[str, dict[str, str]] = {}
    for row in metadata_rows:
        asset_name = _metadata_asset_name(row)
        if not asset_name:
            raise ValueError("metadata row lacks trimmed_adata asset path")
        if asset_name in metadata_by_asset:
            raise ValueError(f"duplicate metadata asset: {asset_name}")
        metadata_by_asset[asset_name] = row

    local_assets = {path.name for path in asset_dir.glob("*.h5ad") if path.is_file()}
    all_assets = sorted(local_assets | set(metadata_by_asset))
    output: list[dict[str, Any]] = []

    for asset_name in all_assets:
        metadata = metadata_by_asset.get(asset_name)
        local_present = asset_name in local_assets
        if metadata is None:
            row = {
                "source_record_id": f"HTAN_VANDERBILT_CRC::{asset_name}",
                "study_id": "HTAN_VANDERBILT_CRC",
                "asset_name": asset_name,
                "local_present": "yes",
                "raw_sample_key": "",
                "raw_patient_id": "",
                "raw_block_id": "",
                "raw_project": "",
                "duplicate_flag": "",
                "identity_conflict": "no",
                "evidence_grade": "E0_unknown",
                "record_status": "LOCAL_ONLY",
                "notes": "No matching row in bundled ST_CRC_cohort_meta2.csv",
            }
        else:
            patient = (metadata.get("patient_name") or "").strip()
            block = (metadata.get("block_name") or "").strip()
            sample = (metadata.get("sample_key") or "").strip()
            required_complete = bool(patient and block and sample)
            status = "MATCHED" if local_present else "METADATA_ONLY"
            row = {
                "source_record_id": f"HTAN_VANDERBILT_CRC::{asset_name}",
                "study_id": "HTAN_VANDERBILT_CRC",
                "asset_name": asset_name,
                "local_present": "yes" if local_present else "no",
                "raw_sample_key": sample,
                "raw_patient_id": patient,
                "raw_block_id": block,
                "raw_project": (metadata.get("project") or "").strip(),
                "duplicate_flag": (metadata.get("duplicated") or "").strip(),
                "identity_conflict": "no" if required_complete else "yes",
                "evidence_grade": "E3_explicit" if required_complete else "EC_conflict",
                "record_status": status,
                "notes": (
                    ""
                    if required_complete
                    else "Bundled metadata lacks patient, block, or sample identity"
                ),
            }
        output.append(row)

    return output


def write_rows(rows: list[dict[str, Any]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(rows, key=lambda row: str(row["source_record_id"]))
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(ordered)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metadata",
        type=Path,
        default=root / "repo/data_meta/ST_CRC_cohort_meta2.csv",
    )
    parser.add_argument(
        "--assets",
        type=Path,
        default=(
            root
            / "data/other_sources/htan/"
            "HTAN_Vanderbilt_CRC_Visium_OSF_hftq2"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "infra/sample-registry/staging/htan_crc_units.tsv",
    )
    args = parser.parse_args()

    rows = extract_htan_rows(args.metadata, args.assets)
    write_rows(rows, args.output)
    counts: dict[str, int] = {}
    for row in rows:
        status = str(row["record_status"])
        counts[status] = counts.get(status, 0) + 1
    summary = " ".join(f"{key}={counts[key]}" for key in sorted(counts))
    print(f"HTAN CRC records={len(rows)} {summary}")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
