#!/usr/bin/env python3
"""Extract explicit 10x block/section identity from bundled metadata."""

from __future__ import annotations

import argparse
import csv
import os
import re
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any


FIELDS = [
    "source_record_id",
    "dataset_id",
    "raw_patient_id",
    "raw_block_id",
    "raw_section_id",
    "atlas_sample_id",
    "local_present",
    "evidence_grade",
    "record_status",
    "source_url",
    "notes",
]
BLOCK_SECTION = re.compile(
    r"_Block_(?P<block>[A-Za-z0-9]+)_Section_(?P<section>\d+)"
)
PAPER_SAMPLE = re.compile(
    r"PAPER SAMPLE.*?:\s*(?:[A-Z]{2,}\s+)?(?P<sample>[A-Za-z]+\d+)\b",
    re.IGNORECASE,
)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def extract_explicit_tenx_units(
    source_manifest: Path,
    atlas_samples: Path,
    local_tenx_root: Path,
) -> list[dict[str, str]]:
    atlas_by_sample = {
        row["sample_id_raw"].strip(): row
        for row in _read_tsv(atlas_samples)
        if row.get("sample_id_raw", "").strip()
    }
    provisional: list[dict[str, Any]] = []
    for row in _read_tsv(source_manifest):
        if row.get("source", "").strip() != "10x Genomics":
            continue
        dataset_id = row.get("dataset_id", "").strip()
        match = BLOCK_SECTION.search(dataset_id)
        if not match:
            continue
        block = f"Block {match.group('block')}"
        section = f"Section {match.group('section')}"
        notes = row.get("notes", "").strip()
        paper_match = PAPER_SAMPLE.search(notes)
        sample_id = paper_match.group("sample") if paper_match else ""
        atlas_row = atlas_by_sample.get(sample_id)
        patient = atlas_row.get("patient_id_raw", "").strip() if atlas_row else ""
        provisional.append(
            {
                "source_record_id": f"TENX::{dataset_id}",
                "dataset_id": dataset_id,
                "raw_patient_id": patient,
                "raw_block_id": block,
                "raw_section_id": section,
                "atlas_sample_id": sample_id,
                "local_present": (
                    "yes" if (local_tenx_root / dataset_id).is_dir() else "no"
                ),
                "evidence_grade": "E2_corroborated" if patient else "E1_weak",
                "record_status": "",
                "source_url": row.get("url", "").strip(),
                "notes": notes,
                "_block_group": dataset_id.rsplit("_Section_", 1)[0],
            }
        )

    patients_by_block: dict[str, set[str]] = defaultdict(set)
    for row in provisional:
        if row["raw_patient_id"]:
            patients_by_block[row["_block_group"]].add(row["raw_patient_id"])
    for row in provisional:
        if not row["raw_patient_id"] and "same block" in row["notes"].lower():
            patients = patients_by_block[row["_block_group"]]
            if len(patients) == 1:
                row["raw_patient_id"] = next(iter(patients))
                row["evidence_grade"] = "E2_corroborated"
                row["notes"] += (
                    " | patient linked through bundled same-block statement "
                    "to the explicitly mapped section"
                )
        if not row["raw_patient_id"]:
            row["record_status"] = "PATIENT_UNKNOWN"
        elif row["local_present"] == "yes":
            row["record_status"] = "MATCHED"
        else:
            row["record_status"] = "METADATA_ONLY"
        del row["_block_group"]
    return sorted(provisional, key=lambda row: row["source_record_id"])


def write_rows(rows: list[dict[str, Any]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = ""
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            newline="",
            encoding="utf-8",
            dir=output.parent,
            prefix=f".{output.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = handle.name
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(sorted(rows, key=lambda row: row["source_record_id"]))
        os.replace(temporary, output)
    finally:
        if temporary:
            Path(temporary).unlink(missing_ok=True)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-manifest",
        type=Path,
        default=root / "data/other_sources/sources_manifest.tsv",
    )
    parser.add_argument(
        "--atlas-samples",
        type=Path,
        default=root / "infra/sample-registry/staging/atlas_samples.tsv",
    )
    parser.add_argument(
        "--local-tenx-root",
        type=Path,
        default=root / "data/other_sources/10x_genomics",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "infra/sample-registry/staging/"
            "tenx_explicit_physical_units.tsv"
        ),
    )
    args = parser.parse_args()
    output = args.output.resolve()
    try:
        output.relative_to(root)
    except ValueError as exc:
        raise ValueError("output escapes project root") from exc
    inputs = {
        args.source_manifest.resolve(),
        args.atlas_samples.resolve(),
    }
    if output in inputs:
        raise ValueError("output collides with input")

    rows = extract_explicit_tenx_units(
        args.source_manifest,
        args.atlas_samples,
        args.local_tenx_root,
    )
    write_rows(rows, output)
    print(
        f"explicit 10x block/section records={len(rows)} "
        f"matched={sum(row['record_status'] == 'MATCHED' for row in rows)}"
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
