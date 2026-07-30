#!/usr/bin/env python3
"""Summarize Atlas metadata into unfrozen R-01 logical-unit candidates."""

from __future__ import annotations

import argparse
import csv
import os
import re
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping


FIELDS = [
    "logical_unit_id",
    "cancer",
    "platforms",
    "preservation",
    "sample_count",
    "patient_value_count",
    "composite_patient_value_count",
    "potential_tls_positive_count",
    "potential_tls_instance_count",
    "references",
    "data_availability",
    "local_data_status",
    "identity_status",
    "gt_acceptance_status",
    "role_freeze_status",
    "notes",
]


def _joined(values: Iterable[str]) -> str:
    return "|".join(sorted({value.strip() for value in values if value.strip()}))


def _local_status(root: Path, availability: str) -> str:
    accessions = sorted(set(re.findall(r"\bGSE\d+\b", availability)))
    if accessions:
        if all((root / "data/GEO" / accession).is_dir() for accession in accessions):
            return "LOCAL_DIRECTORY_PRESENT"
        return "LOCAL_DIRECTORY_MISSING"

    lowered = availability.lower()
    if "humantumoratlas" in lowered or "vanderbilt_crc" in lowered:
        directory = (
            root
            / "data/other_sources/htan/"
            "HTAN_Vanderbilt_CRC_Visium_OSF_hftq2"
        )
        return (
            "LOCAL_DIRECTORY_PRESENT"
            if directory.is_dir()
            else "LOCAL_DIRECTORY_MISSING"
        )
    if "dryad.h70rxwdmj" in lowered:
        directory = root / "data/other_sources/dryad/h70rxwdmj"
        return (
            "LOCAL_DIRECTORY_PRESENT"
            if directory.is_dir()
            else "LOCAL_DIRECTORY_MISSING"
        )
    if "10xgenomics.com/datasets" in lowered:
        directory = root / "data/other_sources/10x_genomics"
        return (
            "LOCAL_BUNDLE_PRESENT"
            if directory.is_dir()
            else "LOCAL_DIRECTORY_MISSING"
        )
    if re.search(r"\b(EGA[SD]|HRA|PRJCA|OEP)\w+", availability, re.I):
        return "HOLD_OR_NOT_LOCAL"
    return "UNKNOWN_OR_NOT_MAPPED"


def _tls_count(raw: str) -> int:
    try:
        return int(float(raw.strip()))
    except (TypeError, ValueError):
        return 0


def summarize_atlas_units(
    sample_rows: Iterable[Mapping[str, str]],
    root: Path,
) -> list[dict[str, str]]:
    grouped: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in sample_rows:
        namespace = (row.get("study_namespace") or "").strip()
        if not namespace:
            raise ValueError("atlas sample lacks study_namespace")
        grouped[namespace].append(row)

    output: list[dict[str, str]] = []
    for namespace in sorted(grouped):
        rows = grouped[namespace]
        patient_values = {
            (row.get("patient_id_raw") or "").strip()
            for row in rows
            if (row.get("patient_id_raw") or "").strip()
        }
        composite = {
            value
            for value in patient_values
            if "\n" in value or "&" in value
        }
        missing_patients = sum(
            not (row.get("patient_id_raw") or "").strip() for row in rows
        )
        if composite:
            identity_status = "PATIENT_CONFLICT_BLOCK_UNKNOWN"
        elif missing_patients:
            identity_status = "PATIENT_INCOMPLETE_BLOCK_UNKNOWN"
        else:
            identity_status = "PATIENT_ONLY_BLOCK_UNKNOWN"

        availability = _joined(
            row.get("data_availability_raw", "") for row in rows
        )
        positive = sum(
            (row.get("tls_presence_raw") or "").strip().lower() == "yes"
            for row in rows
        )
        instances = sum(_tls_count(row.get("tls_count_raw") or "") for row in rows)
        output.append(
            {
                "logical_unit_id": namespace,
                "cancer": _joined(row.get("cancer_raw", "") for row in rows),
                "platforms": _joined(row.get("platform_raw", "") for row in rows),
                "preservation": _joined(
                    row.get("preservation_raw", "") for row in rows
                ),
                "sample_count": str(len(rows)),
                "patient_value_count": str(len(patient_values)),
                "composite_patient_value_count": str(len(composite)),
                "potential_tls_positive_count": str(positive),
                "potential_tls_instance_count": str(instances),
                "references": _joined(
                    row.get("reference_raw", "") for row in rows
                ),
                "data_availability": availability,
                "local_data_status": _local_status(root, availability),
                "identity_status": identity_status,
                "gt_acceptance_status": "NOT_ACCEPTED_R01",
                "role_freeze_status": "UNFROZEN_BLOCK_UNKNOWN",
                "notes": (
                    "TLS fields are recorded for provenance only and were not "
                    "used to accept GT or freeze a role."
                ),
            }
        )
    return output


def write_units(rows: list[dict[str, str]], output: Path) -> None:
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
            writer.writerows(
                sorted(rows, key=lambda row: row["logical_unit_id"])
            )
        os.replace(temporary, output)
    finally:
        if temporary:
            Path(temporary).unlink(missing_ok=True)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--samples",
        type=Path,
        default=root / "infra/sample-registry/staging/atlas_samples.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "infra/sample-registry/staging/atlas_logical_units.tsv",
    )
    args = parser.parse_args()
    output = args.output.resolve()
    try:
        output.relative_to(root)
    except ValueError as exc:
        raise ValueError("output escapes project root") from exc
    if output == args.samples.resolve():
        raise ValueError("output collides with input")

    with args.samples.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    units = summarize_atlas_units(rows, root)
    write_units(units, output)
    local = sum(
        unit["local_data_status"] in {
            "LOCAL_DIRECTORY_PRESENT",
            "LOCAL_BUNDLE_PRESENT",
        }
        for unit in units
    )
    print(f"atlas logical units={len(units)} local_or_bundled={local}")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
