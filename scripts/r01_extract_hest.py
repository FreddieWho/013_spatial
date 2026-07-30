#!/usr/bin/env python3
"""Extract auditable HEST R-01 identity evidence from bundled small metadata."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Iterable, Mapping


RECORD_COLUMNS = [
    "source_record_id",
    "metadata_path",
    "selected_manifest_present",
    "metadata_present",
    "record_status",
    "raw_selected_id",
    "raw_selected_organ",
    "raw_selected_disease_state",
    "raw_selected_oncotree_code",
    "raw_selected_st_technology",
    "raw_selected_preservation_method",
    "raw_id",
    "raw_patient",
    "raw_subseries",
    "raw_study_link",
    "raw_download_page_link1",
    "raw_platform",
    "raw_st_technology",
    "raw_dataset_title",
    "raw_species",
    "raw_organ",
    "raw_tissue",
    "raw_disease_state",
    "raw_disease_comment",
    "raw_oncotree_code",
    "raw_preservation_method",
    "raw_data_publication_date",
    "raw_license",
    "raw_treatment_comment",
    "raw_slide_id",
    "raw_region_name",
    "raw_sample_id",
    "raw_z_step_size",
    "patient_explicit_p_tokens",
    "subseries_explicit_p_tokens",
    "identity_conflict",
    "conflict_reasons",
    "evidence_grade",
    "interpretation",
]

RAW_METADATA_KEYS = {
    "raw_id": "id",
    "raw_patient": "patient",
    "raw_subseries": "subseries",
    "raw_study_link": "study_link",
    "raw_download_page_link1": "download_page_link1",
    "raw_platform": "platform",
    "raw_st_technology": "st_technology",
    "raw_dataset_title": "dataset_title",
    "raw_species": "species",
    "raw_organ": "organ",
    "raw_tissue": "tissue",
    "raw_disease_state": "disease_state",
    "raw_disease_comment": "disease_comment",
    "raw_oncotree_code": "oncotree_code",
    "raw_preservation_method": "preservation_method",
    "raw_data_publication_date": "data_publication_date",
    "raw_license": "license",
    "raw_treatment_comment": "treatment_comment",
    "raw_slide_id": "slide_id",
    "raw_region_name": "region_name",
    "raw_sample_id": "Sample ID",
    "raw_z_step_size": "z_step_size",
}

SELECTED_KEYS = {
    "raw_selected_id": "id",
    "raw_selected_organ": "organ",
    "raw_selected_disease_state": "disease_state",
    "raw_selected_oncotree_code": "oncotree_code",
    "raw_selected_st_technology": "st_technology",
    "raw_selected_preservation_method": "preservation_method",
}

COVERAGE_FIELDS = list(SELECTED_KEYS) + list(RAW_METADATA_KEYS)
CONFLICT_TYPES = [
    "patient_subseries_p_number_conflict",
    "metadata_filename_id_conflict",
]

_EXPLICIT_PATIENT_TOKEN = re.compile(
    r"(?i)(?<![A-Za-z0-9])(?:patient|p)[\s_-]*(\d+)(?!\d)"
)


def _raw(value: object) -> str:
    """Preserve a scalar raw value; represent missing values as unknown/empty."""
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return str(value)


def _read_selected(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "id" not in reader.fieldnames:
            raise ValueError(f"selected manifest lacks id column: {path}")
        records: dict[str, dict[str, str]] = {}
        for row in reader:
            record_id = (row.get("id") or "").strip()
            if not record_id:
                raise ValueError(f"selected manifest contains an empty id: {path}")
            if record_id in records:
                raise ValueError(f"duplicate selected id: {record_id}")
            records[record_id] = {key: _raw(value) for key, value in row.items()}
    return records


def _read_metadata(metadata_dir: Path) -> dict[str, tuple[Path, dict[str, object]]]:
    records: dict[str, tuple[Path, dict[str, object]]] = {}
    for path in sorted(metadata_dir.glob("*.json"), key=lambda item: item.name):
        record_id = path.stem
        if record_id in records:
            raise ValueError(f"duplicate metadata filename id: {record_id}")
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle)
        if not isinstance(payload, dict):
            raise ValueError(f"metadata root is not an object: {path}")
        records[record_id] = (path, payload)
    return records


def _explicit_p_tokens(value: str) -> tuple[int, ...]:
    return tuple(sorted({int(match) for match in _EXPLICIT_PATIENT_TOKEN.findall(value)}))


def _format_tokens(tokens: tuple[int, ...]) -> str:
    return "|".join(f"P{token}" for token in tokens)


def _evidence_grade(
    raw_patient: str,
    patient_tokens: tuple[int, ...],
    subseries_tokens: tuple[int, ...],
    conflict_reasons: list[str],
) -> tuple[str, str]:
    if conflict_reasons:
        return (
            "EC_conflict",
            "Bundled metadata fields conflict; no normalized identity is assigned.",
        )
    if raw_patient and patient_tokens and patient_tokens == subseries_tokens:
        return (
            "E2_corroborated",
            "Bundled patient and subseries fields contain the same explicit P-number.",
        )
    if raw_patient:
        return (
            "E3_explicit",
            "Bundled metadata explicitly supplies a patient value; it remains raw.",
        )
    if subseries_tokens:
        return (
            "E1_weak",
            "Only an explicit P-number inside subseries is available; no patient is assigned.",
        )
    return ("E0_unknown", "Bundled metadata does not establish patient identity.")


def extract_hest_records(
    selected_path: Path | str, metadata_dir: Path | str
) -> list[dict[str, str]]:
    """Return a deterministic union of selected-manifest and JSON metadata records."""
    selected_path = Path(selected_path)
    metadata_dir = Path(metadata_dir)
    selected = _read_selected(selected_path)
    metadata = _read_metadata(metadata_dir)
    rows: list[dict[str, str]] = []

    for record_id in sorted(set(selected) | set(metadata)):
        selected_row = selected.get(record_id)
        metadata_entry = metadata.get(record_id)
        metadata_path: Path | None = metadata_entry[0] if metadata_entry else None
        payload: Mapping[str, object] = metadata_entry[1] if metadata_entry else {}

        row = {column: "" for column in RECORD_COLUMNS}
        row["source_record_id"] = record_id
        row["metadata_path"] = (
            metadata_path.as_posix() if metadata_path is not None else ""
        )
        row["selected_manifest_present"] = "yes" if selected_row is not None else "no"
        row["metadata_present"] = "yes" if metadata_entry is not None else "no"
        if selected_row is not None and metadata_entry is not None:
            row["record_status"] = "MATCHED"
        elif selected_row is not None:
            row["record_status"] = "SELECTED_ONLY"
        else:
            row["record_status"] = "METADATA_ONLY"

        for output_field, input_field in SELECTED_KEYS.items():
            row[output_field] = (
                _raw(selected_row.get(input_field)) if selected_row is not None else ""
            )
        for output_field, input_field in RAW_METADATA_KEYS.items():
            row[output_field] = _raw(payload.get(input_field))

        patient_tokens = _explicit_p_tokens(row["raw_patient"])
        subseries_tokens = _explicit_p_tokens(row["raw_subseries"])
        row["patient_explicit_p_tokens"] = _format_tokens(patient_tokens)
        row["subseries_explicit_p_tokens"] = _format_tokens(subseries_tokens)

        conflict_reasons: list[str] = []
        if patient_tokens and subseries_tokens and patient_tokens != subseries_tokens:
            conflict_reasons.append("patient_subseries_p_number_conflict")
        if row["raw_id"] and row["raw_id"] != record_id:
            conflict_reasons.append("metadata_filename_id_conflict")
        row["identity_conflict"] = "yes" if conflict_reasons else "no"
        row["conflict_reasons"] = "|".join(conflict_reasons)
        row["evidence_grade"], row["interpretation"] = _evidence_grade(
            row["raw_patient"],
            patient_tokens,
            subseries_tokens,
            conflict_reasons,
        )
        rows.append(row)

    return rows


def build_field_coverage(rows: Iterable[Mapping[str, str]]) -> list[dict[str, str]]:
    materialized = list(rows)
    total = len(materialized)
    output: list[dict[str, str]] = []
    for field in COVERAGE_FIELDS:
        count = sum(1 for row in materialized if row.get(field, "") != "")
        output.append(
            {
                "field": field,
                "nonempty_count": str(count),
                "total_records": str(total),
                "coverage_fraction": f"{count / total:.6f}" if total else "0.000000",
            }
        )
    return output


def build_conflict_summary(
    rows: Iterable[Mapping[str, str]]
) -> list[dict[str, str]]:
    materialized = list(rows)
    output: list[dict[str, str]] = []
    for conflict_type in CONFLICT_TYPES:
        record_ids = sorted(
            row["source_record_id"]
            for row in materialized
            if conflict_type in row.get("conflict_reasons", "").split("|")
        )
        output.append(
            {
                "conflict_type": conflict_type,
                "record_count": str(len(record_ids)),
                "record_ids": "|".join(record_ids),
            }
        )
    any_conflicts = sorted(
        row["source_record_id"]
        for row in materialized
        if row.get("evidence_grade") == "EC_conflict"
    )
    output.append(
        {
            "conflict_type": "any_EC_conflict",
            "record_count": str(len(any_conflicts)),
            "record_ids": "|".join(any_conflicts),
        }
    )
    return output


def write_tsv(rows: Iterable[Mapping[str, str]], output: Path | str) -> None:
    materialized = list(rows)
    if not materialized:
        raise ValueError("refusing to write a headerless empty TSV")
    if "source_record_id" in materialized[0]:
        fieldnames = RECORD_COLUMNS
        materialized.sort(key=lambda row: row["source_record_id"])
    else:
        fieldnames = list(materialized[0])
        sort_field = fieldnames[0]
        materialized.sort(key=lambda row: row[sort_field])
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(materialized)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--selected",
        type=Path,
        default=Path("data/other_sources/hest1k/selected_samples.tsv"),
    )
    parser.add_argument(
        "--metadata-dir",
        type=Path,
        default=Path("data/other_sources/hest1k/metadata"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("infra/sample-registry/staging"),
    )
    args = parser.parse_args()

    rows = extract_hest_records(args.selected, args.metadata_dir)
    write_tsv(rows, args.output_dir / "hest_identity_evidence.tsv")
    write_tsv(
        build_field_coverage(rows),
        args.output_dir / "hest_field_coverage.tsv",
    )
    write_tsv(
        build_conflict_summary(rows),
        args.output_dir / "hest_conflict_summary.tsv",
    )
    conflicts = sum(row["evidence_grade"] == "EC_conflict" for row in rows)
    print(f"wrote {len(rows)} HEST records; EC_conflict={conflicts}")


if __name__ == "__main__":
    main()
