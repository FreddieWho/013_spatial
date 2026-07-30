#!/usr/bin/env python3
"""Build the auditable R-01 registry from normalized staging metadata."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping


SOURCE_ASSET_FIELDS = [
    "asset_id",
    "source_namespace",
    "accession",
    "source_sample_id",
    "path",
    "asset_type",
    "format",
    "bytes",
    "checksum_status",
    "checksum",
    "processing_level",
    "parent_asset_id",
    "record_status",
]
PHYSICAL_UNIT_FIELDS = [
    "physical_unit_id",
    "source_namespace",
    "source_record_id",
    "study_id",
    "patient_id",
    "block_id",
    "physical_specimen_id",
    "identity_granularity",
    "block_equivalent_status",
    "block_equivalent_basis",
    "section_id",
    "z_position",
    "section_order",
    "section_thickness",
    "platform",
    "identity_status",
    "evidence_grade",
    "record_status",
]
IDENTITY_EVIDENCE_FIELDS = [
    "evidence_id",
    "entity_type",
    "entity_id",
    "field",
    "source_asset",
    "metadata_path",
    "metadata_key",
    "raw_value",
    "normalized_value",
    "evidence_grade",
    "interpretation",
    "conflict_flag",
]
DUPLICATE_GROUP_FIELDS = [
    "duplicate_group_id",
    "member_type",
    "member_id",
    "relation_type",
    "evidence_grade",
    "resolution",
    "leakage_group_id",
]
ROLE_FREEZE_FIELDS = [
    "logical_unit_id",
    "leakage_group_id",
    "primary_role",
    "special_capability",
    "freeze_version",
    "freeze_date",
    "rationale",
    "allowed_use",
    "forbidden_use",
    "record_status",
]


def _hash(raw: str, length: int = 16) -> str:
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:length]


def _namespaced(kind: str, namespace: str, raw: str) -> str:
    return f"{namespace}::{kind}::{_hash(raw)}" if raw else ""


def _composite_patient(raw: str) -> bool:
    return "\n" in raw or "&" in raw


def merge_external_rows(
    current: Iterable[Mapping[str, str]],
    external: Iterable[Mapping[str, str]],
    fields: list[str],
) -> list[dict[str, str]]:
    """Merge a schema-compatible external staging table fail-closed."""
    allowed = set(fields)
    merged: list[dict[str, str]] = []
    for row in [*current, *external]:
        unexpected = set(row) - allowed
        if unexpected:
            raise ValueError(
                "external staging has unexpected fields: "
                + ", ".join(sorted(unexpected))
            )
        merged.append({field: row.get(field, "") for field in fields})
    return merged


def physical_units_from_atlas(
    rows: Iterable[Mapping[str, str]],
    conflict_studies: Iterable[str] = (),
) -> list[dict[str, str]]:
    conflict_study_ids = set(conflict_studies)
    output = []
    for row in rows:
        record_id = row["record_id"].strip()
        study = row["study_namespace"].strip()
        raw_patient = row.get("patient_id_raw", "").strip()
        conflict = _composite_patient(raw_patient)
        if study in conflict_study_ids:
            patient_id = ""
            grade = "EC_conflict"
            status = "QUARANTINED_METADATA_CONFLICT"
            identity = "OFFICIAL_SOURCE_CONFLICT_PATIENT_BLOCK_UNKNOWN"
        elif conflict:
            patient_id = ""
            grade = "EC_conflict"
            status = "QUARANTINED_METADATA_CONFLICT"
            identity = "PATIENT_CONFLICT_BLOCK_UNKNOWN"
        elif raw_patient:
            patient_id = _namespaced("PATIENT", study, raw_patient)
            grade = "E3_explicit"
            status = "EXCLUDED_MISSING_BLOCK"
            identity = "PATIENT_EXPLICIT_BLOCK_UNKNOWN"
        else:
            patient_id = ""
            grade = "E0_unknown"
            status = "EXCLUDED_MISSING_PATIENT_OR_BLOCK"
            identity = "PATIENT_BLOCK_UNKNOWN"
        output.append(
            {
                "physical_unit_id": f"ATLAS::{record_id}",
                "source_namespace": "ATLAS_TABLE_S2",
                "source_record_id": record_id,
                "study_id": study,
                "patient_id": patient_id,
                "block_id": "",
                "section_id": "",
                "z_position": "",
                "section_order": "",
                "section_thickness": "",
                "platform": row.get("platform_raw", "").strip(),
                "identity_status": identity,
                "evidence_grade": grade,
                "record_status": status,
            }
        )
    return sorted(output, key=lambda row: row["physical_unit_id"])


def physical_units_from_hest(
    rows: Iterable[Mapping[str, str]],
) -> list[dict[str, str]]:
    output = []
    for row in rows:
        record_id = row["source_record_id"].strip()
        raw_study = row.get("raw_study_link", "").strip()
        study = f"HEST_STUDY::{_hash(raw_study)}" if raw_study else ""
        raw_patient = row.get("raw_patient", "").strip()
        conflict = row.get("identity_conflict") == "yes"
        if conflict:
            patient_id = ""
            status = "QUARANTINED_METADATA_CONFLICT"
            identity = "PATIENT_CONFLICT_BLOCK_UNKNOWN"
            grade = "EC_conflict"
        elif raw_patient:
            patient_id = _namespaced("PATIENT", study or "HEST_UNKNOWN", raw_patient)
            status = "EXCLUDED_MISSING_BLOCK"
            identity = "PATIENT_EXPLICIT_BLOCK_UNKNOWN"
            grade = row.get("evidence_grade", "E3_explicit")
        else:
            patient_id = ""
            status = "EXCLUDED_MISSING_PATIENT_OR_BLOCK"
            identity = "PATIENT_BLOCK_UNKNOWN"
            grade = "E0_unknown"

        raw_slide = row.get("raw_slide_id", "").strip()
        section_id = (
            _namespaced("SECTION", study or "HEST_UNKNOWN", raw_slide)
            if raw_slide and raw_slide.lower() not in {"n/a", "na", "-"}
            else ""
        )
        output.append(
            {
                "physical_unit_id": f"HEST::{record_id}",
                "source_namespace": "HEST",
                "source_record_id": record_id,
                "study_id": study,
                "patient_id": patient_id,
                "block_id": "",
                "section_id": section_id,
                "z_position": "",
                "section_order": "",
                "section_thickness": "",
                "platform": row.get("raw_st_technology", "").strip(),
                "identity_status": identity,
                "evidence_grade": grade,
                "record_status": status,
            }
        )
    return sorted(output, key=lambda row: row["physical_unit_id"])


def physical_units_from_htan(
    rows: Iterable[Mapping[str, str]],
) -> list[dict[str, str]]:
    study = "HTAN_VANDERBILT_CRC"
    output = []
    for row in rows:
        record_id = row["source_record_id"].strip()
        patient = row.get("raw_patient_id", "").strip()
        block = row.get("raw_block_id", "").strip()
        sample = row.get("raw_sample_key", "").strip()
        source_status = row.get("record_status", "")
        if source_status == "MATCHED" and patient and block and sample:
            status = "RESOLVED_INCLUDED_CANDIDATE"
            identity = "PATIENT_BLOCK_EXPLICIT_SECTION_UNRESOLVED"
        elif source_status == "METADATA_ONLY":
            status = "RESOLVED_ARCHIVE"
            identity = "METADATA_WITHOUT_LOCAL_ASSET"
        elif row.get("evidence_grade") == "EC_conflict":
            status = "EXCLUDED_METADATA_CONFLICT"
            identity = "IDENTITY_CONFLICT"
        else:
            status = "EXCLUDED_MISSING_PATIENT_OR_BLOCK"
            identity = "LOCAL_ASSET_IDENTITY_UNKNOWN"
        output.append(
            {
                "physical_unit_id": f"HTAN::{row['asset_name'].strip()}",
                "source_namespace": "HTAN_VANDERBILT_CRC",
                "source_record_id": record_id,
                "study_id": study,
                "patient_id": _namespaced("PATIENT", study, patient),
                "block_id": _namespaced("BLOCK", study, f"{patient}|{block}"),
                "section_id": "",
                "z_position": "",
                "section_order": "",
                "section_thickness": "",
                "platform": "Visium",
                "identity_status": identity,
                "evidence_grade": row.get("evidence_grade", "E0_unknown"),
                "record_status": status,
            }
        )
    return sorted(output, key=lambda row: row["physical_unit_id"])


def physical_units_from_tenx(
    rows: Iterable[Mapping[str, str]],
) -> list[dict[str, str]]:
    study = "TENX_V1_BREAST_CANCER_BLOCK_A"
    output = []
    for row in rows:
        dataset_id = row["dataset_id"].strip()
        patient = row.get("raw_patient_id", "").strip()
        block = row.get("raw_block_id", "").strip()
        section = row.get("raw_section_id", "").strip()
        matched = row.get("record_status") == "MATCHED"
        complete = bool(patient and block and section)
        if matched and complete:
            status = "RESOLVED_INCLUDED_CANDIDATE"
            identity = "PATIENT_BLOCK_SECTION_CORROBORATED"
        elif not patient:
            status = "EXCLUDED_MISSING_PATIENT_OR_BLOCK"
            identity = "PATIENT_UNKNOWN"
        elif not matched:
            status = "RESOLVED_ARCHIVE"
            identity = "METADATA_WITHOUT_LOCAL_ASSET"
        else:
            status = "EXCLUDED_MISSING_PATIENT_OR_BLOCK"
            identity = "IDENTITY_INCOMPLETE"
        output.append(
            {
                "physical_unit_id": f"TENX::{dataset_id}",
                "source_namespace": "TENX_GENOMICS",
                "source_record_id": row["source_record_id"].strip(),
                "study_id": study,
                "patient_id": _namespaced("PATIENT", study, patient),
                "block_id": _namespaced("BLOCK", study, f"{patient}|{block}"),
                "section_id": _namespaced(
                    "SECTION",
                    study,
                    f"{patient}|{block}|{section}",
                ),
                "z_position": "",
                "section_order": "",
                "section_thickness": "",
                "platform": "Visium",
                "identity_status": identity,
                "evidence_grade": row.get("evidence_grade", "E0_unknown"),
                "record_status": status,
            }
        )
    return sorted(output, key=lambda row: row["physical_unit_id"])


def build_source_assets(
    inventory_rows: Iterable[Mapping[str, str]],
    htan_rows: Iterable[Mapping[str, str]],
    root: Path,
    tenx_rows: Iterable[Mapping[str, str]] = (),
) -> list[dict[str, str]]:
    output = []
    for row in inventory_rows:
        source_id = row["source_id"]
        metadata_path = root / row["path_pattern"]
        if metadata_path.is_file():
            digest = hashlib.sha256()
            with metadata_path.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            metadata_bytes = str(metadata_path.stat().st_size)
            checksum_status = "sha256_verified"
            checksum = digest.hexdigest()
        else:
            metadata_bytes = row["bytes_total"]
            checksum_status = "AGGREGATE_SHA256"
            checksum = row["sha256"]
        output.append(
            {
                "asset_id": f"META::{source_id}",
                "source_namespace": source_id,
                "accession": "",
                "source_sample_id": "",
                "path": row["path_pattern"],
                "asset_type": "metadata_source",
                "format": row["kind"],
                "bytes": metadata_bytes,
                "checksum_status": checksum_status,
                "checksum": checksum,
                "processing_level": "bundled_metadata",
                "parent_asset_id": "",
                "record_status": row["status"],
            }
        )
    htan_dir = (
        root
        / "data/other_sources/htan/"
        "HTAN_Vanderbilt_CRC_Visium_OSF_hftq2"
    )
    for row in htan_rows:
        asset_name = row["asset_name"]
        path = htan_dir / asset_name
        output.append(
            {
                "asset_id": f"HTAN::{asset_name}",
                "source_namespace": "HTAN_VANDERBILT_CRC",
                "accession": "",
                "source_sample_id": row.get("raw_sample_key", ""),
                "path": str(path.relative_to(root)),
                "asset_type": "processed_expression_container",
                "format": "h5ad",
                "bytes": str(path.stat().st_size) if path.is_file() else "",
                "checksum_status": "NOT_COMPUTED_LARGE_ASSET",
                "checksum": "",
                "processing_level": "processed",
                "parent_asset_id": "META::upstream_crc_cohort_metadata",
                "record_status": row["record_status"],
            }
        )
    tenx_root = root / "data/other_sources/10x_genomics"
    for row in tenx_rows:
        dataset_id = row["dataset_id"]
        path = tenx_root / dataset_id
        output.append(
            {
                "asset_id": f"TENX::{dataset_id}",
                "source_namespace": "TENX_GENOMICS",
                "accession": "",
                "source_sample_id": row.get("atlas_sample_id", ""),
                "path": str(path.relative_to(root)),
                "asset_type": "dataset_directory",
                "format": "directory",
                "bytes": "",
                "checksum_status": "NOT_COMPUTED_DIRECTORY",
                "checksum": "",
                "processing_level": "processed",
                "parent_asset_id": "META::atlas_non_geo_manifest",
                "record_status": row["record_status"],
            }
        )
    return sorted(output, key=lambda row: row["asset_id"])


def build_identity_evidence(
    atlas_rows: Iterable[Mapping[str, str]],
    hest_rows: Iterable[Mapping[str, str]],
    htan_rows: Iterable[Mapping[str, str]],
    tenx_rows: Iterable[Mapping[str, str]] = (),
    conflict_atlas_records: Iterable[str] = (),
) -> list[dict[str, str]]:
    conflict_record_ids = set(conflict_atlas_records)
    output = []
    for row in atlas_rows:
        record_id = row["record_id"]
        field = row["field"]
        official_conflict = record_id in conflict_record_ids
        output.append(
            {
                "evidence_id": row["evidence_id"],
                "entity_type": "physical_unit",
                "entity_id": f"ATLAS::{record_id}",
                "field": field,
                "source_asset": row["source_asset"],
                "metadata_path": (
                    f"{row['source_sheet']}!{row['source_range']}"
                ),
                "metadata_key": field,
                "raw_value": row["raw_value"],
                "normalized_value": "",
                "evidence_grade": (
                    "EC_conflict"
                    if official_conflict
                    else row["evidence_grade"]
                ),
                "interpretation": (
                    "Quarantined because official GEO metadata conflicts "
                    "with bundled Atlas identity/provenance."
                    if official_conflict
                    else row["interpretation"]
                ),
                "conflict_flag": (
                    "yes"
                    if official_conflict
                    or row["evidence_grade"] == "EC_conflict"
                    else "no"
                ),
            }
        )

    hest_fields = [
        "raw_id",
        "raw_patient",
        "raw_subseries",
        "raw_study_link",
        "raw_download_page_link1",
        "raw_slide_id",
        "raw_z_step_size",
        "raw_st_technology",
    ]
    for row in hest_rows:
        record_id = row["source_record_id"]
        conflict = row.get("identity_conflict") == "yes"
        for field in hest_fields:
            raw = row.get(field, "")
            grade = (
                "EC_conflict"
                if conflict and field in {"raw_patient", "raw_subseries"}
                else ("E3_explicit" if raw else "E0_unknown")
            )
            output.append(
                {
                    "evidence_id": f"HEST::{record_id}::{field}",
                    "entity_type": "physical_unit",
                    "entity_id": f"HEST::{record_id}",
                    "field": field.removeprefix("raw_"),
                    "source_asset": "META::hest_metadata_json",
                    "metadata_path": row.get("metadata_path", ""),
                    "metadata_key": field.removeprefix("raw_"),
                    "raw_value": raw,
                    "normalized_value": "",
                    "evidence_grade": grade,
                    "interpretation": (
                        "Bundled HEST metadata; no block or serial inference."
                    ),
                    "conflict_flag": (
                        "yes" if grade == "EC_conflict" else "no"
                    ),
                }
            )

    for row in htan_rows:
        entity_id = f"HTAN::{row['asset_name']}"
        study = "HTAN_VANDERBILT_CRC"
        for field in ["raw_sample_key", "raw_patient_id", "raw_block_id"]:
            raw = row.get(field, "")
            if not raw:
                grade = "E0_unknown"
            elif field == "raw_sample_key":
                grade = "E1_weak"
            else:
                grade = row.get("evidence_grade", "E0_unknown")
            if field == "raw_patient_id" and raw:
                normalized = _namespaced("PATIENT", study, raw)
            elif field == "raw_block_id" and raw:
                normalized = _namespaced(
                    "BLOCK",
                    study,
                    f"{row.get('raw_patient_id', '')}|{raw}",
                )
            else:
                normalized = ""
            output.append(
                {
                    "evidence_id": f"{entity_id}::{field}",
                    "entity_type": "physical_unit",
                    "entity_id": entity_id,
                    "field": field.removeprefix("raw_"),
                    "source_asset": "META::upstream_crc_cohort_metadata",
                    "metadata_path": "repo/data_meta/ST_CRC_cohort_meta2.csv",
                    "metadata_key": field.removeprefix("raw_"),
                    "raw_value": raw,
                    "normalized_value": normalized,
                    "evidence_grade": grade,
                    "interpretation": "Bundled upstream HTAN CRC metadata.",
                    "conflict_flag": (
                        "yes" if grade == "EC_conflict" else "no"
                    ),
                }
            )
    for row in tenx_rows:
        entity_id = f"TENX::{row['dataset_id']}"
        study = "TENX_V1_BREAST_CANCER_BLOCK_A"
        for field in ["raw_patient_id", "raw_block_id", "raw_section_id"]:
            raw = row.get(field, "")
            if field == "raw_patient_id" and raw:
                normalized = _namespaced("PATIENT", study, raw)
            elif field == "raw_block_id" and raw:
                normalized = _namespaced(
                    "BLOCK",
                    study,
                    f"{row.get('raw_patient_id', '')}|{raw}",
                )
            elif field == "raw_section_id" and raw:
                normalized = _namespaced(
                    "SECTION",
                    study,
                    (
                        f"{row.get('raw_patient_id', '')}|"
                        f"{row.get('raw_block_id', '')}|{raw}"
                    ),
                )
            else:
                normalized = ""
            output.append(
                {
                    "evidence_id": f"{entity_id}::{field}",
                    "entity_type": "physical_unit",
                    "entity_id": entity_id,
                    "field": field.removeprefix("raw_"),
                    "source_asset": "META::atlas_non_geo_manifest",
                    "metadata_path": "data/other_sources/sources_manifest.tsv",
                    "metadata_key": field.removeprefix("raw_"),
                    "raw_value": raw,
                    "normalized_value": normalized,
                    "evidence_grade": (
                        row.get("evidence_grade", "E0_unknown")
                        if raw
                        else "E0_unknown"
                    ),
                    "interpretation": (
                        "Bundled 10x block/section description corroborated "
                        "with the Atlas Table S2 patient crosswalk."
                    ),
                    "conflict_flag": "no",
                }
            )
    return sorted(output, key=lambda row: row["evidence_id"])


def build_duplicate_groups(
    cross_rows: Iterable[Mapping[str, str]],
    htan_rows: Iterable[Mapping[str, str]],
    tenx_rows: Iterable[Mapping[str, str]] = (),
) -> list[dict[str, str]]:
    output: dict[tuple[str, str], dict[str, str]] = {}
    for row in cross_rows:
        if row.get("resolution") == "NO_EVIDENCE":
            continue
        left = row["left_record_id"]
        right = row["right_record_id"]
        group = f"EDGE::{_hash('|'.join(sorted([left, right])))}"
        leakage_group = row.get("leakage_group_id") or f"LEAKAGE::{group}"
        for member in [row["left_record_id"], row["right_record_id"]]:
            output[(group, member)] = {
                "duplicate_group_id": group,
                "member_type": "source_record",
                "member_id": member,
                "relation_type": row["relation_type"],
                "evidence_grade": row["evidence_grade"],
                "resolution": row["resolution"],
                "leakage_group_id": leakage_group,
            }

    matched = [
        row
        for row in htan_rows
        if row.get("record_status") == "MATCHED"
        and row.get("raw_patient_id")
        and row.get("raw_block_id")
        and row.get("raw_sample_key")
    ]
    patients: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    blocks: dict[tuple[str, str], list[Mapping[str, str]]] = defaultdict(list)
    sections: dict[tuple[str, str, str], list[Mapping[str, str]]] = defaultdict(list)
    for row in matched:
        patient = row["raw_patient_id"]
        block = row["raw_block_id"]
        sample = row["raw_sample_key"]
        patients[patient].append(row)
        blocks[(patient, block)].append(row)
        sections[(patient, block, sample)].append(row)

    patient_groups: dict[str, str] = {}
    for patient, members in patients.items():
        group = f"HTAN::PATIENT::{_hash(patient)}"
        patient_groups[patient] = group
        for row in members:
            output[(group, row["source_record_id"])] = {
                "duplicate_group_id": group,
                "member_type": "asset",
                "member_id": row["source_record_id"],
                "relation_type": "same_patient_other_block",
                "evidence_grade": "E3_explicit",
                "resolution": "CONFIRMED_METADATA_RELATION",
                "leakage_group_id": group,
            }
    for (patient, block), members in blocks.items():
        group = f"HTAN::BLOCK::{_hash(patient + '|' + block)}"
        for row in members:
            output[(group, row["source_record_id"])] = {
                "duplicate_group_id": group,
                "member_type": "asset",
                "member_id": row["source_record_id"],
                "relation_type": "same_block",
                "evidence_grade": "E3_explicit",
                "resolution": "CONFIRMED_METADATA_RELATION",
                "leakage_group_id": patient_groups[patient],
            }
    for (patient, block, sample), members in sections.items():
        if len(members) < 2:
            continue
        group = f"HTAN::SECTION::{_hash(patient + '|' + block + '|' + sample)}"
        for row in members:
            output[(group, row["source_record_id"])] = {
                "duplicate_group_id": group,
                "member_type": "asset",
                "member_id": row["source_record_id"],
                "relation_type": "same_sample_alias_within_block",
                "evidence_grade": "E1_weak",
                "resolution": "UNRESOLVED_CONSERVATIVE_GROUP",
                "leakage_group_id": patient_groups[patient],
            }

    tenx_by_block: dict[tuple[str, str], list[Mapping[str, str]]] = defaultdict(list)
    for row in tenx_rows:
        if (
            row.get("record_status") == "MATCHED"
            and row.get("raw_patient_id")
            and row.get("raw_block_id")
        ):
            tenx_by_block[
                (row["raw_patient_id"], row["raw_block_id"])
            ].append(row)
    for (patient, block), members in tenx_by_block.items():
        leakage = f"TENX::BLOCK::{_hash(patient + '|' + block)}"
        for row in members:
            member_id = row["source_record_id"]
            output[(leakage, member_id)] = {
                "duplicate_group_id": leakage,
                "member_type": "dataset_directory",
                "member_id": member_id,
                "relation_type": "same_block",
                "evidence_grade": row.get(
                    "evidence_grade",
                    "E2_corroborated",
                ),
                "resolution": "CONFIRMED_METADATA_RELATION",
                "leakage_group_id": leakage,
            }
    if matched:
        group = "SOURCE_LINEAGE::HTAN_VANDERBILT_CRC"
        output[(group, "HTAN_VANDERBILT_CRC")] = {
            "duplicate_group_id": group,
            "member_type": "study",
            "member_id": "HTAN_VANDERBILT_CRC",
            "relation_type": "canonical_source_lineage",
            "evidence_grade": "E3_explicit",
            "resolution": "canonical_study",
            "leakage_group_id": "LINEAGE::HTAN_VANDERBILT_CRC",
        }
    if tenx_by_block:
        group = "SOURCE_LINEAGE::TENX_V1_BREAST_CANCER_BLOCK_A"
        output[(group, "TENX_V1_BREAST_CANCER_BLOCK_A")] = {
            "duplicate_group_id": group,
            "member_type": "study",
            "member_id": "TENX_V1_BREAST_CANCER_BLOCK_A",
            "relation_type": "canonical_source_lineage",
            "evidence_grade": "E2_corroborated",
            "resolution": "canonical_study",
            "leakage_group_id": (
                "LINEAGE::TENX_V1_BREAST_CANCER_BLOCK_A"
            ),
        }
    return sorted(
        output.values(),
        key=lambda row: (row["duplicate_group_id"], row["member_id"]),
    )


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _read_optional_tsv(path: Path) -> list[dict[str, str]]:
    return _read_tsv(path) if path.is_file() else []


def _write_tsv(
    rows: Iterable[Mapping[str, str]],
    fields: list[str],
    output: Path,
) -> None:
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
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, output)
    finally:
        if temporary:
            Path(temporary).unlink(missing_ok=True)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--staging",
        type=Path,
        default=root / "infra/sample-registry/staging",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "infra/sample-registry",
    )
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    try:
        output_dir.relative_to(root)
    except ValueError as exc:
        raise ValueError("output directory escapes project root") from exc

    atlas_samples = _read_tsv(args.staging / "atlas_samples.tsv")
    atlas_evidence = _read_tsv(args.staging / "atlas_evidence.tsv")
    hest = _read_tsv(args.staging / "hest_identity_evidence.tsv")
    htan = _read_tsv(args.staging / "htan_crc_units.tsv")
    tenx = _read_tsv(args.staging / "tenx_explicit_physical_units.tsv")
    cross = _read_tsv(args.staging / "hest_tenx_duplicate_candidates.tsv")
    inventory = _read_tsv(
        root / "infra/sample-registry/source_metadata_inventory.tsv"
    )
    official_conflicts = _read_optional_tsv(
        args.staging / "official_metadata_conflicts.tsv"
    )
    conflict_studies = {
        row["study_id"] for row in official_conflicts if row.get("study_id")
    }
    conflict_atlas_records = {
        row["record_id"]
        for row in atlas_samples
        if row.get("study_namespace") in conflict_studies
    }

    physical = (
        physical_units_from_atlas(atlas_samples, conflict_studies)
        + physical_units_from_hest(hest)
        + physical_units_from_htan(htan)
        + physical_units_from_tenx(tenx)
    )
    physical = merge_external_rows(
        physical,
        _read_optional_tsv(args.staging / "external_geo_physical_units.tsv"),
        PHYSICAL_UNIT_FIELDS,
    )
    physical.sort(key=lambda row: row["physical_unit_id"])
    assets = merge_external_rows(
        build_source_assets(inventory, htan, root, tenx),
        _read_optional_tsv(args.staging / "external_geo_source_assets.tsv"),
        SOURCE_ASSET_FIELDS,
    )
    assets.sort(key=lambda row: row["asset_id"])
    evidence = merge_external_rows(
        build_identity_evidence(
            atlas_evidence,
            hest,
            htan,
            tenx,
            conflict_atlas_records,
        ),
        _read_optional_tsv(args.staging / "external_geo_identity_evidence.tsv"),
        IDENTITY_EVIDENCE_FIELDS,
    )
    evidence.sort(key=lambda row: row["evidence_id"])
    duplicates = merge_external_rows(
        build_duplicate_groups(cross, htan, tenx),
        _read_optional_tsv(args.staging / "external_geo_duplicate_groups.tsv"),
        DUPLICATE_GROUP_FIELDS,
    )
    duplicates.sort(
        key=lambda row: (row["duplicate_group_id"], row["member_id"])
    )

    _write_tsv(assets, SOURCE_ASSET_FIELDS, output_dir / "source_assets.tsv")
    _write_tsv(
        physical,
        PHYSICAL_UNIT_FIELDS,
        output_dir / "physical_units.tsv",
    )
    _write_tsv(
        evidence,
        IDENTITY_EVIDENCE_FIELDS,
        output_dir / "identity_evidence.tsv",
    )
    _write_tsv(
        duplicates,
        DUPLICATE_GROUP_FIELDS,
        output_dir / "duplicate_groups.tsv",
    )
    _write_tsv([], ROLE_FREEZE_FIELDS, output_dir / "role_freeze.tsv")
    print(
        f"registry assets={len(assets)} physical_units={len(physical)} "
        f"evidence={len(evidence)} duplicate_memberships={len(duplicates)} "
        "role_freeze=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
