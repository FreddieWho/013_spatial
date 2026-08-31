#!/usr/bin/env python3
"""Build the R-02 metadata-only structure control plane.

The builder reads tabular metadata, workbook cells, compressed annotation headers,
checksums, and — via ``r02_heiser_gt`` — h5ad per-spot barcode/array-index
coordinate columns only (D-033). It never opens expression matrices or image
pixels.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping

from openpyxl import load_workbook

try:
    from scripts.r02_heiser_gt import VERIFIER_CLASS, build_heiser_records
    from scripts.r02_validation_gt import (
        VERIFIER_KIRC,
        VERIFIER_STCRC,
        VERIFIER_USZ,
        build_validation_records,
    )
except ModuleNotFoundError:
    # Direct ``python scripts/r02_build_registry.py`` execution places the
    # scripts directory, rather than the repository root, on sys.path.
    from r02_heiser_gt import VERIFIER_CLASS, build_heiser_records
    from r02_validation_gt import (
        VERIFIER_KIRC,
        VERIFIER_STCRC,
        VERIFIER_USZ,
        build_validation_records,
    )


ATLAS_PATH = Path("paper/tables/science.adz2742_tables_s1_to_s8.xlsx")
STOMICS_GLOB = (
    "data/other_sources/stomicsdb/STDS0000223/stomics/"
    "STSP*/GSM*_TLS_annotation.csv.gz"
)
EXPECTED_TLS_COUNTS = {
    "HTAN_VANDERBILT_CRC": 44,
    "GEO::GSE175540": 30,
    "GEO::GSE226997": 4,
    "GEO::GSE274103": 8,
    "GEO::GSE274557": 1,
}
HTAN_AVAILABILITY = (
    "https://data.humantumoratlas.org/publications/"
    "vanderbilt_crc_chen_2021?tab=10x-visium#"
)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, rows: Iterable[Mapping[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atlas_tls_candidates(root: Path) -> list[dict[str, str]]:
    workbook = load_workbook(root / ATLAS_PATH, read_only=True, data_only=True)
    summary_sheet = workbook["Table S2"]
    availability_to_unit = {
        HTAN_AVAILABILITY: "HTAN_VANDERBILT_CRC",
        "GSE175540": "GEO::GSE175540",
        "GSE226997": "GEO::GSE226997",
        "GSE274103": "GEO::GSE274103",
        "GSE274557": "GEO::GSE274557",
    }
    sample_scope: dict[str, dict[str, object]] = {}
    summary_counts: dict[str, int] = {key: 0 for key in EXPECTED_TLS_COUNTS}
    current_availability = ""
    for excel_row, values in enumerate(
        summary_sheet.iter_rows(min_row=4, values_only=True), start=4
    ):
        patient = values[2]
        sample = values[3]
        presence = values[8]
        raw_count = values[9]
        if values[11] is not None:
            current_availability = str(values[11]).strip()
        availability = current_availability
        logical_unit_id = availability_to_unit.get(str(availability).strip())
        if logical_unit_id is None:
            continue
        try:
            count = int(raw_count or 0)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"non-integer TLS count at Table S2 row {excel_row}") from exc
        if count < 0:
            raise ValueError(f"negative TLS count at Table S2 row {excel_row}")
        summary_counts[logical_unit_id] += count
        sample_key = str(sample or "")
        if not sample_key:
            raise ValueError(f"missing sample ID at Table S2 row {excel_row}")
        if sample_key in sample_scope:
            previous = sample_scope[sample_key]
            raise ValueError(
                "ambiguous scoped Table S2 sample ID "
                f"{sample_key}: rows {previous['summary_row']} and {excel_row}, "
                f"units {previous['logical_unit_id']} and {logical_unit_id}"
            )
        sample_scope[sample_key] = {
            "logical_unit_id": logical_unit_id,
            "patient": str(patient or ""),
            "presence": str(presence or ""),
            "count": count,
            "summary_row": excel_row,
        }
    if summary_counts != EXPECTED_TLS_COUNTS:
        raise ValueError(
            f"scoped TLS summary mismatch: expected={EXPECTED_TLS_COUNTS}, "
            f"observed={summary_counts}"
        )

    detail_sheet = workbook["Table S4"]
    rows: list[dict[str, str]] = []
    detail_by_sample: dict[str, int] = {key: 0 for key in sample_scope}
    for excel_row, values in enumerate(
        detail_sheet.iter_rows(min_row=4, values_only=True), start=4
    ):
        if values[0] is None:
            continue
        tls_id = str(values[0]).strip()
        sample_key = tls_id.rsplit("_", 1)[0]
        scope = sample_scope.get(sample_key)
        if scope is None:
            continue
        detail_by_sample[sample_key] += 1
        logical_unit_id = str(scope["logical_unit_id"])
        rows.append(
            {
                "candidate_id": f"TLS_SOURCE::{logical_unit_id}::{tls_id}",
                "source_tls_id": tls_id,
                "structure_id": "TLS",
                "logical_unit_id": logical_unit_id,
                "source_record_id": f"Table S4!A{excel_row}:D{excel_row}",
                "source_patient_value": str(scope["patient"]),
                "source_sample_value": sample_key,
                "reported_presence": str(scope["presence"]),
                "reported_count_in_source_row": str(scope["count"]),
                "maturation_state_raw": str(values[1] or ""),
                "spatial_location_raw": str(values[2] or ""),
                "size_um2_raw": str(values[3] or ""),
                "gt_source_id": "ATLAS_TABLE_S4_TLS_IDS",
                "candidate_status": "SOURCE_REPORTED_ID_NONCONFIRMATORY",
                "exclusion_reason": (
                    "source-reported TLS ID has no replayable geometry locator, unique "
                    "R-01 physical link, or closed assay-overlap boundary"
                ),
            }
        )
    workbook.close()
    per_sample_mismatch = {
        sample: (int(scope["count"]), detail_by_sample[sample])
        for sample, scope in sample_scope.items()
        if int(scope["count"]) != detail_by_sample[sample]
    }
    if per_sample_mismatch:
        raise ValueError(f"Table S2/Table S4 scoped TLS mismatch: {per_sample_mismatch}")
    observed = Counter(row["logical_unit_id"] for row in rows)
    if dict(observed) != EXPECTED_TLS_COUNTS:
        raise ValueError(
            f"scoped Table S4 TLS ID mismatch: expected={EXPECTED_TLS_COUNTS}, "
            f"observed={dict(observed)}"
        )
    if len(rows) != 87:
        raise ValueError(f"expected 87 source-reported TLS IDs, observed {len(rows)}")
    return rows


def _stomics_audits(root: Path) -> list[dict[str, str]]:
    paths = sorted(root.glob(STOMICS_GLOB))
    if len(paths) != 19:
        raise ValueError(f"expected 19 STOmics TLS annotation files, observed {len(paths)}")
    rows: list[dict[str, str]] = []
    for index, path in enumerate(paths, start=1):
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            header = next(csv.reader(handle))
        relative = path.relative_to(root).as_posix()
        rows.append(
            {
                "gt_source_id": f"STOMICS_TLS_ANNOTATION::{index:02d}",
                "source_class": "STOMICS_TLS_ANNOTATION_CANDIDATE",
                "structure_id": "TLS",
                "path": relative,
                "sha256": _sha256(path),
                "checksum_status": "VERIFIED",
                "schema_locator": "gzip-csv header:" + ",".join(header),
                "physical_unit_id": "",
                "physical_link_status": "OUTSIDE_R01_FROZEN_CORE",
                "geometry_locator": "",
                "provenance_status": "UNRESOLVED",
                "provenance_locator": "",
                "overlap_status": "UNKNOWN",
                "overlap_locator": "",
                "same_assay_status": "UNKNOWN",
                "audit_status": "NOT_AUDITABLE_GT",
                "notes": (
                    "Header/schema inventory only; file is outside the R-01 frozen core and "
                    "its patient/block mapping, annotation provenance, geometry and overlap "
                    "boundary are not closed. Not included in the 57 scoped summaries."
                ),
            }
        )
    return rows


def _build_outer_splits(root: Path) -> list[dict[str, str]]:
    physical = _read_tsv(root / "infra/sample-registry/physical_units.tsv")
    roles = {
        row["logical_unit_id"]: row
        for row in _read_tsv(root / "infra/sample-registry/role_freeze.tsv")
        if row["record_status"] == "FROZEN"
    }
    rows: list[dict[str, str]] = []
    for row in physical:
        if row["record_status"] != "RESOLVED_INCLUDED_CANDIDATE":
            continue
        logical_unit = row["study_id"]
        role = roles.get(logical_unit)
        if role is None:
            raise ValueError(f"eligible physical unit lacks frozen role: {row['physical_unit_id']}")
        if row["block_id"]:
            granularity = row["identity_granularity"] or "explicit_patient_block"
            block_group = f"block::{row['block_id']}"
            block_level = "yes"
        else:
            if row["identity_granularity"] != "patient_linked_physical_specimen":
                raise ValueError(f"eligible row lacks block and specimen envelope: {row['physical_unit_id']}")
            granularity = row["identity_granularity"]
            block_group = ""
            block_level = "no"
        envelope = f"patient_envelope::{logical_unit}::{row['patient_id']}"
        fold = f"{role['primary_role']}::{hashlib.sha256(envelope.encode()).hexdigest()[:12]}"
        rows.append(
            {
                "physical_unit_id": row["physical_unit_id"],
                "logical_unit_id": logical_unit,
                "patient_id": row["patient_id"],
                "block_id": row["block_id"],
                "physical_specimen_id": row["physical_specimen_id"],
                "identity_granularity": granularity,
                "identity_envelope_id": envelope,
                "block_group_id": block_group,
                "leakage_group_id": role["leakage_group_id"],
                "primary_role": role["primary_role"],
                "outer_fold": fold,
                "block_level_eligible": block_level,
                "split_status": "FROZEN_R02",
            }
        )
    rows.sort(key=lambda row: row["physical_unit_id"])
    if len(rows) != 167:
        raise ValueError(f"expected 167 R-01 eligible physical units, observed {len(rows)}")
    return rows


def _cross_section_links(splits: list[dict[str, str]]) -> list[dict[str, str]]:
    tenx = [
        row
        for row in splits
        if row["logical_unit_id"] == "TENX_V1_BREAST_CANCER_BLOCK_A"
    ]
    if len(tenx) != 2:
        raise ValueError(f"expected two explicit TENX sections, observed {len(tenx)}")
    left, right = sorted(tenx, key=lambda row: row["physical_unit_id"])
    if not left["block_id"] or left["block_id"] != right["block_id"]:
        raise ValueError("TENX serial-section candidates do not share an explicit block")
    return [
        {
            "link_id": "TENX_BLOCK_A_SECTION_PAIR",
            "left_physical_unit_id": left["physical_unit_id"],
            "right_physical_unit_id": right["physical_unit_id"],
            "patient_id": left["patient_id"],
            "block_id": left["block_id"],
            "relation_type": "SAME_EXPLICIT_BLOCK_SECTION_PAIR",
            "evidence_status": "IDENTITY_ONLY_NOT_STRUCTURE_CORRESPONDENCE",
            "outer_fold": left["outer_fold"],
            "cross_fold_allowed": "no",
        }
    ]


def build_registry(root: Path, output: Path) -> None:
    root = root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    candidates = _atlas_tls_candidates(root)
    splits = _build_outer_splits(root)
    atlas_path = root / ATLAS_PATH
    physical_rows = _read_tsv(root / "infra/sample-registry/physical_units.tsv")
    heiser = build_heiser_records(root, physical_rows)
    validation = build_validation_records(root, physical_rows)

    ontology = [
        {
            "structure_id": "TLS",
            "label": "tertiary lymphoid structure",
            "geometry_type": "bounded region or focus; source-specific",
            "minimum_gt_requirement": "auditable instance geometry with independent provenance and overlap boundary",
            "direct_channels_to_exclude": "TLS annotation; defining TLS marker panel; source GT raster/vector",
            "boundary_uncertainty_policy": "source-specific uncertainty required",
            "claim_status": "FROZEN_CLAIM_BEARING",
        },
        {
            "structure_id": "BLOOD_VESSEL",
            "label": "blood vessel",
            "geometry_type": "tubular/linear region",
            "minimum_gt_requirement": "instance-level lumen/wall geometry and independent modality provenance",
            "direct_channels_to_exclude": "vessel annotation; defining endothelial marker panel; source GT raster/vector",
            "boundary_uncertainty_policy": "source-specific uncertainty required",
            "claim_status": "NOT_FROZEN_NO_SCOPED_GT",
        },
        {
            "structure_id": "NECROSIS",
            "label": "necrosis",
            "geometry_type": "bounded region",
            "minimum_gt_requirement": "instance-level pathologist or orthogonal geometry with provenance",
            "direct_channels_to_exclude": "necrosis annotation; source GT raster/vector",
            "boundary_uncertainty_policy": "source-specific uncertainty required",
            "claim_status": "NOT_FROZEN_NO_SCOPED_GT",
        },
        {
            "structure_id": "TUMOR_STROMA_BOUNDARY",
            "label": "tumor-stroma boundary",
            "geometry_type": "interface/curve",
            "minimum_gt_requirement": "auditable tumor and stroma geometry with boundary derivation",
            "direct_channels_to_exclude": "compartment annotation; source boundary raster/vector",
            "boundary_uncertainty_policy": "source-specific uncertainty required",
            "claim_status": "FROZEN_CLAIM_BEARING",
        },
    ]
    _write_tsv(
        output / "structure_ontology.tsv",
        ontology,
        [
            "structure_id",
            "label",
            "geometry_type",
            "minimum_gt_requirement",
            "direct_channels_to_exclude",
            "boundary_uncertainty_policy",
            "claim_status",
        ],
    )

    audits = [
        {
            "gt_source_id": "ATLAS_TABLE_S2_TLS_SUMMARY",
            "source_class": "ATLAS_TLS_COUNT_SUMMARY",
            "structure_id": "TLS",
            "path": ATLAS_PATH.as_posix(),
            "sha256": _sha256(atlas_path),
            "checksum_status": "VERIFIED",
            "schema_locator": "Table S2!I:J (TLS presence/count)",
            "physical_unit_id": "",
            "physical_link_status": "UNRESOLVED_ROW_TO_R01_PHYSICAL_UNIT",
            "geometry_locator": "",
            "provenance_status": "SUMMARY_PROVENANCE_ONLY",
            "provenance_locator": "Table S2!K:L (publication/data source only)",
            "overlap_status": "UNKNOWN",
            "overlap_locator": "",
            "same_assay_status": "UNKNOWN",
            "audit_status": "NOT_AUDITABLE_GT",
            "notes": (
                "Presence/count summaries enumerate candidates but provide no instance geometry, "
                "replayable annotation derivation, or closed physical/overlap mapping."
            ),
        },
        {
            "gt_source_id": "ATLAS_TABLE_S4_TLS_IDS",
            "source_class": "ATLAS_TLS_ID_ATTRIBUTES",
            "structure_id": "TLS",
            "path": ATLAS_PATH.as_posix(),
            "sha256": _sha256(atlas_path),
            "checksum_status": "VERIFIED",
            "schema_locator": "Table S4!A:D (TLS ID, maturation, location, size)",
            "physical_unit_id": "",
            "physical_link_status": "UNRESOLVED_ID_TO_R01_PHYSICAL_UNIT",
            "geometry_locator": "",
            "provenance_status": "SOURCE_REPORTED_ID_ATTRIBUTES",
            "provenance_locator": "Table S4!A:D",
            "overlap_status": "UNKNOWN",
            "overlap_locator": "",
            "same_assay_status": "UNKNOWN",
            "audit_status": "NOT_AUDITABLE_GT",
            "notes": (
                "Verified source-reported IDs replace synthetic count expansion, but "
                "Table S4 supplies attributes rather than replayable structure geometry."
            ),
        },
        {
            "gt_source_id": "UPSTREAM_CRC_TLS_ID_SUMMARY",
            "source_class": "UPSTREAM_CRC_TLS_ID_SUMMARY",
            "structure_id": "TLS",
            "path": "repo/data/ST_CRC_maturation_location.csv",
            "sha256": _sha256(root / "repo/data/ST_CRC_maturation_location.csv"),
            "checksum_status": "VERIFIED",
            "schema_locator": "CSV:ST_ID,sample_key,TLS_ID,Cluster,Location",
            "physical_unit_id": "",
            "physical_link_status": "SAMPLE_KEY_LINK_ONLY_NOT_UNIQUE_PHYSICAL",
            "geometry_locator": "",
            "provenance_status": "CORROBORATING_DERIVED_SUMMARY",
            "provenance_locator": "repo/data/ST_CRC_maturation_location.csv#TLS_ID",
            "overlap_status": "UNKNOWN",
            "overlap_locator": "",
            "same_assay_status": "UNKNOWN",
            "audit_status": "NOT_AUDITABLE_GT",
            "notes": (
                "Corroborates all 44 scoped HTAN Table S4 TLS IDs and sample keys; "
                "contains no replayable geometry and is not independent GT."
            ),
        },
        *_stomics_audits(root),
        *heiser["gt_sources"],
        *validation["gt_sources"],
    ]
    audit_fields = [
        "gt_source_id",
        "source_class",
        "structure_id",
        "path",
        "sha256",
        "checksum_status",
        "schema_locator",
        "physical_unit_id",
        "physical_link_status",
        "geometry_locator",
        "provenance_status",
        "provenance_locator",
        "overlap_status",
        "overlap_locator",
        "same_assay_status",
        "audit_status",
        "notes",
    ]
    _write_tsv(output / "gt_source_audit.tsv", audits, audit_fields)

    candidate_fields = [
        "candidate_id",
        "source_tls_id",
        "structure_id",
        "logical_unit_id",
        "source_record_id",
        "source_patient_value",
        "source_sample_value",
        "reported_presence",
        "reported_count_in_source_row",
        "maturation_state_raw",
        "spatial_location_raw",
        "size_um2_raw",
        "gt_source_id",
        "candidate_status",
        "exclusion_reason",
    ]
    _write_tsv(output / "tls_candidate_summary.tsv", candidates, candidate_fields)
    instances = [
        {
            "instance_id": candidate["candidate_id"],
            "structure_id": "TLS",
            "logical_unit_id": candidate["logical_unit_id"],
            "physical_unit_id": "",
            "patient_id": "",
            "block_id": "",
            "physical_specimen_id": "",
            "gt_source_id": candidate["gt_source_id"],
            "gt_definition_modality": "reported count summary",
            "gt_geometry_locator": "",
            "boundary_uncertainty": "UNKNOWN",
            "confirmation_status": "NONCONFIRMATORY_SOURCE_REPORTED_ID",
            "allowed_use": "inventory and GT acquisition planning only",
            "forbidden_use": "training, validation, localization, distance fields, or claims",
        }
        for candidate in candidates
    ]
    instances.extend(heiser["instances"])
    instances.extend(validation["instances"])
    _write_tsv(
        output / "structure_instances.tsv",
        instances,
        [
            "instance_id",
            "structure_id",
            "logical_unit_id",
            "physical_unit_id",
            "patient_id",
            "block_id",
            "physical_specimen_id",
            "gt_source_id",
            "gt_definition_modality",
            "gt_geometry_locator",
            "boundary_uncertainty",
            "confirmation_status",
            "allowed_use",
            "forbidden_use",
        ],
    )
    replay_rows = heiser["replay_index"] + validation["replay_index"]
    replay_rows.sort(key=lambda row: row["physical_unit_id"])
    _write_tsv(
        output / "h5ad_replay_index.tsv",
        replay_rows,
        [
            "physical_unit_id",
            "path",
            "spot_count",
            "index_fingerprint_sha256",
            "fingerprint_status",
        ],
    )

    input_policy = [
        {
            "channel_class": "data-carried metadata and schema",
            "policy": "ALLOWED_R02_AUDIT",
            "scope": "identity, provenance, annotation schema/header, checksums",
            "reason": "Required to audit GT-input separation without scientific computation.",
        },
        {
            "channel_class": "expression matrix or derived expression",
            "policy": "FORBIDDEN_R02",
            "scope": "X, counts, normalized expression, embeddings, programs",
            "reason": "R-02 is metadata/schema-only and must not inspect model inputs.",
        },
        {
            "channel_class": "image pixels or derived morphology",
            "policy": "FORBIDDEN_R02",
            "scope": "WSI, H&E pixels, masks, embeddings",
            "reason": "Only file/schema provenance may be inventoried.",
        },
        {
            "channel_class": "GT-defining annotation",
            "policy": "FORBIDDEN_AS_MODEL_INPUT",
            "scope": "TLS/structure labels, raster/vector geometry, defining marker panels",
            "reason": "Prevents circular GT and direct-label leakage.",
        },
        {
            "channel_class": "H&E image and image-derived features",
            "policy": "FORBIDDEN_AS_MODEL_INPUT",
            "scope": (
                "WSI pixels, stains and morphology embeddings on tasks whose "
                "ground truth is H&E pathology annotation "
                "(anchor he-annotation-gt-input-exclusion)"
            ),
            "reason": (
                "The auditable Heiser GT is defined by manual annotation on the "
                "H&E image; H&E-derived model inputs on those tasks would be "
                "same-modality circular (KNOWN_OVERLAP)."
            ),
        },
        {
            "channel_class": "outcome or response metadata",
            "policy": "FORBIDDEN_FOR_SELECTION_SPLIT_MODEL",
            "scope": "response tokens in sample titles and any clinical outcome",
            "reason": "Outcome-blind role and fold freeze.",
        },
        {
            "channel_class": "identity, site, treatment, stage, sample and file metadata",
            "policy": "FORBIDDEN_AS_PREDICTIVE_INPUT",
            "scope": "patient/sample IDs, filenames, section labels, tissue site, treatment and TNM/stage",
            "reason": "Prevents identity, acquisition and clinical shortcut leakage.",
        },
    ]
    _write_tsv(
        output / "input_policy.tsv",
        input_policy,
        ["channel_class", "policy", "scope", "reason"],
    )

    _write_tsv(
        output / "outer_splits.tsv",
        splits,
        [
            "physical_unit_id",
            "logical_unit_id",
            "patient_id",
            "block_id",
            "physical_specimen_id",
            "identity_granularity",
            "identity_envelope_id",
            "block_group_id",
            "leakage_group_id",
            "primary_role",
            "outer_fold",
            "block_level_eligible",
            "split_status",
        ],
    )
    _write_tsv(
        output / "cross_section_links.tsv",
        _cross_section_links(splits),
        [
            "link_id",
            "left_physical_unit_id",
            "right_physical_unit_id",
            "patient_id",
            "block_id",
            "relation_type",
            "evidence_status",
            "outer_fold",
            "cross_fold_allowed",
        ],
    )
    leakage = [
        {
            "audit_id": "GSE211956_RESPONSE_TITLE",
            "logical_unit_id": "GEO::GSE211956",
            "risk_type": "OUTCOME_IN_SAMPLE_TITLE",
            "source_path": "infra/bioinf-data-index/raw/geo/GSE211956_gsm_quick.soft",
            "metadata_locator": "!Sample_title",
            "observed_fact": "sample titles include poor/good/partial response tokens",
            "control": "forbidden for role, fold, structure selection, model input, thresholding, or evaluation strata",
            "status": "QUARANTINED_METADATA_FIELD",
        },
        {
            "audit_id": "GSE274103_ACELLULAR_INTERFACE_FILTER",
            "logical_unit_id": "GEO::GSE274103",
            "risk_type": "PREPROCESSING_SELECTION_BIAS",
            "source_path": "infra/bioinf-data-index/raw/geo/GSE274103_gsm_quick.soft",
            "metadata_locator": "!Sample_data_processing",
            "observed_fact": "upstream processing excluded spots in acellular stromal interfaces",
            "control": "cannot support unbiased tumor-stroma boundary prevalence or field claims without raw-space audit",
            "status": "OPEN_LIMITATION",
        },
        {
            "audit_id": "GSE274557_ACELLULAR_INTERFACE_FILTER",
            "logical_unit_id": "GEO::GSE274557",
            "risk_type": "PREPROCESSING_SELECTION_BIAS",
            "source_path": "infra/bioinf-data-index/raw/geo/GSE274557_gsm_quick.soft",
            "metadata_locator": "!Sample_data_processing",
            "observed_fact": "upstream processing excluded spots in acellular stromal interfaces",
            "control": "cannot support unbiased tumor-stroma boundary prevalence or field claims without raw-space audit",
            "status": "OPEN_LIMITATION",
        },
    ]
    _write_tsv(
        output / "leakage_selection_audit.tsv",
        leakage,
        [
            "audit_id",
            "logical_unit_id",
            "risk_type",
            "source_path",
            "metadata_locator",
            "observed_fact",
            "control",
            "status",
        ],
    )

    policy = {
        "phase": "R-02",
        "mode": "metadata_schema_only",
        "network_allowed": False,
        "gpu_allowed": False,
        "expression_matrix_reads_allowed": False,
        "image_pixel_reads_allowed": False,
        "minimum_auditable_gt": {
            "verified_source_checksum": True,
            "verified_r01_physical_link": True,
            "replayable_geometry_locator": True,
            "known_annotation_provenance": True,
            "known_overlap_boundary": True,
            "independent_from_model_input": True,
        },
        "h5ad_obs_index_coordinate_reads_allowed": True,
        "h5ad_obs_ground_truth_label_reads_allowed": True,
        "supported_confirmatory_gt_verifiers": [
            VERIFIER_CLASS,
            VERIFIER_KIRC,
            VERIFIER_USZ,
            VERIFIER_STCRC,
        ],
        "specimen_equivalent_split_rule": "patient_envelope; block_id remains empty; not block-level",
        "claim_rule": "presence/count/summary alone never creates confirmatory structure instances",
        "second_structure_status": "FROZEN:TUMOR_STROMA_BOUNDARY",
    }
    (output / "r02_policy.json").write_text(
        json.dumps(policy, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    try:
        from scripts.r02_validate_gate import evaluate_gate
    except ModuleNotFoundError:
        # Direct ``python scripts/r02_build_registry.py`` execution places the
        # scripts directory, rather than the repository root, on sys.path.
        from r02_validate_gate import evaluate_gate

    gate = evaluate_gate(root, output)
    (output / "r02_gate.json").write_text(
        json.dumps(gate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "README.md").write_text(
        """# R-02 structure registry

This is a metadata/schema-only, fail-closed control plane. It does not read
expression matrices or image pixels, does not use the network, and does not
use a GPU. The only h5ad access is the per-spot barcode/index coordinate
columns (`obs/_index`, `obs/array_row`, `obs/array_col`) needed to replay
annotation geometry (D-033) plus the per-spot human `ground_truth` label
column of the eight USZ Zenodo 14620362 h5ad files (D-038).

The scoped TLS inventory contains 87 source-reported TLS IDs: HTAN
Vanderbilt CRC 44, GSE226997 4, GSE274103 8, GSE274557 1 (Atlas Table S4)
and GSE175540 30 (Meylan 2022 Table S4 R-IDs). These are cross-checked
against Table S2 counts rather than expanded from count ordinals and remain
non-confirmatory inventory rows.

Auditable instance-level GT now exists for the training-lineage HTAN
Vanderbilt CRC unit and for three public-deposit validation lineages
(D-036/D-037/D-038):

- HTAN Vanderbilt CRC: Heiser et al. 2023 (PMC10756562) commit-pinned
  per-spot pathology annotations, crosswalked piece-by-piece through the
  GitHub sample key and the R-01 meta2 evidence asset, replayed to array
  coordinates through the local OSF h5ad spot indexes. 42 piece/rule GT
  sources are auditable; one piece is replay-blocked (missing local h5ad)
  and one annotation capture is unlinked (its R-01 unit is excluded).
- GSE175540 (Meylan et al. 2022, Immunity, PMID 35231421): author-deposited
  per-spot TLS annotation CSVs in the GEO raw data, replayed through the
  local Space Ranger tissue positions. 18 sample-level sources are
  auditable (35 TLS components, 18 patients); three T_agg-only files are
  registered NOT_AUDITABLE, one sample has no deposited annotation file,
  and two zero-TLS samples yield no audit rows.
- TLS_VISIUM_USZ (Zenodo 14620362): per-spot human `ground_truth` labels
  inside the eight deposited h5ad files. 8 sources are auditable (108 TLS
  components, 8 patients).
- ST_CRC_CMS (Zenodo 7760264, Valdeolivas et al. 2024): pathologist
  per-spot category CSVs for 14 sections from 7 CRC patients. 12 sources
  are auditable (551 tumor-stroma boundary components); two replicate
  sections without boundary labels yield no audit rows.

Together these give 889 confirmatory structure instances: 147 TLS (HTAN 4,
GSE175540 35, USZ 108) and 742 tumor-stroma boundary components (HTAN 191,
ST_CRC_CMS 551). Instances in preneoplastic or normal-mucosa context remain
registered as non-confirmatory context only. GSE226997, GSE274103 and
GSE274557 still have no auditable GT; vasculature and necrosis have no
public GT in any examined source.

Nineteen local STOmicsDB TLS annotation files are inventoried by path,
checksum, and compressed CSV header only. They are outside the R-01 frozen
core and are excluded from the 87 scoped summaries. No STOmics row is
accepted as confirmatory GT.

Outer splits use a conservative patient-wide envelope, including all explicit
blocks from the same patient. Explicit blocks remain separately registered for
block-level counting. R-01 patient-linked physical-specimen exceptions stay
block-empty and are not block-level evidence. The one
explicit TENX section pair is registered as an identity link, not a confirmed
cross-section structure correspondence.

The stored gate is expected to be `PARTIAL_GT_READY`, with TLS and
TUMOR_STROMA_BOUNDARY frozen as claim-bearing structures for the training
lineage and the three public-deposit validation lineages. The
`SECOND_CLAIM_BEARING_STRUCTURE_NOT_FROZEN` blocker is cleared; the
remaining validation-lineage GT gaps are tracked in `docs/ISSUES.md`.
""",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--root", type=Path, default=root)
    parser.add_argument(
        "--output", type=Path, default=root / "infra/structure-registry"
    )
    args = parser.parse_args()
    build_registry(args.root, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
