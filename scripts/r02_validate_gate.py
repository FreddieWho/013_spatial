#!/usr/bin/env python3
"""Independently recompute the R-02 fail-closed gate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Mapping

from openpyxl import load_workbook


ATLAS_PATH = Path("paper/tables/science.adz2742_tables_s1_to_s8.xlsx")
EXPECTED_TLS_COUNTS = {
    "HTAN_VANDERBILT_CRC": 44,
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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _within_root(root: Path, relative: str) -> Path | None:
    try:
        path = (root / relative).resolve()
        path.relative_to(root.resolve())
    except (ValueError, OSError):
        return None
    return path


def _official_tls_inventory(
    root: Path,
) -> tuple[
    dict[str, int],
    dict[str, dict[str, str]],
    dict[str, tuple[int, int]],
    list[str],
]:
    workbook = load_workbook(root / ATLAS_PATH, read_only=True, data_only=True)
    summary_sheet = workbook["Table S2"]
    lookup = {
        HTAN_AVAILABILITY: "HTAN_VANDERBILT_CRC",
        "GSE226997": "GEO::GSE226997",
        "GSE274103": "GEO::GSE274103",
        "GSE274557": "GEO::GSE274557",
    }
    counts = {key: 0 for key in EXPECTED_TLS_COUNTS}
    samples: dict[str, tuple[str, int]] = {}
    scope_collisions: list[str] = []
    current_availability = ""
    for values in summary_sheet.iter_rows(min_row=4, values_only=True):
        if values[11] is not None:
            current_availability = str(values[11]).strip()
        unit = lookup.get(current_availability)
        if unit:
            count = int(values[9] or 0)
            counts[unit] += count
            sample_key = str(values[3] or "")
            if sample_key in samples:
                scope_collisions.append(
                    f"{sample_key}:{samples[sample_key][0]}:{unit}"
                )
            else:
                samples[sample_key] = (unit, count)
    details: dict[str, dict[str, str]] = {}
    detail_counts: Counter[str] = Counter()
    for row_number, values in enumerate(
        workbook["Table S4"].iter_rows(min_row=4, values_only=True), start=4
    ):
        if values[0] is None:
            continue
        tls_id = str(values[0]).strip()
        sample_key = tls_id.rsplit("_", 1)[0]
        scope = samples.get(sample_key)
        if scope is None:
            continue
        unit, _ = scope
        detail_counts[sample_key] += 1
        details[tls_id] = {
            "logical_unit_id": unit,
            "source_record_id": f"Table S4!A{row_number}:D{row_number}",
            "source_sample_value": sample_key,
            "gt_source_id": "ATLAS_TABLE_S4_TLS_IDS",
        }
    workbook.close()
    mismatches = {
        sample: (count, detail_counts[sample])
        for sample, (_, count) in samples.items()
        if count != detail_counts[sample]
    }
    return counts, details, mismatches, scope_collisions


def _source_errors(
    root: Path,
    audits: list[dict[str, str]],
    eligible_physical: Mapping[str, dict[str, str]],
    supported_verifiers: set[str],
) -> tuple[list[str], dict[str, str]]:
    errors: list[str] = []
    valid_auditable: dict[str, str] = {}
    seen: set[str] = set()
    for row in audits:
        source_id = row["gt_source_id"]
        if source_id in seen:
            errors.append(f"duplicate_gt_source:{source_id}")
            continue
        seen.add(source_id)
        path = _within_root(root, row["path"])
        if path is None:
            errors.append(f"path_escape:{source_id}")
            checksum_ok = False
        elif not path.is_file():
            errors.append(f"missing_source:{source_id}")
            checksum_ok = False
        else:
            checksum_ok = (
                row["checksum_status"] == "VERIFIED"
                and len(row["sha256"]) == 64
                and _sha256(path) == row["sha256"]
            )
            if not checksum_ok:
                errors.append(f"source_checksum_mismatch:{source_id}")

        if row["audit_status"] != "AUDITABLE_GT":
            continue
        source_claim_errors: list[str] = []
        if row["source_class"] not in supported_verifiers:
            source_claim_errors.append(f"unsupported_gt_verifier:{source_id}")
        if not checksum_ok:
            source_claim_errors.append(f"unverified_source:{source_id}")
        physical_id = row["physical_unit_id"]
        if (
            row["physical_link_status"] != "VERIFIED"
            or not physical_id
            or physical_id not in eligible_physical
        ):
            source_claim_errors.append(f"missing_physical_link:{source_id}")
        if not row["geometry_locator"]:
            source_claim_errors.append(f"missing_geometry:{source_id}")
        elif not row["geometry_locator"].startswith(row["path"] + "#"):
            source_claim_errors.append(f"nonreplayable_geometry:{source_id}")
        if row["provenance_status"] != "KNOWN":
            source_claim_errors.append(f"unknown_provenance:{source_id}")
        elif not row["provenance_locator"]:
            source_claim_errors.append(f"missing_provenance_locator:{source_id}")
        if row["overlap_status"] not in {"KNOWN_NONOVERLAP", "KNOWN_OVERLAP"}:
            source_claim_errors.append(f"unknown_overlap:{source_id}")
        elif not row["overlap_locator"]:
            source_claim_errors.append(f"missing_overlap_locator:{source_id}")
        if row["same_assay_status"] != "INDEPENDENT_ASSAY":
            source_claim_errors.append(f"circular_same_assay_gt:{source_id}")
        if row["source_class"] in {"ATLAS_TLS_COUNT_SUMMARY"}:
            source_claim_errors.append(f"summary_only_not_geometry:{source_id}")
        if source_claim_errors:
            errors.extend(source_claim_errors)
        else:
            valid_auditable[source_id] = physical_id
    return errors, valid_auditable


def _split_errors(
    splits: list[dict[str, str]],
    eligible_physical: Mapping[str, dict[str, str]],
    roles: Mapping[str, dict[str, str]],
    links: list[dict[str, str]],
) -> list[str]:
    errors: list[str] = []
    split_by_physical: dict[str, dict[str, str]] = {}
    envelopes: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for row in splits:
        physical_id = row["physical_unit_id"]
        if physical_id in split_by_physical:
            errors.append(f"duplicate_split_row:{physical_id}")
        split_by_physical[physical_id] = row
        source = eligible_physical.get(physical_id)
        if source is None:
            errors.append(f"unknown_split_physical:{physical_id}")
            continue
        for field in ("patient_id", "block_id", "physical_specimen_id"):
            if row[field] != source[field]:
                errors.append(f"split_identity_mismatch:{physical_id}:{field}")
        role = roles.get(row["logical_unit_id"])
        if (
            source["study_id"] != row["logical_unit_id"]
            or role is None
            or row["primary_role"] != role["primary_role"]
            or row["leakage_group_id"] != role["leakage_group_id"]
        ):
            errors.append(f"split_role_mismatch:{physical_id}")
        if source["identity_granularity"] == "patient_linked_physical_specimen":
            expected = f"patient_envelope::{source['study_id']}::{source['patient_id']}"
            if row["block_id"]:
                errors.append(f"specimen_as_block:{physical_id}")
            if (
                row["identity_envelope_id"] != expected
                or row["block_group_id"]
                or row["block_level_eligible"] != "no"
            ):
                errors.append(f"invalid_patient_envelope:{physical_id}")
        elif source["block_id"]:
            expected = f"patient_envelope::{source['study_id']}::{source['patient_id']}"
            if (
                row["identity_envelope_id"] != expected
                or row["block_group_id"] != f"block::{source['block_id']}"
                or row["block_level_eligible"] != "yes"
            ):
                errors.append(f"invalid_block_envelope:{physical_id}")
        envelopes[row["identity_envelope_id"]].add(
            (row["primary_role"], row["outer_fold"])
        )
    missing = sorted(set(eligible_physical) - set(split_by_physical))
    errors.extend(f"missing_split_row:{item}" for item in missing)
    for envelope, placements in envelopes.items():
        if len(placements) > 1:
            errors.append(f"cross_fold_identity_envelope:{envelope}")
    leakage_placements: dict[str, set[str]] = defaultdict(set)
    for row in splits:
        leakage_placements[row["leakage_group_id"]].add(row["primary_role"])
    for leakage_group, placements in leakage_placements.items():
        if len(placements) > 1:
            errors.append(f"cross_role_leakage_group:{leakage_group}")

    seen_link_ids: set[str] = set()
    seen_link_pairs: set[tuple[str, str]] = set()
    for link in links:
        link_id = link["link_id"]
        pair = tuple(
            sorted((link["left_physical_unit_id"], link["right_physical_unit_id"]))
        )
        if link_id in seen_link_ids or pair in seen_link_pairs:
            errors.append(f"duplicate_cross_section_link:{link_id}")
        seen_link_ids.add(link_id)
        seen_link_pairs.add(pair)
        left = split_by_physical.get(link["left_physical_unit_id"])
        right = split_by_physical.get(link["right_physical_unit_id"])
        if left is None or right is None:
            errors.append(f"cross_section_unknown_physical:{link_id}")
            continue
        left_source = eligible_physical[link["left_physical_unit_id"]]
        right_source = eligible_physical[link["right_physical_unit_id"]]
        if (
            link["patient_id"] != left["patient_id"]
            or link["patient_id"] != right["patient_id"]
            or link["block_id"] != left["block_id"]
            or link["block_id"] != right["block_id"]
            or link["relation_type"] != "SAME_EXPLICIT_BLOCK_SECTION_PAIR"
            or link["evidence_status"]
            != "IDENTITY_ONLY_NOT_STRUCTURE_CORRESPONDENCE"
            or not left_source["section_id"]
            or not right_source["section_id"]
            or left_source["section_id"] == right_source["section_id"]
        ):
            errors.append(f"cross_section_link_mismatch:{link_id}")
        if (
            left["patient_id"] != right["patient_id"]
            or not left["block_id"]
            or left["block_id"] != right["block_id"]
            or left["primary_role"] != right["primary_role"]
            or left["outer_fold"] != right["outer_fold"]
            or link["outer_fold"] != left["outer_fold"]
            or link["cross_fold_allowed"] != "no"
        ):
            errors.append(f"cross_section_leakage:{link_id}")
    return errors


def _duplicate_link_errors(
    root: Path,
    splits: list[dict[str, str]],
    eligible_physical: Mapping[str, dict[str, str]],
) -> list[str]:
    split_by_physical = {row["physical_unit_id"]: row for row in splits}
    source_to_split = {
        source["source_record_id"]: split_by_physical[physical_id]
        for physical_id, source in eligible_physical.items()
        if physical_id in split_by_physical and source["source_record_id"]
    }
    placements: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for member in _read_tsv(root / "infra/sample-registry/duplicate_groups.tsv"):
        split = source_to_split.get(member["member_id"])
        if split:
            placements[member["duplicate_group_id"]].add(
                (split["primary_role"], split["outer_fold"])
            )
    return [
        f"cross_fold_duplicate_group:{group_id}"
        for group_id, group_placements in placements.items()
        if len(group_placements) > 1
    ]


def _metadata_audit_errors(
    root: Path, registry: Path, source_assets: list[dict[str, str]]
) -> list[str]:
    rows = _read_tsv(registry / "leakage_selection_audit.tsv")
    by_id = {row["audit_id"]: row for row in rows}
    asset_by_path = {row["path"]: row for row in source_assets}
    expected = {
        "GSE211956_RESPONSE_TITLE": (
            "infra/bioinf-data-index/raw/geo/GSE211956_gsm_quick.soft",
            "response_",
        ),
        "GSE274103_ACELLULAR_INTERFACE_FILTER": (
            "infra/bioinf-data-index/raw/geo/GSE274103_gsm_quick.soft",
            "acellular stromal interfaces were excluded",
        ),
        "GSE274557_ACELLULAR_INTERFACE_FILTER": (
            "infra/bioinf-data-index/raw/geo/GSE274557_gsm_quick.soft",
            "acellular stromal interfaces were excluded",
        ),
    }
    errors: list[str] = []
    for audit_id, (relative, needle) in expected.items():
        row = by_id.get(audit_id)
        if row is None or row["source_path"] != relative:
            errors.append(f"missing_leakage_selection_audit:{audit_id}")
            continue
        path = _within_root(root, relative)
        asset = asset_by_path.get(relative)
        if (
            path is None
            or not path.is_file()
            or asset is None
            or asset["checksum_status"] != "sha256_verified"
            or _sha256(path) != asset["checksum"]
        ):
            errors.append(f"unverified_leakage_audit_source:{audit_id}")
            continue
        if needle not in path.read_text(encoding="utf-8", errors="replace"):
            errors.append(f"leakage_audit_fact_not_found:{audit_id}")
    return errors


def evaluate_gate(root: Path, registry: Path) -> dict[str, object]:
    root = root.resolve()
    physical_rows = _read_tsv(root / "infra/sample-registry/physical_units.tsv")
    source_assets = _read_tsv(root / "infra/sample-registry/source_assets.tsv")
    eligible = {
        row["physical_unit_id"]: row
        for row in physical_rows
        if row["record_status"] == "RESOLVED_INCLUDED_CANDIDATE"
    }
    roles = {
        row["logical_unit_id"]: row
        for row in _read_tsv(root / "infra/sample-registry/role_freeze.tsv")
        if row["record_status"] == "FROZEN"
    }
    audits = _read_tsv(registry / "gt_source_audit.tsv")
    instances = _read_tsv(registry / "structure_instances.tsv")
    candidates = _read_tsv(registry / "tls_candidate_summary.tsv")
    ontology = _read_tsv(registry / "structure_ontology.tsv")
    splits = _read_tsv(registry / "outer_splits.tsv")
    links = _read_tsv(registry / "cross_section_links.tsv")
    input_policy = _read_tsv(registry / "input_policy.tsv")
    policy = json.loads((registry / "r02_policy.json").read_text(encoding="utf-8"))

    # Verifiers are executable code, not a table assertion. R-02 implements none.
    implemented_verifiers: set[str] = set()
    declared_verifiers = set(policy.get("supported_confirmatory_gt_verifiers", []))
    if declared_verifiers != implemented_verifiers:
        integrity_errors = [
            "policy_declares_unimplemented_gt_verifier:"
            + ",".join(sorted(declared_verifiers))
        ]
    else:
        integrity_errors = []
    supported_verifiers = implemented_verifiers
    source_integrity_errors, auditable_sources = _source_errors(
        root, audits, eligible, supported_verifiers
    )
    integrity_errors.extend(source_integrity_errors)
    integrity_errors.extend(_metadata_audit_errors(root, registry, source_assets))
    split_errors = _split_errors(splits, eligible, roles, links)
    split_errors.extend(_duplicate_link_errors(root, splits, eligible))
    audit_by_id = {row["gt_source_id"]: row for row in audits}

    required_policies = {
        ("expression matrix or derived expression", "FORBIDDEN_R02"),
        ("image pixels or derived morphology", "FORBIDDEN_R02"),
        ("GT-defining annotation", "FORBIDDEN_AS_MODEL_INPUT"),
        ("outcome or response metadata", "FORBIDDEN_FOR_SELECTION_SPLIT_MODEL"),
        (
            "identity, site, treatment, stage, sample and file metadata",
            "FORBIDDEN_AS_PREDICTIVE_INPUT",
        ),
    }
    observed_policies = {
        (row["channel_class"], row["policy"]) for row in input_policy
    }
    for channel, policy in sorted(required_policies - observed_policies):
        integrity_errors.append(f"missing_input_policy:{channel}:{policy}")

    (
        official_counts,
        official_ids,
        source_sample_mismatches,
        source_scope_collisions,
    ) = _official_tls_inventory(root)
    candidate_counts = Counter(row["logical_unit_id"] for row in candidates)
    if official_counts != EXPECTED_TLS_COUNTS:
        integrity_errors.append(
            f"source_tls_count_mismatch:expected={EXPECTED_TLS_COUNTS}:observed={official_counts}"
        )
    if dict(candidate_counts) != EXPECTED_TLS_COUNTS or len(candidates) != 57:
        integrity_errors.append(
            f"candidate_tls_count_mismatch:expected={EXPECTED_TLS_COUNTS}:observed={dict(candidate_counts)}"
        )
    if source_sample_mismatches:
        integrity_errors.append(
            f"table_s2_s4_sample_count_mismatch:{source_sample_mismatches}"
        )
    if source_scope_collisions:
        integrity_errors.append(
            f"ambiguous_table_s2_sample_scope:{source_scope_collisions}"
        )
    candidate_by_tls_id = {row["source_tls_id"]: row for row in candidates}
    if len(candidate_by_tls_id) != len(candidates):
        integrity_errors.append("duplicate_source_tls_id")
    if set(candidate_by_tls_id) != set(official_ids):
        integrity_errors.append("table_s4_candidate_id_set_mismatch")
    else:
        for tls_id, expected in official_ids.items():
            row = candidate_by_tls_id[tls_id]
            for field, value in expected.items():
                if row[field] != value:
                    integrity_errors.append(
                        f"table_s4_candidate_field_mismatch:{tls_id}:{field}"
                    )

    with (root / "repo/data/ST_CRC_maturation_location.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        crc_ids = {row["TLS_ID"] for row in csv.DictReader(handle)}
    htan_ids = {
        tls_id
        for tls_id, row in official_ids.items()
        if row["logical_unit_id"] == "HTAN_VANDERBILT_CRC"
    }
    if crc_ids != htan_ids or len(crc_ids) != 44:
        integrity_errors.append("upstream_crc_tls_id_corroboration_mismatch")

    instance_ids: set[str] = set()
    valid_confirmatory: list[dict[str, str]] = []
    for row in instances:
        instance_id = row["instance_id"]
        if instance_id in instance_ids:
            integrity_errors.append(f"duplicate_instance:{instance_id}")
        instance_ids.add(instance_id)
        if row["confirmation_status"] != "CONFIRMATORY":
            continue
        source_id = row["gt_source_id"]
        if source_id not in auditable_sources:
            integrity_errors.append(f"unverified_gt_source:{instance_id}:{source_id}")
            continue
        source = eligible.get(row["physical_unit_id"])
        if source is None:
            integrity_errors.append(f"unverified_instance_physical:{instance_id}")
            continue
        if auditable_sources[source_id] != row["physical_unit_id"]:
            integrity_errors.append(f"gt_source_instance_physical_mismatch:{instance_id}")
            continue
        source_audit = audit_by_id[source_id]
        if row["gt_geometry_locator"] != source_audit["geometry_locator"]:
            integrity_errors.append(f"instance_geometry_not_source_bound:{instance_id}")
            continue
        if not row["gt_geometry_locator"]:
            integrity_errors.append(f"missing_instance_geometry:{instance_id}")
            continue
        for field in ("patient_id", "block_id", "physical_specimen_id"):
            if row[field] != source[field]:
                integrity_errors.append(f"instance_identity_mismatch:{instance_id}:{field}")
                break
        else:
            valid_confirmatory.append(row)

    candidate_ids = {row["candidate_id"] for row in candidates}
    if instance_ids != candidate_ids:
        integrity_errors.append("candidate_instance_id_set_mismatch")

    second_structure_frozen = any(
        row["structure_id"] != "TLS"
        and row["claim_status"] == "FROZEN_CLAIM_BEARING"
        and any(item["structure_id"] == row["structure_id"] for item in valid_confirmatory)
        for row in ontology
    )
    confirmatory_count = len(valid_confirmatory)
    blockers: list[str] = []
    if integrity_errors:
        status = "HARD_BLOCKED_GT_INTEGRITY"
        blockers.append("GT_SOURCE_OR_INSTANCE_INTEGRITY_FAILED")
    elif split_errors:
        status = "HARD_BLOCKED_SPLIT_INTEGRITY"
        blockers.append("OUTER_SPLIT_OR_IDENTITY_LINK_INTEGRITY_FAILED")
    elif confirmatory_count == 0:
        status = "HARD_BLOCKED_NO_AUDITABLE_GT"
        blockers.append("NO_AUDITABLE_INSTANCE_LEVEL_GT")
    else:
        status = "PARTIAL_GT_READY"
    if not second_structure_frozen:
        blockers.append("SECOND_CLAIM_BEARING_STRUCTURE_NOT_FROZEN")

    stomics_count = sum(
        row["source_class"] == "STOMICS_TLS_ANNOTATION_CANDIDATE" for row in audits
    )
    if stomics_count != 19:
        integrity_errors.append(f"stomics_inventory_count_mismatch:{stomics_count}")
        if status == "HARD_BLOCKED_NO_AUDITABLE_GT":
            status = "HARD_BLOCKED_GT_INTEGRITY"
            blockers.insert(0, "GT_SOURCE_OR_INSTANCE_INTEGRITY_FAILED")

    outer_patient_groups = {
        row["identity_envelope_id"] for row in splits if row["identity_envelope_id"]
    }
    explicit_block_groups = {
        row["block_group_id"] for row in splits if row["block_group_id"]
    }

    return {
        "phase": "R-02",
        "status": status,
        "metadata_schema_only": True,
        "gpu_used": False,
        "network_used": False,
        "expression_or_image_content_read": False,
        "ontology_structure_count": len(ontology),
        "tls_candidate_count": len(candidates),
        "tls_candidate_breakdown": dict(sorted(candidate_counts.items())),
        "stomics_candidate_source_count": stomics_count,
        "r01_eligible_physical_count": len(eligible),
        "outer_split_row_count": len(splits),
        "outer_patient_group_count": len(outer_patient_groups),
        "explicit_block_group_count": len(explicit_block_groups),
        "confirmatory_instance_count": confirmatory_count,
        "auditable_gt_source_count": len(auditable_sources),
        "second_structure_frozen": second_structure_frozen,
        "integrity_errors": sorted(set(integrity_errors + split_errors)),
        "blockers": list(dict.fromkeys(blockers)),
        "conclusion_impact": (
            "R-02 remains incomplete: no structure has auditable instance-level GT, "
            "so no distance field, model training, confirmatory validation, or "
            "second-structure claim may start."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    root = Path(__file__).resolve().parents[1]
    parser.add_argument("--root", type=Path, default=root)
    parser.add_argument(
        "--registry", type=Path, default=root / "infra/structure-registry"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    gate = evaluate_gate(args.root, args.registry)
    payload = json.dumps(gate, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")
    return 0 if gate["status"] != "PARTIAL_GT_READY" else 0


if __name__ == "__main__":
    raise SystemExit(main())
