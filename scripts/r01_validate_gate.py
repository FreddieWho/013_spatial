#!/usr/bin/env python3
"""Evaluate the machine-readable completion gate for roadmap node R-01."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping


def _claim_identity_complete(row: Mapping[str, str]) -> bool:
    if row.get("evidence_grade") not in {"E3_explicit", "E2_corroborated"}:
        return False
    equivalence_marked = any(
        row.get(field)
        for field in (
            "physical_specimen_id",
            "identity_granularity",
            "block_equivalent_status",
            "block_equivalent_basis",
        )
    )
    if equivalence_marked:
        return (
            not row.get("block_id")
            and bool(row.get("physical_specimen_id"))
            and row.get("identity_granularity")
            == "patient_linked_physical_specimen"
            and row.get("block_equivalent_status")
            == "ACCEPTED_BLOCK_EQUIVALENT"
            and bool(row.get("block_equivalent_basis"))
            and row.get("identity_status") == "PATIENT_SPECIMEN_CORROBORATED"
        )
    if row.get("block_id"):
        return True
    return False


def _specimen_integrity_errors(
    row: Mapping[str, str],
    evidence_by_entity: Mapping[str, list[Mapping[str, str]]],
    assets_by_id: Mapping[str, Mapping[str, str]],
    duplicate_groups: list[Mapping[str, str]],
    invalid_asset_ids: set[str],
) -> tuple[list[str], str]:
    unit_id = row.get("physical_unit_id", "")
    study_id = row.get("study_id", "")
    accession = study_id.removeprefix("GEO::")
    errors: list[str] = []
    entity_evidence = evidence_by_entity.get(unit_id, [])
    required = {
        "patient_id": row.get("patient_id", ""),
        "physical_specimen_id": row.get("physical_specimen_id", ""),
        "block_id": "",
    }
    for field, expected in required.items():
        matches = [
            evidence
            for evidence in entity_evidence
            if evidence.get("field") == field
            and evidence.get("evidence_grade")
            in {"E3_explicit", "E2_corroborated"}
            and evidence.get("conflict_flag", "").lower()
            in {"", "no", "false"}
        ]
        if not matches:
            errors.append(f"{unit_id}:missing_{field}_evidence")
            continue
        normalized_matches = [
            evidence
            for evidence in matches
            if evidence.get("normalized_value", "") == expected
        ]
        if not normalized_matches:
            errors.append(f"{unit_id}:{field}_mismatch")
            continue
        if field == "block_id" and not any(
            evidence.get("raw_value", "") == row.get(
                "block_equivalent_basis",
                "",
            )
            for evidence in normalized_matches
        ):
            errors.append(f"{unit_id}:block_equivalent_basis_mismatch")
            continue
        for evidence in normalized_matches:
            if not evidence.get("raw_value") or not evidence.get(
                "metadata_key"
            ):
                errors.append(f"{unit_id}:incomplete_raw_evidence:{field}")
            asset_id = evidence.get("source_asset", "")
            asset = assets_by_id.get(asset_id)
            if not asset:
                errors.append(f"{unit_id}:missing_source_asset:{asset_id}")
                continue
            if asset_id in invalid_asset_ids:
                errors.append(f"{unit_id}:invalid_source_asset:{asset_id}")
            if asset.get("accession") != accession:
                errors.append(f"{unit_id}:source_accession_mismatch:{asset_id}")
            if evidence.get("metadata_path") != asset.get("path"):
                errors.append(f"{unit_id}:source_path_mismatch:{asset_id}")
            if asset.get("record_status") != "active":
                errors.append(f"{unit_id}:inactive_source_asset:{asset_id}")
            if (
                asset.get("checksum_status") != "sha256_verified"
                or len(asset.get("checksum", "")) != 64
            ):
                errors.append(
                    f"{unit_id}:unverified_source_asset:{asset_id}"
                )

    memberships = [
        duplicate
        for duplicate in duplicate_groups
        if duplicate.get("member_type") == "study"
        and duplicate.get("member_id") == study_id
        and duplicate.get("relation_type")
        == "same_accession_same_source_lineage"
        and duplicate.get("resolution")
        == "canonical_geo_study; do_not_count_atlas_separately"
        and duplicate.get("leakage_group_id")
    ]
    valid_memberships = []
    for membership in memberships:
        group_id = membership.get("duplicate_group_id", "")
        peers = [
            duplicate
            for duplicate in duplicate_groups
            if duplicate.get("duplicate_group_id") == group_id
            and duplicate.get("member_type") == "study"
            and duplicate.get("member_id") != study_id
        ]
        if peers:
            valid_memberships.append(membership)
    leakage_groups = {
        membership.get("leakage_group_id", "")
        for membership in valid_memberships
    }
    if len(leakage_groups) != 1:
        errors.append(f"{unit_id}:missing_canonical_lineage")
        return sorted(set(errors)), ""
    return sorted(set(errors)), next(iter(leakage_groups))


def _explicit_identity_integrity_errors(
    row: Mapping[str, str],
    evidence_by_entity: Mapping[str, list[Mapping[str, str]]],
    assets_by_id: Mapping[str, Mapping[str, str]],
    duplicate_groups: list[Mapping[str, str]],
    invalid_asset_ids: set[str],
) -> tuple[list[str], str]:
    unit_id = row.get("physical_unit_id", "")
    study_id = row.get("study_id", "")
    errors: list[str] = []
    entity_evidence = evidence_by_entity.get(unit_id, [])
    for field in ("patient_id", "block_id"):
        matches = [
            evidence
            for evidence in entity_evidence
            if evidence.get("field") == field
            and evidence.get("normalized_value") == row.get(field)
            and evidence.get("evidence_grade")
            in {"E3_explicit", "E2_corroborated"}
            and evidence.get("conflict_flag", "").lower()
            in {"", "no", "false"}
        ]
        if not matches:
            errors.append(f"{unit_id}:missing_or_mismatched_{field}_evidence")
            continue
        for evidence in matches:
            if not evidence.get("raw_value") or not evidence.get(
                "metadata_key"
            ):
                errors.append(f"{unit_id}:incomplete_raw_evidence:{field}")
            asset_id = evidence.get("source_asset", "")
            asset = assets_by_id.get(asset_id)
            if not asset:
                errors.append(f"{unit_id}:missing_source_asset:{asset_id}")
                continue
            if asset_id in invalid_asset_ids:
                errors.append(f"{unit_id}:invalid_source_asset:{asset_id}")
            if evidence.get("metadata_path") != asset.get("path"):
                errors.append(f"{unit_id}:source_path_mismatch:{asset_id}")
            if asset.get("record_status") not in {"active", "PRESENT"}:
                errors.append(f"{unit_id}:inactive_source_asset:{asset_id}")
            if (
                asset.get("checksum_status")
                != "sha256_verified"
                or len(asset.get("checksum", "")) != 64
            ):
                errors.append(
                    f"{unit_id}:unverified_source_asset:{asset_id}"
                )

    memberships = [
        duplicate
        for duplicate in duplicate_groups
        if duplicate.get("member_type") == "study"
        and duplicate.get("member_id") == study_id
        and duplicate.get("relation_type") == "canonical_source_lineage"
        and duplicate.get("resolution") == "canonical_study"
        and duplicate.get("leakage_group_id")
    ]
    leakage_groups = {
        membership.get("leakage_group_id", "")
        for membership in memberships
    }
    if len(leakage_groups) != 1:
        errors.append(f"{unit_id}:missing_canonical_source_lineage")
        return sorted(set(errors)), ""
    return sorted(set(errors)), next(iter(leakage_groups))


def evaluate_gate(
    physical_units: Iterable[Mapping[str, str]],
    role_freeze: Iterable[Mapping[str, str]],
    policy: Mapping[str, Any],
    *,
    identity_evidence: Iterable[Mapping[str, str]] = (),
    source_assets: Iterable[Mapping[str, str]] = (),
    duplicate_groups: Iterable[Mapping[str, str]] = (),
    source_asset_errors: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    physical = list(physical_units)
    roles = list(role_freeze)
    evidence_by_entity: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for evidence in identity_evidence:
        evidence_by_entity[evidence.get("entity_id", "")].append(evidence)
    assets_by_id = {
        asset.get("asset_id", ""): asset for asset in source_assets
    }
    duplicate_rows = list(duplicate_groups)
    invalid_asset_ids = set((source_asset_errors or {}).keys())
    integrity_errors: list[str] = []
    canonical_leakage_groups: dict[str, str] = {}
    integrity_complete_ids: set[str] = set()
    for row in physical:
        if (
            row.get("record_status") != "RESOLVED_INCLUDED_CANDIDATE"
            or not _claim_identity_complete(row)
        ):
            continue
        if row.get("block_equivalent_status") == "ACCEPTED_BLOCK_EQUIVALENT":
            errors, leakage_group = _specimen_integrity_errors(
                row,
                evidence_by_entity,
                assets_by_id,
                duplicate_rows,
                invalid_asset_ids,
            )
        else:
            errors, leakage_group = _explicit_identity_integrity_errors(
                row,
                evidence_by_entity,
                assets_by_id,
                duplicate_rows,
                invalid_asset_ids,
            )
        if errors:
            integrity_errors.extend(errors)
        else:
            integrity_complete_ids.add(row.get("physical_unit_id", ""))
            canonical_leakage_groups[row.get("study_id", "")] = leakage_group
    eligible = [
        row
        for row in physical
        if row.get("record_status") == "RESOLVED_INCLUDED_CANDIDATE"
        and row.get("study_id")
        and row.get("patient_id")
        and _claim_identity_complete(row)
        and row.get("physical_unit_id") in integrity_complete_ids
    ]
    logical_units = sorted({row["study_id"] for row in eligible})
    frozen = [row for row in roles if row.get("record_status") == "FROZEN"]
    frozen_units = sorted({row.get("logical_unit_id", "") for row in frozen if row.get("logical_unit_id")})
    frozen_unit_counts = Counter(
        row.get("logical_unit_id", "") for row in frozen if row.get("logical_unit_id")
    )
    allowed_roles = {
        "discovery",
        "training",
        "internal_validation",
        "external_validation",
        "high_resolution_validation",
        "serial_section_validation",
    }
    non_external_groups = {
        row.get("leakage_group_id", "")
        for row in frozen
        if row.get("primary_role") in allowed_roles
        and row.get("primary_role") != "external_validation"
    }
    external_groups = {
        row.get("leakage_group_id", "")
        for row in frozen
        if row.get("primary_role") == "external_validation"
    }
    independent_external = sorted(
        group for group in external_groups if group and group not in non_external_groups
    )
    role_freeze_errors = sorted(
        {
            *(
                f"frozen_unit_not_claim_eligible:{unit_id}"
                for unit_id in frozen_units
                if unit_id not in logical_units
            ),
            *(
                f"multiple_frozen_roles:{unit_id}"
                for unit_id, count in frozen_unit_counts.items()
                if count != 1
            ),
            *(
                f"leakage_group_role_conflict:{group_id}"
                for group_id in non_external_groups & external_groups
                if group_id
            ),
            *(
                f"invalid_primary_role:{row.get('logical_unit_id', '')}"
                for row in frozen
                if row.get("primary_role") not in allowed_roles
            ),
            *(
                f"missing_leakage_group:{row.get('logical_unit_id', '')}"
                for row in frozen
                if not row.get("leakage_group_id")
            ),
            *(
                f"role_leakage_group_mismatch:{row.get('logical_unit_id', '')}"
                for row in frozen
                if row.get("logical_unit_id") in canonical_leakage_groups
                and row.get("leakage_group_id")
                != canonical_leakage_groups[row.get("logical_unit_id", "")]
            ),
        }
    )

    blockers: list[str] = []
    metadata_exhausted = bool(policy.get("metadata_audit_exhausted", False))
    freeze_attempted = bool(policy.get("role_freeze_attempted", False))
    if len(logical_units) < 6 and metadata_exhausted:
        status = "BLOCKED_IDENTITY"
        blockers.append("fewer_than_6_identity_complete_units")
    elif len(logical_units) >= 6 and freeze_attempted and role_freeze_errors:
        status = "BLOCKED_INDEPENDENCE"
        blockers.extend(role_freeze_errors)
    elif (
        len(logical_units) >= 6
        and freeze_attempted
        and not independent_external
    ):
        status = "BLOCKED_INDEPENDENCE"
        blockers.append("no_independent_external_validation_lineage")
    elif (
        6 <= len(logical_units) <= 10
        and len(frozen_units) >= 6
        and independent_external
        and not role_freeze_errors
    ):
        status = (
            "COMPLETE_WITH_EXCLUSIONS"
            if len(eligible) < len(physical)
            else "COMPLETE"
        )
    else:
        status = "IN_PROGRESS"

    return {
        "node": "R-01",
        "status": status,
        "metadata_audit_exhausted": metadata_exhausted,
        "role_freeze_attempted": freeze_attempted,
        "physical_unit_records": len(physical),
        "claim_eligible_physical_units": len(eligible),
        "claim_eligible_logical_units": len(logical_units),
        "claim_eligible_logical_unit_ids": logical_units,
        "frozen_logical_units": len(frozen_units),
        "independent_external_leakage_groups": independent_external,
        "role_freeze_errors": role_freeze_errors,
        "identity_integrity_errors": sorted(set(integrity_errors)),
        "source_asset_errors": dict(sorted((source_asset_errors or {}).items())),
        "blockers": blockers,
        "completion_requirements": {
            "logical_unit_min": 6,
            "logical_unit_max": 10,
            "patient_block_required": False,
            "patient_block_or_accepted_specimen_required": True,
            "accepted_identity_granularities": [
                "patient_block",
                "patient_linked_physical_specimen",
            ],
            "identity_evidence_grades_allowed": [
                "E3_explicit",
                "E2_corroborated",
            ],
            "section_identity_required_by_r01": False,
            "outcome_blind_role_freeze_required": True,
            "independent_external_lineage_required": True,
            "identity_evidence_referential_integrity_required": True,
            "canonical_source_lineage_required": True,
        },
    }


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _verify_source_assets(
    source_assets: Iterable[Mapping[str, str]],
    root: Path,
) -> dict[str, str]:
    errors: dict[str, str] = {}
    for asset in source_assets:
        asset_id = asset.get("asset_id", "")
        path_value = asset.get("path", "")
        if not asset_id or not path_value:
            continue
        if asset.get("checksum_status") not in {
            "sha256_verified",
            "not_computed_large_local_locator",
        }:
            continue
        path = (root / path_value).resolve()
        try:
            path.relative_to(root.resolve())
        except ValueError:
            errors[asset_id] = "path_escapes_project_root"
            continue
        if not path.is_file():
            errors[asset_id] = "missing_file"
            continue
        expected_bytes = asset.get("bytes", "")
        if expected_bytes and path.stat().st_size != int(expected_bytes):
            errors[asset_id] = "byte_size_mismatch"
            continue
        if asset.get("checksum_status") == "sha256_verified":
            digest = hashlib.sha256()
            with path.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            if digest.hexdigest() != asset.get("checksum"):
                errors[asset_id] = "checksum_mismatch"
    return errors


def _write_json(payload: Mapping[str, Any], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = ""
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=output.parent,
            prefix=f".{output.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = handle.name
            json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, output)
    finally:
        if temporary:
            Path(temporary).unlink(missing_ok=True)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--physical-units",
        type=Path,
        default=root / "infra/sample-registry/physical_units.tsv",
    )
    parser.add_argument(
        "--role-freeze",
        type=Path,
        default=root / "infra/sample-registry/role_freeze.tsv",
    )
    parser.add_argument(
        "--identity-evidence",
        type=Path,
        default=root / "infra/sample-registry/identity_evidence.tsv",
    )
    parser.add_argument(
        "--source-assets",
        type=Path,
        default=root / "infra/sample-registry/source_assets.tsv",
    )
    parser.add_argument(
        "--duplicate-groups",
        type=Path,
        default=root / "infra/sample-registry/duplicate_groups.tsv",
    )
    parser.add_argument(
        "--policy",
        type=Path,
        default=root / "infra/sample-registry/r01_gate_policy.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "infra/sample-registry/r01_gate.json",
    )
    args = parser.parse_args()
    output = args.output.resolve()
    try:
        output.relative_to(root)
    except ValueError as exc:
        raise ValueError("output escapes project root") from exc
    input_paths = {
        args.physical_units.resolve(),
        args.role_freeze.resolve(),
        args.identity_evidence.resolve(),
        args.source_assets.resolve(),
        args.duplicate_groups.resolve(),
        args.policy.resolve(),
    }
    if output in input_paths:
        raise ValueError("output collides with gate input")

    with args.policy.open(encoding="utf-8") as handle:
        policy = json.load(handle)
    physical_units = _read_tsv(args.physical_units)
    role_freeze = _read_tsv(args.role_freeze)
    identity_evidence = _read_tsv(args.identity_evidence)
    source_assets = _read_tsv(args.source_assets)
    duplicate_groups = _read_tsv(args.duplicate_groups)
    result = evaluate_gate(
        physical_units,
        role_freeze,
        policy,
        identity_evidence=identity_evidence,
        source_assets=source_assets,
        duplicate_groups=duplicate_groups,
        source_asset_errors=_verify_source_assets(source_assets, root),
    )
    _write_json(result, output)
    print(
        f"R-01 status={result['status']} "
        f"logical_units={result['claim_eligible_logical_units']} "
        f"frozen={result['frozen_logical_units']}"
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
