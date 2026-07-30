#!/usr/bin/env python3
"""Prepare an outcome-blind request queue for metadata needed to unblock R-01."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import tempfile
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping


FIELDS = [
    "request_candidate_id",
    "priority_tier",
    "source_namespace",
    "logical_unit_id",
    "provenance_group_id",
    "provenance_aliases",
    "provenance_group_size",
    "known_current_eligible_unit",
    "public_locator",
    "reference_or_title",
    "local_metadata_status",
    "physical_records",
    "patient_known_records",
    "distinct_patient_ids",
    "conflict_records",
    "missing_required_fields",
    "requested_artifact",
    "independence_status",
    "request_status",
]


def _hash(value: str, length: int = 16) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()[:length]


def _public(value: str) -> bool:
    normalized = value.strip().lower()
    return bool(normalized) and normalized not in {
        "-",
        "none",
        "none (internal)",
        "in-house",
    }


def _provenance_aliases(*values: str) -> list[str]:
    text = " | ".join(values)
    aliases = {
        match.upper()
        for match in re.findall(
            r"\b(?:GSE\d+|EGAD\d+|PRJCA\d+|PMID\s*:\s*\d+)\b",
            text,
            flags=re.IGNORECASE,
        )
    }
    aliases.update(
        f"PMID:{match}"
        for match in re.findall(
            r"(?:pubmed\.ncbi\.nlm\.nih\.gov|ncbi\.nlm\.nih\.gov/pubmed)/(\d+)",
            text,
            flags=re.IGNORECASE,
        )
    )
    aliases.update(
        f"DOI:{match.rstrip('.,;')}"
        for match in re.findall(
            r"\b10\.\d{4,9}/[^\s?#|]+", text, flags=re.IGNORECASE
        )
    )
    return sorted(re.sub(r"\s+", "", alias.upper()) for alias in aliases)


def build_candidates(
    physical_units: Iterable[Mapping[str, str]],
    atlas_units: Iterable[Mapping[str, str]],
    hest_evidence: Iterable[Mapping[str, str]],
    official_conflicts: Mapping[str, Mapping[str, str]] | None = None,
) -> list[dict[str, str]]:
    """Summarize identity gaps without using structure labels or outcomes."""
    official_conflicts = official_conflicts or {}
    atlas = {row["logical_unit_id"]: row for row in atlas_units}
    hest: dict[str, list[Mapping[str, str]]] = {}
    for row in hest_evidence:
        raw_study = row.get("raw_study_link", "").strip()
        if raw_study:
            hest.setdefault(f"HEST_STUDY::{_hash(raw_study)}", []).append(row)

    grouped: dict[tuple[str, str], list[Mapping[str, str]]] = {}
    for row in physical_units:
        study = row.get("study_id", "").strip()
        if study:
            grouped.setdefault((row.get("source_namespace", ""), study), []).append(row)

    output: list[dict[str, str]] = []
    for (namespace, study), rows in grouped.items():
        official_conflict = official_conflicts.get(study)
        if official_conflict:
            accession = official_conflict.get("accession", "")
            reference = official_conflict.get("official_reference", "")
            output.append(
                {
                    "request_candidate_id": (
                        f"R01META::{_hash(namespace + '|' + study)}"
                    ),
                    "priority_tier": "P5_SCREENED_CONFLICT",
                    "source_namespace": namespace,
                    "logical_unit_id": study,
                    "provenance_group_id": "",
                    "provenance_aliases": ";".join(
                        _provenance_aliases(accession, reference)
                    ),
                    "provenance_group_size": "",
                    "known_current_eligible_unit": "",
                    "public_locator": accession,
                    "reference_or_title": reference,
                    "local_metadata_status": "OFFICIAL_CONFLICT_QUARANTINED",
                    "physical_records": official_conflict.get(
                        "official_record_count",
                        "0",
                    ),
                    "patient_known_records": "0",
                    "distinct_patient_ids": official_conflict.get(
                        "official_patient_count",
                        "",
                    ),
                    "conflict_records": str(len(rows)),
                    "missing_required_fields": "official_metadata_conflict",
                    "requested_artifact": "none; screened conflict excluded",
                    "independence_status": (
                        "OFFICIAL_METADATA_CONFLICT_QUARANTINED"
                    ),
                    "request_status": "SCREENED_CONFLICT_EXCLUDED",
                }
            )
            continue
        patient_rows = [
            row
            for row in rows
            if row.get("patient_id")
            and not row.get("block_id")
            and row.get("block_equivalent_status")
            != "ACCEPTED_BLOCK_EQUIVALENT"
            and row.get("evidence_grade") in {"E3_explicit", "E2_corroborated"}
            and row.get("record_status") != "QUARANTINED_METADATA_CONFLICT"
        ]
        if not patient_rows or any(row.get("block_id") for row in rows):
            continue

        conflicts = sum(
            row.get("evidence_grade") == "EC_conflict"
            or row.get("record_status") == "QUARANTINED_METADATA_CONFLICT"
            for row in rows
        )
        locator = ""
        reference = ""
        local_status = "BUNDLED_REGISTRY_ONLY"
        if namespace == "ATLAS_TABLE_S2" and study in atlas:
            source = atlas[study]
            locator = source.get("data_availability", "").strip()
            reference = source.get("references", "").strip()
            local_status = source.get("local_data_status", "") or local_status
        elif namespace == "HEST" and study in hest:
            source_rows = hest[study]
            locator = source_rows[0].get("raw_study_link", "").strip()
            titles = sorted(
                {
                    row.get("raw_dataset_title", "").strip()
                    for row in source_rows
                    if row.get("raw_dataset_title", "").strip()
                }
            )
            reference = " | ".join(titles)
            local_status = "BUNDLED_HEST_METADATA_PRESENT"

        public_locator = _public(locator) or _public(reference)
        distinct_patients = len({row["patient_id"] for row in patient_rows})
        if public_locator and distinct_patients >= 4 and conflicts == 0:
            priority = "P1_PUBLIC_HIGH_COVERAGE"
        elif public_locator and conflicts == 0:
            priority = "P2_PUBLIC_LOWER_COVERAGE"
        else:
            priority = "P3_ACCESS_OR_CONFLICT_REVIEW"

        output.append(
            {
                "request_candidate_id": f"R01META::{_hash(namespace + '|' + study)}",
                "priority_tier": priority,
                "source_namespace": namespace,
                "logical_unit_id": study,
                "provenance_group_id": "",
                "provenance_aliases": ";".join(
                    _provenance_aliases(locator, reference)
                ),
                "provenance_group_size": "",
                "known_current_eligible_unit": (
                    "HTAN_VANDERBILT_CRC"
                    if "humantumoratlas.org/publications/vanderbilt_crc_chen_2021"
                    in locator.lower()
                    else ""
                ),
                "public_locator": locator,
                "reference_or_title": reference,
                "local_metadata_status": local_status,
                "physical_records": str(len(rows)),
                "patient_known_records": str(len(patient_rows)),
                "distinct_patient_ids": str(distinct_patients),
                "conflict_records": str(conflicts),
                "missing_required_fields": "block_id;sample_to_block_crosswalk",
                "requested_artifact": (
                    "official sample sheet or provenance-preserving "
                    "patient-block-sample crosswalk"
                ),
                "independence_status": "UNRESOLVED_DO_NOT_COUNT",
                "request_status": "PROPOSED_REQUIRES_APPROVAL",
            }
        )

    parents = list(range(len(output)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parents[right_root] = left_root

    alias_owner: dict[str, int] = {}
    for index, row in enumerate(output):
        for alias in filter(None, row["provenance_aliases"].split(";")):
            if alias in alias_owner:
                union(index, alias_owner[alias])
            else:
                alias_owner[alias] = index

    groups: dict[int, list[int]] = {}
    for index in range(len(output)):
        groups.setdefault(find(index), []).append(index)
    for members in groups.values():
        group_key = "|".join(
            sorted(output[index]["request_candidate_id"] for index in members)
        )
        group_id = f"R01PROV::{_hash(group_key)}"
        known_units = sorted(
            {
                output[index]["known_current_eligible_unit"]
                for index in members
                if output[index]["known_current_eligible_unit"]
            }
        )
        for index in members:
            row = output[index]
            row["provenance_group_id"] = group_id
            row["provenance_group_size"] = str(len(members))
            if known_units:
                row["known_current_eligible_unit"] = ";".join(known_units)
                row["priority_tier"] = "P4_KNOWN_CURRENT_LINEAGE"
                row["independence_status"] = (
                    "KNOWN_CURRENT_ELIGIBLE_LINEAGE_NO_INCREMENT"
                )
                row["request_status"] = "NOT_AN_INCREMENTAL_UNIT"
            elif len(members) > 1:
                row["independence_status"] = (
                    "CROSS_AGGREGATOR_OVERLAP_UNRESOLVED_DO_NOT_COUNT"
                )

    tier_order = {
        "P1_PUBLIC_HIGH_COVERAGE": 1,
        "P2_PUBLIC_LOWER_COVERAGE": 2,
        "P3_ACCESS_OR_CONFLICT_REVIEW": 3,
        "P4_KNOWN_CURRENT_LINEAGE": 4,
        "P5_SCREENED_CONFLICT": 5,
    }
    return sorted(
        output,
        key=lambda row: (
            tier_order[row["priority_tier"]],
            -int(row["distinct_patient_ids"] or 0),
            row["source_namespace"],
            row["logical_unit_id"],
        ),
    )


def build_summary(
    candidates: list[Mapping[str, str]],
    gate: Mapping[str, object],
    approval_status: str = "NOT_APPROVED",
) -> dict[str, object]:
    tiers = Counter(row["priority_tier"] for row in candidates)
    provenance_groups = {
        row["provenance_group_id"]
        for row in candidates
        if row.get("provenance_group_id")
    }
    known_groups = {
        row["provenance_group_id"]
        for row in candidates
        if row.get("known_current_eligible_unit")
    }
    screened_conflict_groups = {
        row["provenance_group_id"]
        for row in candidates
        if row.get("request_status") == "SCREENED_CONFLICT_EXCLUDED"
        and row.get("provenance_group_id")
    }
    return {
        "node": "R-01",
        "current_gate_status": gate.get("status"),
        "current_claim_eligible_logical_units": gate.get(
            "claim_eligible_logical_units"
        ),
        "minimum_additional_independent_logical_units": max(
            0, 6 - int(gate.get("claim_eligible_logical_units", 0))
        ),
        "request_candidates": len(candidates),
        "provenance_groups_before_independence_audit": len(provenance_groups),
        "known_current_lineage_groups_no_increment": len(known_groups),
        "screened_conflict_groups_excluded": len(screened_conflict_groups),
        "proposed_provenance_groups_requiring_audit": len(
            provenance_groups - known_groups - screened_conflict_groups
        ),
        "priority_tier_counts": dict(sorted(tiers.items())),
        "request_scope": (
            "official metadata mapping patient, tissue block, sample and section"
        ),
        "excluded_scope": [
            "expression matrices",
            "images",
            "structure labels selected by outcome",
            "model-derived annotations",
        ],
        "ranking_inputs": [
            "bundled patient identity coverage",
            "bundled provenance locator availability",
            "bundled metadata conflict count",
        ],
        "independence_policy": (
            "Every candidate remains unresolved and contributes zero units until "
            "cross-source duplicate and patient/block lineage audit passes."
        ),
        "approval_status": approval_status,
    }


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent, text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--physical-units",
        type=Path,
        default=root / "infra/sample-registry/physical_units.tsv",
    )
    parser.add_argument(
        "--atlas-units",
        type=Path,
        default=root / "infra/sample-registry/staging/atlas_logical_units.tsv",
    )
    parser.add_argument(
        "--hest-evidence",
        type=Path,
        default=root / "infra/sample-registry/staging/hest_identity_evidence.tsv",
    )
    parser.add_argument(
        "--official-conflicts",
        type=Path,
        default=(
            root
            / "infra/sample-registry/staging/official_metadata_conflicts.tsv"
        ),
    )
    parser.add_argument(
        "--gate",
        type=Path,
        default=root / "infra/sample-registry/r01_gate.json",
    )
    parser.add_argument(
        "--output-tsv",
        type=Path,
        default=root / "infra/sample-registry/r01_metadata_request.tsv",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=root / "infra/sample-registry/r01_metadata_request.json",
    )
    parser.add_argument(
        "--approval-status",
        default="NOT_APPROVED",
        choices=[
            "NOT_APPROVED",
            "APPROVED_METADATA_ONLY",
            "APPROVED_METADATA_ONLY_COMPLETED",
        ],
    )
    args = parser.parse_args()

    for output in (args.output_tsv, args.output_json):
        try:
            output.resolve().relative_to(root.resolve())
        except ValueError as exc:
            raise ValueError("output escapes project root") from exc

    candidates = build_candidates(
        _read_tsv(args.physical_units),
        _read_tsv(args.atlas_units),
        _read_tsv(args.hest_evidence),
        official_conflicts={
            row["study_id"]: row
            for row in _read_tsv(args.official_conflicts)
        },
    )
    if args.approval_status != "NOT_APPROVED":
        for candidate in candidates:
            if candidate["request_status"] == "PROPOSED_REQUIRES_APPROVAL":
                candidate["request_status"] = args.approval_status
    with args.gate.open(encoding="utf-8") as handle:
        summary = build_summary(
            candidates,
            json.load(handle),
            approval_status=args.approval_status,
        )

    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=FIELDS, delimiter="\t")
    writer.writeheader()
    writer.writerows(candidates)
    _atomic_write(args.output_tsv, buffer.getvalue())
    _atomic_write(
        args.output_json,
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )
    print(
        f"metadata request candidates={len(candidates)} "
        f"minimum_needed={summary['minimum_additional_independent_logical_units']}"
    )


if __name__ == "__main__":
    main()
