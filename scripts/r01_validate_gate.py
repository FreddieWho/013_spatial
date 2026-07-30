#!/usr/bin/env python3
"""Evaluate the machine-readable completion gate for roadmap node R-01."""

from __future__ import annotations

import argparse
import csv
import json
import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping


def evaluate_gate(
    physical_units: Iterable[Mapping[str, str]],
    role_freeze: Iterable[Mapping[str, str]],
    policy: Mapping[str, Any],
) -> dict[str, Any]:
    physical = list(physical_units)
    roles = list(role_freeze)
    eligible = [
        row
        for row in physical
        if row.get("record_status") == "RESOLVED_INCLUDED_CANDIDATE"
        and row.get("study_id")
        and row.get("patient_id")
        and row.get("block_id")
        and row.get("evidence_grade") in {"E3_explicit", "E2_corroborated"}
    ]
    logical_units = sorted({row["study_id"] for row in eligible})
    frozen = [row for row in roles if row.get("record_status") == "FROZEN"]
    frozen_units = sorted({row.get("logical_unit_id", "") for row in frozen if row.get("logical_unit_id")})
    frozen_unit_counts = Counter(
        row.get("logical_unit_id", "") for row in frozen if row.get("logical_unit_id")
    )
    training_groups = {
        row.get("leakage_group_id", "")
        for row in frozen
        if row.get("primary_role") in {"discovery", "training", "internal_validation"}
    }
    external_groups = {
        row.get("leakage_group_id", "")
        for row in frozen
        if row.get("primary_role") == "external_validation"
    }
    allowed_roles = {
        "discovery",
        "training",
        "internal_validation",
        "external_validation",
    }
    independent_external = sorted(
        group for group in external_groups if group and group not in training_groups
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
                for group_id in training_groups & external_groups
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
        "blockers": blockers,
        "completion_requirements": {
            "logical_unit_min": 6,
            "logical_unit_max": 10,
            "patient_block_required": True,
            "identity_evidence_grades_allowed": [
                "E3_explicit",
                "E2_corroborated",
            ],
            "section_identity_required_by_r01": False,
            "outcome_blind_role_freeze_required": True,
            "independent_external_lineage_required": True,
        },
    }


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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
        args.policy.resolve(),
    }
    if output in input_paths:
        raise ValueError("output collides with gate input")

    with args.policy.open(encoding="utf-8") as handle:
        policy = json.load(handle)
    result = evaluate_gate(
        _read_tsv(args.physical_units),
        _read_tsv(args.role_freeze),
        policy,
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
