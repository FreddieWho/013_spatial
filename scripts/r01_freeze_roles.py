#!/usr/bin/env python3
"""Freeze the outcome-blind R-01 roles after identity eligibility is complete."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from datetime import date
from pathlib import Path
from typing import Iterable, Mapping

try:
    from scripts.r01_build_registry import ROLE_FREEZE_FIELDS
    from scripts.r01_validate_gate import _claim_identity_complete
except ModuleNotFoundError:  # Direct execution via ``python scripts/...``.
    from r01_build_registry import ROLE_FREEZE_FIELDS
    from r01_validate_gate import _claim_identity_complete


ROLE_PLAN = {
    "GEO::GSE274557": {
        "leakage_group_id": "lineage::PRJNA1147209",
        "primary_role": "discovery",
        "special_capability": "multi-site PDAC specimen cohort",
        "rationale": "Largest patient-linked physical-specimen lineage in the approved batch.",
        "allowed_use": "discovery and training-fold method development",
        "forbidden_use": "external confirmation or role reassignment after analysis",
    },
    "HTAN_VANDERBILT_CRC": {
        "leakage_group_id": "LINEAGE::HTAN_VANDERBILT_CRC",
        "primary_role": "training",
        "special_capability": "explicit patient-block crosswalk",
        "rationale": "Largest existing lineage with explicit patient and block identity.",
        "allowed_use": "training and training-fold tuning",
        "forbidden_use": "external confirmation",
    },
    "GEO::GSE274103": {
        "leakage_group_id": "lineage::PRJNA1144961",
        "primary_role": "internal_validation",
        "special_capability": "FFPE PDAC physical-specimen cohort",
        "rationale": "Independent BioProject with five patient-linked FFPE tissues.",
        "allowed_use": "predeclared internal validation",
        "forbidden_use": "external confirmation or discovery-driven replacement",
    },
    "GEO::GSE226997": {
        "leakage_group_id": "lineage::PRJNA942633",
        "primary_role": "internal_validation",
        "special_capability": "four-patient CRC physical-specimen cohort",
        "rationale": "Independent BioProject retained under the strictest E2 fail-closed rule.",
        "allowed_use": "predeclared internal validation",
        "forbidden_use": "external confirmation or identity relaxation",
    },
    "GEO::GSE211956": {
        "leakage_group_id": "lineage::PRJNA873051",
        "primary_role": "external_validation",
        "special_capability": "supplement-backed patient-to-Visium crosswalk",
        "rationale": "Independent ovarian lineage with an official patient-to-GSM crosswalk.",
        "allowed_use": "locked external validation only",
        "forbidden_use": "selection, tuning, threshold choice, or role replacement",
    },
    "TENX_V1_BREAST_CANCER_BLOCK_A": {
        "leakage_group_id": "LINEAGE::TENX_V1_BREAST_CANCER_BLOCK_A",
        "primary_role": "serial_section_validation",
        "special_capability": "two explicit sections from one patient block",
        "rationale": "Only current lineage with explicit same-block section identities.",
        "allowed_use": "serial-section capability validation",
        "forbidden_use": "independent patient-level external confirmation",
    },
}


def freeze_roles(
    physical_units: Iterable[Mapping[str, str]],
    *,
    freeze_date: str,
) -> list[dict[str, str]]:
    eligible = {
        row.get("study_id", "")
        for row in physical_units
        if row.get("record_status") == "RESOLVED_INCLUDED_CANDIDATE"
        and row.get("patient_id")
        and _claim_identity_complete(row)
    }
    missing = sorted(set(ROLE_PLAN) - eligible)
    if missing:
        raise ValueError(
            "planned role unit is not eligible: " + ", ".join(missing)
        )

    return [
        {
            "logical_unit_id": study,
            "leakage_group_id": plan["leakage_group_id"],
            "primary_role": plan["primary_role"],
            "special_capability": plan["special_capability"],
            "freeze_version": "R01-v1",
            "freeze_date": freeze_date,
            "rationale": plan["rationale"],
            "allowed_use": plan["allowed_use"],
            "forbidden_use": plan["forbidden_use"],
            "record_status": "FROZEN",
        }
        for study, plan in ROLE_PLAN.items()
    ]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(rows: list[dict[str, str]], output: Path) -> None:
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
            writer = csv.DictWriter(
                handle,
                fieldnames=ROLE_FREEZE_FIELDS,
                delimiter="\t",
            )
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
        "--physical-units",
        type=Path,
        default=root / "infra/sample-registry/physical_units.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "infra/sample-registry/role_freeze.tsv",
    )
    parser.add_argument("--freeze-date", default=date.today().isoformat())
    args = parser.parse_args()

    output = args.output.resolve()
    try:
        output.relative_to(root)
    except ValueError as exc:
        raise ValueError("output escapes project root") from exc
    if output == args.physical_units.resolve():
        raise ValueError("output collides with physical-unit input")

    rows = freeze_roles(
        _read_tsv(args.physical_units),
        freeze_date=args.freeze_date,
    )
    _write_tsv(rows, output)
    print(f"R-01 role_freeze={len(rows)}")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
