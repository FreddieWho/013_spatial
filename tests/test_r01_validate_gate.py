import unittest

from scripts.r01_validate_gate import evaluate_gate


def candidate(study: str, suffix: str = "1") -> dict[str, str]:
    return {
        "physical_unit_id": f"{study}::{suffix}",
        "study_id": study,
        "patient_id": f"{study}::P",
        "block_id": f"{study}::B",
        "section_id": "",
        "evidence_grade": "E2_corroborated",
        "record_status": "RESOLVED_INCLUDED_CANDIDATE",
    }


def specimen_candidate(study: str, suffix: str = "1") -> dict[str, str]:
    row = candidate(study, suffix)
    row.update(
        {
            "block_id": "",
            "physical_specimen_id": f"{study}::SPECIMEN::{suffix}",
            "identity_granularity": "patient_linked_physical_specimen",
            "block_equivalent_status": "ACCEPTED_BLOCK_EQUIVALENT",
            "block_equivalent_basis": (
                "Official patient mapping, physical tissue description, and "
                "stable specimen locator."
            ),
            "identity_status": "PATIENT_SPECIMEN_CORROBORATED",
        }
    )
    return row


def specimen_support(
    row: dict[str, str],
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    accession = row["study_id"].removeprefix("GEO::")
    gsm_asset = f"asset::{accession}::gsm"
    series_asset = f"asset::{accession}::series"
    evidence = [
        {
            "entity_id": row["physical_unit_id"],
            "field": "patient_id",
            "source_asset": gsm_asset,
            "metadata_path": f"metadata/{accession}_gsm.soft",
            "metadata_key": "patient",
            "normalized_value": row["patient_id"],
            "raw_value": "patient raw",
            "evidence_grade": "E2_corroborated",
            "conflict_flag": "false",
        },
        {
            "entity_id": row["physical_unit_id"],
            "field": "physical_specimen_id",
            "source_asset": gsm_asset,
            "metadata_path": f"metadata/{accession}_gsm.soft",
            "metadata_key": "BioSample",
            "normalized_value": row["physical_specimen_id"],
            "raw_value": row["physical_specimen_id"],
            "evidence_grade": "E2_corroborated",
            "conflict_flag": "false",
        },
        {
            "entity_id": row["physical_unit_id"],
            "field": "block_id",
            "source_asset": series_asset,
            "metadata_path": f"metadata/{accession}_series.soft",
            "metadata_key": "series summary",
            "normalized_value": "",
            "raw_value": row["block_equivalent_basis"],
            "evidence_grade": "E2_corroborated",
            "conflict_flag": "false",
        },
    ]
    assets = [
        {
            "asset_id": asset_id,
            "accession": accession,
            "path": (
                f"metadata/{accession}_gsm.soft"
                if asset_id == gsm_asset
                else f"metadata/{accession}_series.soft"
            ),
            "record_status": "active",
            "checksum_status": "sha256_verified",
            "checksum": "a" * 64,
        }
        for asset_id in (gsm_asset, series_asset)
    ]
    duplicates = [
        {
            "duplicate_group_id": f"same_source::{accession}",
            "member_type": "study",
            "member_id": row["study_id"],
            "relation_type": "same_accession_same_source_lineage",
            "resolution": "canonical_geo_study; do_not_count_atlas_separately",
            "leakage_group_id": f"lineage::{accession}",
        },
        {
            "duplicate_group_id": f"same_source::{accession}",
            "member_type": "study",
            "member_id": f"atlas::{accession}",
            "relation_type": "same_accession_same_source_lineage",
            "resolution": "canonical_geo_study; do_not_count_atlas_separately",
            "leakage_group_id": f"lineage::{accession}",
        },
    ]
    return evidence, assets, duplicates


def explicit_support(
    row: dict[str, str],
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    study = row["study_id"]
    asset_id = f"asset::{study}::metadata"
    path = f"metadata/{study}.tsv"
    evidence = [
        {
            "entity_id": row["physical_unit_id"],
            "field": field,
            "source_asset": asset_id,
            "metadata_path": path,
            "metadata_key": field,
            "raw_value": f"raw::{field}",
            "normalized_value": row[field],
            "evidence_grade": row["evidence_grade"],
            "conflict_flag": "false",
        }
        for field in ("patient_id", "block_id")
    ]
    assets = [
        {
            "asset_id": asset_id,
            "path": path,
            "record_status": "PRESENT",
            "checksum_status": "sha256_verified",
            "checksum": "b" * 64,
        }
    ]
    suffix = study.removeprefix("study")
    leakage = f"leakage{suffix}" if suffix.isdigit() else f"lineage::{study}"
    duplicates = [
        {
            "duplicate_group_id": f"source_lineage::{study}",
            "member_type": "study",
            "member_id": study,
            "relation_type": "canonical_source_lineage",
            "resolution": "canonical_study",
            "leakage_group_id": leakage,
        }
    ]
    return evidence, assets, duplicates


def evaluate_supported(
    physical: list[dict[str, str]],
    roles: list[dict[str, str]],
    policy: dict[str, bool],
) -> dict[str, object]:
    evidence, assets, duplicates = [], [], []
    for row in physical:
        supported = (
            specimen_support(row)
            if row.get("block_equivalent_status") == "ACCEPTED_BLOCK_EQUIVALENT"
            else explicit_support(row)
        )
        evidence.extend(supported[0])
        assets.extend(supported[1])
        duplicates.extend(supported[2])
    return evaluate_gate(
        physical,
        roles,
        policy,
        identity_evidence=evidence,
        source_assets=assets,
        duplicate_groups=duplicates,
    )


class ValidateR01GateTests(unittest.TestCase):
    def test_in_progress_while_metadata_audit_is_not_exhausted(self):
        result = evaluate_supported(
            [candidate("study1")],
            [],
            {"metadata_audit_exhausted": False, "role_freeze_attempted": False},
        )
        self.assertEqual(result["status"], "IN_PROGRESS")
        self.assertEqual(result["claim_eligible_logical_units"], 1)

    def test_hard_blocks_identity_after_exhaustion_below_six_units(self):
        result = evaluate_supported(
            [candidate("study1")],
            [],
            {"metadata_audit_exhausted": True, "role_freeze_attempted": False},
        )
        self.assertEqual(result["status"], "BLOCKED_IDENTITY")
        self.assertIn("fewer_than_6_identity_complete_units", result["blockers"])

    def test_weak_identity_evidence_is_not_claim_eligible(self):
        physical = [candidate(f"study{i}") for i in range(6)]
        for row in physical:
            row["evidence_grade"] = "E1_weak"

        result = evaluate_supported(
            physical,
            [],
            {"metadata_audit_exhausted": True, "role_freeze_attempted": False},
        )

        self.assertEqual(result["status"], "BLOCKED_IDENTITY")
        self.assertEqual(result["claim_eligible_logical_units"], 0)

    def test_accepted_patient_linked_specimen_is_eligible_without_fake_block(self):
        physical = [candidate("explicit0"), candidate("explicit1")]
        specimen_rows = [
            specimen_candidate(f"GEO::GSE{i}") for i in range(4)
        ]
        physical.extend(specimen_rows)
        evidence, assets, duplicates = [], [], []
        for row in physical:
            supported = (
                specimen_support(row)
                if row in specimen_rows
                else explicit_support(row)
            )
            evidence.extend(supported[0])
            assets.extend(supported[1])
            duplicates.extend(supported[2])

        result = evaluate_gate(
            physical,
            [],
            {"metadata_audit_exhausted": True, "role_freeze_attempted": False},
            identity_evidence=evidence,
            source_assets=assets,
            duplicate_groups=duplicates,
        )

        self.assertEqual(result["claim_eligible_logical_units"], 6)
        self.assertEqual(
            result["completion_requirements"]["accepted_identity_granularities"],
            ["patient_block", "patient_linked_physical_specimen"],
        )
        self.assertTrue(all(not row["block_id"] for row in physical[2:]))

    def test_self_asserted_specimen_equivalent_without_evidence_cannot_green(self):
        physical = [candidate("explicit0"), candidate("explicit1")]
        physical.extend(
            specimen_candidate(f"GEO::GSE{i}") for i in range(4)
        )
        evidence, assets, duplicates = [], [], []
        for row in physical[:2]:
            supported = explicit_support(row)
            evidence.extend(supported[0])
            assets.extend(supported[1])
            duplicates.extend(supported[2])

        result = evaluate_gate(
            physical,
            [],
            {"metadata_audit_exhausted": True, "role_freeze_attempted": False},
            identity_evidence=evidence,
            source_assets=assets,
            duplicate_groups=duplicates,
        )

        self.assertEqual(result["claim_eligible_logical_units"], 2)
        self.assertTrue(result["identity_integrity_errors"])

    def test_self_asserted_explicit_patient_block_without_evidence_cannot_green(self):
        physical = [candidate(f"FORGED_EXPLICIT::{i}") for i in range(6)]
        frozen = [
            {
                "logical_unit_id": row["study_id"],
                "leakage_group_id": f"forged::{i}",
                "primary_role": (
                    "external_validation" if i == 5 else "training"
                ),
                "record_status": "FROZEN",
            }
            for i, row in enumerate(physical)
        ]

        result = evaluate_gate(
            physical,
            frozen,
            {"metadata_audit_exhausted": True, "role_freeze_attempted": True},
        )

        self.assertEqual(result["status"], "BLOCKED_IDENTITY")
        self.assertEqual(result["claim_eligible_logical_units"], 0)
        self.assertTrue(result["identity_integrity_errors"])

    def test_specimen_evidence_and_same_source_lineage_are_referentially_checked(self):
        row = specimen_candidate("GEO::GSE1")
        evidence, assets, duplicates = specimen_support(row)
        evidence[0]["normalized_value"] = "GEO::GSE1::WRONG"

        result = evaluate_gate(
            [row],
            [],
            {"metadata_audit_exhausted": True, "role_freeze_attempted": False},
            identity_evidence=evidence,
            source_assets=assets,
            duplicate_groups=duplicates,
        )
        self.assertEqual(result["claim_eligible_logical_units"], 0)
        self.assertTrue(
            any("patient_id_mismatch" in error for error in result["identity_integrity_errors"])
        )

        evidence, assets, _ = specimen_support(row)
        result = evaluate_gate(
            [row],
            [],
            {"metadata_audit_exhausted": True, "role_freeze_attempted": False},
            identity_evidence=evidence,
            source_assets=assets,
            duplicate_groups=[],
        )
        self.assertEqual(result["claim_eligible_logical_units"], 0)
        self.assertTrue(
            any("missing_canonical_lineage" in error for error in result["identity_integrity_errors"])
        )

    def test_specimen_evidence_requires_verified_active_source_assets(self):
        row = specimen_candidate("GEO::GSE1")
        evidence, assets, duplicates = specimen_support(row)
        assets[0]["checksum_status"] = "not_verified"

        result = evaluate_gate(
            [row],
            [],
            {"metadata_audit_exhausted": True, "role_freeze_attempted": False},
            identity_evidence=evidence,
            source_assets=assets,
            duplicate_groups=duplicates,
        )

        self.assertEqual(result["claim_eligible_logical_units"], 0)
        self.assertTrue(
            any("unverified_source_asset" in error for error in result["identity_integrity_errors"])
        )

    def test_specimen_equivalent_fails_closed_if_any_required_field_is_missing(self):
        required = [
            "physical_specimen_id",
            "identity_granularity",
            "block_equivalent_status",
            "block_equivalent_basis",
            "identity_status",
        ]
        for field in required:
            with self.subTest(field=field):
                row = specimen_candidate("specimen")
                row[field] = ""
                result = evaluate_gate(
                    [row],
                    [],
                    {
                        "metadata_audit_exhausted": True,
                        "role_freeze_attempted": False,
                    },
                )
                self.assertEqual(result["claim_eligible_logical_units"], 0)

    def test_specimen_equivalent_cannot_also_smuggle_a_locator_into_block(self):
        row = specimen_candidate("specimen")
        row["block_id"] = "SAMN123"

        result = evaluate_gate(
            [row],
            [],
            {"metadata_audit_exhausted": True, "role_freeze_attempted": False},
        )

        self.assertEqual(result["claim_eligible_logical_units"], 0)

    def test_complete_requires_six_frozen_units_and_external_lineage(self):
        physical = [candidate(f"study{i}") for i in range(6)]
        frozen = [
            {
                "logical_unit_id": f"study{i}",
                "leakage_group_id": f"leakage{i}",
                "primary_role": "external_validation" if i == 5 else "training",
                "record_status": "FROZEN",
            }
            for i in range(6)
        ]
        result = evaluate_supported(
            physical,
            frozen,
            {"metadata_audit_exhausted": True, "role_freeze_attempted": True},
        )
        self.assertEqual(result["status"], "COMPLETE")

    def test_independence_block_requires_attempted_freeze(self):
        physical = [candidate(f"study{i}") for i in range(6)]
        frozen = [
            {
                "logical_unit_id": f"study{i}",
                "leakage_group_id": "same",
                "primary_role": "training",
                "record_status": "FROZEN",
            }
            for i in range(6)
        ]
        result = evaluate_supported(
            physical,
            frozen,
            {"metadata_audit_exhausted": True, "role_freeze_attempted": True},
        )
        self.assertEqual(result["status"], "BLOCKED_INDEPENDENCE")

    def test_freeze_cannot_reference_an_ineligible_logical_unit(self):
        physical = [candidate(f"study{i}") for i in range(6)]
        frozen = [
            {
                "logical_unit_id": "unregistered" if i == 5 else f"study{i}",
                "leakage_group_id": f"leakage{i}",
                "primary_role": "external_validation" if i == 5 else "training",
                "record_status": "FROZEN",
            }
            for i in range(6)
        ]

        result = evaluate_supported(
            physical,
            frozen,
            {"metadata_audit_exhausted": True, "role_freeze_attempted": True},
        )

        self.assertEqual(result["status"], "BLOCKED_INDEPENDENCE")
        self.assertIn("frozen_unit_not_claim_eligible:unregistered", result["blockers"])

    def test_freeze_rejects_duplicate_roles_and_leakage_role_conflict(self):
        physical = [candidate(f"study{i}") for i in range(6)]
        frozen = [
            {
                "logical_unit_id": f"study{i}",
                "leakage_group_id": f"leakage{i}",
                "primary_role": "external_validation" if i == 5 else "training",
                "record_status": "FROZEN",
            }
            for i in range(6)
        ]
        frozen.append(
            {
                "logical_unit_id": "study0",
                "leakage_group_id": "leakage0",
                "primary_role": "external_validation",
                "record_status": "FROZEN",
            }
        )

        result = evaluate_supported(
            physical,
            frozen,
            {"metadata_audit_exhausted": True, "role_freeze_attempted": True},
        )

        self.assertEqual(result["status"], "BLOCKED_INDEPENDENCE")
        self.assertIn("multiple_frozen_roles:study0", result["blockers"])
        self.assertIn("leakage_group_role_conflict:leakage0", result["blockers"])

    def test_freeze_requires_a_valid_primary_role_and_leakage_group(self):
        physical = [candidate(f"study{i}") for i in range(6)]
        frozen = [
            {
                "logical_unit_id": f"study{i}",
                "leakage_group_id": "" if i == 0 else f"leakage{i}",
                "primary_role": "" if i == 0 else (
                    "external_validation" if i == 5 else "training"
                ),
                "record_status": "FROZEN",
            }
            for i in range(6)
        ]

        result = evaluate_supported(
            physical,
            frozen,
            {"metadata_audit_exhausted": True, "role_freeze_attempted": True},
        )

        self.assertEqual(result["status"], "BLOCKED_INDEPENDENCE")
        self.assertIn("invalid_primary_role:study0", result["blockers"])
        self.assertIn("missing_leakage_group:study0", result["blockers"])

    def test_external_lineage_must_be_independent_of_all_claim_bearing_roles(self):
        physical = [candidate(f"study{i}") for i in range(6)]
        frozen = [
            {
                "logical_unit_id": f"study{i}",
                "leakage_group_id": "same" if i in {0, 5} else f"leakage{i}",
                "primary_role": (
                    "serial_section_validation"
                    if i == 0
                    else ("external_validation" if i == 5 else "training")
                ),
                "record_status": "FROZEN",
            }
            for i in range(6)
        ]

        result = evaluate_supported(
            physical,
            frozen,
            {"metadata_audit_exhausted": True, "role_freeze_attempted": True},
        )

        self.assertEqual(result["status"], "BLOCKED_INDEPENDENCE")
        self.assertIn("leakage_group_role_conflict:same", result["blockers"])


if __name__ == "__main__":
    unittest.main()
