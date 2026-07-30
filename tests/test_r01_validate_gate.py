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


class ValidateR01GateTests(unittest.TestCase):
    def test_in_progress_while_metadata_audit_is_not_exhausted(self):
        result = evaluate_gate(
            [candidate("study1")],
            [],
            {"metadata_audit_exhausted": False, "role_freeze_attempted": False},
        )
        self.assertEqual(result["status"], "IN_PROGRESS")
        self.assertEqual(result["claim_eligible_logical_units"], 1)

    def test_hard_blocks_identity_after_exhaustion_below_six_units(self):
        result = evaluate_gate(
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

        result = evaluate_gate(
            physical,
            [],
            {"metadata_audit_exhausted": True, "role_freeze_attempted": False},
        )

        self.assertEqual(result["status"], "BLOCKED_IDENTITY")
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
        result = evaluate_gate(
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
        result = evaluate_gate(
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

        result = evaluate_gate(
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

        result = evaluate_gate(
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

        result = evaluate_gate(
            physical,
            frozen,
            {"metadata_audit_exhausted": True, "role_freeze_attempted": True},
        )

        self.assertEqual(result["status"], "BLOCKED_INDEPENDENCE")
        self.assertIn("invalid_primary_role:study0", result["blockers"])
        self.assertIn("missing_leakage_group:study0", result["blockers"])


if __name__ == "__main__":
    unittest.main()
