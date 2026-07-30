import unittest

from scripts.r01_freeze_roles import ROLE_PLAN, freeze_roles


def eligible(study: str) -> dict[str, str]:
    return {
        "physical_unit_id": f"{study}::1",
        "study_id": study,
        "patient_id": f"{study}::PATIENT::1",
        "block_id": f"{study}::BLOCK::1",
        "evidence_grade": "E2_corroborated",
        "record_status": "RESOLVED_INCLUDED_CANDIDATE",
    }


class FreezeRolesTests(unittest.TestCase):
    def test_freeze_is_outcome_blind_and_covers_exactly_six_units(self):
        rows = [eligible(study) for study in ROLE_PLAN]

        frozen = freeze_roles(rows, freeze_date="2026-07-31")

        self.assertEqual(len(frozen), 6)
        self.assertEqual(
            {row["logical_unit_id"] for row in frozen},
            set(ROLE_PLAN),
        )
        self.assertEqual(
            [
                row["logical_unit_id"]
                for row in frozen
                if row["primary_role"] == "external_validation"
            ],
            ["GEO::GSE211956"],
        )
        forbidden_selection_terms = {"tls", "crs", "model", "performance", "outcome"}
        selection_text = " ".join(
            f"{row['rationale']} {row['allowed_use']} {row['forbidden_use']}"
            for row in frozen
        ).lower()
        self.assertTrue(
            forbidden_selection_terms.isdisjoint(selection_text.split())
        )

    def test_ineligible_planned_unit_fails_closed(self):
        rows = [eligible(study) for study in ROLE_PLAN]
        rows[-1]["record_status"] = "EXCLUDED_MISSING_PATIENT_OR_BLOCK"

        with self.assertRaisesRegex(ValueError, "planned role unit is not eligible"):
            freeze_roles(rows, freeze_date="2026-07-31")


if __name__ == "__main__":
    unittest.main()
