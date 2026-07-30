import unittest

from scripts.r01_prepare_metadata_request import build_candidates, build_summary


def physical(
    study: str,
    patient: str,
    *,
    namespace: str = "ATLAS_TABLE_S2",
    block: str = "",
    grade: str = "E3_explicit",
    status: str = "EXCLUDED_MISSING_BLOCK",
) -> dict[str, str]:
    return {
        "source_namespace": namespace,
        "study_id": study,
        "patient_id": patient,
        "block_id": block,
        "evidence_grade": grade,
        "record_status": status,
    }


class PrepareMetadataRequestTests(unittest.TestCase):
    def test_accepted_specimen_equivalent_is_not_requested_again(self):
        rows = [
            {
                "source_namespace": "external_geo",
                "study_id": "GEO::GSE1",
                "patient_id": "GEO::GSE1::P1",
                "block_id": "",
                "physical_specimen_id": "SAMN1",
                "identity_granularity": "patient_linked_physical_specimen",
                "block_equivalent_status": "ACCEPTED_BLOCK_EQUIVALENT",
                "evidence_grade": "E2_corroborated",
                "record_status": "RESOLVED_INCLUDED_CANDIDATE",
            }
        ]

        self.assertEqual(build_candidates(rows, [], []), [])

    def test_known_official_conflict_is_retained_only_as_screened_exclusion(self):
        rows = [
            {
                "source_namespace": "ATLAS_TABLE_S2",
                "study_id": "atlas::conflict",
                "patient_id": "",
                "block_id": "",
                "evidence_grade": "EC_conflict",
                "record_status": "QUARANTINED_METADATA_CONFLICT",
            }
        ]
        atlas = [
            {
                "logical_unit_id": "atlas::conflict",
                "data_availability": "GSE242311",
                "references": "PMID: 36674951",
                "local_data_status": "LOCAL_DIRECTORY_PRESENT",
            }
        ]
        conflicts = {
            "atlas::conflict": {
                "accession": "GSE242311",
                "official_record_count": "16",
                "official_patient_count": "5",
                "official_reference": "PMID:39456890",
                "conflict_reason": "Atlas patient count and PMID conflict with GEO",
            }
        }

        result = build_candidates(
            rows,
            atlas,
            [],
            official_conflicts=conflicts,
        )

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["request_status"], "SCREENED_CONFLICT_EXCLUDED")
        self.assertEqual(
            result[0]["independence_status"],
            "OFFICIAL_METADATA_CONFLICT_QUARANTINED",
        )
        self.assertEqual(result[0]["distinct_patient_ids"], "5")
        self.assertEqual(result[0]["reference_or_title"], "PMID:39456890")

    def test_only_patient_known_block_missing_lineages_are_requested(self):
        rows = [
            physical("atlas1", f"p{i}") for i in range(4)
        ] + [
            physical("complete", "p", block="b"),
            physical("weak", "p", grade="E1_weak"),
            physical(
                "conflict",
                "p",
                grade="EC_conflict",
                status="QUARANTINED_METADATA_CONFLICT",
            ),
        ]
        atlas = [
            {
                "logical_unit_id": "atlas1",
                "data_availability": "GSE123",
                "references": "PMID: 123",
                "local_data_status": "LOCAL_DIRECTORY_PRESENT",
                "potential_tls_positive_count": "999",
            }
        ]

        result = build_candidates(rows, atlas, [])

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["logical_unit_id"], "atlas1")
        self.assertEqual(result[0]["priority_tier"], "P1_PUBLIC_HIGH_COVERAGE")
        self.assertEqual(result[0]["independence_status"], "UNRESOLVED_DO_NOT_COUNT")
        self.assertNotIn("tls", " ".join(result[0]).lower())

    def test_hest_locator_is_recovered_without_outcome_fields(self):
        link = "https://example.org/study/7"
        from scripts.r01_prepare_metadata_request import _hash

        study = f"HEST_STUDY::{_hash(link)}"
        rows = [
            physical(study, f"p{i}", namespace="HEST") for i in range(2)
        ]
        hest = [
            {
                "raw_study_link": link,
                "raw_dataset_title": "Identity-only title",
            }
        ]

        result = build_candidates(rows, [], hest)

        self.assertEqual(result[0]["public_locator"], link)
        self.assertEqual(result[0]["priority_tier"], "P2_PUBLIC_LOWER_COVERAGE")
        self.assertEqual(
            result[0]["request_status"], "PROPOSED_REQUIRES_APPROVAL"
        )

    def test_shared_pubmed_aliases_form_one_unresolved_provenance_group(self):
        from scripts.r01_prepare_metadata_request import _hash

        link = "https://pubmed.ncbi.nlm.nih.gov/35231421/"
        hest_study = f"HEST_STUDY::{_hash(link)}"
        rows = [
            physical("atlas1", "atlas-p"),
            physical(hest_study, "hest-p", namespace="HEST"),
        ]
        atlas = [
            {
                "logical_unit_id": "atlas1",
                "data_availability": "GSE175540",
                "references": "PMID:\u00a035231421",
                "local_data_status": "LOCAL_DIRECTORY_PRESENT",
            }
        ]
        hest = [{"raw_study_link": link, "raw_dataset_title": "same paper"}]

        result = build_candidates(rows, atlas, hest)

        self.assertEqual(len(result), 2)
        self.assertEqual(
            {row["provenance_group_id"] for row in result},
            {result[0]["provenance_group_id"]},
        )
        self.assertTrue(
            all(
                row["independence_status"]
                == "CROSS_AGGREGATOR_OVERLAP_UNRESOLVED_DO_NOT_COUNT"
                for row in result
            )
        )

    def test_known_current_lineage_is_not_treated_as_incremental(self):
        rows = [physical("atlas-htan", f"p{i}") for i in range(4)]
        atlas = [
            {
                "logical_unit_id": "atlas-htan",
                "data_availability": (
                    "https://data.humantumoratlas.org/publications/"
                    "vanderbilt_crc_chen_2021"
                ),
                "references": "PMID: 38065082",
                "local_data_status": "LOCAL_DIRECTORY_PRESENT",
            }
        ]

        result = build_candidates(rows, atlas, [])

        self.assertEqual(result[0]["priority_tier"], "P4_KNOWN_CURRENT_LINEAGE")
        self.assertEqual(result[0]["known_current_eligible_unit"], "HTAN_VANDERBILT_CRC")
        self.assertEqual(result[0]["request_status"], "NOT_AN_INCREMENTAL_UNIT")

    def test_summary_does_not_precount_unresolved_candidates(self):
        candidates = [
            {
                "priority_tier": "P1_PUBLIC_HIGH_COVERAGE",
            }
            for _ in range(12)
        ]

        summary = build_summary(
            candidates,
            {"status": "BLOCKED_IDENTITY", "claim_eligible_logical_units": 2},
        )

        self.assertEqual(summary["minimum_additional_independent_logical_units"], 4)
        self.assertEqual(summary["approval_status"], "NOT_APPROVED")
        self.assertIn("zero units", summary["independence_policy"])

    def test_summary_can_record_the_user_approved_metadata_only_scope(self):
        summary = build_summary(
            [],
            {"claim_eligible_logical_units": 6, "status": "COMPLETE_WITH_EXCLUSIONS"},
            approval_status="APPROVED_METADATA_ONLY",
        )

        self.assertEqual(
            summary["approval_status"],
            "APPROVED_METADATA_ONLY",
        )
        self.assertEqual(summary["minimum_additional_independent_logical_units"], 0)


if __name__ == "__main__":
    unittest.main()
