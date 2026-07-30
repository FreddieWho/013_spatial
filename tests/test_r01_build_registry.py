import unittest

from scripts.r01_build_registry import (
    build_duplicate_groups,
    physical_units_from_atlas,
    physical_units_from_hest,
    physical_units_from_htan,
    physical_units_from_tenx,
)


class BuildRegistryTests(unittest.TestCase):
    def test_atlas_patient_is_namespaced_but_block_is_not_invented(self):
        rows = [
            {
                "record_id": "A",
                "study_namespace": "study::1",
                "patient_id_raw": "P1",
                "sample_id_raw": "S1",
                "platform_raw": "Visium",
            },
            {
                "record_id": "B",
                "study_namespace": "study::1",
                "patient_id_raw": "P1\nP1 & P2",
                "sample_id_raw": "S2",
                "platform_raw": "Visium",
            },
        ]

        units = physical_units_from_atlas(rows)
        by_record = {row["source_record_id"]: row for row in units}

        self.assertTrue(by_record["A"]["patient_id"].startswith("study::1::PATIENT::"))
        self.assertEqual(by_record["A"]["block_id"], "")
        self.assertEqual(by_record["A"]["section_id"], "")
        self.assertEqual(by_record["A"]["record_status"], "EXCLUDED_MISSING_BLOCK")
        self.assertEqual(by_record["B"]["patient_id"], "")
        self.assertEqual(
            by_record["B"]["record_status"],
            "QUARANTINED_METADATA_CONFLICT",
        )

    def test_hest_conflict_and_unknown_block_remain_excluded(self):
        rows = [
            {
                "source_record_id": "TENX1",
                "raw_patient": "P1",
                "raw_study_link": "https://study/1",
                "raw_slide_id": "",
                "raw_st_technology": "Visium",
                "evidence_grade": "E3_explicit",
                "identity_conflict": "no",
            },
            {
                "source_record_id": "TENX2",
                "raw_patient": "P1",
                "raw_study_link": "https://study/1",
                "raw_slide_id": "",
                "raw_st_technology": "Xenium",
                "evidence_grade": "EC_conflict",
                "identity_conflict": "yes",
            },
        ]
        units = physical_units_from_hest(rows)
        self.assertEqual(units[0]["record_status"], "EXCLUDED_MISSING_BLOCK")
        self.assertEqual(
            units[1]["record_status"],
            "QUARANTINED_METADATA_CONFLICT",
        )

    def test_htan_matched_identity_is_candidate_but_local_only_is_excluded(self):
        rows = [
            {
                "source_record_id": "HTAN::a.h5ad",
                "asset_name": "a.h5ad",
                "raw_sample_key": "S1",
                "raw_patient_id": "P1",
                "raw_block_id": "B1",
                "evidence_grade": "E3_explicit",
                "record_status": "MATCHED",
            },
            {
                "source_record_id": "HTAN::b.h5ad",
                "asset_name": "b.h5ad",
                "raw_sample_key": "",
                "raw_patient_id": "",
                "raw_block_id": "",
                "evidence_grade": "E0_unknown",
                "record_status": "LOCAL_ONLY",
            },
        ]
        units = physical_units_from_htan(rows)
        self.assertEqual(units[0]["record_status"], "RESOLVED_INCLUDED_CANDIDATE")
        self.assertTrue(units[0]["block_id"])
        self.assertEqual(units[0]["section_id"], "")
        self.assertEqual(
            units[1]["record_status"],
            "EXCLUDED_MISSING_PATIENT_OR_BLOCK",
        )

    def test_htan_sample_key_is_not_promoted_to_section_identity(self):
        rows = [
            {
                "source_record_id": "HTAN::a.h5ad",
                "asset_name": "a.h5ad",
                "raw_sample_key": "SAME",
                "raw_patient_id": "P1",
                "raw_block_id": "B1",
                "evidence_grade": "E3_explicit",
                "record_status": "MATCHED",
            },
            {
                "source_record_id": "HTAN::b.h5ad",
                "asset_name": "b.h5ad",
                "raw_sample_key": "SAME",
                "raw_patient_id": "P2",
                "raw_block_id": "B2",
                "evidence_grade": "E3_explicit",
                "record_status": "MATCHED",
            },
        ]
        units = physical_units_from_htan(rows)
        self.assertEqual({unit["section_id"] for unit in units}, {""})
        self.assertNotEqual(units[0]["physical_unit_id"], units[1]["physical_unit_id"])

    def test_explicit_tenx_block_and_section_form_one_candidate_lineage(self):
        rows = [
            {
                "source_record_id": "TENX::section1",
                "dataset_id": "V1_Breast_Cancer_Block_A_Section_1",
                "raw_patient_id": "B_P9",
                "raw_block_id": "Block A",
                "raw_section_id": "Section 1",
                "evidence_grade": "E2_corroborated",
                "record_status": "MATCHED",
            },
            {
                "source_record_id": "TENX::section2",
                "dataset_id": "V1_Breast_Cancer_Block_A_Section_2",
                "raw_patient_id": "B_P9",
                "raw_block_id": "Block A",
                "raw_section_id": "Section 2",
                "evidence_grade": "E2_corroborated",
                "record_status": "MATCHED",
            },
        ]
        units = physical_units_from_tenx(rows)
        self.assertEqual({unit["study_id"] for unit in units}, {
            "TENX_V1_BREAST_CANCER_BLOCK_A"
        })
        self.assertEqual(
            {unit["record_status"] for unit in units},
            {"RESOLVED_INCLUDED_CANDIDATE"},
        )
        self.assertEqual(len({unit["block_id"] for unit in units}), 1)
        self.assertEqual(len({unit["section_id"] for unit in units}), 2)

    def test_duplicate_groups_include_cross_source_and_htan_physical_groups(self):
        cross = [
            {
                "candidate_group_id": "URL::1",
                "left_record_id": "HEST::A",
                "right_record_id": "TENX::A",
                "relation_type": "same_source_page",
                "evidence_grade": "E2_corroborated",
                "resolution": "CONFIRMED_SOURCE_DATASET",
                "leakage_group_id": "LEAKAGE::1",
            },
            {
                "candidate_group_id": "URL::1",
                "left_record_id": "HEST::B",
                "right_record_id": "TENX::A",
                "relation_type": "same_source_page_only",
                "evidence_grade": "E0_unknown",
                "resolution": "NO_EVIDENCE",
                "leakage_group_id": "",
            },
        ]
        htan = [
            {
                "source_record_id": "HTAN::a",
                "raw_sample_key": "S1",
                "raw_patient_id": "P1",
                "raw_block_id": "B1",
                "evidence_grade": "E3_explicit",
                "record_status": "MATCHED",
            },
            {
                "source_record_id": "HTAN::b",
                "raw_sample_key": "S1",
                "raw_patient_id": "P1",
                "raw_block_id": "B1",
                "evidence_grade": "E3_explicit",
                "record_status": "MATCHED",
            },
        ]

        groups = build_duplicate_groups(cross, htan)

        cross_members = [
            row for row in groups if row["leakage_group_id"] == "LEAKAGE::1"
        ]
        self.assertEqual(
            {row["member_id"] for row in cross_members},
            {"HEST::A", "TENX::A"},
        )
        self.assertNotIn("HEST::B", {row["member_id"] for row in groups})
        self.assertTrue(
            any(
                row["relation_type"] == "same_sample_alias_within_block"
                and row["resolution"] == "UNRESOLVED_CONSERVATIVE_GROUP"
                for row in groups
            )
        )
        self.assertTrue(
            any(row["relation_type"] == "same_patient_other_block" for row in groups)
        )


if __name__ == "__main__":
    unittest.main()
