import tempfile
import unittest
from pathlib import Path

from scripts.r01_extract_tenx_explicit import extract_explicit_tenx_units


class ExtractTenxExplicitTests(unittest.TestCase):
    def test_bundled_block_section_and_atlas_patient_mapping_are_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "sources.tsv"
            manifest.write_text(
                "source\tdataset_id\tdescription\turl\tdownloaded\tnotes\n"
                "10x Genomics\tV1_Breast_Cancer_Block_A_Section_1\t"
                "human breast cancer, Block A Section 1\thttps://x/1\tyes\t"
                "PAPER SAMPLE (Table S2): BRCA B37 = this dataset.\n"
                "10x Genomics\tV1_Breast_Cancer_Block_A_Section_2\t"
                "human breast cancer, Block A Section 2\thttps://x/2\tyes\t"
                "Section 2 of same block; NOT used by the paper.\n"
                "10x Genomics\tOther\tno block metadata\thttps://x/3\tyes\t-\n",
                encoding="utf-8",
            )
            atlas = root / "atlas.tsv"
            atlas.write_text(
                "record_id\tcancer_raw\tpatient_id_raw\tsample_id_raw\n"
                "B37\tBRCA\tB_P9\tB37\n",
                encoding="utf-8",
            )
            local = root / "tenx"
            (local / "V1_Breast_Cancer_Block_A_Section_1").mkdir(parents=True)
            (local / "V1_Breast_Cancer_Block_A_Section_2").mkdir(parents=True)

            rows = extract_explicit_tenx_units(manifest, atlas, local)

            self.assertEqual(len(rows), 2)
            self.assertEqual({row["raw_block_id"] for row in rows}, {"Block A"})
            self.assertEqual(
                {row["raw_section_id"] for row in rows},
                {"Section 1", "Section 2"},
            )
            self.assertEqual({row["raw_patient_id"] for row in rows}, {"B_P9"})
            self.assertEqual(
                {row["evidence_grade"] for row in rows},
                {"E2_corroborated"},
            )
            self.assertEqual({row["record_status"] for row in rows}, {"MATCHED"})

    def test_no_patient_mapping_keeps_record_out_of_claim_eligible_state(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "sources.tsv"
            manifest.write_text(
                "source\tdataset_id\tdescription\turl\tdownloaded\tnotes\n"
                "10x Genomics\tX_Block_B_Section_1\t"
                "Block B Section 1\thttps://x/1\tyes\tno paper mapping\n",
                encoding="utf-8",
            )
            atlas = root / "atlas.tsv"
            atlas.write_text(
                "record_id\tcancer_raw\tpatient_id_raw\tsample_id_raw\n",
                encoding="utf-8",
            )
            local = root / "tenx"
            (local / "X_Block_B_Section_1").mkdir(parents=True)

            row = extract_explicit_tenx_units(manifest, atlas, local)[0]

            self.assertEqual(row["raw_patient_id"], "")
            self.assertEqual(row["evidence_grade"], "E1_weak")
            self.assertEqual(row["record_status"], "PATIENT_UNKNOWN")


if __name__ == "__main__":
    unittest.main()
