import csv
import tempfile
import unittest
from pathlib import Path

from scripts.r01_extract_htan import extract_htan_rows, write_rows


class ExtractHtanTests(unittest.TestCase):
    def test_union_preserves_local_only_and_metadata_only_assets(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = root / "meta.csv"
            asset_dir = root / "assets"
            asset_dir.mkdir()
            metadata.write_text(
                "duplicated,sample_key,block_name,patient_name,project,trimmed_adata\n"
                ",S1,B1,P1,HTAN,data/ST/a.h5ad\n"
                "Yes,S1,B1,P1,HTAN,data/ST/b.h5ad\n"
                ",S2,B2,P2,HTAN,data/ST/missing.h5ad\n",
                encoding="utf-8",
            )
            (asset_dir / "a.h5ad").touch()
            (asset_dir / "b.h5ad").touch()
            (asset_dir / "local_only.h5ad").touch()

            rows = extract_htan_rows(metadata, asset_dir)

            self.assertEqual([row["asset_name"] for row in rows], [
                "a.h5ad",
                "b.h5ad",
                "local_only.h5ad",
                "missing.h5ad",
            ])
            by_name = {row["asset_name"]: row for row in rows}
            self.assertEqual(by_name["local_only.h5ad"]["record_status"], "LOCAL_ONLY")
            self.assertEqual(by_name["local_only.h5ad"]["evidence_grade"], "E0_unknown")
            self.assertEqual(by_name["missing.h5ad"]["record_status"], "METADATA_ONLY")
            self.assertEqual(by_name["b.h5ad"]["duplicate_flag"], "Yes")
            self.assertEqual(by_name["a.h5ad"]["raw_patient_id"], "P1")

    def test_duplicate_asset_names_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = root / "meta.csv"
            asset_dir = root / "assets"
            asset_dir.mkdir()
            metadata.write_text(
                "duplicated,sample_key,block_name,patient_name,project,trimmed_adata\n"
                ",S1,B1,P1,HTAN,data/ST/a.h5ad\n"
                ",S2,B2,P2,HTAN,data/ST/a.h5ad\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "duplicate metadata asset"):
                extract_htan_rows(metadata, asset_dir)

    def test_written_output_is_deterministic(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / "out.tsv"
            rows = [
                {
                    "source_record_id": "x",
                    "study_id": "S",
                    "asset_name": "a.h5ad",
                    "local_present": "yes",
                    "raw_sample_key": "A",
                    "raw_patient_id": "P",
                    "raw_block_id": "B",
                    "raw_project": "X",
                    "duplicate_flag": "",
                    "identity_conflict": "no",
                    "evidence_grade": "E3_explicit",
                    "record_status": "MATCHED",
                    "notes": "",
                }
            ]
            write_rows(rows, output)
            first = output.read_bytes()
            write_rows(rows, output)
            self.assertEqual(first, output.read_bytes())
            with output.open(newline="", encoding="utf-8") as handle:
                parsed = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(parsed[0]["source_record_id"], "x")


if __name__ == "__main__":
    unittest.main()
