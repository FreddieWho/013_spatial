import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts.r01_inventory_metadata import inspect_source, run_inventory


class InventoryMetadataTests(unittest.TestCase):
    def test_table_inventory_reports_rows_and_identity_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = root / "samples.tsv"
            table.write_text(
                "sample_id\tpatient_id\tsection_id\tplatform\n"
                "S1\tP1\tA\tVisium\n"
                "S2\tP1\tB\tVisium\n",
                encoding="utf-8",
            )

            record = inspect_source(
                root,
                {
                    "source_id": "samples",
                    "kind": "table",
                    "path": "samples.tsv",
                    "notes": "fixture",
                },
            )

            self.assertEqual(record["status"], "PRESENT")
            self.assertEqual(record["member_count"], 1)
            self.assertEqual(record["row_count"], 2)
            self.assertEqual(
                record["identity_fields"],
                "patient_id;sample_id;section_id",
            )
            self.assertEqual(record["section_fields"], "section_id")

    def test_json_collection_unions_keys_and_flags_derived_gt_risk(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = root / "metadata"
            metadata.mkdir()
            (metadata / "a.json").write_text(
                json.dumps({"id": "A", "patient": "P1", "tls_score": 0.2}),
                encoding="utf-8",
            )
            (metadata / "b.json").write_text(
                json.dumps({"id": "B", "section_order": 2, "annotation": "x"}),
                encoding="utf-8",
            )

            record = inspect_source(
                root,
                {
                    "source_id": "jsons",
                    "kind": "json_collection",
                    "glob": "metadata/*.json",
                    "notes": "fixture",
                },
            )

            self.assertEqual(record["member_count"], 2)
            self.assertEqual(record["row_count"], 2)
            self.assertEqual(
                record["gt_risk_fields"],
                "annotation;tls_score",
            )
            self.assertEqual(record["section_fields"], "section_order")

    def test_inventory_is_deterministic_and_marks_missing_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = root / "catalog.json"
            output = root / "inventory.tsv"
            config.write_text(
                json.dumps(
                    [
                        {
                            "source_id": "missing",
                            "kind": "table",
                            "path": "absent.tsv",
                            "notes": "fixture",
                        }
                    ]
                ),
                encoding="utf-8",
            )

            first = run_inventory(root, config, output)
            first_bytes = output.read_bytes()
            second = run_inventory(root, config, output)

            self.assertEqual(first, second)
            self.assertEqual(first_bytes, output.read_bytes())
            with output.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["status"], "MISSING")


if __name__ == "__main__":
    unittest.main()
