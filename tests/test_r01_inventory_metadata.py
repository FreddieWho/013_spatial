import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

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
                json.dumps(
                    {
                        "id": "A",
                        "patient": "P1",
                        "tls_score": 0.2,
                        "z_step_size": 5,
                    }
                ),
                encoding="utf-8",
            )
            (metadata / "b.json").write_text(
                json.dumps(
                    {
                        "id": "B",
                        "section_order": 2,
                        "annotation": "x",
                        "location": "tumor",
                    }
                ),
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
                "annotation;location;tls_score",
            )
            self.assertEqual(
                record["section_fields"],
                "section_order;z_step_size",
            )

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

    def test_disallowed_metadata_suffix_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            unsafe = root / "expression.h5ad"
            unsafe.write_bytes(b"not opened")

            with self.assertRaisesRegex(ValueError, "disallowed suffix"):
                inspect_source(
                    root,
                    {
                        "source_id": "unsafe",
                        "kind": "file",
                        "path": "expression.h5ad",
                    },
                )

    def test_oversized_metadata_is_rejected_before_hashing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = root / "large.tsv"
            table.write_text("id\nA\n", encoding="utf-8")

            with patch(
                "scripts.r01_inventory_metadata.MAX_METADATA_BYTES",
                1,
            ):
                with self.assertRaisesRegex(ValueError, "metadata file too large"):
                    inspect_source(
                        root,
                        {
                            "source_id": "large",
                            "kind": "table",
                            "path": "large.tsv",
                        },
                    )

    def test_output_must_remain_inside_project_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "root"
            root.mkdir()
            config = root / "catalog.json"
            config.write_text("[]", encoding="utf-8")
            outside = root.parent / "outside.tsv"

            with self.assertRaisesRegex(ValueError, "output escapes project root"):
                run_inventory(root, config, outside)

    def test_output_cannot_overwrite_catalog_or_metadata_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = root / "source.tsv"
            table.write_text("id\nA\n", encoding="utf-8")
            config = root / "catalog.json"
            config.write_text(
                json.dumps(
                    [
                        {
                            "source_id": "source",
                            "kind": "table",
                            "path": "source.tsv",
                        }
                    ]
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "output collides"):
                run_inventory(root, config, table)
            with self.assertRaisesRegex(ValueError, "output collides"):
                run_inventory(root, config, config)

    def test_collection_total_size_is_limited(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = root / "metadata"
            metadata.mkdir()
            (metadata / "a.json").write_text('{"id":"A"}', encoding="utf-8")
            (metadata / "b.json").write_text('{"id":"B"}', encoding="utf-8")

            with patch(
                "scripts.r01_inventory_metadata.MAX_COLLECTION_BYTES",
                1,
            ):
                with self.assertRaisesRegex(ValueError, "metadata collection too large"):
                    inspect_source(
                        root,
                        {
                            "source_id": "collection",
                            "kind": "json_collection",
                            "glob": "metadata/*.json",
                        },
                    )


if __name__ == "__main__":
    unittest.main()
