import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts.r01_extract_hest import (
    build_conflict_summary,
    build_field_coverage,
    extract_hest_records,
    write_tsv,
)


class ExtractHestTests(unittest.TestCase):
    def _write_inputs(self, root: Path) -> tuple[Path, Path]:
        selected = root / "selected_samples.tsv"
        selected.write_text(
            "id\torgan\tdisease_state\toncotree_code\tst_technology\tpreservation_method\n"
            "TENX147\tBowel\tCancer\tCOAD\tXenium\tFFPE\n"
            "ONLY_SELECTED\tBreast\tHealthy\t\tVisium\tFresh Frozen\n",
            encoding="utf-8",
        )
        metadata = root / "metadata"
        metadata.mkdir(exist_ok=True)
        (metadata / "TENX147.json").write_text(
            json.dumps(
                {
                    "id": "TENX147",
                    "patient": "Patient 1",
                    "subseries": "Xenium In Situ, Sample P5 CRC",
                    "study_link": "https://example.test/study",
                    "download_page_link1": "https://example.test/download",
                    "platform": None,
                    "st_technology": "Xenium",
                    "slide_id": float("nan"),
                    "z_step_size": None,
                }
            ),
            encoding="utf-8",
        )
        (metadata / "ONLY_METADATA.json").write_text(
            json.dumps(
                {
                    "id": "ONLY_METADATA",
                    "patient": "P2",
                    "subseries": "patient_2_visit_1",
                    "study_link": "",
                    "download_page_link1": None,
                    "st_technology": "Visium HD",
                }
            ),
            encoding="utf-8",
        )
        return selected, metadata

    def test_preserves_raw_fields_and_flags_only_explicit_p_number_conflict(self):
        with tempfile.TemporaryDirectory() as tmp:
            selected, metadata = self._write_inputs(Path(tmp))

            rows = extract_hest_records(selected, metadata)

            self.assertEqual(
                [row["source_record_id"] for row in rows],
                ["ONLY_METADATA", "ONLY_SELECTED", "TENX147"],
            )
            by_id = {row["source_record_id"]: row for row in rows}
            conflict = by_id["TENX147"]
            self.assertEqual(conflict["raw_patient"], "Patient 1")
            self.assertEqual(
                conflict["raw_subseries"], "Xenium In Situ, Sample P5 CRC"
            )
            self.assertEqual(conflict["raw_study_link"], "https://example.test/study")
            self.assertEqual(
                conflict["raw_download_page_link1"],
                "https://example.test/download",
            )
            self.assertEqual(conflict["raw_id"], "TENX147")
            self.assertEqual(conflict["raw_platform"], "")
            self.assertEqual(conflict["raw_slide_id"], "")
            self.assertEqual(conflict["raw_st_technology"], "Xenium")
            self.assertEqual(conflict["patient_explicit_p_tokens"], "P1")
            self.assertEqual(conflict["subseries_explicit_p_tokens"], "P5")
            self.assertEqual(conflict["identity_conflict"], "yes")
            self.assertEqual(conflict["evidence_grade"], "EC_conflict")
            self.assertEqual(
                conflict["conflict_reasons"],
                "patient_subseries_p_number_conflict",
            )

            corroborated = by_id["ONLY_METADATA"]
            self.assertEqual(corroborated["identity_conflict"], "no")
            self.assertEqual(corroborated["evidence_grade"], "E2_corroborated")
            self.assertEqual(corroborated["record_status"], "METADATA_ONLY")

            selected_only = by_id["ONLY_SELECTED"]
            self.assertEqual(selected_only["record_status"], "SELECTED_ONLY")
            self.assertEqual(selected_only["raw_patient"], "")
            self.assertEqual(selected_only["evidence_grade"], "E0_unknown")

            forbidden_inferences = {
                "patient_id",
                "block_id",
                "section_id",
                "serial_index",
                "z_position",
            }
            self.assertTrue(forbidden_inferences.isdisjoint(conflict))

    def test_coverage_and_conflict_summaries_are_auditable(self):
        with tempfile.TemporaryDirectory() as tmp:
            selected, metadata = self._write_inputs(Path(tmp))
            rows = extract_hest_records(selected, metadata)

            coverage = {
                row["field"]: row for row in build_field_coverage(rows)
            }
            self.assertEqual(coverage["raw_patient"]["nonempty_count"], "2")
            self.assertEqual(coverage["raw_patient"]["total_records"], "3")
            self.assertEqual(coverage["raw_platform"]["nonempty_count"], "0")

            conflicts = {
                row["conflict_type"]: row for row in build_conflict_summary(rows)
            }
            self.assertEqual(
                conflicts["patient_subseries_p_number_conflict"]["record_count"],
                "1",
            )
            self.assertEqual(
                conflicts["patient_subseries_p_number_conflict"]["record_ids"],
                "TENX147",
            )
            self.assertEqual(conflicts["any_EC_conflict"]["record_count"], "1")

    def test_rejects_duplicate_selected_ids_and_is_deterministic(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            selected, metadata = self._write_inputs(root)
            selected.write_text(
                "id\torgan\tdisease_state\toncotree_code\tst_technology\tpreservation_method\n"
                "A\tBowel\tCancer\tCOAD\tVisium\tFFPE\n"
                "A\tBowel\tCancer\tCOAD\tVisium\tFFPE\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "duplicate selected id: A"):
                extract_hest_records(selected, metadata)

            selected, metadata = self._write_inputs(root)
            rows = extract_hest_records(selected, metadata)
            output = root / "records.tsv"
            write_tsv(rows, output)
            first = output.read_bytes()
            write_tsv(list(reversed(rows)), output)
            self.assertEqual(first, output.read_bytes())
            with output.open(newline="", encoding="utf-8") as handle:
                parsed = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(parsed[0]["source_record_id"], "ONLY_METADATA")


if __name__ == "__main__":
    unittest.main()
