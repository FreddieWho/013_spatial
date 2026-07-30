import csv
import tempfile
import unittest
from pathlib import Path

from scripts.r01_summarize_units import summarize_atlas_units, write_units


class SummarizeUnitsTests(unittest.TestCase):
    def test_groups_by_study_namespace_without_accepting_tls_as_gt(self):
        rows = [
            {
                "record_id": "A",
                "study_namespace": "study::1",
                "cancer_raw": "KIRC",
                "patient_id_raw": "P1",
                "platform_raw": "Visium",
                "preservation_raw": "FFPE",
                "reference_raw": "PMID:1",
                "data_availability_raw": "GSE1",
                "tls_presence_raw": "Yes",
                "tls_count_raw": "2",
            },
            {
                "record_id": "B",
                "study_namespace": "study::1",
                "cancer_raw": "KIRC",
                "patient_id_raw": "P2",
                "platform_raw": "Visium",
                "preservation_raw": "FFPE",
                "reference_raw": "PMID:1",
                "data_availability_raw": "GSE1",
                "tls_presence_raw": "No",
                "tls_count_raw": "0",
            },
        ]

        units = summarize_atlas_units(rows, Path("/nonexistent"))

        self.assertEqual(len(units), 1)
        unit = units[0]
        self.assertEqual(unit["sample_count"], "2")
        self.assertEqual(unit["patient_value_count"], "2")
        self.assertEqual(unit["potential_tls_positive_count"], "1")
        self.assertEqual(unit["gt_acceptance_status"], "NOT_ACCEPTED_R01")
        self.assertEqual(unit["role_freeze_status"], "UNFROZEN_BLOCK_UNKNOWN")

    def test_composite_patient_values_are_flagged(self):
        rows = [
            {
                "record_id": "A",
                "study_namespace": "study::1",
                "cancer_raw": "GBM",
                "patient_id_raw": "P1\nP1 & P2",
                "platform_raw": "Visium",
                "preservation_raw": "FF",
                "reference_raw": "",
                "data_availability_raw": "",
                "tls_presence_raw": "",
                "tls_count_raw": "",
            }
        ]
        unit = summarize_atlas_units(rows, Path("/nonexistent"))[0]
        self.assertEqual(unit["composite_patient_value_count"], "1")
        self.assertEqual(unit["identity_status"], "PATIENT_CONFLICT_BLOCK_UNKNOWN")

    def test_local_geo_directory_is_reported_but_not_overstated(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "data/GEO/GSE123").mkdir(parents=True)
            rows = [
                {
                    "record_id": "A",
                    "study_namespace": "study::1",
                    "cancer_raw": "X",
                    "patient_id_raw": "P1",
                    "platform_raw": "Visium",
                    "preservation_raw": "FF",
                    "reference_raw": "",
                    "data_availability_raw": " GSE123 ",
                    "tls_presence_raw": "",
                    "tls_count_raw": "",
                }
            ]
            unit = summarize_atlas_units(rows, root)[0]
            self.assertEqual(unit["local_data_status"], "LOCAL_DIRECTORY_PRESENT")

    def test_output_is_deterministic(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "units.tsv"
            rows = [
                {
                    "logical_unit_id": "x",
                    "cancer": "X",
                    "platforms": "Visium",
                    "preservation": "FF",
                    "sample_count": "1",
                    "patient_value_count": "1",
                    "composite_patient_value_count": "0",
                    "potential_tls_positive_count": "0",
                    "potential_tls_instance_count": "0",
                    "references": "",
                    "data_availability": "",
                    "local_data_status": "UNKNOWN",
                    "identity_status": "PATIENT_ONLY_BLOCK_UNKNOWN",
                    "gt_acceptance_status": "NOT_ACCEPTED_R01",
                    "role_freeze_status": "UNFROZEN_BLOCK_UNKNOWN",
                    "notes": "",
                }
            ]
            write_units(rows, output)
            first = output.read_bytes()
            write_units(rows, output)
            self.assertEqual(first, output.read_bytes())
            with output.open(newline="", encoding="utf-8") as handle:
                parsed = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(parsed[0]["logical_unit_id"], "x")


if __name__ == "__main__":
    unittest.main()
