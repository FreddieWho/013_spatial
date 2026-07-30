import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts.r01_cross_source_duplicates import (
    canonicalize_url,
    find_hest_tenx_candidates,
    write_candidates,
)


class CrossSourceDuplicateTests(unittest.TestCase):
    def test_url_canonicalization_ignores_scheme_www_query_and_fragment(self):
        left = "https://www.example.org/data/sample/?download=1#top"
        right = "http://example.org/data/sample"
        self.assertEqual(canonicalize_url(left), canonicalize_url(right))

    def test_singleton_source_page_join_confirms_source_dataset_not_physical_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = root / "metadata"
            metadata.mkdir()
            (metadata / "TENX1.json").write_text(
                json.dumps(
                    {
                        "id": "TENX1",
                        "patient": "P1",
                        "subseries": "replicate 1",
                        "download_page_link1": "https://www.example.org/datasets/A/",
                    }
                ),
                encoding="utf-8",
            )
            manifest = root / "tenx.tsv"
            manifest.write_text(
                "dataset_dir\tpage_url\n"
                "A\thttp://example.org/datasets/A?source=manifest\n",
                encoding="utf-8",
            )

            rows = find_hest_tenx_candidates(metadata, manifest)

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["left_record_id"], "HEST::TENX1")
            self.assertEqual(rows[0]["right_record_id"], "TENX::A")
            self.assertEqual(
                rows[0]["relation_type"],
                "same_source_dataset_lineage",
            )
            self.assertEqual(rows[0]["evidence_grade"], "E2_corroborated")
            self.assertEqual(
                rows[0]["resolution"],
                "CONFIRMED_SOURCE_DATASET",
            )
            self.assertTrue(rows[0]["leakage_group_id"].startswith("LEAKAGE::"))

    def test_nonmatching_or_missing_urls_are_not_emitted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = root / "metadata"
            metadata.mkdir()
            (metadata / "A.json").write_text(
                json.dumps({"id": "A", "download_page_link1": ""}),
                encoding="utf-8",
            )
            manifest = root / "tenx.tsv"
            manifest.write_text(
                "dataset_dir\tpage_url\n"
                "B\thttps://example.org/B\n",
                encoding="utf-8",
            )
            self.assertEqual(find_hest_tenx_candidates(metadata, manifest), [])

    def test_output_is_deterministic(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "out.tsv"
            rows = [
                {
                    "candidate_group_id": "URL::1",
                    "left_record_id": "HEST::A",
                    "right_record_id": "TENX::A",
                    "relation_type": "same_source_page",
                    "evidence_grade": "E2_corroborated",
                    "resolution": "CONFIRMED_SOURCE_DATASET",
                    "leakage_group_id": "LEAKAGE::1",
                    "canonical_source_url": "example.org/A",
                    "left_raw_patient": "",
                    "left_raw_subseries": "",
                    "notes": "",
                }
            ]
            write_candidates(rows, output)
            first = output.read_bytes()
            write_candidates(rows, output)
            self.assertEqual(first, output.read_bytes())
            with output.open(newline="", encoding="utf-8") as handle:
                parsed = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(parsed[0]["candidate_group_id"], "URL::1")

    def test_audited_bundle_mapping_separates_confirmed_possible_and_no_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = root / "metadata"
            metadata.mkdir()
            common_url = "https://example.org/bundle"
            for record_id, region in [
                ("TENX105", "hKidney_cancer_section"),
                ("TENX106", "hKidney_nondiseased_section"),
            ]:
                (metadata / f"{record_id}.json").write_text(
                    json.dumps(
                        {
                            "id": record_id,
                            "region_name": region,
                            "download_page_link1": common_url,
                        }
                    ),
                    encoding="utf-8",
                )
            manifest = root / "tenx.tsv"
            manifest.write_text(
                "dataset_dir\tpage_url\n"
                "Xenium_V1_hKidney_cancer_section\thttps://example.org/bundle\n",
                encoding="utf-8",
            )

            rows = find_hest_tenx_candidates(metadata, manifest)
            by_hest = {row["left_record_id"]: row for row in rows}

            self.assertEqual(
                by_hest["HEST::TENX105"]["resolution"],
                "CONFIRMED_SOURCE_DATASET",
            )
            self.assertEqual(
                by_hest["HEST::TENX106"]["resolution"],
                "NO_EVIDENCE",
            )
            self.assertEqual(
                by_hest["HEST::TENX106"]["leakage_group_id"],
                "",
            )


if __name__ == "__main__":
    unittest.main()
