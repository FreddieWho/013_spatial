from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "update_bioinf_data_index.py"


def test_generated_index_has_checksums_and_scope_guardrails(tmp_path: Path) -> None:
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--repo-root",
            str(ROOT),
            "--index-root",
            str(tmp_path),
            "--raw-root",
            str(ROOT / "infra" / "bioinf-data-index" / "raw"),
            "--generated-at",
            "2026-07-31T06:16:28+08:00",
        ],
        check=True,
    )

    with (tmp_path / "index.tsv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    assert all(len(row["sha256"]) == 64 for row in rows)
    assert all(row["license"] for row in rows)
    assert all(row["identity_coverage"] for row in rows)
    assert {
        row["license"]
        for row in rows
        if row["accession"] == "PMC10991508"
    } == {"CC BY"}
    assert all(row["retention_scope"] in {"metadata", "identity_crosswalk"} for row in rows)
    assert not any(
        Path(row["local_path"]).suffix.lower()
        in {".h5", ".h5ad", ".jpg", ".jpeg", ".gif", ".png", ".tif", ".tiff"}
        for row in rows
    )

    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["generated_at"] == "2026-07-31T06:16:28+08:00"
    assert summary["retained_expression_files"] == 0
    assert summary["retained_image_files"] == 0
    assert summary["retained_result_files"] == 0
    assert summary["indexed_file_count"] == len(rows)


def test_manifest_records_deleted_scope_contamination_incident(tmp_path: Path) -> None:
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--repo-root",
            str(ROOT),
            "--index-root",
            str(tmp_path),
            "--raw-root",
            str(ROOT / "infra" / "bioinf-data-index" / "raw"),
        ],
        check=True,
    )

    manifest = json.loads((tmp_path / "manifest.json").read_text())
    incident = manifest["scope_contamination_incidents"][0]
    assert incident["source"] == "PMC11508537 supplementaryFiles"
    assert incident["deleted"] is True
    assert incident["retained"] is False
    assert incident["package_sha256"] == (
        "55f7ec7db7baae38b8c5dd8ca5917e4f41701241a05d864e777489ef42d3ba29"
    )
    assert incident["metadata_only_package_sha256"] == (
        "72cdb27990ae78e7606ce00928c1e0feb11c297cd4c847b671c47ddfa43c9378"
    )
    assert incident["deleted_result_file_sha256"] == [
        "62dc15dc772f8e3f42f374a47fe39305a44f9e49bb451ab789d66a4ac769403b",
        "e61bc1511dcaeeb6efdeb5476fbe8048b520f272f390e7258322eb251ebc0b46",
    ]
    full_text_incidents = manifest["scope_contamination_incidents"][1:]
    assert len(full_text_incidents) == 4
    assert {
        incident["deleted_path"] for incident in full_text_incidents
    } == {
        "infra/bioinf-data-index/raw/pmc/PMC9862087_article.xml",
        "infra/bioinf-data-index/raw/pmc/PMC10439131_article.xml",
        "infra/bioinf-data-index/raw/pmc/PMC10991508_article.xml",
        "infra/bioinf-data-index/raw/pmc/PMC11508537_article.xml",
    }
    assert all(
        incident["sha256_status"] == "not_computed_before_deletion"
        and incident["deleted"] is True
        and incident["retained"] is False
        for incident in full_text_incidents
    )
