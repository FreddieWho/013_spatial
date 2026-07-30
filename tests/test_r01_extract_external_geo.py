from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "r01_extract_external_geo.py"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_extracts_four_canonical_studies_without_inventing_blocks(tmp_path: Path) -> None:
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--repo-root",
            str(ROOT),
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
    )

    physical = read_tsv(tmp_path / "external_geo_physical_units.tsv")
    assert len(physical) == 72
    assert {row["study_id"] for row in physical} == {
        "GEO::GSE211956",
        "GEO::GSE226997",
        "GEO::GSE274103",
        "GEO::GSE274557",
    }
    assert all(row["block_id"] == "" for row in physical)
    assert all(row["physical_specimen_id"] for row in physical)
    assert all(
        row["identity_granularity"] == "patient_linked_physical_specimen"
        for row in physical
    )
    assert all(
        row["block_equivalent_status"] == "ACCEPTED_BLOCK_EQUIVALENT"
        for row in physical
    )
    assert all(row["identity_status"] == "PATIENT_SPECIMEN_CORROBORATED" for row in physical)
    assert all(row["record_status"] == "RESOLVED_INCLUDED_CANDIDATE" for row in physical)
    assert all(row["evidence_grade"] == "E2_corroborated" for row in physical)
    assert not any("PDX" in row["source_record_id"].upper() for row in physical)
    assert not any("GSM650610" in row["source_record_id"] for row in physical)

    per_study = {}
    for row in physical:
        per_study.setdefault(row["study_id"], []).append(row)
    assert len(per_study["GEO::GSE274103"]) == 5
    assert len(per_study["GEO::GSE274557"]) == 55
    assert len(per_study["GEO::GSE226997"]) == 4
    assert len(per_study["GEO::GSE211956"]) == 8


def test_records_same_accession_atlas_lineage_and_header_only_assets(
    tmp_path: Path,
) -> None:
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--repo-root",
            str(ROOT),
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
    )

    duplicates = read_tsv(tmp_path / "external_geo_duplicate_groups.tsv")
    assert len(duplicates) == 8
    assert {row["duplicate_group_id"] for row in duplicates} == {
        "same_source::GSE211956",
        "same_source::GSE226997",
        "same_source::GSE274103",
        "same_source::GSE274557",
    }
    assert all(row["relation_type"] == "same_accession_same_source_lineage" for row in duplicates)

    assets = read_tsv(tmp_path / "external_geo_source_assets.tsv")
    tar_rows = [row for row in assets if row["asset_type"] == "local_tar_header_locator"]
    assert len(tar_rows) == 4
    assert all(row["processing_level"] == "header_only" for row in tar_rows)
    assert all(row["checksum_status"] == "not_computed_large_local_locator" for row in tar_rows)
    assert not any(row["format"].lower() in {"jpg", "jpeg", "gif", "png"} for row in assets)


def test_gse226997_fails_closed_when_biosample_evidence_is_missing(
    tmp_path: Path,
) -> None:
    fixture_root = tmp_path / "repo"
    geo_dir = fixture_root / "infra" / "bioinf-data-index" / "raw" / "geo"
    geo_dir.mkdir(parents=True)
    source = (ROOT / "infra/bioinf-data-index/raw/geo/GSE226997_gsm_quick.soft").read_text()
    source = "\n".join(
        line for line in source.splitlines() if "!Sample_relation = BioSample:" not in line
    )
    (geo_dir / "GSE226997_gsm_quick.soft").write_text(source)
    (geo_dir / "GSE226997_series_quick.soft").write_bytes(
        (ROOT / "infra/bioinf-data-index/raw/geo/GSE226997_series_quick.soft").read_bytes()
    )

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--repo-root",
            str(fixture_root),
            "--output-dir",
            str(tmp_path / "out"),
            "--accession",
            "GSE226997",
            "--skip-tar-check",
        ],
        text=True,
        capture_output=True,
    )
    assert result.returncode != 0
    assert "GSE226997" in result.stderr
    assert "BioSample" in result.stderr
