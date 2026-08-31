"""Tests for the three-lineage validation extractor (D-037)."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from scripts.r01_extract_validation_lineages import (
    LINEAGES,
    STCRC_SAMPLES,
    USZ_SAMPLES,
    extract,
)


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.fixture()
def staging(tmp_path: Path) -> Path:
    output = tmp_path / "staging"
    extract(ROOT, output, LINEAGES)
    return output


def test_lineage_unit_counts_and_identity_granularity(staging: Path) -> None:
    units = read_tsv(staging / "validation_lineages_physical_units.tsv")
    assert len(units) == 46
    by_study: dict[str, list[dict[str, str]]] = {}
    for row in units:
        by_study.setdefault(row["study_id"], []).append(row)
    assert set(by_study) == {"GEO::GSE175540", "ST_CRC_CMS", "TLS_VISIUM_USZ"}
    assert len(by_study["GEO::GSE175540"]) == 24
    assert len(by_study["ST_CRC_CMS"]) == 14
    assert len(by_study["TLS_VISIUM_USZ"]) == 8

    for row in by_study["GEO::GSE175540"]:
        assert row["record_status"] == "RESOLVED_INCLUDED_CANDIDATE"
        assert row["block_id"] == ""
        assert row["physical_specimen_id"] == row["source_record_id"]
        assert row["identity_granularity"] == "patient_linked_physical_specimen"
        assert row["block_equivalent_status"] == "ACCEPTED_BLOCK_EQUIVALENT"
        assert row["block_equivalent_basis"]
        assert row["identity_status"] == "PATIENT_SPECIMEN_CORROBORATED"
    # One sample per patient: 24 distinct aliases.
    assert len({row["patient_id"] for row in by_study["GEO::GSE175540"]}) == 24

    for row in by_study["ST_CRC_CMS"]:
        assert row["record_status"] == "RESOLVED_INCLUDED_CANDIDATE"
        assert row["block_id"]
        assert row["section_id"] in {"Rep1", "Rep2"}
        assert not row["physical_specimen_id"]
        assert not row["identity_granularity"]
    assert len({row["patient_id"] for row in by_study["ST_CRC_CMS"]}) == 7
    assert len({row["block_id"] for row in by_study["ST_CRC_CMS"]}) == 7

    for row in by_study["TLS_VISIUM_USZ"]:
        assert row["record_status"] == "RESOLVED_INCLUDED_CANDIDATE"
        assert row["block_id"]
        assert row["evidence_grade"] == "E3_explicit"
        assert not row["identity_granularity"]
    assert len({row["patient_id"] for row in by_study["TLS_VISIUM_USZ"]}) == 8


def test_evidence_covers_gate_required_fields(staging: Path) -> None:
    units = read_tsv(staging / "validation_lineages_physical_units.tsv")
    evidence = read_tsv(staging / "validation_lineages_identity_evidence.tsv")
    by_entity: dict[str, dict[str, str]] = {}
    for row in evidence:
        by_entity.setdefault(row["entity_id"], {})[row["field"]] = row
    for unit in units:
        entity = by_entity[unit["physical_unit_id"]]
        assert entity["patient_id"]["normalized_value"] == unit["patient_id"]
        if unit["study_id"] == "GEO::GSE175540":
            assert entity["physical_specimen_id"]["normalized_value"] == unit[
                "physical_specimen_id"
            ]
            block = entity["block_id"]
            assert block["normalized_value"] == ""
            assert block["raw_value"] == unit["block_equivalent_basis"]
        else:
            assert entity["block_id"]["normalized_value"] == unit["block_id"]
        for row in entity.values():
            assert row["raw_value"]
            assert row["metadata_key"]
            assert row["evidence_grade"] in {"E3_explicit", "E2_corroborated"}
            assert row["conflict_flag"] == "false"


def test_source_assets_verify_and_reference_existing_files(staging: Path) -> None:
    import hashlib

    assets = read_tsv(staging / "validation_lineages_source_assets.tsv")
    # KIRC: family soft + raw tar locator; STCRC: record + readme + annotation
    # zip + 14 sample-zip locators; USZ: record + deposit locator + 8 sidecars.
    assert len(assets) == 2 + 17 + 10
    evidence = read_tsv(staging / "validation_lineages_identity_evidence.tsv")
    asset_ids = {row["asset_id"] for row in assets}
    for row in evidence:
        assert row["source_asset"] in asset_ids
        assert row["metadata_path"]
    for row in assets:
        path = ROOT / row["path"]
        assert path.is_file(), row["path"]
        assert str(path.stat().st_size) == row["bytes"]
        if row["checksum_status"] == "sha256_verified":
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            assert digest == row["checksum"]
        else:
            assert row["checksum_status"] == "not_computed_large_local_locator"


def test_duplicate_groups_match_frozen_role_leakage_groups(staging: Path) -> None:
    duplicates = read_tsv(staging / "validation_lineages_duplicate_groups.tsv")
    assert len(duplicates) == 4
    kirc = {
        row["member_id"]
        for row in duplicates
        if row["duplicate_group_id"] == "same_source::GSE175540"
    }
    assert kirc == {"GEO::GSE175540", "atlas_s2::kirc::1713c35ce741"}
    by_study = {row["member_id"]: row for row in duplicates}
    assert by_study["GEO::GSE175540"]["leakage_group_id"] == "lineage::PRJNA732692"
    assert by_study["ST_CRC_CMS"]["leakage_group_id"] == "LINEAGE::ST_CRC_CMS"
    assert by_study["TLS_VISIUM_USZ"]["leakage_group_id"] == "LINEAGE::TLS_VISIUM_USZ"
    canonical = by_study["ST_CRC_CMS"]
    assert canonical["relation_type"] == "canonical_source_lineage"
    assert canonical["resolution"] == "canonical_study"


def test_missing_official_metadata_fails_closed(tmp_path: Path) -> None:
    empty = tmp_path / "empty_root"
    empty.mkdir()
    with pytest.raises(ValueError, match="family SOFT"):
        extract(empty, tmp_path / "out", ("GSE175540",))
    with pytest.raises((ValueError, OSError), match="ST_CRC_CMS"):
        extract(empty, tmp_path / "out2", ("ST_CRC_CMS",))
    with pytest.raises((ValueError, OSError), match="TLS_VISIUM_USZ"):
        extract(empty, tmp_path / "out3", ("TLS_VISIUM_USZ",))


def test_sample_inventories_are_exact() -> None:
    assert len(STCRC_SAMPLES) == 14
    assert len(USZ_SAMPLES) == 8
    patients = {
        (sample[: -len("_X")] if sample.endswith("_X") else sample).split("_")[1]
        for sample in STCRC_SAMPLES
    }
    assert len(patients) == 7
