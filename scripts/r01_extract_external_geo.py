#!/usr/bin/env python3
"""Extract identity-only R-01 staging rows from approved official GEO metadata.

The extractor reads GEO SOFT metadata, one identity crosswalk supplement, and
tar member names. It never opens a tar member or interprets expression/image
content. A physical specimen is recorded separately from ``block_id``; the
latter remains unknown unless a future, documented policy accepts an explicit
block crosswalk.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
import tarfile
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from openpyxl import load_workbook


SOURCE_ASSET_FIELDS = [
    "asset_id",
    "source_namespace",
    "accession",
    "source_sample_id",
    "path",
    "asset_type",
    "format",
    "bytes",
    "checksum_status",
    "checksum",
    "processing_level",
    "parent_asset_id",
    "record_status",
]
PHYSICAL_UNIT_FIELDS = [
    "physical_unit_id",
    "source_namespace",
    "source_record_id",
    "study_id",
    "patient_id",
    "block_id",
    "section_id",
    "z_position",
    "section_order",
    "section_thickness",
    "platform",
    "identity_status",
    "evidence_grade",
    "record_status",
    "physical_specimen_id",
    "identity_granularity",
    "block_equivalent_status",
    "block_equivalent_basis",
]
IDENTITY_EVIDENCE_FIELDS = [
    "evidence_id",
    "entity_type",
    "entity_id",
    "field",
    "source_asset",
    "metadata_path",
    "metadata_key",
    "raw_value",
    "normalized_value",
    "evidence_grade",
    "interpretation",
    "conflict_flag",
]
DUPLICATE_GROUP_FIELDS = [
    "duplicate_group_id",
    "member_type",
    "member_id",
    "relation_type",
    "evidence_grade",
    "resolution",
    "leakage_group_id",
]

ACCESSIONS = ("GSE274103", "GSE274557", "GSE226997", "GSE211956")
ATLAS_STUDIES = {
    "GSE274103": "atlas_s2::paad::25fc49fa3903",
    "GSE274557": "atlas_s2::paad::37a9a031b038",
    "GSE226997": "atlas_s2::coad::87c0a055d20d",
    "GSE211956": "atlas_s2::ovca::1096a62ee6fb",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def repo_relative(path: Path, repo_root: Path) -> str:
    try:
        return path.resolve().relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def parse_soft_samples(path: Path) -> list[dict[str, object]]:
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("^SAMPLE = "):
            if current is not None:
                samples.append(current)
            current = {"accession": line.split(" = ", 1)[1], "fields": defaultdict(list)}
        elif current is not None and line.startswith("!") and " = " in line:
            key, value = line.split(" = ", 1)
            fields = current["fields"]
            assert isinstance(fields, defaultdict)
            fields[key].append(value)
    if current is not None:
        samples.append(current)
    return samples


def parse_series(path: Path) -> dict[str, list[str]]:
    fields: dict[str, list[str]] = defaultdict(list)
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("!Series_") and " = " in line:
            key, value = line.split(" = ", 1)
            fields[key].append(value)
    return dict(fields)


def one(fields: object, key: str) -> str:
    if not isinstance(fields, dict):
        return ""
    values = fields.get(key, [])
    return str(values[0]) if values else ""


def biosample(fields: object) -> str:
    if not isinstance(fields, dict):
        return ""
    for relation in fields.get("!Sample_relation", []):
        match = re.search(r"BioSample:\s*(?:https?://[^/]+/biosample/)?([A-Z0-9]+)", relation)
        if match:
            return match.group(1)
    return ""


def bioproject(series: dict[str, list[str]]) -> str:
    for relation in series.get("!Series_relation", []):
        match = re.search(r"BioProject:.*?(PRJNA\d+)", relation)
        if match:
            return match.group(1)
    return ""


def make_asset(
    path: Path,
    repo_root: Path,
    accession: str,
    kind: str,
    fmt: str,
    *,
    header_only: bool = False,
) -> dict[str, str]:
    return {
        "asset_id": f"external_geo::{accession}::{kind}",
        "source_namespace": "external_geo",
        "accession": accession,
        "source_sample_id": "",
        "path": repo_relative(path, repo_root),
        "asset_type": "local_tar_header_locator" if header_only else kind,
        "format": fmt,
        "bytes": str(path.stat().st_size),
        "checksum_status": (
            "not_computed_large_local_locator" if header_only else "sha256_verified"
        ),
        "checksum": "" if header_only else sha256(path),
        "processing_level": "header_only" if header_only else "metadata_only",
        "parent_asset_id": "",
        "record_status": "active",
    }


def verify_tar_gsms(tar_path: Path, gsms: Iterable[str]) -> None:
    expected = set(gsms)
    found: set[str] = set()
    with tarfile.open(tar_path, mode="r:") as archive:
        for member in archive:
            prefix = member.name.split("_", 1)[0]
            if prefix in expected:
                found.add(prefix)
    missing = sorted(expected - found)
    if missing:
        raise ValueError(f"{tar_path}: tar header missing GSM locator(s): {', '.join(missing)}")


def supplement_visium_map(path: Path) -> dict[str, str]:
    workbook = load_workbook(path, read_only=True, data_only=True)
    try:
        worksheet = workbook["Sample matrix"]
        rows = worksheet.iter_rows(values_only=True)
        header = next(rows)
        patient_col = header.index(None)
        gsm_col = header.index("Visium sample GEO")
        result: dict[str, str] = {}
        for row in rows:
            patient, gsm = row[patient_col], row[gsm_col]
            if patient and gsm:
                match = re.fullmatch(r"Patient\s+(\d+)", str(patient).strip())
                if not match:
                    raise ValueError(f"{path}: unexpected patient label {patient!r}")
                result[str(gsm).strip()] = f"P{match.group(1)}"
        return result
    finally:
        workbook.close()


def patient_from_title(accession: str, title: str) -> str:
    patterns = {
        "GSE274103": r"PDAC-p(\d+)",
        "GSE274557": r"Pt-(\d+)",
        "GSE226997": r"Visium result of patient\s+(\d+)",
    }
    match = re.search(patterns[accession], title, flags=re.IGNORECASE)
    return f"P{match.group(1)}" if match else ""


def basis_for(accession: str) -> str:
    return {
        "GSE274103": (
            "series: FFPE samples from five patient tissues; "
            "GSM title: patient; GEO relation: BioSample"
        ),
        "GSE274557": (
            "series: 55 individual samples from 13 patients; "
            "GSM title: patient-specimen; GEO relation: BioSample"
        ),
        "GSE226997": (
            "series: Visium of four CRC patients; GSM title: patient; "
            "GSM source: primary colorectal cancer; GEO relation: BioSample"
        ),
        "GSE211956": (
            "supplement: patient-to-Visium-GSM crosswalk; "
            "GSM source: ovarian tumour; GEO relation: BioSample"
        ),
    }[accession]


def validate_series(accession: str, series: dict[str, list[str]]) -> str:
    text = " ".join(
        series.get("!Series_summary", []) + series.get("!Series_overall_design", [])
    )
    project = bioproject(series)
    if not project:
        raise ValueError(f"{accession}: official series metadata lacks BioProject")
    required = {
        "GSE274103": ("FFPE", "five", "patient tissues"),
        "GSE274557": ("55 individual samples", "13 patients"),
        "GSE226997": ("four CRC patients",),
        "GSE211956": ("Visium",),
    }[accession]
    missing = [token for token in required if token.casefold() not in text.casefold()]
    if missing:
        raise ValueError(
            f"{accession}: official series tissue/patient evidence missing: {', '.join(missing)}"
        )
    return project


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", newline="", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    temporary.replace(path)


def extract(
    repo_root: Path,
    output_dir: Path,
    accessions: Iterable[str] = ACCESSIONS,
    *,
    skip_tar_check: bool = False,
) -> None:
    geo_root = repo_root / "infra" / "bioinf-data-index" / "raw" / "geo"
    supplement_path = (
        repo_root
        / "infra"
        / "bioinf-data-index"
        / "raw"
        / "supplements"
        / "PMC10991508_supplementary_data_1.xlsx"
    )
    source_assets: list[dict[str, str]] = []
    physical_units: list[dict[str, str]] = []
    evidence: list[dict[str, str]] = []
    duplicates: list[dict[str, str]] = []

    for accession in accessions:
        if accession not in ACCESSIONS:
            raise ValueError(f"unsupported accession: {accession}")
        gsm_path = geo_root / f"{accession}_gsm_quick.soft"
        series_path = geo_root / f"{accession}_series_quick.soft"
        if not gsm_path.is_file() or not series_path.is_file():
            raise ValueError(f"{accession}: official GEO GSM and series SOFT metadata are required")
        samples = parse_soft_samples(gsm_path)
        series = parse_series(series_path)
        project = validate_series(accession, series)

        supplement_map: dict[str, str] = {}
        if accession == "GSE211956":
            if not supplement_path.is_file():
                raise ValueError(f"{accession}: official patient-to-GSM supplement is required")
            supplement_map = supplement_visium_map(supplement_path)
            source_assets.append(
                make_asset(
                    supplement_path,
                    repo_root,
                    accession,
                    "patient_gsm_crosswalk",
                    "xlsx",
                )
            )

        selected: list[tuple[dict[str, object], str, str]] = []
        for sample in samples:
            gsm = str(sample["accession"])
            fields = sample["fields"]
            title = one(fields, "!Sample_title")
            source = one(fields, "!Sample_source_name_ch1")
            if accession == "GSE274557" and "PDX" in title.upper():
                continue
            if accession == "GSE211956":
                if gsm not in supplement_map:
                    continue
                patient = supplement_map[gsm]
            else:
                patient = patient_from_title(accession, title)
            specimen = biosample(fields)
            if not patient:
                raise ValueError(f"{accession} {gsm}: GSM patient evidence is missing")
            if not specimen:
                raise ValueError(f"{accession} {gsm}: GEO BioSample evidence is missing")
            if accession == "GSE226997" and "colorectal cancer" not in source.casefold():
                raise ValueError(
                    f"{accession} {gsm}: GSM tissue evidence 'Primary colorectal cancer' is missing"
                )
            selected.append((sample, patient, specimen))

        expected_count = {
            "GSE274103": 5,
            "GSE274557": 55,
            "GSE226997": 4,
            "GSE211956": 8,
        }[accession]
        if len(selected) != expected_count:
            raise ValueError(
                f"{accession}: expected {expected_count} qualifying spatial specimens, "
                f"found {len(selected)}"
            )

        tar_path = repo_root / "data" / "GEO" / accession / f"{accession}_RAW.tar"
        if not skip_tar_check:
            if not tar_path.is_file():
                raise ValueError(f"{accession}: local raw tar locator is missing")
            verify_tar_gsms(tar_path, (str(item[0]["accession"]) for item in selected))

        source_assets.extend(
            [
                make_asset(gsm_path, repo_root, accession, "gsm_soft", "soft"),
                make_asset(series_path, repo_root, accession, "series_soft", "soft"),
            ]
        )
        if tar_path.is_file():
            source_assets.append(
                make_asset(
                    tar_path,
                    repo_root,
                    accession,
                    "raw_tar_header",
                    "tar",
                    header_only=True,
                )
            )

        study_id = f"GEO::{accession}"
        group_id = f"same_source::{accession}"
        leakage_id = f"lineage::{project}"
        for member_id in (study_id, ATLAS_STUDIES[accession]):
            duplicates.append(
                {
                    "duplicate_group_id": group_id,
                    "member_type": "study",
                    "member_id": member_id,
                    "relation_type": "same_accession_same_source_lineage",
                    "evidence_grade": "E3_explicit",
                    "resolution": "canonical_geo_study; do_not_count_atlas_separately",
                    "leakage_group_id": leakage_id,
                }
            )

        for sample, patient, specimen in selected:
            gsm = str(sample["accession"])
            fields = sample["fields"]
            title = one(fields, "!Sample_title")
            unit_id = f"external_geo::{gsm}"
            physical_units.append(
                {
                    "physical_unit_id": unit_id,
                    "source_namespace": "external_geo",
                    "source_record_id": gsm,
                    "study_id": study_id,
                    "patient_id": f"{accession}::{patient}",
                    "block_id": "",
                    "section_id": "",
                    "z_position": "",
                    "section_order": "",
                    "section_thickness": "",
                    "platform": "10x_visium",
                    "identity_status": "PATIENT_SPECIMEN_CORROBORATED",
                    "evidence_grade": "E2_corroborated",
                    "record_status": "RESOLVED_INCLUDED_CANDIDATE",
                    "physical_specimen_id": specimen,
                    "identity_granularity": "patient_linked_physical_specimen",
                    "block_equivalent_status": "ACCEPTED_BLOCK_EQUIVALENT",
                    "block_equivalent_basis": basis_for(accession),
                }
            )
            evidence_source = (
                f"external_geo::{accession}::patient_gsm_crosswalk"
                if accession == "GSE211956"
                else f"external_geo::{accession}::gsm_soft"
            )
            evidence.extend(
                [
                    {
                        "evidence_id": f"external_geo::{gsm}::patient",
                        "entity_type": "physical_unit",
                        "entity_id": unit_id,
                        "field": "patient_id",
                        "source_asset": evidence_source,
                        "metadata_path": repo_relative(
                            supplement_path if accession == "GSE211956" else gsm_path,
                            repo_root,
                        ),
                        "metadata_key": (
                            "Sample matrix.Patient + Visium sample GEO"
                            if accession == "GSE211956"
                            else "!Sample_title"
                        ),
                        "raw_value": title,
                        "normalized_value": f"{accession}::{patient}",
                        "evidence_grade": (
                            "E3_explicit"
                            if accession == "GSE211956"
                            else "E2_corroborated"
                        ),
                        "interpretation": "official metadata patient-to-GSM mapping",
                        "conflict_flag": "false",
                    },
                    {
                        "evidence_id": f"external_geo::{gsm}::physical_specimen",
                        "entity_type": "physical_unit",
                        "entity_id": unit_id,
                        "field": "physical_specimen_id",
                        "source_asset": f"external_geo::{accession}::gsm_soft",
                        "metadata_path": repo_relative(gsm_path, repo_root),
                        "metadata_key": "!Sample_relation BioSample",
                        "raw_value": specimen,
                        "normalized_value": specimen,
                        "evidence_grade": "E2_corroborated",
                        "interpretation": "official GEO BioSample identifies physical specimen",
                        "conflict_flag": "false",
                    },
                    {
                        "evidence_id": f"external_geo::{gsm}::block_unknown",
                        "entity_type": "physical_unit",
                        "entity_id": unit_id,
                        "field": "block_id",
                        "source_asset": f"external_geo::{accession}::series_soft",
                        "metadata_path": repo_relative(series_path, repo_root),
                        "metadata_key": "!Series_summary + !Series_overall_design",
                        "raw_value": basis_for(accession),
                        "normalized_value": "",
                        "evidence_grade": "E2_corroborated",
                        "interpretation": (
                            "physical specimen is a candidate block-equivalent; "
                            "no block identifier is asserted"
                        ),
                        "conflict_flag": "false",
                    },
                ]
            )

    write_tsv(output_dir / "external_geo_source_assets.tsv", SOURCE_ASSET_FIELDS, source_assets)
    write_tsv(output_dir / "external_geo_physical_units.tsv", PHYSICAL_UNIT_FIELDS, physical_units)
    write_tsv(
        output_dir / "external_geo_identity_evidence.tsv", IDENTITY_EVIDENCE_FIELDS, evidence
    )
    write_tsv(
        output_dir / "external_geo_duplicate_groups.tsv", DUPLICATE_GROUP_FIELDS, duplicates
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("infra/sample-registry/staging"),
    )
    parser.add_argument("--accession", action="append", choices=ACCESSIONS)
    parser.add_argument("--skip-tar-check", action="store_true")
    args = parser.parse_args()
    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir
    if not output_dir.is_absolute():
        output_dir = repo_root / output_dir
    try:
        extract(
            repo_root,
            output_dir,
            args.accession or ACCESSIONS,
            skip_tar_check=args.skip_tar_check,
        )
    except (OSError, ValueError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
