#!/usr/bin/env python3
"""Extract identity-only R-01 staging rows for the three validation lineages.

Approved under D-037, this extractor registers three public lineages that carry
human-annotated ground truth (D-036: any human annotator counts as
high-confidence at this stage):

- ``GEO::GSE175540`` (Meylan et al. 2022, Immunity, PMID 35231421): 24 ccRCC
  Visium samples. Identity uses the accepted patient-linked physical-specimen
  granularity (one Visium sample per tumor; Atlas Table S2 lists 24 distinct
  patients R_P1-R_P24, one sample each). No BioSample was deposited (patient
  privacy; verified: the family SOFT carries no ``!Sample_relation`` and
  BioProject PRJNA732692 links zero BioSamples), so the GSM accession itself is
  the physical-specimen identifier.
- ``ST_CRC_CMS`` (Valdeolivas et al. 2024, npj Precis Oncol, Zenodo 7760264):
  seven CRC patients, two serial sections per patient (technical replicates),
  with pathologist spot categorization. Explicit patient/block identity.
- ``TLS_VISIUM_USZ`` (Zenodo 14620362): eight FFPE tumor sections (kidney
  KC1-3, lung LC1-5), one Visium section per tumor, expert-annotated spots.
  Explicit identity asserted per tumor; no inter-sample patient linkage is
  stated in the deposit.

The extractor reads official metadata only (GEO SOFT, Zenodo record JSON, the
authors' GitHub README, deposit file listings, and 67-byte per-sample JSON
sidecars). It never opens expression matrices, annotation payloads, or images.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Iterable

try:
    from scripts.r01_build_registry import (
        DUPLICATE_GROUP_FIELDS,
        IDENTITY_EVIDENCE_FIELDS,
        PHYSICAL_UNIT_FIELDS,
        SOURCE_ASSET_FIELDS,
    )
    from scripts.r01_extract_external_geo import (
        parse_soft_samples,
        repo_relative,
        sha256,
        verify_tar_gsms,
        write_tsv,
    )
except ModuleNotFoundError:  # Direct execution via ``python scripts/...``.
    from r01_build_registry import (
        DUPLICATE_GROUP_FIELDS,
        IDENTITY_EVIDENCE_FIELDS,
        PHYSICAL_UNIT_FIELDS,
        SOURCE_ASSET_FIELDS,
    )
    from r01_extract_external_geo import (
        parse_soft_samples,
        repo_relative,
        sha256,
        verify_tar_gsms,
        write_tsv,
    )


LINEAGES = ("GSE175540", "ST_CRC_CMS", "TLS_VISIUM_USZ")

# --- GSE175540 (KIRC) -------------------------------------------------------
KIRC_STUDY = "GEO::GSE175540"
KIRC_ATLAS_PEER = "atlas_s2::kirc::1713c35ce741"
KIRC_PROJECT = "PRJNA732692"
KIRC_LEAKAGE = f"lineage::{KIRC_PROJECT}"
KIRC_TITLE_RE = re.compile(r"^ccRCC tumor\s+((?:ffpe|frozen)_[a-z]_\d+)$")
KIRC_EXPECTED_SAMPLES = 24
KIRC_BLOCK_BASIS = (
    "series: spatial transcriptomics of ccRCC human tumors, one Visium sample "
    "per tumor (12 FFPE + 12 fresh frozen); GSM title: sample alias; atlas "
    "Table S2: 24 distinct patients R_P1-R_P24 with one sample each; no "
    "BioSample deposited (patient privacy), GSM accession is the specimen "
    "identifier"
)

# --- ST_CRC_CMS (Valdeolivas) ----------------------------------------------
STCRC_STUDY = "ST_CRC_CMS"
STCRC_LEAKAGE = "LINEAGE::ST_CRC_CMS"
STCRC_ZENODO = "7760264"
STCRC_SAMPLES = (
    "SN048_A416371_Rep1",
    "SN048_A416371_Rep2",
    "SN048_A121573_Rep1",
    "SN048_A121573_Rep2",
    "SN84_A120838_Rep1",
    "SN84_A120838_Rep2",
    "SN123_A551763_Rep1",
    "SN123_A595688_Rep1",
    "SN123_A798015_Rep1",
    "SN123_A938797_Rep1_X",
    "SN124_A551763_Rep2",
    "SN124_A595688_Rep2",
    "SN124_A798015_Rep2",
    "SN124_A938797_Rep2",
)
STCRC_PATIENTS = (
    "A416371",
    "A121573",
    "A120838",
    "A938797",
    "A595688",
    "A798015",
    "A551763",
)
STCRC_README_QUOTE = (
    "two serial sections per patient to generate technical replicates"
)
STCRC_WSI_GROUPS = {
    "A121573": "(samples A121573_Rep1, A121573_Rep2, A416371_Rep1 and A416371_Rep2)",
    "A416371": "(samples A121573_Rep1, A121573_Rep2, A416371_Rep1 and A416371_Rep2)",
    "A120838": "(samples A120838_Rep1 and A120838_Rep2)",
    "A551763": "(samples A551763_Rep1, A595688_Rep1, A798015_Rep1, A938797_Rep1)",
    "A595688": "(samples A551763_Rep1, A595688_Rep1, A798015_Rep1, A938797_Rep1)",
    "A798015": "(samples A551763_Rep1, A595688_Rep1, A798015_Rep1, A938797_Rep1)",
    "A938797": "(samples A551763_Rep1, A595688_Rep1, A798015_Rep1, A938797_Rep1)",
}

# --- TLS_VISIUM_USZ ---------------------------------------------------------
USZ_STUDY = "TLS_VISIUM_USZ"
USZ_LEAKAGE = "LINEAGE::TLS_VISIUM_USZ"
USZ_ZENODO = "14620362"
USZ_SAMPLES = ("KC1", "KC2", "KC3", "LC1", "LC2", "LC3", "LC4", "LC5")
USZ_TUMOR_QUOTE = "FFPE sections from kidney (3) and lung (5) tumors"
USZ_SAMPLE_QUOTE = (
    "Data samples include KC[1-3] indicating kidney cancer and "
    "LC[1-5] indicating lung cancer"
)


def _hash16(raw: str) -> str:
    import hashlib

    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:16]


def _asset(
    path: Path,
    repo_root: Path,
    namespace: str,
    accession: str,
    kind: str,
    fmt: str,
    *,
    header_only: bool = False,
) -> dict[str, str]:
    return {
        "asset_id": f"{namespace}::{accession}::{kind}",
        "source_namespace": namespace,
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


def _normalized_description(record_path: Path) -> str:
    record = json.loads(record_path.read_text(encoding="utf-8"))
    description = str(record["metadata"]["description"])
    for entity, char in (("&nbsp;", " "), ("&mu;", "μ"), ("&amp;", "&")):
        description = description.replace(entity, char)
    return re.sub(r"<[^>]+>", "", description), record


def _atlas_kirc_patient_count(repo_root: Path) -> int:
    """Atlas Table S2 KIRC rows sourced from GSE175540 (corroboration check)."""
    from openpyxl import load_workbook

    workbook = load_workbook(
        repo_root / "paper/tables/science.adz2742_tables_s1_to_s8.xlsx",
        read_only=True,
        data_only=True,
    )
    try:
        sheet = workbook["Table S2"]
        current_availability = ""
        patients: set[str] = set()
        ffpe = frozen = 0
        for values in sheet.iter_rows(min_row=4, values_only=True):
            if values[11] is not None:
                current_availability = str(values[11]).strip()
            if current_availability != "GSE175540":
                continue
            if values[0] == "BRCA":
                break
            patient = str(values[2] or "").strip()
            sample = str(values[3] or "").strip()
            if not patient or not sample:
                raise ValueError("GSE175540: atlas Table S2 KIRC row lacks patient/sample")
            patients.add(patient)
            if values[7] == "FFPE":
                ffpe += 1
            elif values[7] == "Fresh Frozen":
                frozen += 1
        if len(patients) != 24 or ffpe != 12 or frozen != 12:
            raise ValueError(
                "GSE175540: atlas Table S2 corroboration drifted: "
                f"patients={len(patients)} ffpe={ffpe} frozen={frozen} (expected 24/12/12)"
            )
        return len(patients)
    finally:
        workbook.close()


def extract_kirc(
    repo_root: Path,
) -> tuple[
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    import gzip
    import tempfile

    soft_path = repo_root / "data/GEO/GSE175540/GSE175540_family.soft.gz"
    if not soft_path.is_file():
        raise ValueError("GSE175540: official family SOFT metadata is required")
    # parse_soft_samples expects a plain-text path; decompress to a temp file.
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as src, tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", suffix=".soft", delete=False
    ) as dst:
        dst.write(src.read())
        plain = Path(dst.name)
    try:
        samples = parse_soft_samples(plain)
        series_text = plain.read_text(encoding="utf-8", errors="replace")
    finally:
        plain.unlink()
    if "PRJNA732692" not in series_text:
        raise ValueError("GSE175540: official series metadata lacks BioProject PRJNA732692")
    if re.search(r"!Sample_relation\s*=\s*BioSample", series_text):
        raise ValueError("GSE175540: unexpected BioSample relation; specimen policy needs review")
    _atlas_kirc_patient_count(repo_root)

    selected: list[tuple[str, str, str]] = []  # (gsm, alias, title)
    for sample in samples:
        gsm = str(sample["accession"])
        fields = sample["fields"]
        title_values = fields.get("!Sample_title", [])
        title = str(title_values[0]) if title_values else ""
        match = KIRC_TITLE_RE.match(title)
        if not match:
            raise ValueError(f"GSE175540 {gsm}: unexpected sample title {title!r}")
        selected.append((gsm, match.group(1), title))
    aliases = [alias for _, alias, _ in selected]
    if len(selected) != KIRC_EXPECTED_SAMPLES or len(set(aliases)) != len(aliases):
        raise ValueError(
            f"GSE175540: expected {KIRC_EXPECTED_SAMPLES} uniquely aliased samples, "
            f"observed {len(selected)}"
        )
    if sum(a.startswith("ffpe") for a in aliases) != 12:
        raise ValueError("GSE175540: expected 12 FFPE sample aliases")
    if sum(a.startswith("frozen") for a in aliases) != 12:
        raise ValueError("GSE175540: expected 12 frozen sample aliases")

    assets = [
        _asset(soft_path, repo_root, "external_geo", "GSE175540", "family_soft", "soft")
    ]
    tar_path = repo_root / "data/GEO/GSE175540/GSE175540_RAW.tar"
    if tar_path.is_file():
        verify_tar_gsms(tar_path, (gsm for gsm, _, _ in selected))
        assets.append(
            _asset(
                tar_path,
                repo_root,
                "external_geo",
                "GSE175540",
                "raw_tar_header",
                "tar",
                header_only=True,
            )
        )

    duplicates: list[dict[str, str]] = []
    for member_id in (KIRC_STUDY, KIRC_ATLAS_PEER):
        duplicates.append(
            {
                "duplicate_group_id": "same_source::GSE175540",
                "member_type": "study",
                "member_id": member_id,
                "relation_type": "same_accession_same_source_lineage",
                "evidence_grade": "E3_explicit",
                "resolution": "canonical_geo_study; do_not_count_atlas_separately",
                "leakage_group_id": KIRC_LEAKAGE,
            }
        )

    units: list[dict[str, str]] = []
    evidence: list[dict[str, str]] = []
    soft_rel = repo_relative(soft_path, repo_root)
    soft_asset = "external_geo::GSE175540::family_soft"
    for gsm, alias, title in sorted(selected):
        unit_id = f"external_geo::{gsm}"
        patient_id = f"GSE175540::{alias}"
        units.append(
            {
                "physical_unit_id": unit_id,
                "source_namespace": "external_geo",
                "source_record_id": gsm,
                "study_id": KIRC_STUDY,
                "patient_id": patient_id,
                "block_id": "",
                "section_id": "",
                "z_position": "",
                "section_order": "",
                "section_thickness": "",
                "platform": "10x_visium",
                "identity_status": "PATIENT_SPECIMEN_CORROBORATED",
                "evidence_grade": "E2_corroborated",
                "record_status": "RESOLVED_INCLUDED_CANDIDATE",
                "physical_specimen_id": gsm,
                "identity_granularity": "patient_linked_physical_specimen",
                "block_equivalent_status": "ACCEPTED_BLOCK_EQUIVALENT",
                "block_equivalent_basis": KIRC_BLOCK_BASIS,
            }
        )
        evidence.extend(
            [
                {
                    "evidence_id": f"external_geo::{gsm}::patient",
                    "entity_type": "physical_unit",
                    "entity_id": unit_id,
                    "field": "patient_id",
                    "source_asset": soft_asset,
                    "metadata_path": soft_rel,
                    "metadata_key": "!Sample_title",
                    "raw_value": title,
                    "normalized_value": patient_id,
                    "evidence_grade": "E2_corroborated",
                    "interpretation": (
                        "official GEO sample alias identifies one ccRCC tumor "
                        "sample; atlas Table S2 lists 24 distinct patients "
                        "R_P1-R_P24 with one sample each"
                    ),
                    "conflict_flag": "false",
                },
                {
                    "evidence_id": f"external_geo::{gsm}::physical_specimen",
                    "entity_type": "physical_unit",
                    "entity_id": unit_id,
                    "field": "physical_specimen_id",
                    "source_asset": soft_asset,
                    "metadata_path": soft_rel,
                    "metadata_key": "!Sample_geo_accession",
                    "raw_value": gsm,
                    "normalized_value": gsm,
                    "evidence_grade": "E3_explicit",
                    "interpretation": (
                        "no BioSample was deposited for this series (patient "
                        "privacy; BioProject PRJNA732692 links zero BioSamples), "
                        "so the GSM accession is the only official "
                        "specimen-level identifier"
                    ),
                    "conflict_flag": "false",
                },
                {
                    "evidence_id": f"external_geo::{gsm}::block_unknown",
                    "entity_type": "physical_unit",
                    "entity_id": unit_id,
                    "field": "block_id",
                    "source_asset": soft_asset,
                    "metadata_path": soft_rel,
                    "metadata_key": "!Series_summary + !Series_overall_design + !Sample_title",
                    "raw_value": KIRC_BLOCK_BASIS,
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
    return assets, units, evidence, duplicates


def extract_stcrc(
    repo_root: Path,
) -> tuple[
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    base = repo_root / "data/other_sources/zenodo_st_crc_cms"
    record_path = base / "zenodo_record_7760264.json"
    readme_path = base / "ST_CRC_CMS_README.md"
    annotation_zip = base / "Pathology_SpotAnnotations.zip"
    for path in (record_path, readme_path, annotation_zip):
        if not path.is_file():
            raise ValueError(f"ST_CRC_CMS: required official metadata asset {path} is missing")
    description, record = _normalized_description(record_path)
    readme = readme_path.read_text(encoding="utf-8", errors="replace")
    if STCRC_README_QUOTE not in readme:
        raise ValueError("ST_CRC_CMS: GitHub README no longer states two serial sections per patient")
    for patient, quote in STCRC_WSI_GROUPS.items():
        if quote not in description:
            raise ValueError(f"ST_CRC_CMS: Zenodo record lost the WSI grouping for {patient}")
    file_keys = {str(item.get("key", "")) for item in record.get("files", [])}
    for sample in STCRC_SAMPLES:
        if f"{sample}.zip" not in file_keys:
            raise ValueError(f"ST_CRC_CMS: Zenodo record no longer lists {sample}.zip")
    if "Pathology_SpotAnnotations.zip" not in file_keys:
        raise ValueError("ST_CRC_CMS: Zenodo record no longer lists Pathology_SpotAnnotations.zip")

    assets = [
        _asset(record_path, repo_root, STCRC_STUDY, STCRC_STUDY, "zenodo_record", "json"),
        _asset(readme_path, repo_root, STCRC_STUDY, STCRC_STUDY, "github_readme", "md"),
        _asset(
            annotation_zip, repo_root, STCRC_STUDY, STCRC_STUDY, "pathology_annotations", "zip"
        ),
    ]
    for sample in STCRC_SAMPLES:
        zip_path = base / f"{sample}.zip"
        if not zip_path.is_file():
            raise ValueError(f"ST_CRC_CMS: local sample zip {zip_path} is missing")
        assets.append(
            _asset(
                zip_path,
                repo_root,
                STCRC_STUDY,
                STCRC_STUDY,
                f"sample_zip::{sample}",
                "zip",
                header_only=True,
            )
        )

    duplicates = [
        {
            "duplicate_group_id": f"SOURCE_LINEAGE::{STCRC_STUDY}",
            "member_type": "study",
            "member_id": STCRC_STUDY,
            "relation_type": "canonical_source_lineage",
            "evidence_grade": "E2_corroborated",
            "resolution": "canonical_study",
            "leakage_group_id": STCRC_LEAKAGE,
        }
    ]

    units: list[dict[str, str]] = []
    evidence: list[dict[str, str]] = []
    record_rel = repo_relative(record_path, repo_root)
    readme_rel = repo_relative(readme_path, repo_root)
    record_asset = f"{STCRC_STUDY}::{STCRC_STUDY}::zenodo_record"
    readme_asset = f"{STCRC_STUDY}::{STCRC_STUDY}::github_readme"
    per_patient: dict[str, list[str]] = {}
    for sample in STCRC_SAMPLES:
        stem = sample[: -len("_X")] if sample.endswith("_X") else sample
        _slide, patient, section = stem.split("_", 2)
        del _slide
        if patient not in STCRC_PATIENTS:
            raise ValueError(f"ST_CRC_CMS: unexpected patient alias in {sample}")
        if section not in {"Rep1", "Rep2"}:
            raise ValueError(f"ST_CRC_CMS: unexpected section label in {sample}")
        per_patient.setdefault(patient, []).append(section)
        unit_id = f"{STCRC_STUDY}::{sample}"
        patient_id = f"{STCRC_STUDY}::PATIENT::{_hash16(patient)}"
        block_id = f"{STCRC_STUDY}::BLOCK::{_hash16(patient + '|tumor_block')}"
        units.append(
            {
                "physical_unit_id": unit_id,
                "source_namespace": STCRC_STUDY,
                "source_record_id": f"{sample}.zip",
                "study_id": STCRC_STUDY,
                "patient_id": patient_id,
                "block_id": block_id,
                "section_id": section,
                "z_position": "",
                "section_order": "",
                "section_thickness": "",
                "platform": "10x_visium",
                "identity_status": "PATIENT_BLOCK_SECTION_CORROBORATED",
                "evidence_grade": "E2_corroborated",
                "record_status": "RESOLVED_INCLUDED_CANDIDATE",
                "physical_specimen_id": "",
                "identity_granularity": "",
                "block_equivalent_status": "",
                "block_equivalent_basis": "",
            }
        )
        evidence.extend(
            [
                {
                    "evidence_id": f"{STCRC_STUDY}::{sample}::patient",
                    "entity_type": "physical_unit",
                    "entity_id": unit_id,
                    "field": "patient_id",
                    "source_asset": record_asset,
                    "metadata_path": record_rel,
                    "metadata_key": "files[].key + metadata.description (WSI grouping)",
                    "raw_value": f"{sample}.zip {STCRC_WSI_GROUPS[patient]}",
                    "normalized_value": patient_id,
                    "evidence_grade": "E2_corroborated",
                    "interpretation": (
                        f"deposit file key embeds patient alias {patient}; the "
                        "record description groups this patient's replicate "
                        "sections under one whole-slide image"
                    ),
                    "conflict_flag": "false",
                },
                {
                    "evidence_id": f"{STCRC_STUDY}::{sample}::block",
                    "entity_type": "physical_unit",
                    "entity_id": unit_id,
                    "field": "block_id",
                    "source_asset": readme_asset,
                    "metadata_path": readme_rel,
                    "metadata_key": "README.md introduction",
                    "raw_value": (
                        "We processed fresh-frozen resection samples obtained "
                        "from seven CRC patients ... We considered two serial "
                        "sections per patient to generate technical replicates."
                    ),
                    "normalized_value": block_id,
                    "evidence_grade": "E2_corroborated",
                    "interpretation": (
                        "two serial sections per patient are technical "
                        "replicates of one tumor resection, so the patient "
                        "tumor is the block; section order within a block is "
                        "not asserted"
                    ),
                    "conflict_flag": "false",
                },
            ]
        )
    for patient, sections in per_patient.items():
        if sorted(sections) != ["Rep1", "Rep2"]:
            raise ValueError(
                f"ST_CRC_CMS: patient {patient} does not have exactly Rep1+Rep2 sections"
            )
    return assets, units, evidence, duplicates


def extract_usz(
    repo_root: Path,
) -> tuple[
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    base = repo_root / "data/other_sources/zenodo_usz_tls_visium"
    record_path = base / "zenodo_record_14620362.json"
    big_zip = base / "TLS_VISIUM_USZ.zip"
    if not record_path.is_file() or not big_zip.is_file():
        raise ValueError("TLS_VISIUM_USZ: Zenodo record JSON and deposit zip are required")
    description, record = _normalized_description(record_path)
    for quote in (USZ_TUMOR_QUOTE, USZ_SAMPLE_QUOTE):
        if quote not in description:
            raise ValueError(f"TLS_VISIUM_USZ: Zenodo record no longer states: {quote}")
    file_keys = {str(item.get("key", "")) for item in record.get("files", [])}
    if file_keys != {"TLS_VISIUM_USZ.zip"}:
        raise ValueError(f"TLS_VISIUM_USZ: unexpected record file set {sorted(file_keys)}")

    sidecar_dir = base / "TLS_VISIUM_USZ/h5ad_preprocessed"
    assets = [
        _asset(record_path, repo_root, USZ_STUDY, USZ_STUDY, "zenodo_record", "json"),
        _asset(
            big_zip,
            repo_root,
            USZ_STUDY,
            USZ_STUDY,
            "deposit_zip_header",
            "zip",
            header_only=True,
        ),
    ]
    units: list[dict[str, str]] = []
    evidence: list[dict[str, str]] = []
    record_rel = repo_relative(record_path, repo_root)
    record_asset = f"{USZ_STUDY}::{USZ_STUDY}::zenodo_record"
    for sample in USZ_SAMPLES:
        sidecar = sidecar_dir / f"{sample}.json"
        if not sidecar.is_file():
            raise ValueError(f"TLS_VISIUM_USZ: missing deposit sidecar {sidecar}")
        payload = json.loads(sidecar.read_text(encoding="utf-8"))
        if payload.get("SAMPLE") != sample:
            raise ValueError(f"TLS_VISIUM_USZ: sidecar {sidecar} SAMPLE mismatch")
        assets.append(
            _asset(sidecar, repo_root, USZ_STUDY, USZ_STUDY, f"sample_sidecar::{sample}", "json")
        )
        sidecar_rel = repo_relative(sidecar, repo_root)
        sidecar_asset = f"{USZ_STUDY}::{USZ_STUDY}::sample_sidecar::{sample}"
        unit_id = f"{USZ_STUDY}::{sample}"
        patient_id = f"{USZ_STUDY}::PATIENT::{_hash16(sample)}"
        block_id = f"{USZ_STUDY}::BLOCK::{_hash16(sample + '|tumor_section')}"
        units.append(
            {
                "physical_unit_id": unit_id,
                "source_namespace": USZ_STUDY,
                "source_record_id": f"TLS_VISIUM_USZ/h5ad_preprocessed/{sample}.h5ad",
                "study_id": USZ_STUDY,
                "patient_id": patient_id,
                "block_id": block_id,
                "section_id": sample,
                "z_position": "",
                "section_order": "",
                "section_thickness": "5 um",
                "platform": "10x_visium",
                "identity_status": "PATIENT_BLOCK_EXPLICIT",
                "evidence_grade": "E3_explicit",
                "record_status": "RESOLVED_INCLUDED_CANDIDATE",
                "physical_specimen_id": "",
                "identity_granularity": "",
                "block_equivalent_status": "",
                "block_equivalent_basis": "",
            }
        )
        evidence.extend(
            [
                {
                    "evidence_id": f"{USZ_STUDY}::{sample}::patient",
                    "entity_type": "physical_unit",
                    "entity_id": unit_id,
                    "field": "patient_id",
                    "source_asset": sidecar_asset,
                    "metadata_path": sidecar_rel,
                    "metadata_key": "SAMPLE",
                    "raw_value": sample,
                    "normalized_value": patient_id,
                    "evidence_grade": "E3_explicit",
                    "interpretation": (
                        "deposit sidecar identifies the sample; the record "
                        "states KC[1-3] are kidney tumors and LC[1-5] are lung "
                        "tumors with no inter-sample patient linkage, so "
                        "identity is asserted per tumor"
                    ),
                    "conflict_flag": "false",
                },
                {
                    "evidence_id": f"{USZ_STUDY}::{sample}::block",
                    "entity_type": "physical_unit",
                    "entity_id": unit_id,
                    "field": "block_id",
                    "source_asset": record_asset,
                    "metadata_path": record_rel,
                    "metadata_key": "metadata.description",
                    "raw_value": (
                        "consists of 5 μm thick FFPE sections from kidney (3) "
                        "and lung (5) tumors ... Data samples include KC[1-3] "
                        "indicating kidney cancer and LC[1-5] indicating lung "
                        "cancer"
                    ),
                    "normalized_value": block_id,
                    "evidence_grade": "E3_explicit",
                    "interpretation": (
                        "one 5 um FFPE Visium section per tumor; the tumor is "
                        "the block-equivalent unit of identity"
                    ),
                    "conflict_flag": "false",
                },
            ]
        )
    duplicates = [
        {
            "duplicate_group_id": f"SOURCE_LINEAGE::{USZ_STUDY}",
            "member_type": "study",
            "member_id": USZ_STUDY,
            "relation_type": "canonical_source_lineage",
            "evidence_grade": "E3_explicit",
            "resolution": "canonical_study",
            "leakage_group_id": USZ_LEAKAGE,
        }
    ]
    return assets, units, evidence, duplicates


def extract(repo_root: Path, output_dir: Path, lineages: Iterable[str]) -> None:
    source_assets: list[dict[str, str]] = []
    physical_units: list[dict[str, str]] = []
    evidence: list[dict[str, str]] = []
    duplicates: list[dict[str, str]] = []
    builders = {
        "GSE175540": extract_kirc,
        "ST_CRC_CMS": extract_stcrc,
        "TLS_VISIUM_USZ": extract_usz,
    }
    for lineage in lineages:
        if lineage not in builders:
            raise ValueError(f"unsupported validation lineage: {lineage}")
        assets, units, rows, groups = builders[lineage](repo_root)
        source_assets.extend(assets)
        physical_units.extend(units)
        evidence.extend(rows)
        duplicates.extend(groups)

    write_tsv(
        output_dir / "validation_lineages_source_assets.tsv",
        SOURCE_ASSET_FIELDS,
        source_assets,
    )
    write_tsv(
        output_dir / "validation_lineages_physical_units.tsv",
        PHYSICAL_UNIT_FIELDS,
        physical_units,
    )
    write_tsv(
        output_dir / "validation_lineages_identity_evidence.tsv",
        IDENTITY_EVIDENCE_FIELDS,
        evidence,
    )
    write_tsv(
        output_dir / "validation_lineages_duplicate_groups.tsv",
        DUPLICATE_GROUP_FIELDS,
        duplicates,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("infra/sample-registry/staging"),
    )
    parser.add_argument("--lineage", action="append", choices=LINEAGES)
    args = parser.parse_args()
    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir
    if not output_dir.is_absolute():
        output_dir = repo_root / output_dir
    try:
        extract(repo_root, output_dir, args.lineage or LINEAGES)
    except (OSError, ValueError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
