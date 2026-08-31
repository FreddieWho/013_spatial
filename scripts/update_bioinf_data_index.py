#!/usr/bin/env python3
"""Generate the provenance index for approved metadata-only external inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import tempfile
from datetime import datetime, timezone
from pathlib import Path


INDEX_FIELDS = [
    "source_id",
    "accession",
    "official_url",
    "local_path",
    "bytes",
    "sha256",
    "license",
    "identity_coverage",
    "retention_scope",
    "status",
    "notes",
]
FORBIDDEN_SUFFIXES = {
    ".h5",
    ".h5ad",
    ".jpg",
    ".jpeg",
    ".gif",
    ".png",
    ".tif",
    ".tiff",
    ".svs",
}
INCIDENT = {
    "source": "PMC11508537 supplementaryFiles",
    "package_sha256": "55f7ec7db7baae38b8c5dd8ca5917e4f41701241a05d864e777489ef42d3ba29",
    "metadata_only_package_sha256": (
        "72cdb27990ae78e7606ce00928c1e0feb11c297cd4c847b671c47ddfa43c9378"
    ),
    "deleted_result_file_sha256": [
        "62dc15dc772f8e3f42f374a47fe39305a44f9e49bb451ab789d66a4ac769403b",
        "e61bc1511dcaeeb6efdeb5476fbe8048b520f272f390e7258322eb251ebc0b46",
    ],
    "observed_members": {"jpg": 7, "gif": 7, "nested_zip": 1},
    "content_inspection": "filenames and archive member types only; images were not opened",
    "retained": False,
    "deleted": True,
}
FULL_TEXT_INCIDENTS = [
    {
        "source": "PMC article XML scope cleanup",
        "deleted_path": (
            f"infra/bioinf-data-index/raw/pmc/{accession}_article.xml"
        ),
        "sha256": None,
        "sha256_status": "not_computed_before_deletion",
        "retained": False,
        "deleted": True,
    }
    for accession in (
        "PMC9862087",
        "PMC10439131",
        "PMC10991508",
        "PMC11508537",
    )
]
ACCEPTED_GEO_ACCESSIONS = {
    "GSE211956",
    "GSE226997",
    "GSE274103",
    "GSE274557",
}
R04_REFERENCE_SOURCES = (
    {
        "source_id": "local:005_preCan:GSE132465_counts",
        "accession": "GSE132465",
        "local_path": "/home/huyudi/005_preCan/data/downloaded/CIT_CRC_005/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz",
        "identity_coverage": "donor_cell_type",
        "notes": "R04 CRC reference; read in place, capped and never copied into this repository",
    },
    {
        "source_id": "local:005_preCan:GSE132465_annotation",
        "accession": "GSE132465",
        "local_path": "/home/huyudi/005_preCan/data/downloaded/CIT_CRC_005/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz",
        "identity_coverage": "donor_cell_type",
        "notes": "R04 annotation; only Index/Patient/Cell_type are allowlisted",
    },
    {
        "source_id": "local:006:gse115978",
        "accession": "GSE115978",
        "local_path": "/home/huyudi/006/data/processed/srt/raw/gse115978.h5ad",
        "identity_coverage": "donor_cell_type",
        "notes": "R04 h5ad reference; counts and patient_id/anno_orig allowlist only",
    },
    {
        "source_id": "local:006:gse232240",
        "accession": "GSE232240",
        "local_path": "/home/huyudi/006/data/processed/srt/raw/gse232240.h5ad",
        "identity_coverage": "donor_cell_type",
        "notes": "R04 h5ad reference; X and patient_id/cell_type allowlist only",
    },
)
R04_SOURCE_METADATA = (
    {
        "source_id": "zenodo:7760264:record",
        "accession": "ZENODO:7760264",
        "local_path": "data/other_sources/zenodo_st_crc_cms/zenodo_record_7760264.json",
        "identity_coverage": "spatial_cohort_section",
        "notes": (
            "R04 ST-CRC-CMS source metadata; only the official record JSON is indexed, "
            "while pre-existing local matrices remain outside this metadata index"
        ),
    },
    {
        "source_id": "zenodo:14620362:record",
        "accession": "ZENODO:14620362",
        "local_path": "data/other_sources/zenodo_usz_tls_visium/zenodo_record_14620362.json",
        "identity_coverage": "spatial_cohort_section",
        "notes": (
            "R04 USZ TLS Visium source metadata; only the official record JSON is indexed, "
            "while pre-existing local matrices remain outside this metadata index"
        ),
    },
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def accession_for(path: Path) -> str:
    name = path.name
    for pattern in (r"GSE\d+", r"PMC\d+", r"PMID\d+"):
        import re

        match = re.search(pattern, name)
        if match:
            return match.group(0)
    return "MULTI"


def official_url(accession: str) -> str:
    if accession.startswith("GSE"):
        return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
    if accession.startswith("PMC"):
        return f"https://pmc.ncbi.nlm.nih.gov/articles/{accession}/"
    if accession.startswith("PMID"):
        return f"https://pubmed.ncbi.nlm.nih.gov/{accession.removeprefix('PMID')}/"
    if accession.startswith("ZENODO:"):
        return f"https://zenodo.org/records/{accession.removeprefix('ZENODO:')}"
    return "https://pubmed.ncbi.nlm.nih.gov/"


def license_for(accession: str, relative_raw: Path) -> str:
    if accession.startswith("PMC"):
        return "CC BY"
    if relative_raw.as_posix() == (
        "supplements/PMC10991508_supplementary_data_1.xlsx"
    ):
        return "CC BY"
    return "NOT_ASSERTED_IN_RETRIEVED_METADATA"


def identity_coverage_for(accession: str, relative_raw: Path) -> str:
    relative = relative_raw.as_posix()
    if relative == "supplements/PMC10991508_supplementary_data_1.xlsx":
        return "patient_to_gsm_crosswalk"
    if accession in ACCEPTED_GEO_ACCESSIONS and relative.endswith(
        "_series_quick.soft"
    ):
        return "cohort_tissue_patient_bioproject"
    if accession in ACCEPTED_GEO_ACCESSIONS and relative.endswith(
        "_gsm_quick.soft"
    ):
        return "patient_tissue_biosample"
    if accession.startswith("GSE"):
        return "screened_not_accepted"
    return "provenance_only"


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        handle.write(text)
        temporary = Path(handle.name)
    temporary.replace(path)


def build_index(
    repo_root: Path,
    index_root: Path,
    raw_root: Path,
    *,
    generated_at: str | None = None,
    include_r04_references: bool = False,
) -> list[dict[str, object]]:
    rows: list[dict[str, str]] = []
    for path in sorted(raw_root.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix.casefold() in FORBIDDEN_SUFFIXES:
            raise ValueError(f"forbidden expression/image file retained in metadata index: {path}")
        relative_raw = path.relative_to(raw_root)
        accession = accession_for(path)
        scope = (
            "identity_crosswalk"
            if relative_raw.as_posix()
            == "supplements/PMC10991508_supplementary_data_1.xlsx"
            else "metadata"
        )
        try:
            local_path = path.resolve().relative_to(repo_root.resolve()).as_posix()
        except ValueError:
            local_path = path.resolve().as_posix()
        rows.append(
            {
                "source_id": f"{relative_raw.parent.name}:{path.stem}",
                "accession": accession,
                "official_url": official_url(accession),
                "local_path": local_path,
                "bytes": str(path.stat().st_size),
                "sha256": sha256(path),
                "license": license_for(accession, relative_raw),
                "identity_coverage": identity_coverage_for(
                    accession,
                    relative_raw,
                ),
                "retention_scope": scope,
                "status": "retained",
                "notes": (
                    "official patient-to-GSM identity crosswalk; biological result "
                    "columns were not used"
                    if scope == "identity_crosswalk"
                    else "official metadata/document record"
                ),
            }
        )
    if include_r04_references:
        for source in R04_SOURCE_METADATA:
            path = repo_root / source["local_path"]
            if not path.is_file():
                raise FileNotFoundError(f"R04 source metadata is unavailable: {path}")
            rows.append(
                {
                    "source_id": source["source_id"],
                    "accession": source["accession"],
                    "official_url": official_url(source["accession"]),
                    "local_path": source["local_path"],
                    "bytes": str(path.stat().st_size),
                    "sha256": sha256(path),
                    "license": "NOT_ASSERTED_IN_RETRIEVED_METADATA",
                    "identity_coverage": source["identity_coverage"],
                    "retention_scope": "r04_source_metadata_pointer",
                    "status": "retained_local_metadata",
                    "notes": source["notes"],
                }
            )
        for source in R04_REFERENCE_SOURCES:
            path = Path(source["local_path"])
            if not path.is_file():
                raise FileNotFoundError(f"R04 reference pointer is unavailable: {path}")
            rows.append(
                {
                    "source_id": source["source_id"],
                    "accession": source["accession"],
                    "official_url": official_url(source["accession"]),
                    "local_path": str(path),
                    "bytes": str(path.stat().st_size),
                    "sha256": sha256(path),
                    "license": "LOCAL_USER_APPROVED_SOURCE",
                    "identity_coverage": source["identity_coverage"],
                    "retention_scope": "approved_local_reference_pointer",
                    "status": "not_copied_external",
                    "notes": source["notes"],
                }
            )

    index_root.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", newline="", encoding="utf-8", dir=index_root, delete=False
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=INDEX_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    temporary.replace(index_root / "index.tsv")

    generated_at = generated_at or datetime.now(timezone.utc).isoformat()
    manifest: dict[str, object] = {
        "schema_version": 3 if include_r04_references else 1,
        "generated_at": generated_at,
        "policy": {
            "allowed": [
                "official metadata",
                "official identity crosswalk",
                *( ["R-04 official source metadata"] if include_r04_references else []),
                *( ["approved local reference pointer read without repository copy"] if include_r04_references else []),
            ],
            "forbidden": ["expression matrix", "image", "biological result table"],
            "image_content_opened": False,
        },
        "sources": rows,
        "scope_contamination_incidents": [INCIDENT, *FULL_TEXT_INCIDENTS],
    }
    summary: dict[str, object] = {
        "schema_version": 3 if include_r04_references else 1,
        "generated_at": generated_at,
        "indexed_file_count": len(rows),
        "indexed_bytes": sum(int(row["bytes"]) for row in rows),
        "retained_metadata_files": sum(row["retention_scope"] == "metadata" for row in rows),
        "retained_identity_crosswalk_files": sum(
            row["retention_scope"] == "identity_crosswalk" for row in rows
        ),
        "retained_expression_files": 0,
        "retained_image_files": 0,
        "retained_result_files": 0,
        "external_reference_pointer_count": sum(
            row["retention_scope"] == "approved_local_reference_pointer" for row in rows
        ),
        "external_reference_pointer_bytes": sum(
            int(row["bytes"])
            for row in rows
            if row["retention_scope"] == "approved_local_reference_pointer"
        ),
        "r04_source_metadata_pointer_count": sum(
            row["retention_scope"] == "r04_source_metadata_pointer" for row in rows
        ),
        "r04_source_metadata_pointer_bytes": sum(
            int(row["bytes"])
            for row in rows
            if row["retention_scope"] == "r04_source_metadata_pointer"
        ),
        "deleted_scope_contamination_incidents": 1 + len(
            FULL_TEXT_INCIDENTS
        ),
    }
    atomic_text(index_root / "manifest.json", json.dumps(manifest, indent=2) + "\n")
    atomic_text(index_root / "summary.json", json.dumps(summary, indent=2) + "\n")
    atomic_text(
        index_root / "README.md",
        """# External bioinformatics metadata index

This directory records the approved metadata-only external inputs used for
R-01 identity resolution. `index.tsv` is the checksum inventory,
`manifest.json` records provenance and scope incidents, and `summary.json`
provides machine-readable retention totals.

Only official metadata and an official patient-to-GSM identity crosswalk are
retained here. When explicitly requested, R-04 also records official Zenodo
source metadata JSON and hashed pointers to existing local reference files;
the latter are read in place with an allowlist and are never copied into this
repository. No image, expression matrix or biological result file is retained.
Local GEO raw archives outside this directory are used only as header-level
GSM locators.
""",
    )
    return list(manifest["sources"])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--index-root", type=Path)
    parser.add_argument("--raw-root", type=Path)
    parser.add_argument(
        "--generated-at",
        help="Stable ISO-8601 acquisition timestamp for versioned rebuilds.",
    )
    parser.add_argument(
        "--include-r04-references",
        action="store_true",
        help="Record hashed pointers to the approved local 005/006 reference files.",
    )
    args = parser.parse_args()
    repo_root = args.repo_root.resolve()
    index_root = args.index_root or repo_root / "infra" / "bioinf-data-index"
    raw_root = args.raw_root or index_root / "raw"
    try:
        build_index(
            repo_root,
            index_root.resolve(),
            raw_root.resolve(),
            generated_at=args.generated_at,
            include_r04_references=args.include_r04_references,
        )
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
