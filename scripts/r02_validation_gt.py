#!/usr/bin/env python3
"""Shared deterministic validation-lineage GT logic for R-02 (D-037).

Three public lineages carry human-annotated ground truth (D-036: any human
annotator is high-confidence at this stage):

- ``GEO::GSE175540`` (Meylan et al. 2022, Immunity, PMID 35231421): per-spot
  TLS annotations deposited at GEO (``*_TLS_annotation.csv.gz``), replayed to
  array coordinates through the same deposit's ``tissue_positions_list`` CSVs.
  Files with the ``Barcode,TLS`` header carry only ``T_agg`` (T-cell aggregate)
  labels; T_agg is not TLS and those files are fail-closed NOT_AUDITABLE rows.
  Empty annotation cells are unknown, never negative.
- ``ST_CRC_CMS`` (Valdeolivas et al. 2024, npj Precis Oncol, Zenodo 7760264):
  pathologist spot categorization CSVs (``Barcode,Pathologist*``), replayed
  through ``tissue_positions_list.csv`` members extracted verbatim from the
  deposit sample zips. Only the explicit ``tumor&stroma*`` label family is
  scoped as TUMOR_STROMA_BOUNDARY; every other label (including all
  misspelling variants) is frozen deposit vocabulary but unscoped. ``IC
  aggregate*`` labels are immune-cell aggregates, not verified TLS, and are
  fail-closed excluded from the TLS scope.
- ``TLS_VISIUM_USZ`` (Zenodo 14620362): expert-annotated labels stored inside
  the deposit h5ad ``obs/ground_truth`` categorical column, replayed through
  ``obs/_index`` + ``obs/x_array``/``obs/y_array`` from the same files (D-038
  authorizes these per-deposit label/coordinate column reads; no expression or
  image data is read anywhere).

Instances are deterministic hex-grid connected components of the scoped spots,
the same construction used for the Heiser training lineage. No expression
values and no image pixels are read anywhere in this module. The verifier
never trusts registry TSV self-reports; every check recomputes from the raw
files.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import re
from functools import lru_cache
from pathlib import Path

import h5py

try:
    from scripts.r02_heiser_gt import connected_components, spot_index_fingerprint
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from r02_heiser_gt import connected_components, spot_index_fingerprint


VERIFIER_KIRC = "SPOT_BARCODE_TLS_ANNOTATION_CSV"
VERIFIER_USZ = "H5AD_OBS_GROUND_TRUTH_LABELS"
VERIFIER_STCRC = "SPOT_BARCODE_PATHOLOGY_CATEGORY_CSV"
VALIDATION_VERIFIER_CLASSES = (VERIFIER_KIRC, VERIFIER_USZ, VERIFIER_STCRC)

OVERLAP_LOCATOR = (
    "infra/structure-registry/input_policy.tsv#he-annotation-gt-input-exclusion"
)

# --- GSE175540 (KIRC) -------------------------------------------------------
KIRC_RAW_DIR = Path("data/GEO/GSE175540/raw")
KIRC_LOGICAL_UNIT = "GEO::GSE175540"
KIRC_ANNOTATION_RE = re.compile(r"^(GSM\d+)_((?:ffpe|frozen)_[a-z]_\d+)_TLS_annotation\.csv\.gz$")
KIRC_EXPECTED_ANNOTATION_FILES = 23
KIRC_TLS_CAT_VALUES = {"TLS", "NO_TLS", ""}
KIRC_TAGG_VALUES = {"T_agg", ""}
KIRC_PROVENANCE_NOTE = (
    "infra/structure-registry/provenance/MEYLAN_2022_IMMUNITY_PMID35231421.md"
)
KIRC_GT_DEFINITION_MODALITY = (
    "author TLS spot annotation guided by CD3/CD20 immunostaining "
    "(Meylan et al. 2022, Immunity, PMID 35231421), deposited at GEO"
)
KIRC_BOUNDARY_UNCERTAINTY = (
    "spot-resolution (55 um) annotation; TLS instance = hex-connected "
    "component of TLS-labeled spots; atlas Table S4 histological structure "
    "counts differ by definition granularity and are not used as instance GT"
)

# --- TLS_VISIUM_USZ ---------------------------------------------------------
USZ_H5AD_DIR = Path("data/other_sources/zenodo_usz_tls_visium/TLS_VISIUM_USZ/h5ad_preprocessed")
USZ_LOGICAL_UNIT = "TLS_VISIUM_USZ"
USZ_SAMPLES = ("KC1", "KC2", "KC3", "LC1", "LC2", "LC3", "LC4", "LC5")
USZ_LABEL_VOCABULARY = {"TLS", "INFL", "TUM", "NOR", "UNASSIGNED", "LN"}
USZ_PROVENANCE_NOTE = (
    "infra/structure-registry/provenance/USZ_ZENODO_14620362_TLS_VISIUM.md"
)
USZ_GT_DEFINITION_MODALITY = (
    "manual expert annotation of Visium spots on H&E images "
    "(Zenodo 14620362; annotators K.S. and S.D., expert researchers)"
)
USZ_BOUNDARY_UNCERTAINTY = (
    "spot-resolution (55 um) manual annotation; TLS instance = hex-connected "
    "component of ground_truth=TLS spots; UNASSIGNED cells are unknown, never "
    "negative"
)

# --- ST_CRC_CMS (Valdeolivas) ----------------------------------------------
STCRC_DIR = Path("data/other_sources/zenodo_st_crc_cms")
STCRC_LOGICAL_UNIT = "ST_CRC_CMS"
STCRC_EXPECTED_ANNOTATION_FILES = 14
STCRC_ANNOTATION_RE = re.compile(r"^Pathologist_Annotations_(SN\d+_A\d+_Rep[12](?:_X)?)\.csv$")
STCRC_HEADERS = {
    "Pathologist annotation",
    "Pathologist Annotation",
    "Pathologist Annotations",
    "Pathologist_KH",
}
# Verbatim frozen deposit vocabulary (misspellings included, fail-closed).
STCRC_LABEL_VOCABULARY = {
    "",
    "IC aggragate_stroma or muscularis",
    "IC aggreagate_connective tissue",
    "IC aggregate connective tissue",
    "IC aggregate submucosa",
    "IC aggregate_muscularis or stroma",
    "IC aggregate_stroma or muscularis",
    "IC aggregate_submucosa",
    "IC aggregregate_submucosa",
    "IC aggreates_stroma or muscularis",
    "connective tissue_1_edema",
    "connective tissue_2_fibroblastic_IC low",
    "connective tissue_3_fibroblastic_IC med",
    "connective tissue_4_muscularis_IC low",
    "connective tissue_6_hemosiderin?",
    "epithelium&lam propria",
    "epithelium&submucosa",
    "epitehlium&submucosa",
    "exclude",
    "glandular tissue",
    "lamina propria",
    "muscularis_IC med to high",
    "non neo epithelium",
    "squamous epithelium",
    "stroma desmoplastic_IC low",
    "stroma desmoplastic_IC med to high",
    "stroma_desmoplastic_IC low",
    "stroma_desmoplastic_IC med to high",
    "stroma_fibroblastic_IC high",
    "stroma_fibroblastic_IC low",
    "stroma_fibroblastic_IC med",
    "stroma_fibroblastic_IC_high",
    "stroma_fibroblastic_IC_med",
    "submucosa",
    "tumor",
    "tumor&stroma",
    "tumor&stroma IC med to high",
    "tumor&stroma_IC low",
    "tumor&stroma_IC med to high",
}
STCRC_TSB_LABELS = (
    "tumor&stroma",
    "tumor&stroma IC med to high",
    "tumor&stroma_IC low",
    "tumor&stroma_IC med to high",
)
STCRC_PROVENANCE_NOTE = (
    "infra/structure-registry/provenance/"
    "VALDEOLIVAS_2024_NPJ_PRECIS_ONCOL_DOI_10.1038_s41698-023-00488-4.md"
)
STCRC_GT_DEFINITION_MODALITY = (
    "pathologist spot categorization on H&E (QuPath/Loupe; Valdeolivas et "
    "al. 2024, npj Precis Oncol, doi:10.1038/s41698-023-00488-4), deposited "
    "at Zenodo 7760264"
)
STCRC_BOUNDARY_UNCERTAINTY = (
    "spot-resolution (55 um) pathologist categorization; boundary instance = "
    "hex-connected component of tumor&stroma-family spots; per-file annotator "
    "identity is attributed only for A938797 (Pathologist_KH header)"
)

BARCODE_RE = re.compile(r"^[ACGT]{16}-1$")

ROLE_USAGE = {
    "external_validation": (
        "locked external validation within frozen outer splits; "
        "molecular-only input budgets",
        "selection, tuning, threshold choice, training, or discovery-driven "
        "replacement; as model input feature or label; any use with H&E- or "
        "image-derived inputs",
    ),
    "internal_validation": (
        "predeclared internal validation within frozen outer splits; "
        "molecular-only input budgets",
        "training, tuning, external confirmation, or role replacement; as "
        "model input feature or label; any use with H&E- or image-derived "
        "inputs",
    ),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _hash16(raw: str) -> str:
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:16]


def _read_tissue_positions_csv(rows: csv.reader) -> dict[str, tuple[int, int, int]]:
    """barcode -> (in_tissue, array_row, array_col); optional header row."""
    positions: dict[str, tuple[int, int, int]] = {}
    for row_number, row in enumerate(rows, start=1):
        if row_number == 1 and row and row[0].strip().lower() == "barcode":
            continue
        if len(row) < 4:
            raise ValueError(f"malformed tissue-positions row {row_number}")
        barcode = row[0].strip()
        if not BARCODE_RE.match(barcode):
            raise ValueError(f"invalid tissue-positions barcode {barcode!r} row {row_number}")
        if barcode in positions:
            raise ValueError(f"duplicate tissue-positions barcode {barcode!r}")
        positions[barcode] = (int(row[1]), int(row[2]), int(row[3]))
    return positions


def _spot_geometry(
    positions: dict[str, tuple[int, int, int]]
) -> dict[str, tuple[int, int]]:
    """in-tissue barcode -> (array_row, array_col)."""
    return {
        barcode: (row, col)
        for barcode, (in_tissue, row, col) in positions.items()
        if in_tissue == 1
    }


def _role_usage(logical_unit: str) -> tuple[str, str]:
    role = {
        KIRC_LOGICAL_UNIT: "external_validation",
        USZ_LOGICAL_UNIT: "external_validation",
        STCRC_LOGICAL_UNIT: "internal_validation",
    }[logical_unit]
    return ROLE_USAGE[role]


def _audit_row(
    *,
    source_id: str,
    source_class: str,
    structure_id: str,
    path: Path,
    root: Path,
    schema_locator: str,
    physical_unit_id: str,
    geometry_locator: str,
    provenance_note: str,
    audit_status: str,
    notes: str,
    physical_link_status: str = "VERIFIED",
) -> dict[str, str]:
    return {
        "gt_source_id": source_id,
        "source_class": source_class,
        "structure_id": structure_id,
        "path": path.relative_to(root).as_posix(),
        "sha256": _sha256(path),
        "checksum_status": "VERIFIED",
        "schema_locator": schema_locator,
        "physical_unit_id": physical_unit_id,
        "physical_link_status": physical_link_status,
        "geometry_locator": geometry_locator,
        "provenance_status": "KNOWN",
        "provenance_locator": provenance_note,
        "overlap_status": "KNOWN_OVERLAP",
        "overlap_locator": OVERLAP_LOCATOR,
        "same_assay_status": "INDEPENDENT_ASSAY",
        "audit_status": audit_status,
        "notes": notes,
    }


def _instance_row(
    *,
    instance_id: str,
    structure_id: str,
    logical_unit: str,
    r01_row: dict[str, str],
    source_id: str,
    modality: str,
    geometry_locator: str,
    boundary_uncertainty: str,
) -> dict[str, str]:
    allowed, forbidden = _role_usage(logical_unit)
    return {
        "instance_id": instance_id,
        "structure_id": structure_id,
        "logical_unit_id": logical_unit,
        "physical_unit_id": r01_row["physical_unit_id"],
        "patient_id": r01_row["patient_id"],
        "block_id": r01_row["block_id"],
        "physical_specimen_id": r01_row["physical_specimen_id"],
        "gt_source_id": source_id,
        "gt_definition_modality": modality,
        "gt_geometry_locator": geometry_locator,
        "boundary_uncertainty": boundary_uncertainty,
        "confirmation_status": "CONFIRMATORY",
        "allowed_use": allowed,
        "forbidden_use": forbidden,
    }


def _require_r01_unit(
    r01_units: dict[str, dict[str, str]], physical_id: str, lineage: str
) -> dict[str, str]:
    row = r01_units.get(physical_id)
    if row is None:
        raise ValueError(f"{lineage}: no R-01 physical unit for {physical_id}")
    if row["record_status"] != "RESOLVED_INCLUDED_CANDIDATE":
        raise ValueError(
            f"{lineage}: R-01 unit {physical_id} status "
            f"{row['record_status']} is not eligible"
        )
    return row


def _replay_row(
    physical_id: str, path: Path, root: Path, geometry: dict[str, tuple[int, int]]
) -> dict[str, str]:
    return {
        "physical_unit_id": physical_id,
        "path": path.relative_to(root).as_posix(),
        "spot_count": str(len(geometry)),
        "index_fingerprint_sha256": spot_index_fingerprint(geometry),
        "fingerprint_status": "INDEX_VALUE_PINNED",
    }


# ---------------------------------------------------------------------------
# GSE175540 (KIRC)


def kirc_annotation_files(root: Path) -> list[tuple[str, str, Path]]:
    paths = sorted((root / KIRC_RAW_DIR).glob("*_TLS_annotation.csv.gz"))
    if len(paths) != KIRC_EXPECTED_ANNOTATION_FILES:
        raise ValueError(
            f"GSE175540: expected {KIRC_EXPECTED_ANNOTATION_FILES} TLS "
            f"annotation files, observed {len(paths)}"
        )
    found: list[tuple[str, str, Path]] = []
    for path in paths:
        match = KIRC_ANNOTATION_RE.match(path.name)
        if not match:
            raise ValueError(f"GSE175540: unexpected annotation file name {path.name}")
        found.append((match.group(1), match.group(2), path))
    return found


def read_kirc_annotation(path: Path) -> tuple[str, dict[str, list[str]]]:
    """(header_kind, label -> barcodes). Empty cells are unknown, never negative."""
    labels: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if header == ["Barcode", "TLS_2_cat"]:
            kind, vocabulary = "TLS_2_cat", KIRC_TLS_CAT_VALUES
        elif header == ["Barcode", "TLS"]:
            kind, vocabulary = "TLS", KIRC_TAGG_VALUES
        else:
            raise ValueError(f"GSE175540: unexpected annotation header in {path}: {header}")
        seen: set[str] = set()
        for row_number, row in enumerate(reader, start=2):
            if len(row) != 2:
                raise ValueError(f"GSE175540: malformed row {row_number} in {path}")
            barcode, label = row[0].strip(), row[1].strip()
            if not BARCODE_RE.match(barcode):
                raise ValueError(f"GSE175540: invalid barcode {barcode!r} in {path}:{row_number}")
            if barcode in seen:
                raise ValueError(f"GSE175540: duplicate barcode {barcode!r} in {path}:{row_number}")
            seen.add(barcode)
            if label not in vocabulary:
                raise ValueError(f"GSE175540: unknown label {label!r} in {path}:{row_number}")
            if label:
                labels.setdefault(label, []).append(barcode)
    if not seen:
        raise ValueError(f"GSE175540: empty annotation file {path}")
    return kind, labels


@lru_cache(maxsize=None)
def load_kirc_spot_geometry(root_str: str, gsm: str, alias: str) -> dict[str, tuple[int, int]]:
    path = Path(root_str) / KIRC_RAW_DIR / f"{gsm}_{alias}_tissue_positions_list.csv.gz"
    if not path.is_file():
        raise FileNotFoundError(f"GSE175540: missing tissue positions {path}")
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        positions = _read_tissue_positions_csv(csv.reader(handle))
    return _spot_geometry(positions)


def _kirc_records(
    root: Path, r01_units: dict[str, dict[str, str]]
) -> dict[str, list[dict[str, str]]]:
    gt_sources: list[dict[str, str]] = []
    instances: list[dict[str, str]] = []
    replay_index: list[dict[str, str]] = []
    for gsm, alias, path in kirc_annotation_files(root):
        physical_id = f"external_geo::{gsm}"
        r01_row = _require_r01_unit(r01_units, physical_id, "GSE175540")
        if r01_row["patient_id"] != f"GSE175540::{alias}":
            raise ValueError(
                f"GSE175540: annotation alias {alias} conflicts with R-01 "
                f"patient {r01_row['patient_id']} for {gsm}"
            )
        geometry = load_kirc_spot_geometry(str(root), gsm, alias)
        replay_index.append(_replay_row(physical_id, path.parent / path.name.replace(
            "_TLS_annotation.csv.gz", "_tissue_positions_list.csv.gz"), root, geometry))
        kind, labels = read_kirc_annotation(path)
        relative = path.relative_to(root).as_posix()
        if kind == "TLS":
            gt_sources.append(
                _audit_row(
                    source_id=f"GSE175540_TLS_ANNOTATION::{alias}::T_AGG_ONLY",
                    source_class=VERIFIER_KIRC,
                    structure_id="",
                    path=path,
                    root=root,
                    schema_locator="gzip-csv header:Barcode,TLS",
                    physical_unit_id=physical_id,
                    geometry_locator="",
                    provenance_note=KIRC_PROVENANCE_NOTE,
                    audit_status="NOT_AUDITABLE_GT",
                    notes=(
                        "file carries only T_agg (T-cell aggregate) labels; "
                        "T_agg is not TLS and is fail-closed excluded from the "
                        "scoped TLS ground truth"
                    ),
                )
            )
            continue
        tls_barcodes = [b for b in labels.get("TLS", []) if b in geometry]
        if not tls_barcodes:
            # Audited negative (all-NO_TLS or empty-of-TLS file): no scoped
            # evidence, no row (same convention as the Heiser lineage).
            continue
        locator = f"{relative}#TLS_2_cat=TLS"
        source_id = f"GSE175540_TLS_ANNOTATION::{alias}"
        gt_sources.append(
            _audit_row(
                source_id=source_id,
                source_class=VERIFIER_KIRC,
                structure_id="TLS",
                path=path,
                root=root,
                schema_locator="gzip-csv header:Barcode,TLS_2_cat",
                physical_unit_id=physical_id,
                geometry_locator=locator,
                provenance_note=KIRC_PROVENANCE_NOTE,
                audit_status="AUDITABLE_GT",
                notes=(
                    f"{len(tls_barcodes)} in-tissue TLS-labeled spots replayed "
                    "to array coordinates; empty annotation cells are unknown, "
                    "never negative"
                ),
            )
        )
        components = connected_components({geometry[b] for b in tls_barcodes})
        for index, component in enumerate(components, start=1):
            instances.append(
                _instance_row(
                    instance_id=f"TLS::GSE175540::{alias}::comp{index:02d}",
                    structure_id="TLS",
                    logical_unit=KIRC_LOGICAL_UNIT,
                    r01_row=r01_row,
                    source_id=source_id,
                    modality=KIRC_GT_DEFINITION_MODALITY,
                    geometry_locator=locator,
                    boundary_uncertainty=KIRC_BOUNDARY_UNCERTAINTY,
                )
            )
    return {"gt_sources": gt_sources, "instances": instances, "replay_index": replay_index}


# ---------------------------------------------------------------------------
# TLS_VISIUM_USZ


@lru_cache(maxsize=None)
def load_usz_obs(root_str: str, sample: str) -> dict[str, tuple[int, int, str]]:
    """barcode -> (x_array, y_array, ground_truth label). D-038 scoped reads."""
    path = Path(root_str) / USZ_H5AD_DIR / f"{sample}.h5ad"
    if not path.is_file():
        raise FileNotFoundError(f"TLS_VISIUM_USZ: missing h5ad {path}")
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        raw_index = obs["_index"][:]
        x_array = obs["x_array"][:]
        y_array = obs["y_array"][:]
        group = obs["ground_truth"]
        categories = [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in group["categories"][:]
        ]
        codes = group["codes"][:]
    result: dict[str, tuple[int, int, str]] = {}
    for raw, x, y, code in zip(raw_index, x_array, y_array, codes):
        barcode = raw.decode() if isinstance(raw, bytes) else str(raw)
        if code < 0:
            raise ValueError(f"TLS_VISIUM_USZ: missing ground_truth code in {path}")
        label = categories[int(code)]
        if label not in USZ_LABEL_VOCABULARY:
            raise ValueError(f"TLS_VISIUM_USZ: unknown label {label!r} in {path}")
        if barcode in result:
            raise ValueError(f"TLS_VISIUM_USZ: duplicate barcode {barcode!r} in {path}")
        result[barcode] = (int(x), int(y), label)
    return result


def _usz_records(
    root: Path, r01_units: dict[str, dict[str, str]]
) -> dict[str, list[dict[str, str]]]:
    gt_sources: list[dict[str, str]] = []
    instances: list[dict[str, str]] = []
    replay_index: list[dict[str, str]] = []
    for sample in USZ_SAMPLES:
        physical_id = f"{USZ_LOGICAL_UNIT}::{sample}"
        r01_row = _require_r01_unit(r01_units, physical_id, USZ_LOGICAL_UNIT)
        path = root / USZ_H5AD_DIR / f"{sample}.h5ad"
        if r01_row["source_record_id"] != (
            f"TLS_VISIUM_USZ/h5ad_preprocessed/{sample}.h5ad"
        ):
            raise ValueError(
                f"TLS_VISIUM_USZ: R-01 source record drifted for {sample}: "
                f"{r01_row['source_record_id']}"
            )
        obs = load_usz_obs(str(root), sample)
        geometry = {barcode: (x, y) for barcode, (x, y, _label) in obs.items()}
        replay_index.append(_replay_row(physical_id, path, root, geometry))
        tls_barcodes = [b for b, (_x, _y, label) in obs.items() if label == "TLS"]
        if not tls_barcodes:
            raise ValueError(f"TLS_VISIUM_USZ: no TLS spots in {sample}; scope review needed")
        relative = path.relative_to(root).as_posix()
        locator = f"{relative}#obs/ground_truth=TLS"
        source_id = f"USZ_H5AD_GROUND_TRUTH::{sample}"
        gt_sources.append(
            _audit_row(
                source_id=source_id,
                source_class=VERIFIER_USZ,
                structure_id="TLS",
                path=path,
                root=root,
                schema_locator="h5ad obs: _index + x_array/y_array + ground_truth categorical",
                physical_unit_id=physical_id,
                geometry_locator=locator,
                provenance_note=USZ_PROVENANCE_NOTE,
                audit_status="AUDITABLE_GT",
                notes=(
                    f"{len(tls_barcodes)} ground_truth=TLS spots replayed from "
                    "the deposit h5ad obs; UNASSIGNED cells are unknown, never "
                    "negative"
                ),
            )
        )
        components = connected_components({geometry[b] for b in tls_barcodes})
        for index, component in enumerate(components, start=1):
            instances.append(
                _instance_row(
                    instance_id=f"TLS::TLS_VISIUM_USZ::{sample}::comp{index:02d}",
                    structure_id="TLS",
                    logical_unit=USZ_LOGICAL_UNIT,
                    r01_row=r01_row,
                    source_id=source_id,
                    modality=USZ_GT_DEFINITION_MODALITY,
                    geometry_locator=locator,
                    boundary_uncertainty=USZ_BOUNDARY_UNCERTAINTY,
                )
            )
    return {"gt_sources": gt_sources, "instances": instances, "replay_index": replay_index}


# ---------------------------------------------------------------------------
# ST_CRC_CMS (Valdeolivas)


def stcrc_annotation_files(root: Path) -> list[tuple[str, Path]]:
    paths = sorted((root / STCRC_DIR / "Pathology_SpotAnnotations").glob("*.csv"))
    if len(paths) != STCRC_EXPECTED_ANNOTATION_FILES:
        raise ValueError(
            f"ST_CRC_CMS: expected {STCRC_EXPECTED_ANNOTATION_FILES} pathology "
            f"annotation CSVs, observed {len(paths)}"
        )
    found: list[tuple[str, Path]] = []
    for path in paths:
        match = STCRC_ANNOTATION_RE.match(path.name)
        if not match:
            raise ValueError(f"ST_CRC_CMS: unexpected annotation file name {path.name}")
        found.append((match.group(1), path))
    return found


def read_stcrc_annotation(path: Path) -> dict[str, list[str]]:
    """label -> barcodes. Empty annotation cells are unknown, never negative."""
    labels: dict[str, list[str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if not header or len(header) != 2 or header[0] != "Barcode" or header[1] not in STCRC_HEADERS:
            raise ValueError(f"ST_CRC_CMS: unexpected annotation header in {path}: {header}")
        seen: set[str] = set()
        for row_number, row in enumerate(reader, start=2):
            if len(row) != 2:
                raise ValueError(f"ST_CRC_CMS: malformed row {row_number} in {path}")
            barcode, label = row[0].strip(), row[1].strip()
            if not BARCODE_RE.match(barcode):
                raise ValueError(f"ST_CRC_CMS: invalid barcode {barcode!r} in {path}:{row_number}")
            if barcode in seen:
                raise ValueError(f"ST_CRC_CMS: duplicate barcode {barcode!r} in {path}:{row_number}")
            seen.add(barcode)
            if label not in STCRC_LABEL_VOCABULARY:
                raise ValueError(f"ST_CRC_CMS: unknown label {label!r} in {path}:{row_number}")
            if label:
                labels.setdefault(label, []).append(barcode)
    if not seen:
        raise ValueError(f"ST_CRC_CMS: empty annotation file {path}")
    return labels


@lru_cache(maxsize=None)
def load_stcrc_spot_geometry(root_str: str, sample: str) -> dict[str, tuple[int, int]]:
    path = (
        Path(root_str)
        / STCRC_DIR
        / "tissue_positions"
        / f"{sample}_tissue_positions_list.csv"
    )
    if not path.is_file():
        raise FileNotFoundError(f"ST_CRC_CMS: missing tissue positions {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        positions = _read_tissue_positions_csv(csv.reader(handle))
    return _spot_geometry(positions)


def _stcrc_records(
    root: Path, r01_units: dict[str, dict[str, str]]
) -> dict[str, list[dict[str, str]]]:
    gt_sources: list[dict[str, str]] = []
    instances: list[dict[str, str]] = []
    replay_index: list[dict[str, str]] = []
    locator_labels = ",".join(STCRC_TSB_LABELS)
    for sample, path in stcrc_annotation_files(root):
        physical_id = f"{STCRC_LOGICAL_UNIT}::{sample}"
        r01_row = _require_r01_unit(r01_units, physical_id, STCRC_LOGICAL_UNIT)
        stem = sample[: -len("_X")] if sample.endswith("_X") else sample
        patient_raw = stem.split("_")[1]
        if r01_row["patient_id"] != f"{STCRC_LOGICAL_UNIT}::PATIENT::{_hash16(patient_raw)}":
            raise ValueError(
                f"ST_CRC_CMS: annotation patient {patient_raw} conflicts with "
                f"R-01 patient {r01_row['patient_id']} for {sample}"
            )
        if r01_row["block_id"] != (
            f"{STCRC_LOGICAL_UNIT}::BLOCK::{_hash16(patient_raw + '|tumor_block')}"
        ):
            raise ValueError(f"ST_CRC_CMS: R-01 block drifted for {sample}")
        geometry = load_stcrc_spot_geometry(str(root), sample)
        replay_index.append(
            _replay_row(
                physical_id,
                root / STCRC_DIR / "tissue_positions" / f"{sample}_tissue_positions_list.csv",
                root,
                geometry,
            )
        )
        labels = read_stcrc_annotation(path)
        tsb_barcodes = [
            barcode
            for label in STCRC_TSB_LABELS
            for barcode in labels.get(label, [])
            if barcode in geometry
        ]
        if not tsb_barcodes:
            # Audited negative for the TSB family (e.g. A798015 Rep1/Rep2): no
            # scoped evidence, no row (Heiser-lineage convention).
            continue
        relative = path.relative_to(root).as_posix()
        locator = f"{relative}#tsb_label_family={locator_labels}"
        source_id = f"STCRC_PATHOLOGY::{sample}"
        gt_sources.append(
            _audit_row(
                source_id=source_id,
                source_class=VERIFIER_STCRC,
                structure_id="TUMOR_STROMA_BOUNDARY",
                path=path,
                root=root,
                schema_locator="csv header:Barcode,Pathologist*",
                physical_unit_id=physical_id,
                geometry_locator=locator,
                provenance_note=STCRC_PROVENANCE_NOTE,
                audit_status="AUDITABLE_GT",
                notes=(
                    f"{len(tsb_barcodes)} in-tissue tumor&stroma-family spots "
                    "replayed to array coordinates; empty annotation cells are "
                    "unknown, never negative"
                ),
            )
        )
        components = connected_components({geometry[b] for b in tsb_barcodes})
        for index, component in enumerate(components, start=1):
            instances.append(
                _instance_row(
                    instance_id=f"TSB::ST_CRC_CMS::{sample}::comp{index:02d}",
                    structure_id="TUMOR_STROMA_BOUNDARY",
                    logical_unit=STCRC_LOGICAL_UNIT,
                    r01_row=r01_row,
                    source_id=source_id,
                    modality=STCRC_GT_DEFINITION_MODALITY,
                    geometry_locator=locator,
                    boundary_uncertainty=STCRC_BOUNDARY_UNCERTAINTY,
                )
            )
    return {"gt_sources": gt_sources, "instances": instances, "replay_index": replay_index}


# ---------------------------------------------------------------------------


def build_validation_records(
    root: Path, physical_rows: list[dict[str, str]]
) -> dict[str, list[dict[str, str]]]:
    """Derive validation-lineage gt_source, instance, and replay-index rows."""
    root = Path(root).resolve()
    r01_units = {
        row["physical_unit_id"]: row
        for row in physical_rows
        if row["study_id"] in {KIRC_LOGICAL_UNIT, USZ_LOGICAL_UNIT, STCRC_LOGICAL_UNIT}
    }
    kirc = _kirc_records(root, r01_units)
    usz = _usz_records(root, r01_units)
    stcrc = _stcrc_records(root, r01_units)
    gt_sources = kirc["gt_sources"] + usz["gt_sources"] + stcrc["gt_sources"]
    instances = kirc["instances"] + usz["instances"] + stcrc["instances"]
    replay = kirc["replay_index"] + usz["replay_index"] + stcrc["replay_index"]
    gt_sources.sort(key=lambda row: row["gt_source_id"])
    instances.sort(key=lambda row: row["instance_id"])
    replay.sort(key=lambda row: row["physical_unit_id"])
    return {"gt_sources": gt_sources, "instances": instances, "replay_index": replay}


# ---------------------------------------------------------------------------
# Executable gate verifiers (never trust registry self-reports).


def verify_kirc_source(
    root: Path, row: dict[str, str], eligible_physical: dict[str, dict[str, str]]
) -> list[str]:
    """Executable verifier for SPOT_BARCODE_TLS_ANNOTATION_CSV rows."""
    root = Path(root).resolve()
    errors: list[str] = []
    source_id = row["gt_source_id"]
    path = (root / row["path"]).resolve()
    try:
        path.relative_to(root)
    except ValueError:
        return [f"verifier_path_escape:{source_id}"]
    if not path.is_file():
        return [f"verifier_missing_file:{source_id}"]
    try:
        kind, labels = read_kirc_annotation(path)
    except ValueError as exc:
        return [f"verifier_annotation_schema:{source_id}:{exc}"]
    match = re.fullmatch(r"GSE175540_TLS_ANNOTATION::((?:ffpe|frozen)_[a-z]_\d+)", source_id)
    if match is None:
        errors.append(f"verifier_source_id_form:{source_id}")
        return errors
    alias = match.group(1)
    name_match = KIRC_ANNOTATION_RE.match(path.name)
    if name_match is None or name_match.group(2) != alias:
        errors.append(f"verifier_source_id_form:{source_id}")
        return errors
    gsm = name_match.group(1)
    expected_physical = f"external_geo::{gsm}"
    if row["physical_unit_id"] != expected_physical:
        errors.append(f"verifier_physical_link_mismatch:{source_id}")
    physical = eligible_physical.get(row["physical_unit_id"])
    if physical is not None and physical["patient_id"] != f"GSE175540::{alias}":
        errors.append(f"verifier_patient_raw_value_mismatch:{source_id}")
    if row["structure_id"] != "TLS":
        errors.append(f"verifier_structure_mismatch:{source_id}")
    if row["geometry_locator"] != f"{row['path']}#TLS_2_cat=TLS":
        errors.append(f"verifier_geometry_locator_mismatch:{source_id}")
    if kind != "TLS_2_cat":
        errors.append(f"verifier_header_kind_mismatch:{source_id}")
        return errors
    try:
        geometry = load_kirc_spot_geometry(str(root), gsm, alias)
    except FileNotFoundError:
        errors.append(f"verifier_missing_replay_asset:{source_id}")
        return errors
    if not [b for b in labels.get("TLS", []) if b in geometry]:
        errors.append(f"verifier_label_without_spots:{source_id}")
    return errors


def verify_usz_source(
    root: Path, row: dict[str, str], eligible_physical: dict[str, dict[str, str]]
) -> list[str]:
    """Executable verifier for H5AD_OBS_GROUND_TRUTH_LABELS rows."""
    root = Path(root).resolve()
    errors: list[str] = []
    source_id = row["gt_source_id"]
    path = (root / row["path"]).resolve()
    try:
        path.relative_to(root)
    except ValueError:
        return [f"verifier_path_escape:{source_id}"]
    if not path.is_file():
        return [f"verifier_missing_file:{source_id}"]
    match = re.fullmatch(r"USZ_H5AD_GROUND_TRUTH::((?:KC|LC)\d)", source_id)
    if match is None:
        errors.append(f"verifier_source_id_form:{source_id}")
        return errors
    sample = match.group(1)
    if path.name != f"{sample}.h5ad":
        errors.append(f"verifier_source_id_form:{source_id}")
        return errors
    if row["physical_unit_id"] != f"{USZ_LOGICAL_UNIT}::{sample}":
        errors.append(f"verifier_physical_link_mismatch:{source_id}")
    if row["structure_id"] != "TLS":
        errors.append(f"verifier_structure_mismatch:{source_id}")
    if row["geometry_locator"] != f"{row['path']}#obs/ground_truth=TLS":
        errors.append(f"verifier_geometry_locator_mismatch:{source_id}")
    try:
        obs = load_usz_obs(str(root), sample)
    except (FileNotFoundError, ValueError, KeyError, OSError) as exc:
        return errors + [f"verifier_h5ad_schema:{source_id}:{exc}"]
    if not [b for b, (_x, _y, label) in obs.items() if label == "TLS"]:
        errors.append(f"verifier_label_without_spots:{source_id}")
    return errors


def verify_stcrc_source(
    root: Path, row: dict[str, str], eligible_physical: dict[str, dict[str, str]]
) -> list[str]:
    """Executable verifier for SPOT_BARCODE_PATHOLOGY_CATEGORY_CSV rows."""
    root = Path(root).resolve()
    errors: list[str] = []
    source_id = row["gt_source_id"]
    path = (root / row["path"]).resolve()
    try:
        path.relative_to(root)
    except ValueError:
        return [f"verifier_path_escape:{source_id}"]
    if not path.is_file():
        return [f"verifier_missing_file:{source_id}"]
    try:
        labels = read_stcrc_annotation(path)
    except ValueError as exc:
        return [f"verifier_annotation_schema:{source_id}:{exc}"]
    match = re.fullmatch(r"STCRC_PATHOLOGY::(SN\d+_A\d+_Rep[12](?:_X)?)", source_id)
    if match is None:
        errors.append(f"verifier_source_id_form:{source_id}")
        return errors
    sample = match.group(1)
    if path.name != f"Pathologist_Annotations_{sample}.csv":
        errors.append(f"verifier_source_id_form:{source_id}")
        return errors
    if row["physical_unit_id"] != f"{STCRC_LOGICAL_UNIT}::{sample}":
        errors.append(f"verifier_physical_link_mismatch:{source_id}")
    physical = eligible_physical.get(row["physical_unit_id"])
    if physical is not None:
        stem = sample[: -len("_X")] if sample.endswith("_X") else sample
        patient_raw = stem.split("_")[1]
        if physical["patient_id"] != (
            f"{STCRC_LOGICAL_UNIT}::PATIENT::{_hash16(patient_raw)}"
        ):
            errors.append(f"verifier_patient_raw_value_mismatch:{source_id}")
    if row["structure_id"] != "TUMOR_STROMA_BOUNDARY":
        errors.append(f"verifier_structure_mismatch:{source_id}")
    locator_prefix = f"{row['path']}#tsb_label_family="
    if not row["geometry_locator"].startswith(locator_prefix):
        errors.append(f"verifier_geometry_locator_mismatch:{source_id}")
    elif tuple(row["geometry_locator"][len(locator_prefix):].split(",")) != STCRC_TSB_LABELS:
        errors.append(f"verifier_geometry_labels_mismatch:{source_id}")
    try:
        geometry = load_stcrc_spot_geometry(str(root), sample)
    except FileNotFoundError:
        errors.append(f"verifier_missing_replay_asset:{source_id}")
        return errors
    tsb_barcodes = [
        barcode
        for label in STCRC_TSB_LABELS
        for barcode in labels.get(label, [])
        if barcode in geometry
    ]
    if not tsb_barcodes:
        errors.append(f"verifier_label_without_spots:{source_id}")
    return errors


VERIFIER_FUNCTIONS = {
    VERIFIER_KIRC: verify_kirc_source,
    VERIFIER_USZ: verify_usz_source,
    VERIFIER_STCRC: verify_stcrc_source,
}
