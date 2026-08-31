#!/usr/bin/env python3
"""Shared deterministic Heiser 2023 pathology-annotation GT logic for R-02.

The registry builder and the independent gate validator both derive the Heiser
GT sources, structure instances, and h5ad replay index from the raw inputs
through this module:

- raw annotation CSVs in ``data/other_sources/htan/spatial_CRC_atlas_repo``
  (commit-pinned clone, see ``REPO_COMMIT``);
- the GitHub ``visium_sample_key.csv`` cross-checked against the verified R-01
  evidence asset ``repo/data_meta/ST_CRC_cohort_meta2.csv``;
- R-01 ``physical_units.tsv`` rows (patient/block identity is copied from, and
  independently re-derived against, the raw sample-key values);
- local OSF h5ad ``obs/_index``, ``obs/array_row``, ``obs/array_col`` only.

Identity granularity is the *piece*, not the capture area: one Visium capture
area can host several physical pieces (TMA cores), each with its own h5ad and
its own patient/block identity (e.g. capture 8578_AS_1 hosts pieces
WD33469/WD33473/WD33474 from three different patients). GT sources and
structure instances are therefore keyed by piece; the annotation CSVs are
capture-level and are partitioned by piece through the per-piece h5ad barcode
sets.

No expression values and no image pixels are read anywhere in this module
(D-033). The verifier never trusts registry TSV self-reports; every check here
recomputes from the raw files.
"""

from __future__ import annotations

import csv
import hashlib
import re
from functools import lru_cache
from pathlib import Path

import h5py


ANNOTATION_DIR = Path("data/other_sources/htan/spatial_CRC_atlas_repo/resources/ST")
H5AD_DIR = Path("data/other_sources/htan/HTAN_Vanderbilt_CRC_Visium_OSF_hftq2")
SAMPLE_KEY_PATH = ANNOTATION_DIR / "visium_sample_key.csv"
META2_PATH = Path("repo/data_meta/ST_CRC_cohort_meta2.csv")
PROVENANCE_NOTE = "infra/structure-registry/provenance/HEISER_2023_CELL_PMC10756562.md"
REPO_COMMIT = "64e585453514801a3b730f86aec90bbc0f0595da"
VERIFIER_CLASS = "SPOT_BARCODE_PATHOLOGY_ANNOTATION_CSV"
EXPECTED_ANNOTATION_CSV_COUNT = 42
EXPECTED_PIECE_COUNT = 48
LOGICAL_UNIT = "HTAN_VANDERBILT_CRC"

BARCODE_RE = re.compile(r"^[ACGT]{16}-1$")
LABEL_VOCABULARY = (
    "adenoma",
    "adenoma_border",
    "carcinoma",
    "carcinoma_border",
    "carcinoma_edge",
    "lymphoid_follicle",
    "normal_mucosa",
    "smooth_muscle",
)

CARCINOMA_TUMOR_TYPES = ("MSI-H", "MSS")
PRENEOPLASTIC_TUMOR_TYPES = ("SSL/HP", "TA/TVA")
NORMAL_TUMOR_TYPES = ("NL",)
CONTEXT_BY_TUMOR_TYPE = {
    **{kind: "CARCINOMA" for kind in CARCINOMA_TUMOR_TYPES},
    **{kind: "PRENEOPLASTIC" for kind in PRENEOPLASTIC_TUMOR_TYPES},
    **{kind: "NORMAL_MUCOSA" for kind in NORMAL_TUMOR_TYPES},
}

# (rule_id, structure_id, labels). A piece supports a rule when at least one
# of its own spots carries one of the rule labels.
STRUCTURE_RULES = (
    ("TLS_FOLLICLE", "TLS", ("lymphoid_follicle",)),
    (
        "CARCINOMA_BOUNDARY",
        "TUMOR_STROMA_BOUNDARY",
        ("carcinoma_border", "carcinoma_edge"),
    ),
    ("ADENOMA_BOUNDARY", "TUMOR_STROMA_BOUNDARY", ("adenoma_border",)),
)
RULES_BY_ID = {rule_id: (structure, labels) for rule_id, structure, labels in STRUCTURE_RULES}

INSTANCE_PREFIX = {"TLS": "TLS", "TUMOR_STROMA_BOUNDARY": "TSB"}
NONCONFIRMATORY_STATUS = {
    "PRENEOPLASTIC": "NONCONFIRMATORY_CONTEXT_PRENEOPLASTIC",
    "NORMAL_MUCOSA": "NONCONFIRMATORY_CONTEXT_NORMAL_MUCOSA",
}

GT_DEFINITION_MODALITY = (
    "manual pathology annotation on H&E in 10X Loupe Browser "
    "(Heiser et al. 2023, PMC10756562)"
)
BOUNDARY_UNCERTAINTY = (
    "spot-resolution (55 um microwell) manual annotation; per-file annotator "
    "identity is not individually attributed; see provenance note"
)
CONFIRMATORY_ALLOWED_USE = (
    "training-fold distance fields, core masking and localization benchmarks "
    "within frozen outer splits; molecular-only input budgets"
)
CONFIRMATORY_FORBIDDEN_USE = (
    "as model input feature or training label; any use with H&E- or "
    "image-derived inputs; internal/external validation claims"
)
CONTEXT_ALLOWED_USE = "discovery and context reference only"
CONTEXT_FORBIDDEN_USE = (
    "confirmatory claims, training, validation, localization benchmarks"
)

_NEIGHBOR_DELTAS = ((0, -2), (0, 2), (-1, -1), (-1, 1), (1, -1), (1, 1))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _hash16(raw: str) -> str:
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:16]


def expected_patient_id(patient_raw: str) -> str:
    return f"{LOGICAL_UNIT}::PATIENT::{_hash16(patient_raw)}"


def expected_block_id(patient_raw: str, block_raw: str) -> str:
    return f"{LOGICAL_UNIT}::BLOCK::{_hash16(patient_raw + '|' + block_raw)}"


@lru_cache(maxsize=1)
def load_piece_index(root_str: str) -> dict[str, dict[str, str]]:
    """piece key -> {capture, annotation_capture, h5ad, patient, block, tumor_type}.

    The piece key is the leading index column shared by both key files
    (``6723_1_WD86056`` style). Both files carry one row per piece; the GitHub
    ``visium_sample_key.csv`` and the verified R-01 evidence asset
    ``ST_CRC_cohort_meta2.csv`` must agree on every piece; any drift raises
    (fail-closed).
    """
    root = Path(root_str)
    per_file: list[dict[str, dict[str, str]]] = []
    for path in (META2_PATH, SAMPLE_KEY_PATH):
        records: dict[str, dict[str, str]] = {}
        with (root / path).open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                piece_key = (row.get("") or "").strip()
                if not piece_key:
                    raise ValueError(f"{path}: row without piece key")
                if piece_key in records:
                    raise ValueError(f"{path}: duplicate piece key {piece_key}")
                h5ad = Path(row["trimmed_adata"].strip()).name
                if not h5ad:
                    raise ValueError(
                        f"{path}:{piece_key} lacks trimmed_adata h5ad name"
                    )
                records[piece_key] = {
                    "capture": row["sample_key"].strip(),
                    "annotation_capture": row["sample_key_short"].strip(),
                    "block": row["block_name"].strip(),
                    "patient": row["patient_name"].strip(),
                    "tumor_type": row["tumor_type"].strip(),
                    "h5ad": h5ad,
                }
        per_file.append(records)
    meta2, github = per_file
    if set(meta2) != set(github):
        raise ValueError(
            "piece-key sets differ between meta2 and GitHub sample key: "
            f"meta2-only={sorted(set(meta2) - set(github))} "
            f"github-only={sorted(set(github) - set(meta2))}"
        )
    pieces: dict[str, dict[str, str]] = {}
    for piece_key in sorted(meta2):
        m_row, g_row = meta2[piece_key], github[piece_key]
        for field in (
            "capture",
            "annotation_capture",
            "block",
            "patient",
            "tumor_type",
            "h5ad",
        ):
            if m_row[field] != g_row[field]:
                raise ValueError(
                    f"sample-key/meta2 conflict for {piece_key}:{field}: "
                    f"{g_row[field]!r} vs {m_row[field]!r}"
                )
        pieces[piece_key] = dict(m_row)
    if len(pieces) != EXPECTED_PIECE_COUNT:
        raise ValueError(
            f"expected {EXPECTED_PIECE_COUNT} keyed pieces, observed {len(pieces)}"
        )
    return pieces


def annotation_csvs(root: Path) -> list[tuple[str, Path]]:
    paths = sorted((root / ANNOTATION_DIR).glob("*_pathology_annotation.csv"))
    if len(paths) != EXPECTED_ANNOTATION_CSV_COUNT:
        raise ValueError(
            f"expected {EXPECTED_ANNOTATION_CSV_COUNT} Heiser annotation CSVs, "
            f"observed {len(paths)}"
        )
    return [
        (path.name[: -len("_pathology_annotation.csv")], path) for path in paths
    ]


def read_annotation(path: Path) -> dict[str, list[str]]:
    """label -> barcodes. Empty annotation cells are unknown, never negative."""
    labels: dict[str, list[str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if header != ["Barcode", "pathology_annotation"]:
            raise ValueError(f"unexpected annotation header in {path}: {header}")
        seen: set[str] = set()
        for row_number, row in enumerate(reader, start=2):
            if len(row) != 2:
                raise ValueError(f"malformed row {row_number} in {path}")
            barcode, label = row[0].strip(), row[1].strip()
            if not BARCODE_RE.match(barcode):
                raise ValueError(f"invalid barcode {barcode!r} in {path}:{row_number}")
            if barcode in seen:
                raise ValueError(f"duplicate barcode {barcode!r} in {path}:{row_number}")
            seen.add(barcode)
            if not label:
                continue
            if label not in LABEL_VOCABULARY:
                raise ValueError(f"unknown label {label!r} in {path}:{row_number}")
            labels.setdefault(label, []).append(barcode)
    if not seen:
        raise ValueError(f"empty annotation file {path}")
    return labels


@lru_cache(maxsize=None)
def load_spot_index(root_str: str, h5ad_name: str) -> dict[str, tuple[int, int]]:
    """barcode -> (array_row, array_col). Identity/coordinate reads only."""
    path = Path(root_str) / H5AD_DIR / h5ad_name
    if not path.is_file():
        raise FileNotFoundError(f"missing local h5ad replay asset {path}")
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        raw_index = obs["_index"][:]
        array_row = obs["array_row"][:]
        array_col = obs["array_col"][:]
    index: dict[str, tuple[int, int]] = {}
    for raw, row, col in zip(raw_index, array_row, array_col):
        barcode = raw.decode() if isinstance(raw, bytes) else str(raw)
        if barcode in index:
            raise ValueError(f"duplicate barcode {barcode!r} in {path} obs index")
        index[barcode] = (int(row), int(col))
    return index


def spot_index_fingerprint(index: dict[str, tuple[int, int]]) -> str:
    digest = hashlib.sha256()
    digest.update(f"spots={len(index)};".encode())
    for barcode in sorted(index):
        row, col = index[barcode]
        digest.update(f"{barcode}\x00{row},{col};".encode())
    return digest.hexdigest()


def connected_components(points: set[tuple[int, int]]) -> list[list[tuple[int, int]]]:
    remaining = set(points)
    components: list[list[tuple[int, int]]] = []
    while remaining:
        stack = [remaining.pop()]
        component: list[tuple[int, int]] = []
        while stack:
            spot = stack.pop()
            component.append(spot)
            for d_row, d_col in _NEIGHBOR_DELTAS:
                neighbor = (spot[0] + d_row, spot[1] + d_col)
                if neighbor in remaining:
                    remaining.discard(neighbor)
                    stack.append(neighbor)
        components.append(sorted(component))
    components.sort(key=lambda comp: (min(comp), -len(comp)))
    return components


def _context_status(context: str) -> str:
    if context == "CARCINOMA":
        return "CONFIRMATORY"
    return NONCONFIRMATORY_STATUS[context]


def build_heiser_records(
    root: Path, physical_rows: list[dict[str, str]]
) -> dict[str, list[dict[str, str]]]:
    """Derive Heiser gt_source, instance, and h5ad replay-index rows."""
    root = Path(root).resolve()
    pieces = load_piece_index(str(root))
    by_capture: dict[str, list[str]] = {}
    for piece_key, record in pieces.items():
        by_capture.setdefault(record["annotation_capture"], []).append(piece_key)
    r01_htan = {
        row["physical_unit_id"]: row
        for row in physical_rows
        if row["study_id"] == LOGICAL_UNIT
    }
    gt_sources: list[dict[str, str]] = []
    instances: list[dict[str, str]] = []
    replay_index: dict[str, dict[str, str]] = {}

    # Classify every keyed piece once: eligible pieces replay to coordinates;
    # blocked pieces are fail-closed NOT_AUDITABLE rows.
    eligible_pieces: dict[str, tuple[dict[str, str], dict[str, str]]] = {}
    blocked_pieces: dict[str, tuple[dict[str, str], dict[str, str] | None, list[str]]] = {}
    for piece_key, record in pieces.items():
        physical_id = f"HTAN::{record['h5ad']}"
        r01_row = r01_htan.get(physical_id)
        local_asset = (root / H5AD_DIR / record["h5ad"]).is_file()
        reasons: list[str] = []
        if r01_row is None:
            reasons.append("no R-01 physical unit for keyed h5ad")
        elif r01_row["record_status"] != "RESOLVED_INCLUDED_CANDIDATE":
            reasons.append(
                f"R-01 unit status {r01_row['record_status']} is not eligible"
            )
        if not local_asset:
            reasons.append("local OSF h5ad replay asset missing")
        if reasons:
            blocked_pieces[piece_key] = (record, r01_row, reasons)
        else:
            assert r01_row is not None
            eligible_pieces[piece_key] = (record, r01_row)

    # Barcode ownership is scoped per capture: Visium barcode sequences are
    # reused across capture areas, so they are unique only within one capture.
    # Within a capture, a barcode in two pieces is an identity conflict
    # (raise). Labeled barcodes owned by no piece of their capture were lost
    # to upstream QC/trimming and are unusable; that is expected for the
    # filtered_trimmed h5ads and is not an error.
    def capture_owner(piece_keys: list[str]) -> dict[str, str]:
        scoped: dict[str, str] = {}
        for piece_key in piece_keys:
            if piece_key not in eligible_pieces:
                continue
            record, _r01_row = eligible_pieces[piece_key]
            for barcode in load_spot_index(str(root), record["h5ad"]):
                if barcode in scoped:
                    raise ValueError(
                        f"barcode {barcode} appears in two pieces of capture "
                        f"{record['annotation_capture']}: "
                        f"{scoped[barcode]} and {piece_key}"
                    )
                scoped[barcode] = piece_key
        return scoped

    for annotation_capture, csv_path in annotation_csvs(root):
        relative = csv_path.relative_to(root).as_posix()
        labels = read_annotation(csv_path)
        sha = _sha256(csv_path)
        piece_keys = sorted(by_capture.get(annotation_capture, []))

        if not piece_keys:
            r01_hint = r01_htan.get(
                f"HTAN::{annotation_capture.replace('_', '_AS_')}"
                "_filtered_trimmed.h5ad",
                {},
            )
            gt_sources.append(
                {
                    "gt_source_id": f"HEISER_PATHOLOGY::{annotation_capture}::UNLINKED",
                    "source_class": VERIFIER_CLASS,
                    "structure_id": "",
                    "path": relative,
                    "sha256": sha,
                    "checksum_status": "VERIFIED",
                    "schema_locator": "csv header:Barcode,pathology_annotation",
                    "physical_unit_id": "",
                    "physical_link_status": "UNRESOLVED_NO_SAMPLE_KEY_ROW",
                    "geometry_locator": "",
                    "provenance_status": "KNOWN",
                    "provenance_locator": PROVENANCE_NOTE,
                    "overlap_status": "UNKNOWN",
                    "overlap_locator": "",
                    "same_assay_status": "UNKNOWN",
                    "audit_status": "NOT_AUDITABLE_GT",
                    "notes": (
                        "No visium_sample_key or meta2 row links this capture area "
                        "to an R-01 physical unit; fail-closed. R-01 unit for the "
                        f"matching h5ad name has status "
                        f"{r01_hint.get('record_status', 'ABSENT')}."
                    ),
                }
            )
            continue

        # Enforce within-capture piece disjointness (identity integrity).
        # Per-piece filtering below uses each piece's own spot index, which
        # is unambiguous because barcodes are unique within a capture.
        capture_owner(piece_keys)

        for piece_key in piece_keys:
            if piece_key not in blocked_pieces:
                continue
            record, r01_row, reasons = blocked_pieces[piece_key]
            physical_id = f"HTAN::{record['h5ad']}"
            gt_sources.append(
                {
                    "gt_source_id": f"HEISER_PATHOLOGY::{piece_key}::REPLAY_BLOCKED",
                    "source_class": VERIFIER_CLASS,
                    "structure_id": "",
                    "path": relative,
                    "sha256": sha,
                    "checksum_status": "VERIFIED",
                    "schema_locator": "csv header:Barcode,pathology_annotation",
                    "physical_unit_id": physical_id if r01_row else "",
                    "physical_link_status": "R01_UNIT_NOT_ELIGIBLE_OR_ASSET_MISSING",
                    "geometry_locator": "",
                    "provenance_status": "KNOWN",
                    "provenance_locator": PROVENANCE_NOTE,
                    "overlap_status": "UNKNOWN",
                    "overlap_locator": "",
                    "same_assay_status": "UNKNOWN",
                    "audit_status": "NOT_AUDITABLE_GT",
                    "notes": (
                        f"Keyed piece {piece_key} cannot be replayed to "
                        "coordinates: "
                        + "; ".join(reasons)
                        + ". Fail-closed pending a single OSF h5ad fetch and/or "
                        "R-01 status review."
                    ),
                }
            )

        for piece_key in piece_keys:
            if piece_key not in eligible_pieces:
                continue
            record, r01_row = eligible_pieces[piece_key]
            h5ad_name = record["h5ad"]
            physical_id = f"HTAN::{h5ad_name}"
            spot_index = load_spot_index(str(root), h5ad_name)
            replay_index[physical_id] = {
                "physical_unit_id": physical_id,
                "path": (H5AD_DIR / h5ad_name).as_posix(),
                "spot_count": str(len(spot_index)),
                "index_fingerprint_sha256": spot_index_fingerprint(spot_index),
                "fingerprint_status": "INDEX_VALUE_PINNED",
            }
            context = CONTEXT_BY_TUMOR_TYPE.get(record["tumor_type"])
            if context is None:
                raise ValueError(
                    f"unknown tumor_type {record['tumor_type']!r} for piece {piece_key}"
                )

            for rule_id, structure_id, rule_labels in STRUCTURE_RULES:
                barcodes = [
                    barcode
                    for label in rule_labels
                    for barcode in labels.get(label, [])
                    if barcode in spot_index
                ]
                if not barcodes:
                    continue
                geometry_locator = (
                    f"{relative}#pathology_annotation={','.join(rule_labels)}"
                )
                status = _context_status(context)
                source_id = f"HEISER_PATHOLOGY::{piece_key}::{rule_id}"
                gt_sources.append(
                    {
                        "gt_source_id": source_id,
                        "source_class": VERIFIER_CLASS,
                        "structure_id": structure_id,
                        "path": relative,
                        "sha256": sha,
                        "checksum_status": "VERIFIED",
                        "schema_locator": "csv header:Barcode,pathology_annotation",
                        "physical_unit_id": physical_id,
                        "physical_link_status": "VERIFIED",
                        "geometry_locator": geometry_locator,
                        "provenance_status": "KNOWN",
                        "provenance_locator": PROVENANCE_NOTE,
                        "overlap_status": "KNOWN_OVERLAP",
                        "overlap_locator": (
                            "infra/structure-registry/input_policy.tsv"
                            "#he-annotation-gt-input-exclusion"
                        ),
                        "same_assay_status": "INDEPENDENT_ASSAY",
                        "audit_status": "AUDITABLE_GT",
                        "notes": (
                            f"spot-level manual pathology annotation (Heiser 2023, "
                            f"commit {REPO_COMMIT}); {len(barcodes)} label-group "
                            f"spots in piece {piece_key}; context {context} from "
                            f"tumor_type {record['tumor_type']}; per-file "
                            "annotator identity unattributed"
                        ),
                    }
                )
                components = connected_components(
                    {spot_index[barcode] for barcode in barcodes}
                )
                for index, component in enumerate(components, start=1):
                    instance_id = (
                        f"{INSTANCE_PREFIX[structure_id]}::HEISER::{piece_key}"
                        f"::comp{index:02d}"
                    )
                    confirmatory = status == "CONFIRMATORY"
                    instances.append(
                        {
                            "instance_id": instance_id,
                            "structure_id": structure_id,
                            "logical_unit_id": LOGICAL_UNIT,
                            "physical_unit_id": physical_id,
                            "patient_id": r01_row["patient_id"],
                            "block_id": r01_row["block_id"],
                            "physical_specimen_id": r01_row["physical_specimen_id"],
                            "gt_source_id": source_id,
                            "gt_definition_modality": GT_DEFINITION_MODALITY,
                            "gt_geometry_locator": geometry_locator,
                            "boundary_uncertainty": BOUNDARY_UNCERTAINTY,
                            "confirmation_status": status,
                            "allowed_use": (
                                CONFIRMATORY_ALLOWED_USE
                                if confirmatory
                                else CONTEXT_ALLOWED_USE
                            ),
                            "forbidden_use": (
                                CONFIRMATORY_FORBIDDEN_USE
                                if confirmatory
                                else CONTEXT_FORBIDDEN_USE
                            ),
                        }
                    )

    gt_sources.sort(key=lambda row: row["gt_source_id"])
    instances.sort(key=lambda row: row["instance_id"])
    replay = [replay_index[key] for key in sorted(replay_index)]
    return {"gt_sources": gt_sources, "instances": instances, "replay_index": replay}


def verify_annotation_source(
    root: Path,
    row: dict[str, str],
    eligible_physical: dict[str, dict[str, str]],
) -> list[str]:
    """Executable verifier for SPOT_BARCODE_PATHOLOGY_ANNOTATION_CSV rows."""
    root = Path(root).resolve()
    errors: list[str] = []
    source_id = row["gt_source_id"]
    relative = row["path"]
    path = (root / relative).resolve()
    try:
        path.relative_to(root)
    except ValueError:
        return [f"verifier_path_escape:{source_id}"]
    if not path.is_file():
        return [f"verifier_missing_file:{source_id}"]

    try:
        labels = read_annotation(path)
    except ValueError as exc:
        return [f"verifier_annotation_schema:{source_id}:{exc}"]

    locator_prefix = f"{relative}#pathology_annotation="
    if not row["geometry_locator"].startswith(locator_prefix):
        errors.append(f"verifier_geometry_locator_mismatch:{source_id}")
        locator_labels: tuple[str, ...] = ()
    else:
        locator_labels = tuple(row["geometry_locator"][len(locator_prefix):].split(","))

    match = re.fullmatch(
        r"HEISER_PATHOLOGY::(.+)::(TLS_FOLLICLE|CARCINOMA_BOUNDARY|ADENOMA_BOUNDARY)",
        source_id,
    )
    if match is None:
        errors.append(f"verifier_source_id_form:{source_id}")
        return errors
    piece_key, rule_id = match.groups()
    structure_id, rule_labels = RULES_BY_ID[rule_id]
    if row["structure_id"] != structure_id:
        errors.append(f"verifier_structure_mismatch:{source_id}")
    if locator_labels != rule_labels:
        errors.append(f"verifier_geometry_labels_mismatch:{source_id}")

    try:
        pieces = load_piece_index(str(root))
    except ValueError as exc:
        return errors + [f"verifier_crosswalk_conflict:{source_id}:{exc}"]
    piece = pieces.get(piece_key)
    if piece is None:
        errors.append(f"verifier_unkeyed_piece:{source_id}")
        return errors

    expected_physical = f"HTAN::{piece['h5ad']}"
    if row["physical_unit_id"] != expected_physical:
        errors.append(f"verifier_physical_link_mismatch:{source_id}")
    physical = eligible_physical.get(row["physical_unit_id"])
    if physical is not None:
        if physical["patient_id"] != expected_patient_id(piece["patient"]):
            errors.append(f"verifier_patient_raw_value_mismatch:{source_id}")
        if physical["block_id"] != expected_block_id(piece["patient"], piece["block"]):
            errors.append(f"verifier_block_raw_value_mismatch:{source_id}")

    try:
        spot_index = load_spot_index(str(root), piece["h5ad"])
    except FileNotFoundError:
        errors.append(f"verifier_missing_h5ad_replay_asset:{source_id}")
        return errors
    group_barcodes = [
        barcode
        for label in rule_labels
        for barcode in labels.get(label, [])
        if barcode in spot_index
    ]
    if not group_barcodes:
        errors.append(f"verifier_label_without_spots:{source_id}")
    return errors
