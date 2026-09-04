"""Strict alignment of training-role structure annotations.

The registry, rather than a filename or expression pattern, decides whether a
spot annotation may be used.  Missing barcodes remain ``NaN`` and are never
silently converted to structural negatives.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np


STRUCTURE_GT_SCHEMA = "r04.training_structure_gt.v1"


def _read_json(path: Path) -> dict[str, object]:
    import json

    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _fragment_values(locator: str) -> tuple[Path, tuple[str, ...]]:
    path_text, _, fragment = locator.partition("#")
    if not path_text or not fragment:
        raise ValueError("GT geometry locator must contain a path and semantic fragment")
    match = re.search(r"(?:^|,)pathology_annotation=([^,]+(?:,[^,]+)*)", fragment)
    if match is None:
        match = re.search(r"pathology_annotation=([^,]+)", fragment)
    if match is None:
        raise ValueError("GT locator has no pathology annotation value declaration")
    values = tuple(
        value.strip().lower()
        for value in match.group(1).split(",")
        if value.strip()
    )
    if not values:
        raise ValueError("GT locator declares no positive annotation values")
    return Path(path_text), values


def _read_annotation_csv(path: Path) -> dict[str, str]:
    if not path.is_file():
        raise ValueError(f"training GT annotation is missing: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"training GT annotation has no header: {path}")
        columns = {str(name).strip().lower(): name for name in reader.fieldnames if name is not None}
        barcode_column = columns.get("barcode")
        value_column = columns.get("pathology_annotation")
        if barcode_column is None or value_column is None:
            raise ValueError(f"training GT annotation columns are incomplete: {path}")
        result: dict[str, str] = {}
        for row in reader:
            barcode = str(row.get(barcode_column, "")).strip()
            if not barcode:
                continue
            if barcode in result:
                raise ValueError(f"duplicate training GT barcode: {path}::{barcode}")
            result[barcode] = str(row.get(value_column, "")).strip().lower()
    if not result:
        raise ValueError(f"training GT annotation is empty: {path}")
    return result


def load_training_gt_catalog(
    registry_path: Path,
    training_manifest_path: Path,
    *,
    structures: Sequence[str] = ("TLS", "TUMOR_STROMA_BOUNDARY"),
) -> dict[str, object]:
    """Load only confirmatory training-fold GT from the structure registry."""
    manifest = _read_json(training_manifest_path)
    if manifest.get("schema") != "r04.role_manifest.v1" or manifest.get("status") != "READY" or manifest.get("role") != "training":
        raise ValueError("structure GT loader requires a READY training manifest")
    rows = manifest.get("rows")
    if not isinstance(rows, list) or not rows:
        raise ValueError("training manifest has no rows")
    training_sections = {str(row["section_id"]) for row in rows if isinstance(row, dict) and "section_id" in row}
    if not training_sections:
        raise ValueError("training manifest has no section IDs")

    requested = tuple(str(value) for value in structures)
    if len(set(requested)) != len(requested) or not requested:
        raise ValueError("structure names must be unique and non-empty")
    with registry_path.open(newline="", encoding="utf-8-sig") as handle:
        registry_rows = list(csv.DictReader(handle, delimiter="\t"))
    catalog: dict[str, dict[str, dict[str, object]]] = {name: {} for name in requested}
    skipped = {name: 0 for name in requested}
    loaded_files: dict[str, dict[str, object]] = {}
    for row in registry_rows:
        structure = str(row.get("structure_id", ""))
        section_id = str(row.get("physical_unit_id", ""))
        if structure not in catalog or section_id not in training_sections:
            continue
        if row.get("confirmation_status") != "CONFIRMATORY" or "training-fold" not in str(row.get("allowed_use", "")):
            skipped[structure] += 1
            continue
        locator = str(row.get("gt_geometry_locator", ""))
        relative_path, positive_values = _fragment_values(locator)
        path = relative_path if relative_path.is_absolute() else registry_path.parent.parent.parent / relative_path
        # The project registry stores paths relative to the project root.  The
        # fallback above is only for an explicitly project-relative locator.
        if not path.is_file():
            path = registry_path.parent.parent.parent / relative_path
        key = str(path.resolve())
        if key not in loaded_files:
            loaded_files[key] = {
                "path": str(path),
                "values": _read_annotation_csv(path),
                "positive_values": set(positive_values),
            }
        else:
            loaded_files[key]["positive_values"] = set(loaded_files[key]["positive_values"]) | set(positive_values)
        existing = catalog[structure].get(section_id)
        if existing is None:
            catalog[structure][section_id] = {
                "path": str(path),
                "positive_values": sorted(set(positive_values)),
                "annotation": loaded_files[key]["values"],
            }
        else:
            if existing["path"] != str(path):
                raise ValueError(f"multiple GT files disagree for {structure}/{section_id}")
            existing["positive_values"] = sorted(
                set(existing["positive_values"]) | set(positive_values)
            )
    return {
        "schema": STRUCTURE_GT_SCHEMA,
        "registry_path": str(registry_path),
        "training_manifest_path": str(training_manifest_path),
        "training_section_count": len(training_sections),
        "structures": {
            name: {
                "section_count": len(catalog[name]),
                "sections": catalog[name],
                "skipped_ineligible_rows": skipped[name],
            }
            for name in requested
        },
        "loaded_annotation_files": loaded_files,
        "role_boundary": "training_role_confirmatory_training_fold_localization_readout_only",
        "internal_external_validation_gt": "SEALED_NOT_READ",
    }


def labels_for_spots(
    catalog: Mapping[str, object],
    structure: str,
    section_ids: Sequence[object],
    barcodes: Sequence[object],
) -> np.ndarray:
    """Return tri-state labels aligned to exact ``(section_id, barcode)``."""
    if len(section_ids) != len(barcodes):
        raise ValueError("section IDs and barcodes are not aligned")
    structures = catalog.get("structures")
    if not isinstance(structures, Mapping) or structure not in structures:
        raise ValueError(f"structure is absent from training GT catalog: {structure}")
    entry = structures[structure]
    if not isinstance(entry, Mapping) or not isinstance(entry.get("sections"), Mapping):
        raise ValueError("training GT catalog section map is invalid")
    labels = np.full(len(section_ids), np.nan, dtype=float)
    for index, (section_id, barcode) in enumerate(zip(section_ids, barcodes)):
        section = entry["sections"].get(str(section_id))
        if not isinstance(section, Mapping):
            continue
        annotation = section.get("annotation")
        positive_values = {str(value).lower() for value in section.get("positive_values", [])}
        if not isinstance(annotation, Mapping) or not positive_values:
            raise ValueError("training GT section annotation is invalid")
        value = annotation.get(str(barcode))
        if value is None or not str(value).strip():
            continue
        labels[index] = 1.0 if str(value).lower() in positive_values else 0.0
    return labels
