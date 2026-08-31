"""Post-freeze known-structure anchor evaluation.

No discovery module imports this module.  The candidate registry hash is
checked before any GT file is read, and unknown labels are never negatives.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from typing import Iterable

import numpy as np

from .types import GateResult


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_confirmatory_gt(path: Path) -> list[dict[str, str]]:
    """Read only after the caller has verified the frozen candidate hash."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"section_id", "barcode", "structure_id", "label"}
    if not rows or not required <= set(rows[0]):
        raise ValueError("GT anchor table requires section_id, barcode, structure_id and label")
    return [row for row in rows if row["label"].lower() in {"positive", "confirmatory"}]


def evaluate_anchor_correspondence(
    *,
    candidate_registry_hash: str,
    expected_registry_hash: str,
    field_values: np.ndarray,
    anchor_values: np.ndarray,
    groups: Iterable[str],
) -> GateResult:
    if candidate_registry_hash != expected_registry_hash:
        return GateResult("BLOCKED_CANDIDATE_HASH_MISMATCH", ("candidate_registry_hash_mismatch",))
    field = np.asarray(field_values, dtype=float)
    anchor = np.asarray(anchor_values, dtype=float)
    if field.shape != anchor.shape:
        return GateResult("SCIENTIFIC_FAIL_ANCHOR_SHAPE", ("field_anchor_shape_mismatch",))
    groups = np.asarray(list(groups))
    if len(np.unique(groups)) < 3:
        return GateResult("NOT_TESTABLE", ("fewer_than_three_outer_groups",))
    if field.std() == 0 or anchor.std() == 0:
        return GateResult("SCIENTIFIC_FAIL_NO_ANCHOR_VARIATION", ("constant_anchor_or_field",))
    correlation = float(np.corrcoef(field, anchor)[0, 1])
    return GateResult("ANCHOR_EVALUATED", metrics={"field_anchor_correlation": correlation})
