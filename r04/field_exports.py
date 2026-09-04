"""Auditable persistence for held-out continuous spatial effects.

The exporter stores the held-out field values and the gene loading needed to
reconstruct a rotation-stable, centered linear spot-by-gene effect.  It does
not assign biological names to factors and it does not read ground truth.
"""

from __future__ import annotations

import hashlib
import tempfile
from pathlib import Path
from typing import Sequence

import numpy as np

from .diagnostics import field_effect_matrix, mnsf_spatial_effect_matrix
from .runtime import atomic_json
from .types import FieldFit, SectionData


SCHEMA = "r04.heldout_field_export.v1"
EFFECT_DEFINITION = (
    "centered_linear_spot_by_gene_effect=section_center(field_mean)@evaluation_loading.T"
)
RATE_DEFINITION = (
    "section_center(exp(field_mean)@"
    "(evaluation_loading*factor_amplitude).T)"
)


def _as_string_array(values: Sequence[object]) -> np.ndarray:
    return np.asarray([str(value) for value in values], dtype="U")


def _atomic_npz(path: Path, arrays: dict[str, np.ndarray]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="wb", suffix=".npz", dir=path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        np.savez_compressed(handle, **arrays)
    try:
        temporary.replace(path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _field_matrix(
    sections: Sequence[SectionData], groups: Sequence[Sequence[FieldFit]]
) -> tuple[np.ndarray, np.ndarray]:
    if len(sections) != len(groups):
        raise ValueError("field groups and sections must have equal length")
    expected_k: int | None = None
    means: list[np.ndarray] = []
    sds: list[np.ndarray] = []
    for section, group in zip(sections, groups):
        if any(field.section_id != section.section_id for field in group):
            raise ValueError("field section IDs do not match validation sections")
        k = len(group)
        if expected_k is None:
            expected_k = k
        if k != expected_k:
            raise ValueError("field factor counts differ between sections")
        n_spots = len(section.barcode)
        if k:
            mean = np.column_stack([np.asarray(field.field_mean, dtype=float) for field in group])
            sd = np.column_stack([np.asarray(field.field_sd, dtype=float) for field in group])
        else:
            mean = np.empty((n_spots, 0), dtype=float)
            sd = np.empty((n_spots, 0), dtype=float)
        if mean.shape != (n_spots, k) or sd.shape != (n_spots, k):
            raise ValueError("field arrays are not spot-aligned")
        means.append(mean)
        sds.append(sd)
    if expected_k is None or expected_k < 1:
        raise ValueError("held-out field export requires at least one factor")
    return np.concatenate(means, axis=0), np.concatenate(sds, axis=0)


def _section_center(matrix: np.ndarray, section_ids: np.ndarray) -> np.ndarray:
    centered = np.asarray(matrix, dtype=float).copy()
    for section_id in np.unique(section_ids):
        indices = np.flatnonzero(section_ids == section_id)
        centered[indices] -= centered[indices].mean(axis=0, keepdims=True)
    return centered


def write_heldout_field_export(
    path: Path,
    *,
    sections: Sequence[SectionData],
    field_groups: Sequence[Sequence[FieldFit]],
    evaluation_loading: np.ndarray,
    factor_amplitude: np.ndarray,
    frozen_gene_ids: Sequence[str],
    adaptation_gene_ids: Sequence[str],
    evaluation_gene_ids: Sequence[str],
    adaptation_indices: np.ndarray,
    evaluation_indices: np.ndarray,
    provenance: dict[str, object],
) -> dict[str, object]:
    """Persist one held-out or training field export.

    ``field_groups`` are returned by ``MNSFEstimator.infer`` in the exact
    section order supplied here.  The exported loading is restricted to the
    evaluation genes, so a downstream readout cannot accidentally use the
    genes that adapted the held-out field.
    ``adaptation_gene_ids`` may be empty when the export intentionally uses the
    complete frozen gene panel for a post-fit structure readout.
    """
    field_mean, field_sd = _field_matrix(sections, field_groups)
    loading = np.asarray(evaluation_loading, dtype=float)
    amplitude = np.asarray(factor_amplitude, dtype=float)
    frozen = tuple(str(value) for value in frozen_gene_ids)
    adaptation = tuple(str(value) for value in adaptation_gene_ids)
    evaluation = tuple(str(value) for value in evaluation_gene_ids)
    adaptation_idx = np.asarray(adaptation_indices, dtype=int)
    evaluation_idx = np.asarray(evaluation_indices, dtype=int)
    if loading.ndim != 2 or loading.shape[1] != field_mean.shape[1]:
        raise ValueError("evaluation loading and field factor dimensions differ")
    if loading.shape[0] != len(evaluation) or amplitude.shape != (field_mean.shape[1],):
        raise ValueError("evaluation loading or amplitude shape is invalid")
    if len(set(frozen)) != len(frozen):
        raise ValueError("frozen gene universe contains duplicates")
    if len(frozen) != len(adaptation) + len(evaluation):
        raise ValueError("gene split does not cover the frozen gene universe")
    if set(adaptation).intersection(evaluation) or set(adaptation) | set(evaluation) != set(frozen):
        raise ValueError("adaptation and evaluation genes must partition the frozen universe")
    if adaptation_idx.shape != (len(adaptation),) or evaluation_idx.shape != (len(evaluation),):
        raise ValueError("gene split index shapes are invalid")
    if not np.array_equal(np.sort(np.concatenate([adaptation_idx, evaluation_idx])), np.arange(len(frozen))):
        raise ValueError("gene split indices must partition frozen gene positions")
    if np.any(~np.isfinite(field_mean)) or np.any(~np.isfinite(field_sd)):
        raise ValueError("field export contains non-finite values")
    if np.any(field_sd < 0) or np.any(~np.isfinite(loading)) or np.any(loading < 0):
        raise ValueError("field uncertainty or loading is invalid")
    if np.any(~np.isfinite(amplitude)) or np.any(amplitude <= 0):
        raise ValueError("factor amplitude is invalid")

    barcodes = _as_string_array([barcode for section in sections for barcode in section.barcode])
    section_ids = _as_string_array([
        section.section_id for section in sections for _ in section.barcode
    ])
    patient_ids = _as_string_array([
        section.patient_id for section in sections for _ in section.barcode
    ])
    block_ids = _as_string_array([
        section.block_id or "" for section in sections for _ in section.barcode
    ])
    lineages = _as_string_array([
        section.lineage for section in sections for _ in section.barcode
    ])
    coords = np.concatenate([np.asarray(section.coords, dtype=float) for section in sections], axis=0)
    section_offsets = np.cumsum(
        np.asarray([0, *[len(section.barcode) for section in sections]], dtype=int)
    )
    if len(barcodes) != len(field_mean) or coords.shape != (len(field_mean), 2):
        raise ValueError("held-out metadata is not spot-aligned")
    rate_effect = mnsf_spatial_effect_matrix(field_mean, loading, amplitude)
    rate_effect = _section_center(rate_effect, section_ids)

    arrays = {
        "field_mean": field_mean.astype(np.float32),
        "field_sd": field_sd.astype(np.float32),
        "spatial_rate_effect_centered": rate_effect.astype(np.float32),
        "evaluation_loading": loading.astype(np.float32),
        "factor_amplitude": amplitude.astype(np.float32),
        "barcodes": barcodes,
        "section_ids": section_ids,
        "patient_ids": patient_ids,
        "block_ids": block_ids,
        "lineages": lineages,
        "coords": coords.astype(np.float32),
        "section_offsets": section_offsets,
        "adaptation_indices": adaptation_idx,
        "evaluation_indices": evaluation_idx,
        "frozen_gene_ids": _as_string_array(frozen),
        "adaptation_gene_ids": _as_string_array(adaptation),
        "evaluation_gene_ids": _as_string_array(evaluation),
    }
    digest = _atomic_npz(path, arrays)
    array_manifest = {
        key: {
            "shape": list(value.shape),
            "dtype": str(value.dtype),
            "sha256": hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest(),
        }
        for key, value in arrays.items()
    }
    metadata = {
        "schema": SCHEMA,
        "npz_path": str(path),
        "npz_sha256": digest,
        "n_spots": int(field_mean.shape[0]),
        "k_model": int(field_mean.shape[1]),
        "evaluation_gene_count": int(len(evaluation)),
        "field_mean_centering": "within_section_zero_mean",
        "effect_representation": "centered_mnsf_spatial_rate_contribution_primary",
        "effect_definition": RATE_DEFINITION,
        "rate_contribution_definition": RATE_DEFINITION,
        "linearized_effect_definition": "section_center(field_mean)@evaluation_loading.T",
        "linearized_effect_role": "secondary_rotation_stable_diagnostic",
        "array_manifest": array_manifest,
        "single_factor_biological_naming": "forbidden",
        "nested_gene_crossfit": bool(adaptation),
        "gene_panel_role": (
            "nested_gene_crossfit_evaluation_genes"
            if adaptation else "full_frozen_model_panel"
        ),
        "provenance": provenance,
    }
    metadata_path = path.with_suffix(".json")
    atomic_json(metadata_path, metadata)
    metadata["metadata_path"] = str(metadata_path)
    atomic_json(metadata_path, metadata)
    return metadata


def read_heldout_field_export(path: Path) -> dict[str, np.ndarray | dict[str, object]]:
    """Read and fail closed on a tampered held-out field export."""
    metadata_path = path.with_suffix(".json")
    if not path.is_file() or not metadata_path.is_file():
        raise ValueError("held-out field export or metadata is missing")
    import json

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("schema") != SCHEMA or metadata.get("npz_path") != str(path):
        raise ValueError("held-out field export schema/path mismatch")
    if hashlib.sha256(path.read_bytes()).hexdigest() != metadata.get("npz_sha256"):
        raise ValueError("held-out field export hash mismatch")
    with np.load(path, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    required = {
        "field_mean", "field_sd", "spatial_rate_effect_centered", "evaluation_loading", "factor_amplitude",
        "barcodes", "section_ids", "patient_ids", "block_ids", "lineages", "coords", "section_offsets", "adaptation_indices",
        "evaluation_indices", "frozen_gene_ids", "adaptation_gene_ids", "evaluation_gene_ids",
    }
    if set(arrays) != required:
        raise ValueError("held-out field export array schema mismatch")
    if arrays["field_mean"].ndim != 2 or arrays["field_sd"].shape != arrays["field_mean"].shape:
        raise ValueError("held-out field export field shape mismatch")
    if arrays["spatial_rate_effect_centered"].shape != (
        len(arrays["field_mean"]), len(arrays["evaluation_gene_ids"])
    ):
        raise ValueError("held-out field export rate-effect shape mismatch")
    if arrays["evaluation_loading"].shape != (
        len(arrays["evaluation_gene_ids"]), arrays["field_mean"].shape[1]
    ):
        raise ValueError("held-out field export loading shape mismatch")
    if arrays["coords"].shape != (len(arrays["field_mean"]), 2):
        raise ValueError("held-out field export coordinate shape mismatch")
    if len(arrays["barcodes"]) != len(arrays["field_mean"]):
        raise ValueError("held-out field export barcode alignment mismatch")
    if (
        not np.isfinite(arrays["field_mean"]).all()
        or not np.isfinite(arrays["field_sd"]).all()
        or np.any(arrays["field_sd"] < 0)
        or not np.isfinite(arrays["spatial_rate_effect_centered"]).all()
    ):
        raise ValueError("held-out field export contains invalid numeric values")
    offsets = arrays["section_offsets"]
    if offsets.ndim != 1 or len(offsets) < 2 or int(offsets[0]) != 0 or int(offsets[-1]) != len(arrays["field_mean"]):
        raise ValueError("held-out field export section offsets are invalid")
    if np.any(np.diff(offsets) <= 0):
        raise ValueError("held-out field export section offsets are not increasing")
    array_manifest = metadata.get("array_manifest")
    if not isinstance(array_manifest, dict):
        raise ValueError("held-out field export array manifest is missing")
    for key, value in arrays.items():
        expected = array_manifest.get(key)
        if not isinstance(expected, dict) or expected.get("sha256") != hashlib.sha256(
            np.ascontiguousarray(value).tobytes()
        ).hexdigest():
            raise ValueError(f"held-out field export array hash mismatch: {key}")
    return {**arrays, "metadata": metadata}


def reconstructed_linear_effect(export: dict[str, object]) -> np.ndarray:
    """Reconstruct the rotation-stable centered linear spot-by-gene effect."""
    field_mean = np.asarray(export["field_mean"], dtype=float)
    loading = np.asarray(export["evaluation_loading"], dtype=float)
    return field_effect_matrix(field_mean, loading)


def reconstructed_rate_effect(export: dict[str, object]) -> np.ndarray:
    """Return the primary section-centered mNSF spatial-rate contribution."""
    if "spatial_rate_effect_centered" in export:
        return np.asarray(export["spatial_rate_effect_centered"], dtype=float)
    field_mean = np.asarray(export["field_mean"], dtype=float)
    loading = np.asarray(export["evaluation_loading"], dtype=float)
    amplitude = np.asarray(export["factor_amplitude"], dtype=float)
    section_ids = np.asarray(export["section_ids"])
    return _section_center(
        mnsf_spatial_effect_matrix(field_mean, loading, amplitude), section_ids
    )
