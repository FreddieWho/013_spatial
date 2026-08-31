"""Small, serialisable R-04 data contracts."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np
from scipy import sparse


class R04ContractError(ValueError):
    """Raised when a data or candidate contract is violated."""


@dataclass(frozen=True)
class SectionData:
    """Counts and coordinates for one physical section.

    ``counts`` is spot-by-gene and must contain non-negative integer counts.
    ``barcode`` remains an opaque string; no suffix stripping or cross-capture
    deduplication is permitted.
    """

    section_id: str
    patient_id: str
    block_id: str | None
    lineage: str
    barcode: tuple[str, ...]
    coords: np.ndarray
    counts: Any
    gene_id: tuple[str, ...]
    library_size: np.ndarray | None = None
    metadata: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.section_id or not self.patient_id or not self.lineage:
            raise R04ContractError("section, patient and lineage IDs are required")
        coords = np.asarray(self.coords)
        if coords.ndim != 2 or coords.shape[1] != 2:
            raise R04ContractError("coords must have shape (n_spots, 2)")
        if len(self.barcode) != coords.shape[0]:
            raise R04ContractError("barcode and coordinate lengths differ")
        if len(self.gene_id) == 0:
            raise R04ContractError("gene_id cannot be empty")
        shape = getattr(self.counts, "shape", None)
        if shape != (coords.shape[0], len(self.gene_id)):
            raise R04ContractError("counts must be spot-by-gene")
        values = self.counts.data if sparse.issparse(self.counts) else np.asarray(self.counts)
        if np.any(~np.isfinite(values)) or np.any(values < 0) or np.any(values != np.floor(values)):
            raise R04ContractError("counts must be finite, non-negative integer counts")
        if not np.isfinite(coords).all():
            raise R04ContractError("coordinates must be finite")
        if self.library_size is not None and len(self.library_size) != len(self.barcode):
            raise R04ContractError("library_size length differs from barcode length")


@dataclass(frozen=True)
class FieldFit:
    """Posterior-like output for one continuous field.

    Arrays are spot-aligned to the input section.  ``field_mean`` is a value
    of a continuous function at observed coordinates, not a class assignment.
    """

    model_id: str
    factor_id: str
    section_id: str
    loading: np.ndarray
    field_mean: np.ndarray
    field_sd: np.ndarray
    lengthscale: float
    spatial_variance: float
    input_hash: str
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.model_id or not self.factor_id or not self.section_id:
            raise R04ContractError("field identifiers are required")
        if self.lengthscale <= 0 or self.spatial_variance < 0:
            raise R04ContractError("invalid field scale or variance")
        if self.loading.ndim != 1 or not np.isfinite(self.loading).all():
            raise R04ContractError("field loading must be a finite one-dimensional array")
        if self.field_mean.ndim != 1 or self.field_sd.shape != self.field_mean.shape:
            raise R04ContractError("field posterior arrays must be one-dimensional")
        if not np.isfinite(self.field_mean).all() or not np.isfinite(self.field_sd).all():
            raise R04ContractError("field posterior contains non-finite values")
        if (self.field_sd < 0).any():
            raise R04ContractError("field SD cannot be negative")


@dataclass(frozen=True)
class CandidateField:
    """A cross-run equivalence class of continuous fields."""

    candidate_id: str
    tier: str
    member_factor_ids: tuple[str, ...]
    loading_cosine: float
    field_correlation: float
    restart_fraction: float
    patient_fraction: float
    raw_status: str = "UNTESTED"
    composition_status: str = "NOT_TESTED"
    validation_status: str = "NOT_TESTED"
    diagnostics: Mapping[str, Any] = field(default_factory=dict)
    input_hash: str = ""
    field_hash: str = ""

    def __post_init__(self) -> None:
        if self.tier not in {"BOTH", "MNSF_ONLY", "SIGNED_ONLY"}:
            raise R04ContractError(f"invalid candidate tier: {self.tier}")
        for value in (self.loading_cosine, self.field_correlation):
            if not -1 <= value <= 1:
                raise R04ContractError("candidate similarities must be in [-1, 1]")
        for value in (self.restart_fraction, self.patient_fraction):
            if not 0 <= value <= 1:
                raise R04ContractError("candidate stability fractions must be in [0, 1]")


@dataclass(frozen=True)
class GateResult:
    """Machine-readable scientific or operational result."""

    status: str
    reasons: tuple[str, ...] = ()
    metrics: Mapping[str, float] = field(default_factory=dict)
    input_hash: str = ""
    candidate_hash: str = ""


def as_float_array(value: Any, *, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.isfinite(array).all():
        raise R04ContractError(f"{name} contains non-finite values")
    return array
