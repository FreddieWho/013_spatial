"""Observation bundle, gene roles and mask-enforcing statistics.

Design rules (FA_CONTRACT_ADDENDUM_20261020.md):
  * masked values are NaN and live in a separate file from the model-visible bundle;
  * every statistic is a function of (values, visible_mask) and refuses to run unmasked;
  * a witness (B) gene may never appear on any input layer's gene axis.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from scipy import sparse

ROLE_A_FEATURE = "A_FEATURE"
ROLE_B_WITNESS = "B_WITNESS"
ROLE_BACKGROUND = "BACKGROUND"
ROLE_EXCLUDED = "EXCLUDED"
ROLES = (ROLE_A_FEATURE, ROLE_B_WITNESS, ROLE_BACKGROUND, ROLE_EXCLUDED)

INPUT_ROLES = (ROLE_A_FEATURE, ROLE_BACKGROUND)


class MaskViolation(AssertionError):
    """Raised when hidden information could reach a model input."""


class ContractError(ValueError):
    """Raised when an observation bundle violates the observation contract."""


# --------------------------------------------------------------------------- #
# gene roles
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class GeneRoles:
    """Explicit gene -> role registry.

    A boolean "witness mode" switch is not enough to keep the B genes out of the
    input layers, so the registry is a required object, not an option.
    """

    roles: Mapping[str, str]

    def __post_init__(self) -> None:
        bad = {g: r for g, r in self.roles.items() if r not in ROLES}
        if bad:
            raise ContractError(f"unknown gene roles: {sorted(bad.items())[:5]}")

    @classmethod
    def from_mapping(cls, mapping: Mapping[str, str]) -> "GeneRoles":
        return cls(dict(mapping))

    @property
    def witness_genes(self) -> set[str]:
        return {g for g, r in self.roles.items() if r == ROLE_B_WITNESS}

    def role_of(self, gene: str, default: str = ROLE_EXCLUDED) -> str:
        return self.roles.get(str(gene), default)

    def mask(self, genes: Sequence[str], allowed: Iterable[str] = INPUT_ROLES) -> np.ndarray:
        allowed = set(allowed)
        return np.array([self.role_of(str(g)) in allowed for g in genes], dtype=bool)

    def assert_no_witness(self, genes: Sequence[str], context: str = "") -> None:
        present = sorted(self.witness_genes.intersection(map(str, genes)))
        if present:
            where = f" ({context})" if context else ""
            raise MaskViolation(
                f"{len(present)} witness gene(s) reached a model input{where}: {present[:8]}"
            )


# --------------------------------------------------------------------------- #
# arrays / fingerprints
# --------------------------------------------------------------------------- #
def _as_2d_float(values: Any, name: str = "values") -> np.ndarray:
    arr = np.asarray(values)
    if arr.ndim != 2:
        raise ContractError(f"{name} must be 2-D, got shape {arr.shape}")
    return arr


def fingerprint_np(arr: Any) -> str:
    """Content hash of an array: dtype, shape and raw bytes."""
    a = np.ascontiguousarray(np.asarray(arr))
    digest = hashlib.sha256()
    digest.update(str(a.dtype).encode())
    digest.update(repr(a.shape).encode())
    digest.update(a.tobytes())
    return digest.hexdigest()


def fingerprint_arrays(mapping: Mapping[str, Any]) -> dict[str, str]:
    out: dict[str, str] = {}
    for key in sorted(mapping):
        value = mapping[key]
        out[key] = fingerprint_np(value)
    return out


def max_abs_delta(a: Any, b: Any) -> float:
    """max|a-b| with NaN treated as a value (NaN != NaN counts as a difference)."""
    x = np.asarray(a)
    y = np.asarray(b)
    if x.shape != y.shape:
        return float("inf")
    if x.dtype == bool or y.dtype == bool:
        return 0.0 if np.array_equal(x, y) else float("inf")
    if x.dtype != y.dtype:
        x = x.astype(np.float64)
        y = y.astype(np.float64)
    nan_x = np.isnan(x) if np.issubdtype(x.dtype, np.floating) else np.zeros(x.shape, bool)
    nan_y = np.isnan(y) if np.issubdtype(y.dtype, np.floating) else np.zeros(y.shape, bool)
    if not np.array_equal(nan_x, nan_y):
        return float("inf")
    both = ~nan_x
    if not both.any():
        return 0.0
    return float(np.max(np.abs(x[both] - y[both])))


# --------------------------------------------------------------------------- #
# masked statistics (the only sanctioned way to summarise observations)
# --------------------------------------------------------------------------- #
def _check_mask(mask: Any, n_obs: int) -> np.ndarray:
    m = np.asarray(mask, dtype=bool)
    if m.shape != (n_obs,):
        raise ContractError(f"visible_mask must have shape ({n_obs},), got {m.shape}")
    if not m.any():
        raise ContractError("visible_mask is empty; nothing is observable")
    return m


def visible_mean_std(values: Any, visible_mask: Any) -> tuple[np.ndarray, np.ndarray]:
    """Per-gene mean/std over *visible rows only*. Hidden rows never contribute."""
    if sparse.issparse(values):
        values = values.toarray()
    x = _as_2d_float(values)
    m = _check_mask(visible_mask, x.shape[0])
    sub = x[m]
    if not np.all(np.isfinite(sub)):
        raise MaskViolation("visible rows contain non-finite values; mask was not applied first")
    return np.nanmean(sub, axis=0), np.nanstd(sub, axis=0)


def visible_library_size(
    counts: Any, visible_mask: Any, gene_mask: Sequence[bool] | None = None
) -> np.ndarray:
    """Library size per spot over *feature genes only*, and only on visible spots.

    Hidden spots get NaN. The denominator must never include witness genes
    (otherwise B enters A algebraically through the shared library).
    """
    m = np.asarray(visible_mask, dtype=bool)
    n_obs = counts.shape[0]
    m = _check_mask(m, n_obs)
    if gene_mask is not None:
        gm = np.asarray(gene_mask, dtype=bool)
        if gm.shape != (counts.shape[1],):
            raise ContractError("gene_mask length does not match gene axis")
    else:
        gm = None
    out = np.full(n_obs, np.nan, dtype=np.float64)
    if sparse.issparse(counts):
        sub = counts[:, gm].tocsr() if gm is not None else sparse.csr_matrix(counts)
        values = np.asarray(sub[m].sum(axis=1)).ravel()
        if not np.all(np.isfinite(values)):
            raise ContractError("non-finite library size on visible rows")
        out[m] = values
        return out
    c = _as_2d_float(counts)
    if gm is not None:
        c = c[:, gm]
    if not np.all(np.isfinite(c[m])):
        raise MaskViolation("visible rows contain non-finite values; mask was not applied first")
    out[m] = c[m].sum(axis=1)
    return out


def visible_detection_rate(
    counts: Any, visible_mask: Any, gene_mask: Sequence[bool] | None = None
) -> np.ndarray:
    """Fraction of visible feature genes detected per visible spot; NaN elsewhere."""
    m = np.asarray(visible_mask, dtype=bool)
    n_obs = counts.shape[0]
    m = _check_mask(m, n_obs)
    gm = None
    if gene_mask is not None:
        gm = np.asarray(gene_mask, dtype=bool)
        if gm.shape != (counts.shape[1],):
            raise ContractError("gene_mask length does not match gene axis")
    out = np.full(n_obs, np.nan, dtype=np.float64)
    if sparse.issparse(counts):
        sub = counts[:, gm].tocsr() if gm is not None else sparse.csr_matrix(counts)
        denom = max(sub.shape[1], 1)
        out[m] = np.asarray((sub[m] > 0).sum(axis=1)).ravel() / denom
        return out
    c = _as_2d_float(counts)
    if gm is not None:
        c = c[:, gm]
    denom = max(c.shape[1], 1)
    out[m] = (c[m] > 0).sum(axis=1) / denom
    return out


def log_normalize(counts: Any, library: Any, scale: float = 1e4) -> np.ndarray:
    """log1p(counts * scale / library) with NaN for hidden rows."""
    c = np.asarray(counts, dtype=np.float64)
    lib = np.asarray(library, dtype=np.float64).reshape(-1, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        factor = np.where(np.isfinite(lib) & (lib > 0), scale / np.maximum(lib, 1.0), np.nan)
    return np.log1p(c * factor)


# --------------------------------------------------------------------------- #
# bundle
# --------------------------------------------------------------------------- #
@dataclass
class ObservationBundle:
    """Everything a model is allowed to see about one section."""

    sample_id: str
    unit_id: str
    source_id: str
    patient_status: str
    barcodes: np.ndarray
    source_rows: np.ndarray
    xy: np.ndarray
    coordinate_unit: str
    scale_provenance: Mapping[str, Any]
    tissue_support: np.ndarray
    visible_mask: np.ndarray
    genes: np.ndarray
    gene_roles: GeneRoles
    values_visible: np.ndarray
    field_ids: Sequence[str] = field(default_factory=tuple)
    gene_sets: Mapping[str, Sequence[str]] = field(default_factory=dict)
    score_definition: Mapping[str, Any] = field(default_factory=dict)
    auxiliary_visible: Mapping[str, np.ndarray] = field(default_factory=dict)
    input_provenance: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        n_obs = len(self.barcodes)
        self.barcodes = np.asarray(self.barcodes).astype(str)
        self.genes = np.asarray(self.genes).astype(str)
        self.source_rows = np.asarray(self.source_rows)
        self.xy = np.asarray(self.xy, dtype=np.float64)
        self.tissue_support = np.asarray(self.tissue_support, dtype=bool)
        self.visible_mask = _check_mask(self.visible_mask, n_obs)
        self.values_visible = _as_2d_float(self.values_visible, "values_visible")

        if self.xy.shape != (n_obs, 2):
            raise ContractError(f"xy must be ({n_obs}, 2), got {self.xy.shape}")
        if self.source_rows.shape != (n_obs,):
            raise ContractError("source_rows length does not match barcodes")
        if self.tissue_support.shape != (n_obs,):
            raise ContractError("tissue_support length does not match barcodes")
        if self.values_visible.shape[0] != n_obs:
            raise ContractError("values_visible row count does not match barcodes")
        if self.values_visible.shape[1] != len(self.genes):
            raise ContractError("values_visible column count does not match genes")
        if not np.all(np.isfinite(self.xy)):
            raise ContractError("coordinates must be finite")

        hidden = ~self.visible_mask
        if hidden.any() and not np.all(np.isnan(self.values_visible[hidden])):
            raise MaskViolation("hidden rows must be NaN in the model-visible bundle")
        if not np.all(np.isfinite(self.values_visible[self.visible_mask])):
            raise ContractError("visible rows must be finite (no NaN inside the observed region)")

        for key, aux in self.auxiliary_visible.items():
            arr = np.asarray(aux)
            if arr.shape != (n_obs,):
                raise ContractError(f"auxiliary '{key}' must have shape ({n_obs},), got {arr.shape}")
            if hidden.any() and not np.all(np.isnan(arr[hidden])):
                raise MaskViolation(f"auxiliary '{key}' carries hidden values")

        self.gene_roles.assert_no_witness(self.genes, "ObservationBundle.genes")

    # -- convenience ------------------------------------------------------- #
    @property
    def n_obs(self) -> int:
        return len(self.barcodes)

    @property
    def n_genes(self) -> int:
        return len(self.genes)

    @property
    def n_visible(self) -> int:
        return int(self.visible_mask.sum())

    def hidden_rows(self) -> np.ndarray:
        return np.flatnonzero(~self.visible_mask)

    def summary(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "unit_id": self.unit_id,
            "source_id": self.source_id,
            "patient_status": self.patient_status,
            "n_obs": self.n_obs,
            "n_visible": self.n_visible,
            "n_genes": self.n_genes,
            "coordinate_unit": self.coordinate_unit,
            "field_ids": list(self.field_ids),
        }


# --------------------------------------------------------------------------- #
# split / persistence: hidden truth never enters the bundle
# --------------------------------------------------------------------------- #
@dataclass
class HiddenTruth:
    """Evaluation-only container. Trainers must never receive this object."""

    values: np.ndarray
    rows: np.ndarray
    genes: np.ndarray
    evaluation_mask: np.ndarray
    meta: Mapping[str, Any] = field(default_factory=dict)

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            path,
            values=self.values,
            rows=self.rows,
            genes=self.genes,
            evaluation_mask=self.evaluation_mask,
        )
        meta_path = path.with_suffix(".meta.json")
        meta_path.write_text(json.dumps(dict(self.meta), indent=2, ensure_ascii=False))
        return path


def split_for_evaluation(
    values_full: Any, visible_mask: Any, *, keep="hidden"
) -> tuple[np.ndarray, HiddenTruth]:
    """Split full values into (visible-with-NaN, HiddenTruth).

    ``values_full`` is the *raw* matrix; this is the only function allowed to see
    the hidden rows, and it hands them to a separate object that is written to a
    separate file.
    """
    x = _as_2d_float(values_full, "values_full")
    m = _check_mask(visible_mask, x.shape[0])
    visible = np.array(x, dtype=np.float64, copy=True)
    hidden_rows = np.flatnonzero(~m)
    visible[~m] = np.nan
    if keep not in ("hidden", "all"):
        raise ContractError("keep must be 'hidden' or 'all'")
    rows = hidden_rows if keep == "hidden" else np.arange(x.shape[0])
    truth = HiddenTruth(
        values=np.array(x[rows], dtype=np.float64, copy=True),
        rows=rows.astype(np.int64),
        genes=np.array([], dtype=str),
        evaluation_mask=np.ones(len(rows), dtype=bool),
        meta={"keep": keep, "n_rows": int(len(rows))},
    )
    return visible, truth


def _encode_strings(arr: Any) -> np.ndarray:
    return np.asarray(arr).astype(str)


def save_observation(bundle: ObservationBundle, directory: str | Path, stem: str = "observation") -> Path:
    """Write the model-visible bundle. No hidden array is ever written here."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{stem}.npz"
    payload = {
        "barcodes": _encode_strings(bundle.barcodes),
        "source_rows": np.asarray(bundle.source_rows),
        "xy": np.asarray(bundle.xy, dtype=np.float64),
        "tissue_support": np.asarray(bundle.tissue_support, dtype=bool),
        "visible_mask": np.asarray(bundle.visible_mask, dtype=bool),
        "genes": _encode_strings(bundle.genes),
        "values_visible": np.asarray(bundle.values_visible, dtype=np.float64),
        "coordinate_unit": np.array(bundle.coordinate_unit),
    }
    for key, aux in bundle.auxiliary_visible.items():
        payload[f"aux::{key}"] = np.asarray(aux, dtype=np.float64)
    np.savez_compressed(path, **payload)

    sidecar = {
        "summary": bundle.summary(),
        "scale_provenance": dict(bundle.scale_provenance),
        "scoring": {
            "field_ids": list(bundle.field_ids),
            "gene_sets": {k: list(v) for k, v in bundle.gene_sets.items()},
            "score_definition": dict(bundle.score_definition),
        },
        "input_provenance": dict(bundle.input_provenance),
        "gene_roles": dict(bundle.gene_roles.roles),
        "fingerprints": fingerprint_arrays(
            {
                "values_visible": np.asarray(bundle.values_visible, dtype=np.float64),
                "visible_mask": np.asarray(bundle.visible_mask, dtype=bool),
                "xy": np.asarray(bundle.xy, dtype=np.float64),
                "tissue_support": np.asarray(bundle.tissue_support, dtype=bool),
            }
        ),
        "contains_hidden_values": False,
    }
    (directory / f"{stem}.json").write_text(json.dumps(sidecar, indent=2, ensure_ascii=False))
    return path


def load_observation(directory: str | Path, stem: str = "observation") -> tuple[ObservationBundle, dict]:
    """Load a model-visible bundle, refusing any payload that smells hidden."""
    directory = Path(directory)
    npz = np.load(directory / f"{stem}.npz", allow_pickle=False)
    for key in npz.files:
        if key.startswith("hidden") or key.startswith("eval::"):
            raise MaskViolation(f"bundle file contains evaluation-only key '{key}'")
    sidecar = json.loads((directory / f"{stem}.json").read_text())
    aux = {k[len("aux::"):]: npz[k] for k in npz.files if k.startswith("aux::")}
    bundle = ObservationBundle(
        sample_id=sidecar["summary"]["sample_id"],
        unit_id=sidecar["summary"]["unit_id"],
        source_id=sidecar["summary"]["source_id"],
        patient_status=sidecar["summary"]["patient_status"],
        barcodes=npz["barcodes"],
        source_rows=npz["source_rows"],
        xy=npz["xy"],
        coordinate_unit=str(npz["coordinate_unit"]),
        scale_provenance=sidecar["scale_provenance"],
        tissue_support=npz["tissue_support"],
        visible_mask=npz["visible_mask"],
        genes=npz["genes"],
        gene_roles=GeneRoles.from_mapping(sidecar["gene_roles"]),
        values_visible=npz["values_visible"],
        field_ids=sidecar["scoring"]["field_ids"],
        gene_sets=sidecar["scoring"]["gene_sets"],
        score_definition=sidecar["scoring"]["score_definition"],
        auxiliary_visible=aux,
        input_provenance=sidecar["input_provenance"],
    )
    return bundle, sidecar
