"""Readers for local count matrices with an explicit, fail-fast axis contract.

The legacy reader in ``scripts/field_discovery_20261018.py`` silently transposed a
matrix whenever the shape looked transposed.  A square (or coincidentally matching)
matrix then permutes spots and genes; the addendum forbids that.  Here the axis
orientation must be declared, and every mismatch is an error.
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import h5py
import numpy as np
from scipy import sparse

from .observation import ContractError, MaskViolation

OBS_BY_GENE = "obs_by_gene"
GENE_BY_OBS = "gene_by_obs"
AXIS_ORDERS = (OBS_BY_GENE, GENE_BY_OBS)

# obs columns that are derived from the full expression matrix; they would leak
# hidden-region statistics if a model received them.
EXPRESSION_DERIVED_OBS_COLUMNS = (
    "n_genes_by_counts",
    "total_counts",
    "log1p_n_genes_by_counts",
    "log1p_total_counts",
    "log1p_total_counts_mito",
    "pct_counts_mito",
    "pct_counts_in_top_50_genes",
    "pct_counts_in_top_100_genes",
    "pct_counts_in_top_200_genes",
    "pct_counts_in_top_500_genes",
    "mean_counts",
    "n_cells_by_counts",
    "pct_dropout_by_counts",
)

DEFAULT_OBS_AUX = ("array_row", "array_col", "in_tissue")


@dataclass(frozen=True)
class AxisSpec:
    """Explicit orientation of the stored matrix."""

    n_obs: int | None = None
    n_genes: int | None = None
    axis_order: str = OBS_BY_GENE

    def __post_init__(self) -> None:
        if self.axis_order not in AXIS_ORDERS:
            raise ContractError(f"axis_order must be one of {AXIS_ORDERS}")


@dataclass
class ReadReport:
    path: str
    bytes: int
    encoding: str
    axis_order: str
    n_obs: int
    n_genes: int
    n_duplicate_symbols: int
    duplicate_policy: str
    duplicate_examples: list[str] = field(default_factory=list)
    zero_depth_spots: int = 0
    zero_depth_examples: list[str] = field(default_factory=list)
    aux_kept: list[str] = field(default_factory=list)
    aux_dropped: list[str] = field(default_factory=list)
    spatial_key: str | None = None
    negative_values: int = 0
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class CountTable:
    counts: sparse.csr_matrix
    genes: np.ndarray
    barcodes: np.ndarray
    xy: np.ndarray
    aux: dict[str, np.ndarray]
    report: ReadReport

    @property
    def n_obs(self) -> int:
        return self.counts.shape[0]

    @property
    def n_genes(self) -> int:
        return self.counts.shape[1]


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _decode(values: Sequence[Any]) -> np.ndarray:
    return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in values])


def _attr_str(node: h5py.Group | h5py.Dataset, key: str, default: str = "") -> str:
    if key not in node.attrs:
        return default
    value = node.attrs[key]
    return value.decode() if isinstance(value, bytes) else str(value)


def _sparse_from_group(node: h5py.Group, encoding: str) -> sparse.spmatrix:
    shape = tuple(int(v) for v in node.attrs.get("shape", ()))
    if len(shape) != 2:
        raise ContractError(f"sparse group without usable 'shape' attribute: {shape}")
    data = node["data"][:]
    indices = node["indices"][:]
    indptr = node["indptr"][:]
    if encoding == "csc_matrix":
        return sparse.csc_matrix((data, indices, indptr), shape=shape)
    if encoding == "csr_matrix":
        return sparse.csr_matrix((data, indices, indptr), shape=shape)
    # no encoding: derive orientation, then still validate against the declared axis
    n_rows = len(indptr) - 1
    if shape[0] == shape[1] and n_rows in (shape[0], shape[1]):
        raise ContractError(
            "ambiguous square sparse matrix without 'encoding-type'; "
            "declare the axis explicitly at the reader level"
        )
    if n_rows == shape[0]:
        return sparse.csr_matrix((data, indices, indptr), shape=shape)
    if n_rows == shape[1]:
        return sparse.csc_matrix((data, indices, indptr), shape=shape)
    raise ContractError("sparse indptr length matches neither axis")


def _load_matrix(node: h5py.Group | h5py.Dataset) -> tuple[sparse.csr_matrix, str]:
    if isinstance(node, h5py.Dataset):
        arr = np.asarray(node[:])
        return sparse.csr_matrix(arr), "dense"
    encoding = _attr_str(node, "encoding-type")
    matrix = _sparse_from_group(node, encoding)
    return matrix.tocsr(), encoding or "sparse_unlabelled"


def _axis_of(indptr_len: int, shape: tuple[int, int], encoding: str) -> str:
    n_rows = indptr_len
    if shape[0] == shape[1]:
        return OBS_BY_GENE if encoding == "csr_matrix" else GENE_BY_OBS
    if n_rows == shape[0]:
        return OBS_BY_GENE
    if n_rows == shape[1]:
        return GENE_BY_OBS
    raise ContractError("cannot determine axis from indptr length")


def scan_h5ad_header(path: str | Path) -> dict[str, Any]:
    """Metadata-only probe (no matrix payload) used for asset inventory."""
    path = Path(path)
    info: dict[str, Any] = {"path": str(path), "exists": path.exists()}
    if not path.exists():
        return info
    info["bytes"] = path.stat().st_size
    with h5py.File(path, "r") as f:
        x = f["X"]
        if isinstance(x, h5py.Dataset):
            info["encoding"] = "dense"
            info["shape"] = tuple(int(v) for v in x.shape)
        else:
            info["encoding"] = _attr_str(x, "encoding-type") or "unknown"
            info["shape"] = tuple(int(v) for v in x.attrs.get("shape", ()))
        info["n_barcodes"] = int(len(f["obs"]["_index"]))
        info["n_genes"] = int(len(f["var"]["_index"]))
        info["has_spatial"] = "spatial" in f.get("obsm", {})
        if info["has_spatial"]:
            info["spatial_shape"] = tuple(int(v) for v in f["obsm"]["spatial"].shape)
        info["obs_columns"] = sorted(f["obs"].keys())
        info["uns_keys"] = sorted(f["uns"].keys()) if "uns" in f else []
    shapes_match_obs = info.get("shape", (None,))[0] == info.get("n_barcodes")
    info["declared_axis_consistent"] = bool(
        shapes_match_obs and info.get("shape", (None, None))[-1] == info.get("n_genes")
    )
    return info


# --------------------------------------------------------------------------- #
# main reader
# --------------------------------------------------------------------------- #
def read_h5ad_counts(
    path: str | Path,
    *,
    axis: AxisSpec | None = None,
    on_duplicate: str = "error",
    obs_aux_columns: Sequence[str] = DEFAULT_OBS_AUX,
    spatial_key: str = "spatial",
    load_matrix: bool = True,
) -> CountTable:
    """Read a local .h5ad count matrix under the explicit-axis contract."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if on_duplicate not in ("error", "sum"):
        raise ContractError("on_duplicate must be 'error' or 'sum'")
    axis = axis or AxisSpec()

    with h5py.File(path, "r") as f:
        barcodes = _decode(f["obs"]["_index"][:])
        genes = _decode(f["var"]["_index"][:])
        n_obs_stored = len(barcodes)
        n_genes_stored = len(genes)

        if axis.n_obs is not None and axis.n_obs != n_obs_stored:
            raise ContractError(
                f"declared n_obs={axis.n_obs} but obs/_index has {n_obs_stored} rows"
            )
        if axis.n_genes is not None and axis.n_genes != n_genes_stored:
            raise ContractError(
                f"declared n_genes={axis.n_genes} but var/_index has {n_genes_stored} columns"
            )

        if "spatial" not in f.get("obsm", {}):
            raise MaskViolation("dataset has no obsm['spatial']; coordinates are required")
        xy = np.asarray(f["obsm"][spatial_key][:], dtype=np.float64)
        if xy.shape != (n_obs_stored, 2):
            raise ContractError(f"obsm['{spatial_key}'] shape {xy.shape} != ({n_obs_stored}, 2)")

        aux: dict[str, np.ndarray] = {}
        obs_keys = set(f["obs"].keys())
        for column in obs_aux_columns:
            if column in obs_keys:
                aux[column] = np.asarray(f["obs"][column][:])
        dropped = sorted(
            k for k in obs_keys if k != "_index" and k not in aux and k in EXPRESSION_DERIVED_OBS_COLUMNS
        )

        encoding = _attr_str(f["X"], "encoding-type") if isinstance(f["X"], h5py.Group) else "dense"
        matrix = None
        if load_matrix:
            matrix, encoding = _load_matrix(f["X"])
            stored_shape = matrix.shape
        else:
            stored_shape = (
                tuple(int(v) for v in f["X"].attrs.get("shape", ()))
                if isinstance(f["X"], h5py.Group)
                else tuple(int(v) for v in f["X"].shape)
            )

    expected = (n_obs_stored, n_genes_stored)
    stored = tuple(stored_shape)
    if axis.axis_order == OBS_BY_GENE:
        if stored != expected:
            if stored == expected[::-1]:
                raise ContractError(
                    f"matrix is stored transposed: shape {stored} vs obs/var {expected}. "
                    "Refusing to guess; re-read with AxisSpec(axis_order='gene_by_obs') if that is intended."
                )
            raise ContractError(f"stored shape {stored} does not match obs/var {expected}")
    else:  # declared gene_by_obs
        if stored != expected[::-1]:
            raise ContractError(
                f"declared gene_by_obs but stored shape {stored} != {expected[::-1]}"
            )

    report = ReadReport(
        path=str(path),
        bytes=path.stat().st_size,
        encoding=encoding,
        axis_order=axis.axis_order,
        n_obs=n_obs_stored,
        n_genes=n_genes_stored,
        n_duplicate_symbols=int(len(genes) - len(set(genes))),
        duplicate_policy=on_duplicate,
        duplicate_examples=[
            g for g, c in Counter(genes.tolist()).most_common() if c > 1
        ][:5],
        aux_kept=sorted(aux),
        aux_dropped=dropped,
        spatial_key=spatial_key,
        notes=[],
    )

    if matrix is None:
        return CountTable(
            counts=sparse.csr_matrix((n_obs_stored, n_genes_stored)),
            genes=genes,
            barcodes=barcodes,
            xy=xy,
            aux=aux,
            report=report,
        )

    if axis.axis_order == GENE_BY_OBS:
        matrix = matrix.T.tocsr()

    if not np.all(np.isfinite(matrix.data)):
        raise ContractError("count matrix contains non-finite values")
    negative = int((matrix.data < 0).sum())
    if negative:
        raise ContractError(f"count matrix contains {negative} negative values")
    report.negative_values = negative

    if report.n_duplicate_symbols:
        if on_duplicate == "error":
            raise ContractError(
                f"{report.n_duplicate_symbols} duplicated gene symbols, e.g. "
                f"{report.duplicate_examples}; pass on_duplicate='sum' to aggregate explicitly"
            )
        genes, matrix = _sum_duplicates(genes, matrix)
        report.notes.append(f"duplicates summed; gene axis reduced to {len(genes)}")

    depth = np.asarray(matrix.sum(axis=1)).ravel()
    zero = np.flatnonzero(depth <= 0)
    report.zero_depth_spots = int(len(zero))
    report.zero_depth_examples = [str(barcodes[i]) for i in zero[:5]]
    if len(zero):
        report.notes.append("zero-depth spots retained (not dropped); they are not biological negatives")

    return CountTable(
        counts=matrix,
        genes=genes,
        barcodes=barcodes,
        xy=xy,
        aux=aux,
        report=report,
    )


def _sum_duplicates(genes: np.ndarray, matrix: sparse.csr_matrix) -> tuple[np.ndarray, sparse.csr_matrix]:
    order: list[str] = []
    index: dict[str, int] = {}
    for g in genes:
        if g not in index:
            index[g] = len(order)
            order.append(g)
    mapping = np.array([index[g] for g in genes], dtype=np.int64)
    selection = sparse.csr_matrix(
        (np.ones(len(genes), dtype=np.float64), (mapping, np.arange(len(genes)))),
        shape=(len(order), len(genes)),
    )
    merged = (matrix @ selection.T).tocsr()
    return np.array(order), merged


def write_read_report(report: ReadReport, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report.to_dict(), indent=2, ensure_ascii=False))
    return path
