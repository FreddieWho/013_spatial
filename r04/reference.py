"""Bounded local single-cell reference construction.

Only explicitly allowlisted counts, donor and cell-type fields are read.  No
outcome, response or treatment column is materialised by these adapters.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
from scipy import sparse


@dataclass(frozen=True)
class ReferenceSource:
    source_id: str
    count_locator: Path
    annotation_locator: Path | None
    donor_column: str
    cell_type_column: str
    provenance: str


@dataclass(frozen=True)
class ReferencePanel:
    counts: sparse.csr_matrix
    gene_id: tuple[str, ...]
    donor_id: tuple[str, ...]
    broad_type: tuple[str, ...]
    source_id: tuple[str, ...]
    uncertainty_status: str = "POINT_REFERENCE_NO_POSTERIOR"


TYPE_MAP = {
    "mal": "epithelial_malignant", "tumor": "epithelial_malignant", "epithelial": "epithelial_malignant",
    "caf": "fibroblast", "fibroblast": "fibroblast", "endo": "endothelial", "endothelial": "endothelial",
    "macrophage": "myeloid", "myeloid": "myeloid", "monocyte": "myeloid",
    "t.cd4": "T_CD4", "cd4": "T_CD4", "t.cd8": "T_CD8", "cd8": "T_CD8",
    "nk": "NK", "b.cell": "B_plasma", "b": "B_plasma", "plasma": "B_plasma",
}


def harmonize_cell_type(value: object) -> str:
    text = str(value).strip().lower()
    for token, broad in TYPE_MAP.items():
        if token in text:
            return broad
    return "other_unknown"


def deterministic_cap(
    donor: Iterable[str], broad_type: Iterable[str], *, per_donor_type: int = 200,
    max_total: int = 50_000,
) -> np.ndarray:
    donor = np.asarray(list(donor), dtype=str)
    broad_type = np.asarray(list(broad_type), dtype=str)
    if len(donor) != len(broad_type):
        raise ValueError("donor and cell type lengths differ")
    selected: list[int] = []
    for d in sorted(np.unique(donor)):
        for cell_type in sorted(np.unique(broad_type[donor == d])):
            indices = np.flatnonzero((donor == d) & (broad_type == cell_type))
            selected.extend(indices[:per_donor_type].tolist())
    return np.asarray(selected[:max_total], dtype=int)


def build_reference_panel(
    matrices: Iterable[sparse.spmatrix],
    genes: Iterable[str] | Iterable[Iterable[str]],
    donors: Iterable[Iterable[str]],
    cell_types: Iterable[Iterable[str]],
    source_ids: Iterable[str],
    *,
    spatial_genes: Iterable[str] | None = None,
    per_donor_type: int = 200,
    max_total: int = 50_000,
) -> ReferencePanel:
    matrices = list(matrices)
    gene_values = list(genes)
    if gene_values and isinstance(gene_values[0], str):
        gene_lists = [tuple(map(str, gene_values))] * len(matrices)
    else:
        gene_lists = [tuple(map(str, value)) for value in gene_values]
    donor_arrays = [np.asarray(list(value), dtype=str) for value in donors]
    type_arrays = [np.asarray([harmonize_cell_type(item) for item in value], dtype=str) for value in cell_types]
    source_ids = list(map(str, source_ids))
    if not (len(matrices) == len(donor_arrays) == len(type_arrays) == len(source_ids) == len(gene_lists)):
        raise ValueError("reference source arrays differ in number")
    common = set(gene_lists[0])
    for source_genes in gene_lists[1:]:
        common &= set(source_genes)
    if spatial_genes is not None:
        common &= set(map(str, spatial_genes))
    gene_tuple = tuple(gene for gene in gene_lists[0] if gene in common)
    if not gene_tuple:
        raise ValueError("reference sources have no common genes")
    chunks: list[sparse.csr_matrix] = []
    donor_out: list[str] = []
    type_out: list[str] = []
    source_out: list[str] = []
    for matrix, source_genes, donor, broad, source in zip(matrices, gene_lists, donor_arrays, type_arrays, source_ids):
        if matrix.shape[0] != len(donor) or len(donor) != len(broad) or matrix.shape[1] != len(source_genes):
            raise ValueError(f"reference dimensions differ for {source}")
        selected = deterministic_cap(donor, broad, per_donor_type=per_donor_type, max_total=max_total)
        columns = [source_genes.index(gene) for gene in gene_tuple]
        chunks.append(sparse.csr_matrix(matrix)[selected][:, columns])
        donor_out.extend(donor[selected].tolist())
        type_out.extend(broad[selected].tolist())
        source_out.extend([source] * len(selected))
    counts = sparse.vstack(chunks, format="csr")
    return ReferencePanel(counts, gene_tuple, tuple(donor_out), tuple(type_out), tuple(source_out))


def read_gse132465_annotation(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Read only the six annotation columns needed by the CRC adapter."""
    import pandas as pd
    frame = pd.read_csv(path, sep="\t", compression="infer", dtype=str)
    required = {"Index", "Patient", "Cell_type"}
    if not required <= set(frame.columns):
        raise ValueError("GSE132465 annotation lacks Index, Patient or Cell_type")
    donor = dict(zip(frame["Index"], frame["Patient"]))
    cell_type = dict(zip(frame["Index"], frame["Cell_type"]))
    return donor, cell_type


def read_gse132465_counts(
    path: Path,
    annotation_path: Path,
    *,
    spatial_genes: Iterable[str] | None = None,
    per_donor_type: int = 200,
    max_total: int = 50_000,
) -> tuple[sparse.csr_matrix, tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    """Select a bounded cell subset before materialising the count matrix."""
    import pandas as pd
    donor_map, type_map = read_gse132465_annotation(annotation_path)
    header = pd.read_csv(path, sep="\t", compression="infer", nrows=0)
    columns = list(header.columns)
    gene_column = columns[0]
    cells = [cell for cell in columns[1:] if cell in donor_map and cell in type_map]
    donor = np.asarray([donor_map[cell] for cell in cells], dtype=str)
    broad = np.asarray([harmonize_cell_type(type_map[cell]) for cell in cells], dtype=str)
    selected = deterministic_cap(donor, broad, per_donor_type=per_donor_type, max_total=max_total)
    selected_cells = [cells[index] for index in selected]
    wanted_genes = set(map(str, spatial_genes)) if spatial_genes is not None else None
    usecols = [gene_column, *selected_cells]
    frame = pd.read_csv(path, sep="\t", compression="infer", usecols=usecols)
    genes = np.asarray(frame[gene_column].astype(str))
    keep = np.ones(len(genes), dtype=bool) if wanted_genes is None else np.asarray([gene in wanted_genes for gene in genes])
    counts = sparse.csr_matrix(frame.loc[keep, selected_cells].to_numpy(dtype=np.int32).T)
    return counts, tuple(genes[keep]), tuple(donor[selected].tolist()), tuple(broad[selected].tolist())


def _decode_h5_values(values: object) -> np.ndarray:
    array = np.asarray(values)
    if array.dtype.kind == "S":
        return np.asarray([value.decode() for value in array], dtype=str)
    if array.dtype.kind == "O":
        return np.asarray([value.decode() if isinstance(value, bytes) else str(value) for value in array], dtype=str)
    return array


def _read_h5ad_obs_column(obs: object, name: str) -> np.ndarray:
    import h5py
    if name not in obs:
        raise ValueError(f"allowlisted h5ad observation column is absent: {name}")
    value = obs[name]
    if isinstance(value, h5py.Dataset):
        return _decode_h5_values(value[()])
    if "categories" not in value or "codes" not in value:
        raise ValueError(f"unsupported h5ad categorical encoding for {name}")
    categories = _decode_h5_values(value["categories"][()])
    codes = np.asarray(value["codes"][()], dtype=int)
    return np.asarray([categories[code] if code >= 0 else "other_unknown" for code in codes], dtype=str)


def _read_h5ad_matrix(handle: object, key: str, rows: np.ndarray, columns: np.ndarray) -> sparse.csr_matrix:
    import h5py
    matrix = handle[key]
    n_obs = len(handle["obs"]["_index"])
    n_var = len(handle["var"]["_index"])
    if isinstance(matrix, h5py.Group):
        data = matrix["data"][()]
        indices = matrix["indices"][()]
        indptr = matrix["indptr"][()]
        full = sparse.csr_matrix((data, indices, indptr), shape=(n_obs, n_var))
        return full[rows][:, columns].tocsr()
    chunks: list[sparse.csr_matrix] = []
    for start in range(0, len(rows), 256):
        selected_rows = rows[start:start + 256]
        row_order = np.argsort(selected_rows)
        sorted_rows = selected_rows[row_order]
        dense = np.asarray(matrix[sorted_rows, :], dtype=np.int32)[:, columns]
        dense = dense[np.argsort(row_order)]
        chunks.append(sparse.csr_matrix(dense))
    return sparse.vstack(chunks, format="csr")


def read_allowlisted_h5ad(
    path: Path,
    *,
    count_key: str,
    donor_column: str,
    cell_type_column: str,
    source_id: str,
    spatial_genes: Iterable[str] | None = None,
    per_donor_type: int = 200,
    max_total: int = 50_000,
) -> tuple[sparse.csr_matrix, tuple[str, ...], tuple[str, ...], tuple[str, ...], str]:
    """Read a compact panel from an h5ad using a strict field allowlist.

    ``count_key`` may be ``X`` or ``layers/counts`` only.  The adapter does
    not touch any other obs field, including response and treatment columns.
    """
    import h5py
    allowed_count_keys = {"X", "layers/counts"}
    if count_key not in allowed_count_keys:
        raise ValueError(f"count key is outside the R04 allowlist: {count_key}")
    key = count_key.replace("/", "/")
    with h5py.File(path, "r") as handle:
        var_names = _decode_h5_values(handle["var"]["_index"][()])
        donors = _read_h5ad_obs_column(handle["obs"], donor_column)
        raw_types = _read_h5ad_obs_column(handle["obs"], cell_type_column)
        broad = np.asarray([harmonize_cell_type(value) for value in raw_types], dtype=str)
        selected = deterministic_cap(donors, broad, per_donor_type=per_donor_type, max_total=max_total)
        if spatial_genes is None:
            columns = np.arange(len(var_names), dtype=int)
        else:
            wanted = set(map(str, spatial_genes))
            columns = np.asarray([index for index, name in enumerate(var_names) if str(name) in wanted], dtype=int)
        if len(columns) == 0:
            raise ValueError(f"no spatial genes overlap h5ad reference: {path}")
        matrix = _read_h5ad_matrix(handle, key, selected, columns)
        genes = tuple(str(var_names[index]) for index in columns)
    return matrix, genes, tuple(donors[selected]), tuple(broad[selected]), source_id
