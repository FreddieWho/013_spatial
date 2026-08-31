"""Materialise one section from an R-04 molecule-only manifest row."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Mapping

import numpy as np
from scipy import sparse

from .io_contract import build_section, read_10x_h5, read_10x_mtx, read_h5ad_counts
from .types import SectionData


def _read_text_table(path: Path) -> list[dict[str, str]]:
    opener = __import__("gzip").open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def read_10x_positions(path: Path) -> tuple[tuple[str, ...], np.ndarray]:
    rows = _read_text_table(path)
    if not rows:
        raise ValueError(f"empty coordinate table: {path}")
    fields = set(rows[0])
    barcode_key = next((key for key in ("barcode", "barcodes", "Barcode", "_index") if key in fields), None)
    coordinate_pair = next((pair for pair in (("array_row", "array_col"), ("row", "col"), ("x_array", "y_array")) if set(pair) <= fields), None)
    if barcode_key is None or coordinate_pair is None:
        # 10x tissue_positions_list.csv is a standard six-column, headerless
        # table: barcode, in_tissue, array_row, array_col, pixel_row,
        # pixel_col.  DictReader consumes the first data row as a header, so
        # reread it positionally when no named columns are present.
        opener = __import__("gzip").open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            positional = list(csv.reader(handle))
        if not positional or any(len(row) < 4 for row in positional):
            raise ValueError(f"coordinate table lacks barcode and array coordinate columns: {path}")
        selected = [row for row in positional if row[1] in {"1", "1.0", "True", "true"}]
        if not selected:
            raise ValueError(f"coordinate table has no in-tissue rows: {path}")
        return tuple(row[0] for row in selected), np.asarray([[float(row[2]), float(row[3])] for row in selected], dtype=float)
    in_tissue_key = next((key for key in ("in_tissue", "tissue") if key in fields), None)
    selected = [row for row in rows if in_tissue_key is None or row[in_tissue_key] in {"1", "1.0", "True", "true"}]
    barcodes = tuple(row[barcode_key] for row in selected)
    coords = np.asarray([[float(row[coordinate_pair[0]]), float(row[coordinate_pair[1]])] for row in selected], dtype=float)
    return barcodes, coords


def load_section_from_row(row: Mapping[str, str]) -> SectionData:
    matrix_path = Path(row["matrix_locator"])
    h5ad_coords = None
    cached_coords = None
    cached_library_size = None
    if row.get("matrix_kind") == "r04_panel_npz":
        metadata_locator = row.get("matrix_metadata_locator")
        if not metadata_locator:
            raise ValueError("r04_panel_npz rows require matrix_metadata_locator")
        metadata = json.loads(Path(metadata_locator).read_text(encoding="utf-8"))
        counts = sparse.load_npz(matrix_path).tocsr()
        barcodes = tuple(str(value) for value in metadata["barcodes"])
        genes = tuple(str(value) for value in metadata["gene_id"])
        if "coords" in metadata:
            cached_coords = np.asarray(metadata["coords"], dtype=float)
        if metadata.get("library_size") is not None:
            cached_library_size = np.asarray(metadata["library_size"], dtype=float)
        if counts.shape != (len(barcodes), len(genes)):
            raise ValueError(f"cached panel shape does not match metadata: {matrix_path}")
    elif matrix_path.is_dir():
        counts, barcodes, genes = read_10x_mtx(matrix_path)
    elif matrix_path.suffix in {".h5", ".hdf5"}:
        counts, barcodes, genes = read_10x_h5(matrix_path)
    elif matrix_path.suffix == ".h5ad":
        counts, barcodes, genes, h5ad_coords = read_h5ad_counts(matrix_path, count_key=row.get("count_key", "X"))
    else:
        raise ValueError(f"R04 loader accepts 10x H5/Matrix Market or allowlisted h5ad counts: {matrix_path}")
    if cached_coords is not None:
        coordinate_barcodes, coords = barcodes, cached_coords
    elif h5ad_coords is not None:
        coordinate_barcodes, coords = barcodes, h5ad_coords
    else:
        coordinate_barcodes, coords = read_10x_positions(Path(row["coordinate_locator"]))
    positions = {barcode: index for index, barcode in enumerate(coordinate_barcodes)}
    keep = [index for index, barcode in enumerate(barcodes) if barcode in positions]
    if not keep:
        raise ValueError(f"no barcode overlap between matrix and coordinates: {row['section_id']}")
    coordinate_index = np.asarray([positions[barcodes[index]] for index in keep], dtype=int)
    counts_kept = counts[keep]
    libraries = (
        cached_library_size
        if cached_library_size is not None
        else np.asarray(counts_kept.sum(axis=1)).ravel()
    )
    if len(libraries) != counts_kept.shape[0]:
        raise ValueError(f"library_size metadata does not match cached panel: {matrix_path}")
    nonzero = libraries > 0
    if not nonzero.any():
        raise ValueError(f"all barcode-overlap spots have zero library size: {row['section_id']}")
    excluded = int((~nonzero).sum())
    return build_section(
        section_id=row["section_id"], patient_id=row["patient_id"], block_id=row.get("block_id") or None,
        lineage=row["lineage"], counts=counts_kept[nonzero], barcode=[barcodes[index] for index in np.asarray(keep)[nonzero]],
        coords=coords[coordinate_index[nonzero]], gene_id=genes,
        library_size=libraries[nonzero],
        metadata={"zero_library_spots_excluded": str(excluded)},
    )
