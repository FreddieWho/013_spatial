"""Load one molecule-only section for R-16 replication (symbols, not ENSG).

Does not use the Vanderbilt HVG-10k cache. Gene symbols come from the source
10x feature table so marker modules can be scored on each lineage independently.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Mapping

import h5py
import numpy as np
from scipy import sparse

from r04.io_contract import _decode_strings
from r04.loaders import load_section_from_row
from r04.types import SectionData


def _symbols_from_10x_h5(path: Path, gene_ids: tuple[str, ...]) -> tuple[str, ...]:
    with h5py.File(path, "r") as handle:
        feats = handle["matrix"]["features"]
        ids = _decode_strings(feats["id"][()] if "id" in feats else feats["name"][()])
        names = _decode_strings(feats["name"][()] if "name" in feats else feats["id"][()])
    lut = {}
    for i, n in zip(ids, names):
        lut.setdefault(i, n)
        lut.setdefault(n, n)
    return tuple(lut.get(g, g) for g in gene_ids)


def _symbols_from_10x_mtx(directory: Path, gene_ids: tuple[str, ...]) -> tuple[str, ...]:
    feature_path = next(
        (p for p in (
            directory / "features.tsv.gz", directory / "features.tsv",
            directory / "genes.tsv.gz", directory / "genes.tsv",
        ) if p.exists()),
        None,
    )
    if feature_path is None:
        return gene_ids
    opener = gzip.open if feature_path.suffix == ".gz" else open
    rows = []
    with opener(feature_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            rows.append(line.rstrip("\n").split("\t"))
    lut = {}
    for row in rows:
        ens = row[0]
        sym = row[1] if len(row) > 1 else row[0]
        lut.setdefault(ens, sym)
        lut.setdefault(sym, sym)
    return tuple(lut.get(g, g) for g in gene_ids)


def load_section_symbols(row: Mapping[str, str]) -> SectionData:
    """Like ``load_section_from_row`` but ``gene_id`` is gene symbols."""
    section = load_section_from_row(row)
    matrix_path = Path(row["matrix_locator"])
    kind = row.get("matrix_kind", "")
    if kind == "10x_h5" or matrix_path.suffix in {".h5", ".hdf5"}:
        symbols = _symbols_from_10x_h5(matrix_path, section.gene_id)
    elif kind == "10x_mtx" or matrix_path.is_dir():
        symbols = _symbols_from_10x_mtx(matrix_path, section.gene_id)
    else:
        symbols = section.gene_id
    return SectionData(
        section_id=section.section_id,
        patient_id=section.patient_id,
        block_id=section.block_id,
        lineage=section.lineage,
        barcode=section.barcode,
        coords=section.coords,
        counts=section.counts,
        gene_id=symbols,
        library_size=section.library_size,
        metadata=section.metadata,
    )


def usz_tls_labels(section_id: str, barcodes: tuple[str, ...], root: Path) -> np.ndarray | None:
    """Return TLS 0/1/nan aligned to barcodes from USZ preprocessed h5ad, or None."""
    alias = section_id.split("::")[-1]
    path = root / "data/other_sources/zenodo_usz_tls_visium/TLS_VISIUM_USZ/h5ad_preprocessed" / f"{alias}.h5ad"
    if not path.exists():
        return None
    with h5py.File(path, "r") as handle:
        idx = _decode_strings(handle["obs"]["_index"][()])
        cats = _decode_strings(handle["obs"]["ground_truth"]["categories"][()])
        codes = np.asarray(handle["obs"]["ground_truth"]["codes"][:], dtype=int)
    lut = {b: int(codes[i]) for i, b in enumerate(idx)}
    tls_code = next((i for i, c in enumerate(cats) if c.upper() == "TLS"), None)
    if tls_code is None:
        return None
    out = np.full(len(barcodes), np.nan)
    for i, b in enumerate(barcodes):
        if b not in lut:
            continue
        code = lut[b]
        label = cats[code] if 0 <= code < len(cats) else "UNASSIGNED"
        if label.upper() == "UNASSIGNED":
            continue
        out[i] = 1.0 if code == tls_code else 0.0
    return out
