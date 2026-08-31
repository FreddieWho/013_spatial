"""Fail-closed molecule-only input readers for R-04."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
from scipy import sparse
from scipy.io import mmread

from .types import R04ContractError, SectionData


FORBIDDEN_INPUT_TOKENS = (
    "ground_truth", "gt", "tls", "tumor_stroma", "boundary", "vessel",
    "blood_vessel", "necrosis", "mask", "distance", "outcome", "response",
    "treatment", "stage", "site", "sample_name", "file_name", "hande",
    "survival", "pfs", "overall_survival", "efs", "dfs", "recurrence", "relapse",
    "clinical", "diagnosis", "tnm", "response",
)


@dataclass(frozen=True)
class ManifestRow:
    section_id: str
    patient_id: str
    block_id: str
    lineage: str
    matrix_locator: str
    coordinate_locator: str
    matrix_kind: str
    read_authorization: str = "MOLECULES_AND_COORDINATES"
    count_key: str = "X"
    outer_fold: str = ""
    leakage_group_id: str = ""
    primary_role: str = ""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_json_hash(value: object) -> str:
    encoded = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def section_content_hash(section: SectionData) -> str:
    """Hash identities, coordinates and count payload—not only shape metadata."""
    digest = hashlib.sha256()
    for value in (section.section_id, section.patient_id, section.block_id or "", section.lineage):
        digest.update(str(value).encode("utf-8"))
        digest.update(b"\0")
    for barcode in section.barcode:
        digest.update(barcode.encode("utf-8"))
        digest.update(b"\0")
    for gene in section.gene_id:
        digest.update(gene.encode("utf-8"))
        digest.update(b"\0")
    digest.update(np.asarray(section.coords, dtype=np.float64, order="C").tobytes())
    if sparse.issparse(section.counts):
        matrix = section.counts.tocsr()
        for value in (matrix.data, matrix.indices, matrix.indptr):
            digest.update(np.asarray(value).tobytes())
    else:
        digest.update(np.asarray(section.counts).tobytes())
    return digest.hexdigest()


def sections_content_hash(sections: Iterable[SectionData]) -> str:
    digest = hashlib.sha256()
    for section in sections:
        digest.update(section_content_hash(section).encode("ascii"))
    return digest.hexdigest()


def read_manifest(path: Path) -> list[ManifestRow]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"section_id", "patient_id", "block_id", "lineage", "matrix_locator", "coordinate_locator", "matrix_kind"}
    if not rows or not required <= set(rows[0]):
        raise R04ContractError("R04 manifest is missing required molecule-only columns")
    result = [ManifestRow(**{key: row.get(key, "") for key in ManifestRow.__dataclass_fields__}) for row in rows]
    if len({row.section_id for row in result}) != len(result):
        raise R04ContractError("section_id must be unique in R04 manifest")
    for row in result:
        if row.read_authorization not in {"MOLECULES_AND_COORDINATES", "HTAN_COUNTS_AND_COORDINATES"}:
            raise R04ContractError(f"unsupported read authorization for {row.section_id}")
        if row.matrix_kind == "h5ad_counts" and row.read_authorization != "HTAN_COUNTS_AND_COORDINATES":
            raise R04ContractError("h5ad expression requires explicit HTAN count authorization")
        if row.matrix_kind != "h5ad_counts" and row.read_authorization != "MOLECULES_AND_COORDINATES":
            raise R04ContractError("non-h5ad input requires molecule-only authorization")
        # Locator paths are provenance, not model covariates; inspect the
        # semantic identity/role fields but do not reject a harmless directory
        # name containing a target token.
        for key in ("section_id", "patient_id", "block_id", "lineage", "matrix_kind", "read_authorization"):
            value = getattr(row, key)
            if any(token in str(value).lower() for token in FORBIDDEN_INPUT_TOKENS):
                raise R04ContractError(f"forbidden target-like token in manifest row {row.section_id}")
    return result


def validate_count_matrix(counts: object) -> None:
    values = counts.data if sparse.issparse(counts) else np.asarray(counts)
    if np.any(~np.isfinite(values)) or np.any(values < 0):
        raise R04ContractError("counts must be finite and non-negative")
    if np.any(values != np.floor(values)):
        raise R04ContractError("R04 discovery requires integer raw counts")


def _decode_strings(values: object) -> tuple[str, ...]:
    return tuple(value.decode() if isinstance(value, bytes) else str(value) for value in values)


def read_10x_h5(path: Path) -> tuple[sparse.csr_matrix, tuple[str, ...], tuple[str, ...]]:
    """Read a 10x HDF5 matrix without loading images or annotations."""
    import h5py

    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        shape = tuple(int(v) for v in matrix["shape"][()])
        values = matrix["data"][()]
        indices = matrix["indices"][()]
        indptr = matrix["indptr"][()]
        gene_group = matrix["features"]
        gene_key = "id" if "id" in gene_group else "name"
        genes = _decode_strings(gene_group[gene_key][()])
        barcodes = _decode_strings(matrix["barcodes"][()])
    # 10x stores gene-by-barcode CSC; R04 uses spot-by-gene CSR.
    values = sparse.csc_matrix((values, indices, indptr), shape=shape)
    counts = values.T.tocsr()
    validate_count_matrix(counts)
    return counts, barcodes, genes


def read_10x_mtx(directory: Path) -> tuple[sparse.csr_matrix, tuple[str, ...], tuple[str, ...]]:
    matrix_path = directory / "matrix.mtx"
    if not matrix_path.exists():
        matrix_path = directory / "matrix.mtx.gz"
    feature_path = next((p for p in (directory / "features.tsv.gz", directory / "features.tsv", directory / "genes.tsv.gz", directory / "genes.tsv") if p.exists()), None)
    barcode_path = next((p for p in (directory / "barcodes.tsv.gz", directory / "barcodes.tsv") if p.exists()), None)
    if not matrix_path or not feature_path or not barcode_path:
        raise FileNotFoundError(f"incomplete 10x Matrix Market directory: {directory}")
    matrix = sparse.csc_matrix(mmread(matrix_path))
    feature_open = gzip.open if feature_path.suffix == ".gz" else open
    barcode_open = gzip.open if barcode_path.suffix == ".gz" else open
    with feature_open(feature_path, "rb") as handle:
        feature_rows = [line.decode().rstrip("\n").split("\t") for line in handle]
    with barcode_open(barcode_path, "rb") as handle:
        barcodes = tuple(line.decode().rstrip("\n") for line in handle)
    genes = tuple(row[0] for row in feature_rows)
    if matrix.shape != (len(genes), len(barcodes)):
        raise R04ContractError("10x matrix dimensions do not match feature/barcode files")
    counts = matrix.T.tocsr()
    validate_count_matrix(counts)
    return counts, barcodes, genes


def read_h5ad_counts(path: Path, *, count_key: str = "X") -> tuple[sparse.csr_matrix, tuple[str, ...], tuple[str, ...], np.ndarray]:
    """Read HTAN-style count h5ad content without loading obs annotations."""
    import h5py
    if count_key not in {"X", "layers/counts"}:
        raise R04ContractError(f"h5ad count key is outside the R04 allowlist: {count_key}")
    with h5py.File(path, "r") as handle:
        barcodes = _decode_strings(handle["obs"]["_index"][()])
        gene_key = "gene_ids" if "gene_ids" in handle["var"] else "_index"
        genes = _decode_strings(handle["var"][gene_key][()])
        matrix = handle[count_key]
        if isinstance(matrix, h5py.Group):
            values = matrix["data"][()]
            indices = matrix["indices"][()]
            indptr = matrix["indptr"][()]
            counts = sparse.csr_matrix((values, indices, indptr), shape=(len(barcodes), len(genes)))
        else:
            # Validate before casting.  Dense h5ad matrices are frequently
            # stored as float32 even when they contain integer counts; an
            # early int64 cast would silently truncate normalized values.
            dense = np.asarray(matrix[()])
            validate_count_matrix(dense)
            counts = sparse.csr_matrix(dense.astype(np.int64, copy=False))
        if "spatial" not in handle["obsm"]:
            raise R04ContractError(f"h5ad has no explicit obsm/spatial coordinates: {path}")
        coords = np.asarray(handle["obsm"]["spatial"][()], dtype=float)
    if coords.ndim != 2 or coords.shape[0] != len(barcodes) or coords.shape[1] < 2:
        raise R04ContractError("h5ad obsm/spatial dimensions are invalid")
    coords = coords[:, :2]
    validate_count_matrix(counts)
    counts = counts.astype(np.int64, copy=False)
    return counts, barcodes, genes, coords


def build_section(
    *,
    section_id: str,
    patient_id: str,
    block_id: str | None,
    lineage: str,
    counts: object,
    barcode: Iterable[str],
    coords: np.ndarray,
    gene_id: Iterable[str],
    library_size: np.ndarray | None = None,
    metadata: Mapping[str, str] | None = None,
) -> SectionData:
    validate_count_matrix(counts)
    return SectionData(
        section_id=section_id, patient_id=patient_id, block_id=block_id,
        lineage=lineage, barcode=tuple(barcode), coords=np.asarray(coords, dtype=float),
        counts=counts, gene_id=tuple(gene_id), library_size=library_size,
        metadata=metadata or {},
    )


def reject_target_metadata(metadata: Mapping[str, object]) -> None:
    for key in metadata:
        if any(token in key.lower() for token in FORBIDDEN_INPUT_TOKENS):
            raise R04ContractError(f"target-like metadata is not allowed in model input: {key}")
