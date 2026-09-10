#!/usr/bin/env python3
"""Build a downsampled CRC scRNA reference from GSE236581 (original annotations only).

Reads the author-provided count matrix + metadata, maps panel ENSG IDs to
symbols via the frozen gencode v49 table, stratified-caps cells per
(MajorCellType x Patient), and writes a compact reference bundle:
counts CSR (panel genes x sampled cells), cell labels (Major + Sub, verbatim
author strings), patient/tissue vectors, and full provenance (hashes, seed,
cap rule, excluded inputs).

No re-annotation, no harmonization, no clustering. Gene symbols that map
ambiguously (one symbol <- multiple ENSG) are resolved by preferring the
panel ENSG's own symbol; unmapped panel genes are recorded and dropped.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_symbol_map(mapping_csv: Path, panel_ids: list[str]) -> tuple[dict[str, str], list[str]]:
    """Map panel ENSG IDs -> symbols. Returns (ensg->symbol, ordered panel symbols)."""
    table: dict[str, list[str]] = {}
    with open(mapping_csv, newline="") as f:
        for row in csv.reader(f):
            if len(row) >= 2:
                table.setdefault(row[0], []).append(row[1])
    mapping: dict[str, str] = {}
    dropped: list[str] = []
    for g in panel_ids:
        base = g.replace("DEPRECATED_", "").split(".")[0]
        syms = table.get(base, [])
        if not syms:
            dropped.append(g)
            continue
        mapping[g] = syms[0]  # panel ENSG's own symbol; duplicates resolved identically
    ordered = [mapping[g] for g in panel_ids if g not in dropped]
    return mapping, dropped


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mtx", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--barcodes", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--symbol-map", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260910)
    parser.add_argument("--cap-per-class-patient", type=int, default=1000)
    parser.add_argument("--tissues", default="Tumor",
                        help="comma-separated Tissue values to include")
    args = parser.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    panel_ids = [l.strip() for l in Path(args.panel).read_text().splitlines() if l.strip()]
    ensg_to_sym, dropped = load_symbol_map(args.symbol_map, panel_ids)
    panel_syms = [ensg_to_sym[g] for g in panel_ids if g not in dropped]

    with gzip.open(args.features, "rt") as f:
        ref_genes = [l.split("\t")[0] for l in f]
    ref_index = {g: i for i, g in enumerate(ref_genes)}
    keep_rows = sorted({ref_index[s] for s in panel_syms if s in ref_index})
    missing_syms = sorted(set(panel_syms) - set(ref_index))
    print(f"panel={len(panel_ids)} mapped={len(panel_syms)} in_ref={len(keep_rows)} "
          f"missing_syms={len(missing_syms)} dropped_ensg={len(dropped)}", flush=True)

    with gzip.open(args.barcodes, "rt") as f:
        barcodes = [l.strip() for l in f]
    meta_major, meta_sub, meta_pat, meta_tis = {}, {}, {}, {}
    with gzip.open(args.metadata, "rt") as f:
        rdr = csv.reader(f, delimiter=" ", skipinitialspace=True)
        next(rdr)
        for row in rdr:
            bc = row[0].strip('"')
            meta_major[bc] = row[8].strip('"')
            meta_sub[bc] = row[9].strip('"')
            meta_pat[bc] = row[5].strip('"')
            meta_tis[bc] = row[7].strip('"')  # row[6] is Treatment/stage, NOT tissue
    assert len(meta_major) == len(barcodes), "metadata/barcode count mismatch"

    tissues = {t.strip() for t in args.tissues.split(",")}
    rng = np.random.default_rng(args.seed)
    groups: dict[tuple[str, str], list[int]] = {}
    for i, bc in enumerate(barcodes):
        if meta_tis[bc] not in tissues:
            continue
        groups.setdefault((meta_major[bc], meta_pat[bc]), []).append(i)
    sampled: list[int] = []
    for key in sorted(groups):
        idx = groups[key]
        if len(idx) > args.cap_per_class_patient:
            idx = sorted(rng.choice(idx, args.cap_per_class_patient, replace=False).tolist())
        sampled.extend(idx)
    sampled = sorted(sampled)
    keep_cols = set(sampled)
    print(f"groups={len(groups)} sampled_cells={len(sampled)}", flush=True)

    # Stream the mtx once, keeping panel rows x sampled cols.
    row_set = set(keep_rows)
    row_remap = {r: k for k, r in enumerate(keep_rows)}
    col_remap = {c: k for k, c in enumerate(sampled)}
    data, indices, indptr = [], [], [0]
    n_rows_file = n_cols_file = 0
    with gzip.open(args.mtx, "rt") as f:
        for line in f:
            if line.startswith("%"):
                continue
            parts = line.split()
            if n_rows_file == 0 and len(parts) == 3:
                n_rows_file, n_cols_file, _ = map(int, parts)
                # mtx is genes(rows) x cells(cols); build CSR row by row is
                # awkward streaming column-major; collect COO triplets instead.
                coo_r, coo_c, coo_v = [], [], []
                continue
            r, c = int(parts[0]) - 1, int(parts[1]) - 1
            if r in row_set and c in keep_cols:
                coo_r.append(row_remap[r])
                coo_c.append(col_remap[c])
                coo_v.append(float(parts[2]))
    mat = sparse.csr_matrix((coo_v, (coo_r, coo_c)),
                            shape=(len(keep_rows), len(sampled)), dtype=np.float32)
    print(f"kept nnz={mat.nnz} shape={mat.shape}", flush=True)

    sym_to_ensg: dict[str, str] = {}
    for g in panel_ids:
        if g not in dropped:
            sym_to_ensg.setdefault(ensg_to_sym[g], g)
    kept_symbols = [ref_genes[r] for r in keep_rows]
    np.savez_compressed(
        out_dir / "reference_counts.npz",
        counts_data=mat.data, counts_indices=mat.indices,
        counts_indptr=mat.indptr, counts_shape=np.asarray(mat.shape),
        panel_symbols=np.asarray(kept_symbols),
        panel_ensg=np.asarray([sym_to_ensg[s] for s in kept_symbols]),
        barcodes=np.asarray([barcodes[c] for c in sampled]),
        major=np.asarray([meta_major[barcodes[c]] for c in sampled]),
        sub=np.asarray([meta_sub[barcodes[c]] for c in sampled]),
        patient=np.asarray([meta_pat[barcodes[c]] for c in sampled]),
        tissue=np.asarray([meta_tis[barcodes[c]] for c in sampled]),
    )
    provenance = {
        "schema": "r04.composition_reference.v1",
        "source": "GSE236581 CRC-ICB (author counts + author metadata, unmodified)",
        "annotations": "verbatim author MajorCellType/SubCellType; no re-annotation, no harmonization",
        "inputs": {k: {"path": str(getattr(args, k)), "sha256": sha256_file(getattr(args, k))}
                   for k in ("mtx", "features", "barcodes", "metadata")},
        "panel": str(args.panel),
        "symbol_map": str(args.symbol_map),
        "seed": args.seed,
        "cap_per_class_patient": args.cap_per_class_patient,
        "tissues": sorted(tissues),
        "dropped_ensg_unmapped": dropped,
        "missing_symbols_not_in_ref": missing_syms,
        "n_cells_sampled": len(sampled),
        "n_genes_kept": len(keep_rows),
        "matrix_dims_file": [n_rows_file, n_cols_file],
    }
    (out_dir / "provenance.json").write_text(json.dumps(provenance, indent=1))
    print("REFERENCE_OK", out_dir, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
