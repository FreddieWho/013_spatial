#!/usr/bin/env python3
"""Materialize the R-16 census cache (HVG top-10000, D-114).

Reads each training-row source h5ad, keeps union genes with pilot final rank
< top_n, writes CSR counts + coords + barcodes per section into
infra/r16/census_cache_hvg10k/sections/. Deterministic; R-04 assets untouched.
Duplicate symbols in a section's var keep the first occurrence (disclosed).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parent.parent


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--training-manifest", type=Path,
                    default=ROOT / "infra/r04/role_manifests/training_manifest.json")
    ap.add_argument("--rank-file", type=Path,
                    default=ROOT / "infra/r16/hvg_final_rank_20260914.npy")
    ap.add_argument("--union-genes", type=Path, default=Path("/tmp/det_union_genes.json"))
    ap.add_argument("--top-n", type=int, default=10000)
    ap.add_argument("--out", type=Path, default=ROOT / "infra/r16/census_cache_hvg10k/sections")
    args = ap.parse_args()

    final_rank = np.load(args.rank_file)
    union = json.load(open(args.union_genes))
    union_index = {s: i for i, s in enumerate(union)}
    keep = [int(i) for i in np.flatnonzero(final_rank < args.top_n)]
    keep_set = set(keep)
    genes_out = [union[i] for i in keep]
    assert len(genes_out) == args.top_n, f"expected {args.top_n}, got {len(genes_out)}"
    out_pos_of_union = {u: k for k, u in enumerate(keep)}

    args.out.mkdir(parents=True, exist_ok=True)
    man = json.load(open(args.training_manifest))
    rows_out = []
    n_dup = 0
    for r in man["rows"]:
        with h5py.File(r["matrix_locator"]) as f:
            X = f["X"]
            if isinstance(X, h5py.Group):
                mat = sparse.csr_matrix(
                    (X["data"][:], X["indices"][:], X["indptr"][:]),
                    shape=tuple(X.attrs["shape"]))
            else:
                mat = sparse.csr_matrix(X[:].astype(np.float32))
            syms = [s.decode() if isinstance(s, bytes) else s for s in f["var"]["_index"][:]]
            barcodes = [b.decode() if isinstance(b, bytes) else b for b in f["obs"]["_index"][:]]
            coords = f["obsm"]["spatial"][:].astype(np.float32)
        src_cols, out_pos, seen = [], [], set()
        for j, s in enumerate(syms):
            u = union_index.get(s)
            if u is None or u not in keep_set or u in seen:
                if u in seen:
                    n_dup += 1
                continue
            seen.add(u)
            src_cols.append(j)
            out_pos.append(out_pos_of_union[u])
        # sparse column-selection matrix: source_col -> out_pos (missing genes
        # become all-zero columns automatically, order always matches genes_out)
        sel = sparse.csr_matrix(
            (np.ones(len(src_cols)), (src_cols, out_pos)),
            shape=(mat.shape[1], args.top_n))
        sub = (mat @ sel).tocsr()
        missing = args.top_n - len(src_cols)
        assert sub.shape == (len(barcodes), args.top_n), sub.shape
        assert len(coords) == len(barcodes)
        stem = f"{len(rows_out):04d}_" + r["section_id"].replace("::", "_").replace("/", "_")
        np.savez_compressed(args.out / f"{stem}.npz",
                            data=sub.data, indices=sub.indices, indptr=sub.indptr,
                            shape=np.array(sub.shape), coords=coords)
        meta = {
            "schema": "r16.census_cache_section.v1",
            "section_id": r["section_id"],
            "patient_id": r["patient_id"],
            "lineage": r["lineage"],
            "barcodes": barcodes,
            "genes": genes_out,
            "n_genes_missing_in_source": int(missing),
            "source_matrix_locator": r["matrix_locator"],
        }
        (args.out / f"{stem}.json").write_text(json.dumps(meta))
        rows_out.append({"stem": stem, "section_id": r["section_id"],
                         "patient_id": r["patient_id"], "n_spots": len(barcodes),
                         "n_genes_missing_in_source": int(missing)})
        print(f"{stem}: spots={len(barcodes)} missing_genes={missing}")

    manifest = {
        "schema": "r16.census_cache_manifest.v1",
        "gene_universe": "integration_hvg_top10000_D114",
        "n_genes": args.top_n,
        "n_sections": len(rows_out),
        "duplicate_symbols_skipped": n_dup,
        "rank_file_sha256": sha256(args.rank_file),
        "rows": rows_out,
    }
    (args.out.parent / "cache_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"materialized {len(rows_out)} sections, {args.top_n} genes, dup_skipped={n_dup}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
