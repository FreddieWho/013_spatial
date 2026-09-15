#!/usr/bin/env python3
"""Export R-16 census cache sections as per-section AnnData (.h5ad) for GraphST.

Bridge for the GraphST GPU track (D-121). Reads the HVG-10k cache (raw
integer counts, gene symbols) and writes, per section, <stem>.h5ad with
X=raw counts CSR, obs[section_id, patient_id, barcode],
var_names=gene symbols, obsm['spatial']=coords. Plus manifest.json.
GraphST runs its NATIVE per-section preprocess on these files
(seurat_v3 top-3000 HVG); the HVG-10k universe is the input pool, disclosed
in registry notes. Deterministic; read-only over the cache.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from r16 import census as C

ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=ROOT / "infra/r16/census_cache_hvg10k/sections")
    ap.add_argument("--out", type=Path, default=ROOT / "infra/r16/graphst_bridge")
    ap.add_argument("--sections", nargs="*", default=None,
                    help="subset of stems (default: all)")
    args = ap.parse_args()

    import anndata as ad

    stems = C.list_stems(args.cache)
    if args.sections:
        want = set(args.sections)
        stems = [s for s in stems if s in want]
    if not stems:
        raise SystemExit(f"empty selection in cache: {args.cache}")
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = []
    for stem in stems:
        mat, barcodes, genes, coords, meta = C.load_section(args.cache, stem)
        # integrity: raw integer counts required downstream
        if not np.all(np.asarray(mat.data) == np.floor(np.asarray(mat.data))):
            raise SystemExit(f"non-integer counts in {stem}; abort")
        if len(set(genes)) != len(genes):
            raise SystemExit(f"duplicate gene symbols in {stem}; abort")
        a = ad.AnnData(X=mat.tocsr())
        a.obs["section_id"] = meta["section_id"]
        a.obs["patient_id"] = meta["patient_id"]
        a.obs["barcode"] = list(barcodes)
        a.obs_names = list(barcodes)
        a.var_names = list(genes)
        a.obsm["spatial"] = np.asarray(coords, dtype=np.float32)
        a.write(args.out / f"{stem}.h5ad")
        manifest.append({"stem": stem, "section_id": meta["section_id"],
                         "patient_id": meta["patient_id"],
                         "n_spots": len(barcodes), "n_genes": len(genes)})
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"exported {len(manifest)} sections -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
