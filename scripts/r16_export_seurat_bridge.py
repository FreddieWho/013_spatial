#!/usr/bin/env python3
"""Export R-16 census cache sections as 10x-style Matrix Market triples + coords.

Bridge for the R integration track (Seurat anchors / Harmony / PRECAST).
Reads the HVG-10k cache (raw integer counts, gene symbols) and writes, per
section, matrix.mtx.gz (genes x spots), features.tsv.gz (symbol x2),
barcodes.tsv.gz, coords.csv (barcode,x,y). Plus manifest.json.
Deterministic; read-only over the cache.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmwrite

from r16 import census as C

ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=ROOT / "infra/r16/census_cache_hvg10k/sections")
    ap.add_argument("--out", type=Path, default=ROOT / "infra/r16/seurat_bridge")
    args = ap.parse_args()

    from r16 import census as C

    stems = C.list_stems(args.cache)
    if not stems:
        raise SystemExit(f"empty cache: {args.cache}")
    manifest = []
    for stem in stems:
        mat, barcodes, genes, coords, meta = C.load_section(args.cache, stem)
        # integrity: raw integer counts required downstream (SCTransform)
        if not np.all(np.asarray(mat.data) == np.floor(np.asarray(mat.data))):
            raise SystemExit(f"non-integer counts in {stem}; abort (SCT contract)")
        d = args.out / stem
        d.mkdir(parents=True, exist_ok=True)
        mmwrite(str(d / "matrix.mtx"), mat.T.tocoo())
        with gzip.open(d / "matrix.mtx.gz", "wb") as f:
            f.write((d / "matrix.mtx").read_bytes())
        (d / "matrix.mtx").unlink()
        with gzip.open(d / "features.tsv.gz", "wt") as f:
            for g in genes:
                f.write(f"{g}\t{g}\tGene Expression\n")
        with gzip.open(d / "barcodes.tsv.gz", "wt") as f:
            for b in barcodes:
                f.write(f"{b}\n")
        with open(d / "coords.csv", "w") as f:
            f.write("barcode,x,y\n")
            for b, (x, y) in zip(barcodes, coords):
                f.write(f"{b},{x},{y}\n")
        manifest.append({"stem": stem, "section_id": meta["section_id"],
                         "patient_id": meta["patient_id"],
                         "n_spots": len(barcodes), "n_genes": len(genes)})
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"exported {len(manifest)} sections -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
