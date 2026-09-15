#!/usr/bin/env python3
"""GraphST GPU run driver, arm G (D-121).

Per-section training (NOT merged: GraphST builds a DENSE O(n^2) distance
matrix in construct_interaction; 113k pooled spots ~= 100GB, infeasible on
24/48GB cards; disclosed design deviation). Cross-section recurrence is
tested downstream by the standard registry machinery (matching/purity/
Moran/split-half), never forced by joint training.

Native GraphST behavior otherwise: per-section seurat_v3 top-3000 HVG
preprocess on the HVG-10k input pool (disclosed approximation), default
loss weights, full epochs (NO reduction per user requirement).
Deviations from upstream defaults, all disclosed: random_seed=20260914
(upstream 41); Leiden clustering instead of mclust (no R bridge on rental;
matches Tier-2 registry convention); resolutions 0.25/0.5/1.0 + min 20.

Outputs per section (incremental, preemption-safe): labels_armG.csv rows
(section, orig_barcode, GRAPHST_025/05/10); params.json; run.log.
Scored downstream by scripts/r16_score_r_labels.py (arm G).
Exploratory grade; no claims.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from r16 import census as C

ROOT = Path(__file__).resolve().parent.parent
RES_COLS = ["GRAPHST_025", "GRAPHST_05", "GRAPHST_10"]


def resolve_device(want: str) -> str:
    import torch

    if want == "cpu":
        return "cpu"
    if want == "cuda":
        if not torch.cuda.is_available():
            raise SystemExit("device=cuda requested but unavailable")
        return "cuda:0"
    return "cuda:0" if torch.cuda.is_available() else "cpu"


def leiden_on_embedding(emb: np.ndarray):
    """Leiden at registry resolutions on embedding kNN graph + min-size filter."""
    g = C.knn_igraph(np.asarray(emb, dtype=np.float64), C.KNN_EXPR)
    out = {}
    for res, col in zip(C.RESOLUTIONS, RES_COLS):
        out[col] = C.filter_min_size(C.leiden_labels(g, res), C.MIN_CLUSTER_SIZE)
    return out


def run_one(h5ad_path: Path, device: str, epochs: int):
    """Train GraphST on one section file.
    Returns (adata_with_emb, emb, seconds, hvg_method).
    Fallback (disclosed): scanpy seurat_v3 HVG uses loess, which goes
    singular on degenerate sections (hit 2026-09-16, full run aborted at
    section 19/47). On ValueError/singularity we replicate GraphST's own
    preprocess steps with flavor='seurat' (dispersion, no loess) and retry
    once; the section is tagged hvg-seurat-fallback in logs/params.
    """
    import anndata as ad
    # NOTE: GraphST/__init__.py does NOT re-export the class, so
    # `from GraphST import GraphST` binds the SUBMODULE (TypeError).
    # Correct: import the class from the submodule (verified 2026-09-16).
    from GraphST.GraphST import GraphST

    t0 = time.time()
    adata = ad.read_h5ad(h5ad_path)
    hvg_method = "seurat_v3"
    try:
        model = GraphST(
            adata, device=device, epochs=epochs,
            dim_input=3000, dim_output=64, random_seed=C.SEED,
        )
    except ValueError as e:
        if "singularit" not in str(e).lower():
            raise
        import scanpy as sc
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.scale(adata, zero_center=False, max_value=10)
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=3000)
        hvg_method = "seurat-fallback"
        model = GraphST(
            adata, device=device, epochs=epochs,
            dim_input=3000, dim_output=64, random_seed=C.SEED,
        )
    adata = model.train()
    emb = np.asarray(adata.obsm["emb"])
    return adata, emb, time.time() - t0, hvg_method


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, required=True,
                    help="graphst_bridge dir with <stem>.h5ad + manifest.json")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--epochs", type=int, default=600)
    ap.add_argument("--device", default="auto", choices=["auto", "cuda", "cpu"])
    ap.add_argument("--dry-run", action="store_true",
                    help="first section only, 3 epochs, full path incl. clustering")
    ap.add_argument("--sections", nargs="*", default=None)
    args = ap.parse_args()

    import torch

    dev = resolve_device(args.device)
    manifest = json.loads((args.data / "manifest.json").read_text())
    stems = [r["stem"] for r in manifest]
    if args.sections:
        want = set(args.sections)
        stems = [s for s in stems if s in want]
    if args.dry_run:
        stems = stems[:1]
        args.epochs = 3
    if not stems:
        raise SystemExit("empty section selection")
    args.out.mkdir(parents=True, exist_ok=True)
    gver_file = args.out / "graphst_version.txt"
    graphst_version = gver_file.read_text().strip() if gver_file.exists() else "unpinned"
    logf = open(args.out / "run.log", "a")
    params = {
        "schema": "r16.graphst_run.v1", "seed": C.SEED,
        "graphst_commit": graphst_version,
        "torch": torch.__version__, "cuda_available": torch.cuda.is_available(),
        "cuda_device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        "epochs": args.epochs, "resolutions": list(C.RESOLUTIONS),
        "min_cluster_size": C.MIN_CLUSTER_SIZE, "knn_expr": C.KNN_EXPR,
        "per_section_training": True,
        "merged_training_rejected": "dense O(n^2) dist matrix ~100GB at 113k spots",
        "hvg": "native seurat_v3 top-3000 per section on HVG-10k pool",
        "clustering": "Leiden (not upstream mclust): no R bridge, Tier-2 convention",
        "upstream_seed_default_overridden": 41,
    }
    (args.out / "params.json").write_text(json.dumps(params, indent=2))
    labels_path = args.out / "labels_armG.csv"
    header_written = labels_path.exists()
    # resume: skip sections already fully written (incremental safe rerun)
    done_stems = set()
    if header_written:
        import pandas as _pd
        done_stems = set(_pd.read_csv(labels_path, usecols=["section"])["section"].tolist())
        print(f"resume: {len(done_stems)} sections already done, skipping", flush=True)
    fallbacks = []
    total_t = 0.0
    for i, stem in enumerate(stems):
        if stem in done_stems:
            print(f"[{i+1}/{len(stems)}] {stem} SKIPPED (already done)", flush=True)
            continue
        adata, emb, secs, hvg_method = run_one(args.data / f"{stem}.h5ad", dev, args.epochs)
        if hvg_method != "seurat_v3":
            fallbacks.append(stem)
        total_t += secs
        labs = leiden_on_embedding(emb)
        df = pd.DataFrame({
            "section": stem,
            "orig_barcode": list(adata.obs["barcode"]),
            RES_COLS[0]: labs[RES_COLS[0]],
            RES_COLS[1]: labs[RES_COLS[1]],
            RES_COLS[2]: labs[RES_COLS[2]],
        })
        df.to_csv(labels_path, mode="a", header=not header_written, index=False)
        header_written = True
        n_cl = {c: int((df[c] >= 0).sum()) for c in RES_COLS}
        msg = (f"[{i+1}/{len(stems)}] {stem} n={len(df)} "
               f"clusters={ {c: int(df[c].nunique()) for c in RES_COLS} } "
               f"hvg={hvg_method} {secs:.0f}s elapsed={total_t:.0f}s")
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()
    logf.close()
    print(f"DONE sections={len(stems)} hvg_fallbacks={fallbacks} total_s={total_t:.0f}")
    (args.out / "hvg_fallbacks.json").write_text(json.dumps(fallbacks, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
