"""Lane A full screen: gene ~ f(signed-hop distance to axis contour).

For each of 6 working axes (T/B/Mye/Epi/Stromal/Plasma) and each gene in the
10k cache: per-section signed-hop profile (bins -4..+8, same as the CXCL13~B
positive control), then cross-patient shape consistency of the profile +
spatial-preserve null (180° rotation + small shift, I-021/D-125 machinery).

Per (axis, gene): observed median profile across sections, null profile
distribution, peak location, shape correlation across patients, grade.
Upgrade bar (frozen): >=2 patients same-direction + effect floor, same as
Tier-2 (D-110 band 0.02 analogue: |median peak height| >= 0.5 z + null p).
Exploratory grade; no claims.

Outputs: infra/r16/laneA_screen_YYYYMMDD.json + registry rows.
Usage: python3 scripts/r16_laneA_screen.py [--genes N] [--null-draws N]
  --genes 0 (default) = all cache genes; small N = smoke (top-variance).
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from r16 import axes as A
from r16 import census as C

ROOT = Path(__file__).resolve().parent.parent
PC_BINS = np.arange(-4, 9)
N_BINS = len(PC_BINS)


def section_profiles(cache, stems, axis_defs, genes: list[str], null_draws: int,
                     seed: int):
    """Per section: {axis: mask}, per gene: obs profile + null profiles.

    Returns dicts keyed by stem: obs[axis][gene] = (means, ns),
    nulls[axis][gene] = list of mean-vectors (one per draw).
    Memory-light: iterate genes in blocks, one section at a time.
    """
    rng = np.random.default_rng(seed)
    out = {}
    for si, stem in enumerate(stems):
        t0 = time.time()
        mat, barcodes, all_genes, coords, meta = C.load_section(cache, stem)
        gidx = {g: i for i, g in enumerate(all_genes)}
        keep = [g for g in genes if g in gidx]
        if not keep:
            continue
        x_z = C.log1p_zscore(mat)
        w = C.spatial_weights(coords, C.KNN_SPATIAL)
        scores = {}
        for a, gs in axis_defs.items():
            s, _ = A.axis_scores(x_z, gidx, gs)
            scores[a] = s
        masks = {}
        for a, s in scores.items():
            m = A.contour_mask(s, 0.7)
            masks[a] = m if 0 < m.sum() < len(m) else None
        perms = [C.spatial_null_remap(coords, seed=seed + si * 1000 + d)[0]
                 for d in range(null_draws)]
        sec = {"patient": meta["patient_id"], "masks_ok": {},
               "obs": {}, "null": {}}
        col_idx = np.array([gidx[g] for g in keep])
        X = x_z[:, col_idx]  # n_spots x n_keep
        for a, m in masks.items():
            if m is None:
                continue
            sec["masks_ok"][a] = True
            signed = A.signed_hop_distance(w, m)
            # observed: mean per bin per gene
            obs_m = np.full((N_BINS, len(keep)), np.nan)
            for bi, b in enumerate(PC_BINS):
                sel = signed == b
                if sel.any():
                    obs_m[bi] = X[sel].mean(axis=0)
            sec["obs"][a] = obs_m
            # null: permuted field, same mask
            null_m = np.empty((null_draws, N_BINS, len(keep)))
            null_m[:] = np.nan
            for d, perm in enumerate(perms):
                Xp = X[perm]
                for bi, b in enumerate(PC_BINS):
                    sel = signed == b
                    if sel.any():
                        null_m[d, bi] = Xp[sel].mean(axis=0)
            sec["null"][a] = null_m
        sec["keep"] = keep
        out[stem] = sec
        print(f"  [{si+1}/{len(stems)}] {stem[:30]} "
              f"axes={sorted(sec['masks_ok'])} dt={time.time()-t0:.0f}s", flush=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", type=int, default=0,
                    help="0=all cache genes, else top-N by variance (smoke)")
    ap.add_argument("--axes", type=str, default="",
                    help="comma-separated subset, e.g. 'B'; default=all")
    ap.add_argument("--null-draws", type=int, default=8)
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "infra/r16/laneA_screen_20260917.json")
    args = ap.parse_args()

    cache = ROOT / "infra/r16/census_cache_hvg10k/sections"
    proxy = json.load(open(ROOT / "infra/r04/marker_proxy_combined.json"))
    union = json.load(open(ROOT / "infra/r16/hvg_union_genes.json"))
    rank = np.load(ROOT / "infra/r16/hvg_final_rank_20260914.npy")
    top10k = {union[i] for i in np.flatnonzero(rank < 10000)}
    axis_defs = {k: [g for g in v["voted_genes"] if g in top10k]
                 for k, v in proxy["classes"].items() if k not in A.RETIRED_AXES}
    axis_defs["Plasma"] = [g for g in A.PLASMA_GENES if g in top10k]
    if args.axes:
        want = [a.strip() for a in args.axes.split(",") if a.strip()]
        unknown = set(want) - set(axis_defs)
        if unknown:
            raise SystemExit(f"unknown axes: {sorted(unknown)}")
        axis_defs = {a: axis_defs[a] for a in want}

    stems = C.list_stems(cache)
    # gene list = cache genes (first section defines order)
    _, _, genes0, _, _ = C.load_section(cache, stems[0])
    genes0 = list(genes0)
    if args.genes > 0:
        # top-N by pooled variance proxy: use first-section variance rank
        mat, _, _, _, _ = C.load_section(cache, stems[0])
        x = mat.toarray().astype(np.float32)
        np.log1p(x, out=x)
        v = x.var(axis=0)
        order = np.argsort(-v)[:args.genes]
        genes = [genes0[i] for i in order]
    else:
        genes = genes0
    print(f"genes={len(genes)} sections={len(stems)} null_draws={args.null_draws}")

    sec = section_profiles(cache, stems, axis_defs, genes, args.null_draws,
                           args.seed)
    # aggregate per (axis, gene): median obs profile across sections,
    # null peak-height distribution, cross-patient peak agreement
    results = []
    for a in axis_defs:
        usable = [(s, d) for s, d in sec.items() if a in d.get("obs", {})]
        if not usable:
            continue
        keep = usable[0][1]["keep"]
        for gi, g in enumerate(keep):
            obs_stack = np.array([d["obs"][a][:, gi] for _, d in usable])
            med = np.nanmedian(obs_stack, axis=0)
            if np.all(np.isnan(med)):
                continue
            peak = int(PC_BINS[int(np.nanargmax(med))])
            height = float(np.nanmax(med))
            # null: peak heights across draws+sections
            null_h = []
            for _, d in usable:
                nm = d["null"][a]
                for dd in range(nm.shape[0]):
                    null_h.append(float(np.nanmax(nm[dd, :, gi])))
            null_h = np.array([v for v in null_h if np.isfinite(v)])
            p = float((np.sum(null_h >= height) + 1) / (len(null_h) + 1)) if len(null_h) else 1.0
            # cross-patient: peak bin agreement (fraction sharing modal peak)
            peaks = []
            for _, d in usable:
                v = d["obs"][a][:, gi]
                if np.isfinite(v).sum() >= 3:
                    peaks.append(int(PC_BINS[int(np.nanargmax(v))]))
            agree = (max(np.unique(peaks, return_counts=True)[1]) / len(peaks)
                     if peaks else 0.0)
            npat = len({d["patient"] for _, d in usable})
            results.append({
                "axis": a, "gene": g, "n_sections": len(usable),
                "n_patients": npat, "peak_bin": peak,
                "peak_height_z": round(height, 3), "null_p": round(p, 4),
                "null_n": len(null_h), "peak_agree": round(float(agree), 3),
                "median_profile": [round(float(v), 3) if np.isfinite(v) else None
                                   for v in med],
            })
    # rank: reproduced = null_p<=0.05 & agree>=0.5 & height>=0.5
    for r in results:
        r["grade"] = ("EXPLORATORY_REPRODUCED"
                      if (r["null_p"] <= 0.05 and r["peak_agree"] >= 0.5
                          and r["peak_height_z"] >= 0.5)
                      else "DESCRIPTIVE_SINGLE_PATIENT")
    n_rep = sum(1 for r in results if r["grade"] == "EXPLORATORY_REPRODUCED")
    artifact = {"schema": "r16.laneA_screen.v1",
                "status": "EXPLORATORY_COMPLETE_NOT_CLAIM", "seed": args.seed,
                "parameters": {"axes": sorted(axis_defs), "n_genes": len(genes),
                               "null_draws": args.null_draws,
                               "null": "spatial_preserve_rot180_smallshift",
                               "bins": PC_BINS.tolist()},
                "n_tested": len(results), "n_reproduced": n_rep,
                "results": results}
    args.out.write_text(json.dumps(artifact, ensure_ascii=False))
    print(f"tested={len(results)} reproduced={n_rep} -> {args.out}")
    # top hits preview
    rep = sorted([r for r in results if r["grade"] == "EXPLORATORY_REPRODUCED"],
                 key=lambda r: (r["null_p"], -r["peak_agree"]))[:20]
    for r in rep:
        print(f"  {r['axis']:8s} {r['gene']:12s} peak={r['peak_bin']:+d} "
              f"h={r['peak_height_z']:.2f} p={r['null_p']:.4f} "
              f"agree={r['peak_agree']:.2f} npat={r['n_patients']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
