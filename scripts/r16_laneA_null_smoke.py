"""Smoke-test the spatial-preserve null on real cache geometry (I-021, Lane A).

1. Every section: remap acceptance diagnostics (frac within 0.5 hop).
2. Positive control: CXCL13~B signed-hop profile shape must survive the
   OBSERVED curve; null curves (field rotated/translated) must flatten it.
   Pass bar: observed peak <= +1 hop AND null peak distribution centered
   away (median null peak |.| > 1 or null peak height < 50% of observed).
   Exploratory grade; no claims.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from r16 import axes as A
from r16 import census as C

ROOT = Path(__file__).resolve().parent.parent
PC_BINS = np.arange(-4, 9)


def main() -> int:
    cache = ROOT / "infra/r16/census_cache_hvg10k/sections"
    proxy = json.load(open(ROOT / "infra/r04/marker_proxy_combined.json"))
    union = json.load(open(ROOT / "infra/r16/hvg_union_genes.json"))
    rank = np.load(ROOT / "infra/r16/hvg_final_rank_20260914.npy")
    top = {union[i] for i in np.flatnonzero(rank < 10000)}
    b_genes = [g for g in proxy["classes"]["B"]["voted_genes"] if g in top]

    stems = C.list_stems(cache)
    acc_ok, acc_all, null_peaks, obs_peaks = 0, 0, [], []
    rng_seed = 0
    for stem in stems:
        mat, barcodes, genes, coords, meta = C.load_section(cache, stem)
        gidx = {g: i for i, g in enumerate(genes)}
        if "CXCL13" not in gidx:
            continue
        x_z = C.log1p_zscore(mat)
        w = C.spatial_weights(coords, C.KNN_SPATIAL)
        bscore, _ = A.axis_scores(x_z, gidx, b_genes)
        bmask = A.contour_mask(bscore, 0.7)
        if not (0 < bmask.sum() < len(bmask)):
            continue
        signed = A.signed_hop_distance(w, bmask)
        obs = A.binned_profile(x_z[:, gidx["CXCL13"]], signed, PC_BINS)
        ov = np.array([v for v, n in obs])
        if np.all(np.isnan(ov)):
            continue
        obs_peaks.append(int(PC_BINS[int(np.nanargmax(ov))]))
        # 8 null draws: rotate/translate the FIELD, keep mask fixed
        for d in range(8):
            perm, info = C.spatial_null_remap(coords, seed=C.SEED + rng_seed)
            rng_seed += 1
            acc_all += 1
            acc_ok += int(info["accepted"])
            moved_vals = x_z[perm][:, gidx["CXCL13"]]
            null_prof = A.binned_profile(moved_vals, signed, PC_BINS)
            nv = np.array([v for v, n in null_prof])
            if not np.all(np.isnan(nv)):
                null_peaks.append(int(PC_BINS[int(np.nanargmax(nv))]))
    obs = np.array(obs_peaks)
    nul = np.array(null_peaks)
    print(f"sections_with_PC={len(obs)} null_draws={len(nul)}")
    print(f"remap_accept_rate={acc_ok}/{acc_all}={acc_ok/max(acc_all,1):.3f}")
    print(f"obs peaks<=+1: {(obs <= 1).sum()}/{len(obs)} "
          f"dist={dict(zip(*np.unique(obs, return_counts=True)))}")
    print(f"null peaks median={np.median(nul):.1f} "
          f"dist={dict(zip(*np.unique(nul, return_counts=True)))}")
    print(f"null |peak|>1 frac={(np.abs(nul) > 1).mean():.3f}")
    ok = (obs <= 1).mean() >= 0.6 and (np.abs(nul) > 1).mean() >= 0.5
    print("SMOKE", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
