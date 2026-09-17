"""Lane B: residual spatial structure after removing axes + distance.

Per section: regress each gene (log1p z-scored) on the 6 axis scores +
signed-hop distance shells to each axis contour (design §7: 回归掉全部
轴＋距离协变量). On the residual field: Moran's I + value-permutation null
(200 draws, valid for autocorrelation per I-021 refined) + LISA tile check
(BH-FDR 999-draw LISA would cost too much per gene; use Moran screen first,
LISA only on screen hits).

Per (gene): median residual-Moran across sections, null p, tile rate,
cross-patient consistency. Upgrade bar: >=2 patients residual-coherent +
effect floor, same family as Tier-2. Exploratory; no claims.

Outputs: infra/r16/laneB_residual_YYYYMMDD.json.
Usage: python3 scripts/r16_laneB_screen.py [--genes N] [--null-draws N]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from r16 import axes as A
from r16 import axis_factory as F
from r16 import census as C

ROOT = Path(__file__).resolve().parent.parent
SHELL_EDGES = (-4, -2, 0, 2, 4, 8)  # distance shells (hops) per axis


def shell_dummies(signed: np.ndarray) -> np.ndarray:
    """One-hot distance shells; returns n×(len-1) design block."""
    cols = []
    for lo, hi in zip(SHELL_EDGES[:-1], SHELL_EDGES[1:]):
        cols.append(((signed >= lo) & (signed < hi)).astype(float))
    D = np.column_stack(cols)
    D = D[:, D.mean(axis=0) > 0.01]  # drop near-empty shells
    return D


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", type=int, default=0,
                    help="0=all cache genes, else top-N by variance (smoke)")
    ap.add_argument("--null-draws", type=int, default=C.NULL_DRAWS)
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "infra/r16/laneB_residual_20260917.json")
    args = ap.parse_args()

    cache = ROOT / "infra/r16/census_cache_hvg10k/sections"
    proxy = json.load(open(ROOT / "infra/r04/marker_proxy_combined.json"))
    union = json.load(open(ROOT / "infra/r16/hvg_union_genes.json"))
    rank = np.load(ROOT / "infra/r16/hvg_final_rank_20260914.npy")
    top10k = {union[i] for i in np.flatnonzero(rank < 10000)}
    axis_defs = {k: [g for g in v["voted_genes"] if g in top10k]
                 for k, v in proxy["classes"].items() if k not in A.RETIRED_AXES}
    axis_defs["Plasma"] = [g for g in A.PLASMA_GENES if g in top10k]

    stems = C.list_stems(cache)
    _, _, genes0, _, _ = C.load_section(cache, stems[0])
    genes0 = list(genes0)
    if args.genes > 0:
        mat, _, _, _, _ = C.load_section(cache, stems[0])
        x = mat.toarray().astype(np.float32)
        np.log1p(x, out=x)
        v = x.var(axis=0)
        genes = [genes0[i] for i in np.argsort(-v)[:args.genes]]
    else:
        genes = genes0
    print(f"genes={len(genes)} sections={len(stems)} null_draws={args.null_draws}",
          flush=True)

    rng = np.random.default_rng(args.seed)
    # per gene: (patient, resid Moran, null p) per section.
    # NOTE (bug caught 2026-09-17): coherence bar must respect the
    # permutation floor 1/(draws+1). With 50 draws the floor is 0.0196,
    # so a p<=0.01 bar is UNREACHABLE by construction. Bar = p<=0.05
    # (reachable: floor 0.0196 < 0.05), counted BY PATIENT (design §7).
    res: dict[str, list] = {g: [] for g in genes}
    n_sec = 0
    for si, stem in enumerate(stems):
        t0 = time.time()
        mat, barcodes, all_genes, coords, meta = C.load_section(cache, stem)
        gidx = {g: i for i, g in enumerate(all_genes)}
        keep = [g for g in genes if g in gidx]
        if not keep:
            continue
        n_sec += 1
        x_z = C.log1p_zscore(mat)
        w = C.spatial_weights(coords, C.KNN_SPATIAL)
        # design: 6 axes + distance shells per axis
        blocks = [np.ones((len(coords), 1))]
        for a, gs in axis_defs.items():
            s, _ = A.axis_scores(x_z, gidx, gs)
            blocks.append(s[:, None])
            m = A.contour_mask(s, 0.7)
            if 0 < m.sum() < len(m):
                signed = A.signed_hop_distance(w, m)
                blocks.append(shell_dummies(signed))
        D = np.column_stack(blocks)
        # OLS residualize all kept genes at once (lstsq, disclosed linear)
        col_idx = np.array([gidx[g] for g in keep])
        X = x_z[:, col_idx]
        beta, _, _, _ = np.linalg.lstsq(D, X, rcond=None)
        R = X - D @ beta
        r2 = 1 - (R ** 2).sum(axis=0) / np.maximum((X ** 2).sum(axis=0), 1e-12)
        # Moran screen on residuals, per gene
        pat = meta["patient_id"]
        for gi, g in enumerate(keep):
            I, p = C.moran_permutation_p(R[:, gi], w, args.null_draws, rng)
            res[g].append((pat, I, p))
        print(f"  [{si+1}/{len(stems)}] {stem[:30]} "
              f"design={D.shape[1]} R2med={np.median(r2):.3f} dt={time.time()-t0:.0f}s",
              flush=True)

    results = []
    for g in genes:
        if not res[g]:
            continue
        Is = [I for _, I, _ in res[g]]
        I_med = float(np.median(Is))
        p_med = float(np.median([p for _, _, p in res[g]]))
        # patient coherence: median p across the patient's sections <= 0.05
        by_pat: dict[str, list] = {}
        for pat, I, p in res[g]:
            by_pat.setdefault(pat, []).append(p)
        coh_pats = [pat for pat, ps in by_pat.items()
                    if float(np.median(ps)) <= 0.05]
        results.append({"gene": g, "n_sections": len(res[g]),
                        "n_patients": len(by_pat),
                        "resid_moran_median": round(I_med, 3),
                        "null_p_median": round(p_med, 4),
                        "n_patients_coherent": len(coh_pats),
                        "coherent_patients": sorted(coh_pats)[:10],
                        "grade": ("EXPLORATORY_REPRODUCED"
                                  if (len(coh_pats) >= 2 and I_med >= 0.1)
                                  else "DESCRIPTIVE_SINGLE_PATIENT")})
    n_rep = sum(1 for r in results if r["grade"] == "EXPLORATORY_REPRODUCED")
    artifact = {"schema": "r16.laneB_residual.v1",
                "status": "EXPLORATORY_COMPLETE_NOT_CLAIM", "seed": args.seed,
                "parameters": {"axes": sorted(axis_defs),
                               "shells": list(zip(SHELL_EDGES[:-1], SHELL_EDGES[1:])),
                               "n_genes": len(genes), "null_draws": args.null_draws,
                               "null": "value_permutation_on_residuals",
                               "regression": "OLS_on_axes_plus_shells"},
                "n_tested": len(results), "n_reproduced": n_rep,
                "results": results}
    args.out.write_text(json.dumps(artifact, ensure_ascii=False))
    print(f"tested={len(results)} reproduced={n_rep} -> {args.out}")
    for r in sorted([x for x in results if x["grade"] == "EXPLORATORY_REPRODUCED"],
                    key=lambda x: -x["resid_moran_median"])[:20]:
        print(f"  {r['gene']:12s} I={r['resid_moran_median']:.3f} "
              f"pmed={r['null_p_median']:.4f} npat_coh={r['n_patients_coherent']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
