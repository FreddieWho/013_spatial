#!/usr/bin/env python3
"""R-16 Tier-1 axis ruler factory (design v2.0, D-115).

Per section: build composition axis scores (T/B/Mye/ILC/Epi/Stromal/Plasma),
test spatial coherence (Moran's I, value-permutation null, 200 draws), measure
contour-mask stability across disclosed quantiles, and run the end-to-end
positive control (CXCL13 profile vs B-cell contour, signed hop distance).

Outputs axis_rulers JSON artifact + tier-1 registry rows JSON for
scripts/r16_build_registry.py. Exploratory grade; no claims.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from r16 import axes as A
from r16 import census as C

ROOT = Path(__file__).resolve().parent.parent
PC_BINS = np.arange(-4, 9)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=ROOT / "infra/r16/census_cache_hvg10k/sections")
    ap.add_argument("--marker-proxy", type=Path, default=ROOT / "infra/r04/marker_proxy_combined.json")
    ap.add_argument("--union-genes", type=Path, default=ROOT / "infra/r16/hvg_union_genes.json")
    ap.add_argument("--rank-file", type=Path, default=ROOT / "infra/r16/hvg_final_rank_20260914.npy")
    ap.add_argument("--top-n", type=int, default=10000)
    ap.add_argument("--null-draws", type=int, default=C.NULL_DRAWS)
    ap.add_argument("--out", type=Path, default=ROOT / "infra/r16/axis_rulers_20260914.json")
    ap.add_argument("--out-rows", type=Path, default=ROOT / "infra/r16/registry_rows_tier1_20260914.json")
    args = ap.parse_args()

    proxy = json.load(open(args.marker_proxy))
    union = json.load(open(args.union_genes))
    rank = np.load(args.rank_file)
    top = {union[i] for i in np.flatnonzero(rank < args.top_n)}
    axis_defs = {k: [g for g in v["voted_genes"] if g in top]
                 for k, v in proxy["classes"].items()}
    axis_defs["Plasma"] = [g for g in A.PLASMA_GENES if g in top]
    rng = np.random.default_rng(C.SEED)
    stems = C.list_stems(args.cache)

    coh = {a: [] for a in axis_defs}
    stab = {a: [] for a in axis_defs}
    n_genes_used = {a: [] for a in axis_defs}
    pat_profiles: dict[str, np.ndarray] = {}
    n_sections = 0
    for stem in stems:
        mat, barcodes, genes, coords, meta = C.load_section(args.cache, stem)
        n_sections += 1
        gidx = {g: i for i, g in enumerate(genes)}
        x_z = C.log1p_zscore(mat)
        w = C.spatial_weights(coords, C.KNN_SPATIAL)
        scores = {}
        for a, gs in axis_defs.items():
            s, n_used = A.axis_scores(x_z, gidx, gs)
            scores[a] = s
            n_genes_used[a].append(n_used)
            I, p = C.moran_permutation_p(s, w, args.null_draws, rng)
            coh[a].append({"patient_id": meta["patient_id"], "section_id": meta["section_id"],
                           "moran_i": I, "moran_p": p})
            stab[a].append(A.contour_stability(s)[0])
        if "CXCL13" in gidx:
            bmask = A.contour_mask(scores["B"], 0.7)
            if 0 < bmask.sum() < len(bmask):
                signed = A.signed_hop_distance(w, bmask)
                prof = A.binned_profile(x_z[:, gidx["CXCL13"]], signed, PC_BINS)
                acc = pat_profiles.setdefault(meta["patient_id"], np.zeros((len(PC_BINS), 2)))
                for i, (v, n) in enumerate(prof):
                    if n and np.isfinite(v):
                        acc[i, 0] += v * n
                        acc[i, 1] += n

    axes_summary = {}
    for a in axis_defs:
        per_pat: dict[str, list[float]] = {}
        for r in coh[a]:
            per_pat.setdefault(r["patient_id"], []).append(r["moran_p"])
        n_coh = sum(1 for v in per_pat.values() if float(np.median(v)) <= 0.01)
        axes_summary[a] = {
            "n_genes": int(np.median(n_genes_used[a])),
            "patients_coherent": n_coh,
            "n_patients": len(per_pat),
            "moran_i_median": float(np.median([r["moran_i"] for r in coh[a]])),
            "moran_p_median": float(np.median([r["moran_p"] for r in coh[a]])),
            "contour_stability_mean": float(np.mean(stab[a])),
        }

    pc_peaks = {}
    for pat, acc in pat_profiles.items():
        prof = np.where(acc[:, 1] > 0, acc[:, 0] / np.maximum(acc[:, 1], 1), np.nan)
        pc_peaks[pat] = int(PC_BINS[int(np.nanargmax(prof))])
    peaks = np.array(list(pc_peaks.values()))
    positive_control = {
        "gene": "CXCL13", "axis": "B", "contour_quantile": 0.7,
        "bins": PC_BINS.tolist(),
        "n_patients": int(len(peaks)),
        "peak_hop_distribution": {str(int(k)): int(v) for k, v in
                                  zip(*np.unique(peaks, return_counts=True))},
        "n_peak_inside_or_edge_le1": int((peaks <= 1).sum()),
    }

    artifact = {
        "schema": "r16.axis_rulers.v1",
        "status": "EXPLORATORY_COMPLETE_NOT_CLAIM",
        "seed": C.SEED,
        "parameters": {"quantiles": list(A.QUANTILES), "null_type": "value_permutation",
                       "null_draws": args.null_draws, "knn_spatial": C.KNN_SPATIAL},
        "n_sections": n_sections,
        "axes": axes_summary,
        "positive_control": positive_control,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, ensure_ascii=False))

    rows = []
    for a, s in axes_summary.items():
        grade = "EXPLORATORY_REPRODUCED" if s["patients_coherent"] >= 2 else "DESCRIPTIVE_SINGLE_PATIENT"
        rows.append({
            "pattern_id": f"R16A_{a}", "tier": "1", "kind": "composition_axis",
            "axis_name": a, "lineage": "HTAN_VANDERBILT_CRC",
            "n_sections": n_sections, "n_patients": s["n_patients"],
            "resolution_support": "NA", "top_markers": f"module:{s['n_genes']}genes",
            "geometry": f"contour@{A.QUANTILES}",
            "effect_size": f"{s['moran_i_median']:.3f}",
            "null_type": "value_permutation", "null_draws": args.null_draws,
            "null_p": f"{s['moran_p_median']:.4f}", "fdr_q": "NA",
            "cross_patient_shape_corr": "NA", "residual_dimension": "NA",
            "uncertainty_ci95": "NA", "evidence_grade": grade,
            "notes": f"coherent_patients={s['patients_coherent']}/{s['n_patients']};"
                     f"contour_jaccard={s['contour_stability_mean']:.2f}",
        })
    args.out_rows.write_text(json.dumps(rows, indent=2, ensure_ascii=False))

    print(f"sections={n_sections}")
    for a, s in axes_summary.items():
        print(f"  {a:8s} coherent={s['patients_coherent']}/{s['n_patients']} "
              f"I={s['moran_i_median']:.2f} stab={s['contour_stability_mean']:.2f}")
    pc = positive_control
    print(f"  PC CXCL13~B: n={pc['n_patients']} peak<=+1: {pc['n_peak_inside_or_edge_le1']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
