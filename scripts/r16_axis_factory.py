#!/usr/bin/env python3
"""R-16 Tier-1 v2 ruler factory driver (design v2.1, D-118).

Per section per lambda: neighbor-augment -> PCA top-10 candidate axes ->
LISA tile coherence (999-draw value-permutation + BH-FDR) -> novelty vs the
6 working axes (|Pearson r| >= 0.7 KNOWN). Then cross-patient matching on
|cosine| of loading vectors (average linkage, cut 0.25) with recurrence
counted over tile_ok members only.

Outputs: artifact JSON (tracked) + loadings NPZ (local-only) + tier1b rows
JSON (tracked, merged by scripts/r16_build_registry.py).
Exploratory grade; no claims.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from r16 import axes as A
from r16 import axis_factory as F
from r16 import census as C

ROOT = Path(__file__).resolve().parent.parent
N_PC = F.N_PC
KNOW_R_THRESHOLD = 0.7


def load_axis_defs(marker_proxy: Path, union_genes: Path, rank_file: Path, top_n: int):
    proxy = json.load(open(marker_proxy))
    union = json.load(open(union_genes))
    rank = np.load(rank_file)
    top = {union[i] for i in np.flatnonzero(rank < top_n)}
    defs = {k: [g for g in v["voted_genes"] if g in top]
            for k, v in proxy["classes"].items()
            if k not in A.RETIRED_AXES}
    defs["Plasma"] = [g for g in A.PLASMA_GENES if g in top]
    return defs


def pick_pilot_stems(cache, stems: list[str]) -> list[str]:
    """5 stems spanning the spot-count range (deterministic)."""
    sizes = []
    for stem in stems:
        js = json.loads((cache / f"{stem}.json").read_text())
        z = np.load(cache / f"{stem}.npz")
        n = int(tuple(z["shape"])[0])
        sizes.append((n, stem, js))
    sizes.sort()
    picks = [sizes[i][1] for i in
             sorted({0, len(sizes) // 4, len(sizes) // 2, 3 * len(sizes) // 4, len(sizes) - 1})]
    return picks


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=ROOT / "infra/r16/census_cache_hvg10k/sections")
    ap.add_argument("--marker-proxy", type=Path, default=ROOT / "infra/r04/marker_proxy_combined.json")
    ap.add_argument("--union-genes", type=Path, default=ROOT / "infra/r16/hvg_union_genes.json")
    ap.add_argument("--rank-file", type=Path, default=ROOT / "infra/r16/hvg_final_rank_20260914.npy")
    ap.add_argument("--top-n", type=int, default=10000)
    ap.add_argument("--lambdas", type=float, nargs="+", default=list(F.LAMBDAS))
    ap.add_argument("--k-spatial", type=int, default=F.K_PRIMARY)
    ap.add_argument("--radius-check", action="store_true",
                    help="also run k=30 on processed sections and report agreement")
    ap.add_argument("--pilot", action="store_true", help="5 size-spanning sections only")
    ap.add_argument("--sections", type=str, nargs="*", default=None)
    ap.add_argument("--tag", type=str, default="20260915",
                    help="output filename tag (pilot uses pilot_ prefix)")
    args = ap.parse_args()

    rng = np.random.default_rng(C.SEED)
    stems = C.list_stems(args.cache)
    if args.sections:
        stems = [s for s in stems if s in set(args.sections)]
    if args.pilot:
        stems = pick_pilot_stems(args.cache, stems)
        tag = f"pilot_{args.tag}"
    else:
        tag = args.tag
    axis_defs = load_axis_defs(args.marker_proxy, args.union_genes, args.rank_file, args.top_n)

    out_dir = ROOT / "infra/r16"
    loadings_path = out_dir / f"axis_factory_loadings_{tag}.npz"
    artifact_path = out_dir / f"axis_factory_{tag}.json"

    t0 = time.time()
    gene_order = None
    candidates = []   # one row per (section, lambda, pc)
    load_list, load_keys = [], []
    for si, stem in enumerate(stems):
        mat, barcodes, genes, coords, meta = C.load_section(args.cache, stem)
        if gene_order is None:
            gene_order = list(genes)
        elif list(genes) != gene_order:
            raise SystemExit(f"gene order drift in {stem}; abort (cache contract)")
        x_log = np.log1p(mat.toarray().astype(np.float32))
        x_z = C.log1p_zscore(mat)
        gidx = {g: i for i, g in enumerate(genes)}
        w = C.spatial_weights(coords, C.KNN_SPATIAL)
        # working-axis fields for novelty (computed once on raw x_z)
        ax_fields = {}
        for aname, agenes in axis_defs.items():
            s, _ = A.axis_scores(x_z, gidx, agenes)
            ax_fields[aname] = s
        # neighborhood mean on log1p, then z-scored like X
        kw = F.spatial_mean_weights(coords, args.k_spatial)
        m_z = F.zscore_columns(np.asarray((kw @ x_log)))
        for lam in args.lambdas:
            z_aug = F.augment(x_z, m_z.astype(np.float32), lam)
            scores, loadings, vratio = F.candidate_pcs(z_aug, N_PC)
            for pc in range(loadings.shape[1]):
                s = scores[:, pc]
                lisa = F.lisa_tile_stats(s, w, F.LISA_DRAWS, rng)
                rs = {an: float(abs(np.corrcoef(s, f)[0, 1]))
                      for an, f in ax_fields.items()}
                best = max(rs, key=rs.get)
                candidates.append({
                    "stem": stem, "section_id": meta["section_id"],
                    "patient_id": meta["patient_id"], "lambda": lam, "pc": pc,
                    "var_ratio": float(vratio[pc]),
                    "tile_ok": lisa["tile_ok"], "n_sig": lisa["n_sig"],
                    "frac_sig": lisa["frac_sig"], "n_hotspots": lisa["n_hotspots"],
                    "max_hotspot": lisa["max_hotspot"],
                    "best_axis": best, "best_r": rs[best],
                    "known": bool(rs[best] >= KNOW_R_THRESHOLD),
                })
                load_list.append(loadings[:, pc].astype(np.float32))
                load_keys.append(f"{stem}::lam{lam}::pc{pc}")
        print(f"[{si + 1}/{len(stems)}] {stem} done ({time.time() - t0:.0f}s)", flush=True)

    np.savez_compressed(loadings_path,
                        loadings=np.vstack(load_list),
                        keys=np.array(load_keys),
                        genes=np.array(gene_order))
    print(f"saved loadings {loadings_path} ({time.time() - t0:.0f}s)")

    # ---- cross-patient matching on |cosine| of loadings
    L = np.vstack(load_list).astype(np.float64)
    nrm = L / np.linalg.norm(L, axis=1, keepdims=True)
    sim = nrm @ nrm.T
    np.fill_diagonal(sim, 1.0)
    Z = linkage(squareform(1 - np.clip(sim, -1, 1), checks=False), method="average")
    lab = fcluster(Z, t=1 - C.MATCH_COSINE, criterion="distance")
    groups: dict[int, list[int]] = {}
    for i, g in enumerate(lab):
        groups.setdefault(int(g), []).append(i)

    cand_by_idx = candidates
    rows, group_recs = [], []
    gid = 0
    for g in sorted(groups, key=lambda k: -len(groups[k])):
        mem = groups[g]
        gid += 1
        pats = sorted({cand_by_idx[i]["patient_id"] for i in mem})
        rec_pats = sorted({cand_by_idx[i]["patient_id"] for i in mem
                           if cand_by_idx[i]["tile_ok"]})
        lams = sorted({cand_by_idx[i]["lambda"] for i in mem})
        sub = sim[np.ix_(mem, mem)].copy()
        np.fill_diagonal(sub, 1.0)
        purity = float(sub.min()) if len(mem) > 1 else 1.0
        new_frac = float(np.mean([not cand_by_idx[i]["known"] for i in mem]))
        axes_hit = [cand_by_idx[i]["best_axis"] for i in mem if cand_by_idx[i]["known"]]
        top_axis = max(set(axes_hit), key=axes_hit.count) if axes_hit else "NONE"
        mean_load = np.mean([L[i] for i in mem], axis=0)
        markers = F.top_loading_genes(mean_load, gene_order, 10)
        var_med = float(np.median([cand_by_idx[i]["var_ratio"] for i in mem]))
        tile_frac = float(np.mean([cand_by_idx[i]["tile_ok"] for i in mem]))
        if purity < 0.5:
            grade = "DESCRIPTIVE_SINGLE_PATIENT"
            note_extra = "LOW_PURITY_MIXED_GROUP"
        else:
            grade = ("EXPLORATORY_REPRODUCED" if len(rec_pats) >= 2
                     else "DESCRIPTIVE_SINGLE_PATIENT")
            note_extra = f"purity={purity:.2f}"
        pid = f"R16B_{gid:03d}"
        novelty = "NEW" if new_frac >= 0.5 else f"KNOWN:{top_axis}"
        group_recs.append({
            "pattern_id": pid, "n_members": len(mem),
            "n_sections": len({cand_by_idx[i]["section_id"] for i in mem}),
            "n_patients": len(pats), "recurrence_patients": len(rec_pats),
            "patients": pats, "lambda_support": lams,
            "purity": purity, "novelty": novelty, "new_frac": new_frac,
            "top_axis": top_axis,
            "var_ratio_median": var_med, "tile_frac": tile_frac,
            "top_markers": markers, "evidence_grade": grade,
        })
        rows.append({
            "pattern_id": pid, "tier": "1b", "kind": "discovered_axis",
            "axis_name": novelty, "lineage": "HTAN_VANDERBILT_CRC",
            "n_sections": len({cand_by_idx[i]["section_id"] for i in mem}),
            "n_patients": len(pats),
            "resolution_support": f"lambda:{len(lams)}/4",
            "top_markers": ";".join(markers), "geometry": "NA",
            "effect_size": f"{var_med:.4f}",
            "null_type": "value_permutation_LISA", "null_draws": F.LISA_DRAWS,
            "null_p": "NA", "fdr_q": f"{F.FDR_ALPHA}",
            "cross_patient_shape_corr": f"{purity:.3f}",
            "residual_dimension": f"new_frac={new_frac:.2f}",
            "uncertainty_ci95": "NA", "evidence_grade": grade,
            "notes": f"{note_extra};tile_frac={tile_frac:.2f};"
                     f"recurrence_patients={len(rec_pats)}",
        })

    artifact = {
        "schema": "r16.axis_factory.v1",
        "status": "EXPLORATORY_COMPLETE_NOT_CLAIM",
        "seed": C.SEED,
        "parameters": {"lambdas": list(args.lambdas), "k_spatial": args.k_spatial,
                       "n_pc": N_PC, "lisa_draws": F.LISA_DRAWS,
                       "fdr_alpha": F.FDR_ALPHA,
                       "min_hotspot": F.MIN_HOTSPOT_SIZE,
                       "match": "average_linkage_|cos|_0.75",
                       "novelty_r": KNOW_R_THRESHOLD,
                       "radius_check": bool(args.radius_check)},
        "n_sections": len(stems),
        "n_candidates": len(candidates),
        "n_groups": len(group_recs),
        "candidates": candidates,
        "groups": group_recs,
    }
    artifact_path.write_text(json.dumps(artifact, indent=2, ensure_ascii=False))
    (out_dir / f"registry_rows_tier1b_{tag}.json").write_text(
        json.dumps(rows, indent=2, ensure_ascii=False))

    n_rep = sum(1 for g in group_recs if g["evidence_grade"] == "EXPLORATORY_REPRODUCED")
    n_new = sum(1 for g in group_recs if g["novelty"] == "NEW" and
                g["evidence_grade"] == "EXPLORATORY_REPRODUCED")
    print(f"candidates={len(candidates)} groups={len(group_recs)} "
          f"reproduced={n_rep} new_reproduced={n_new} ({time.time() - t0:.0f}s)")

    if args.radius_check:
        agree = radius_check_report(args, stems, axis_defs)
        artifact["radius_check"] = agree
        artifact_path.write_text(json.dumps(artifact, indent=2, ensure_ascii=False))
        print("radius_check:", json.dumps(agree, ensure_ascii=False))
    return 0


def radius_check_report(args, stems, axis_defs) -> dict:
    """Rerun the same sections at k=30; report loading agreement per (section, lambda, pc).

    Escalation rule (disclosed): if median best-match |cos| < 0.9, do NOT
    silently proceed — report only.
    """
    rng = np.random.default_rng(C.SEED + 1)
    cosines = []
    for stem in stems:
        mat, barcodes, genes, coords, meta = C.load_section(args.cache, stem)
        x_log = np.log1p(mat.toarray().astype(np.float32))
        x_z = C.log1p_zscore(mat)
        m15 = F.zscore_columns(np.asarray((F.spatial_mean_weights(coords, 15) @ x_log)))
        m30 = F.zscore_columns(np.asarray((F.spatial_mean_weights(coords, 30) @ x_log)))
        for lam in args.lambdas:
            _, l15, _ = F.candidate_pcs(F.augment(x_z, m15.astype(np.float32), lam), N_PC)
            _, l30, _ = F.candidate_pcs(F.augment(x_z, m30.astype(np.float32), lam), N_PC)
            for pc in range(l15.shape[1]):
                best = max(F.abs_cosine(l15[:, pc], l30[:, q]) for q in range(l30.shape[1]))
                cosines.append(best)
    cosines = np.array(cosines)
    return {"median_best_match_cos": float(np.median(cosines)),
            "frac_above_0.9": float((cosines >= 0.9).mean()),
            "n_compared": int(len(cosines))}


if __name__ == "__main__":
    raise SystemExit(main())
