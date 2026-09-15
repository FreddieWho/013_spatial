#!/usr/bin/env python3
"""R-16 step 1: composition pattern census (design doc v1.3 §3).

For each cached section: log1p -> z-score -> PCA30 -> kNN15 -> Leiden at three
resolutions (min cluster 20). Then: cluster signatures, cross-resolution
Jaccard stability, Moran's I with value-permutation null (200 draws), and
cross-section cosine matching (>=0.75 connected components) into pattern
groups. Universality counted by PATIENT. Outputs census JSON + registry TSV.

Exploratory grade. No claims. Deterministic (seed 20260914).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from r16 import axes as A
from r16 import census as C

ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=ROOT / "infra/r16/census_cache_hvg10k/sections")
    ap.add_argument("--out-json", type=Path,
                    default=ROOT / "infra/r16/composition_pattern_census_20260914.json")
    ap.add_argument("--out-clusters", type=Path,
                    default=ROOT / "infra/r16/census_clusters_20260914.npz")
    ap.add_argument("--null-draws", type=int, default=C.NULL_DRAWS)
    ap.add_argument("--marker-proxy", type=Path,
                    default=ROOT / "infra/r04/marker_proxy_combined.json")
    ap.add_argument("--out-rows", type=Path,
                    default=ROOT / "infra/r16/registry_rows_tier2_20260914.json")
    args = ap.parse_args()

    rng = np.random.default_rng(C.SEED)
    stems = C.list_stems(args.cache)
    if not stems:
        raise SystemExit(f"empty cache: {args.cache}")

    clusters = []  # every kept cluster across sections/resolutions
    n_fail = 0
    for stem in stems:
        try:
            mat, barcodes, genes, coords, meta = C.load_section(args.cache, stem)
            x_z = C.log1p_zscore(mat)
            x_log = np.log1p(mat.toarray().astype(np.float32))
            scores = C.pca_scores(x_z)
            g = C.knn_igraph(scores, C.KNN_EXPR)
            w_sp = C.spatial_weights(coords, C.KNN_SPATIAL)
            labels_by_res = {}
            for res in C.RESOLUTIONS:
                labels_by_res[res] = C.filter_min_size(
                    C.leiden_labels(g, res), C.MIN_CLUSTER_SIZE)
            for res in C.RESOLUTIONS:
                labels = labels_by_res[res]
                sigs = C.cluster_signatures(x_z, labels)
                for c, sig in sigs.items():
                    mask = labels == c
                    # stability: best Jaccard against adjacent resolutions
                    jac = []
                    for r2 in C.RESOLUTIONS:
                        if r2 == res:
                            continue
                        lab2 = labels_by_res[r2]
                        best = 0.0
                        for c2 in np.unique(lab2):
                            if c2 < 0:
                                continue
                            best = max(best, C.jaccard(mask, lab2 == c2))
                        jac.append(best)
                    I, p = C.moran_permutation_p(mask.astype(float), w_sp,
                                                 args.null_draws, rng)
                    clusters.append({
                        "stem": stem,
                        "section_id": meta["section_id"],
                        "patient_id": meta["patient_id"],
                        "resolution": res,
                        "cluster_id": int(c),
                        "n_spots": int(mask.sum()),
                        "signature": sig,
                        "markers": C.top_markers(x_log, labels, c, genes),
                        "stability_jaccard_min": float(min(jac)),
                        "moran_i": I,
                        "moran_p": p,
                        "genes_ref": genes,
                    })
        except Exception as e:  # fail-closed per section, disclosed
            clusters.append({"stem": stem, "error": f"{type(e).__name__}: {e}"})
            n_fail += 1

    ok = [c for c in clusters if "error" not in c]
    errs = [c for c in clusters if "error" in c]

    # cross-section matching: average-linkage hierarchical on cosine distance
    # (v2.0: kills single-linkage chaining; cut at 1 - MATCH_COSINE)
    n = len(ok)
    if n:
        sig_mat = np.vstack([c["signature"] for c in ok]).astype(np.float32)
        norm = sig_mat / np.linalg.norm(sig_mat, axis=1, keepdims=True)
        sim = norm @ norm.T
        np.fill_diagonal(sim, 1.0)
        dist = 1 - np.clip(sim, -1, 1)
        Z = linkage(squareform(dist, checks=False), method="average")
        lab = fcluster(Z, t=1 - C.MATCH_COSINE, criterion="distance")
        groups: dict[int, list[int]] = {}
        for i, g in enumerate(lab):
            groups.setdefault(int(g), []).append(i)
    else:
        sim = np.zeros((0, 0))
        groups = {}

    # Tier-1 axis alignment per group (v2.0): mean signature scored on each axis
    proxy = json.load(open(args.marker_proxy))
    axis_defs = {k: v["voted_genes"] for k, v in proxy["classes"].items()
                 if k not in A.RETIRED_AXES}
    axis_defs["Plasma"] = A.PLASMA_GENES
    genes0 = ok[0]["genes_ref"] if ok else []
    gidx0 = {g: i for i, g in enumerate(genes0)}

    registry_rows = []
    row_dicts = []
    group_summaries = []
    gid = 0
    for root, members in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        gid += 1
        ms = [ok[i] for i in members]
        patients = sorted({m["patient_id"] for m in ms})
        sections = sorted({m["section_id"] for m in ms})
        resolutions = sorted({m["resolution"] for m in ms})
        markers: list[str] = []
        for m in ms:
            for gname in m["markers"]:
                if gname not in markers:
                    markers.append(gname)
                if len(markers) >= 10:
                    break
            if len(markers) >= 10:
                break
        stable_frac = float(np.mean([m["stability_jaccard_min"] >= C.STABILITY_JACCARD
                                     for m in ms]))
        moran_med = float(np.median([m["moran_i"] for m in ms]))
        moran_p_med = float(np.median([m["moran_p"] for m in ms]))
        # purity: min pairwise cosine within group (v2.0)
        if len(members) > 1:
            sub = sim[np.ix_(members, members)].copy()
            np.fill_diagonal(sub, 1.0)
            purity = float(sub.min())
        else:
            purity = 1.0
        # Tier-1 axis alignment: mean signature scored on each axis
        mean_sig = np.mean([m["signature"] for m in ms], axis=0)
        axis_sc = {}
        for aname, agenes in axis_defs.items():
            cols = [gidx0[g] for g in agenes if g in gidx0]
            axis_sc[aname] = float(mean_sig[cols].mean()) if cols else 0.0
        best_axis = max(axis_sc, key=axis_sc.get)
        best_axis_score = axis_sc[best_axis]
        explained = best_axis_score >= 0.5
        if purity < 0.5:
            grade = "DESCRIPTIVE_SINGLE_PATIENT"
            purity_note = "LOW_PURITY_MIXED_GROUP"
        else:
            grade = ("EXPLORATORY_REPRODUCED" if len(patients) >= 2
                     else "DESCRIPTIVE_SINGLE_PATIENT")
            purity_note = f"purity={purity:.2f}"
        pid = f"R16P{gid:03d}"
        group_summaries.append({
            "pattern_id": pid,
            "n_clusters": len(ms),
            "n_sections": len(sections),
            "n_patients": len(patients),
            "patients": patients,
            "resolutions_present": resolutions,
            "resolution_support": len(resolutions),
            "stable_cluster_fraction": stable_frac,
            "moran_i_median": moran_med,
            "moran_p_median": moran_p_med,
            "purity": purity,
            "best_axis": best_axis,
            "best_axis_score": best_axis_score,
            "axis_explained": explained,
            "top_markers": markers,
            "evidence_grade": grade,
        })
        notes = (f"stable_frac={stable_frac:.2f};{purity_note};"
                 f"best_axis={best_axis}({best_axis_score:.2f})"
                 f"{'_explained' if explained else '_UNEXPLAINED'}")
        registry_rows.append([
            pid, "2", "pattern_group", "NA", "HTAN_VANDERBILT_CRC",
            str(len(sections)), str(len(patients)), f"{len(resolutions)}/3",
            ";".join(markers), "NA", f"{moran_med:.3f}",
            "value_permutation", str(args.null_draws), f"{moran_p_med:.4f}", "NA",
            "NA", "NA", "NA", grade, notes,
        ])
        row_dicts.append({
            "pattern_id": pid, "tier": "2", "kind": "pattern_group",
            "axis_name": "NA", "lineage": "HTAN_VANDERBILT_CRC",
            "n_sections": len(sections), "n_patients": len(patients),
            "resolution_support": f"{len(resolutions)}/3",
            "top_markers": ";".join(markers), "geometry": "NA",
            "effect_size": f"{moran_med:.3f}", "null_type": "value_permutation",
            "null_draws": args.null_draws, "null_p": f"{moran_p_med:.4f}",
            "fdr_q": "NA", "cross_patient_shape_corr": "NA",
            "residual_dimension": f"axis_explained={explained}",
            "uncertainty_ci95": "NA", "evidence_grade": grade, "notes": notes,
        })

    # persist cluster signatures + metadata for downstream diagnosis/reuse
    if ok:
        sig_mat = np.vstack([c["signature"] for c in ok]).astype(np.float32)
        np.savez_compressed(
            args.out_clusters,
            signatures=sig_mat,
            section_idx=np.array([stems.index(c["stem"]) for c in ok], dtype=np.int32),
            resolutions=np.array([c["resolution"] for c in ok], dtype=np.float32),
            n_spots=np.array([c["n_spots"] for c in ok], dtype=np.int32),
            moran_i=np.array([c["moran_i"] for c in ok], dtype=np.float32),
            moran_p=np.array([c["moran_p"] for c in ok], dtype=np.float32),
            stability=np.array([c["stability_jaccard_min"] for c in ok], dtype=np.float32),
        )
        (args.out_clusters.parent / (args.out_clusters.stem + ".json")).write_text(
            json.dumps([{"section_id": c["section_id"], "patient_id": c["patient_id"],
                         "resolution": c["resolution"], "cluster_id": c["cluster_id"],
                         "markers": c["markers"]} for c in ok], ensure_ascii=False))

    artifact = {
        "schema": "r16.composition_pattern_census.v2",
        "status": "EXPLORATORY_COMPLETE_NOT_CLAIM",
        "seed": C.SEED,
        "parameters": {
            "resolutions": list(C.RESOLUTIONS), "min_cluster_size": C.MIN_CLUSTER_SIZE,
            "pca_dims": C.PCA_DIMS, "knn_expr": C.KNN_EXPR, "knn_spatial": C.KNN_SPATIAL,
            "match_cosine": C.MATCH_COSINE, "stability_jaccard": C.STABILITY_JACCARD,
            "null_draws": args.null_draws, "null_type": "value_permutation",
            "matching_method": "average_linkage_hierarchical",
            "purity_floor": 0.5, "axis_explained_floor": 0.5,
        },
        "n_sections": len(stems),
        "n_sections_failed": n_fail,
        "section_errors": errs,
        "n_clusters_total": n,
        "n_pattern_groups": len(group_summaries),
        "pattern_groups": group_summaries,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, ensure_ascii=False))
    args.out_rows.write_text(json.dumps(row_dicts, indent=2, ensure_ascii=False))

    print(f"sections={len(stems)} failed={n_fail} clusters={n} groups={len(group_summaries)}")
    for g in group_summaries[:20]:
        print(f"  {g['pattern_id']}: pat={g['n_patients']:2d} sec={g['n_sections']:2d} "
              f"res={g['resolution_support']}/3 I={g['moran_i_median']:.2f} "
              f"purity={g['purity']:.2f} axis={g['best_axis']}({g['best_axis_score']:.2f}) "
              f"{g['evidence_grade']} markers={','.join(g['top_markers'][:5])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
