#!/usr/bin/env python3
"""Score R-track joint clusters into registry rows (D-120).

Reads Seurat/Harmony label CSVs (global cluster ids across sections by
construction: no cross-section matching needed) + census cache (counts,
coords, genes). Per candidate: composition, pooled-log1p markers, per-section
Moran (200 draws), cross-arm agreement (best Jaccard vs every other arm/res:
numpy A/B + R S/H). Grade bar identical to Tier-2 (>=2 patients).
No split-half refit for R arms (cost): judge = cross-arm agreement +
later external replication. Disclosed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from r16 import census as C

ROOT = Path(__file__).resolve().parent.parent
ARMS = {
    "S": ("labels_armS.csv", ["SCT_025", "SCT_05", "SCT_10"]),
    "H": ("labels_armH.csv", ["HARM_025", "HARM_05", "HARM_10"]),
    "G": ("labels_armG.csv", ["GRAPHST_025", "GRAPHST_05", "GRAPHST_10"]),
}


def load_cache_index(cache_dir: Path):
    """section stem -> (counts csr, barcodes, genes, coords, meta)."""
    out = {}
    for stem in C.list_stems(cache_dir):
        mat, barcodes, genes, coords, meta = C.load_section(cache_dir, stem)
        out[stem] = {"mat": mat, "barcodes": list(barcodes), "genes": list(genes),
                     "coords": coords, "meta": meta}
    genes0 = next(iter(out.values()))["mat"].shape[1]
    return out, genes0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels-dir", type=Path, required=True)
    ap.add_argument("--cache", type=Path, default=ROOT / "infra/r16/census_cache_hvg10k/sections")
    ap.add_argument("--bridge-manifest", type=Path, default=ROOT / "infra/r16/seurat_bridge/manifest.json")
    ap.add_argument("--null-draws", type=int, default=C.NULL_DRAWS)
    ap.add_argument("--out-json", type=Path, default=ROOT / "infra/r16/joint_embed_R_20260915.json")
    ap.add_argument("--out-rows", type=Path, default=ROOT / "infra/r16/registry_rows_tier2r_20260915.json")
    args = ap.parse_args()
    rng = np.random.default_rng(C.SEED)

    sec_of_patient = {}
    for row in json.load(open(args.bridge_manifest)):
        sec_of_patient[row["stem"]] = (row["patient_id"], row["section_id"])
    cache, _ = load_cache_index(args.cache)
    # pooled log1p for markers (transient ~4.5GB)
    pooled_log, pooled_stem = [], []
    for stem, cd in cache.items():
        x = cd["mat"].toarray().astype(np.float32)
        np.log1p(x, out=x)
        pooled_log.append(x)
        pooled_stem.extend([stem] * x.shape[0])
    pooled_log = np.vstack(pooled_log)
    pooled_stem = np.array(pooled_stem)
    genes = cache[next(iter(cache))]["genes"]

    # all arm/res labelings, spot-aligned to pooled order
    labelings = {}  # (arm, col) -> (labels array, stems array)
    for arm, (fname, cols) in ARMS.items():
        if not (args.labels_dir / fname).exists():
            print(f"arm {arm}: {fname} absent, skipped", flush=True)
            continue
        df = pd.read_csv(args.labels_dir / fname, dtype=str)
        key = list(zip(df["section"].tolist(), df["orig_barcode"].tolist()))
        for col in cols:
            labmap = dict(zip(key, df[col].tolist()))
            labs = []
            ok = True
            for s, cd in cache.items():
                for b in cd["barcodes"]:
                    v = labmap.get((s, b))
                    if v is None:
                        ok = False
                        break
                    labs.append(v)
                if not ok:
                    break
            if not ok:
                raise SystemExit(f"barcode join failed for {arm}/{col}")
            labelings[(arm, col)] = (np.array(labs), pooled_stem.copy())
    print(f"arms loaded: {sorted(labelings)}", flush=True)

    # per-candidate characterization
    cands = []
    for (arm, col), (labs, stems_arr) in labelings.items():
        for c in sorted(set(labs.tolist())):
            idx = np.flatnonzero(labs == c)
            if len(idx) < C.MIN_CLUSTER_SIZE:
                continue
            secs = sorted(set(stems_arr[idx].tolist()))
            pats = sorted({sec_of_patient[s][0] for s in secs})
            sig = pooled_log[idx].mean(axis=0)
            # markers vs rest (pooled t-stat, reuse)
            rest = np.ones(len(labs), bool)
            rest[idx] = False
            markers = C.top_markers(pooled_log, np.where(np.isin(np.arange(len(labs)), idx), 1, 0),
                                    1, genes) if rest.any() else []
            morans, moranps, per_sec_n = [], [], {}
            for s in secs:
                m = idx[stems_arr[idx] == s]
                per_sec_n[s] = int(len(m))
                if len(m) < C.MIN_CLUSTER_SIZE:
                    continue
                cd = cache[s]
                # map global pooled positions back: pooled order == cache order
                mask = np.zeros(len(cd["barcodes"]))
                base = 0
                # find offset of section s in pooled array
                for stem2, cd2 in cache.items():
                    if stem2 == s:
                        break
                    base += len(cd2["barcodes"])
                local = m - base
                mask[local] = 1.0
                w = C.spatial_weights(cd["coords"], C.KNN_SPATIAL)
                I, p = C.moran_permutation_p(mask, w, args.null_draws, rng)
                morans.append(I)
                moranps.append(p)
            cands.append({
                "arm": arm, "col": col, "cluster_id": str(c),
                "n_spots": int(len(idx)), "n_sections": len(secs),
                "n_patients": len(pats), "patients": pats, "sections": secs,
                "per_section_n": per_sec_n, "sig": sig.astype(np.float32),
                "markers": markers,
                "moran_i_median": float(np.median(morans)) if morans else None,
                "moran_p_median": float(np.median(moranps)) if moranps else None,
                "moran_sections_tested": int(len(morans)),
            })
    print(f"candidates: {len(cands)}", flush=True)

    # cross-arm agreement: best Jaccard of member-spot sets vs every other labeling
    memb = {}
    for k, (labs, _) in labelings.items():
        d = {}
        for c in set(labs.tolist()):
            d[c] = set(np.flatnonzero(labs == c).tolist())
        memb[k] = d
    for r in cands:
        key = (r["arm"], {"SCT_025": "SCT_025", "SCT_05": "SCT_05", "SCT_10": "SCT_10",
                          "HARM_025": "HARM_025", "HARM_05": "HARM_05", "HARM_10": "HARM_10",
                          "GRAPHST_025": "GRAPHST_025", "GRAPHST_05": "GRAPHST_05",
                          "GRAPHST_10": "GRAPHST_10"}[r["col"]])
        mine = memb[key][r["cluster_id"]]
        best = {}
        for k2, d2 in memb.items():
            if k2 == key:
                continue
            bj = 0.0
            for c2, s2 in d2.items():
                u = len(mine | s2)
                j = len(mine & s2) / u if u else 0.0
                bj = max(bj, j)
            best[f"{k2[0]}:{k2[1]}"] = round(float(bj), 3)
        r["xarm_best_jaccard"] = best

    rows = []
    for r in cands:
        grade = ("EXPLORATORY_REPRODUCED" if r["n_patients"] >= 2
                 else "DESCRIPTIVE_SINGLE_PATIENT")
        pid = f"R16J{r['arm']}{r['col'].split('_')[1]}-{r['cluster_id']}"
        agree = ",".join(f"{k}={v}" for k, v in sorted(r["xarm_best_jaccard"].items()))
        rows.append({
            "pattern_id": pid, "tier": "2", "kind": "joint_pattern_group",
            "axis_name": "NA", "lineage": "HTAN_VANDERBILT_CRC",
            "n_sections": r["n_sections"], "n_patients": r["n_patients"],
            "resolution_support": r["col"],
            "top_markers": ";".join(r["markers"]), "geometry": "NA",
            "effect_size": (f"{r['moran_i_median']:.3f}" if r["moran_i_median"] is not None else "NA"),
            "null_type": "value_permutation", "null_draws": args.null_draws,
            "null_p": (f"{r['moran_p_median']:.4f}" if r["moran_p_median"] is not None else "NA"),
            "fdr_q": "NA", "cross_patient_shape_corr": "NA",
            "residual_dimension": "NA", "uncertainty_ci95": "NA",
            "evidence_grade": grade,
            "notes": f"arm={r['arm']};no_split_half_refit(cost);xarm_agree:{agree}",
        })
    slim = []
    for r in cands:
        d = {k: v for k, v in r.items() if k != "sig"}
        slim.append(d)
    artifact = {"schema": "r16.joint_embed_R.v1",
                "status": "EXPLORATORY_COMPLETE_NOT_CLAIM", "seed": C.SEED,
                "parameters": {"arms": sorted({a for a, _ in labelings}), "resolutions": [0.25, 0.5, 1.0],
                               "min_cluster_size": C.MIN_CLUSTER_SIZE,
                               "null_draws": args.null_draws,
                               "split_half": "not_refit_for_R_arms_judge_is_xarm_agreement_plus_external"},
                "n_candidates": len(slim), "candidates": slim}
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, ensure_ascii=False))
    args.out_rows.write_text(json.dumps(rows, indent=2, ensure_ascii=False))
    n_rep = sum(1 for r in rows if r["evidence_grade"] == "EXPLORATORY_REPRODUCED")
    print(f"candidates={len(rows)} reproduced(>=2pat)={n_rep}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
