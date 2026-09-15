#!/usr/bin/env python3
"""R-16 joint-embedding discovery track (D-119).

Pooled single-cell-style integration over HTAN Vanderbilt sections ONLY.
  Arm A: pool log1p -> ComBat(batch=section) -> z-score -> PCA30 -> kNN15
         -> Leiden x3 (min 20).
  Arm B: same on BANKSY-style [own, neighbor-mean] concat (spatial kNN k=6,
         self excluded, row-normalized; equal-weight mixing disclosed).
Internal judge: split-half patients (15/15, seeded) rerun independently per
arm (ComBat REFIT inside each half: no leakage); full-run clusters annotated
by best cosine to EACH half (>=0.75 both halves = split-supported).
Cross-run matching uses RAW log1p-space signatures (common space; corrected
spaces are run-specific and not comparable).
Validation lineages (ST-CRC/USZ) never touched. Exploratory; no claims.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from r16 import census as C
from r16 import joint as J

ROOT = Path(__file__).resolve().parent.parent
MATCH_COS = C.MATCH_COSINE


def load_all(cache_dir: Path):
    stems = C.list_stems(cache_dir)
    mats, coords_l, sec_ids, pat_ids, genes0 = [], [], [], [], None
    for stem in stems:
        mat, barcodes, genes, coords, meta = C.load_section(cache_dir, stem)
        if genes0 is None:
            genes0 = list(genes)
        elif list(genes) != genes0:
            raise SystemExit(f"gene order drift in {stem}; abort (cache contract)")
        mats.append(mat.tocsr())
        coords_l.append(coords)
        sec_ids.append(meta["section_id"])
        pat_ids.append(meta["patient_id"])
    return stems, mats, coords_l, sec_ids, pat_ids, genes0


def cluster_corrected(x_corr_z: np.ndarray, rng_seed: int = C.SEED):
    scores = C.pca_scores(x_corr_z)
    g = C.knn_igraph(scores, C.KNN_EXPR)
    out = {}
    for res in C.RESOLUTIONS:
        out[res] = C.filter_min_size(C.leiden_labels(g, res), C.MIN_CLUSTER_SIZE)
    return out


def describe_clusters(labels_by_res, x_raw_log, spot_sec, spot_pat, sec_ids,
                      coords_l, genes, arm, draws, rng):
    """Characterize clusters: signatures (raw+corrected space handled by caller
    via sig_z), markers, composition, per-section Moran.
    spot_sec/spot_pat are per-SPOT section/patient id arrays."""
    recs = []
    for res, labels in labels_by_res.items():
        for c in sorted(set(labels.tolist())):
            if c < 0:
                continue
            idx = np.flatnonzero(labels == c)
            if len(idx) < C.MIN_CLUSTER_SIZE:
                continue
            secs = sorted(set(spot_sec[idx].tolist()))
            pats = sorted(set(spot_pat[idx].tolist()))
            sig_raw = x_raw_log[idx].mean(axis=0)
            markers = C.top_markers(x_raw_log, labels, c, genes)
            # Moran per contributing section (>=20 spots of this cluster)
            morans, moranps = [], []
            per_sec_n = {}
            for s in secs:
                m = idx[spot_sec[idx] == s]
                per_sec_n[s] = int(len(m))
                if len(m) < C.MIN_CLUSTER_SIZE:
                    continue
                si = sec_ids.index(s)
                w = C.spatial_weights(coords_l[si], C.KNN_SPATIAL)
                # map global idx to local: positions of section spots
                loc = np.flatnonzero(spot_sec == s)
                pos = {g: k for k, g in enumerate(loc)}
                mask = np.zeros(len(loc))
                mask[[pos[i] for i in m]] = 1.0
                I, p = C.moran_permutation_p(mask, w, draws, rng)
                morans.append(I)
                moranps.append(p)
            recs.append({
                "arm": arm, "resolution": res, "cluster_id": int(c),
                "n_spots": int(len(idx)), "n_sections": len(secs),
                "n_patients": len(pats), "patients": pats,
                "sections": secs, "per_section_n": per_sec_n,
                "sig_raw": np.asarray(sig_raw, dtype=np.float32),
                "markers": markers,
                "moran_i_median": float(np.median(morans)) if morans else None,
                "moran_p_median": float(np.median(moranps)) if moranps else None,
                "moran_sections_tested": int(len(morans)),
            })
    return recs


def run_arm(x_feat_log: np.ndarray, batches: np.ndarray, sec_ids: list,
            pat_ids: list, coords_l: list, genes: list, arm: str,
            draws: int, rng: np.random.Generator, t0: float, log_prefix: str):
    x_corr = J.combat_correct(x_feat_log, batches)
    x_z = J.zscore_pooled(x_corr)
    del x_corr
    labels_by_res = cluster_corrected(x_z, C.SEED)
    # raw-space matrix for signatures/markers: own-log only (first block if augmented)
    return labels_by_res, x_z


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=ROOT / "infra/r16/census_cache_hvg10k/sections")
    ap.add_argument("--arms", nargs="+", default=["A", "B"], choices=["A", "B"])
    ap.add_argument("--tag", type=str, default="20260915")
    ap.add_argument("--null-draws", type=int, default=C.NULL_DRAWS)
    ap.add_argument("--max-sections", type=int, default=None)
    ap.add_argument("--skip-halves", action="store_true")
    ap.add_argument("--out-json", type=Path, default=None)
    ap.add_argument("--out-rows", type=Path, default=None)
    args = ap.parse_args()
    t0 = time.time()
    out_json = args.out_json or (ROOT / f"infra/r16/joint_embed_{args.tag}.json")
    out_rows = args.out_rows or (ROOT / f"infra/r16/registry_rows_tier2j_{args.tag}.json")
    rng = np.random.default_rng(C.SEED)

    stems, mats, coords_l, sec_ids, pat_ids, genes = load_all(args.cache)
    assert len(set(sec_ids)) == len(sec_ids), "section_id must be unique"
    if args.max_sections:
        keep = sorted(range(len(stems)))[:args.max_sections]
        mats = [mats[i] for i in keep]
        coords_l = [coords_l[i] for i in keep]
        sec_ids = [sec_ids[i] for i in keep]
        pat_ids = [pat_ids[i] for i in keep]
        stems = [stems[i] for i in keep]
    print(f"loaded {len(stems)} sections, {sum(m.shape[0] for m in mats)} spots", flush=True)

    ns = [int(m.shape[0]) for m in mats]
    X_log = J.pool_log1p(mats)
    del mats
    batches = np.repeat(np.array(sec_ids), ns)
    spot_pat_all = np.repeat(np.array(pat_ids), ns)
    # section row offsets from true spot counts
    bounds, off = [], 0
    for n in ns:
        bounds.append((off, off + n))
        off += n
    del ns

    # Arm B neighbor means (per section, then stack)
    N_log = None
    if "B" in args.arms:
        parts = []
        for s in sec_ids:
            a, b = bounds[sec_ids.index(s)]
            parts.append(J.neighbor_mean_log1p(X_log[a:b], coords_l[sec_ids.index(s)]))
        N_log = np.vstack(parts).astype(np.float32)
        print(f"neighbor means stacked {N_log.shape} ({time.time()-t0:.0f}s)", flush=True)

    results = {}
    for arm in args.arms:
        Xf = X_log if arm == "A" else np.hstack([X_log, N_log]).astype(np.float32)
        labels_by_res, x_z = run_arm(Xf, batches, sec_ids, pat_ids, coords_l,
                                     genes, arm, args.null_draws, rng, t0, arm)
        spot_sec = batches
        recs = describe_clusters(labels_by_res, X_log, spot_sec, spot_pat_all, sec_ids,
                                 coords_l, genes, arm, args.null_draws, rng)
        results[arm] = {"labels": labels_by_res, "recs": recs, "x_z": x_z}
        ncl = {res: int(np.unique(lab[lab >= 0]).size) for res, lab in labels_by_res.items()}
        print(f"arm {arm}: clusters_per_res={ncl} kept_recs={len(recs)} "
              f"({time.time()-t0:.0f}s)", flush=True)
        del x_z

    # ---- split-half internal judge (ComBat refit inside each half)
    half_info = {"halves": None, "done": False}
    if not args.skip_halves:
        uniq_pats = sorted(set(pat_ids))
        half_a, half_b = J.patient_halves(uniq_pats, C.SEED)
        half_info["halves"] = {"A": sorted(half_a), "B": sorted(half_b)}
        for arm in args.arms:
            Xf = X_log if arm == "A" else np.hstack([X_log, N_log]).astype(np.float32)
            half_sigs = {}
            for hname, hset in (("A", half_a), ("B", half_b)):
                keep = np.array([p in hset for p in spot_pat_all])
                Xh = Xf[keep]
                bh = batches[keep]
                xcorr = J.combat_correct(Xh, bh)
                xz = J.zscore_pooled(xcorr)
                del xcorr
                lab = cluster_corrected(xz, C.SEED)
                del xz
                # raw-space signatures for cross-half matching (common space)
                Xraw = X_log[keep]
                sigs = []
                for res, labels in lab.items():
                    for c in sorted(set(labels.tolist())):
                        if c < 0:
                            continue
                        idx = np.flatnonzero(labels == c)
                        if len(idx) >= C.MIN_CLUSTER_SIZE:
                            sigs.append(Xraw[idx].mean(axis=0))
                half_sigs[hname] = np.vstack(sigs).astype(np.float32) if sigs else np.zeros((0, X_log.shape[1]), np.float32)
                print(f"arm {arm} half {hname}: {len(sigs)} clusters ({time.time()-t0:.0f}s)", flush=True)
            results[arm]["half_sigs"] = half_sigs
            # mutual best-match between halves
            if len(half_sigs["A"]) and len(half_sigs["B"]):
                ba = J.best_match_table(half_sigs["A"], half_sigs["B"])
                bb = J.best_match_table(half_sigs["B"], half_sigs["A"])
                results[arm]["half_match"] = {
                    "A_to_B_median": float(np.median(ba)),
                    "B_to_A_median": float(np.median(bb)),
                    "A_frac_ge075": float((ba >= MATCH_COS).mean()),
                    "B_frac_ge075": float((bb >= MATCH_COS).mean()),
                }
            else:
                results[arm]["half_match"] = None
        half_info["done"] = True

    # ---- annotate full-run clusters with split-half support + build rows
    all_rows, all_cands = [], []
    for arm in args.arms:
        hs = results[arm].get("half_sigs")
        for r in results[arm]["recs"]:
            sup = {"A": None, "B": None}
            if hs is not None and len(hs["A"]) and len(hs["B"]):
                s = r["sig_raw"].astype(np.float64)
                for hname in ("A", "B"):
                    H = hs[hname].astype(np.float64)
                    Hn = H / np.linalg.norm(H, axis=1, keepdims=True)
                    sn = s / (np.linalg.norm(s) or 1.0)
                    sup[hname] = float((Hn @ sn).max())
            r["split_support"] = sup
            both = (sup["A"] is not None and sup["B"] is not None
                    and sup["A"] >= MATCH_COS and sup["B"] >= MATCH_COS)
            pats = r["patients"]
            grade = ("EXPLORATORY_REPRODUCED" if len(pats) >= 2
                     else "DESCRIPTIVE_SINGLE_PATIENT")
            pid = f"R16J{arm}{r['resolution']}-{r['cluster_id']:03d}"
            notes = (f"arm={arm};split_half={'yes' if both else 'no'}"
                     + (f"(A={sup['A']:.2f},B={sup['B']:.2f})" if sup["A"] is not None else "( halves_skipped)")
                     + f";res={r['resolution']}")
            all_rows.append({
                "pattern_id": pid, "tier": "2", "kind": "joint_pattern_group",
                "axis_name": "NA", "lineage": "HTAN_VANDERBILT_CRC",
                "n_sections": r["n_sections"], "n_patients": len(pats),
                "resolution_support": str(r["resolution"]),
                "top_markers": ";".join(r["markers"]), "geometry": "NA",
                "effect_size": (f"{r['moran_i_median']:.3f}" if r["moran_i_median"] is not None else "NA"),
                "null_type": "value_permutation", "null_draws": args.null_draws,
                "null_p": (f"{r['moran_p_median']:.4f}" if r["moran_p_median"] is not None else "NA"),
                "fdr_q": "NA", "cross_patient_shape_corr": "NA",
                "residual_dimension": "NA", "uncertainty_ci95": "NA",
                "evidence_grade": grade, "notes": notes,
            })
            all_cands.append({k: (v.tolist() if isinstance(v, np.ndarray) else v)
                              for k, v in r.items() if k != "sig_raw"})
    artifact = {
        "schema": "r16.joint_embed.v1",
        "status": "EXPLORATORY_COMPLETE_NOT_CLAIM",
        "seed": C.SEED,
        "parameters": {"arms": args.arms, "batch_key": "section",
                       "pca_dims": C.PCA_DIMS, "knn_expr": C.KNN_EXPR,
                       "resolutions": list(C.RESOLUTIONS),
                       "min_cluster_size": C.MIN_CLUSTER_SIZE,
                       "match_cosine": MATCH_COS, "null_draws": args.null_draws,
                       "augment_lambda": J.AUGMENT_LAMBDA, "augment_knn": J.KNN_NB,
                       "halves": half_info["halves"]},
        "half_match": {a: results[a].get("half_match") for a in args.arms},
        "n_candidates": len(all_cands),
        "candidates": all_cands,
    }
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, ensure_ascii=False))
    out_rows.write_text(json.dumps(all_rows, indent=2, ensure_ascii=False))
    n_rep = sum(1 for r in all_rows if r["evidence_grade"] == "EXPLORATORY_REPRODUCED")
    n_sup = sum(1 for r in all_cands if r["split_support"].get("A") is not None
                and r["split_support"]["A"] >= MATCH_COS and r["split_support"]["B"] >= MATCH_COS)
    print(f"candidates={len(all_cands)} reproduced(>=2pat)={n_rep} split_supported={n_sup} "
          f"({time.time()-t0:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
