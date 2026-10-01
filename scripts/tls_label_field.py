#!/usr/bin/env python3
"""v8.1 field test with independent (pathology) TLS anchors.

Primary statistic: SECTION-level edge-distance profile (whole TLS mask), with a
matched-autocorrelation smooth-surrogate null using the same mask.
Secondary: per-anchor profile for sections with enough evaluable anchors (descriptive).

Contract: infra/tls_label_field_20261001/run_contract_v8.json (v8.1, frozen pre-scoring).
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/tls_label_field_20261001"
CONTRACT = OUT / "run_contract_v8.json"
V7_TARGETS = ROOT / "infra/gobp_halo_20261001/target_sets_v7.tsv"
HALLMARK = ROOT / "data/geneset/h.all.v2026.1.Hs.symbols.gmt"
SINGLE_GENES = ["CXCL13", "CCL19", "CCL21", "LTB", "CD74", "HLA-DRA", "IFNG",
                "CXCL9", "GZMB", "CD8A", "MKI67", "COL1A1"]
BIN_EDGES = np.arange(0, 1001, 100)
NBINS = len(BIN_EDGES) - 1
NEAR_I, FAR_I = 2, 7
MIN_BIN_SPOTS = 30
MIN_SECTIONS = 8
RING_FACTOR = 1.3


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_gmt(path):
    out = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3:
            out[p[0]] = [g for g in p[2:] if g]
    return out


def distance_bins(xy, tls_mask, pitch):
    """Distance to nearest TLS spot, ring exclusion, bin incidence and counts."""
    d = cKDTree(xy[tls_mask]).query(xy)[0]
    valid = d >= RING_FACTOR * pitch
    bin_id = np.digitize(d, BIN_EDGES) - 1
    keep = valid & (bin_id >= 0) & (bin_id < NBINS)
    rows = np.where(keep)[0]
    keys = bin_id[keep]
    B = sparse.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, keys)),
                          shape=(len(xy), NBINS))
    cnt = np.asarray(B.sum(axis=0)).ravel()
    return B, cnt, d, keep


def smooth_gaussian(xy, values, sigma, mult=3.0):
    tree = cKDTree(xy)
    pr = tree.query_pairs(mult * sigma, output_type="ndarray")
    n = len(xy)
    rows = np.concatenate([pr[:, 0], pr[:, 1], np.arange(n)])
    cols = np.concatenate([pr[:, 1], pr[:, 0], np.arange(n)])
    dd = np.linalg.norm(xy[rows] - xy[cols], axis=1)
    w = np.exp(-(dd ** 2) / (2 * sigma ** 2))
    W = sparse.csr_matrix((w, (rows, cols)), shape=(n, n))
    rs = np.asarray(W.sum(axis=1)).ravel()
    rs[rs == 0] = 1.0
    return np.asarray((sparse.diags(1.0 / rs) @ W) @ values)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--sections", default="")
    ap.add_argument("--tag", default="")
    ap.add_argument("--limit-readouts", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(20261001)

    v8 = json.loads(CONTRACT.read_text())
    assert v8["status"] == "FROZEN_BEFORE_SCORING" and str(v8["version"]).startswith("8.1")
    screen = load_module("screen", ROOT / "scripts/gobp_halo_screen.py")
    nullmod = load_module("nullmod", ROOT / "scripts/gobp_halo_null.py")
    pool = load_module("pool", ROOT / "scripts/tls_pool_expand.py")
    hallmark = load_module("hallmark", ROOT / "scripts/tls_field_hallmark.py")
    assert screen.NBINS == NBINS and (screen.NEAR_I, screen.FAR_I) == (NEAR_I, FAR_I)

    targets = list(csv.DictReader(open(V7_TARGETS), delimiter="\t"))
    readouts = [(t["set_id"], t["genes"].split(",")) for t in targets]
    hm = read_gmt(HALLMARK)
    for key in ("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                "HALLMARK_G2M_CHECKPOINT", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"):
        readouts.append((key, hm[key]))
    for g in SINGLE_GENES:
        readouts.append((f"GENE_{g.replace('-', '_')}", [g]))
    if args.limit_readouts:
        readouts = readouts[:args.limit_readouts]

    diag = list(csv.DictReader(open(ROOT / "infra/tls_pool_expand_20260930/pool_diagnostics.tsv"),
                               delimiter="\t"))
    ok = [r for r in diag if r["status"] == "ok" and r["cohort"] in ("USZ", "GSE175540")]
    if args.sections:
        keep = set(args.sections.split(","))
        ok = [r for r in ok if r["section_id"] in keep]

    rows, null_rows, secondary = [], [], []
    t_start = time.time()
    for row in ok:
        sid, cohort = row["section_id"], row["cohort"]
        t0 = time.time()
        loaded = hallmark.load_section_by_id(pool, cohort, sid)
        if loaded is None or loaded[0] is None:
            continue
        _, X, genes, xy, pitch, labels, index = loaded
        y = np.asarray(labels.get("y"), dtype=object)
        tls_mask = (y == "TLS")
        if tls_mask.sum() < 20:
            print(f"{sid}: too few TLS spots", flush=True)
            continue
        gidx = dict(index)

        def member(gl):
            idxs = []
            for g in gl:
                v = gidx.get(g)
                if v is None:
                    continue
                idxs.extend(v if isinstance(v, list) else [v])
            return np.array(sorted(idxs), dtype=np.int32)

        members = [member(gl) for _, gl in readouts] + [member(pool.signature_symbols())]
        S, k = screen.score_sets(X, members)
        s = S[:, -1].astype(np.float64)
        if not np.isfinite(s).all():
            continue
        B, cnt, dist, keep = distance_bins(xy, tls_mask, pitch)
        if min(cnt[NEAR_I], cnt[FAR_I]) < MIN_BIN_SPOTS:
            print(f"{sid}: bins too small (near {int(cnt[NEAR_I])}, far {int(cnt[FAR_I])})",
                  flush=True)
            continue
        Bn = B @ sparse.diags(1.0 / np.maximum(cnt, 1))
        Q_own = screen.residual_operator(screen.design(s, xy))
        s_smooth = smooth_gaussian(xy, s, 300.0)
        Q_loc = screen.residual_operator(screen.design(s_smooth, xy))
        # secondary per-anchor profile (descriptive)
        comps = pool.components(xy, tls_mask, pitch)
        n_anchor_ok = 0
        if len(comps):
            Ca, cnta = screen.anchor_operators(xy, comps, pitch)
            if Ca.shape[1]:
                cma = cnta.reshape(-1, NBINS)
                n_anchor_ok = int(((cma[:, NEAR_I] >= 5) & (cma[:, FAR_I] >= 5)).sum())
        pairs, pd = nullmod.pair_list(xy, rng)
        cache = {}
        for j, (set_id, _) in enumerate(readouts):
            if len(members[j]) == 0:
                continue
            r = S[:, j].astype(np.float64)
            if not np.isfinite(r).all():
                continue
            R_own = screen.residualize(r, Q_own)
            R_loc = screen.residualize(r, Q_loc)
            def contrast_1d(vec):
                m = np.asarray(Bn.T @ vec)
                return float(m[NEAR_I] - m[FAR_I])
            obs_raw = contrast_1d(r)
            obs_own = contrast_1d(R_own)
            obs_loc = contrast_1d(R_loc)
            rows.append({"section_id": sid, "cohort": cohort, "set_id": set_id,
                         "raw_contrast": round(obs_raw, 6), "own_contrast": round(obs_own, 6),
                         "local_contrast": round(obs_loc, 6), "n_near": int(cnt[NEAR_I]),
                         "n_far": int(cnt[FAR_I]), "k_used": int(k[j])})
            if n_anchor_ok >= 8:
                secondary.append({"section_id": sid, "set_id": set_id, "set": set_id,
                                  "n_anchor_ok": n_anchor_ok})
            corr = nullmod.band_corr(R_own, pairs, pd) if pairs is not None else {}
            sigma = nullmod.sigma_from_bands(corr)
            if not np.isfinite(sigma):
                # v8.2 fallback: the conditioned residual shows no detectable
                # spatial autocorrelation, so the matched surrogate is
                # uncorrelated noise with the same marginal distribution.
                S_sim = nullmod.quantile_map(rng.standard_normal((len(xy), args.draws)), R_own)
                means = np.asarray(Bn.T @ S_sim)
                null_con = means[NEAR_I, :] - means[FAR_I, :]
                p = float((1 + np.sum(null_con >= obs_own)) / (args.draws + 1))
                null_rows.append({"section_id": sid, "cohort": cohort, "set_id": set_id,
                                  "obs": round(obs_own, 6),
                                  "null_median": round(float(np.median(null_con)), 6),
                                  "p_one_sided": round(p, 4), "sigma": "",
                                  "note": "white_fallback"})
                continue
            key = round(sigma, 1)
            if key not in cache:
                cache[key] = nullmod.smoothed_noise(xy, sigma, args.draws, rng)
            S_sim = nullmod.quantile_map(cache[key], R_own)
            means = np.asarray(Bn.T @ S_sim)          # (NBINS, draws)
            null_con = means[NEAR_I, :] - means[FAR_I, :]
            p = float((1 + np.sum(null_con >= obs_own)) / (args.draws + 1))
            null_rows.append({"section_id": sid, "cohort": cohort, "set_id": set_id,
                              "obs": round(obs_own, 6),
                              "null_median": round(float(np.median(null_con)), 6),
                              "p_one_sided": round(p, 4), "sigma": round(sigma, 1),
                              "note": ""})
        print(f"{sid}: tls={int(tls_mask.sum())} near={int(cnt[NEAR_I])} far={int(cnt[FAR_I])} "
              f"anchor_ok={n_anchor_ok} elapsed={time.time()-t0:.1f}s", flush=True)

    with (OUT / f"label_field_sections{args.tag}.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    with (OUT / f"label_field_null{args.tag}.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(null_rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(null_rows)
    if secondary:
        with (OUT / f"label_field_secondary{args.tag}.tsv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(secondary[0]), delimiter="\t")
            w.writeheader(); w.writerows(secondary)

    out = []
    for set_id, _ in readouts:
        sub = [r for r in null_rows if r["set_id"] == set_id]
        ps = np.array([float(r["p_one_sided"]) for r in sub if r["p_one_sided"] != ""])
        obs = np.array([float(r["obs"]) for r in sub if r["obs"] != ""])
        nl = np.array([float(r["null_median"]) for r in sub if r["null_median"] != ""])
        rsub = [r for r in rows if r["set_id"] == set_id]
        raw = np.array([float(r["raw_contrast"]) for r in rsub])
        own = np.array([float(r["own_contrast"]) for r in rsub])
        frac = float((ps <= 0.05).mean()) if len(ps) else float("nan")
        med_obs = float(np.median(obs)) if len(obs) else float("nan")
        med_null = float(np.median(nl)) if len(nl) else float("nan")
        out.append({"set_id": set_id, "n_sections": len(rsub), "n_evaluable": len(ps),
                    "raw_pos_frac": round(float((raw > 0).mean()), 3) if len(raw) else "",
                    "own_pos_frac": round(float((own > 0).mean()), 3) if len(own) else "",
                    "raw_median": round(float(np.median(raw)), 6) if len(raw) else "",
                    "own_median": round(float(np.median(own)), 6) if len(own) else "",
                    "null_frac_p05": round(frac, 3) if frac == frac else "",
                    "median_null": round(med_null, 6) if med_null == med_null else "",
                    "median_p": round(float(np.median(ps)), 4) if len(ps) else "",
                    "claim": int(len(ps) >= MIN_SECTIONS and frac >= 0.5 and med_obs > med_null)})
    out.sort(key=lambda r: (-(r["null_frac_p05"] if isinstance(r["null_frac_p05"], float) else -1),
                            -(r["own_median"] if isinstance(r["own_median"], float) else -9)))
    with (OUT / f"label_field_summary{args.tag}.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader(); w.writerows(out)
    print(f"wrote sections={len(rows)} null={len(null_rows)} elapsed={time.time()-t_start:.0f}s")
    for r in out[:12]:
        print(f"  {r['set_id'][:50]:50} n={r['n_evaluable']:2} own_pos={r['own_pos_frac']} "
              f"own_med={r['own_median']} frac_p05={r['null_frac_p05']} claim={r['claim']}")


if __name__ == "__main__":
    main()
