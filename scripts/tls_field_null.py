#!/usr/bin/env python3
"""v5 halo calibration: edge distance + selection/autocorrelation null.

Contract: infra/tls_field_20260930/run_contract_v5.json (frozen before running).
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/tls_field_20260930"
V5 = OUT / "run_contract_v5.json"
BIN_EDGES = np.arange(0, 1001, 100)
NBINS = len(BIN_EDGES) - 1
NEAR_I = 2   # 200-300
FAR_I = 7    # 700-800
MIN_BIN_SPOTS = 5
NEAR_FLOOR = 0.005
MIN_ANCHORS = 8
MIN_DRAW_ANCHORS = 8
LAG_BANDS = [(100, 200), (200, 400), (400, 700), (700, 1000)]
MAX_PAIRS = 200_000


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_gmt(path):
    sets = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3:
            sets[p[0]] = [g for g in p[2:] if g]
    return sets


# ---------------------------------------------------------------- geometry
def estimate_sigma(xy, values, rng):
    """Fit rho(d) = exp(-d^2/(4 sigma^2)) on lag-band correlations."""
    v = np.asarray(values, float)
    v = (v - np.nanmean(v)) / (np.nanstd(v) + 1e-12)
    pairs = cKDTree(xy).query_pairs(1000.0, output_type="ndarray")
    if len(pairs) > MAX_PAIRS:
        pairs = pairs[rng.choice(len(pairs), MAX_PAIRS, replace=False)]
    if len(pairs) < 100:
        return float("nan"), {}
    d = np.linalg.norm(xy[pairs[:, 0]] - xy[pairs[:, 1]], axis=1)
    prod = v[pairs[:, 0]] * v[pairs[:, 1]]
    xs, ys, cal = [], [], []
    for lo, hi in LAG_BANDS:
        sel = (d >= lo) & (d < hi)
        if sel.sum() < 50:
            continue
        rho = float(np.nanmean(prod[sel]))
        cal.append((f"{lo}-{hi}", round(rho, 4)))
        if rho > 0.05:
            xs.append(((lo + hi) / 2.0) ** 2 / 4.0)
            ys.append(-np.log(rho))
    xs, ys = np.array(xs), np.array(ys)
    if len(xs) < 2 or (xs ** 2).sum() <= 0:
        return float("nan"), dict(cal)
    a = float((xs * ys).sum() / (xs ** 2).sum())   # a = 1/(4 sigma^2)
    sigma = float(np.sqrt(1.0 / (4.0 * a))) if a > 0 else float("nan")
    sigma = float(np.clip(sigma, 40.0, 800.0))
    return sigma, dict(cal)


def smooth_op(xy, sigma, mult=3.0):
    tree = cKDTree(xy)
    pairs = tree.query_pairs(mult * sigma, output_type="ndarray")
    if len(pairs) == 0:
        return sparse.identity(len(xy), format="csr")
    i = np.concatenate([pairs[:, 0], pairs[:, 1], np.arange(len(xy))])
    j = np.concatenate([pairs[:, 1], pairs[:, 0], np.arange(len(xy))])
    d = np.linalg.norm(xy[i] - xy[j], axis=1)
    w = np.exp(-(d ** 2) / (2.0 * sigma ** 2))
    W = sparse.csr_matrix((w, (i, j)), shape=(len(xy), len(xy)))
    rs = np.asarray(W.sum(axis=1)).ravel()
    rs[rs == 0] = 1.0
    return sparse.diags(1.0 / rs) @ W


def quantile_map(z, observed):
    """Monotone map of each column of z to the empirical quantiles of observed."""
    obs = np.sort(np.asarray(observed, float))
    ranks = np.argsort(np.argsort(z, axis=0), axis=0)
    return obs[ranks]


def per_anchor_bin_means(xy, comps, pitch_um, values):
    """Voronoi edge-distance bin means per component. Returns means, counts (ncomp x NBINS)."""
    ncomp = len(comps)
    means = np.full((ncomp, NBINS), np.nan)
    counts = np.zeros((ncomp, NBINS), dtype=int)
    if ncomp == 0:
        return means, counts
    anchor_idx = np.concatenate([np.asarray(c) for c in comps])
    comp_of_anchor = np.repeat(np.arange(ncomp), [len(c) for c in comps])
    d, idx = cKDTree(xy[anchor_idx]).query(xy)
    cid = comp_of_anchor[idx]
    bin_id = np.digitize(d, BIN_EDGES) - 1
    vals = np.asarray(values, float)
    valid = (bin_id >= 0) & (bin_id < NBINS) & (d >= 1.3 * pitch_um) & np.isfinite(vals)
    if not valid.any():
        return means, counts
    key = cid[valid] * NBINS + bin_id[valid]
    counts = np.bincount(key, minlength=ncomp * NBINS).reshape(ncomp, NBINS)
    sums = np.bincount(key, weights=vals[valid], minlength=ncomp * NBINS).reshape(ncomp, NBINS)
    with np.errstate(invalid="ignore", divide="ignore"):
        means = sums / np.where(counts > 0, counts, 1)
    means[counts == 0] = np.nan
    return means, counts


def section_statistic(means, counts):
    if len(means) == 0:
        return float("nan"), 0
    ok = ((counts[:, NEAR_I] >= MIN_BIN_SPOTS) & (counts[:, FAR_I] >= MIN_BIN_SPOTS)
          & np.isfinite(means[:, NEAR_I]) & np.isfinite(means[:, FAR_I])
          & (means[:, NEAR_I] >= NEAR_FLOOR))
    if ok.sum() < MIN_ANCHORS:
        return float("nan"), int(ok.sum())
    rel = 1.0 - means[ok, FAR_I] / means[ok, NEAR_I]
    return float(np.median(rel)), int(ok.sum())


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20261001)
    ap.add_argument("--sections", default="", help="comma list to subset (smoke)")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    v5 = json.loads(V5.read_text())
    assert v5["status"] == "FROZEN_BEFORE_SCORING"
    primary = v5["operator"]["primary_readouts"]
    pool = load_module("pool", ROOT / "scripts/tls_pool_expand.py")
    hallmark = load_module("hallmark", ROOT / "scripts/tls_field_hallmark.py")
    gmt = read_gmt(ROOT / "data/geneset/h.all.v2026.1.Hs.symbols.gmt")
    proxy = json.load(open(ROOT / "infra/r04/marker_proxy_combined.json"))
    b_genes = sorted(set(proxy["classes"]["B"]["voted_genes"]))
    readouts = {k: v for k, v in gmt.items() if k in primary}
    readouts["B_AXIS"] = b_genes

    diag = list(csv.DictReader(open(ROOT / "infra/tls_pool_expand_20260930/pool_diagnostics.tsv"),
                               delimiter="\t"))
    ok = [r for r in diag if r["status"] == "ok" and r["cohort"] != "Xenium"]
    if args.sections:
        keep = set(args.sections.split(","))
        ok = [r for r in ok if r["section_id"] in keep]

    curve_rows, null_rows = [], []
    for row in ok:
        sid, cohort = row["section_id"], row["cohort"]
        loaded = hallmark.load_section_by_id(pool, cohort, sid)
        if loaded is None or loaded[0] is None:
            print(f"{sid}: load failed", flush=True)
            continue
        _, X, genes, xy, pitch, labels, index = loaded
        gidx = dict(index)
        cols_of = {}
        for tag, glist in list(readouts.items()) + [("__sig__", pool.signature_symbols())]:
            cols = []
            for g in glist:
                if g in gidx:
                    v = gidx[g]
                    cols.extend(v if isinstance(v, list) else [v])
            cols_of[tag] = cols
        scores = {t: pool.aucell(X, c, 0.05) for t, c in cols_of.items()}
        s = scores.pop("__sig__")
        if not np.isfinite(s).all():
            print(f"{sid}: detection score has non-finite values, skipped", flush=True)
            continue
        finite = np.isfinite(s)
        thr = float(np.quantile(s[finite], 0.9))
        mask = finite & (s >= thr)
        comps = pool.components(xy, mask, pitch)
        sigma_s, cal_s = estimate_sigma(xy, s, rng)
        # observed edge curves + statistics
        obs = {}
        for tag, v in scores.items():
            if not np.isfinite(v).all():
                print(f"{sid}: {tag} has non-finite scores, dropped", flush=True)
                continue
            means, counts = per_anchor_bin_means(xy, comps, pitch, v)
            stat, n_eff = section_statistic(means, counts)
            obs[tag] = (stat, n_eff)
            med = np.nanmedian(means, axis=0)
            for bi in range(NBINS):
                curve_rows.append({"section_id": sid, "cohort": cohort, "set": tag,
                                   "bin_lo": int(BIN_EDGES[bi]), "median": round(float(med[bi]), 5),
                                   "n_anchors": int(np.isfinite(means[:, bi]).sum())})
        scores = {t: v for t, v in scores.items() if np.isfinite(v).all()}
        # surrogate null
        if not np.isfinite(sigma_s):
            print(f"{sid}: sigma_s unfit, null skipped", flush=True)
            for tag, (stat, n_eff) in obs.items():
                null_rows.append({"section_id": sid, "cohort": cohort, "set": tag,
                                  "obs": stat, "n_effective": n_eff, "null_median": "",
                                  "p_one_sided": "", "draws_used": 0,
                                  "sigma_s": "", "sigma_r": "", "rho_s": "", "note": "sigma_s unfit"})
            continue
        Ws = smooth_op(xy, sigma_s)
        N = args.draws
        Z1 = rng.standard_normal((len(xy), N))
        S_sim = Ws @ Z1
        S_sim = (S_sim - S_sim.mean(axis=0)) / (S_sim.std(axis=0) + 1e-12)
        for tag, r in scores.items():
            sigma_r, cal_r = estimate_sigma(xy, r, rng)
            if not np.isfinite(sigma_r):
                null_rows.append({"section_id": sid, "cohort": cohort, "set": tag,
                                  "obs": obs[tag][0], "n_effective": obs[tag][1],
                                  "null_median": "", "p_one_sided": "", "draws_used": 0,
                                  "sigma_s": round(sigma_s, 1), "sigma_r": "", "rho_s": "",
                                  "note": "sigma_r unfit"})
                continue
            rho = float(np.corrcoef(np.argsort(np.argsort(s)), np.argsort(np.argsort(r)))[0, 1])
            Wr = smooth_op(xy, sigma_r)
            Z2 = rng.standard_normal((len(xy), N))
            R_sim = Wr @ Z2
            R_sim = (R_sim - R_sim.mean(axis=0)) / (R_sim.std(axis=0) + 1e-12)
            mix = rho * S_sim + np.sqrt(max(0.0, 1.0 - rho ** 2)) * R_sim
            S_q = quantile_map(S_sim, s)
            R_q = quantile_map(mix, r)
            null_stats = []
            comp_counts = []
            for j in range(N):
                m = S_q[:, j] >= float(np.quantile(S_q[:, j], 0.9))
                cj = pool.components(xy, m, pitch)
                if len(cj) == 0:
                    continue
                comp_counts.append(len(cj))
                mm, cc = per_anchor_bin_means(xy, cj, pitch, R_q[:, j])
                st, ne = section_statistic(mm, cc)
                if np.isfinite(st) and ne >= MIN_DRAW_ANCHORS:
                    null_stats.append(st)
            if len(null_stats) < max(20, int(0.5 * N)):
                null_rows.append({"section_id": sid, "cohort": cohort, "set": tag,
                                  "obs": obs[tag][0], "n_effective": obs[tag][1],
                                  "null_median": "", "p_one_sided": "", "draws_used": len(null_stats),
                                  "sigma_s": round(sigma_s, 1), "sigma_r": round(sigma_r, 1),
                                  "rho_s": round(rho, 4), "note": "too few valid draws"})
                continue
            ns = np.array(null_stats)
            stat = obs[tag][0]
            p = float((1 + np.sum(ns >= stat)) / (len(ns) + 1)) if np.isfinite(stat) else float("nan")
            null_rows.append({"section_id": sid, "cohort": cohort, "set": tag,
                              "obs": stat, "n_effective": obs[tag][1],
                              "null_median": round(float(np.median(ns)), 5),
                              "p_one_sided": round(p, 4) if p == p else "",
                              "draws_used": len(ns), "sigma_s": round(sigma_s, 1),
                              "sigma_r": round(sigma_r, 1), "rho_s": round(rho, 4),
                              "note": f"ncomp_obs={len(comps)};ncomp_null_med={int(np.median(comp_counts)) if comp_counts else -1};cal_s={cal_s}"})
        print(f"{sid}: anchors={len(comps)} sigma_s={sigma_s:.0f} cal={cal_s}", flush=True)

    with (OUT / "field_edge_curves_v5.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["section_id", "cohort", "set", "bin_lo", "median", "n_anchors"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(curve_rows)
    with (OUT / "field_null_v5.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["section_id", "cohort", "set", "obs", "n_effective",
                                          "null_median", "p_one_sided", "draws_used",
                                          "sigma_s", "sigma_r", "rho_s", "note"], delimiter="\t")
        w.writeheader()
        w.writerows(null_rows)
    print(f"wrote curves={len(curve_rows)} null={len(null_rows)}")


if __name__ == "__main__":
    main()
