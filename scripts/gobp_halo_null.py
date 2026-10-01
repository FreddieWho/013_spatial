#!/usr/bin/env python3
"""v7 halo null for the top-30 GOBP sets + 3 controls (contract v7).

Fixed-anchor matched-autocorrelation null on the own-conditioned residual.
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
OUT = ROOT / "infra/gobp_halo_20261001"
CONTRACT = OUT / "run_contract_v7.json"
TARGETS = OUT / "target_sets_v7.tsv"
GMT = Path("/home/huyudi/006/data/pathway/c5.go.v2025.1.Hs.symbols.gmt")
BIN_EDGES = np.arange(0, 1001, 100)
NBINS = len(BIN_EDGES) - 1
NEAR_I, FAR_I = 2, 7
MIN_BIN_SPOTS = 5
MIN_ANCHORS = 8
LAG_BANDS = [(100, 200), (200, 400), (400, 700), (700, 1000)]
MAX_PAIRS = 200_000


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pair_list(xy, rng, max_pairs=MAX_PAIRS):
    pairs = cKDTree(xy).query_pairs(1000.0, output_type="ndarray")
    if len(pairs) > max_pairs:
        pairs = pairs[rng.choice(len(pairs), max_pairs, replace=False)]
    if len(pairs) < 100:
        return None, None
    d = np.linalg.norm(xy[pairs[:, 0]] - xy[pairs[:, 1]], axis=1)
    return pairs, d


def band_corr(v, pairs, d):
    sd = v.std()
    if not np.isfinite(sd) or sd <= 0:
        return {}
    z = (v - v.mean()) / sd
    prod = z[pairs[:, 0]] * z[pairs[:, 1]]
    out = {}
    for lo, hi in LAG_BANDS:
        sel = (d >= lo) & (d < hi)
        if sel.sum() >= 50:
            out[(lo, hi)] = float(np.nanmean(prod[sel]))
    return out


def sigma_from_bands(corr):
    xs, ys = [], []
    for (lo, hi), rho in corr.items():
        if rho > 0.05:
            xs.append(((lo + hi) / 2.0) ** 2 / 4.0)
            ys.append(-np.log(rho))
    if len(xs) < 2:
        return float("nan")
    xs, ys = np.array(xs), np.array(ys)
    a = float((xs * ys).sum() / (xs ** 2).sum())
    if a <= 0:
        return float("nan")
    return float(np.clip(np.sqrt(1.0 / (4.0 * a)), 40.0, 800.0))


def smoothed_noise(xy, sigma, n_draws, rng, mult=3.0):
    pairs = cKDTree(xy).query_pairs(mult * sigma, output_type="ndarray")
    n = len(xy)
    rows = np.concatenate([pairs[:, 0], pairs[:, 1], np.arange(n)])
    cols = np.concatenate([pairs[:, 1], pairs[:, 0], np.arange(n)])
    dd = np.linalg.norm(xy[rows] - xy[cols], axis=1)
    w = np.exp(-(dd ** 2) / (2.0 * sigma ** 2))
    W = sparse.csr_matrix((w, (rows, cols)), shape=(n, n))
    rs = np.asarray(W.sum(axis=1)).ravel()
    rs[rs == 0] = 1.0
    W = sparse.diags(1.0 / rs) @ W
    Z = rng.standard_normal((n, n_draws))
    S = W @ Z
    S = (S - S.mean(axis=0)) / (S.std(axis=0) + 1e-12)
    return S


def quantile_map(z, observed):
    obs = np.sort(np.asarray(observed, float))
    ranks = np.argsort(np.argsort(z, axis=0), axis=0)
    return obs[ranks]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--sections", default="")
    ap.add_argument("--tag", default="")
    ap.add_argument("--limit-sets", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(20261001)

    v7 = json.loads(CONTRACT.read_text())
    assert v7["status"] == "FROZEN_BEFORE_SCORING"
    screen = load_module("screen", ROOT / "scripts/gobp_halo_screen.py")
    pool = load_module("pool", ROOT / "scripts/tls_pool_expand.py")
    hallmark = load_module("hallmark", ROOT / "scripts/tls_field_hallmark.py")

    targets = list(csv.DictReader(open(TARGETS), delimiter="\t"))
    if args.limit_sets:
        targets = targets[:args.limit_sets]
    set_ids = [t["set_id"] for t in targets]
    set_genes = [t["genes"].split(",") for t in targets]

    diag = list(csv.DictReader(open(ROOT / "infra/tls_pool_expand_20260930/pool_diagnostics.tsv"),
                               delimiter="\t"))
    ok = [r for r in diag if r["status"] == "ok" and r["cohort"] != "Xenium"]
    if args.sections:
        keep = set(args.sections.split(","))
        ok = [r for r in ok if r["section_id"] in keep]

    rows = []
    t_start = time.time()
    for row in ok:
        sid, cohort = row["section_id"], row["cohort"]
        t0 = time.time()
        loaded = hallmark.load_section_by_id(pool, cohort, sid)
        if loaded is None or loaded[0] is None:
            continue
        _, X, genes, xy, pitch, labels, index = loaded
        gidx = dict(index)

        def member(glist):
            idxs = []
            for g in glist:
                v = gidx.get(g)
                if v is None:
                    continue
                idxs.extend(v if isinstance(v, list) else [v])
            return np.array(sorted(idxs), dtype=np.int32)

        members = [member(g) for g in set_genes] + [member(pool.signature_symbols())]
        S, k = screen.score_sets(X, members)
        s = S[:, -1].astype(np.float64)
        if not np.isfinite(s).all():
            continue
        thr = float(np.quantile(s, 0.9))
        comps = pool.components(xy, s >= thr, pitch)
        if len(comps) == 0:
            continue
        C, cnt = screen.anchor_operators(xy, comps, pitch)
        if C.shape[1] == 0:
            continue
        A = C @ sparse.diags(1.0 / np.maximum(cnt, 1))
        cnt_m = cnt.reshape(-1, NBINS)
        pairs, d = pair_list(xy, rng)
        Q = screen.residual_operator(screen.design(s, xy))
        W_cache = {}
        for j, sid_set in enumerate(set_ids):
            if len(members[j]) == 0:
                continue
            r = S[:, j].astype(np.float64)
            if not np.isfinite(r).all():
                continue
            R = screen.residualize(r, Q)
            obs_arr, n_eff = screen.anchor_contrasts(np.asarray(A.T @ R), cnt)
            obs = float(np.ravel(obs_arr)[0])
            if not np.isfinite(obs):
                rows.append({"section_id": sid, "cohort": cohort, "set_id": sid_set,
                             "obs": "", "null_median": "", "p_one_sided": "", "n_anchors": n_eff,
                             "sigma": "", "rho_s": "", "note": "not evaluable"})
                continue
            corr = band_corr(R, pairs, d) if pairs is not None else {}
            sigma = sigma_from_bands(corr)
            rho = float(np.corrcoef(np.argsort(np.argsort(s)), np.argsort(np.argsort(r)))[0, 1])
            if not np.isfinite(sigma):
                rows.append({"section_id": sid, "cohort": cohort, "set_id": sid_set,
                             "obs": round(obs, 6), "null_median": "", "p_one_sided": "",
                             "n_anchors": n_eff, "sigma": "", "rho_s": round(rho, 4),
                             "note": "sigma unfit"})
                continue
            key = round(sigma, 1)
            if key not in W_cache:
                W_cache[key] = smoothed_noise(xy, sigma, args.draws, rng)
            S_sim = quantile_map(W_cache[key], R)
            means = np.asarray(A.T @ S_sim)                     # (K, draws)
            ncomp = len(cnt) // NBINS
            M = means.reshape(ncomp, NBINS, args.draws)
            ok_anchor = ((cnt_m[:, NEAR_I] >= MIN_BIN_SPOTS) & (cnt_m[:, FAR_I] >= MIN_BIN_SPOTS))
            if ok_anchor.sum() < MIN_ANCHORS:
                rows.append({"section_id": sid, "cohort": cohort, "set_id": sid_set,
                             "obs": round(obs, 6), "null_median": "", "p_one_sided": "",
                             "n_anchors": int(ok_anchor.sum()), "sigma": round(sigma, 1),
                             "rho_s": round(rho, 4), "note": "too few anchors for null"})
                continue
            deltas = M[ok_anchor, NEAR_I, :] - M[ok_anchor, FAR_I, :]
            null_med = np.median(deltas, axis=0)
            p = float((1 + np.sum(null_med >= obs)) / (args.draws + 1))
            rows.append({"section_id": sid, "cohort": cohort, "set_id": sid_set,
                         "obs": round(obs, 6),
                         "null_median": round(float(np.median(null_med)), 6),
                         "p_one_sided": round(p, 4), "n_anchors": int(ok_anchor.sum()),
                         "sigma": round(sigma, 1), "rho_s": round(rho, 4), "note": ""})
        print(f"{sid}: anchors={len(comps)} sets={len(set_ids)} "
              f"elapsed={time.time()-t0:.1f}s", flush=True)

    with (OUT / f"gobp_halo_null_v7{args.tag}.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["section_id", "cohort", "set_id", "obs", "null_median",
                                          "p_one_sided", "n_anchors", "sigma", "rho_s", "note"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    # per-set summary under the pre-declared claim rule
    out = []
    for sid_set in set_ids:
        sub = [r for r in rows if r["set_id"] == sid_set]
        ps = np.array([float(r["p_one_sided"]) for r in sub if r["p_one_sided"] != ""])
        obs = np.array([float(r["obs"]) for r in sub if r["obs"] != ""])
        nulls = np.array([float(r["null_median"]) for r in sub if r["null_median"] != ""])
        frac = float((ps <= 0.05).mean()) if len(ps) else float("nan")
        med_obs = float(np.median(obs)) if len(obs) else float("nan")
        med_null = float(np.median(nulls)) if len(nulls) else float("nan")
        claim = int(len(ps) >= 8 and frac >= 0.5 and med_obs > med_null)
        out.append({"set_id": sid_set, "n_evaluable": len(ps), "frac_p_le_0.05": round(frac, 3)
                    if frac == frac else "", "median_obs": round(med_obs, 6) if med_obs == med_obs else "",
                    "median_null": round(med_null, 6) if med_null == med_null else "",
                    "median_p": round(float(np.median(ps)), 4) if len(ps) else "",
                    "min_p": round(float(ps.min()), 4) if len(ps) else "",
                    "claim": claim})
    out.sort(key=lambda r: (-(r["frac_p_le_0.05"] if isinstance(r["frac_p_le_0.05"], float) else -1),
                            -(r["median_obs"] if isinstance(r["median_obs"], float) else -9)))
    with (OUT / f"gobp_halo_null_summary_v7{args.tag}.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    print(f"wrote rows={len(rows)} summary={len(out)} elapsed={time.time()-t_start:.0f}s")
    for r in out[:12]:
        print(f"  {r['set_id'][:52]:52} n={r['n_evaluable']:3} frac={r['frac_p_le_0.05']} "
              f"med_obs={r['median_obs']} med_null={r['median_null']} claim={r['claim']}")


if __name__ == "__main__":
    main()
