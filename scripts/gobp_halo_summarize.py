#!/usr/bin/env python3
"""Summarise the v6 GOBP halo screen: ranking, conditioning loss, circularity,
sign-consistency against a binomial reference."""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/gobp_halo_20261001"
PANEL = ROOT / "infra/gobp_compress_20261001/survivors_specific_floor20_jac70.tsv"
GMT = Path("/home/huyudi/006/data/pathway/c5.go.v2025.1.Hs.symbols.gmt")
OBO = Path("/home/huyudi/012_conference/aaai2027/data/raw/ontology/go-basic.obo")


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    summ = list(csv.DictReader(open(OUT / f"gobp_halo_summary{args.tag}.tsv"), delimiter="\t"))
    secs = list(csv.DictReader(open(OUT / f"gobp_halo_sections{args.tag}.tsv"), delimiter="\t"))
    print(f"sets={len(summ)} section-rows={len(secs)}")

    def f(v):
        return float(v) if v not in ("", None) else np.nan

    pool = load_module("pool", ROOT / "scripts/tls_pool_expand.py")
    study = load_module("study", ROOT / "scripts/gobp_compress_study.py")
    gmt = {n: g for n, g in study.read_gmt(GMT).items() if n.startswith("GOBP_")}
    sig = set(pool.signature_symbols())

    # circularity + family label
    terms, parents = study.parse_obo(OBO)
    depth, anc = study.depths_and_ancestors(parents)
    by_norm = defaultdict(list)
    for tid, d in terms.items():
        by_norm[study.norm(d["name"])].append(tid)

    def fam(pid):
        cand = by_norm.get(study.norm(pid.replace("GOBP_", "", 1)))
        if not cand:
            return ""
        tid = sorted(cand)[0]
        a = [(depth.get(x, 99), x) for x in anc.get(tid, set())]
        top = sorted(a)[:3]
        return " < ".join(terms[x]["name"] for _, x in top[:2]) if top else "root"

    rows = []
    for s in summ:
        pid = s["program_id"]
        g = set(gmt.get(pid, []))
        ov = len(g & sig)
        rows.append({
            "program_id": pid,
            "n_sections": int(s["n_sections"]),
            "raw_pos": f(s["raw_pos_frac"]), "raw_med": f(s["raw_median"]),
            "own_pos": f(s["own_pos_frac"]), "own_med": f(s["own_median"]),
            "local_pos": f(s["local_pos_frac"]), "local_med": f(s["local_median"]),
            "k_used": f(s["median_k_used"]), "raw_size": f(s["raw_size"]),
            "covered": f(s["covered_size"]), "depth": f(s["depth"]),
            "sig_overlap": ov, "sig_frac": round(ov / max(1, len(g)), 3),
            "family": fam(pid),
        })
    arr = {k: np.array([r[k] for r in rows], dtype=float)
           for k in ("raw_pos", "raw_med", "own_pos", "own_med", "local_pos", "local_med")}
    ok = np.isfinite(arr["own_med"]) & np.isfinite(arr["raw_med"])
    print(f"evaluable sets (finite own+raw): {int(ok.sum())}")

    # conditioning loss
    loss = arr["raw_med"] - arr["own_med"]
    print(f"raw_med > 0: {np.nanmean(arr['raw_med'] > 0):.3f}   own_med > 0: {np.nanmean(arr['own_med'] > 0):.3f}")
    print(f"median raw_med {np.nanmedian(arr['raw_med']):+.5f}  median own_med {np.nanmedian(arr['own_med']):+.5f}")
    print(f"corr(raw_med, own_med) = {np.corrcoef(arr['raw_med'][ok], arr['own_med'][ok])[0,1]:.3f}")
    print(f"sets with raw_pos>=0.6: {np.nanmean(arr['raw_pos'] >= 0.6):.3f}   "
          f"own_pos>=0.6: {np.nanmean(arr['own_pos'] >= 0.6):.3f}")

    # binomial reference for sign consistency (53 sections)
    from scipy.stats import binom
    n = 53
    print("\nbinomial reference (n=53, p=0.5): expected sets by chance")
    for frac in (0.6, 0.7, 0.8):
        k = int(np.ceil(frac * n))
        p = float(binom.sf(k - 1, n, 0.5))
        print(f"  >= {frac:.0%} positive: p={p:.3g} -> {p * len(rows):.1f} of {len(rows)} sets")

    order = np.argsort(-np.nan_to_num(arr["own_med"], nan=-9))
    top = [rows[i] for i in order[:40]]
    with (OUT / f"gobp_halo_top{args.tag}.tsv").open("w", newline="") as fh:
        keys = list(rows[0])
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader()
        for r in top:
            w.writerow(r)
    print(f"\ntop 25 by own_med (own_pos/raw_pos, size, sig overlap, family):")
    for r in top[:25]:
        print(f"  {r['program_id'][:58]:58} own={r['own_med']:+.5f} ({r['own_pos']:.2f}) "
              f"raw={r['raw_med']:+.5f} ({r['raw_pos']:.2f}) n={int(r['raw_size']):4d} "
              f"sig={r['sig_overlap']:2d} d={int(r['depth']) if np.isfinite(r['depth']) else 0} "
              f"{r['family'][:38]}")
    circ = [r for r in rows if r["sig_frac"] >= 0.2]
    print(f"\nsets with >=20% detection-signature overlap: {len(circ)}")

    # ---- label-permutation multiplicity reference (preserves the per-section
    # marginal distribution of the statistic, destroys set identity) ----
    secs_list = sorted({r["section_id"] for r in secs})
    pids_list = [r["program_id"] for r in rows]
    idx_s = {s: i for i, s in enumerate(secs_list)}
    idx_p = {p: i for i, p in enumerate(pids_list)}
    Ms = {}
    for arm in ("own", "raw"):
        M = np.full((len(secs_list), len(pids_list)), np.nan)
        for r in secs:
            v = r[f"{arm}_delta"]
            if v not in ("", None):
                M[idx_s[r["section_id"]], idx_p[r["program_id"]]] = float(v)
        Ms[arm] = M

    def count_sets(M, thr, min_sec=8):
        f = np.isfinite(M)
        n = f.sum(axis=0)
        frac = ((M > 0) & f).sum(axis=0) / np.maximum(n, 1)
        return int(((frac >= thr) & (n >= min_sec)).sum())

    rng = np.random.default_rng(20261001)
    mult_rows = []
    for arm, M in Ms.items():
        for thr in (0.6, 0.65, 0.7):
            obs = count_sets(M, thr)
            null = []
            for _ in range(200):
                Mp = np.full_like(M, np.nan)
                for i in range(M.shape[0]):
                    idx = np.where(np.isfinite(M[i]))[0]
                    v = M[i][idx].copy()
                    rng.shuffle(v)
                    Mp[i, idx] = v
                null.append(count_sets(Mp, thr))
            null = np.array(null)
            mult_rows.append({"arm": arm, "threshold": thr, "observed": obs,
                              "null_median": int(np.median(null)),
                              "null_q95": int(np.quantile(null, 0.95)),
                              "null_max": int(null.max()),
                              "ratio_vs_null_median": round(obs / max(1, np.median(null)), 2)})
            print(f"  multiplicity [{arm}] >= {thr:.2f}: observed {obs} vs permutation null "
                  f"median {np.median(null):.0f} (q95 {np.quantile(null,0.95):.0f})")
    with (OUT / f"gobp_halo_multiplicity{args.tag}.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(mult_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(mult_rows)


if __name__ == "__main__":
    main()
