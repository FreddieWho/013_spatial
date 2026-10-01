#!/usr/bin/env python3
"""v5.1 conditioned-residual halo test (primary arm).

Selection effect is removed by conditioning each readout on the detection-field
level (deciles of the spot's own signature score, or of a 300 um-smoothed
signature score) plus a quadratic tissue trend. The statistic is the per-anchor
near-vs-far contrast of the residual, aggregated per section and then across
sections.

Contract: infra/tls_field_20260930/run_contract_v5.json (v5.1 addendum).
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from collections import defaultdict  # noqa: F401 (kept for parity with sibling scripts)
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/tls_field_20260930"
BIN_EDGES = np.arange(0, 1001, 100)
NBINS = len(BIN_EDGES) - 1
NEAR_I, FAR_I = 2, 7          # 200-300 vs 700-800
MIN_BIN_SPOTS = 5
MIN_ANCHORS = 8
SMOOTH_UM = 300.0
DECILES = 10


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


def anchor_bin_members(xy, comps, pitch_um):
    ncomp = len(comps)
    members = [[np.array([], dtype=int) for _ in range(NBINS)] for _ in range(ncomp)]
    if ncomp == 0:
        return members
    anchor_idx = np.concatenate([np.asarray(c) for c in comps])
    comp_of = np.repeat(np.arange(ncomp), [len(c) for c in comps])
    d, idx = cKDTree(xy[anchor_idx]).query(xy)
    cid = comp_of[idx]
    bin_id = np.digitize(d, BIN_EDGES) - 1
    valid = (bin_id >= 0) & (bin_id < NBINS) & (d >= 1.3 * pitch_um)
    buckets = [[] for _ in range(ncomp * NBINS)]
    for i in np.where(valid)[0]:
        buckets[cid[i] * NBINS + bin_id[i]].append(i)
    for c in range(ncomp):
        for b in range(NBINS):
            v = buckets[c * NBINS + b]
            if v:
                members[c][b] = np.array(v, dtype=int)
    return members


def design(s, xy, deciles=DECILES):
    """[decile dummies of s] + [1,u,v,u2,uv,v2]."""
    q = np.quantile(s, np.linspace(0, 1, deciles + 1)[1:-1])
    dec = np.digitize(s, q)
    u = (xy[:, 0] - xy[:, 0].mean()) / (np.ptp(xy[:, 0]) / 2 + 1e-9)
    v = (xy[:, 1] - xy[:, 1].mean()) / (np.ptp(xy[:, 1]) / 2 + 1e-9)
    D = np.zeros((len(s), deciles))
    D[np.arange(len(s)), dec] = 1.0
    return np.column_stack([D, np.ones(len(s)), u, v, u * u, u * v, v * v])


def residual(r, X):
    beta, *_ = np.linalg.lstsq(X, r, rcond=None)
    return r - X @ beta


def contrast(res, members):
    out = []
    for c in members:
        if len(c[NEAR_I]) >= MIN_BIN_SPOTS and len(c[FAR_I]) >= MIN_BIN_SPOTS:
            out.append(float(res[c[NEAR_I]].mean() - res[c[FAR_I]].mean()))
    if len(out) < MIN_ANCHORS:
        return float("nan"), len(out)
    return float(np.median(out)), len(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sections", default="")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    v5 = json.loads((OUT / "run_contract_v5.json").read_text())
    primary = v5["operator"]["primary_readouts"]
    pool = load_module("pool", ROOT / "scripts/tls_pool_expand.py")
    hallmark = load_module("hallmark", ROOT / "scripts/tls_field_hallmark.py")
    nullmod = load_module("nullmod", ROOT / "scripts/tls_field_null.py")
    gmt = read_gmt(ROOT / "data/geneset/h.all.v2026.1.Hs.symbols.gmt")
    proxy = json.load(open(ROOT / "infra/r04/marker_proxy_combined.json"))
    readouts = {k: v for k, v in gmt.items() if k in primary}
    readouts["B_AXIS"] = sorted(set(proxy["classes"]["B"]["voted_genes"]))

    diag = list(csv.DictReader(open(ROOT / "infra/tls_pool_expand_20260930/pool_diagnostics.tsv"),
                               delimiter="\t"))
    ok = [r for r in diag if r["status"] == "ok" and r["cohort"] != "Xenium"]
    if args.sections:
        keep = set(args.sections.split(","))
        ok = [r for r in ok if r["section_id"] in keep]

    rows = []
    for row in ok:
        sid, cohort = row["section_id"], row["cohort"]
        loaded = hallmark.load_section_by_id(pool, cohort, sid)
        if loaded is None or loaded[0] is None:
            print(f"{sid}: load failed", flush=True)
            continue
        _, X, genes, xy, pitch, labels, index = loaded
        gidx = dict(index)
        def cols_for(gl):
            cols = []
            for g in gl:
                if g in gidx:
                    v = gidx[g]
                    cols.extend(v if isinstance(v, list) else [v])
            return cols
        s = pool.aucell(X, cols_for(pool.signature_symbols()), 0.05)
        if not np.isfinite(s).all():
            print(f"{sid}: detection non-finite, skipped", flush=True)
            continue
        thr = float(np.quantile(s, 0.9))
        comps = pool.components(xy, s >= thr, pitch)
        members = anchor_bin_members(xy, comps, pitch)
        s_smooth = nullmod.smooth_op(xy, SMOOTH_UM) @ s
        X_own = design(s, xy)
        X_loc = design(s_smooth, xy)
        for tag, gl in readouts.items():
            r = pool.aucell(X, cols_for(gl), 0.05)
            if not np.isfinite(r).all():
                print(f"{sid}: {tag} non-finite, dropped", flush=True)
                continue
            raw, n_eff = contrast(r, members)
            own, _ = contrast(residual(r, X_own), members)
            loc, _ = contrast(residual(r, X_loc), members)
            rows.append({"section_id": sid, "cohort": cohort, "set": tag,
                         "raw_delta": raw, "own_cond_delta": own, "local_cond_delta": loc,
                         "n_anchors": n_eff, "n_components": len(comps)})
        print(f"{sid}: comps={len(comps)}", flush=True)

    def fmt(v):
        return "" if (v is None or not np.isfinite(v)) else round(float(v), 6)

    with (OUT / f"field_conditioned{args.tag}_v5.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["section_id", "cohort", "set", "raw_delta",
                                          "own_cond_delta", "local_cond_delta",
                                          "n_anchors", "n_components"], delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({**r, "raw_delta": fmt(r["raw_delta"]),
                        "own_cond_delta": fmt(r["own_cond_delta"]),
                        "local_cond_delta": fmt(r["local_cond_delta"])})

    summary = []
    for tag in sorted({r["set"] for r in rows}):
        sub = [r for r in rows if r["set"] == tag]
        item = {"set": tag, "n_sections": len(sub)}
        for arm in ("raw_delta", "own_cond_delta", "local_cond_delta"):
            v = np.array([r[arm] for r in sub if np.isfinite(r[arm])])
            item[f"{arm}_n"] = len(v)
            item[f"{arm}_pos_frac"] = round(float((v > 0).mean()), 3) if len(v) else ""
            item[f"{arm}_median"] = round(float(np.median(v)), 5) if len(v) else ""
        summary.append(item)
    with (OUT / f"field_conditioned_summary{args.tag}_v5.tsv").open("w", newline="") as f:
        keys = ["set", "n_sections", "raw_delta_n", "raw_delta_pos_frac", "raw_delta_median",
                "own_cond_delta_n", "own_cond_delta_pos_frac", "own_cond_delta_median",
                "local_cond_delta_n", "local_cond_delta_pos_frac", "local_cond_delta_median"]
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
        w.writeheader()
        w.writerows(summary)
    print(f"wrote rows={len(rows)} summary={len(summary)}")


if __name__ == "__main__":
    main()
