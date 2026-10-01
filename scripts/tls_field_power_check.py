#!/usr/bin/env python3
"""Power check for the conditioned halo statistic.

Injects a synthetic halo (constant bump in the 200-300 um ring around real
anchors) into a readout and verifies that the conditioned statistic detects it.
Purpose: show the negative v5 result is not a dead pipeline.

Usage: python3 scripts/tls_field_power_check.py --sections a,b,c --bump 0.5 --bump 1.0
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/tls_field_20260930"
BIN_EDGES = np.arange(0, 1001, 100)
NEAR_I, FAR_I = 2, 7


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sections", required=True)
    ap.add_argument("--bumps", default="0.25,0.5,1.0")
    ap.add_argument("--readouts", default="HALLMARK_INTERFERON_GAMMA_RESPONSE,B_AXIS")
    args = ap.parse_args()
    bumps = [float(b) for b in args.bumps.split(",")]
    wanted = set(args.readouts.split(","))

    pool = load_module("pool", ROOT / "scripts/tls_pool_expand.py")
    hallmark = load_module("hallmark", ROOT / "scripts/tls_field_hallmark.py")
    cond = load_module("cond", ROOT / "scripts/tls_field_conditioned.py")
    gmt = cond.read_gmt(ROOT / "data/geneset/h.all.v2026.1.Hs.symbols.gmt")
    proxy = json.load(open(ROOT / "infra/r04/marker_proxy_combined.json"))
    readouts = {k: v for k, v in gmt.items() if k in wanted}
    if "B_AXIS" in wanted:
        readouts["B_AXIS"] = sorted(set(proxy["classes"]["B"]["voted_genes"]))

    diag = list(csv.DictReader(open(ROOT / "infra/tls_pool_expand_20260930/pool_diagnostics.tsv"),
                               delimiter="\t"))
    want = set(args.sections.split(","))
    rows = []
    for row in [r for r in diag if r["section_id"] in want and r["status"] == "ok"]:
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
        comps = pool.components(xy, s >= float(np.quantile(s, 0.9)), pitch)
        members = cond.anchor_bin_members(xy, comps, pitch)
        Xd = cond.design(s, xy)
        # distance to nearest anchor spot, for the injection
        d_near = cKDTree(xy[np.concatenate(comps)]).query(xy)[0]
        for tag, gl in readouts.items():
            r = pool.aucell(X, cols_for(gl), 0.05)
            base, n_eff = cond.contrast(cond.residual(r, Xd), members)
            for b in bumps:
                inj = r + b * np.std(r) * ((d_near > 200) & (d_near <= 300)).astype(float)
                stat, _ = cond.contrast(cond.residual(inj, Xd), members)
                rows.append({"section_id": sid, "cohort": cohort, "set": tag,
                             "bump_sd": b, "base": base, "injected": stat,
                             "detected": int(np.isfinite(stat) and stat > 0.002),
                             "n_anchors": n_eff})
        print(f"{sid}: comps={len(comps)}", flush=True)

    with (OUT / "field_power_check_v5.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["section_id", "cohort", "set", "bump_sd", "base",
                                          "injected", "detected", "n_anchors"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for tag in sorted({r["set"] for r in rows}):
        for b in bumps:
            sub = [r for r in rows if r["set"] == tag and r["bump_sd"] == b]
            print(f"{tag[:40]:40} bump={b:<5} detected {sum(r['detected'] for r in sub)}/{len(sub)} "
                  f"median_injected={np.median([r['injected'] for r in sub]):+.5f} "
                  f"median_base={np.median([r['base'] for r in sub]):+.5f}")


if __name__ == "__main__":
    main()
