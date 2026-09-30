#!/usr/bin/env python3
"""Field screening v4: Hallmark 50 + B axis around pool anchors.

Frozen operator: infra/tls_field_20260930/run_contract_v4.json.
P0 angular bins + P1 corridor control on raw AUCell. Descriptive only:
no p-values, no inversion. P2/P3 gated on P0/P1 positive (separate run).
Xenium sections skipped (frozen).

Mask: anchor component spots + one-ring neighbors excluded from readout.
Distance: absolute micron bins (edges 0..1000 step 100), min 5 spots per bin.
Angular: 4 sectors per radial band; sector 0 centered on the direction to the
nearest tumor-labeled spot where labels exist, else +x axis (flagged
arbitrary_zero; only heterogeneity interpretable there).
"""
from __future__ import annotations

import csv
import importlib.util
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/tls_field_20260930"
V4 = OUT / "run_contract_v4.json"
BIN_EDGES = np.arange(0, 1001, 100)
MIN_BIN_SPOTS = 5
RBANDS = [(100, 200), (200, 400), (400, 700), (700, 1000)]
NSECT = 4


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_gmt(path):
    sets = {}
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def main():
    v4 = json.loads(V4.read_text())
    assert v4["status"] == "FROZEN_BEFORE_SCORING"
    pool = load_module("pool", ROOT / "scripts/tls_pool_expand.py")
    gmt = read_gmt(ROOT / "data/geneset/h.all.v2026.1.Hs.symbols.gmt")
    assert len(gmt) == 50, f"expected 50 Hallmark sets, got {len(gmt)}"
    proxy = json.load(open(ROOT / "infra/r04/marker_proxy_combined.json"))
    b_genes = sorted(set(proxy["classes"]["B"]["voted_genes"]))
    readouts = {f"HALLMARK_{k.replace('HALLMARK_', '')}": v for k, v in gmt.items()}
    readouts["B_AXIS"] = b_genes

    diag = list(csv.DictReader(open(ROOT / "infra/tls_pool_expand_20260930/pool_diagnostics.tsv"),
                               delimiter="\t"))
    ok = {r["section_id"]: r for r in diag if r["status"] == "ok"}
    # section_id -> loader call
    curves, angular, summ = [], [], []
    for sid in sorted(ok):
        cohort = ok[sid]["cohort"]
        if cohort == "Xenium":
            continue  # frozen
        out = load_section_by_id(pool, cohort, sid)
        if out is None:
            print(f"{sid}: loader unresolved, skipped", flush=True)
            continue
        _, X, genes, xy, pitch, labels, index = out
        gidx = {}
        for name, rows in index.items():
            gidx[name] = rows
        max_rank = max(1, int(np.floor(0.05 * X.shape[1])))
        scores = {}
        for tag, glist in readouts.items():
            cols = []
            for g in glist:
                if g in gidx:
                    v = gidx[g]
                    cols.extend(v if isinstance(v, list) else [v])
            scores[tag] = pool.aucell(X, cols, 0.05)
        # recompute anchors identically to pool run
        sig = pool.signature_symbols()
        scols = []
        for g in sig:
            if g in gidx:
                v = gidx[g]
                scols.extend(v if isinstance(v, list) else [v])
        sscores = pool.aucell(X, scols, 0.05)
        finite = np.isfinite(sscores)
        thr = float(np.quantile(sscores[finite], 0.9))
        mask = finite & (sscores >= thr)
        comps = pool.components(xy, mask, pitch)
        # tumor direction where labels exist
        lab = labels.get("y") if labels else None
        tum = None
        if lab is not None:
            tum = np.where(np.asarray(lab) == "TUM")[0]
            if len(tum) == 0:
                tum = None
        tree = cKDTree(xy)
        ring = 1.3 * pitch
        anchor_curves = defaultdict(list)   # (tag, bin) -> [means]
        anchor_ang = defaultdict(list)      # (tag, rband, sector) -> [means]
        n_anchors = 0
        for c in comps:
            c = np.asarray(c)
            # mask: members + one-ring neighbors
            neigh = set()
            for i in c:
                neigh.update(tree.query_ball_point(xy[i], ring))
            keep = np.ones(len(xy), dtype=bool)
            keep[list(neigh)] = False
            d = np.linalg.norm(xy - xy[c].mean(axis=0), axis=1)
            rmask = keep & np.isfinite(d)
            if rmask.sum() < 30:
                continue
            n_anchors += 1
            # sector zero: tumor direction if available
            if tum is not None and len(tum):
                t0 = xy[tum[np.argmin(np.linalg.norm(xy[tum] - xy[c].mean(axis=0), axis=1))]]
                ang0 = float(np.arctan2(t0[1] - xy[c, 1].mean(), t0[0] - xy[c, 0].mean()))
                arb = 0
            else:
                ang0, arb = 0.0, 1
            ang = np.arctan2(xy[:, 1] - xy[c, 1].mean(), xy[:, 0] - xy[c, 0].mean()) - ang0
            ang = (ang + np.pi) % (2 * np.pi) - np.pi
            for tag, v in scores.items():
                vv = v[rmask]
                dd = d[rmask]
                aa = ang[rmask]
                for bi in range(len(BIN_EDGES) - 1):
                    m = (dd >= BIN_EDGES[bi]) & (dd < BIN_EDGES[bi + 1])
                    if m.sum() >= MIN_BIN_SPOTS and np.isfinite(vv[m]).sum() >= MIN_BIN_SPOTS:
                        anchor_curves[(tag, bi)].append(float(np.nanmean(vv[m])))
                for rbi, (lo, hi) in enumerate(RBANDS):
                    m = (dd >= lo) & (dd < hi)
                    if m.sum() < MIN_BIN_SPOTS:
                        continue
                    for s_ in range(NSECT):
                        sm = m & (aa >= -np.pi + s_ * np.pi / 2) & (aa < -np.pi + (s_ + 1) * np.pi / 2)
                        if sm.sum() >= 3 and np.isfinite(vv[sm]).sum() >= 3:
                            anchor_ang[(tag, rbi, s_)].append(float(np.nanmean(vv[sm])))
                _ = arb
        for (tag, bi), vals in sorted(anchor_curves.items()):
            curves.append({"section_id": sid, "cohort": cohort, "set": tag,
                           "bin_lo": int(BIN_EDGES[bi]), "bin_hi": int(BIN_EDGES[bi + 1]),
                           "median": round(float(np.median(vals)), 5),
                           "n_anchors": len(vals)})
        for (tag, rbi, s_), vals in sorted(anchor_ang.items()):
            angular.append({"section_id": sid, "cohort": cohort, "set": tag,
                            "rband": f"{RBANDS[rbi][0]}-{RBANDS[rbi][1]}",
                            "sector": s_, "median": round(float(np.median(vals)), 5),
                            "n_anchors": len(vals)})
        print(f"{sid}: anchors={n_anchors} curves={sum(1 for k in anchor_curves if k)}", flush=True)
    with (OUT / "field_curves_v4.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["section_id", "cohort", "set", "bin_lo", "bin_hi",
                                          "median", "n_anchors"], delimiter="\t")
        w.writeheader(); w.writerows(curves)
    with (OUT / "field_angular_v4.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["section_id", "cohort", "set", "rband", "sector",
                                          "median", "n_anchors"], delimiter="\t")
        w.writeheader(); w.writerows(angular)
    # summary: radial drop (0-100 vs 700-1000) and angular spread per section/set
    by = defaultdict(list)
    for r in curves:
        by[(r["section_id"], r["set"])].append(r)
    for (sid, tag), rs in sorted(by.items()):
        near = [r["median"] for r in rs if r["bin_lo"] == 0]
        far = [r["median"] for r in rs if r["bin_lo"] == 700]
        drop = (near[0] - far[0]) if (near and far) else float("nan")
        an = defaultdict(list)
        for r in angular:
            if r["section_id"] == sid and r["set"] == tag and r["rband"] in ("200-400", "400-700"):
                an[r["rband"]].append(r["median"])
        spread = float(np.nanmax([np.nanmax(v) - np.nanmin(v) for v in an.values()])) if an else float("nan")
        summ.append({"section_id": sid, "set": tag,
                     "radial_drop_0_700": round(drop, 5) if drop == drop else "",
                     "angular_spread": round(spread, 5) if spread == spread else ""})
    with (OUT / "field_summary_v4.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["section_id", "set", "radial_drop_0_700",
                                          "angular_spread"], delimiter="\t")
        w.writeheader(); w.writerows(summ)
    print(f"wrote curves={len(curves)} angular={len(angular)} summary={len(summ)}")


def load_section_by_id(pool, cohort, sid):
    if cohort == "USZ":
        return pool.load_usz(sid.split("-", 1)[1])
    if cohort == "GSE175540":
        return pool.load_gse175540(sid.split("-", 1)[1])
    if cohort == "STCRC":
        return pool.load_stcrc(sid.split("-", 1)[1])
    if cohort == "HTAN":
        cap = sid.split("-", 1)[1]
        key_map = {"6723_KL_4": "6723_4", "8578_AS_4": "8578_4",
                   "7319_AS_2": "7319_2", "8899_AS_5": "8899_5",
                   "8899_AS_6": "8899_6", "8270_AS_2": "8270_2",
                   "7003_AS_5": "7003_5", "8899_AS_8": "8899_8"}
        if cap not in key_map:
            return None, f"HTAN capture not in v3 roster: {cap}"
        return pool.attach_htan_labels(pool.load_htan(cap), key_map[cap])
    if cohort == "Cervilla":
        if sid.endswith("Visium-v1"):
            return pool.load_cervilla_visium()
        return pool.load_cervilla_cytassist()
    if cohort == "Parent":
        if "CytAssist" in sid or "cytassist" in sid:
            return pool.load_cytassist11_crc()
        return pool.load_parent_crc()
    if cohort == "VisiumHD":
        return pool.load_visiumhd_6p5()
    return None, f"no dispatch for cohort {cohort}"


if __name__ == "__main__":
    main()
