#!/usr/bin/env python3
"""Fig-4 style gradients anchored on COMPUTATIONAL TLS components.

Anchors: expression-only components (AUCell B+plasma+chemokine score above the
per-section quantile Q, hex-connected, >=2 spots) that touch a TUM spot within
130 um. Anchor spots are EXCLUDED from the readout set (input/readout disjoint).

Readout: all TUM spots of the section minus the anchor's own spots. Distance
from anchor centroid -> rescale 0-1 -> AUCell -> sliding bins. Section curve is
the median across anchors; per-anchor endpoint deltas are reported separately so
a single anchor cannot carry the result.

Circularity: detection uses B/plasma/chemokines, so those readouts are circular.
Tumor-intrinsic readouts (G2M, EMT, MYC, MKI67, COL1A1) are the clean test.
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
OUT = ROOT / "infra/tls_repro_20260930"
H5 = ROOT / "data/other_sources/zenodo_usz_tls_visium/TLS_VISIUM_USZ/h5ad_preprocessed"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    repro = load_module("repro", ROOT / "scripts/tls_pancancer_repro.py")
    detect = load_module("detect", ROOT / "scripts/tls_computational_detect.py")
    contract = json.loads((OUT / "detect_contract.json").read_text())
    Q = contract["region_rule"]["Q_primary"]
    sig = detect.signature_genes()
    sets = repro.gene_sets_from_local()

    per_comp, per_bin = [], []
    for path in sorted(H5.glob("*.h5ad")):
        alias = path.stem
        X, genes, lab, xy, _ = repro.load_section(alias)
        gidx = {}
        for i, g in enumerate(genes):
            gidx.setdefault(g, i)
        bcols = [gidx[g] for g in sig if g in gidx]
        bscore = detect.aucell_score(X, bcols)
        thr = float(np.quantile(bscore, Q))
        comps = [c for c in detect.components(xy, bscore >= thr) if len(c) >= 2]
        tum = np.where(lab == "TUM")[0]
        tls = np.where(lab == "TLS")[0]
        if not len(tum):
            continue
        max_rank = max(1, int(np.floor(repro.MAX_RANK_FRAC * X.shape[1])))
        readout_scores = {}
        for tag, glist in sets.items():
            cols = [gidx[g] for g in glist if g in gidx]
            readout_scores[tag] = repro.aucell(X, cols, max_rank)
        for g in repro.SINGLE:
            if g in gidx:
                readout_scores[g] = repro.aucell(X, [gidx[g]], max_rank)
        kept = 0
        for ci, c in enumerate(comps):
            touch = cKDTree(xy[tum]).query(xy[c], k=1)[0].min()
            if touch > detect.TOUCH:
                continue
            readout = np.setdiff1d(tum, c)  # disjoint from the anchor
            if len(readout) < 30:
                continue
            centroid = xy[c].mean(axis=0)
            dist = np.linalg.norm(xy[readout] - centroid, axis=1)
            if np.ptp(dist) <= 0:
                continue
            s = (dist - dist.min()) / np.ptp(dist)
            known = float(np.isin(c, tls).mean()) if len(tls) else 0.0
            kept += 1
            for tag, v in readout_scores.items():
                vv = v[readout]
                vals = np.array([
                    float(np.nanmean(vv[np.abs(s - b) <= repro.BIN_HALF]))
                    if (np.abs(s - b) <= repro.BIN_HALF).sum() >= repro.MIN_SPOTS_BIN else np.nan
                    for b in repro.BIN_CENTERS
                ])
                if np.isfinite(vals).sum() < 5:
                    continue
                z = (vals - np.nanmean(vals)) / (np.nanstd(vals) + 1e-12)
                centers = np.array(repro.BIN_CENTERS)
                lab_shape, delta, rho, paper = repro.classify(z, centers)
                per_comp.append({
                    "section": alias, "comp": ci, "set": tag, "n_spots": len(c),
                    "known_tls_frac": known, "is_new": known == 0.0,
                    "um1": float(dist.max()), "delta_z": delta,
                    "class": lab_shape, "paper_class": paper,
                })
                for b, zz in zip(centers, z):
                    per_bin.append({"section": alias, "comp": ci, "set": tag,
                                    "bin": b, "z": zz})
        print(f"{alias}: anchors={kept} (of {len(comps)} comps)", flush=True)

    cur = defaultdict(list)
    for r in per_bin:
        cur[(r["section"], r["set"], r["bin"])].append(r["z"])
    per_curve = []
    for (sec, tag, b), zs in sorted(cur.items()):
        per_curve.append({"section": sec, "set": tag, "bin": b,
                          "z_median": float(np.median(zs)), "n_comp": len(zs)})
    summ = []
    by = defaultdict(list)
    for r in per_curve:
        by[(r["section"], r["set"])].append(r)
    for (sec, tag), rs in sorted(by.items()):
        rs = sorted(rs, key=lambda r: r["bin"])
        z = np.array([r["z_median"] for r in rs])
        centers = np.array([r["bin"] for r in rs])
        lab_shape, delta, rho, paper = repro.classify(z, centers)
        pc = [r for r in per_comp if r["section"] == sec and r["set"] == tag]
        d = np.array([r["delta_z"] for r in pc], dtype=float)
        summ.append({
            "section": sec, "set": tag, "n_anchors": len(pc),
            "delta_z": delta, "class": lab_shape, "paper_class": paper,
            "anchor_delta_median": float(np.median(d)) if len(d) else np.nan,
            "anchor_pos_frac": float((d > 0).mean()) if len(d) else np.nan,
        })
    with (OUT / "comp_anchor_percomp.tsv").open("w", newline="") as f:
        cols = ["section", "comp", "set", "n_spots", "known_tls_frac", "is_new",
                "um1", "delta_z", "class", "paper_class"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(per_comp)
    with (OUT / "comp_anchor_curves.tsv").open("w", newline="") as f:
        cols = ["section", "set", "bin", "z_median", "n_comp"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(per_curve)
    with (OUT / "comp_anchor_summary.tsv").open("w", newline="") as f:
        cols = ["section", "set", "n_anchors", "delta_z", "class", "paper_class",
                "anchor_delta_median", "anchor_pos_frac"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(summ)
    print("per_comp", len(per_comp), "curves", len(per_curve), "summary", len(summ))


if __name__ == "__main__":
    main()
