#!/usr/bin/env python3
"""Reproduce Cho et al. Science 2026 (pan-cancer TLS atlas) Fig-4 style distance
gradients on the USZ Visium cohort.

Paper operator, as written in its methods:
  - radial distance from the TLS CENTER to every spot (Semla RadialDistance)
  - retain tumor spots only (IT TLS) / non-tumor spots only (DT TLS)
  - rescale those distances to 0-1 per TLS
  - AUCell score per spot for each pathway (rank-based)
  - Monocle natural spline df=3 -> smooth curve -> z-score per pathway
  - classify "high-to-low" if it declines from distance 0 to 0.75

This script, in addition to the paper's operator, reports SLIDING-BIN means
(bin centers 0.025..0.975, half-width 0.075). Bins are a user-requested
visual/robust readout; classification still uses the paper's endpoints
(nearest bins to 0.05 and 0.75) plus Spearman monotonicity.

Modes
  IT      a TLS component that touches a TUM spot in its first ring (<=130 um);
          retained spots = all ground_truth==TUM spots of the section.
  DTlike  a TLS component with no TUM spot within 400 um; retained spots =
          all ground_truth==NOR spots (paper's histologically-normal branch).

Not a confirmation run: no p-values, no patient split. Reproduction only.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parent.parent
H5 = ROOT / "data/other_sources/zenodo_usz_tls_visium/TLS_VISIUM_USZ/h5ad_preprocessed"
OUT = ROOT / "infra/tls_repro_20260930"
PITCH = 100.0
TOUCH = 1.3 * PITCH
FAR = 400.0
BIN_CENTERS = np.linspace(0.025, 0.975, 20)
BIN_HALF = 0.075
MIN_SPOTS_BIN = 3
MAX_RANK_FRAC = 0.05

SETS = {
    "IFN_I": "GOBP_RESPONSE_TO_TYPE_I_INTERFERON",
    "IFN_II": "GOBP_RESPONSE_TO_TYPE_II_INTERFERON",
    "MHC_II": "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II",
    "ANTIGEN_PRES": "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
    "INFLAMMATION": "GOBP_INFLAMMATORY_RESPONSE",
    "G2M": "GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION",
    "EMT": "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
}
SINGLE = ["CXCL13", "CXCL10", "CXCL11", "CXCL16", "TNFSF14", "CCL19", "CCL21",
          "MKI67", "MYC", "COL1A1", "CD74", "HLA-DRA"]
B_AXIS = "B"


def decode(x):
    return [v.decode() if isinstance(v, bytes) else v for v in x]


def load_section(alias):
    path = H5 / f"{alias}.h5ad"
    with h5py.File(path, "r") as h:
        X = h["X"][:]
        genes = decode(h["var"]["gene_name"][:])
        obs = h["obs"]
        cats = decode(obs["ground_truth"]["categories"][:])
        codes = np.asarray(obs["ground_truth"]["codes"][:])
        lab = np.array([cats[c] if c >= 0 else "" for c in codes], dtype=object)
        xa = np.asarray(obs["x_array"][:], float)
        ya = np.asarray(obs["y_array"][:], float)
        barcodes = decode(obs["_index"][:])
    xy = np.column_stack([xa + 0.5 * (ya % 2), ya * (3 ** 0.5) / 2]) * PITCH
    return X, genes, lab, xy, barcodes


def foci(xy, radius):
    parent = np.arange(len(xy))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in cKDTree(xy).query_pairs(radius):
        a, b = find(i), find(j)
        if a != b:
            parent[b] = a
    groups = {}
    for i in range(len(xy)):
        groups.setdefault(find(i), []).append(i)
    return list(groups.values())


def aucell(X, cols, max_rank, chunk=800):
    n, k = X.shape[0], len(cols)
    out = np.full(n, np.nan)
    if k == 0:
        return out
    for a in range(0, n, chunk):
        b = min(n, a + chunk)
        r = rankdata(-X[a:b], axis=1, method="average")  # rank 1 = highest expression (AUCell convention)
        pos = np.clip(max_rank - r[:, cols] + 1.0, 0.0, None)
        out[a:b] = pos.sum(axis=1) / (k * max_rank)
    return out


def gene_sets_from_local():
    defs = json.loads((ROOT / "infra/r16/recovery_20260918/programs_go/definitions.json").read_text())
    sets = {}
    for tag, rep in SETS.items():
        d = defs[rep]
        sets[tag] = sorted(set(d.get("input_genes") or []))
    proxy = json.loads((ROOT / "infra/r04/marker_proxy_combined.json").read_text())
    sets[B_AXIS] = sorted(set(proxy["classes"]["B"]["voted_genes"]))
    return sets


def classify(curve_z, centers):
    ok = np.isfinite(curve_z)
    if ok.sum() < 5:
        return "not_evaluable", np.nan, np.nan, "not_evaluable"
    i0 = int(np.nanargmin(np.abs(centers - 0.05)))
    i75 = int(np.nanargmin(np.abs(centers - 0.75)))
    delta = float(curve_z[i0] - curve_z[i75])
    rho, _ = spearmanr(centers[ok], curve_z[ok])
    rho = float(rho)
    paper = "high-to-low" if delta > 0 else "low-to-high"
    if delta > 0 and rho <= -0.5:
        lab = "high-to-low"
    elif delta < 0 and rho >= 0.5:
        lab = "low-to-high"
    elif abs(delta) <= 0.2:
        lab = "flat"
    else:
        lab = "nonmonotonic"
    return lab, delta, rho, paper


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sets = gene_sets_from_local()
    per_tls, per_curve, summ = [], [], []
    for path in sorted(H5.glob("*.h5ad")):
        alias = path.stem
        X, genes, lab, xy, _ = load_section(alias)
        gidx = {}
        for i, g in enumerate(genes):
            gidx.setdefault(g, i)
        max_rank = max(1, int(np.floor(MAX_RANK_FRAC * X.shape[1])))
        scores, coverage = {}, {}
        for tag, glist in list(sets.items()) + [("B", sets["B"])]:
            cols = [gidx[g] for g in glist if g in gidx]
            coverage[tag] = (len(cols), len(glist))
            scores[tag] = aucell(X, cols, max_rank)
        for g in SINGLE:
            if g in gidx:
                scores[g] = aucell(X, [gidx[g]], max_rank)
        tls = np.where(lab == "TLS")[0]
        tum = np.where(lab == "TUM")[0]
        nor = np.where(lab == "NOR")[0]
        tree = cKDTree(xy)
        groups = foci(xy[tls], TOUCH)
        for gi, grp in enumerate(groups):
            mem = tls[np.array(grp)]
            centroid = xy[mem].mean(axis=0)
            neigh = set()
            for i in mem:
                neigh.update(tree.query_ball_point(xy[i], TOUCH))
            named = [j for j in neigh if j not in set(mem.tolist())]
            labj = [lab[j] for j in named]
            touch_tum = "TUM" in labj
            d_tum = cKDTree(xy[tum]).query(centroid, k=1)[0] if len(tum) else np.inf
            if touch_tum:
                mode, spots = "IT", tum
            elif d_tum > FAR and len(nor) >= 30:
                mode, spots = "DTlike", nor
            else:
                continue
            if len(spots) < 30:
                continue
            dist = np.linalg.norm(xy[spots] - centroid, axis=1)
            if np.ptp(dist) <= 0:
                continue
            s = (dist - dist.min()) / np.ptp(dist)
            for tag in list(scores):
                v = scores[tag][spots]
                for c in BIN_CENTERS:
                    m = np.abs(s - c) <= BIN_HALF
                    if m.sum() < MIN_SPOTS_BIN:
                        continue
                    per_tls.append({"mode": mode, "section": alias, "tls": gi,
                                    "set": tag, "bin": round(float(c), 4),
                                    "n_spots": int(m.sum()),
                                    "mean_score": float(np.nanmean(v[m]))})
            print(f"{alias} tls={gi} {mode} um0={dist.min():.0f} um1={dist.max():.0f} n={len(spots)}", flush=True)
        meta = {"coverage": coverage, "n_tls": int((lab == "TLS").sum()),
                "n_tum": int(len(tum)), "n_nor": int(len(nor))}
        print(alias, "meta", meta, flush=True)

    rows = per_tls
    by = {}
    for r in rows:
        by.setdefault((r["mode"], r["section"], r["set"]), []).append(r)
    for (mode, section, tag), rs in sorted(by.items()):
        counts = {}
        for r in rs:
            counts.setdefault(r["bin"], []).append(r["mean_score"])
        centers = np.array(sorted(counts))
        med = np.array([np.nanmedian(counts[c]) for c in centers])
        if np.nanstd(med) > 0:
            z = (med - np.nanmean(med)) / np.nanstd(med)
        else:
            z = np.zeros_like(med)
        cancer = "KC" if section.upper().startswith("KC") else "LC"
        n_tls = len({r["tls"] for r in rs})
        for c, zz, mm in zip(centers, z, med):
            per_curve.append({"mode": mode, "cancer": cancer, "section": section,
                              "set": tag, "bin": c, "z": float(zz), "median_raw": float(mm),
                              "n_tls": n_tls})
        lab_shape, delta, rho, paper = classify(z, centers)
        summ.append({"mode": mode, "cancer": cancer, "section": section, "set": tag,
                     "n_tls": n_tls, "delta_z": delta, "spearman": rho,
                     "class": lab_shape, "paper_class": paper})

    with (OUT / "repro_tls_bins.tsv").open("w", newline="") as f:
        cols = ["mode", "section", "tls", "set", "bin", "n_spots", "mean_score"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with (OUT / "repro_curves.tsv").open("w", newline="") as f:
        cols = ["mode", "cancer", "section", "set", "bin", "z", "median_raw", "n_tls"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(per_curve)
    with (OUT / "repro_summary.tsv").open("w", newline="") as f:
        cols = ["mode", "cancer", "section", "set", "n_tls", "delta_z", "spearman", "class", "paper_class"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(summ)
    print("rows", len(rows), "curves", len(per_curve), "summary", len(summ))


if __name__ == "__main__":
    main()
