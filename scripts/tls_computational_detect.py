#!/usr/bin/env python3
"""Computational TLS detection on USZ (emulate Cho et al. 2026).

Paper: iStar -> TLS signature score -> unsupervised clustering -> segmentation.
No H&E for USZ, so iStar is skipped; expression-only detection here.

Definition is frozen in detect_contract.json. This script only measures yield
and overlap with the pathologist labels; it fits no gradient and computes no
p-value.
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

import h5py
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
H5 = ROOT / "data/other_sources/zenodo_usz_tls_visium/TLS_VISIUM_USZ/h5ad_preprocessed"
OUT = ROOT / "infra/tls_repro_20260930"
PITCH = 100.0
TOUCH = 1.3 * PITCH


def decode(x):
    return [v.decode() if isinstance(v, bytes) else v for v in x]


def signature_genes():
    from r16 import axes as A
    proxy = json.loads((ROOT / "infra/r04/marker_proxy_combined.json").read_text())
    genes = set(proxy["classes"]["B"]["voted_genes"]) | set(A.PLASMA_GENES)
    genes |= {"CXCL13", "CCL19", "CCL21", "LTB"}
    return sorted(genes)


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
    xy = np.column_stack([xa + 0.5 * (ya % 2), ya * (3 ** 0.5) / 2]) * PITCH
    return X, genes, lab, xy


def zscore_log1p(X):
    Z = np.log1p(X.astype(np.float32))
    mu = Z.mean(axis=0)
    sd = Z.std(axis=0)
    sd[sd < 1e-8] = 1.0
    return (Z - mu) / sd


def aucell_score(X, cols, max_frac=0.05):
    """Rank-based AUCell (paper convention): rank 1 = highest expression."""
    from scipy.stats import rankdata
    out = np.full(X.shape[0], np.nan)
    k = len(cols)
    if k == 0:
        return out
    max_rank = max(1, int(np.floor(max_frac * X.shape[1])))
    chunk = 800
    for a in range(0, X.shape[0], chunk):
        b = min(X.shape[0], a + chunk)
        r = rankdata(-X[a:b], axis=1, method="average")
        pos = np.clip(max_rank - r[:, cols] + 1.0, 0.0, None)
        out[a:b] = pos.sum(axis=1) / (k * max_rank)
    return out


def components(xy, mask):
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return []
    parent = np.arange(len(idx))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in cKDTree(xy[idx]).query_pairs(TOUCH):
        a, b = find(i), find(j)
        if a != b:
            parent[b] = a
    groups = {}
    for i in range(len(idx)):
        groups.setdefault(find(i), []).append(idx[i])
    return list(groups.values())


def main():
    contract = json.loads((OUT / "detect_contract.json").read_text())
    assert contract["status"] == "FROZEN_BEFORE_DETECTION"
    sig = signature_genes()
    q_primary = contract["region_rule"]["Q_primary"]
    q_list = [q_primary] + list(contract["region_rule"]["Q_sensitivity"])
    rows, sweep = [], []
    for path in sorted(H5.glob("*.h5ad")):
        alias = path.stem
        X, genes, lab, xy = load_section(alias)
        gidx = {}
        for i, g in enumerate(genes):
            gidx.setdefault(g, i)
        cols = [gidx[g] for g in sig if g in gidx]
        coverage = len(cols), len(sig)
        score = aucell_score(X, cols)
        # depth check kept as a diagnostic (v1 was rejected for dependence on depth)
        Z = zscore_log1p(X[:, cols])
        from scipy.stats import spearmanr
        lib = X.sum(axis=1)
        rho_depth = float(spearmanr(score, lib)[0])
        tls = np.where(lab == "TLS")[0]
        tum = np.where(lab == "TUM")[0]
        nor = np.where(lab == "NOR")[0]
        # detector quality against pathologist labels
        from scipy.stats import rankdata
        ref = np.concatenate([tls, nor])
        r = rankdata(score[ref])
        n1, n2 = len(tls), len(nor)
        auc = (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n2) if n1 and n2 else np.nan
        for q in q_list:
            thr = float(np.quantile(score, q))
            comps = components(xy, score >= thr)
            comps = [c for c in comps if len(c) >= contract["region_rule"]["min_spots"]]
            modes = Counter()
            new_it = 0
            for c in comps:
                d = cKDTree(xy[tum]).query(xy[c].mean(axis=0), k=1)[0] if len(tum) else np.inf
                touch = cKDTree(xy[tum]).query(xy[c], k=1)[0].min() if len(tum) else np.inf
                mode = "IT" if touch <= TOUCH else ("near" if d <= 400 else "distal")
                modes[mode] += 1
                if mode == "IT":
                    known = float(np.isin(c, tls).mean())
                    if known == 0.0:
                        new_it += 1
                    if q == q_primary:
                        rows.append({
                            "section": alias, "comp": len(rows), "n_spots": len(c),
                            "mean_score": float(score[c].mean()),
                            "cx_um": float(xy[c, 0].mean()), "cy_um": float(xy[c, 1].mean()),
                            "min_dist_tum_um": float(touch), "centroid_dist_tum_um": float(d),
                            "mode": mode, "known_tls_frac": known, "is_new": known == 0.0,
                        })
            sweep.append({"section": alias, "q": q, "threshold": thr, "n_components": len(comps),
                          "IT": modes.get("IT", 0), "near": modes.get("near", 0),
                          "distal": modes.get("distal", 0), "IT_new": new_it})
        print(f"{alias} sig={coverage[0]}/{coverage[1]} auc_tls_vs_nor={auc:.3f} "
              f"rho_score_lib={rho_depth:+.2f} tls_spots={len(tls)} tum_spots={len(tum)}", flush=True)

    with (OUT / "detect_sweep.tsv").open("w", newline="") as f:
        cols = ["section", "q", "threshold", "n_components", "IT", "near", "distal", "IT_new"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(sweep)
    with (OUT / "detect_components_it.tsv").open("w", newline="") as f:
        if rows:
            cols = list(rows[0])
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    print(json.dumps({"Q_primary": q_primary, "IT_components": len(rows),
                      "IT_new": sum(1 for r in rows if r["is_new"])}, indent=2))
    for r in sweep:
        if r["q"] == q_primary:
            print(r)


if __name__ == "__main__":
    main()
