#!/usr/bin/env python3
"""Discovery-patient TLS field curves. Rules are read from the frozen contract.

Does not refit the split, gene sets, or shape cut. Does not pool scaled
distances. Does not compute a spatial p-value.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse
from scipy.interpolate import BSpline
from scipy.spatial import cKDTree

from r16.section_io import load_section_symbols

ROOT = Path(__file__).resolve().parent.parent
CONTRACT = ROOT / "infra/tls_field_20260930/run_contract.json"
OUT = ROOT / "infra/tls_field_20260930"


def hex_um(xa, ya, pitch):
    x = xa + 0.5 * (ya % 2)
    y = ya * (3 ** 0.5) / 2
    return np.column_stack([x, y]) * pitch


def foci(xy, radius):
    parent = np.arange(len(xy))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    if len(xy):
        for i, j in cKDTree(xy).query_pairs(radius):
            a, b = find(i), find(j)
            if a != b:
                parent[b] = a
    groups = {}
    for i in range(len(xy)):
        groups.setdefault(find(i), []).append(i)
    return list(groups.values())


def labels(path, barcodes):
    with h5py.File(path, "r") as h:
        cats = [c.decode() if isinstance(c, bytes) else c
                for c in h["obs"]["ground_truth"]["categories"][:]]
        codes = np.asarray(h["obs"]["ground_truth"]["codes"][:])
        idx = [b.decode() if isinstance(b, bytes) else b for b in h["obs"]["_index"][:]]
        xa = np.asarray(h["obs"]["x_array"][:], float)
        ya = np.asarray(h["obs"]["y_array"][:], float)
    lut = {b: (cats[codes[i]] if codes[i] >= 0 else "", xa[i], ya[i])
           for i, b in enumerate(idx)}
    lab = []
    xy = np.zeros((len(barcodes), 2))
    for i, b in enumerate(barcodes):
        key = b.decode() if isinstance(b, bytes) else str(b)
        item = lut.get(key)
        if item is None:
            lab.append("")
        else:
            lab.append(item[0])
            xy[i] = (item[1], item[2])
    return np.array(lab), xy


def aucell(counts, cols, max_rank):
    dense = counts.toarray() if sparse.issparse(counts) else np.asarray(counts)
    n_genes = dense.shape[1]
    ranks = np.argsort(np.argsort(-dense, axis=1), axis=1) + 1
    k = len(cols)
    if k == 0:
        return np.full(dense.shape[0], np.nan)
    in_top = ranks[:, cols] <= max_rank
    # recovery curve area inside the top window
    pos = np.clip(max_rank - ranks[:, cols] + 1, 0, None)
    u = pos.sum(axis=1)
    return u / (k * max_rank)


def ns_basis(x):
    x = np.asarray(x, float)
    a, b = float(np.min(x)), float(np.max(x))
    if b <= a:
        b = a + 1e-6
    k1, k2 = np.quantile(x, [1/3, 2/3])
    knots = np.array([a, a, a, a, k1, k2, b, b, b, b])
    eye = np.eye(len(knots) - 4)
    B = np.column_stack([BSpline(knots, eye[i], 3)(x) for i in range(eye.shape[0])])
    # natural-style reduction: intercept plus two interior contrasts
    return np.column_stack([np.ones(len(x)), B[:, 1:3]])


def shape_of(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 30:
        return "not_evaluable", [np.nan] * 4
    X = ns_basis(x[ok])
    beta, _, _, _ = np.linalg.lstsq(X, y[ok], rcond=None)
    grid = np.array([0.0, 0.25, 0.5, 0.75])
    # evaluate with the same knots as the fit
    a, b = float(np.min(x[ok])), float(np.max(x[ok]))
    k1, k2 = np.quantile(x[ok], [1/3, 2/3])
    knots = np.array([a, a, a, a, k1, k2, b, b, b, b])
    eye = np.eye(len(knots) - 4)
    Bg = np.column_stack([BSpline(knots, eye[i], 3)(grid) for i in range(eye.shape[0])])
    Xg = np.column_stack([np.ones(4), Bg[:, 1:3]])
    fit = Xg @ beta
    drop = float(fit[0] - fit[3])
    if np.all(np.diff(fit) <= 1e-6) and drop >= 0.02:
        return "high-to-low", fit
    if np.all(np.diff(fit) >= -1e-6) and -drop >= 0.02:
        return "low-to-high", fit
    return "flat/nonmonotonic", fit


def main():
    contract = json.loads(CONTRACT.read_text())
    man = json.loads((ROOT / "infra/r04/role_manifests/external_validation_manifest.json").read_text())
    by_patient = {r["patient_id"]: r for r in man["rows"]}
    rows = []
    for pid in contract["discovery_patients"]:
        sec = load_section_symbols(by_patient[pid])
        alias = sec.section_id.split("::")[-1]
        path = ROOT / "data/other_sources/zenodo_usz_tls_visium/TLS_VISIUM_USZ/h5ad_preprocessed" / f"{alias}.h5ad"
        lab, array_xy = labels(path, sec.barcode)
        xy = hex_um(array_xy[:, 0], array_xy[:, 1], contract["pitch_um"])
        tls_ix = np.where(lab == "TLS")[0]
        groups = foci(xy[tls_ix], 1.3 * contract["pitch_um"])
        counts = sec.counts.tocsr() if sparse.issparse(sec.counts) else sparse.csr_matrix(sec.counts)
        genes = [g if isinstance(g, str) else g.decode() for g in sec.gene_id]
        gidx = {}
        for i, g in enumerate(genes):
            gidx.setdefault(g, i)
        max_rank = max(1, int(np.floor(0.05 * counts.shape[1])))
        scores = {}
        for name, glist in {**contract["gene_sets"], "B": contract["b_axis_genes"]}.items():
            cols = [gidx[g] for g in glist if g in gidx]
            scores[name] = aucell(counts, cols, max_rank)
        for gi, group in enumerate(groups):
            members = tls_ix[group]
            neigh = set()
            tree = cKDTree(xy)
            for i in members:
                neigh.update(tree.query_ball_point(xy[i], 1.3 * contract["pitch_um"]))
            neigh = [j for j in neigh if j not in set(members)]
            named = [j for j in neigh if lab[j] in ("TUM", "NOR", "INFL", "TLS")]
            touches_tum = any(lab[j] == "TUM" for j in named)
            has_nor = any(lab[j] == "NOR" for j in named)
            if not (touches_tum and not has_nor):
                continue
            tumor = np.where(lab == "TUM")[0]
            dist = cKDTree(xy[members]).query(xy[tumor], k=1)[0]
            if len(dist) < 30 or np.ptp(dist) <= 0:
                continue
            scaled = (dist - dist.min()) / np.ptp(dist)
            if int((scaled >= 0.5).sum()) < 8:
                continue
            um0 = float(dist.min())
            um1 = float(dist.max())
            um75 = um0 + 0.75 * (um1 - um0)
            for name, score in scores.items():
                kind, fit = shape_of(scaled, score[tumor])
                rows.append({
                    "patient": pid, "section": alias, "tls": gi, "set": name,
                    "shape": kind, "n_tumor": int(len(tumor)),
                    "um_0": round(um0, 1), "um_0.75": round(um75, 1), "um_1": round(um1, 1),
                    "fit_0": round(float(fit[0]), 4), "fit_0.25": round(float(fit[1]), 4),
                    "fit_0.5": round(float(fit[2]), 4), "fit_0.75": round(float(fit[3]), 4),
                })
        print(alias, "it_curves", sum(1 for r in rows if r["section"] == alias), flush=True)
    out = OUT / "discovery_curves.tsv"
    with out.open("w", newline="") as f:
        cols = ["patient", "section", "tls", "set", "shape", "n_tumor",
                "um_0", "um_0.75", "um_1", "fit_0", "fit_0.25", "fit_0.5", "fit_0.75"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("wrote", out, "rows", len(rows))


if __name__ == "__main__":
    main()
