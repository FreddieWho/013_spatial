#!/usr/bin/env python3
"""Supervised expansion (D-212/D-213): score 5 new GSE sections with the
frozen 1742 panel + geometry, and label maturation (E/P/S) on
independently-defined TLS instances across all sections.

No ST-detection labels are produced or used for training.
"""
import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.tls_peripheral import source_metadata, counts  # noqa: E402
from scripts.tls_multiscale_gobp import build_endpoints  # noqa: E402
from scripts.gobp_halo_screen import score_sets  # noqa: E402

OLD = ROOT / "infra/tls_multiscale_gobp_20261002"
OUT = ROOT / "infra/tls_expand_20261006"
# D-214: only TLS_2_cat-schema sections enter the supervised core.
NEW_GSM = ["GSM5924036_ffpe_c_21", "GSM5924048_frozen_b_13"]
NEW_SID = ["gse175540-" + g for g in NEW_GSM]

FDC = ["CR2", "FCER2"]
GCM = ["MKI67", "BCL6"]


def score_new():
    (OUT / "scores").mkdir(parents=True, exist_ok=True)
    (OUT / "geometry").mkdir(parents=True, exist_ok=True)
    panel = json.loads((OLD / "panel_gobp.json").read_text())
    control = json.loads((OLD / "panel_controls.json").read_text())
    allsets = panel | control
    assert len(allsets) == 1742, len(allsets)
    for sid in NEW_SID:
        if (OUT / "scores" / f"{sid}.npz").exists() and (OUT / "geometry" / f"{sid}.npz").exists():
            print("skip completed", sid, flush=True)
            continue
        print("scoring", sid, flush=True)
        m = source_metadata("GSE175540", sid)
        X, genes = counts("GSE175540", sid, m)
        index = {g: i for i, g in enumerate(genes)}
        members = [np.array([index[g] for g in gl if g in index], dtype=int) for gl in allsets.values()]
        S, k = score_sets(X, members)
        fraction = k / np.array([len(v) for v in allsets.values()])
        valid = (fraction >= .8) & (k >= 10) & (np.nanstd(S, axis=0) > 1e-12)
        lib = np.asarray(X.sum(axis=1)).ravel().astype(np.float32)
        np.savez_compressed(OUT / "scores" / f"{sid}.npz", scores=S.astype(np.float32),
                            set_ids=np.array(list(allsets)), coverage=fraction, n_genes=k,
                            valid=valid, libsize=lib, detected=np.diff(X.indptr))
        graph, ends, members_c, edge = build_endpoints(m["xy"], m["labels"])
        sparse.save_npz(OUT / "geometry" / f"{sid}_graph.npz", graph)
        np.savez_compressed(OUT / "geometry" / f"{sid}.npz", xy=m["xy"], labels=m["labels"],
                            order=m["order"], edge=edge,
                            endpoint_members=np.concatenate(members_c) if members_c else np.array([], int),
                            endpoint_indptr=np.r_[0, np.cumsum([len(v) for v in members_c])])
        print("done", sid, S.shape, flush=True)


def tls_components(xy, labels):
    from scipy.spatial import cKDTree
    pts = np.flatnonzero(labels == "TLS")
    if len(pts) == 0:
        return []
    parent = {int(i): int(i) for i in pts}

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    if len(pts) > 1:
        for a, b in cKDTree(xy[pts]).query_pairs(130.01):
            ia, ib = int(pts[a]), int(pts[b])
            ra, rb = find(ia), find(ib)
            if ra != rb:
                parent[rb] = ra
    groups = {}
    for i in pts:
        groups.setdefault(find(int(i)), []).append(int(i))
    return [np.array(v) for v in groups.values() if len(v) >= 2]


def mature_all():
    import csv
    import json as js
    idrows = list(csv.DictReader(open(ROOT / "infra/bioinf-data-index/raw/user_gse175540_20261003/current_26_source_identity_confirmed.tsv"), delimiter="\t"))
    sections = [(r["section_id"], r["cohort"]) for r in idrows] + [(s, "GSE175540") for s in NEW_SID]
    from sklearn.cluster import KMeans
    rows, feats_fdc, feats_gc, keys = [], [], [], []
    cache = {}
    for sid, cohort in sections:
        geo = OLD / "geometry" if (OLD / "geometry" / f"{sid}.npz").exists() else OUT / "geometry"
        g = np.load(geo / f"{sid}.npz")
        xy, labels = np.asarray(g["xy"]), np.array([str(x) for x in g["labels"]])
        comps = tls_components(xy, labels)
        m = source_metadata(cohort, sid)
        X, genes = counts(cohort, sid, m)
        gidx = {g_: i for i, g_ in enumerate(genes)}
        for ci, c in enumerate(comps):
            vals = {}
            for name, markers in (("FDC", FDC), ("GC", GCM)):
                sub = [np.asarray(X[c][:, gidx[g_]].todense()).ravel() for g_ in markers if g_ in gidx]
                vals[name] = float(np.percentile(np.log1p(np.mean(sub, axis=0)), 95)) if sub else float("nan")
            rows.append({"section_id": sid, "cohort": cohort, "comp": ci, "n_spots": len(c),
                         "cx": float(xy[c][:, 0].mean()), "cy": float(xy[c][:, 1].mean()),
                         "fdc95": vals["FDC"], "gc95": vals["GC"]})
            feats_fdc.append(vals["FDC"])
            feats_gc.append(vals["GC"])
            keys.append(len(rows) - 1)
    feats_fdc = np.array(feats_fdc)
    feats_gc = np.array(feats_gc)
    ok = np.isfinite(feats_fdc)
    lab_fdc = np.full(len(rows), -1)
    lab_fdc[ok] = KMeans(n_clusters=2, n_init=20, random_state=0).fit(feats_fdc[ok].reshape(-1, 1)).labels_
    # orient: higher FDC = mature
    if feats_fdc[lab_fdc == 0].mean() > feats_fdc[lab_fdc == 1].mean():
        lab_fdc = 1 - lab_fdc
    mature_mask = lab_fdc == 1
    lab_gc = np.full(len(rows), -1)
    if mature_mask.sum() >= 2 and np.isfinite(feats_gc[mature_mask]).sum() >= 2:
        sub = np.where(mature_mask)[0]
        vv = feats_gc[sub]
        ok2 = np.isfinite(vv)
        l2 = KMeans(n_clusters=2, n_init=20, random_state=0).fit(vv[ok2].reshape(-1, 1)).labels_
        if vv[ok2][l2 == 0].mean() > vv[ok2][l2 == 1].mean():
            l2 = 1 - l2
        lab_gc[sub[ok2]] = l2
    for i, r in enumerate(rows):
        if lab_fdc[i] == 0:
            r["maturation"] = "E-TLS"
        elif lab_gc[i] == 0:
            r["maturation"] = "P-TLS"
        elif lab_gc[i] == 1:
            r["maturation"] = "S-TLS"
        else:
            r["maturation"] = "mTLS-unresolved"
    with open(OUT / "tls_maturation.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    print("instances:", len(rows), Counter(r["maturation"] for r in rows), flush=True)
    (OUT / "STATUS").write_text("EXPAND_COMPLETE\n")


if __name__ == "__main__":
    t0 = time.time()
    score_new()
    mature_all()
    print(f"elapsed {time.time()-t0:.0f}s", flush=True)
