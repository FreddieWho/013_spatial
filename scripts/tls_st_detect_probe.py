#!/usr/bin/env python3
"""Probe of citation ST-based TLS detection on labeled sections.

Pipeline (adapted; deviations recorded):
- TLS score per spot = mean log1p(count+1) of repo TLS_signature genes
  (paper uses iStar-imputed pixel expression; we use spot level, no NN).
- Smooth over spatial graph neighbors (paper: 32x32 px moving average).
- k-means k=5; TLS cluster = highest-mean-score cluster (frozen rule,
  chosen before seeing labels).
- Cross-check with spot-level CD3D/MS4A1 (no imputation), as in paper.
- Connected components on selected mask; discard <3 spots.
- Maturation: per-component 95th pct CR2/FCER2 -> kmeans k=2 (E vs m);
  MKI67/BCL6 -> P vs S.
Compares against independent author labels. Weak-label use only.
"""
import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.metrics import confusion_matrix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.tls_peripheral import source_metadata, counts  # noqa: E402

SIG = [l.strip() for l in (ROOT / "repo/data/TLS_signature/genes.txt").read_text().splitlines() if l.strip()]


def tls_score(X, genes):
    idx = {g: i for i, g in enumerate(genes)}
    keep = [idx[g] for g in SIG if g in idx]
    print(f"signature genes present: {len(keep)}/{len(SIG)}", flush=True)
    sub = X[:, keep].toarray() if sparse.issparse(X) else np.asarray(X[:, keep])
    return np.log1p(sub).mean(axis=1)


def smooth(y, xy, radius=200.0):
    from scipy.spatial import cKDTree
    tree = cKDTree(xy)
    nbrs = tree.query_ball_point(xy, radius)
    out = np.array([y[n].mean() for n in nbrs])
    return out


def components(xy, mask, radius=130.01):
    from scipy.spatial import cKDTree
    pts = np.flatnonzero(mask)
    if len(pts) == 0:
        return []
    parent = {int(i): int(i) for i in pts}

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    if len(pts) > 1:
        for a, b in cKDTree(xy[points := pts]).query_pairs(radius):
            ia, ib = int(pts[a]), int(pts[b])
            ra, rb = find(ia), find(ib)
            if ra != rb:
                parent[rb] = ra
    groups = {}
    for i in pts:
        groups.setdefault(find(int(i)), []).append(int(i))
    return [np.array(v) for v in groups.values() if len(v) >= 3]


def run(gsm, cohort="GSE175540"):
    sid = f"gse175540-{gsm}" if cohort == "GSE175540" else gsm
    m = source_metadata(cohort, sid)
    X, genes = counts(cohort, sid, m)
    score = tls_score(X, np.array(genes))
    sm = smooth(score, m["xy"])
    km = KMeans(n_clusters=5, n_init=10, random_state=0).fit(sm.reshape(-1, 1))
    order = np.argsort([sm[km.labels_ == k].mean() for k in range(5)])
    tls_cluster = int(order[-1])  # frozen rule: top-score cluster
    mask = km.labels_ == tls_cluster
    # cross-check markers (expression only, no labels)
    gidx = {g: i for i, g in enumerate(genes)}
    cd3 = np.asarray(X[:, gidx["CD3D"]].todense()).ravel() if "CD3D" in gidx else np.zeros(len(mask))
    ms = np.asarray(X[:, gidx["MS4A1"]].todense()).ravel() if "MS4A1" in gidx else np.zeros(len(mask))
    print(f"{gsm}: n={len(mask)} frac_selected={mask.mean():.3f} "
          f"CD3D_in={np.log1p(cd3[mask]).mean():.2f} CD3D_out={np.log1p(cd3[~mask]).mean():.2f} "
          f"MS4A1_in={np.log1p(ms[mask]).mean():.2f} MS4A1_out={np.log1p(ms[~mask]).mean():.2f}", flush=True)
    comps = components(m["xy"], mask)
    print(f"{gsm}: components={len(comps)} sizes={sorted([len(c) for c in comps], reverse=True)[:8]}", flush=True)
    # audit vs independent labels
    lab = m["labels"]
    known = np.isin(lab, ["TLS", "NO_TLS"])
    y = (lab[known] == "TLS").astype(int)
    p = mask[known].astype(int)
    tn, fp, fn, tp = confusion_matrix(y, p, labels=[0, 1]).ravel()
    sens = tp / max(1, tp + fn)
    spec = tn / max(1, tn + fp)
    print(f"{gsm}: AUDIT sens={sens:.3f} spec={spec:.3f} tp={tp} fp={fp} fn={fn} tn={tn}", flush=True)
    # component hit rate: fraction of true TLS spots covered
    cov = mask[lab == "TLS"].mean() if (lab == "TLS").any() else float("nan")
    print(f"{gsm}: TLS-spot coverage={cov:.3f}", flush=True)
    # maturation on detected components
    for name, markers in (("FDC", ["CR2", "FCER2"]), ("GC", ["MKI67", "BCL6"])):
        vals = []
        for c in comps:
            sub = []
            for g in markers:
                if g in gidx:
                    sub.append(np.asarray(X[c][:, gidx[g]].todense()).ravel())
            if sub:
                vals.append(np.percentile(np.log1p(np.mean(sub, axis=0)), 95))
        vals = np.array(vals)
        if len(vals) >= 2:
            lab2 = KMeans(n_clusters=2, n_init=10, random_state=0).fit(vals.reshape(-1, 1)).labels_
            print(f"{gsm}: {name} split={np.bincount(lab2).tolist()}", flush=True)
    return {"gsm": gsm, "sens": sens, "spec": spec, "n_comp": len(comps), "coverage": float(cov)}


if __name__ == "__main__":
    out = [run(g) for g in sys.argv[1:]]
    print(json.dumps(out, indent=1))
