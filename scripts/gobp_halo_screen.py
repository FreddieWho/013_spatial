#!/usr/bin/env python3
"""GOBP halo screen (contract v6).

Scores the compressed keep-specific GOBP panel (1691 sets) around TLS-pool
anchors and measures the conditioned near-far contrast, with the same operator
as the v5 conditioned arm. Screen only: no p-values, no claims.

Speed design: per section the gene ranks are computed once and all sets are
scored by a single sparse product, so 1691 sets cost about as much as 1-2 sets
in the old per-set AUCell loop. Equivalence with tls_pool_expand.aucell is
asserted on a sample before any scan output is written.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/gobp_halo_20261001"
PANEL = ROOT / "infra/gobp_compress_20261001/survivors_specific_floor20_jac70.tsv"
CONTRACT = OUT / "run_contract_v6.json"
BIN_EDGES = np.arange(0, 1001, 100)
NBINS = len(BIN_EDGES) - 1
NEAR_I, FAR_I = 2, 7
MIN_BIN_SPOTS = 5
MIN_ANCHORS = 8
DECILES = 10
SMOOTH_UM = 300.0
RING_FACTOR = 1.3


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_panel(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    return rows


# ------------------------------------------------------------------ scoring
def score_sets(X, members, max_frac=0.05, chunk=600):
    """One rank pass per chunk; all sets scored by sparse product.

    Mirrors tls_pool_expand.aucell exactly: rank 1 = highest expression,
    contribution clip(max_rank - rank + 1, 0) / (k * max_rank).
    """
    n, G = X.shape
    max_rank = max(1, int(np.floor(max_frac * G)))
    nsets = len(members)
    cols = np.concatenate(members) if nsets else np.array([], dtype=np.int32)
    set_ids = np.repeat(np.arange(nsets), [len(m) for m in members])
    M = sparse.csr_matrix((np.ones(len(cols), dtype=np.float32), (cols, set_ids)),
                          shape=(G, nsets))
    k = np.bincount(set_ids, minlength=nsets).astype(np.float32)
    out = np.zeros((n, nsets), dtype=np.float32)
    for a in range(0, n, chunk):
        b = min(n, a + chunk)
        block = X[a:b].toarray() if sparse.issparse(X) else np.asarray(X[a:b])
        r = rankdata(-block, axis=1, method="average")
        pos = np.clip(max_rank - r + 1.0, 0.0, None).astype(np.float32)
        out[a:b] = (pos @ M)
    denom = np.where(k > 0, k * max_rank, np.nan)
    return out / denom[None, :], k


# --------------------------------------------------------------- geometry
def anchor_operators(xy, comps, pitch_um):
    """Voronoi edge-distance bin incidence (n x ncomp*NBINS) + counts."""
    ncomp = len(comps)
    K = ncomp * NBINS
    n = len(xy)
    if ncomp == 0:
        return sparse.csr_matrix((n, 0)), np.zeros(0)
    anchor_idx = np.concatenate([np.asarray(c) for c in comps])
    comp_of = np.repeat(np.arange(ncomp), [len(c) for c in comps])
    d, idx = cKDTree(xy[anchor_idx]).query(xy)
    cid = comp_of[idx]
    bin_id = np.digitize(d, BIN_EDGES) - 1
    valid = (bin_id >= 0) & (bin_id < NBINS) & (d >= RING_FACTOR * pitch_um)
    rows = np.where(valid)[0]
    keys = cid[valid] * NBINS + bin_id[valid]
    C = sparse.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, keys)), shape=(n, K))
    cnt = np.asarray(C.sum(axis=0)).ravel()
    return C, cnt


def design(s, xy, deciles=DECILES):
    q = np.quantile(s, np.linspace(0, 1, deciles + 1)[1:-1])
    dec = np.digitize(s, q)
    u = (xy[:, 0] - xy[:, 0].mean()) / (np.ptp(xy[:, 0]) / 2 + 1e-9)
    v = (xy[:, 1] - xy[:, 1].mean()) / (np.ptp(xy[:, 1]) / 2 + 1e-9)
    D = np.zeros((len(s), deciles), dtype=np.float64)
    D[np.arange(len(s)), dec] = 1.0
    return np.column_stack([D, np.ones(len(s)), u, v, u * u, u * v, v * v])


def residual_operator(X_design):
    # Decile indicators plus an intercept are necessarily rank deficient.
    # Unpivoted QR adds arbitrary directions at dependent/empty columns.
    X_design = np.asarray(X_design, dtype=float)
    if X_design.ndim != 2 or not np.isfinite(X_design).all():
        raise ValueError("design must be a finite matrix")
    U, s, _ = np.linalg.svd(X_design, full_matrices=False)
    tol = max(X_design.shape) * np.finfo(float).eps * (s[0] if len(s) else 0)
    return U[:, s > tol]


def residualize(R, Q):
    return R - Q @ (Q.T @ R)


def anchor_contrasts(means_kn, cnt):
    """means_kn: (ncomp*NBINS, nsets) -> per-anchor delta (nsets,) with counts guard."""
    ncomp = len(cnt) // NBINS
    M = means_kn.reshape(ncomp, NBINS, -1)
    cnt_m = cnt.reshape(ncomp, NBINS)
    ok = (cnt_m[:, NEAR_I] >= MIN_BIN_SPOTS) & (cnt_m[:, FAR_I] >= MIN_BIN_SPOTS)
    if ok.sum() < MIN_ANCHORS:
        return np.full(M.shape[2], np.nan), int(ok.sum())
    delta = np.full((ok.sum(), M.shape[2]), np.nan)
    delta = M[ok, NEAR_I, :] - M[ok, FAR_I, :]
    return np.nanmedian(delta, axis=0), int(ok.sum())


def main():
    global OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--sections", default="", help="comma list for smoke tests")
    ap.add_argument("--limit-sets", type=int, default=0)
    ap.add_argument("--tag", default="")
    ap.add_argument("--check-sets", type=int, default=5,
                    help="number of sets to verify against tls_pool_expand.aucell")
    ap.add_argument("--output-dir", type=Path, required=True, help="New audit output directory (D-167)")
    args = ap.parse_args()
    OUT = args.output_dir
    OUT.mkdir(parents=True, exist_ok=False)

    OUT.mkdir(parents=True, exist_ok=True)
    v6 = json.loads(CONTRACT.read_text())
    assert v6["status"] == "FROZEN_BEFORE_SCORING"
    pool = load_module("pool", ROOT / "scripts/tls_pool_expand.py")
    hallmark = load_module("hallmark", ROOT / "scripts/tls_field_hallmark.py")
    nullmod = load_module("nullmod", ROOT / "scripts/tls_field_null.py")
    gmt = {n: g for n, g in nullmod.read_gmt(
        Path("/home/huyudi/006/data/pathway/c5.go.v2025.1.Hs.symbols.gmt")).items()
        if n.startswith("GOBP_")}
    panel_rows = read_panel(PANEL)
    if args.limit_sets:
        panel_rows = panel_rows[:args.limit_sets]
    panel_ids = [r["program_id"] for r in panel_rows]
    panel_genes = [gmt[p] for p in panel_ids]
    sig = pool.signature_symbols()
    all_sets = panel_genes + [sig]
    rng = np.random.default_rng(20261001)

    diag = list(csv.DictReader(open(ROOT / "infra/tls_pool_expand_20260930/pool_diagnostics.tsv"),
                               delimiter="\t"))
    ok = [r for r in diag if r["status"] == "ok" and r["cohort"] != "Xenium"]
    if args.sections:
        keep = set(args.sections.split(","))
        ok = [r for r in ok if r["section_id"] in keep]

    long_rows, checked = [], False
    t_start = time.time()
    for row in ok:
        sid, cohort = row["section_id"], row["cohort"]
        t0 = time.time()
        loaded = hallmark.load_section_by_id(pool, cohort, sid)
        if loaded is None or loaded[0] is None:
            print(f"{sid}: load failed", flush=True)
            continue
        _, X, genes, xy, pitch, labels, index = loaded
        gidx = dict(index)
        members = []
        for glist in all_sets:
            idxs = []
            for g in glist:
                v = gidx.get(g)
                if v is None:
                    continue
                idxs.extend(v if isinstance(v, list) else [v])
            members.append(np.array(sorted(idxs), dtype=np.int32))
        t_load = time.time() - t0
        t1 = time.time()
        S, k = score_sets(X, members)
        t_score = time.time() - t1
        s = S[:, -1].astype(np.float64)
        R_all = S[:, :-1].astype(np.float64)
        if not np.isfinite(s).all():
            print(f"{sid}: detection score non-finite, skipped", flush=True)
            continue
        if not checked and args.check_sets > 0:
            bad = []
            for j in rng.choice(len(panel_ids), size=min(args.check_sets, len(panel_ids)),
                                replace=False):
                ref = pool.aucell(X, members[j], 0.05)
                d = float(np.nanmax(np.abs(ref - S[:, j])))
                if not np.isfinite(d) or d > 1e-6:
                    bad.append((panel_ids[j], d))
            assert not bad, f"fast scorer mismatch vs aucell: {bad}"
            print(f"equivalence check passed on {min(args.check_sets, len(panel_ids))} sets "
                  f"(max |diff| <= 1e-6)", flush=True)
            checked = True
        thr = float(np.quantile(s, 0.9))
        comps = pool.components(xy, s >= thr, pitch)
        C, cnt = anchor_operators(xy, comps, pitch)
        ncomp = len(comps)
        if ncomp == 0 or C.shape[1] == 0:
            print(f"{sid}: no anchors", flush=True)
            continue
        s_smooth = nullmod.smooth_op(xy, SMOOTH_UM) @ s
        Q_own = residual_operator(design(s, xy))
        Q_loc = residual_operator(design(s_smooth, xy))
        arms = {
            "raw": R_all,
            "own": residualize(R_all, Q_own),
            "local": residualize(R_all, Q_loc),
        }
        counts_norm = sparse.diags(1.0 / np.maximum(cnt, 1))
        A = C @ counts_norm
        deltas = {}
        for arm, R in arms.items():
            means = np.asarray((A.T @ R))
            deltas[arm], n_eff = anchor_contrasts(means, cnt)
        for j, pid in enumerate(panel_ids):
            long_rows.append((sid, cohort, pid, deltas["raw"][j], deltas["own"][j],
                              deltas["local"][j], n_eff, int(k[j])))
        print(f"{sid}: comps={ncomp} k_med={int(np.median(k))} "
              f"load={t_load:.1f}s score={t_score:.1f}s total={time.time()-t0:.1f}s", flush=True)

    tag = args.tag
    with (OUT / f"gobp_halo_sections{tag}.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["section_id", "cohort", "program_id", "raw_delta", "own_delta",
                    "local_delta", "n_anchors", "k_used"])
        for r in long_rows:
            w.writerow([r[0], r[1], r[2]] + ["" if v is None or not np.isfinite(v)
                                             else round(float(v), 6) for v in r[3:6]]
                       + [r[6], r[7]])
    # summary per set
    by_set = {}
    for r in long_rows:
        by_set.setdefault(r[2], []).append(r)
    meta = {r["program_id"]: r for r in panel_rows}
    with (OUT / f"gobp_halo_summary{tag}.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["program_id", "n_sections", "raw_pos_frac", "raw_median",
                    "own_pos_frac", "own_median", "local_pos_frac", "local_median",
                    "median_k_used", "raw_size", "covered_size", "depth"])
        for pid in panel_ids:
            rs = by_set.get(pid, [])
            line = [pid, len(rs)]
            for arm_i in (3, 4, 5):
                v = np.array([r[arm_i] for r in rs if np.isfinite(r[arm_i])], dtype=float)
                if len(v):
                    line += [round(float((v > 0).mean()), 3), round(float(np.median(v)), 6)]
                else:
                    line += ["", ""]
            kk = [r[7] for r in rs if r[7]]
            m = meta.get(pid, {})
            line += [int(np.median(kk)) if kk else "", m.get("raw_size", ""),
                     m.get("covered_size", ""), m.get("depth", "")]
            w.writerow(line)
    print(f"wrote sections={len(long_rows)} sets={len(panel_ids)} "
          f"elapsed={time.time()-t_start:.0f}s")


if __name__ == "__main__":
    main()
