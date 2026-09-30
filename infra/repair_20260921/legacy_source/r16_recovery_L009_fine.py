#!/usr/bin/env python3
"""L-009: P-mhc2 TLS指向验证 — TLS-foci自适应细窗 + 多种子稳定性。

设计（相对 T5B 的两处改变）：
1. 细窗：以每个 TLS focus 为中心 titration 三档半窗 (3/5/8 hops)，
   阴性窗=同切片等大小随机窗（GT无关生成，事后套GT定阴阳，03 §7）。
2. 多种子：3 seeds × 3档 = 每focus 9次测量；head-to-head: mhc2 vs B轴，
   同窗同seed配对比较，符号检验 + 中位差。
3. 终点：隐藏窗均值分排序阳性/阴性窗的 AUC（与T5B同口径，可比）；
   新增：focus级命中率（阳性窗分>阴性窗分位数）。
输出：masked_structure_fine.tsv, L009_verdict.md
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from r16 import axes as A
from r16 import census as C
from r16.recovery import mask_tasks as T
from r16.section_io import load_section_symbols, usz_tls_labels
from r04.structure_readout import _auc

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/r16/recovery_20260918"


def prog_score(counts, gidx: dict, genes: list[str]) -> np.ndarray:
    cols = [gidx[g] for g in genes if g in gidx]
    if not cols:
        return np.full(counts.shape[0], np.nan)
    X = counts[:, cols]
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    return np.log1p(X.astype(float)).mean(axis=1)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--half-hops", type=float, nargs="+", default=[3.0, 5.0, 8.0])
    ap.add_argument("--seeds", type=int, nargs="+", default=[11, 22, 33])
    ap.add_argument("--n-neg-per-focus", type=int, default=4)
    args = ap.parse_args()
    t0 = time.time()
    defs = json.load(open(OUT / "programs/program_definitions.json"))
    proxy = json.load(open(ROOT / "infra/r04/marker_proxy_combined.json"))
    axis_full = {k: list(v["voted_genes"]) for k, v in proxy["classes"].items()
                 if k not in A.RETIRED_AXES}
    axis_full["Plasma"] = list(A.PLASMA_GENES)
    MHC2 = defs["P-mhc2"]["input_genes"]
    BGEN = axis_full["B"]

    man = json.load(open(ROOT / "infra/r04/role_manifests/external_validation_manifest.json"))
    rows_out = []
    for row in man["rows"]:
        sec = load_section_symbols(row)
        genes = [g if isinstance(g, str) else g.decode() for g in sec.gene_id]
        gidx = {}
        for i, g in enumerate(genes):
            if g not in gidx:
                gidx[g] = i
        counts = sec.counts.tocsr() if sparse.issparse(sec.counts) \
            else sparse.csr_matrix(sec.counts)
        coords = np.asarray(sec.coords, dtype=float)
        y = usz_tls_labels(sec.section_id, tuple(sec.barcode), ROOT)
        gt = np.where(np.isfinite(y), y, 0).astype(int)
        # hop单位
        d, _ = cKDTree(coords).query(coords, k=2)
        hop = float(np.median(d[:, 1]))
        # TLS foci（图连通分支，kNN图）
        tls = np.where(gt == 1)[0]
        w = C.spatial_weights(coords, C.KNN_SPATIAL)
        nc, lab = connected_components(w[tls][:, tls], directed=False)
        foci = [tls[lab == i] for i in range(nc)]
        ps_m = prog_score(counts, gidx, MHC2)
        ps_b = prog_score(counts, gidx, BGEN)
        sname = sec.section_id.split("::")[-1]
        for fi, fidx in enumerate(foci):
            c = coords[fidx].mean(axis=0)
            for hh in args.half_hops:
                h = hh * hop
                win = {"center": (float(c[0]), float(c[1])), "half": h}
                m = T.window_mask(coords, win)
                if m.sum() == 0:
                    continue
                cov = gt[m].mean()  # focus窗的TLS覆盖率（事后描述，非输入）
                for seed in args.seeds:
                    rng = np.random.default_rng(seed * 1000 + fi * 97 + int(hh * 13))
                    # 阴性窗：同大小随机窗，零覆盖才要
                    negs = []
                    tries = 0
                    while len(negs) < args.n_neg_per_focus and tries < 60:
                        tries += 1
                        rc = coords[rng.integers(len(coords))]
                        nw = {"center": (float(rc[0]), float(rc[1])), "half": h}
                        nm = T.window_mask(coords, nw)
                        if nm.sum() == 0 or gt[nm].mean() > 0:
                            continue
                        negs.append(nw)
                    hm_m = float(np.nanmean(ps_m[m]))
                    hb_m = float(np.nanmean(ps_b[m]))
                    hm_n = [float(np.nanmean(ps_m[T.window_mask(coords, nw)])) for nw in negs]
                    hb_n = [float(np.nanmean(ps_b[T.window_mask(coords, nw)])) for nw in negs]
                    rows_out.append({
                        "section": sec.section_id, "patient": sec.patient_id,
                        "focus": fi, "focus_size": len(fidx),
                        "half_hops": hh, "seed": seed,
                        "pos_cov": round(float(cov), 3),
                        "mhc2_pos": round(hm_m, 4), "b_pos": round(hb_m, 4),
                        "mhc2_neg_med": round(float(np.median(hm_n)), 4) if hm_n else "NA",
                        "b_neg_med": round(float(np.median(hb_n)), 4) if hb_n else "NA",
                        "mhc2_hit": bool(hm_m > np.median(hm_n)) if hm_n else None,
                        "b_hit": bool(hb_m > np.median(hb_n)) if hb_n else None,
                        "n_neg": len(negs)})
        print(f"  {sname}: {nc} foci done", flush=True)
    with open(OUT / "masked_structure_fine.tsv", "w") as f:
        cols = ["section", "patient", "focus", "focus_size", "half_hops", "seed",
                "pos_cov", "mhc2_pos", "b_pos", "mhc2_neg_med", "b_neg_med",
                "mhc2_hit", "b_hit", "n_neg"]
        f.write("\t".join(cols) + "\n")
        for r in rows_out:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"L-009 fine windows done ({time.time()-t0:.0f}s): {len(rows_out)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
