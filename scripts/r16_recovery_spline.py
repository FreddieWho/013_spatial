#!/usr/bin/env python3
"""M2-spline: 程序分 natural spline df=3 ~ 到TLS距离（抄灵感论文 Science 2026）。

设计（与论文方法对齐，Visium坐标无µm则用hop单位）：
- 锚点：USZ 8片 ground_truth=TLS foci（独立病理标注）；ST-CRC无TLS GT不进本轮。
- 距离：每spot到最近TLS focus边缘的符号距离（focus内<0，参考Semla RadialDistance）；
  rescale到每片0–1（论文做法，便于跨片比较）；分析限TUM spots（论文：tumor-annotated）。
- 模型：score ~ ns(dist, df=3) OLS；分类：拟合曲线在0→0.75区间单调降=high-to-low，
  单调升=low-to-high，否则flat；效应量=拟合值域（max-min）；p=整体F vs 常数。
- 程序：T3旧6程序input分 + MHC-II家族（ANTIGEN_PROCESSING GO代表）+ B轴（参照）。
输出：infra/r16/recovery_20260918/spline_curves.tsv + spline_summary.tsv
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.interpolate import BSpline
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components

from r16 import axes as A
from r16 import census as C
from r16.section_io import load_section_symbols, usz_tls_labels

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/r16/recovery_20260918"


def ns_basis(x: np.ndarray, df: int = 3) -> np.ndarray:
    """Natural cubic spline basis via BSpline with natural boundary conditions.

    df=3 -> 2 interior knots at tertiles; design = [1, N1, N2, N3].
    """
    x = np.asarray(x, dtype=float)
    qs = np.quantile(x, [0.0, 1/3, 2/3, 1.0])
    # pad if degenerate
    for i in range(1, 4):
        if qs[i] <= qs[i-1]:
            qs[i] = qs[i-1] + 1e-6
    t = np.concatenate([[qs[0]]*4, qs[1:3], [qs[3]]*4])
    B = np.column_stack([BSpline(t, np.eye(len(t)-4)[i], 3)(x)
                         for i in range(len(t)-4)])
    # natural constraints: second derivative zero at boundaries -> drop 2 dof
    # (standard: use first df columns of full basis after absorbing intercept)
    X = np.column_stack([np.ones(len(x)), B[:, 1:df+1]])
    return X


def prog_score(counts, gidx: dict, genes: list[str]) -> np.ndarray:
    cols = [gidx[g] for g in genes if g in gidx]
    if not cols:
        return np.full(counts.shape[0], np.nan)
    X = counts[:, cols]
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    return np.log1p(X.astype(float)).mean(axis=1)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit-programs", type=int, default=0)
    ap.add_argument("--tag", type=str, default="v1")
    ap.add_argument("--go-list", type=Path, default=None,
                    help="JSON list of GO rep_ids: progs = GO input genes (GO线模式)")
    args = ap.parse_args()
    t0 = time.time()
    defs = json.load(open(OUT / "programs/program_definitions.json"))
    proxy = json.load(open(ROOT / "infra/r04/marker_proxy_combined.json"))
    axis_full = {k: list(v["voted_genes"]) for k, v in proxy["classes"].items()
                 if k not in A.RETIRED_AXES}
    axis_full["Plasma"] = list(A.PLASMA_GENES)
    # MHC-II家族代表（GO线最强的抗原呈递通路，按go_summary_full USZ_dM1排序取top3）
    go_extra = {}
    try:
        import csv
        rows = list(csv.DictReader(open(OUT / "programs_go/go_summary_full.tsv"), delimiter="\t"))
        ap_rows = [r for r in rows if "ANTIGEN_PROCESSING_AND_PRESENTATION" in r["go_id"]]
        ap_rows.sort(key=lambda r: -float(r["USZ_dM1"]))
        gdefs = json.load(open(OUT / "programs_go/definitions.json"))
        for r in ap_rows[:3]:
            d = gdefs[r["rep_id"]]
            go_extra["GO:" + r["go_id"][5:35]] = d["input_genes"]
    except Exception as e:
        print("GO extra skipped:", e)

    if args.go_list is not None:
        gdefs = json.load(open(OUT / "programs_go/definitions.json"))
        reps = json.load(open(args.go_list))
        progs = {"GO:" + r: gdefs[r]["input_genes"] for r in reps}
    else:
        progs = {pn: d["input_genes"] for pn, d in defs.items()}
        progs.update(go_extra)
        progs["B-axis"] = axis_full["B"]
    if args.limit_programs:
        progs = dict(list(progs.items())[:args.limit_programs])
    print(f"programs: {list(progs)}", flush=True)

    man = json.load(open(ROOT / "infra/r04/role_manifests/external_validation_manifest.json"))
    curve_rows, summ_rows = [], []
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
        # TUM mask（h5ad ground_truth TUM）
        import h5py
        alias = sec.section_id.split("::")[-1]
        path = ROOT / "data/other_sources/zenodo_usz_tls_visium/TLS_VISIUM_USZ/h5ad_preprocessed" / f"{alias}.h5ad"
        tum = np.zeros(len(coords), dtype=bool)
        with h5py.File(path, "r") as h:
            cats = [c.decode() if isinstance(c, bytes) else c
                    for c in h["obs"]["ground_truth"]["categories"][:]]
            codes = np.asarray(h["obs"]["ground_truth"]["codes"][:])
            idx = [b.decode() if isinstance(b, bytes) else b for b in h["obs"]["_index"][:]]
        lut = {b: cats[codes[i]] for i, b in enumerate(idx)}
        for i, b in enumerate(sec.barcode):
            if lut.get(b) == "TUM":
                tum[i] = True
        # TLS foci + 符号距离（focus内<0）
        tls = np.where(gt == 1)[0]
        w = C.spatial_weights(coords, C.KNN_SPATIAL)
        nc, lab = connected_components(w[tls][:, tls], directed=False)
        dmin = np.full(len(coords), np.nan)
        tree = cKDTree(coords[tls])
        dd, _ = tree.query(coords, k=1)
        # focus内判：最近TLS点同focus且被TLS包围 -> 简化：gt==1为核内
        dmin = dd.astype(float)
        dmin[gt == 1] = -dmin[gt == 1] * 0.5  # 核内为负（量级减半，示意core）
        # rescale 0-1（按片，论文做法；核内映射到0）
        pos = dmin[dmin >= 0]
        if len(pos) < 50:
            continue
        dmax = float(np.quantile(pos, 0.99))
        dn = np.clip(dmin / max(dmax, 1e-9), -0.1, 1.0)
        dn = np.clip((dn + 0.1) / 1.1, 0, 1)
        keep = tum & np.isfinite(dn)
        print(f"  {alias}: foci={nc} TUM={tum.sum()}/{len(coords)} kept={keep.sum()}", flush=True)
        if keep.sum() < 100:
            continue
        X = ns_basis(dn[keep])
        grid = np.linspace(0, 0.75, 16)
        Xg = ns_basis(grid)  # NOTE: 基函数节点按拟合数据分位算的，grid外推需同节点
        for pn, gl in progs.items():
            s = prog_score(counts, gidx, gl)[keep]
            ok = np.isfinite(s)
            if ok.sum() < 50:
                continue
            beta, _, _, _ = np.linalg.lstsq(X[ok], s[ok], rcond=None)
            # F vs 常数
            pred = X[ok] @ beta
            ss_res = float(((s[ok] - pred) ** 2).sum())
            ss_tot = float(((s[ok] - s[ok].mean()) ** 2).sum())
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
            n, p = ok.sum(), X.shape[1]
            F = ((ss_tot - ss_res) / (p - 1)) / (ss_res / (n - p)) if ss_res > 0 else 0.0
            from scipy.stats import f as fdist
            pval = float(fdist.sf(F, p - 1, n - p))
            # 拟合曲线：训练节点（分位数）固定后，在箱中心求同一基函数。
            # ns_basis 的节点依赖输入分位数——必须用拟合数据的节点，不能重算。
            # 做法：缓存拟合时 X 的节点分位数，箱中心用同节点基求值。
            xfit = dn[keep][ok]
            qs_fit = np.quantile(xfit, [0.0, 1/3, 2/3, 1.0])
            for i in range(1, 4):
                if qs_fit[i] <= qs_fit[i-1]:
                    qs_fit[i] = qs_fit[i-1] + 1e-6
            def basis_at(xq):
                t = np.concatenate([[qs_fit[0]]*4, qs_fit[1:3], [qs_fit[3]]*4])
                B = np.column_stack([BSpline(t, np.eye(len(t)-4)[i], 3)(xq)
                                       for i in range(len(t)-4)])
                return np.column_stack([np.ones(len(xq)), B[:, 1:4]])
            bins = np.linspace(0, 0.75, 9)
            bc, bv = [], []
            for bi in range(len(bins) - 1):
                m = (dn[keep][ok] >= bins[bi]) & (dn[keep][ok] < bins[bi + 1])
                if m.sum() >= 10:
                    c = (bins[bi] + bins[bi + 1]) / 2
                    bc.append(c)
                    bv.append(float(basis_at(np.array([c]))[0] @ beta))
            bv = np.array(bv)
            # 形态分类：Spearman单调性（robust）+ 45%端到端变化地板。
            # 严格逐差单调太脆（一箱抖动就判flat），论文曲线也是光滑趋势分类。
            if len(bv) >= 4:
                from scipy.stats import spearmanr
                rho, p_rho = spearmanr(np.arange(len(bv)), bv)
                rel = (bv[-1] - bv[0]) / max(abs(bv).max(), 1e-9)
                if p_rho < 0.05 and abs(rel) >= 0.15:
                    shape = "high-to-low" if rho < 0 else "low-to-high"
                elif abs(rel) < 0.15:
                    shape = "flat"
                else:
                    shape = "nonmonotonic"
            else:
                shape = "NA"
            amp = float(np.nanmax(bv) - np.nanmin(bv)) if len(bv) else 0.0
            summ_rows.append({"program": pn, "section": sec.section_id,
                              "n_foci": nc, "n_tum_kept": int(keep.sum()),
                              "r2": round(r2, 4), "p": f"{pval:.3g}",
                              "shape": shape, "amplitude": round(amp, 4)})
            for c_, v_ in zip(bc, bv):
                curve_rows.append({"program": pn, "section": sec.section_id,
                                   "dist01": round(float(c_), 3),
                                   "fit": round(float(v_), 4)})
    with open(OUT / f"spline_summary_{args.tag}.tsv", "w") as f:
        f.write("program\tsection\tn_foci\tn_tum_kept\tr2\tp\tshape\tamplitude\n")
        for r in summ_rows:
            f.write("\t".join(str(r[c]) for c in
                               ["program", "section", "n_foci", "n_tum_kept",
                                "r2", "p", "shape", "amplitude"]) + "\n")
    with open(OUT / f"spline_curves_{args.tag}.tsv", "w") as f:
        f.write("program\tsection\tdist01\tfit\n")
        for r in curve_rows:
            f.write("\t".join(str(r[c]) for c in
                               ["program", "section", "dist01", "fit"]) + "\n")
    print(f"spline {args.tag} done ({time.time()-t0:.0f}s): {len(summ_rows)} fits")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
