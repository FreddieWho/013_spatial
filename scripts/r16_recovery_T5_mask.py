#!/usr/bin/env python3
"""T5: 遮蔽任务 — A. 场状态预测（全程序） + B. USZ-TLS 隐藏核心检测（03 §7）。

A（场预测）：ST-CRC + USZ，每片规则窗遮蔽，预测冻结 readout 均值；
  基线：M0(Q+C) / 组成+marker / 最近邻插值；完整模型=M1/M2特征。
  终点：MSE + 分程序分队列留一患者汇总。主张上限=场预测。
B（结构定位）：USZ 8 片，GT=ground_truth=TLS（独立病理标注，非表达派生）；
  通用窗事后套 GT 定阴阳；输出存在概率（窗内 TLS 覆盖率预测）+ 位置质量；
  基线含 B 轴（已知最强 TLS 组成代理）——程序必须超 B 轴才算增量。
ST-CRC IC aggregate 不进 B（D-033：非已验证 TLS）。
输出：mask_task_manifest.tsv, leakage_tests.json,
      masked_field_results.tsv, masked_structure_results.tsv
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse

from r16 import axes as A
from r16 import census as C
from r16.recovery import mask_tasks as T
from r16.section_io import load_section_symbols

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
    ap.add_argument("--n-windows", type=int, default=12)
    ap.add_argument("--window-size", type=int, default=8)
    ap.add_argument("--seed", type=int, default=C.SEED)
    args = ap.parse_args()
    t0 = time.time()
    defs = json.load(open(OUT / "programs/program_definitions.json"))
    proxy = json.load(open(ROOT / "infra/r04/marker_proxy_combined.json"))
    axis_full = {k: list(v["voted_genes"]) for k, v in proxy["classes"].items()
                 if k not in A.RETIRED_AXES}
    axis_full["Plasma"] = list(A.PLASMA_GENES)

    manifest, leak, field_rows, struct_rows = [], [], [], []
    # ---- B 前置：USZ TLS GT 读取（section_io 同口径） ----
    from r16.section_io import usz_tls_labels

    for mf, cohort in [("infra/r04/role_manifests/internal_validation_manifest.json", "ST-CRC"),
                       ("infra/r04/role_manifests/external_validation_manifest.json", "USZ")]:
        man = json.load(open(ROOT / mf))
        for row in man["rows"]:
            sec = load_section_symbols(row)
            genes = [g if isinstance(g, str) else g.decode() for g in sec.gene_id]
            gidx = {}
            for i, g in enumerate(genes):
                if g not in gidx:
                    gidx[g] = i
            counts = sec.counts.tocsr() if sparse.issparse(sec.counts) \
                else sparse.csr_matrix(sec.counts)
            lib = np.asarray(counts.sum(axis=1)).ravel()
            ng = np.asarray((counts > 0).sum(axis=1)).ravel()
            coords = np.asarray(sec.coords, dtype=float)
            wins = T.generate_windows(coords, args.n_windows, args.window_size,
                                      args.seed)
            gt = None
            if cohort == "USZ":
                y = usz_tls_labels(sec.section_id, tuple(sec.barcode), ROOT)
                if y is not None:
                    gt = np.where(np.isfinite(y), y, 0).astype(int)
            manifest.append({"cohort": cohort, "section": sec.section_id,
                             "patient": sec.patient_id, "n_spots": len(coords),
                             "n_windows": len(wins), "has_struct_gt": gt is not None})
            Xd = counts.toarray().astype(np.float32)
            np.log1p(Xd, out=Xd)
            mu, sd = Xd.mean(0), Xd.std(0)
            sd[sd < 1e-8] = 1.0
            Xz = (Xd - mu) / sd
            C6 = np.column_stack([A.axis_scores(Xz, gidx, gs)[0][:, None]
                                  for gs in axis_full.values()])
            Q = np.column_stack([np.ones(len(coords)), np.log10(lib + 1),
                                 (ng - ng.mean()) / (ng.std() + 1e-9)])
            # ---- A. 场预测：每窗遮蔽，预测 readout；基线 NN 插值 vs M1 vs M2 ----
            for pn in defs:
                inp = prog_score(counts, gidx, defs[pn]["input_genes"])
                ro = prog_score(counts, gidx, defs[pn]["readout_genes"])
                errs = {"nn": [], "m1": [], "m2": []}
                for wi, w in enumerate(wins):
                    m = T.window_mask(coords, w)
                    if m.sum() == 0 or m.sum() == len(coords):
                        continue
                    vis = ~m
                    # NN插值基线：可见点 readout 的最近邻
                    from scipy.spatial import cKDTree
                    tree = cKDTree(coords[vis])
                    _, idx = tree.query(coords[m], k=1)
                    pred_nn = ro[vis][idx]
                    # M1/M2：可见点拟合（Q+C+inp / +邻域），系数窗内估计=局部校准口径
                    D0 = np.column_stack([Q[vis], C6[vis]])
                    yv = ro[vis]
                    b0, _, _, _ = np.linalg.lstsq(
                        np.column_stack([D0, inp[vis]]), yv, rcond=None)
                    # M2 邻域项仅可见点消息
                    feat, _ = T.visible_features(coords, inp, w)
                    b2, _, _, _ = np.linalg.lstsq(
                        np.column_stack([D0, inp[vis], feat[vis]]), yv, rcond=None)
                    pred_m1 = np.column_stack([Q[m], C6[m], inp[m]]) @ b0 \
                        if len(b0) == D0.shape[1] + 1 else np.full(m.sum(), np.nan)
                    row2 = np.column_stack([Q[m], C6[m], inp[m], feat[m]])
                    pred_m2 = row2 @ b2 if len(b2) == row2.shape[1] \
                        else np.full(m.sum(), np.nan)
                    for k, pr in [("nn", pred_nn), ("m1", pred_m1), ("m2", pred_m2)]:
                        errs[k].append(float(((pr - ro[m]) ** 2).mean()))
                for k in errs:
                    field_rows.append({"program": pn, "cohort": cohort,
                                       "section": sec.section_id,
                                       "patient": sec.patient_id,
                                       "baseline": k,
                                       "mse_med": round(float(np.median(errs[k])), 5)
                                       if errs[k] else "NA",
                                       "n_windows": len(errs[k])})
            # ---- 泄漏探针：隐藏值改任意数，输入特征不变 ----
            probe = {"section": sec.section_id, "n_windows": len(wins)}
            w0 = wins[0]
            f1, _ = T.visible_features(coords, prog_score(
                counts, gidx, defs["P-epi"]["input_genes"]), w0)
            Xd2 = counts.toarray().astype(float)
            m0 = T.window_mask(coords, w0)
            Xd2[m0] = 999.0
            f2, _ = T.visible_features(coords, np.log1p(Xd2).mean(axis=1), w0) \
                if False else (None, None)
            # 直接探针：同一函数，隐藏表达改值，特征必须逐位相等
            e1 = prog_score(counts, gidx, defs["P-epi"]["input_genes"])
            e2 = e1.copy()
            e2[m0] = 999.0
            g1, _ = T.visible_features(coords, e1, w0)
            g2, _ = T.visible_features(coords, e2, w0)
            probe["leak_free"] = bool(np.array_equal(g1, g2, equal_nan=True))
            leak.append(probe)
            # ---- B. USZ-TLS 隐藏核心检测 ----
            if gt is not None:
                pos, neg = T.label_windows(wins, coords, gt)
                # 存在概率：窗内 TLS 覆盖率的程序分预测（可见点拟合 logistic 近似→线性分）
                bscore = prog_score(counts, gidx, axis_full["B"])
                for pn in list(defs) + ["B-axis"]:
                    ps = prog_score(counts, gidx,
                                    defs[pn]["input_genes"] if pn in defs else axis_full["B"])
                    aucs = {}
                    for wi in pos + neg:
                        m = T.window_mask(coords, wins[wi])
                        vis = ~m
                        # 可见点上：TLS标签~程序分 的 AUC（训练可见，测隐藏窗均值分）
                        from r04.structure_readout import _auc
                        yv = gt[vis]
                        pv = ps[vis]
                        fin = np.isfinite(yv) & np.isfinite(pv)
                        if fin.sum() < 20 or len(np.unique(yv[fin])) < 2:
                            continue
                        aucs[wi] = float(_auc(yv[fin], pv[fin],
                                             np.ones(fin.sum())))
                    # 窗排序质量：阳性窗的隐藏窗均值分是否高于阴性窗
                    hmeans = {}
                    for wi in pos + neg:
                        m = T.window_mask(coords, wins[wi])
                        hmeans[wi] = float(np.nanmean(ps[m]))
                    from r04.structure_readout import _auc as auc_fn
                    yy = np.array([1 if i in pos else 0 for i in hmeans])
                    pp = np.array([hmeans[i] for i in hmeans])
                    fin = np.isfinite(pp)
                    auc_h = float(auc_fn(yy[fin], pp[fin], np.ones(fin.sum()))) \
                        if fin.sum() >= 4 and len(np.unique(yy[fin])) == 2 else None
                    struct_rows.append({"program": pn, "section": sec.section_id,
                                        "patient": sec.patient_id,
                                        "n_pos_windows": len(pos), "n_neg_windows": len(neg),
                                        "hidden_mean_auc": round(auc_h, 4)
                                        if auc_h is not None else "NA"})
    with open(OUT / "mask_task_manifest.tsv", "w") as f:
        f.write("cohort\tsection\tpatient\tn_spots\tn_windows\thas_struct_gt\n")
        for r in manifest:
            f.write(f"{r['cohort']}\t{r['section']}\t{r['patient']}\t"
                    f"{r['n_spots']}\t{r['n_windows']}\t{r['has_struct_gt']}\n")
    json.dump(leak, open(OUT / "leakage_tests.json", "w"), indent=2)
    with open(OUT / "masked_field_results.tsv", "w") as f:
        f.write("program\tcohort\tsection\tpatient\tbaseline\tmse_med\tn_windows\n")
        for r in field_rows:
            f.write("\t".join(str(r[c]) for c in
                               ["program", "cohort", "section", "patient",
                                "baseline", "mse_med", "n_windows"]) + "\n")
    with open(OUT / "masked_structure_results.tsv", "w") as f:
        f.write("program\tsection\tpatient\tn_pos_windows\tn_neg_windows\thidden_mean_auc\n")
        for r in struct_rows:
            f.write("\t".join(str(r[c]) for c in
                               ["program", "section", "patient",
                                "n_pos_windows", "n_neg_windows",
                                "hidden_mean_auc"]) + "\n")
    print(f"T5 done ({time.time()-t0:.0f}s): field={len(field_rows)} "
          f"struct={len(struct_rows)} leak_free={all(r['leak_free'] for r in leak)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
