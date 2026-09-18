#!/usr/bin/env python3
"""T2: Lane A 修正版小规模诊断（02 T2，不重跑 60000 组合）。

修正点（A05/A06/A07/A09）：
- 患者先汇总（median_s），再 cohort median；观察与每次 null 同一算子。
- 联合 draw：每次 surrogate 对全队列产生一个统计量；B=联合重复次数。
- remap 全部 info 保留：accepted/frac/唯一源比例/重复率/Moran 改变。
- 主指标：固定对比（区内-区外）+ max 峰高（诊断）；峰位仅辅助。
- 诊断集：6 患者 x 先一片（按组织大小/形状/覆盖选，不按效应挑）；
  基因：组成正对照（CXCL13~B）+ LaneB 191 代表 + 模拟阴性/注入信号。
输出：infra/r16/recovery_20260918/laneA/{laneA_corrected_small.json,
      surrogate_diagnostics.tsv, calibration_curves.tsv, legacy_vs_corrected.md}
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from r16 import axes as A
from r16 import census as C
from r16.recovery import laneA_stats as S

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/r16/recovery_20260918/laneA"
PC_BINS = np.arange(-4, 9)

# 6 诊断患者：组织大小/形状/覆盖 spanning，不按效应挑选
# （全 ID，T0 时简写未对上 manifest，已修正）
DIAG_PATIENTS = [
    "HTAN_VANDERBILT_CRC::PATIENT::4aa571d63ab1aa20",  # 4 切片，大组织
    "HTAN_VANDERBILT_CRC::PATIENT::4f81822a6f2ecd63",  # 3 切片
    "HTAN_VANDERBILT_CRC::PATIENT::0bd3490036cbd4de",  # 2 切片
    "HTAN_VANDERBILT_CRC::PATIENT::286667cab598bb4a",  # 2 切片
    "HTAN_VANDERBILT_CRC::PATIENT::0def2e0b01043545",  # 2 切片
    "HTAN_VANDERBILT_CRC::PATIENT::07dc676d92f9d153",  # 小组织 665 spots
]
DIAG_GENES = ["CXCL13", "EPCAM", "MUC2", "COL1A1", "FN1", "FOS",
              "MKI67", "PTPRC", "SYNTH_NEG", "SYNTH_BUMP"]


def binned_profile(values: np.ndarray, signed: np.ndarray):
    m = np.full(len(PC_BINS), np.nan)
    ns = np.zeros(len(PC_BINS), dtype=int)
    for bi, b in enumerate(PC_BINS):
        sel = signed == b
        if sel.any():
            m[bi] = values[sel].mean()
            ns[bi] = int(sel.sum())
    return m, ns


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-draws", type=int, default=199)
    ap.add_argument("--seed", type=int, default=C.SEED)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    cache = ROOT / "infra/r16/census_cache_hvg10k/sections"
    man = json.load(open(ROOT / "infra/r16/census_cache_hvg10k/cache_manifest.json"))
    # 每患者先一片：该患者第一片（manifest 顺序，效应盲选）
    pat_first = {}
    for r in man["rows"]:
        pat_first.setdefault(r["patient_id"], r["stem"])
    stems = [pat_first[p] for p in DIAG_PATIENTS]
    print("diag stems:", [s[:12] for s in stems], flush=True)

    proxy = json.load(open(ROOT / "infra/r04/marker_proxy_combined.json"))
    union = json.load(open(ROOT / "infra/r16/hvg_union_genes.json"))
    rank = np.load(ROOT / "infra/r16/hvg_final_rank_20260914.npy")
    top10k = {union[i] for i in np.flatnonzero(rank < 10000)}
    axis_defs = {k: [g for g in v["voted_genes"] if g in top10k]
                 for k, v in proxy["classes"].items() if k not in A.RETIRED_AXES}
    axis_defs["Plasma"] = [g for g in A.PLASMA_GENES if g in top10k]

    rng = np.random.default_rng(args.seed)
    # per (axis, gene): by_patient profiles + section ns
    data = {}
    diag_rows = []
    for si, stem in enumerate(stems):
        mat, bcs, genes, coords, meta = C.load_section(cache, stem)
        gidx = {g: i for i, g in enumerate(genes)}
        x_z = C.log1p_zscore(mat)
        # depth proxy: 源 counts 行和（cache 存 counts，可直接算）
        libsize = np.asarray(mat.sum(axis=1)).ravel()
        w = C.spatial_weights(coords, C.KNN_SPATIAL)
        sB, _ = A.axis_scores(x_z, gidx, axis_defs["B"])
        mB = A.contour_mask(sB, 0.7)
        signed = A.signed_hop_distance(w, mB)
        # 合成对照：阴性（噪声）与注入 bump（区内 +1.5z）
        synth_neg = rng.normal(0, 1, len(bcs))
        synth_bump = rng.normal(0, 1, len(bcs)) + np.where(signed <= -1, 1.5, 0.0)
        X = {"SYNTH_NEG": synth_neg, "SYNTH_BUMP": synth_bump}
        for g in DIAG_GENES:
            if g in gidx:
                X[g] = x_z[:, gidx[g]]
        for g, vals in X.items():
            prof, ns = binned_profile(vals, signed)
            data.setdefault(("B", g), {})[meta["patient_id"]] = [prof]
            diag_rows.append({"stem": stem, "patient": meta["patient_id"][:18],
                              "gene": g, "n_spots": len(bcs),
                              "libsize_med": round(float(np.median(libsize)), 1)})
        print(f"  [{si+1}/{len(stems)}] {stem[:24]} n={len(bcs)}", flush=True)

    # remap 诊断：每片每次 draw 的 info + 唯一源比例/重复率/Moran 改变
    sur_rows = []
    for si, stem in enumerate(stems):
        mat, bcs, genes, coords, meta = C.load_section(cache, stem)
        gidx = {g: i for i, g in enumerate(genes)}
        x_z = C.log1p_zscore(mat)
        w = C.spatial_weights(coords, C.KNN_SPATIAL)
        cx = x_z[:, gidx["CXCL13"]] if "CXCL13" in gidx else x_z[:, 0]
        moran_obs = C.moran_i(cx, w)
        for d in range(args.n_draws):
            perm, info = C.spatial_null_remap(coords, seed=args.seed + si * 1000 + d)
            uniq = len(np.unique(perm)) / len(perm)
            # 重复率：被多个目标复用的源点比例
            rep = 1.0 - uniq
            moran_sur = C.moran_i(cx[perm], w)
            sur_rows.append({"stem": stem[:24], "draw": d,
                             "accepted": info["accepted"],
                             "frac_within_tol": round(info["frac_within_tol"], 3),
                             "unique_src_frac": round(uniq, 3),
                             "dup_rate": round(rep, 3),
                             "moran_obs": round(moran_obs, 3),
                             "moran_sur": round(moran_sur, 3),
                             "moran_ratio": round(moran_sur / moran_obs, 3)
                             if moran_obs > 1e-6 else None})
    with open(OUT / "surrogate_diagnostics.tsv", "w") as f:
        cols = list(sur_rows[0])
        f.write("\t".join(cols) + "\n")
        for r in sur_rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    import statistics
    acc = sum(1 for r in sur_rows if r["accepted"]) / len(sur_rows)
    print(f"surrogate: {len(sur_rows)} draws, accept_rate={acc:.3f}, "
          f"uniq_med={statistics.median(r['unique_src_frac'] for r in sur_rows):.3f}, "
          f"moran_ratio_med={statistics.median(r['moran_ratio'] for r in sur_rows if r['moran_ratio']):.3f}")

    # 修正统计量：每 (axis,gene) 联合 null
    results = []
    for (a, g), by_pat in data.items():
        pp = S.summarize_patient_profiles(by_pat)
        t_max = S.cohort_max_statistic(pp)
        t_con = S.cohort_contrast_statistic(pp)
        # 联合 surrogate：每次 draw 全队列同一机制
        null_max, null_con = [], []
        for d in range(args.n_draws):
            pp_null = {}
            for p, secs in by_pat.items():
                # surrogate: 该患者的 profile 在 bin 间置换（保边际，破距离关联）
                # + 注：这是诊断集的占位关联 null；主方法变异函数匹配见 T2 后续
                pp_null[p] = [rng.permutation(s) for s in secs]
            ppn = S.summarize_patient_profiles(pp_null)
            null_max.append(S.cohort_max_statistic(ppn))
            null_con.append(abs(S.cohort_contrast_statistic(ppn)))
        p_max = S.matched_null_p(t_max, null_max)
        p_con = S.matched_null_p(abs(t_con), null_con)
        results.append({"axis": a, "gene": g, "n_patients": len(by_pat),
                        "t_max": round(t_max, 3), "p_max": round(p_max, 4),
                        "t_contrast": round(t_con, 3), "p_contrast": round(p_con, 4),
                        "null_draws": args.n_draws,
                        "null_kind": "PLACEHOLDER_bin_permutation_pending_variogram"})
    (OUT / "laneA_corrected_small.json").write_text(
        json.dumps({"schema": "013_spatial.recovery_T2_small.v1",
                    "diag_patients": 6, "diag_sections": 6,
                    "null_draws": args.n_draws, "seed": args.seed,
                    "results": results}, ensure_ascii=False, indent=2))
    # 校准曲线：合成阴性/注入在两种统计量下的 p
    with open(OUT / "calibration_curves.tsv", "w") as f:
        f.write("gene\tt_max\tp_max\tt_contrast\tp_contrast\texpectation\n")
        for r in results:
            exp = ("null_uniform" if r["gene"] == "SYNTH_NEG"
                   else ("small_p" if r["gene"] == "SYNTH_BUMP" else "unknown"))
            f.write(f"{r['gene']}\t{r['t_max']}\t{r['p_max']}\t"
                    f"{r['t_contrast']}\t{r['p_contrast']}\t{exp}\n")
    for r in sorted(results, key=lambda x: x["p_contrast"]):
        print(f"  {r['gene']:10s} tmax={r['t_max']:.2f} pmax={r['p_max']:.3f} "
              f"tcon={r['t_contrast']:+.2f} pcon={r['p_contrast']:.3f}")
    print(f"done ({time.time()-t0:.0f}s) -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
