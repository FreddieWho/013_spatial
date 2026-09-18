#!/usr/bin/env python3
"""T1: PC matching 重算 — abs-cosine 分组 + 定向平均 + 新旧映射 + 联合解释力。

只重算 matching/registry 派生（02 T1-2）。不重训 PCA，不重跑 LISA。
输入：infra/r16/axis_factory_20260915.json + axis_factory_loadings_20260915.npz
输出：infra/r16/recovery_20260918/audit/{pc_matching_reanalysis.json,
      old_to_new_groups.tsv, subspace_diagnostic.json, geometry_diagnostic.json}
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.linalg import qr, svd

from r16 import axes as A
from r16 import axis_factory as F
from r16 import census as C
from r16.recovery import pc_matching as M

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/r16/recovery_20260918/audit"


def joint_r2(y: np.ndarray, B: np.ndarray) -> float:
    """ loading y 被已知轴 loading 矩阵 B 联合解释的比例（A03 的正确新奇性）。"""
    Q, _ = qr(B, mode="economic")
    return float((Q @ (Q.T @ y)) @ (Q @ (Q.T @ y)) / max(y @ y, 1e-12))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--artifact", type=Path,
                    default=ROOT / "infra/r16/axis_factory_20260915.json")
    ap.add_argument("--loadings", type=Path,
                    default=ROOT / "infra/r16/axis_factory_loadings_20260915.npz")
    ap.add_argument("--marker-proxy", type=Path,
                    default=ROOT / "infra/r04/marker_proxy_combined.json")
    ap.add_argument("--union-genes", type=Path, default=ROOT / "infra/r16/hvg_union_genes.json")
    ap.add_argument("--rank-file", type=Path,
                    default=ROOT / "infra/r16/hvg_final_rank_20260914.npy")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    art = json.load(open(args.artifact))
    z = np.load(args.loadings)
    L = z["loadings"].astype(np.float64)
    genes = list(z["genes"])
    # 旧分组（复现 legacy，供映射）：signed cosine
    nrm = L / np.linalg.norm(L, axis=1, keepdims=True)
    sim_old = nrm @ nrm.T
    np.fill_diagonal(sim_old, 1.0)
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    Z_old = linkage(squareform(1 - np.clip(sim_old, -1, 1), checks=False), method="average")
    lab_old = fcluster(Z_old, t=1 - C.MATCH_COSINE, criterion="distance")
    # 新分组：|cosine|
    lab_new = M.match_groups(L, cutoff=C.MATCH_COSINE)

    cands = art["candidates"]
    assert len(cands) == len(L) == len(lab_old) == len(lab_new)
    old_groups: dict[int, list[int]] = {}
    new_groups: dict[int, list[int]] = {}
    for i, g in enumerate(lab_old):
        old_groups.setdefault(int(g), []).append(i)
    for i, g in enumerate(lab_new):
        new_groups.setdefault(int(g), []).append(i)

    # 旧->新映射：旧组内最大份额的新组
    from collections import Counter
    old2new = {}
    for g, mem in old_groups.items():
        c = Counter(int(lab_new[i]) for i in mem)
        top, n = c.most_common(1)[0]
        old2new[g] = {"new_group": top, "frac": n / len(mem), "size": len(mem)}
    with open(OUT / "old_to_new_groups.tsv", "w") as f:
        f.write("old_group\tsize\tnew_group\tfrac\n")
        for g in sorted(old2new):
            r = old2new[g]
            f.write(f"{g}\t{r['size']}\t{r['new_group']}\t{r['frac']:.3f}\n")

    # 新组记录：定向平均 + conflict + 联合解释力
    # 已知轴 loading：在同一基因坐标下用轴 marker 基因的单位向量近似
    gidx = {g: i for i, g in enumerate(genes)}
    proxy = json.load(open(args.marker_proxy))
    union = json.load(open(args.union_genes))
    rank = np.load(args.rank_file)
    top = {union[i] for i in np.flatnonzero(rank < 10000)}
    axis_defs = {k: [g for g in v["voted_genes"] if g in top]
                 for k, v in proxy["classes"].items() if k not in A.RETIRED_AXES}
    axis_defs["Plasma"] = [g for g in A.PLASMA_GENES if g in top]
    Bcols = []
    for an in sorted(axis_defs):
        v = np.zeros(len(genes))
        cols = [gidx[g] for g in axis_defs[an] if g in gidx]
        if cols:
            v[cols] = 1.0 / np.sqrt(len(cols))
            Bcols.append(v)
    B = np.stack(Bcols, axis=1)

    new_recs = []
    split_count = sum(1 for r in old2new.values() if r["frac"] < 1.0)
    merge_info = []
    new2old: dict[int, Counter] = {}
    for i, g in enumerate(lab_new):
        new2old.setdefault(int(g), Counter())[int(lab_old[i])] += 1
    for g, c in new2old.items():
        if len(c) > 1:
            merge_info.append({"new_group": g, "n_old_groups": len(c),
                               "old_groups": sorted(c)})
    for g in sorted(new_groups):
        mem = new_groups[g]
        mean, flips, conflict = M.oriented_mean(L, mem)
        markers = F.top_loading_genes(mean, genes, 10)
        r2 = joint_r2(mean, B)
        # 旧 novelty：成员单轴 max|r|>=0.7 比例；新：联合 R²
        pats = sorted({cands[i]["patient_id"] for i in mem})
        rec_pats = sorted({cands[i]["patient_id"] for i in mem if cands[i]["tile_ok"]})
        new_recs.append({"new_group": g, "n_members": len(mem),
                         "n_patients": len(pats), "recurrence_patients": len(rec_pats),
                         "sign_conflict": conflict, "n_flipped": flips.count(-1),
                         "joint_axis_R2": round(r2, 4), "top_markers": markers,
                         "merged_from_old_groups": len(new2old[g])})

    # 复现判定（与 legacy 同规则，但用新组）：purity 用 |cos| min，>=0.5 + rec>=2
    sim_new = M.abs_cosine_matrix(L)
    n_rep_new = n_new_rep = 0
    for rec in new_recs:
        mem = new_groups[rec["new_group"]]
        sub = sim_new[np.ix_(mem, mem)].copy()
        np.fill_diagonal(sub, 1.0)
        purity = float(sub.min()) if len(mem) > 1 else 1.0
        rec["abs_purity"] = round(purity, 4)
        grade = ("EXPLORATORY_REPRODUCED"
                 if (purity >= 0.5 and rec["recurrence_patients"] >= 2)
                 else "DESCRIPTIVE_SINGLE_PATIENT")
        rec["grade_new_rule"] = grade
        known_by_joint = rec["joint_axis_R2"] >= 0.49  # ~= 单轴0.7²，联合口径
        rec["known_by_joint_R2"] = bool(known_by_joint)
        if grade == "EXPLORATORY_REPRODUCED":
            n_rep_new += 1
            if not known_by_joint:
                n_new_rep += 1

    reanalysis = {
        "schema": "013_spatial.recovery_T1.v1",
        "n_candidates": len(L),
        "n_old_groups": len(old_groups), "n_new_groups": len(new_groups),
        "old_groups_split": split_count,
        "new_groups_merged_from_multiple_old": merge_info,
        "legacy": {"n_groups": art["n_groups"],
                   "reproduced": sum(1 for g in art["groups"]
                                     if g["evidence_grade"] == "EXPLORATORY_REPRODUCED"),
                   "new_reproduced": sum(1 for g in art["groups"]
                                         if g["novelty"] == "NEW" and g["evidence_grade"] == "EXPLORATORY_REPRODUCED")},
        "corrected": {"n_groups": len(new_groups), "reproduced": n_rep_new,
                      "new_reproduced_joint_R2": n_new_rep},
        "groups": new_recs,
        "elapsed_s": round(time.time() - t0, 1),
    }
    (OUT / "pc_matching_reanalysis.json").write_text(
        json.dumps(reanalysis, ensure_ascii=False, indent=2))
    print(f"old_groups={len(old_groups)} new_groups={len(new_groups)} "
          f"split={split_count} merged={len(merge_info)} "
          f"rep: legacy->corrected, new_rep_joint={n_new_rep} "
          f"({time.time()-t0:.0f}s)")
    for m in merge_info[:10]:
        print("  merged:", m)
    top_new = sorted([r for r in new_recs if r["grade_new_rule"] == "EXPLORATORY_REPRODUCED"
                      and not r["known_by_joint_R2"]],
                     key=lambda r: -r["recurrence_patients"])[:10]
    for r in top_new:
        print(f"  NEWcandidate g{r['new_group']}: nmem={r['n_members']} "
              f"rec={r['recurrence_patients']} R2={r['joint_axis_R2']} "
              f"conflict={r['sign_conflict']} {r['top_markers'][:6]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
