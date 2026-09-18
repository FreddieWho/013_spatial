#!/usr/bin/env python3
"""T3-2: 患者内共表达模块 — union 基因的患者等权相关 → 平均链接模块。

对 556 union 基因：每片 log1p（注意：共表达用 log1p counts，不用 z，
避免片内标准化吃掉幅度信息）→ 患者内 median 相关 → 30 患者平均 →
1-|r| 平均链接，切 0.5（|r|>=0.5 同模块）。输出模块表 + 患者支持。
流式：每片只读一次，只取 union 基因列。
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from r16 import census as C

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/r16/recovery_20260918/programs"


def main() -> int:
    cache = ROOT / "infra/r16/census_cache_hvg10k/sections"
    genes = [l.split("\t")[0] for l in
             open(OUT / "candidate_union.tsv").read().strip().split("\n")[1:]]
    print(f"union genes: {len(genes)}")
    stems = C.list_stems(cache)
    _, _, genes0, _, _ = C.load_section(cache, stems[0])
    gidx = {g: i for i, g in enumerate(genes0)}
    miss = [g for g in genes if g not in gidx]
    print(f"missing in cache: {len(miss)} {miss[:10]}")
    keep = [g for g in genes if g in gidx]
    cols = np.array([gidx[g] for g in keep])

    # 患者内 median 相关：先按患者聚合该患者各片的相关矩阵
    pat_mats: dict[str, list] = defaultdict(list)
    for si, stem in enumerate(stems):
        mat, _, _, _, meta = C.load_section(cache, stem)
        x = np.log1p(mat.toarray().astype(np.float32))[:, cols]
        # 过滤零方差列（该片不表达的基因）
        sd = x.std(axis=0)
        ok = sd > 1e-9
        r = np.full((len(keep), len(keep)), np.nan)
        if ok.sum() >= 2:
            cc = np.corrcoef(x[:, ok].T)
            r[np.ix_(ok, ok)] = cc
        pat_mats[meta["patient_id"]].append(r)
        if (si + 1) % 10 == 0:
            print(f"  [{si+1}/{len(stems)}]", flush=True)
    # 患者内 median → 跨患者 mean（患者等权；NaN 感知）
    pat_med = {p: np.nanmedian(np.stack(v), axis=0) for p, v in pat_mats.items()}
    R = np.nanmean(np.stack(list(pat_med.values())), axis=0)
    n_obs = np.sum([~np.isnan(np.stack(v)).all(axis=0) for v in pat_mats.values()], axis=0)
    np.savez_compressed(OUT / "coexpression_R.npz", R=R, genes=np.array(keep))
    # 平均链接 1-|r|，切 0.5；NaN→距离 1（无共表达证据）
    D = 1 - np.abs(np.nan_to_num(R, nan=0.0))
    np.fill_diagonal(D, 0.0)
    Z = linkage(squareform(D, checks=False), method="average")
    lab = fcluster(Z, t=0.5, criterion="distance")
    mods: dict[int, list[int]] = defaultdict(list)
    for i, g in enumerate(lab):
        mods[int(g)].append(i)
    with open(OUT / "coexpression_modules.tsv", "w") as fh:
        fh.write("module\tsize\tmean_abs_r\tmin_abs_r\tgenes\n")
        recs = []
        for m, mem in sorted(mods.items(), key=lambda kv: -len(kv[1])):
            sub = np.abs(R[np.ix_(mem, mem)])
            # ignore diagonal
            off = sub[~np.eye(len(mem), dtype=bool)]
            recs.append((m, len(mem), float(np.nanmean(off)) if len(off) else 1.0,
                         float(np.nanmin(off)) if len(off) else 1.0, mem))
        for m, sz, me, mi, mem in sorted(recs, key=lambda r: -r[1]):
            fh.write(f"M{m:02d}\t{sz}\t{me:.3f}\t{mi:.3f}\t"
                     f"{','.join(keep[i] for i in mem)}\n")
    print(f"modules: {len(mods)}; sizes: {sorted([len(v) for v in mods.values()], reverse=True)[:15]}")
    # 大模块预览
    for line in open(OUT / "coexpression_modules.tsv").read().strip().split("\n")[:8]:
        print("  " + line[:200])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
