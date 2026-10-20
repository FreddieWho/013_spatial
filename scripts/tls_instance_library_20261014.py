"""74 张 HookNet 图的失败模式分层 + 候选实例库（标签前冻结的模型侧事实）。

标签到位后本库作为：评估协变量（focal/diffuse 分层）、标注分诊依据、
GC-TLS 嵌套结构清单。与标注无关，不降级。

分层规则（冻结）：
  diffuse  最大连通块面积 >= 0.79 mm²（TLS 罕见超过 ~1mm 直径，π·0.5²）
  focal    有 >=1 个 >=0.03 mm² 的块且非 diffuse
  sparse   其余（信号稀少/无）

输出：
  instance_library.jsonl  每实例一行（slide, channel, thr, area_mm2, 质心 fullres px）
  slide_stratification.tsv 每片一行（计数、focality、嵌套率、分层标签、fold/unit）
"""
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parent.parent
MAP_DIR = ROOT / "infra/hooknet_gpu_20261008/maps"
SIDECAR = ROOT / "infra/hooknet_gpu_20261008"
SECTIONS = ROOT / "infra/tls_four_cancer_labels_20261008/eligibility/sections.tsv"
FOLDS = ROOT / "infra/tls_hooknet_rna_triage_20261014/split_folds_v1.tsv"
OUT = ROOT / "infra/tls_instance_library_20261014"

UM = 4.0                      # mosaic µm/px
CELL_MM2 = UM * UM / 1e6
TLS_THRS = [0.3, 0.5, 0.7]
GC_THR = 0.5
MIN_AREA_MM2 = 0.005          # 聚集体下限
DIFFUSE_BLOB_MM2 = 0.79       # 单块超此即非 TLS 尺度
FOCAL_MIN_MM2 = 0.03          # "像 TLS" 的面积下限


def components(mask, thr_tag, channel, sid, ps_um_px):
    lab, n = ndimage.label(mask)
    out = []
    if n == 0:
        return out
    sizes = np.bincount(lab.ravel())
    for i in range(1, n + 1):
        nc = int(sizes[i])
        area = nc * CELL_MM2
        if area < MIN_AREA_MM2:
            continue
        ys, xs = np.nonzero(lab == i)
        out.append({
            "slide": sid, "channel": channel, "thr": thr_tag,
            "area_mm2": round(area, 4), "n_cells": nc,
            "cx_fullres_px": round(float(xs.mean()) * UM / ps_um_px, 1),
            "cy_fullres_px": round(float(ys.mean()) * UM / ps_um_px, 1),
            "bbox_cells": [int(ys.min()), int(xs.min()), int(ys.max()), int(xs.max())],
        })
    return out


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    with open(SECTIONS) as f:
        sec = {r["section_id"]: r for r in csv.DictReader(f, delimiter="\t")}
    folds = {}
    if FOLDS.exists():
        with open(FOLDS) as f:
            folds = {r["section_id"]: r for r in csv.DictReader(f, delimiter="\t")}

    inst_f = open(OUT / "instance_library.jsonl", "w")
    rows = []
    slides = sorted(p.stem for p in MAP_DIR.glob("*.npz"))
    for i, sid in enumerate(slides):
        z = np.load(MAP_DIR / f"{sid}.npz")
        prob = z["prob"].astype(np.float32) / 255.0
        ps = float(sec[sid]["pixel_size_um_estimated"])
        tls, gc = prob[:, :, 2], prob[:, :, 1]

        row = {"section_id": sid, "organ": sec[sid]["organ"],
               "patient_status": sec[sid]["patient_status"]}
        if sid in folds:
            row["unit_id"], row["fold"] = folds[sid]["unit_id"], folds[sid]["fold"]

        per_thr = {}
        for thr in TLS_THRS:
            insts = components(tls >= thr, thr, "tls", sid, ps)
            for e in insts:
                inst_f.write(json.dumps(e) + "\n")
            areas = sorted((e["area_mm2"] for e in insts), reverse=True)
            tot = sum(areas)
            per_thr[thr] = (len(insts), areas[0] if areas else 0.0, tot)
            row[f"n_tls_t{int(thr*10)}"] = len(insts)
            row[f"max_mm2_t{int(thr*10)}"] = round(areas[0], 4) if areas else 0.0
        gc_insts = components(gc >= GC_THR, GC_THR, "gc", sid, ps)
        for e in gc_insts:
            inst_f.write(json.dumps(e) + "\n")
        row["n_gc"] = len(gc_insts)

        # focality = 最大块占比；嵌套率 = gc 质心落在 tls(0.5) 膨胀邻域的比例
        n5, big5, tot5 = per_thr[0.5]
        row["focality_t5"] = round(big5 / tot5, 3) if tot5 > 0 else 0.0
        if gc_insts and n5:
            tls_mask = ndimage.binary_dilation(tls >= 0.5, iterations=10)
            inside = sum(1 for e in gc_insts
                         if tls_mask[min(int(e["cy_fullres_px"] * ps / UM), tls_mask.shape[0] - 1),
                                   min(int(e["cx_fullres_px"] * ps / UM), tls_mask.shape[1] - 1)])
            row["gc_nesting"] = round(inside / len(gc_insts), 3)
        else:
            row["gc_nesting"] = ""

        if big5 >= DIFFUSE_BLOB_MM2:
            row["stratum"] = "diffuse"
        elif big5 >= FOCAL_MIN_MM2:
            row["stratum"] = "focal"
        else:
            row["stratum"] = "sparse"
        rows.append(row)
        if (i + 1) % 20 == 0:
            print(f"[{i+1}/{len(slides)}] {time.time()-t0:.0f}s", flush=True)
    inst_f.close()

    keys = sorted({k for r in rows for k in r})
    with open(OUT / "slide_stratification.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    print("strata:", dict(Counter(r["stratum"] for r in rows)))
    print("by organ×stratum:")
    c = Counter((r["organ"], r["stratum"]) for r in rows)
    for k in sorted(c):
        print(f"  {k}: {c[k]}")
    print(f"elapsed {time.time()-t0:.0f}s; instances in library:",
          sum(1 for _ in open(OUT / 'instance_library.jsonl')))


if __name__ == "__main__":
    main()
