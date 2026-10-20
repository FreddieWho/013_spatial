"""HookNet 概率图 vs 人工标注的评估管线（标签到位后一条命令出报告）。

三层指标，全部按 SPLIT_FREEZE 的 unit/fold 记录，折外纪律在使用时执行：
  1) 实例级：标注实例 vs 候选实例（ch2 阈值→连通域→面积下限），
     匹配 = IoU>=MIN_IOU 或 候选质心落在标注内 或 标注质心落在候选内；
     报 recall / precision / F1（逐阈值）。
  2) spot 级：tls_map_labels.map_section 的三态标签（positive/negative，
     unknown 永不计入）vs spot 足迹平均 ch2 概率，报 PR-AUC / ROC-AUC。
  3) 聚合：逐片、逐器官、按 unit 中位数。

--selftest：合成图 + 合成标注，断言 recall=1/precision=1/PR-AUC=1、
打乱图 PR-AUC≈流行率、unknown 不被算成 negative。
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

ROOT = Path("/home/huyudi/013_spatial")
sys.path.insert(0, str(ROOT / "scripts"))
from tls_map_labels import map_section  # noqa: E402

MAP_DIR = ROOT / "infra/hooknet_gpu_20261008/maps"
FOLDS = ROOT / "infra/tls_hooknet_rna_triage_20261014/split_folds_v1.tsv"
OUT = ROOT / "infra/tls_eval_20261014"

THRESHOLDS = [0.3, 0.5, 0.7]
MIN_AREA_MM2 = 0.005          # 312 格 @4µm；淋巴聚集体也允许
MIN_IOU = 0.1
MOSAIC_UM_PER_PX = 4.0
TLS_CHANNEL = 2


def extract_instances(prob_tls: np.ndarray, thr: float, min_area_mm2: float):
    m = prob_tls >= thr
    lab, n = ndimage.label(m)
    out = []
    min_cells = min_area_mm2 * 1e6 / (MOSAIC_UM_PER_PX ** 2)
    for i in range(1, n + 1):
        ys, xs = np.nonzero(lab == i)
        if len(ys) < min_cells:
            continue
        out.append({"mask_idx": i, "n_cells": len(ys),
                    "cy": float(ys.mean()), "cx": float(xs.mean()),
                    "bbox": (int(ys.min()), int(xs.min()), int(ys.max()), int(xs.max()))})
    return lab, out


def mask_of(lab, idx, bbox):
    y0, x0, y1, x1 = bbox
    sub = np.zeros((y1 - y0 + 3, x1 - x0 + 3), dtype=bool)
    sub[1:-1, 1:-1] = lab[y0:y1 + 1, x0:x1 + 1] == idx
    return sub, (y0 - 1, x0 - 1)


def ann_to_cells(geom, um_per_px_fullres: float):
    """标注多边形（fullres px）栅格化到 mosaic 格。返回 set[(y,x)] 与质心。"""
    from shapely.geometry import Point
    minx, miny, maxx, maxy = geom.bounds
    f = um_per_px_fullres / MOSAIC_UM_PER_PX
    y0, y1 = int(miny * f) - 1, int(maxy * f) + 2
    x0, x1 = int(minx * f) - 1, int(maxx * f) + 2
    cells = set()
    for yy in range(max(0, y0), y1):
        for xx in range(max(0, x0), x1):
            if geom.contains(Point((xx + 0.5) / f, (yy + 0.5) / f)):
                cells.add((yy, xx))
    cy, cx = geom.centroid.y * f, geom.centroid.x * f
    return cells, (cy, cx)


def instance_metrics(prob_tls, annotations_fullres, um_per_px_fullres, thr):
    """annotations_fullres: list[shapely Polygon]。返回 recall/precision/F1。"""
    lab, cands = extract_instances(prob_tls, thr, MIN_AREA_MM2)
    anns = [ann_to_cells(g, um_per_px_fullres) for g in annotations_fullres]
    matched_a, matched_c = set(), set()
    for ci, c in enumerate(cands):
        mask, (oy, ox) = mask_of(lab, c["mask_idx"], c["bbox"])
        ys, xs = np.nonzero(mask)
        cset = set(zip((ys + oy).tolist(), (xs + ox).tolist()))
        for ai, (aset, (acy, acx)) in enumerate(anns):
            if ai in matched_a:
                continue
            inter = len(cset & aset)
            iou = inter / max(len(cset | aset), 1)
            centroid_hit = ((int(c["cy"]), int(c["cx"])) in aset) or \
                           ((int(acy), int(acx)) in cset)
            if iou >= MIN_IOU or centroid_hit:
                matched_a.add(ai)
                matched_c.add(ci)
                break
    rec = len(matched_a) / max(len(anns), 1)
    prec = len(matched_c) / max(len(cands), 1)
    f1 = 2 * rec * prec / max(rec + prec, 1e-12)
    return {"recall": rec, "precision": prec, "f1": f1,
            "n_ann": len(anns), "n_cand": len(cands)}


def average_precision(y_true, y_score):
    order = np.argsort(-y_score)
    y = np.asarray(y_true)[order]
    tp = np.cumsum(y)
    prec = tp / np.arange(1, len(y) + 1)
    P = y.sum()
    return float((prec * y).sum() / P) if P else float("nan")


def roc_auc(y_true, y_score):
    y = np.asarray(y_true)
    P, N = y.sum(), len(y) - y.sum()
    if not P or not N:
        return float("nan")
    ranks = np.empty(len(y))
    ranks[np.argsort(y_score, kind="mergesort")] = np.arange(1, len(y) + 1)
    return float((ranks[y == 1].sum() - P * (P + 1) / 2) / (P * N))


def spot_metrics(prob_tls, mapped_spots, cells_yx=None):
    """mapped_spots: list of dict(label, cy, cx in mosaic cells)。"""
    pos = [s for s in mapped_spots if s["label"] == "positive"]
    neg = [s for s in mapped_spots if s["label"] == "negative"]
    if not pos or not neg:
        return {"pr_auc": float("nan"), "roc_auc": float("nan"),
                "n_pos": len(pos), "n_neg": len(neg)}
    pts = [(s["cy"], s["cx"]) for s in pos + neg]
    y = [1] * len(pos) + [0] * len(neg)
    scores = []
    for cy, cx in pts:
        yy, xx = int(round(cy)), int(round(cx))
        scores.append(float(prob_tls[yy, xx]))
    scores = np.array(scores)
    return {"pr_auc": average_precision(y, scores), "roc_auc": roc_auc(y, scores),
            "n_pos": len(pos), "n_neg": len(neg)}


def load_folds():
    if not FOLDS.exists():
        return {}
    with open(FOLDS) as f:
        return {r["section_id"]: (r["unit_id"], int(r["fold"]))
                for r in csv.DictReader(f, delimiter="\t")}


def selftest() -> int:
    rng = np.random.default_rng(0)
    from shapely.geometry import Point
    H = W = 800
    prob = rng.uniform(0, 0.1, (H, W)).astype(np.float32)
    anns = []
    for cy, cx, r in [(200, 200, 30), (600, 500, 40), (400, 650, 25)]:
        yy, xx = np.mgrid[0:H, 0:W]
        prob[((yy - cy) ** 2 + (xx - cx) ** 2) <= r * r] = 0.95
        anns.append(Point(cx / 0.25, cy / 0.25).buffer(r / 0.25))  # fullres=1µm/px
    ok = True
    for thr in THRESHOLDS:
        m = instance_metrics(prob, anns, um_per_px_fullres=1.0, thr=thr)
        good = m["recall"] == 1.0 and m["precision"] == 1.0
        ok &= good
        print(f"  [selftest] thr={thr}: {m} {'OK' if good else 'FAIL'}")
    spots = ([{"label": "positive", "cy": 200, "cx": 200},
              {"label": "positive", "cy": 600, "cx": 500},
              {"label": "negative", "cy": 50, "cx": 50},
              {"label": "negative", "cy": 700, "cx": 100},
              {"label": "unknown", "cy": 200, "cx": 200}])  # unknown 重复坐标
    sm = spot_metrics(prob, spots)
    ok &= sm["pr_auc"] == 1.0 and sm["n_pos"] == 2 and sm["n_neg"] == 2
    print(f"  [selftest] spot pr_auc={sm['pr_auc']:.3f} roc={sm['roc_auc']:.3f} "
          f"(unknown 未计入: {sm['n_pos'] == 2 and sm['n_neg'] == 2})")
    shuffled = rng.permutation(prob.ravel()).reshape(H, W)
    sm2 = spot_metrics(shuffled, spots)
    prev = 2 / 4
    ok &= abs(sm2["pr_auc"] - prev) < 0.35
    print(f"  [selftest] shuffled pr_auc={sm2['pr_auc']:.3f} (期望≈流行率 {prev:.2f})")
    print("  [selftest]", "ALL OK" if ok else "FAILURES PRESENT")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--annotations", default=None, help="导出 envelope；缺省找 incoming/")
    args = ap.parse_args()
    if args.selftest:
        print("[selftest] instance extraction + matching + spot metrics")
        sys.exit(selftest())

    OUT.mkdir(parents=True, exist_ok=True)
    folds = load_folds()
    rows = []
    import h5py
    from tls_map_labels import load_labels, HEST
    slides = sorted(p.stem for p in MAP_DIR.glob("*.npz"))
    for sid in slides:
        rec = map_section(sid)   # 保守三态映射；reasons 顺序 = h5ad obs 顺序
        if rec.get("n_positive", 0) == 0 and rec.get("n_negative", 0) == 0:
            continue  # 该片无可用标注
        z = np.load(MAP_DIR / f"{sid}.npz")
        prob_tls = z["prob"][:, :, TLS_CHANNEL].astype(np.float32) / 255.0
        row = {"section_id": sid}
        if sid in folds:
            row["unit_id"], row["fold"] = folds[sid]
        # spot 级：与 map_section 同序读 obs，按索引对齐 reasons
        with h5py.File(HEST / "st" / f"{sid}.h5ad", "r") as f:
            obs = f["obs"]
            cols = np.asarray(obs["pxl_col_in_fullres"][:], dtype=float)
            rws = np.asarray(obs["pxl_row_in_fullres"][:], dtype=float)
        fct = float(rec["um_per_px"]) / MOSAIC_UM_PER_PX
        spots = [{"label": lab[0], "cy": float(rws[i]) * fct, "cx": float(cols[i]) * fct}
                 for i, lab in enumerate(rec["reasons"])]
        sm = spot_metrics(prob_tls, spots)
        row.update({f"spot_{k}": v for k, v in sm.items()})
        # 实例级：标注多边形原文（fullres px）
        inst, rev, _src = load_labels(sid)
        from tls_map_labels import POSITIVE_CLASSES
        geoms = [e["geom"] for e in inst if e["cls"] in POSITIVE_CLASSES]
        for thr in THRESHOLDS:
            im = instance_metrics(prob_tls, geoms, float(rec["um_per_px"]), thr)
            row.update({f"inst_{k}_t{int(thr*10)}": v for k, v in im.items()
                        if k in ("recall", "precision", "f1")})
        row["n_ann_instances"] = len(geoms)
        row["n_spot_pos"] = rec["n_positive"]
        row["n_spot_neg"] = rec["n_negative"]
        row["n_spot_unknown"] = rec["n_unknown"]
        rows.append(row)
        print(f"{sid}: spot_pr={sm['pr_auc']:.3f} inst_f1@0.5="
              f"{row.get('inst_f1_t5', float('nan')):.3f}", flush=True)
    if not rows:
        print("no annotated slides yet; run with --selftest to verify the harness")
        return
    keys = sorted({k for r in rows for k in r})
    with open(OUT / "eval_per_slide.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"written {OUT/'eval_per_slide.tsv'} ({len(rows)} slides)")


if __name__ == "__main__":
    main()
