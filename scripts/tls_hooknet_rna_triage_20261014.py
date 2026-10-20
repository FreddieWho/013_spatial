"""HookNet 图像概率 x RNA TLS 签名的空间一致性（跨模态三角验证，D 记录见 decisions.md）。

为什么是合法的互证：图像模型输出与转录签名是两个独立模态，不存在
"用同一 RNA 签名验证同一 RNA 场"的循环。产出不作真值，只做三件事：
  1) 每片 Spearman 一致性 + 环面平移零假设 p 值（TLS 通道 ch2 为主，
     ch1=gc、ch0=tumor/其他 作特异对照——ch2 应 >= ch1 > ch0）；
  2) 患者级聚合（已知患者单元）；
  3) 图像高/RNA 低的冲突区域清单 → 专家标注分诊单。

坐标链（逐环验证过）：
  h5ad obs pxl_col/row_in_fullres（metadata fullres 帧）
    --x pixel_size_um_estimated--> 物理 µm
    --/ 4.0-----------------------> mosaic 格（sidecar: output_um_per_px=4.0）
  注意：infra/tls_fourcancer_20261006/geometry/*.npz 的 xy 是
  hex_xy(array_row, array_col) 的阵列布局坐标，不能用于图像映射。

spot 足迹：Visium 直径 55µm → 半径 27.5µm ≈ 6.875 格。
映射有效：足迹内 coverage>0 的格占比 >= 0.5，否则该 spot 记 unmapped。
"""

import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
MAP_DIR = ROOT / "infra/hooknet_gpu_20261008/maps"
SIDECAR_DIR = ROOT / "infra/hooknet_gpu_20261008"
SCORE_DIR = ROOT / "infra/tls_fourcancer_20261006/scores"
ST_DIR = ROOT / "data/other_sources/hest1k/st"
SECTIONS = ROOT / "infra/tls_four_cancer_labels_20261008/eligibility/sections.tsv"
OUT = ROOT / "infra/tls_hooknet_rna_triage_20261014"

TLS_CANONICAL = [
    "CXCL13", "MS4A1", "CD3D", "CD3E", "CD19", "CCL19", "CCL21", "LTB",
    "CXCR5", "CCR7", "CD79A", "BANK1", "SELL", "CR2",
]
MIN_PANEL_GENES = 10
SPOT_RADIUS_UM = 27.5          # Visium 55µm 直径
MOSAIC_UM_PER_PX = 4.0         # sidecar output_um_per_px
MIN_DISK_COVERAGE = 0.5        # 足迹内 coverage>0 占比下限
N_NULL = 100                   # 环面平移零假设次数
RNG_SEED = 20261014


def load_h5ad_spots(sid: str):
    """返回 (X csr in_tissue, genes, pxl_xy fullres, array_rc, barcodes)。"""
    with h5py.File(ST_DIR / f"{sid}.h5ad", "r") as f:
        g = f["X"]
        if isinstance(g, h5py.Group):
            data = g["data"][:]
            indices = g["indices"][:]
            indptr = g["indptr"][:]
            shape = tuple(g.attrs["shape"]) if "shape" in g.attrs else g.attrs["h5sparse_shape"]
            X = sparse.csr_matrix((data, indices, indptr), shape=shape)
            if X.shape[1] != len(f["var"]["_index"]):
                X = X.T.tocsr()
        else:
            X = sparse.csr_matrix(g[:])
        genes = np.array([b.decode() if isinstance(b, bytes) else str(b)
                          for b in f["var"]["_index"][:]])
        obs = f["obs"]
        it = np.asarray(obs["in_tissue"][:]).astype(int) == 1
        px = np.asarray(obs["pxl_col_in_fullres"][:], dtype=float)
        py = np.asarray(obs["pxl_row_in_fullres"][:], dtype=float)
        bc = np.array([b.decode() if isinstance(b, bytes) else str(b)
                       for b in obs["_index"][:]])
    return X[it].tocsr(), genes, np.stack([px[it], py[it]], axis=1), bc[it]


def tls_signature(X, genes):
    """log1p(CPM) 逐基因 z 后取面板均值。返回 (score, n_genes_present)。"""
    idx = {g: i for i, g in enumerate(genes)}
    present = [g for g in TLS_CANONICAL if g in idx]
    if len(present) < MIN_PANEL_GENES:
        return None, len(present)
    lib = np.asarray(X.sum(axis=1)).ravel()
    lib = np.where(lib > 0, lib, 1.0)
    cols = []
    for gname in present:
        x = np.asarray(X[:, idx[gname]].todense()).ravel()
        v = np.log1p(1e4 * x / lib)
        sd = v.std()
        cols.append((v - v.mean()) / sd if sd > 1e-12 else v * 0.0)
    return np.mean(cols, axis=0), len(present)


def disk_offsets(radius_cells: float):
    r = int(np.ceil(radius_cells))
    offs = [(dy, dx) for dy in range(-r, r + 1) for dx in range(-r, r + 1)
            if dy * dy + dx * dx <= radius_cells ** 2]
    return np.array(offs, dtype=int)


def spot_disk_means(prob, coverage, cells_yx, offs):
    """cells_yx: (n,2) float mosaic 格坐标。返回 (means (n,3), frac_valid (n,))。"""
    h, w = coverage.shape
    cy = np.rint(cells_yx[:, 0]).astype(int)
    cx = np.rint(cells_yx[:, 1]).astype(int)
    yy = cy[:, None] + offs[None, :, 0]
    xx = cx[:, None] + offs[None, :, 1]
    inside = (yy >= 0) & (yy < h) & (xx >= 0) & (xx < w)
    yyc = np.clip(yy, 0, h - 1)
    xxc = np.clip(xx, 0, w - 1)
    cov = coverage[yyc, xxc] > 0
    valid = inside & cov
    frac = valid.sum(axis=1) / offs.shape[0]
    means = np.zeros((len(cells_yx), prob.shape[2]), dtype=np.float32)
    for c in range(prob.shape[2]):
        vals = prob[yyc, xxc, c].astype(np.float32)
        vals[~valid] = 0.0
        means[:, c] = vals.sum(axis=1) / np.maximum(valid.sum(axis=1), 1) / 255.0
    return means, frac


def null_shift_p(prob, coverage, cells_yx, offs, score, ok, rho_obs, rng, n=N_NULL):
    """环面平移零假设：固定 spot，平移概率图重算 rho。"""
    h, w = coverage.shape
    better = 0
    m1 = np.ones_like(coverage, dtype=np.uint8)  # 平移后全部"有覆盖"
    for _ in range(n):
        sy = int(rng.integers(1, h))
        sx = int(rng.integers(1, w))
        shifted = np.roll(np.roll(prob, sy, axis=0), sx, axis=1)
        means, _ = spot_disk_means(shifted, m1, cells_yx, offs)
        r = spearmanr(means[ok, 2], score[ok]).statistic
        if np.isfinite(r) and r >= rho_obs:
            better += 1
    return (1 + better) / (n + 1)


def main():
    only = sys.argv[1].split(",") if len(sys.argv) > 1 else None
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)
    offs = disk_offsets(SPOT_RADIUS_UM / MOSAIC_UM_PER_PX)

    import csv
    with open(SECTIONS) as f:
        sec = {r["section_id"]: r for r in csv.DictReader(f, delimiter="\t")}

    maps = {p.stem for p in MAP_DIR.glob("*.npz")}
    scored = {p.stem[5:] for p in SCORE_DIR.glob("hest-*.npz")}
    todo = sorted(maps & scored)
    if only:
        todo = [s for s in todo if s in only]
    print(f"slides with map+scores: {len(maps & scored)}, running {len(todo)}", flush=True)

    rows = []
    triage_spots = {}
    for i, sid in enumerate(todo):
        t1 = time.time()
        row = {"section_id": sid, "organ": sec[sid]["organ"],
               "patient_key": sec[sid]["patient_key"],
               "patient_status": sec[sid]["patient_status"]}
        try:
            ps = float(sec[sid]["pixel_size_um_estimated"])
            X, genes, pxl, bc = load_h5ad_spots(sid)
            score, npanel = tls_signature(X, genes)
            row["n_panel_genes"] = npanel
            if score is None:
                row["status"] = "PANEL_TOO_THIN"
                rows.append(row)
                print(f"[{i+1}/{len(todo)}] {sid}: panel thin ({npanel})", flush=True)
                continue
            z = np.load(MAP_DIR / f"{sid}.npz")
            prob, coverage = z["prob"], z["coverage"]
            # fullres px -> µm -> mosaic 格
            cells = pxl * ps / MOSAIC_UM_PER_PX
            cells_yx = np.stack([cells[:, 1], cells[:, 0]], axis=1)
            means, frac = spot_disk_means(prob, coverage, cells_yx, offs)
            ok = frac >= MIN_DISK_COVERAGE
            row["n_spots"] = len(cells_yx)
            row["n_mapped"] = int(ok.sum())
            if ok.sum() < 50:
                row["status"] = "TOO_FEW_MAPPED"
                rows.append(row)
                print(f"[{i+1}/{len(todo)}] {sid}: mapped {ok.sum()} < 50", flush=True)
                continue
            rho_tls = spearmanr(means[ok, 2], score[ok]).statistic
            rho_gc = spearmanr(means[ok, 1], score[ok]).statistic
            rho_other = spearmanr(means[ok, 0], score[ok]).statistic
            if not np.isfinite(rho_tls):
                # 常数通道（如全零 TLS 图）spearman 得 NaN；NaN 参与零假设比较
                # 会假性显著（所有 null >= NaN 为 False → p=0.0099），必须单列
                row.update(status="CONSTANT_CHANNEL", rho_tls=rho_tls,
                           rho_gc=rho_gc, rho_other=rho_other, p_shift="")
                rows.append(row)
                print(f"[{i+1}/{len(todo)}] {sid}: rho NaN (constant channel), "
                      f"mapped {ok.sum()}", flush=True)
                continue
            p_null = null_shift_p(prob, coverage, cells_yx, offs, score, ok, rho_tls, rng)
            row.update(status="OK", rho_tls=rho_tls, rho_gc=rho_gc, rho_other=rho_other,
                       p_shift=p_null)
            # 冲突区：图像前 10% 但 RNA 后半（标注优先）
            thr_img = np.quantile(means[ok, 2], 0.9)
            thr_rna = np.quantile(score[ok], 0.5)
            conf = np.flatnonzero(ok)[(means[ok, 2] >= thr_img) & (score[ok] <= thr_rna)]
            order = conf[np.argsort(-means[conf, 2])[:10]]
            triage_spots[sid] = [
                {"barcode": bc[j], "fullres_x": round(float(pxl[j, 0]), 1),
                 "fullres_y": round(float(pxl[j, 1]), 1),
                 "img_tls": round(float(means[j, 2]), 3),
                 "rna_tls": round(float(score[j]), 3)}
                for j in order]
            print(f"[{i+1}/{len(todo)}] {sid}: rho_tls={rho_tls:.3f} "
                  f"(gc {rho_gc:.3f}, other {rho_other:.3f}) p={p_null:.3f} "
                  f"mapped {ok.sum()}/{len(ok)} [{time.time()-t1:.0f}s]", flush=True)
        except Exception as e:
            row["status"] = f"FAIL: {type(e).__name__}: {e}"
            print(f"[{i+1}/{len(todo)}] {sid}: {row['status']}", flush=True)
        rows.append(row)

    # 写表
    keys = sorted({k for r in rows for k in r})
    with open(OUT / "per_slide.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    (OUT / "triage_spots.json").write_text(json.dumps(triage_spots, indent=1) + "\n")

    okrows = [r for r in rows if r.get("status") == "OK"]
    summary = {
        "n_ok": len(okrows), "n_total": len(rows),
        "median_rho_tls": float(np.nanmedian([r["rho_tls"] for r in okrows])) if okrows else None,
        "median_rho_gc": float(np.nanmedian([r["rho_gc"] for r in okrows])) if okrows else None,
        "median_rho_other": float(np.nanmedian([r["rho_other"] for r in okrows])) if okrows else None,
        "frac_p_lt_0.05": float(np.mean([r["p_shift"] < 0.05 for r in okrows])) if okrows else None,
        "channel_order_ok": int(sum(r["rho_tls"] >= r["rho_gc"] >= r["rho_other"]
                                    for r in okrows)),
        "elapsed_s": round(time.time() - t0, 1),
        "params": {"panel": TLS_CANONICAL, "min_panel": MIN_PANEL_GENES,
                   "spot_radius_um": SPOT_RADIUS_UM, "mosaic_um_per_px": MOSAIC_UM_PER_PX,
                   "min_disk_coverage": MIN_DISK_COVERAGE, "n_null": N_NULL,
                   "seed": RNG_SEED},
        "status_counts": {},
    }
    from collections import Counter
    summary["status_counts"] = dict(Counter(r["status"] for r in rows))
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "params"}, indent=1), flush=True)


if __name__ == "__main__":
    main()
