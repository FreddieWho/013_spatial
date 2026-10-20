#!/usr/bin/env python3
"""iStar probe input prep (D-220): 1 HEST (TENX13, true hi-res) + 1 GSE (c_57, hires).

HEST arm: register wsi tif to hires space via FFT phase correlation, map
pxl_fullres coords into tif pixels. GSE arm: hires PNG directly.
Probe genes: repo TLS_signature/genes.txt intersected with matrix genes.
"""
import csv
import gzip
import json
import sys
from pathlib import Path

import h5py
import numpy as np
from PIL import Image
from scipy import sparse

Image.MAX_IMAGE_PIXELS = None
ROOT = Path("/home/huyudi/013_spatial")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
OUT = ROOT / "infra/istar_probe_20261007"
SIG = [l.strip() for l in (ROOT / "repo/data/TLS_signature/genes.txt").read_text().splitlines()
       if l.strip()]


def load_h5_counts(path):
    with h5py.File(path, "r") as f:
        X = f["X"]
        if isinstance(X, h5py.Group):
            n_obs = len(f["obs"]["_index"])
            data, indices, indptr = X["data"][:], X["indices"][:], X["indptr"][:]
            if len(indptr) - 1 == n_obs:
                X = sparse.csr_matrix((data, indices, indptr))
            else:
                X = sparse.csc_matrix((data, indices, indptr)).tocsr()
        else:
            X = sparse.csr_matrix(X[:])
        genes = np.array([g.decode() if isinstance(g, bytes) else str(g)
                          for g in f["var"]["_index"][:]])
    return X, genes


def tissue_bbox(mask):
    ys, xs = np.where(mask)
    return xs.min(), ys.min(), xs.max(), ys.max()


def register_tif_to_hires(tif_path, hires_arr, spot_xy_hires=None):
    """Return (scale, dx, dy): tif_px = hires_px*scale + (dx,dy), via thumbnails."""
    tif = Image.open(tif_path).convert("RGB")
    tw, th = tif.size
    hw, hh = hires_arr.shape[1], hires_arr.shape[0]
    # tissue masks via Otsu on grayscale
    import cv2
    g_tif = np.asarray(tif.resize((hw, hh), Image.BILINEAR).convert("L"))
    _, m_tif = cv2.threshold(g_tif, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    tissue_tif = m_tif < 128  # tissue darker than background
    if tissue_tif.mean() > 0.9:  # inverted guess
        tissue_tif = ~tissue_tif
    g_h = hires_arr.mean(axis=2).astype(np.uint8) if hires_arr.ndim == 3 else hires_arr
    _, m_h = cv2.threshold(g_h, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    tissue_h = m_h < 128
    if tissue_h.mean() > 0.9:
        tissue_h = ~tissue_h
    x0t, y0t, x1t, y1t = tissue_bbox(tissue_tif)
    x0h, y0h, x1h, y1h = tissue_bbox(tissue_h)
    # affine from hires frame to full-tif frame via tissue bboxes
    scx = (x1t - x0t) * (tw / hw) / max(1, (x1h - x0h))
    scy = (y1t - y0t) * (th / hh) / max(1, (y1h - y0h))
    if abs(scx - scy) / max(scx, scy) > 0.05:
        print(f"WARN anisotropic bbox scale x={scx:.3f} y={scy:.3f}; possible rotation/padding", flush=True)
    s = float((scx + scy) / 2)
    ox = x0t * (tw / hw) - x0h * scx
    oy = y0t * (th / hh) - y0h * scy
    # validate with FFT residual shift at fitted scale
    w2, h2 = round(tw / s), round(th / s)
    r = np.asarray(tif.resize((w2, h2), Image.BILINEAR)).astype(float)
    r = r.mean(axis=2) if r.ndim == 3 else r
    ref = hires_arr.astype(float).mean(axis=2) if hires_arr.ndim == 3 else hires_arr.astype(float)
    y0, x0 = max(0, (h2 - hh) // 2), max(0, (w2 - hw) // 2)
    r = r[y0:y0 + hh, x0:x0 + hw]
    F = np.fft.fft2(r - r.mean())
    F2 = np.fft.fft2(ref - ref.mean())
    cc = np.abs(np.fft.fftshift(np.fft.ifft2(F * np.conj(F2))))
    ddy, ddx = np.unravel_index(np.argmax(cc), cc.shape)
    peak = float(cc[ddy, ddx] / (cc.mean() + 1e-9))
    print(f"bbox tif={(x0t,y0t,x1t,y1t)} hires={(x0h,y0h,x1h,y1h)} scale={s:.4f} off=({ox:.1f},{oy:.1f})"
          f" fft_resid=({ddx-hw//2},{ddy-hh//2}) peak_ratio={peak:.1f}", flush=True)
    # fold FFT residual (in hires px) into offset
    ox = float(ox + (ddx - hw // 2) * s)
    oy = float(oy + (ddy - hh // 2) * s)
    return s, ox, oy


def prep_hest(sid="TENX13", um_px=0.27):
    d = OUT / f"HEST_{sid}"
    d.mkdir(parents=True, exist_ok=True)
    X, genes = load_h5_counts(ROOT / f"data/other_sources/hest1k/st/{sid}.h5ad")
    keep = [g for g in SIG if g in set(genes)]
    gi = {g: i for i, g in enumerate(genes)}
    cols = [gi[g] for g in keep]
    with h5py.File(ROOT / f"data/other_sources/hest1k/st/{sid}.h5ad", "r") as f:
        o = f["obs"]
        bc = np.array([b.decode() if isinstance(b, bytes) else str(b) for b in o["_index"][:]])
        pr = np.asarray(o["pxl_row_in_fullres"][:]).astype(float)
        pc = np.asarray(o["pxl_col_in_fullres"][:]).astype(float)
    sp = ROOT / f"data/other_sources/10x_genomics/V1_Breast_Cancer_Block_A_Section_1/spatial"
    hires = np.asarray(Image.open(sp / "tissue_hires_image.png").convert("RGB"))
    tif_path = ROOT / f"data/other_sources/hest1k/wsis/{sid}.tif"
    s, ox, oy = register_tif_to_hires(tif_path, hires)
    print(f"HEST {sid}: reg scale={s:.4f} off=({ox:.1f},{oy:.1f})", flush=True)
    # he-raw.jpg = full tif converted (keep native resolution)
    jpg = d / "he-raw.jpg"
    if not jpg.exists():
        Image.open(tif_path).convert("RGB").save(jpg, quality=92)
    # coords in tif px: fullres->hires (known scalef) -> tif via registration
    import json as js
    sf = js.load(open(sp / "scalefactors_json.json"))
    hsf = float(sf["tissue_hires_scalef"])
    hx = pc * hsf
    hy = pr * hsf
    tx = hx * s + ox
    ty = hy * s + oy
    sub = X[:, cols].toarray()
    with open(d / "cnts.tsv", "w") as fh:
        fh.write("spot\t" + "\t".join(keep) + "\n")
        for b, row in zip(bc, sub):
            fh.write(b + "\t" + "\t".join(str(int(v)) for v in row) + "\n")
    with open(d / "locs-raw.tsv", "w") as fh:
        fh.write("spot\tx\ty\n")
        for b, x, y in zip(bc, tx, ty):
            fh.write(f"{b}\t{x:.1f}\t{y:.1f}\n")
    (d / "pixel-size-raw.txt").write_text(f"{um_px}\n")
    (d / "radius-raw.txt").write_text(f"{27.5 / um_px:.1f}\n")
    (d / "gene-names.txt").write_text("\n".join(keep) + "\n")
    (d / "REGISTRATION.txt").write_text(
        f"scale={s} off=({ox},{oy}) um_px={um_px} n_genes={len(keep)}\n")
    print(f"HEST {sid}: spots={len(bc)} genes={len(keep)}", flush=True)


def prep_gse(gsm="GSM5924053_frozen_c_57"):
    from scripts.tls_peripheral import source_metadata, counts  # noqa
    d = OUT / f"GSE_{gsm}"
    d.mkdir(parents=True, exist_ok=True)
    sid = "gse175540-" + gsm
    m = source_metadata("GSE175540", sid)
    X, genes = counts("GSE175540", sid, m)
    keep = [g for g in SIG if g in set(genes)]
    gi = {g: i for i, g in enumerate(genes)}
    cols = [gi[g] for g in keep]
    with gzip.open(ROOT / f"data/GEO/GSE175540/raw/{gsm}_tissue_hires_image.png.gz", "rb") as fh:
        img = Image.open(fh).convert("RGB")
        img.load()
    img.save(d / "he-raw.jpg", quality=92)
    with gzip.open(ROOT / f"data/GEO/GSE175540/raw/{gsm}_scalefactors_json.json.gz", "rt") as fh:
        sf = json.load(fh)
    hsf = float(sf["tissue_hires_scalef"])
    sdf = float(sf["spot_diameter_fullres"])
    # hires um/px = (55 um / sdf fullres-px) / hires_scalef
    um_px = (55.0 / sdf) / hsf
    print(f"GSE {gsm}: hires {img.size} um_px={um_px:.2f}", flush=True)
    sub = X[:, cols].toarray() if sparse.issparse(X) else np.asarray(X[:, cols])
    # barcodes from matrix order
    with h5py.File(ROOT / f"data/GEO/GSE175540/raw/{gsm}_filtered_feature_bc_matrix.h5", "r") as f:
        bc_all = [b.decode() if isinstance(b, bytes) else str(b) for b in f["matrix/barcodes"][:]]
    bc = [bc_all[i] for i in m["order"]]
    with gzip.open(ROOT / f"data/GEO/GSE175540/raw/{gsm}_tissue_positions_list.csv.gz", "rt") as fh:
        import csv as cv
        rows = list(cv.reader(fh))
    # cols: barcode,in_tissue,array_row,array_col,pxl_row_fullres,pxl_col_fullres -> hires x=col,y=row
    pos = {r[0]: (float(r[5]) * hsf, float(r[4]) * hsf) for r in rows if r[1] == "1"}
    with open(d / "cnts.tsv", "w") as fh:
        fh.write("spot\t" + "\t".join(keep) + "\n")
        for b, row in zip(bc, sub):
            fh.write(b + "\t" + "\t".join(str(int(v)) for v in row) + "\n")
    with open(d / "locs-raw.tsv", "w") as fh:
        fh.write("spot\tx\ty\n")
        for b in bc:
            x, y = pos[b]
            fh.write(f"{b}\t{x:.1f}\t{y:.1f}\n")
    (d / "pixel-size-raw.txt").write_text(f"{um_px}\n")
    (d / "radius-raw.txt").write_text(f"{27.5 / um_px:.1f}\n")
    (d / "gene-names.txt").write_text("\n".join(keep) + "\n")
    print(f"GSE {gsm}: spots={len(bc)} genes={len(keep)}", flush=True)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    if which in ("hest", "both"):
        prep_hest()
    if which in ("gse", "both"):
        prep_gse()
