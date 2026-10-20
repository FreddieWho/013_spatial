#!/usr/bin/env python3
"""HookNet-TLS crop probe (hooknet conda env, CPU).

Builds HookNet from the installed package, loads the downloaded weights,
runs target+context patches for one TLS centroid and one negative spot.
"""
import gzip
import json
import time
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path("/home/huyudi/013_spatial")
WORK = ROOT / "infra/tls_weaklabel_probe_20261006" / "hooknet_crop"
GSM = "GSM5924053_frozen_c_57"
SPOT_UM = 0.49647431285076987  # measured fullres um/px from array pitch (100 um)


def load_fullres():
    with gzip.open(ROOT / f"data/GEO/GSE175540/raw/{GSM}_tissue_hires_image.png.gz", "rb") as f:
        img = Image.open(f).convert("RGB")
        img.load()
    with gzip.open(ROOT / f"data/GEO/GSE175540/raw/{GSM}_scalefactors_json.json.gz", "rt") as f:
        scalef = float(json.load(f)["tissue_hires_scalef"])
    full = img.resize((round(img.size[0] / scalef), round(img.size[1] / scalef)), Image.BILINEAR)
    return np.asarray(full).astype(np.float32) / 255.0


def patch_at(full, cx_um, cy_um, spacing, size=1244):
    factor = SPOT_UM / spacing
    small = Image.fromarray((full * 255).astype(np.uint8)).resize(
        (round(full.shape[1] * factor), round(full.shape[0] * factor)), Image.BILINEAR)
    arr = np.asarray(small).astype(np.float32) / 255.0
    cx, cy = int(cx_um / spacing), int(cy_um / spacing)
    x0, y0 = cx - size // 2, cy - size // 2
    canvas = np.ones((size, size, 3), np.float32)
    sx0, sy0 = max(0, x0), max(0, y0)
    sx1, sy1 = min(arr.shape[1], x0 + size), min(arr.shape[0], y0 + size)
    dx0, dy0 = sx0 - x0, sy0 - y0
    canvas[dy0:dy0 + (sy1 - sy0), dx0:dx0 + (sx1 - sx0)] = arr[sy0:sy1, sx0:sx1]
    return canvas


def main():
    import tensorflow as tf
    from tensorflow.keras.layers import LeakyReLU
    from hooknet.models.tensorflowmodel import HookNet
    WORK.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    full = load_fullres()
    print(f"fullres-equiv shape: {full.shape}", flush=True)
    geo = np.load(ROOT / "infra/tls_multiscale_gobp_20261002/geometry" / f"gse175540-{GSM}.npz")
    xy, labels = np.asarray(geo["xy"], float), np.array([str(x) for x in geo["labels"]])
    tls_pts = xy[labels == "TLS"]
    neg_pts = xy[labels == "NO_TLS"]
    tls_c = tls_pts.mean(axis=0)
    d = np.linalg.norm(neg_pts - tls_c, axis=1)
    neg_c = neg_pts[np.argmax(d)]
    print(f"TLS centroid um={tls_c}, neg um={neg_c}", flush=True)

    model = HookNet(input_shape=[[1244, 1244, 3], [1244, 1244, 3]],
                    hook_indices=[2, 4], n_classes=3, batch_norm=True,
                    n_filters=32, activation=LeakyReLU(alpha=0.1),
                    model_weights=str(ROOT / "infra/tls_weaklabel_probe_20261006/weights/weights.h5"))
    print(f"model built in {time.time()-t0:.0f}s", flush=True)
    for name, c in (("tls", tls_c), ("neg", neg_c)):
        t1 = time.time()
        pt = patch_at(full, c[0], c[1], 0.5)[None]
        pc = patch_at(full, c[0], c[1], 2.0)[None]
        pred = model.predict([pt, pc], verbose=0)[0]
        print(f"{name}: {time.time()-t1:.1f}s out_shape={pred.shape} "
              f"mean={pred.reshape(-1, pred.shape[-1]).mean(axis=0).round(4)} "
              f"max={pred.reshape(-1, pred.shape[-1]).max(axis=0).round(4)}", flush=True)
        np.save(WORK / f"{name}_pred.npy", pred)


if __name__ == "__main__":
    main()
