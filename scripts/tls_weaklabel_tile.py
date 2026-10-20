#!/usr/bin/env python3
"""Route A probe: tile paired H&E around each spot, classify tiles with local
YOLO steps, map back to nearest spots. Weak labels only.

Usage: python3 scripts/tls_weaklabel_tile.py GSM5924053_frozen_c_57
"""
import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "infra/tls_weaklabel_probe_20261006" / "tiles"
RAW = ROOT / "data/GEO/GSE175540/raw"
CROP = 500  # hires px window around each spot


def main():
    gsm = sys.argv[1]
    WORK.mkdir(parents=True, exist_ok=True)
    with gzip.open(RAW / f"{gsm}_tissue_hires_image.png.gz", "rb") as f:
        img = Image.open(f).convert("RGB")
        img.load()
    with gzip.open(RAW / f"{gsm}_scalefactors_json.json.gz", "rt") as f:
        sf = json.load(f)
    scalef = float(sf["tissue_hires_scalef"])
    with gzip.open(RAW / f"{gsm}_tissue_positions_list.csv.gz", "rt") as f:
        rows = [l.strip().split(",") for l in f]
    header = rows[0]
    spots = [r for r in rows[1:] if r[1] == "1"]
    print(f"spots in tissue: {len(spots)}, hires={img.size}, scalef={scalef}", flush=True)

    from ultralytics import YOLO
    m1 = YOLO(str(ROOT / "repo/data/TLS_classification_HE_YOLO_step1.pt"))
    m2 = YOLO(str(ROOT / "repo/data/TLS_classification_HE_YOLO_step2.pt"))

    half = CROP // 2
    w, h = img.size
    out = []
    t0 = time.time()
    for i, r in enumerate(spots):
        barcode, _, arr_row, arr_col, py_full, px_full = r[0], r[1], r[2], r[3], float(r[4]), float(r[5])
        cx, cy = int(px_full * scalef), int(py_full * scalef)
        crop = img.crop((max(0, cx - half), max(0, cy - half),
                         min(w, cx + half), min(h, cy + half)))
        r1 = m1.predict(crop, imgsz=640, verbose=False)[0]
        r2 = m2.predict(crop, imgsz=640, verbose=False)[0]
        top1 = r1.names[int(r1.probs.top1)] if r1.probs is not None else "none"
        top2 = r2.names[int(r2.probs.top1)] if r2.probs is not None else "none"
        out.append({"barcode": barcode, "maturity": top1,
                    "conf1": float(r1.probs.top1conf) if r1.probs is not None else 0.0,
                    "subtype": top2,
                    "conf2": float(r2.probs.top1conf) if r2.probs is not None else 0.0})
        if (i + 1) % 200 == 0:
            print(f"{i+1}/{len(spots)} {time.time()-t0:.0f}s", flush=True)
    (WORK / f"{gsm}_tile_labels.json").write_text(json.dumps(out, indent=1))
    print(f"done {len(out)} tiles in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
