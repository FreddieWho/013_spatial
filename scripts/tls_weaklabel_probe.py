#!/usr/bin/env python3
"""Timing probe: local YOLO TLS detection on one paired H&E image (CPU only).

Weak labels only. Does not touch independent annotations.
"""
import gzip
import shutil
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "infra/tls_weaklabel_probe_20261006"
RAW = ROOT / "data/GEO/GSE175540/raw"


def main():
    gsm = sys.argv[1] if len(sys.argv) > 1 else "GSM5924053_frozen_c_57"
    WORK.mkdir(parents=True, exist_ok=True)
    src = RAW / f"{gsm}_tissue_hires_image.png.gz"
    img = WORK / f"{gsm}_hires.png"
    if not img.exists():
        with gzip.open(src, "rb") as f, open(img, "wb") as o:
            shutil.copyfileobj(f, o)
    print("image bytes:", img.stat().st_size, flush=True)

    from ultralytics import YOLO
    for step in (1, 2):
        w = ROOT / f"repo/data/TLS_classification_HE_YOLO_step{step}.pt"
        t0 = time.time()
        model = YOLO(str(w))
        out = WORK / f"{gsm}_yolo_step{step}"
        model(source=str(img), project=str(out), save_txt=True, imgsz=640, verbose=False)
        dt = time.time() - t0
        txts = sorted(out.rglob("*.txt"))
        n_boxes = sum(1 for t in txts for _ in t.read_text().splitlines())
        print(f"step{step}: {dt:.1f}s, txt_files={len(txts)}, boxes={n_boxes}", flush=True)
    (WORK / f"{gsm}_DONE").write_text("probe done\n")


if __name__ == "__main__":
    main()
