#!/usr/bin/env python3
"""Build overlay heat tiles from HookNet TLS probability mosaics.

The mosaics live at 4 um/px over the full staged frame. The web pack's base
image lives at base_um_per_px (0.9-2.2 um/px) over the same physical frame.
Heat tiles are built in *base coordinates* so the viewer draws them with the
exact same tile math as the H&E layer -- no second coordinate system.

Heat L0 = mosaic native resolution (4 um/px); each further level halves.
Nothing is upsampled: the heat layer never pretends to have more detail.
Only the TLS channel (index 2) is shipped.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

ROOT = Path("/home/huyudi/013_spatial")
PACK = ROOT / "dist/tls_annotation_pack"
TILE = 256
TLS_CH = 2


def load_gate(path: Path) -> dict:
    out = {}
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["section_id"]] = r
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sections", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    with (PACK / "data/index.json").open() as fh:
        index = json.load(fh)
    slides = index["slides"]
    if args.sections:
        slides = [s for s in slides if s["section_id"] in args.sections]
    if args.limit:
        slides = slides[: args.limit]

    ok, skipped = 0, []
    for i, s in enumerate(slides, 1):
        sid = s["section_id"]
        mp = PACK / "data" / sid
        npz = ROOT / "infra/hooknet_gpu_20261008/maps" / f"{sid}.npz"
        if not npz.exists():
            skipped.append((sid, "no mosaic"))
            continue
        z = np.load(npz)
        prob = z["prob"]
        if prob.shape[2] <= TLS_CH:
            skipped.append((sid, "mosaic has no TLS channel"))
            continue
        tls = prob[:, :, TLS_CH].astype(np.uint8)          # 0..255
        mh, mw = tls.shape
        base_um = float(s["base_um_per_px"])
        # Heat L0 stays at mosaic native resolution (4 um/px). Upsampling it to
        # base resolution would multiply bytes ~19x with zero new information;
        # the viewer maps heat px -> base px with one recorded factor instead.
        img = Image.fromarray(tls, mode="L")
        levels = []
        L = 0
        cur = img
        hdir = mp / "heat"
        while True:
            w, h = cur.size
            levels.append({"level": L, "w": w, "h": h,
                           "um_per_px": round(4.0 * (mw / w), 5)})
            ld = hdir / str(L)
            ld.mkdir(parents=True, exist_ok=True)
            for row in range((h + TILE - 1) // TILE):
                for col in range((w + TILE - 1) // TILE):
                    box = (col * TILE, row * TILE,
                           min((col + 1) * TILE, w), min((row + 1) * TILE, h))
                    t = cur.crop(box)
                    if t.size[0] < 8 or t.size[1] < 8:
                        t = t.resize((max(t.size[0], 8), max(t.size[1], 8)))
                    t.save(ld / f"{col}_{row}.png")
            if max(w, h) <= TILE * 4:
                break
            L += 1
            cur = cur.resize((max(1, w // 2), max(1, h // 2)), Image.BILINEAR)
        s["heat"] = {"levels": levels, "tile": TILE,
                     "aligned_by": "heat px * (4.0*2^L/base_um) = base px",
                     "mosaic_shape": [mh, mw],
                     "mosaic_um_per_px": 4.0,
                     "channel": TLS_CH,
                     "note": "HookNet TLS probability 0..255; blind-pass rule: overlay use re-tiers labels to reviewed_candidate"}
        ok += 1
        print(f"[{i}/{len(slides)}] {sid:<9} heat L0 {mw}x{mh} "
              f"({len(levels)} levels)", flush=True)

    with (PACK / "data/index.json").open("w") as fh:
        json.dump(index, fh, indent=1)
    print(f"\nheat tiles for {ok} slides; skipped: {skipped or 'none'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
