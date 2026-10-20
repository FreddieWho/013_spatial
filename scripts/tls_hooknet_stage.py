#!/usr/bin/env python3
"""Stage the minimal input a HookNet GPU run needs, so a rented box gets ~4 GB
instead of the 649 GB of source pyramids.

Per section:
  <sid>.jpg     the source level at or below 0.6 um/px, resampled to EXACTLY
                0.5 um/px so the target branch matches the model's training
                spacing; the context branch is a 4x downsample of this file
  <sid>.json    dims, spacing, tissue bbox, source provenance
  positions.json  the patch grid on a declared stride, inside the tissue mask

Nothing about the model or the labels is decided here; this only moves pixels.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import time
import warnings
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

warnings.filterwarnings("ignore")
Image.MAX_IMAGE_PIXELS = None

ROOT = Path("/home/huyudi/013_spatial")
H = ROOT / "data/other_sources/hest1k"
GATE = ROOT / "infra/tls_four_cancer_labels_20261008/eligibility/annotatable_gate.tsv"
OUT = ROOT / "dist/hooknet_gpu_input"

TARGET_UM = 0.5          # model target-branch spacing
SOURCE_MAX_UM = 0.6      # coarsest source level we accept
JPEG_Q = 90


def pick_level(sid: str):
    meta = json.loads((H / "metadata" / f"{sid}.json").read_text())
    fw = int(meta["fullres_px_width"])
    px = float(meta["pixel_size_um_estimated"])
    with Image.open(H / "wsis" / f"{sid}.tif") as im:
        n = getattr(im, "n_frames", 1)
        sizes = []
        for i in range(n):
            im.seek(i)
            sizes.append(im.size)
    lvl = 0
    for i in range(n):
        if (px / (sizes[i][0] / fw)) <= SOURCE_MAX_UM:
            lvl = i
    um = px / (sizes[lvl][0] / fw)
    return lvl, sizes[lvl], um, fw


def tissue_mask(sid: str, w: int, h: int, fw: int) -> np.ndarray:
    scale = w / fw
    m = Image.new("L", (w, h), 0)
    d = ImageDraw.Draw(m)
    for feat in json.loads((H / "tissue_seg" / f"{sid}_contours.geojson").read_text())["features"]:
        g = feat["geometry"]
        c = g.get("coordinates") or []
        rings = c if g.get("type") == "Polygon" else [r for p in c for r in p]
        for i, ring in enumerate(rings):
            pts = [(p[0] * scale, p[1] * scale) for p in ring]
            if len(pts) >= 3:
                d.polygon(pts, fill=255 if i == 0 else 0)
    return np.asarray(m) > 127


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--stride-um", type=float, default=311.0,
                    help="50%% overlap for HookNet's 622 um field")
    ap.add_argument("--sections", nargs="*", default=None)
    ap.add_argument("--quality", type=int, default=JPEG_Q)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "img").mkdir(exist_ok=True)

    rows = [r for r in csv.DictReader(GATE.open(), delimiter="\t")
            if r["gate_status_human"] == "annotatable_candidate"]
    if args.sections:
        rows = [r for r in rows if r["section_id"] in args.sections]
    rows.sort(key=lambda r: (r["organ"], r["section_id"]))

    manifest = []
    t_all = time.time()
    for i, r in enumerate(rows, 1):
        sid = r["section_id"]
        t0 = time.time()
        lvl, (lw, lh), source_um, fw = pick_level(sid)
        with Image.open(H / "wsis" / f"{sid}.tif") as im:
            im.seek(lvl)
            img = im.convert("RGB")
        # resample so one pixel is exactly TARGET_UM
        f = source_um / TARGET_UM
        tw, th = int(round(lw * f)), int(round(lh * f))
        if abs(f - 1.0) > 1e-6:
            img = img.resize((tw, th), Image.BILINEAR)
        p = out / "img" / f"{sid}.jpg"
        img.save(p, quality=args.quality, optimize=False, subsampling=0)

        mask = tissue_mask(sid, tw, th, fw)
        ys, xs = np.nonzero(mask)
        if len(xs) == 0:
            bbox = (0, 0, 1, 1)
        else:
            bbox = (int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1)

        # patch grid on the declared stride, only where the patch centre is tissue
        step = args.stride_um / TARGET_UM
        gx = np.arange(step / 2, tw, step)
        gy = np.arange(step / 2, th, step)
        pos = []
        for yy in gy:
            iy = int(round(yy))
            for xx in gx:
                ix = int(round(xx))
                if mask[iy, ix]:
                    pos.append([ix, iy])

        rec = {
            "section_id": sid, "organ": r["organ"],
            "disease_state": r["disease_state"],
            "image": f"img/{sid}.jpg",
            "width": tw, "height": th,
            "um_per_px": TARGET_UM,
            "source_level": lvl, "source_um_per_px": round(source_um, 5),
            "resample_factor": round(f, 5),
            "fullres_px": [fw, int(json.loads((H / 'metadata' / f'{sid}.json').read_text())["fullres_px_height"])],
            "tissue_bbox_px": list(bbox),
            "n_positions": len(pos),
            "jpeg_bytes": p.stat().st_size,
            "patient_status": r["patient_status"],
            "specimen_key": r["specimen_key"],
            "seconds": round(time.time() - t0, 1),
        }
        (out / f"{sid}.json").write_text(json.dumps(rec, indent=1) + "\n")
        (out / f"{sid}.positions.json").write_text(json.dumps(pos))
        manifest.append(rec)
        print(f'[{i}/{len(rows)}] {sid:<9} {tw}x{th} from L{lvl}@{source_um:.3f} '
              f'-> {p.stat().st_size/1e6:6.1f} MB  positions {len(pos):>5}  '
              f'{rec["seconds"]}s', flush=True)

    doc = {
        "generated": time.strftime("%Y-%m-%d"),
        "purpose": "minimal HookNet input for a rented GPU; no labels here",
        "target_um_per_px": TARGET_UM,
        "context_factor": 4,
        "context_um_per_px": TARGET_UM * 4,
        "patch_px": 1244,
        "stride_um": args.stride_um,
        "overlap": "50%",
        "field_of_view_um": 1244 * TARGET_UM,
        "jpeg_quality": args.quality,
        "n_slides": len(manifest),
        "total_positions": sum(m["n_positions"] for m in manifest),
        "total_image_bytes": sum(m["jpeg_bytes"] for m in manifest),
        "slides": manifest,
    }
    (out / "manifest.json").write_text(json.dumps(doc, indent=1) + "\n")
    tot = sum(p.stat().st_size for p in out.rglob("*") if p.is_file())
    print(f'\npack {out}: {tot/1e9:.2f} GB, {len(manifest)} slides, '
          f'{doc["total_positions"]} positions, {time.time()-t_all:.0f}s')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())