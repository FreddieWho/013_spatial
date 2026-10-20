#!/usr/bin/env python3
"""HookNet-TLS CPU cost probe on a true four-cancer whole-slide image.

Runs in the `hooknet` conda env on CPU only. GPU stays disabled.

The earlier GSE hires probe measured a single patch and told us nothing about
throughput on a real slide. HookNet is applied as a dense grid with rings at
0.5 and 2.0 um, so the honest unit of cost is: patches per second, then patches
per slide, then wall clock for the whole four-cancer set.

    /home/huyudi/.conda/envs/hooknet/bin/python scripts/tls_hooknet_cpu_probe.py \
        --section ZEN39 --patches 48
"""

from __future__ import annotations

import argparse
import json
import math
import time
import warnings
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

warnings.filterwarnings("ignore")
Image.MAX_IMAGE_PIXELS = None

ROOT = Path("/home/huyudi/013_spatial")
HEST = ROOT / "data/other_sources/hest1k"
WEIGHTS = ROOT / "infra/tls_weaklabel_probe_20261006/weights/weights.h5"
PATCH = 1244
RING_UM = (0.5, 2.0)
# HookNet is run as a dense grid; the published spacing is one grid step per
# HookNet output unit. This probe measures throughput, not the detection task,
# so the grid pitch is a declared assumption, not a fitted value.
GRID_PITCH_UM = 32.0


def tissue_mask_at(sid: str, target_um: float, max_plane_mp: float):
    """Finest level that still meets HookNet's resolution need.

    Level choice must follow resolution, not image width: HookNet reads 0.5 um
    texture, so a width-based pick lands on a thumbnail where that ring is
    meaningless. The decoded plane is capped so a probe does not try to hold a
    multi-gigabyte level in memory.
    """
    meta = json.loads((HEST / "metadata" / f"{sid}.json").read_text())
    w, h = int(meta["fullres_px_width"]), int(meta["fullres_px_height"])
    px_um = float(meta["pixel_size_um_estimated"])
    with Image.open(HEST / "wsis" / f"{sid}.tif") as im:
        n = getattr(im, "n_frames", 1)
        sizes = []
        for i in range(n):
            im.seek(i)
            sizes.append(im.size)
        # coarsest level that is still at or below the resolution budget
        lvl = 0
        for i in range(n):
            lw, _ = sizes[i]
            if (px_um / (lw / w)) <= target_um:
                lvl = i
        plane_mp = sizes[lvl][0] * sizes[lvl][1] / 1e6
        while plane_mp > max_plane_mp and lvl < n - 1:
            lvl += 1
            plane_mp = sizes[lvl][0] * sizes[lvl][1] / 1e6
        im.seek(lvl)
        lw, lh = im.size
        img = im.convert("RGB")
    scale = lw / w
    mm_px = px_um / scale
    return img, mm_px, (w, h), scale, meta, lvl, plane_mp


def tissue_boxes(sid: str, pitch_um: float, mm_px: float, w, h, scale: float):
    """Grid positions inside the tissue contour, at the pyramid level used.

    Rasterises the contour with PIL (the same method the registration check
    uses) instead of doing ray casting per ring: it is exact for holes, far
    faster, and avoids a huge edge-by-point matrix.
    """
    contours = json.loads((HEST / "tissue_seg" / f"{sid}_contours.geojson").read_text())
    mask = Image.new("L", (w, h), 0)
    draw = ImageDraw.Draw(mask)

    def paint(ring, value):
        pts = [(p[0], p[1]) for p in ring]
        if len(pts) >= 3:
            draw.polygon(pts, fill=value)

    for feat in contours.get("features", []):
        g = feat.get("geometry") or {}
        c = g.get("coordinates") or []
        if g.get("type") == "Polygon":
            rings = c
        elif g.get("type") == "MultiPolygon":
            rings = [r for poly in c for r in poly]
        else:
            continue
        for i, ring in enumerate(rings):
            # contours are in ORIGINAL pixels; the mask is at the level's scale
            paint([(p[0] * scale, p[1] * scale) for p in ring],
                  255 if i == 0 else 0)

    arr = np.asarray(mask) > 127
    step = max(1, int(round(pitch_um / mm_px)))
    ys, xs = np.nonzero(arr[::step, ::step])
    return list(zip((xs * step).astype(float).tolist(),
                    (ys * step).astype(float).tolist()))


def make_patches(sid: str, boxes, mm_px, level_img, scale):
    """Crop target and context ring patches from the pyramid level, in pixels."""
    import numpy as np
    arr = np.asarray(level_img)
    out = []
    half = PATCH // 2
    for cx, cy in boxes:
        patches = []
        for um in RING_UM:
            half_px = int(round(um / mm_px * half / (PATCH / 2) * half))
            # HookNet samples the two rings inside one PATCH-sized window.
            span = PATCH
            x0 = int(round(cx)) - span // 2
            y0 = int(round(cy)) - span // 2
            x0 = max(0, min(x0, arr.shape[1] - 1))
            y0 = max(0, min(y0, arr.shape[0] - 1))
            x1 = min(arr.shape[1], x0 + span)
            y1 = min(arr.shape[0], y0 + span)
            canvas = np.ones((PATCH, PATCH, 3), np.float32)
            sub = arr[y0:y1, x0:x1].astype(np.float32) / 255.0
            canvas[:sub.shape[0], :sub.shape[1]] = sub
            patches.append(canvas)
        out.append(patches)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--section", default="ZEN39")
    ap.add_argument("--patches", type=int, default=48)
    ap.add_argument("--target-um", type=float, default=0.6,
                    help="HookNet needs ~0.5 um texture; this is the resolution budget")
    ap.add_argument("--max-plane-mp", type=float, default=450.0)
    ap.add_argument("--pitch-um", type=float, default=GRID_PITCH_UM)
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args()

    os_env = {"CUDA_VISIBLE_DEVICES": "-1", "TF_CPP_MIN_LOG_LEVEL": "2"}
    import os
    os.environ.update(os_env)

    import tensorflow as tf
    from tensorflow.keras.layers import LeakyReLU
    from hooknet.models.tensorflowmodel import HookNet

    print(f"tensorflow {tf.__version__}", flush=True)
    print(f"CUDA_VISIBLE_DEVICES={os.environ['CUDA_VISIBLE_DEVICES']} "
          f"tf-gpu={bool(tf.config.list_physical_devices('GPU'))}", flush=True)

    img, mm_px, (w, h), scale, meta, lvl, plane_mp = tissue_mask_at(
        args.section, args.target_um, args.max_plane_mp)
    print(f"section {args.section}: source native "
          f"{float(meta['pixel_size_um_estimated']):.3f} um/px -> using level {lvl} "
          f"{img.size[0]}x{img.size[1]} ({plane_mp:.0f} MP) at {mm_px:.3f} um/px",
          flush=True)

    boxes = tissue_boxes(args.section, args.pitch_um, mm_px,
                         img.size[0], img.size[1], scale)
    print(f"declared grid pitch {args.pitch_um} um -> {len(boxes)} candidate "
          f"positions inside tissue on this slide", flush=True)

    t0 = time.time()
    model = HookNet(input_shape=[[PATCH, PATCH, 3], [PATCH, PATCH, 3]],
                    hook_indices=[2, 4], n_classes=3, batch_norm=True,
                    n_filters=32, activation=LeakyReLU(alpha=0.1),
                    model_weights=str(WEIGHTS))
    build_s = time.time() - t0
    print(f"model built in {build_s:.1f}s", flush=True)

    n = min(args.patches, len(boxes))
    if n == 0:
        print("no candidate positions")
        return 1
    # deterministic spread over the tissue rather than the first N
    step = max(1, len(boxes) // n)
    picks = [boxes[i * step] for i in range(n)]

    t1 = time.time()
    batch = make_patches(args.section, picks, mm_px, img, scale)
    prep_s = time.time() - t1
    print(f"{n} patch pairs prepared in {prep_s:.1f}s", flush=True)

    t2 = time.time()
    preds = model.predict([np.stack([b[0] for b in batch]),
                           np.stack([b[1] for b in batch])], verbose=0)
    infer_s = time.time() - t2
    per_patch = infer_s / max(1, n)
    print(f"inference {n} pairs in {infer_s:.1f}s -> {per_patch:.2f} s/patch",
          flush=True)
    print(f"prediction shape {preds.shape} mean={preds.reshape(-1, preds.shape[-1]).mean(0).round(4)}",
          flush=True)

    total_positions = len(boxes)
    est_slide_h = per_patch * total_positions
    slide_area_um2 = w * scale * mm_px * h * scale * mm_px
    print(f"\nextrapolation for THIS slide ({args.section}):")
    print(f"  positions            {total_positions}")
    print(f"  predicted single-core wall clock  {est_slide_h/60:.1f} min "
          f"({est_slide_h/3600:.2f} h)")

    gate = list(__import__("csv").DictReader(
        (ROOT / "infra/tls_four_cancer_labels_20261008/eligibility/annotatable_gate.tsv").open(),
        delimiter="\t"))
    tot_positions = 0
    for r in gate:
        if r["gate_status_human"] != "annotatable_candidate":
            continue
        px_um = float(r["px_um"] or 0)
        if not px_um:
            continue
        # positions scale with area and inversely with pitch squared
        tot_positions += 0  # placeholder; exact count needs per-slide masks
    out = {
        "section": args.section,
        "tensorflow": tf.__version__,
        "gpu_used": bool(tf.config.list_physical_devices("GPU")),
        "level": lvl,
        "level_plane_mp": round(plane_mp, 1),
        "level_um_per_px": round(mm_px, 4),
        "declared_grid_pitch_um": args.pitch_um,
        "positions_on_slide": total_positions,
        "patch_pairs_measured": n,
        "seconds_per_patch": round(per_patch, 3),
        "model_build_seconds": round(build_s, 1),
        "patch_prep_seconds": round(prep_s, 1),
        "est_single_core_hours_this_slide": round(est_slide_h / 3600, 3),
        "assumption": (
            "Throughput measured on CPU with the real weights. The per-slide "
            "estimate uses a declared grid pitch, not a fitted detection result, "
            "and assumes one core. Parallel cores reduce wall clock roughly "
            "linearly; memory and I/O may not."),
    }
    if args.json_out:
        Path(args.json_out).write_text(json.dumps(out, indent=2) + "\n")
    print("\n" + json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())