#!/usr/bin/env python3
"""Test whether near-threshold sections are recoverable by a translation fix.

A section whose spots sit slightly off its tissue contour may reflect a global
frame offset rather than a bad section. A single translation is only believable
if the optimum is a broad basin: a spike means we are fitting the contour, not
measuring an offset.

Acceptance rule, fixed before looking at any result:
  R1  best inside-fraction >= 0.95
  R2  gain over identity >= 0.05
  R3  the optimum is a plateau: >= 2% of the evaluated grid lies within 0.5
      percentage points of the best value
  R4  |dx|, |dy| <= 1000 px and <= 5% of image width/height
Nothing is applied automatically; this only reports what could be recovered.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import warnings
from pathlib import Path

import h5py
import numpy as np

warnings.filterwarnings("ignore")
from PIL import Image, ImageDraw

Image.MAX_IMAGE_PIXELS = None

ROOT = Path("/home/huyudi/013_spatial")
HEST = ROOT / "data/other_sources/hest1k"
ELIG = ROOT / "infra/tls_four_cancer_labels_20261008/eligibility"
OUT = ELIG / "registration_repair"

LEVEL_TARGET_W = 2000
MAX_DXY = 1000
COARSE_STEP = 40
FINE_STEP = 8
PLATEAU_FRAC = 0.02
PLATEAU_TOL = 0.005


def load_sections() -> list[dict]:
    path = ELIG / "annotatable_gate.tsv"
    out = []
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["gate_status"] != "blocked_registration":
                continue
            v = float(r["spot_inside_tissue"] or 0)
            if v >= 0.80:                      # near-threshold only
                out.append(r)
    return sorted(out, key=lambda r: float(r["spot_inside_tissue"]))


def tissue_mask(sid: str, px_um: float) -> tuple[np.ndarray, float, int, int]:
    """Boolean tissue mask at a reduced level, plus its downsample factor."""
    meta = json.loads((HEST / "metadata" / f"{sid}.json").read_text())
    w, h = int(meta["fullres_px_width"]), int(meta["fullres_px_height"])
    with Image.open(HEST / "wsis" / f"{sid}.tif") as im:
        n = getattr(im, "n_frames", 1)
        lvl = 0
        for i in range(n):
            im.seek(i)
            if im.size[0] <= LEVEL_TARGET_W:
                lvl = i
        im.seek(lvl)
        lw, lh = im.size
    scale = lw / w                                   # <= 1
    mask = Image.new("L", (lw, lh), 0)
    draw = ImageDraw.Draw(mask)
    contours = json.loads((HEST / "tissue_seg" / f"{sid}_contours.geojson").read_text())

    def paint(ring, value):
        pts = [(p[0] * scale, p[1] * scale) for p in ring]
        if len(pts) >= 3:
            draw.polygon(pts, fill=value)

    for feat in contours.get("features", []):
        geom = feat.get("geometry") or {}
        coords = geom.get("coordinates") or []
        if geom.get("type") == "Polygon":
            rings = coords
        elif geom.get("type") == "MultiPolygon":
            rings = [r for poly in coords for r in poly]
        else:
            continue
        for i, ring in enumerate(rings):
            paint(ring, 255 if i == 0 else 0)
    return np.asarray(mask) > 127, scale, w, h


def spot_xy(sid: str, scale: float):
    with h5py.File(HEST / "st" / f"{sid}.h5ad", "r") as f:
        obs = f["obs"]
        c = np.asarray(obs["pxl_col_in_fullres"][:], dtype=np.int64)
        r = np.asarray(obs["pxl_row_in_fullres"][:], dtype=np.int64)
        it = (np.asarray(obs["in_tissue"][:]).astype(bool) if "in_tissue" in obs
              else np.ones(len(c), bool))
    return c[it] * scale, r[it] * scale


def score(mask: np.ndarray, xs: np.ndarray, ys: np.ndarray, dx: int, dy: int) -> float:
    """Fraction of spots inside the mask after shifting spots by (dx, dy)."""
    ix = np.rint(xs + dx).astype(np.int64)
    iy = np.rint(ys + dy).astype(np.int64)
    ok = (ix >= 0) & (iy >= 0) & (ix < mask.shape[1]) & (iy < mask.shape[0])
    if not ok.any():
        return 0.0
    return float(mask[iy[ok], ix[ok]].mean())


def search(mask, xs, ys, step, limit):
    best = (-1.0, 0, 0)
    grid = []
    for dy in range(-limit, limit + 1, step):
        for dx in range(-limit, limit + 1, step):
            v = score(mask, xs, ys, dx, dy)
            grid.append(v)
            if v > best[0]:
                best = (v, dx, dy)
    return best, grid


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    rows = load_sections()
    if args.only:
        rows = [r for r in rows if r["section_id"] in args.only]
    print(f"near-threshold registration-blocked sections: {len(rows)}", flush=True)

    results = []
    for i, r in enumerate(rows, 1):
        sid = r["section_id"]
        t0 = time.time()
        mask, scale, w, h = tissue_mask(sid, 1.0)
        xs, ys = spot_xy(sid, scale)
        ident = score(mask, xs, ys, 0, 0)
        (bv, bdx, bdy), coarse = search(mask, xs, ys, COARSE_STEP, MAX_DXY)
        (fv, fdx, fdy), fine = search(mask, xs, ys, FINE_STEP,
                                      COARSE_STEP)
        allv = coarse + fine
        tol = PLATEAU_TOL
        plateau = sum(1 for v in allv if v >= fv - tol) / len(allv)
        dx_px, dy_px = fdx / scale, fdy / scale
        r1 = fv >= 0.95
        r2 = (fv - ident) >= 0.05
        r3 = plateau >= PLATEAU_FRAC
        r4 = (abs(dx_px) <= MAX_DXY and abs(dy_px) <= MAX_DXY
              and abs(dx_px) <= 0.05 * w and abs(dy_px) <= 0.05 * h)
        rec = {
            "section_id": sid, "organ": r["organ"],
            "identity_inside": round(ident, 4),
            "best_inside": round(fv, 4),
            "gain": round(fv - ident, 4),
            "dx_px": round(dx_px, 1), "dy_px": round(dy_px, 1),
            "plateau_frac": round(plateau, 4),
            "R1_best>=0.95": "yes" if r1 else "no",
            "R2_gain>=0.05": "yes" if r2 else "no",
            "R3_plateau": "yes" if r3 else "no",
            "R4_translation_sane": "yes" if r4 else "no",
            "recoverable": "yes" if (r1 and r2 and r3 and r4) else "no",
            "seconds": round(time.time() - t0, 1),
        }
        results.append(rec)
        print(f'  [{i}/{len(rows)}] {sid} ident={ident:.4f} best={fv:.4f} '
              f'gain={fv-ident:+.4f} shift=({dx_px:+.0f},{dy_px:+.0f})px '
              f'plateau={plateau:.3f} -> {rec["recoverable"]}', flush=True)

    with (OUT / "repair_report.tsv").open("w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(results[0]), delimiter="\t")
        wtr.writeheader()
        wtr.writerows(results)
    ok = [r for r in results if r["recoverable"] == "yes"]
    print(f"\nrecoverable under the declared rule: {len(ok)}/{len(results)}")
    for r in ok:
        print(f'   {r["section_id"]} gain={r["gain"]} shift=({r["dx_px"]},{r["dy_px"]})')
    print(f"report: {OUT / 'repair_report.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())