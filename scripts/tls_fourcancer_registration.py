#!/usr/bin/env python3
"""Visual/quantitative registration check for the four-cancer HEST material.

RNA-free. For each section it asks one question: do the spot coordinates land on
the tissue in the whole-slide image, and under which flip/rotation?

Outputs per section:
  - inside-tissue fraction for all 8 dihedral transforms (identity included)
  - the best transform and its margin over identity
  - centroid offset in micrometres
  - 5 overlay crops at spread-out positions for human confirmation

Operational rule, fixed before reading the numbers:
  misregistered        best inside-fraction < 0.90
  transform_needed     best transform is not identity and beats identity by > 0.05
  ok                   otherwise
"""

from __future__ import annotations

import csv
import json
import warnings
from pathlib import Path

import h5py
import numpy as np
from PIL import Image, ImageDraw

Image.MAX_IMAGE_PIXELS = None
warnings.simplefilter("ignore")

ROOT = Path("/home/huyudi/013_spatial")
HEST = ROOT / "data/other_sources/hest1k"
ELIG = ROOT / "infra/tls_four_cancer_labels_20261008/eligibility"
OUT = ELIG / "registration"
ORGANS = {"Bowel", "Lung", "Breast", "Liver"}
TARGET_LEVEL_WIDTH = 2600
CROP = 760
N_PREVIEW = 5
SPOT_DIAMETER_UM = 55.0

TRANSFORMS = {
    "identity": lambda x, y, w, h: (x, y),
    "flip_lr": lambda x, y, w, h: (w - x, y),
    "flip_ud": lambda x, y, w, h: (x, h - y),
    "rot180": lambda x, y, w, h: (w - x, h - y),
    "transpose": lambda x, y, w, h: (y, x),
    "rot90cw": lambda x, y, w, h: (y, w - x),
    "rot90ccw": lambda x, y, w, h: (h - y, x),
    "antitranspose": lambda x, y, w, h: (h - y, w - x),
}


def load_obs(path: Path):
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        cols = np.asarray(obs["pxl_col_in_fullres"][:], dtype=float)
        rows = np.asarray(obs["pxl_row_in_fullres"][:], dtype=float)
        if "in_tissue" in obs:
            it = np.asarray(obs["in_tissue"][:]).astype(bool)
        else:
            it = np.ones(len(cols), dtype=bool)
    return cols, rows, it


def pick_level(im: Image.Image, width: int, height: int) -> tuple[int, float]:
    n = getattr(im, "n_frames", 1)
    best = 0
    for i in range(n):
        im.seek(i)
        if im.size[0] >= TARGET_LEVEL_WIDTH:
            best = i
    im.seek(best)
    # Downsample factor: fullres pixel -> level pixel. Always <= 1.
    return best, im.size[0] / width


def build_mask(contours: dict, lw: int, lh: int, scale: float) -> np.ndarray:
    """Rasterise tissue contours at the chosen pyramid level.

    Even-odd is approximated by filling outer rings then punching holes black,
    which matches the GeoJSON ring order for these files.
    """
    mask = Image.new("L", (lw, lh), 0)
    draw = ImageDraw.Draw(mask)

    def rings(coords):
        if coords and isinstance(coords[0][0], (int, float)):
            return [coords]
        if coords and isinstance(coords[0][0][0], (int, float)):
            return coords
        return []

    def paint(ring, value):
        pts = [(p[0] * scale, p[1] * scale) for p in ring]
        if len(pts) >= 3:
            draw.polygon(pts, fill=value)

    for feat in contours.get("features", []):
        geom = feat.get("geometry") or {}
        gtype = geom.get("type")
        coords = geom.get("coordinates") or []
        if gtype == "Polygon":
            for i, ring in enumerate(rings(coords)):
                paint(ring, 255 if i == 0 else 0)
        elif gtype == "MultiPolygon":
            for poly in coords:
                for i, ring in enumerate(rings(poly)):
                    paint(ring, 255 if i == 0 else 0)
    return np.asarray(mask) > 127


def inside_fraction(mask: np.ndarray, xs: np.ndarray, ys: np.ndarray,
                    w: int, h: int, name: str, scale: float) -> tuple[float, np.ndarray]:
    fn = TRANSFORMS[name]
    tx, ty = fn(xs, ys, w, h)
    ix = np.rint(tx * scale).astype(np.int64)
    iy = np.rint(ty * scale).astype(np.int64)
    ok = (ix >= 0) & (iy >= 0) & (ix < mask.shape[1]) & (iy < mask.shape[0])
    hit = np.zeros(len(ix), dtype=bool)
    if ok.any():
        hit[ok] = mask[iy[ok], ix[ok]]
    return float(hit.mean()), hit


def make_previews(sid: str, wsi: Path, level: int, scale: float,
                  xs: np.ndarray, ys: np.ndarray, name: str,
                  px_per_um: float, tissue: np.ndarray):
    fn = TRANSFORMS[name]
    out_dir = OUT / "previews" / sid
    out_dir.mkdir(parents=True, exist_ok=True)
    with Image.open(wsi) as im:
        n = getattr(im, "n_frames", 1)
        if level >= n:
            level = n - 1
        im.seek(level)
        lw, lh = im.size
        q = [5, 27, 50, 73, 95][:N_PREVIEW]
        picks = []
        for p in q:
            tx, ty = fn(xs, ys, lw / scale, lh / scale)
            idx = int(np.argsort(tx + ty)[int(len(tx) * p / 100)])
            picks.append((float(tx[idx]), float(ty[idx])))
        written = []
        for k, (cx, cy) in enumerate(picks):
            cx_l, cy_l = cx * scale, cy * scale
            left, top = cx_l - CROP / 2, cy_l - CROP / 2
            box = (int(round(left)), int(round(top)),
                   int(round(left + CROP)), int(round(top + CROP)))
            box = (max(0, box[0]), max(0, box[1]),
                   min(lw, box[2]), min(lh, box[3]))
            if box[2] - box[0] < 32 or box[3] - box[1] < 32:
                continue
            crop = im.crop(box).convert("RGB")
            d = ImageDraw.Draw(crop)
            r = max(2, int(round(SPOT_DIAMETER_UM * px_per_um * scale / 2)))
            for sx, sy in zip(xs, ys):
                tx, ty = fn(np.array([sx]), np.array([sy]), lw / scale, lh / scale)
                px, py = tx[0] * scale - box[0], ty[0] * scale - box[1]
                if r <= px < crop.size[0] and r <= py < crop.size[1]:
                    d.ellipse([px - r, py - r, px + r, py + r],
                              outline=(255, 40, 40), width=2)
            d.rectangle([0, 0, crop.size[0] - 1, 18], fill=(0, 0, 0))
            d.text((4, 4), f"{sid} {name} pos{k + 1} lvl{level}", fill=(255, 255, 255))
            p = out_dir / f"pos{k + 1}.jpg"
            crop.save(p, quality=82, optimize=True)
            written.append(str(p.relative_to(ROOT)))
    return written


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with (HEST / "selected_samples.tsv").open() as fh:
        selected = list(csv.DictReader(fh, delimiter="\t"))
    sections = [r for r in selected
                if r["organ"] in ORGANS and "Visium" in r["st_technology"]]
    sections.sort(key=lambda r: r["id"])
    audit = {r["section_id"]: r for r in
             csv.DictReader((ELIG / "data_audit_sections.tsv").open(), delimiter="\t")}

    results = []
    for row in sections:
        sid = row["id"]
        meta = json.loads((HEST / "metadata" / f"{sid}.json").read_text())
        w = float(meta["fullres_px_width"])
        h = float(meta["fullres_px_height"])
        est = float(meta["pixel_size_um_estimated"])
        cols, rows_, it = load_obs(HEST / "st" / f"{sid}.h5ad")
        keep = it if it.sum() > 50 else np.ones(len(cols), dtype=bool)
        xs, ys = cols[keep], rows_[keep]
        contours = json.loads((HEST / "tissue_seg" / f"{sid}_contours.geojson").read_text())

        wsi = HEST / "wsis" / f"{sid}.tif"
        with Image.open(wsi) as im:
            level, scale = pick_level(im, w, h)
            im.seek(level)
            lw, lh = im.size
            mask = build_mask(contours, lw, lh, scale)
            tissue_area_frac = float(mask.mean())
            scores = {}
            hits = {}
            for name in TRANSFORMS:
                frac, hit = inside_fraction(mask, xs, ys, w, h, name, scale)
                scores[name] = round(frac, 4)
                hits[name] = hit
            best = max(scores, key=scores.get)
            identity = scores["identity"]
            margin = round(scores[best] - identity, 4)

            if best == "identity" or margin <= 0.05:
                verdict = "ok" if scores[best] >= 0.90 else "misregistered"
            else:
                verdict = "transform_needed"

            ys_i, xs_i = np.nonzero(mask)
            if len(xs_i) and scores[best] > 0:
                tx, ty = TRANSFORMS[best](xs, ys, w, h)
                hit = hits[best]
                off_x = (tx[hit].mean() - xs_i.mean()) * est
                off_y = (ty[hit].mean() - ys_i.mean()) * est
            else:
                off_x = off_y = float("nan")

            previews = []
            if scores[best] >= 0.90:
                px_per_um = 1.0 / est
                previews = make_previews(sid, wsi, level, scale, xs, ys, best,
                                         px_per_um, mask)

        rec = {
            "section_id": sid,
            "organ": row["organ"],
            "platform": row["st_technology"],
            "disease_state": row["disease_state"],
            "n_spots_used": int(len(xs)),
            "image_w": int(w), "image_h": int(h),
            "level": level, "level_w": lw, "level_h": lh,
            "tissue_area_frac": round(tissue_area_frac, 5),
            "identity_inside": identity,
            "best_transform": best,
            "best_inside": scores[best],
            "margin_over_identity": margin,
            "offset_um_x": round(float(off_x), 1),
            "offset_um_y": round(float(off_y), 1),
            "verdict": verdict,
            "n_previews": len(previews),
            "preview_paths": ";".join(previews),
            **{f"inside_{k}": v for k, v in scores.items()},
        }
        results.append(rec)
        print(f'{sid} {row["organ"]:<7} best={best:<12} {scores[best]:.4f} '
              f'id={identity:.4f} margin={margin:+.4f} -> {verdict}', flush=True)

    with (OUT / "registration_check.tsv").open("w", newline="") as fh:
        keys = list(results[0])
        wtr = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        wtr.writeheader()
        wtr.writerows(results)

    from collections import Counter
    summary = {
        "n_sections": len(results),
        "verdicts": dict(Counter(r["verdict"] for r in results)),
        "n_previews_written": sum(r["n_previews"] for r in results),
        "transform_needed": [r["section_id"] for r in results
                             if r["verdict"] == "transform_needed"],
        "misregistered": [r["section_id"] for r in results
                          if r["verdict"] == "misregistered"],
        "ok_sections": [r["section_id"] for r in results if r["verdict"] == "ok"],
        "rule": "misregistered if best<0.90; transform_needed if best!=identity and margin>0.05",
        "rna_used": "no",
        "labels_created": False,
    }
    (OUT / "registration_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items()
                      if k not in ("ok_sections",)}, indent=2))


if __name__ == "__main__":
    main()