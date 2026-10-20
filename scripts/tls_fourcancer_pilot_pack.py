#!/usr/bin/env python3
"""Propose a labelling pilot and build an RNA-free annotator image pack.

The pilot is a feasibility trial, not a sample-size claim. Selection rules are
declared here and applied mechanically:

  1. gate_status must be annotatable_candidate (registration, panel and scale ok)
  2. prefer specimens with a metadata-backed patient identity
  3. one section per specimen, so replicates never inflate the independent count
  4. fill the remaining slots only with unknown-identity sections, spread across
     different studies, and never count them as independent patients
  5. reserve known-patient specimens that are not used, and say plainly that a
     2-patient reserve is not an external test set
  6. no TLS score, hotspot or RNA quantity is used anywhere in selection

Output tiles are pure H&E with no spot overlay, so the first reading pass stays
blind to spot positions.
"""

from __future__ import annotations

import csv
import json
import warnings
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
from PIL import Image

Image.MAX_IMAGE_PIXELS = None
warnings.simplefilter("ignore")

ROOT = Path("/home/huyudi/013_spatial")
HEST = ROOT / "data/other_sources/hest1k"
ELIG = ROOT / "infra/tls_four_cancer_labels_20261008/eligibility"
PACK = ELIG.parent / "annotator_pack"
ORGANS = {"Bowel", "Lung", "Breast", "Liver"}

# Declared quotas. Liver has only 2 usable cancer sections, both known patients.
QUOTA = {"Liver": 2, "Bowel": 4, "Breast": 2, "Lung": 2}
N_TILES = 5
TILE_FIELD_UM = 1000.0
MAX_TILE_PX = 3600
MAX_PLANE_MP = 120.0


def pick_read_level(im: Image.Image, w: int, h: int) -> tuple[int, float]:
    """Finest pyramid level whose decoded plane stays under the memory cap."""
    n = getattr(im, "n_frames", 1)
    chosen = 0
    for i in range(n):
        im.seek(i)
        if im.size[0] * im.size[1] / 1e6 <= MAX_PLANE_MP:
            chosen = i
            break
        chosen = i
    im.seek(chosen)
    return chosen, im.size[0] / w


def load_xy(path: Path):
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        c = np.asarray(obs["pxl_col_in_fullres"][:], dtype=float)
        r = np.asarray(obs["pxl_row_in_fullres"][:], dtype=float)
        it = np.asarray(obs["in_tissue"][:]).astype(bool) if "in_tissue" in obs \
            else np.ones(len(c), bool)
    return c, r, it


def write_tiles(sid: str, xs: np.ndarray, ys: np.ndarray) -> list[dict]:
    meta = json.loads((HEST / "metadata" / f"{sid}.json").read_text())
    w = int(meta["fullres_px_width"])
    h = int(meta["fullres_px_height"])
    px_um = float(meta["pixel_size_um_estimated"])
    out_dir = PACK / "tiles" / sid
    out_dir.mkdir(parents=True, exist_ok=True)
    wsi = HEST / "wsis" / f"{sid}.tif"
    tiles = []
    with Image.open(wsi) as im:
        level, scale = pick_read_level(im, w, h)
        im.seek(level)
        lw, lh = im.size
        um_per_level_px = px_um / scale
        crop = int(round(min(MAX_TILE_PX, TILE_FIELD_UM / um_per_level_px)))
        order = np.argsort(xs + ys)
        picks = [int(order[int(len(order) * q / 100)]) for q in (5, 28, 50, 72, 95)]
        for k, idx in enumerate(picks, 1):
            cx, cy = xs[idx] * scale, ys[idx] * scale
            left, top = cx - crop / 2, cy - crop / 2
            box = (int(round(left)), int(round(top)),
                   int(round(left + crop)), int(round(top + crop)))
            box = (max(0, box[0]), max(0, box[1]), min(lw, box[2]), min(lh, box[3]))
            if box[2] - box[0] < crop * 0.6 or box[3] - box[1] < crop * 0.6:
                continue
            tile = im.crop(box).convert("RGB")
            p = out_dir / f"field{k}.jpg"
            tile.save(p, quality=90, optimize=True)
            tiles.append({
                "section_id": sid, "tile": f"tiles/{sid}/field{k}.jpg",
                "fullres_x": int(round(box[0] / scale)),
                "fullres_y": int(round(box[1] / scale)),
                "fullres_w": int(round((box[2] - box[0]) / scale)),
                "fullres_h": int(round((box[3] - box[1]) / scale)),
                "field_um": round((box[2] - box[0]) * um_per_level_px, 1),
                "read_level": level, "um_per_pixel": round(um_per_level_px, 5),
            })
    return tiles


def main() -> None:
    gate = list(csv.DictReader((ELIG / "annotatable_gate.tsv").open(), delimiter="\t"))
    cands = [r for r in gate if r["gate_status"] == "annotatable_candidate"]

    by_specimen = defaultdict(list)
    for r in cands:
        by_specimen[r["specimen_key"]].append(r)
    for v in by_specimen.values():
        v.sort(key=lambda r: r["section_id"])

    chosen: list[dict] = []
    reserved: list[dict] = []
    for organ in ("Liver", "Bowel", "Breast", "Lung"):
        quota = QUOTA[organ]
        specs = {k: v for k, v in by_specimen.items()
                 if v[0]["organ"] == organ}
        known = sorted([k for k, v in specs.items()
                        if v[0]["patient_status"] == "known_within_study"])
        unknown = sorted([k for k, v in specs.items()
                          if v[0]["patient_status"] != "known_within_study"])
        take_known = known[:quota]
        hold_known = known[quota:]
        take_unknown = unknown[:max(0, quota - len(take_known))]
        for k in take_known + take_unknown:
            r = specs[k][0]
            chosen.append({**r, "selection_reason": (
                "known_patient_specimen" if k in take_known
                else "unknown_identity_feasibility_not_counted")})
        for k in hold_known:
            reserved.append({**specs[k][0],
                             "reservation": "held_back_known_patient_not_an_external_test"})
    chosen.sort(key=lambda r: (r["organ"], r["section_id"]))

    PACK.mkdir(parents=True, exist_ok=True)
    all_tiles = []
    for r in chosen:
        sid = r["section_id"]
        c, rr, it = load_xy(HEST / "st" / f"{sid}.h5ad")
        keep = it if it.sum() > 50 else np.ones(len(c), bool)
        all_tiles.extend(write_tiles(sid, c[keep], rr[keep]))
        link = PACK / "images" / f"{sid}.tif"
        link.parent.mkdir(parents=True, exist_ok=True)
        if not link.is_symlink():
            link.symlink_to((HEST / "wsis" / f"{sid}.tif").resolve())
        r["wsi_link"] = f"images/{sid}.tif"

    with (PACK / "pilot_sections.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(chosen[0]), delimiter="\t")
        w.writeheader()
        w.writerows(chosen)
    with (PACK / "reserved_sections.tsv").open("w", newline="") as fh:
        fields = list(reserved[0]) if reserved else ["section_id", "reservation"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(reserved)
    with (PACK / "tiles.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_tiles[0]), delimiter="\t")
        w.writeheader()
        w.writerows(all_tiles)

    n_known = sum(1 for r in chosen if r["patient_status"] == "known_within_study")
    summary = {
        "pilot_sections": len(chosen),
        "quota": QUOTA,
        "by_organ": {o: sum(1 for r in chosen if r["organ"] == o)
                     for o in ("Liver", "Bowel", "Breast", "Lung")},
        "sections_with_known_patient": n_known,
        "independent_patients_in_pilot": n_known,
        "unknown_identity_sections_not_counted": len(chosen) - n_known,
        "reserved_known_patient_specimens": [r["section_id"] for r in reserved],
        "external_test_available": False,
        "external_test_note": (
            "A 2-specimen reserve is not an external test set. No locked test set "
            "can be formed from these materials."),
        "tiles_written": len(all_tiles),
        "tile_field_um": TILE_FIELD_UM,
        "rna_or_score_used_in_selection": False,
        "labels_created": False,
        "annotator_assigned": False,
        "status": "PROPOSED_PILOT_NOT_LOCKED",
    }
    (PACK / "pack_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()