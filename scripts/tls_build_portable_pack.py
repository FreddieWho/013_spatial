#!/usr/bin/env python3
"""Build the portable, offline TLS annotation pack.

Renders a JPEG tile pyramid per section from the HEST pyramid TIFF, copies the
viewer and the frozen rules, and writes an index. No annotation data is created.

    python3 scripts/tls_build_portable_pack.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

ROOT = Path("/home/huyudi/013_spatial")
HEST = ROOT / "data/other_sources/hest1k"
SRC = ROOT / "infra/tls_annotation_pack_src"
PACK = ROOT / "infra/tls_four_cancer_labels_20261008/annotator_pack"
DIST = ROOT / "dist"
TILE = 256
MIN_TILE = 4          # stop descending when the level is this small


def build_pyramid(sid: str, out_dir: Path, target_um: float, quality: int,
                 tol: float = 1.08) -> dict:
    meta = json.loads((HEST / "metadata" / f"{sid}.json").read_text())
    full_w = int(meta["fullres_px_width"])
    full_h = int(meta["fullres_px_height"])
    px_um = float(meta["pixel_size_um_estimated"])

    t0 = time.time()
    im = Image.open(HEST / "wsis" / f"{sid}.tif")
    n = getattr(im, "n_frames", 1)
    # Levels get coarser with rising index; take the coarsest one that is still
    # within tolerance of the target. A hard <= comparison traps slides whose
    # levels land a hair above the target (0.2757 um/px slides step
    # 0.551 -> 1.103, and 1.103 > 1.10), which costs 4x the pixels.
    base = 0
    for i in range(n):
        im.seek(i)
        lw, _ = im.size
        if (px_um / (lw / full_w)) <= target_um * tol:
            base = i
    im.seek(base)
    img = im.convert("RGB")
    base_w, base_h = img.size
    base_um = px_um / (base_w / full_w)
    im.close()

    levels = []
    L = 0
    cur = img
    while True:
        w, h = cur.size
        d = 1 << L
        levels.append({
            "level": L, "w": w, "h": h,
            "um_per_px": round(base_um * d, 5),
        })
        ld = out_dir / str(L)
        ld.mkdir(parents=True, exist_ok=True)
        for row in range((h + TILE - 1) // TILE):
            for col in range((w + TILE - 1) // TILE):
                box = (col * TILE, row * TILE,
                       min((col + 1) * TILE, w), min((row + 1) * TILE, h))
                tile = cur.crop(box)
                if tile.size[0] < 8 or tile.size[1] < 8:
                    tile = tile.resize((max(tile.size[0], 8), max(tile.size[1], 8)))
                tile.save(ld / f"{col}_{row}.jpg", quality=quality,
                          optimize=False, subsampling=0)
        if max(w, h) <= TILE * MIN_TILE:
            break
        L += 1
        cur = cur.reduce(2)
        if max(cur.size) > 6000:            # guard against a runaway pyramid
            cur = cur.resize((cur.size[0] // 2, cur.size[1] // 2))

    out = {
        "section_id": sid,
        "organ": meta.get("organ"),
        "disease_state": meta.get("disease_state"),
        "full_w": full_w, "full_h": full_h,
        "native_um_per_px": round(px_um, 5),
        "base_level": base, "base_w": base_w, "base_h": base_h,
        "base_um_per_px": round(base_um, 5),
        "tile": TILE,
        "levels": levels,
        "coords": "base image pixels; multiply by full_w/base_w for original pixels",
        "build_seconds": round(time.time() - t0, 1),
    }
    (out_dir / "meta.json").write_text(json.dumps(out, indent=1) + "\n")

    tissue = HEST / "tissue_seg" / f"{sid}_contours.geojson"
    if tissue.exists():
        shutil.copy2(tissue, out_dir / "tissue.geojson")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-um", type=float, default=1.10,
                    help="finest delivered resolution in micrometres per pixel")
    ap.add_argument("--quality", type=int, default=82)
    ap.add_argument("--out", default=str(DIST / "tls_annotation_pack"))
    ap.add_argument("--sections", nargs="*", default=None)
    ap.add_argument("--include-coarse", action="store_true",
                    help="also include sections that fail only the 0.6 um/px "
                         "resolution gate (still readable for human annotation)")
    ap.add_argument("--two-tier", dest="two_tier", action="store_true", default=True,
                    help="fine resolution for the pilot and every known-patient "
                         "section, coarser for the rest")
    ap.add_argument("--no-two-tier", dest="two_tier", action="store_false")
    ap.add_argument("--add-only", action="store_true",
                    help="build only --sections into an existing pack and merge "
                         "the index, instead of rebuilding everything")
    ap.add_argument("--human-gate", action="store_true",
                    help="select by the human-legibility gate instead of the "
                         "automated-pipeline gate")
    ap.add_argument("--keep-duplicates", action="store_true",
                    help="keep sections whose whole-slide image is byte-identical "
                         "to an already-selected one")
    ap.add_argument("--fine-um", type=float, default=1.10)
    ap.add_argument("--bulk-um", type=float, default=2.20)
    ap.add_argument("--quality-bulk", type=int, default=None)
    args = ap.parse_args()
    if args.quality_bulk is None:
        args.quality_bulk = args.quality

    import hashlib

    def wsi_md5(sid: str) -> str:
        h = hashlib.md5()
        with (HEST / "wsis" / f"{sid}.tif").open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 22), b""):
                h.update(chunk)
        return h.hexdigest()
    out_root = Path(args.out)
    if out_root.exists() and not args.add_only:
        shutil.rmtree(out_root)
    (out_root / "data").mkdir(parents=True, exist_ok=True)

    # viewer + scripts + docs
    shutil.copytree(SRC / "viewer", out_root / "viewer", dirs_exist_ok=True)
    for f in ("run.sh", "run.bat", "README.md", "CHECKLIST.md"):
        shutil.copy2(SRC / f, out_root / f)
    (out_root / "run.sh").chmod(0o755)
    (out_root / "docs").mkdir(exist_ok=True)
    guide_md = (ROOT / "docs/ANNOTATION_GUIDE_20261008.md").read_text()
    guide_md = guide_md.replace(
        "scripts/tls_annotate_gui.py", "../viewer/index.html（离线版）")
    # Ship HTML, not raw markdown: http.server sends no charset for .md, which
    # makes browsers decode the Chinese as something else (mojibake).
    subprocess.run(
        [sys.executable, str(ROOT / "scripts/md_to_html.py"),
         str(ROOT / "docs/ANNOTATION_GUIDE_20261008.md"),
         str(out_root / "docs/ANNOTATION_GUIDE.html"),
         "--title", "TLS 判读规则"], check=True)
    subprocess.run(
        [sys.executable, str(ROOT / "scripts/md_to_html.py"),
         str(SRC / "CHECKLIST.md"), str(out_root / "docs/CHECKLIST.html"),
         "--title", "每张自查表"], check=True)
    shutil.copy2(ROOT / "docs/ANNOTATION_GUIDE_20261008.md",
                 out_root / "docs/ANNOTATION_GUIDE.md")

    # The annotation set is every section that passed the annotatable gate, not
    # just the pilot. Role is recorded so the pilot stays distinguishable.
    gate_path = (ROOT / "infra/tls_four_cancer_labels_20261008/eligibility"
                 / "annotatable_gate.tsv")
    pilot_ids = {r["section_id"] for r in
                 csv.DictReader((PACK / "pilot_sections.tsv").open(), delimiter="\t")}
    gate_col = "gate_status_human" if args.human_gate else "gate_status"
    wanted = {"annotatable_candidate"}
    if args.human_gate:
        # the human gate already folds coarse-only and panel-blocked sections in
        wanted = {"annotatable_candidate"}
    rows = [r for r in csv.DictReader(gate_path.open(), delimiter="\t")
            if r.get(gate_col) in wanted]
    if args.sections:
        rows = [r for r in rows if r["section_id"] in args.sections]
    rows.sort(key=lambda r: (r["organ"], r["section_id"]))

    index = []
    if args.add_only:
        prev = json.loads((out_root / "data/index.json").read_text())
        index = [s for s in prev["slides"] if s["section_id"] not in
                 {r["section_id"] for r in rows}]

    seen_md5: dict[str, str] = {}
    for s in index:
        seen_md5[s.get("wsi_md5", "")] = s["section_id"]
    for s in index:
        pass
    for p in out_root.glob("data/*/meta.json"):
        pass

    skipped = []
    for i, r in enumerate(rows, 1):
        sid = r["section_id"]
        digest = wsi_md5(sid)
        dup = seen_md5.get(digest)
        if dup and not args.keep_duplicates:
            skipped.append({"section_id": sid, "duplicate_of": dup,
                            "wsi_md5": digest})
            print(f'[{i}/{len(rows)}] {sid} SKIPPED - byte-identical to {dup}', flush=True)
            continue
        seen_md5[digest] = sid
        fine = args.two_tier and (
            sid in pilot_ids or r["patient_status"] == "known_within_study")
        tgt = args.fine_um if fine else args.bulk_um
        qual = args.quality if fine else args.quality_bulk
        tier = "fine" if fine else "bulk"
        print(f'[{i}/{len(rows)}] {sid} {tier} @<={tgt}um ...', flush=True)
        m = build_pyramid(sid, out_root / "data" / sid, tgt, qual)
        meta = json.loads((HEST / "metadata" / f"{sid}.json").read_text())
        m["gate_status"] = r["gate_status"]
        m["gate_status_human"] = r.get("gate_status_human", "")
        m["gate_reason_human"] = r.get("gate_reason_human", "")
        m["role"] = "pilot" if sid in pilot_ids else "expanded"
        m["tier"] = tier
        m["target_um"] = tgt
        m["organ"] = r["organ"]
        m["disease_state"] = r["disease_state"]
        m["disease_comment"] = str(meta.get("disease_comment") or "")
        m["preservation_method"] = meta.get("preservation_method")
        m["wsi_md5"] = digest
        m["duplicate_of"] = ""
        m["patient_status"] = r["patient_status"]
        m["patient_key"] = r["patient_key"]
        m["specimen_key"] = r["specimen_key"]
        m["spot_inside_tissue"] = r["spot_inside_tissue"]
        m["registration_verdict"] = r["registration_verdict"]
        m["tenx_possible_same_release"] = r["tenx_possible_same_release"]
        index.append(m)
        print(f'    {m["base_w"]}x{m["base_h"]} @ {m["base_um_per_px"]} µm/px, '
              f'{len(m["levels"])} levels, {m["build_seconds"]}s', flush=True)

    doc = {
        "generated": time.strftime("%Y-%m-%d"),
        "target_um_per_px": args.target_um,
        "jpeg_quality": args.quality,
        "tile": TILE,
        "selection": ("annotatable_candidates" + ("+coarse" if args.include_coarse else "")),
        "resolution_tiers": ({"fine_um": args.fine_um, "bulk_um": args.bulk_um,
                              "two_tier": bool(args.two_tier)} if args.two_tier
                             else {"single_um": args.target_um, "two_tier": False}),
        "n_pilot": sum(1 for m in index if m["role"] == "pilot"),
        "n_expanded": sum(1 for m in index if m["role"] == "expanded"),
        "n_fine": sum(1 for m in index if m.get("tier") == "fine"),
        "reviewer_credential": "expert_human_review",
        "double_read": False,
        "annotation_protocol": "blind_first_pass_single_reader",
        "contains_labels": False,
        "skipped_duplicates": skipped,
        "slides": index,
    }
    (out_root / "data" / "index.json").write_text(json.dumps(doc, indent=1) + "\n")

    total = sum(p.stat().st_size for p in out_root.rglob("*") if p.is_file())
    print(f"\npack: {out_root}")
    print(f"size: {total/1e6:.0f} MB  sections: {len(index)}")


if __name__ == "__main__":
    sys.exit(main())