#!/usr/bin/env python3
"""Validate annotations submitted for the four-cancer TLS pilot.

Checks format, coordinate frame and coverage completeness. It does not judge
whether a structure is really a TLS, and it does not map anything to spots yet.

    python3 scripts/tls_annotation_check.py
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

from shapely.geometry import Polygon, shape
from shapely.ops import unary_union

ROOT = Path("/home/huyudi/013_spatial")
HEST = ROOT / "data/other_sources/hest1k"
PACK = ROOT / "infra/tls_four_cancer_labels_20261008/annotator_pack"
ANN = ROOT / "infra/tls_four_cancer_labels_20261008/annotation/incoming"
OUT = ROOT / "infra/tls_four_cancer_labels_20261008/annotation"
ALLOWED_CLASS = {"TLS", "lymphoid_aggregate", "non_TLS", "unresolved"}
ALLOWED_CONF = {"high", "medium", "low", "unknown"}
ALLOWED_MAT = {"GC_present", "GC_not_identified", "unknown"}
BORDER_TOL_UM = 200.0


def tissue_area_um2(sid: str, px_um: float) -> float:
    p = HEST / "tissue_seg" / f"{sid}_contours.geojson"
    if not p.exists():
        return float("nan")
    geoms = []
    for f in json.loads(p.read_text()).get("features", []):
        g = shape(f["geometry"])
        if not g.is_valid:
            g = g.buffer(0)
        geoms.append(g)
    if not geoms:
        return float("nan")
    return unary_union(geoms).area * px_um * px_um


def check_section(sid: str) -> dict:
    meta = json.loads((HEST / "metadata" / f"{sid}.json").read_text())
    w = float(meta["fullres_px_width"])
    h = float(meta["fullres_px_height"])
    px_um = float(meta["pixel_size_um_estimated"])
    rec = {"section_id": sid, "image_w": int(w), "image_h": int(h),
           "um_per_px": round(px_um, 4), "errors": [], "warnings": []}

    for kind, fname in (("instances", "instances.geojson"),
                        ("reviewed", "reviewed_regions.geojson")):
        p = ANN / sid / fname
        rec[f"{kind}_file"] = "present" if p.exists() else "MISSING"
        if not p.exists():
            rec[f"n_{kind}"] = 0
            rec[f"area_{kind}_mm2"] = 0.0
            rec["errors"].append(f"{fname} missing")
            continue
        try:
            doc = json.loads(p.read_text())
        except json.JSONDecodeError as exc:
            rec[f"n_{kind}"] = 0
            rec[f"area_{kind}_mm2"] = 0.0
            rec["errors"].append(f"{fname} not valid JSON: {exc}")
            continue
        feats = doc.get("features") or []
        rec[f"n_{kind}"] = len(feats)
        polys = []
        ids = []
        for i, f in enumerate(feats):
            geom = f.get("geometry") or {}
            if geom.get("type") != "Polygon":
                rec["errors"].append(f"{fname}[{i}] geometry is {geom.get('type')}, not Polygon")
                continue
            ring = geom["coordinates"][0]
            if len(ring) < 3:
                rec["errors"].append(f"{fname}[{i}] has {len(ring)} points")
                continue
            g = shape(geom)
            if not g.is_valid:
                rec["warnings"].append(f"{fname}[{i}] self-intersecting")
                g = g.buffer(0)
            polys.append(g)
            xs = [pt[0] for pt in ring]
            ys = [pt[1] for pt in ring]
            tol = BORDER_TOL_UM / px_um
            if min(xs) < -tol or min(ys) < -tol or max(xs) > w + tol or max(ys) > h + tol:
                rec["warnings"].append(
                    f"{fname}[{i}] outside image (x {min(xs):.0f}-{max(xs):.0f}, "
                    f"y {min(ys):.0f}-{max(ys):.0f}, image {int(w)}x{int(h)})")
            props = f.get("properties") or {}
            if kind == "instances":
                ids.append(props.get("instance_id"))
                cls = props.get("structure_class")
                if cls not in ALLOWED_CLASS:
                    rec["errors"].append(f"instances[{i}] structure_class={cls!r}")
                if props.get("confidence") not in ALLOWED_CONF:
                    rec["errors"].append(f"instances[{i}] confidence={props.get('confidence')!r}")
                if props.get("maturity") not in ALLOWED_MAT:
                    rec["errors"].append(f"instances[{i}] maturity={props.get('maturity')!r}")
                if props.get("evidence_level") != "independent_reference":
                    rec["errors"].append(f"instances[{i}] evidence_level="
                                         f"{props.get('evidence_level')!r}")
                if props.get("double_read") != "false":
                    rec["warnings"].append(f"instances[{i}] double_read="
                                           f"{props.get('double_read')!r}")
        if kind == "instances":
            dup = [k for k, v in Counter(ids).items() if v > 1]
            if dup:
                rec["errors"].append(f"duplicate instance_id: {dup[:3]}")
        area = unary_union(polys).area * px_um * px_um / 1e6 if polys else 0.0
        rec[f"area_{kind}_mm2"] = round(area, 3)

    ta = tissue_area_um2(sid, px_um) / 1e6
    rec["tissue_area_mm2"] = round(ta, 3) if ta == ta else ""
    if rec["tissue_area_mm2"] and rec["area_reviewed_mm2"]:
        rec["reviewed_over_tissue"] = round(
            rec["area_reviewed_mm2"] / rec["tissue_area_mm2"], 3)
        if rec["reviewed_over_tissue"] < 0.5:
            rec["warnings"].append(
                f"reviewed covers only {rec['reviewed_over_tissue']:.0%} of the "
                f"tissue contour; check the coverage pass is finished")
        if rec["reviewed_over_tissue"] > 1.6:
            rec["warnings"].append(
                f"reviewed area is {rec['reviewed_over_tissue']:.0%} of the tissue "
                f"contour; some reviewed polygons may sit on background")
    if rec["n_instances"] == 0 and rec["area_reviewed_mm2"] > 0:
        rec["warnings"].append("coverage drawn but no instance marked")
    rec["errors"] = "; ".join(rec["errors"])[:400]
    rec["warnings"] = "; ".join(rec["warnings"])[:400]
    rec["status"] = ("ERROR" if rec["errors"]
                     else "INCOMPLETE" if rec["warnings"] or rec["n_instances"] == 0
                     else "OK")
    return rec


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with (PACK / "pilot_sections.tsv").open() as fh:
        sections = [r["section_id"] for r in csv.DictReader(fh, delimiter="\t")]
    rows = [check_section(s) for s in sections]
    fields = ["section_id", "status", "n_instances", "n_reviewed",
              "area_instances_mm2", "area_reviewed_mm2", "tissue_area_mm2",
              "reviewed_over_tissue", "image_w", "image_h", "um_per_px",
              "errors", "warnings"]
    with (OUT / "annotation_check.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    counts = Counter(r["status"] for r in rows)
    print(f"sections {len(rows)}  {dict(counts)}")
    missing = [r for r in rows if r["instances_file"] == "MISSING"]
    for r in rows:
        if r["status"] != "OK":
            print(f'  {r["section_id"]:<9} {r["status"]:<11} '
                  f'inst={r["n_instances"]:<4} rev={r["n_reviewed"]:<4} '
                  f'{(r["errors"] or r["warnings"])[:110]}')
    if len(missing) == len(rows):
        print("no annotations submitted yet; run the viewer to start")
    print(f"report: {OUT / 'annotation_check.tsv'}")


if __name__ == "__main__":
    main()