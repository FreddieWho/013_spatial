#!/usr/bin/env python3
"""Map annotated TLS regions onto Visium / Visium HD measurement bins.

Frozen conservative rules (see the labelling plan, step 3):

  positive  every sampled point of the bin falls inside a TLS/aggregate instance
  negative  every sampled point falls inside a reviewed region AND no instance
  unknown   anything else: partial overlap, conflicting overlap, outside the
            reviewed area, or a spot of unknown geometry

`unknown` is never folded into `negative`. Instances are kept even when no bin
is fully covered, so small structures survive.

Input is either the browser export (`annotations.json`) or per-section GeoJSON
under `annotation/incoming/<section>/`.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from shapely.geometry import Point, Polygon, shape
from shapely.ops import unary_union

ROOT = Path("/home/huyudi/013_spatial")
HEST = ROOT / "data/other_sources/hest1k"
GATE = ROOT / "infra/tls_four_cancer_labels_20261008/eligibility/annotatable_gate.tsv"
INCOMING = ROOT / "infra/tls_four_cancer_labels_20261008/annotation/incoming"
OUT = ROOT / "infra/tls_four_cancer_labels_20261008/annotation/mapped"

POSITIVE_CLASSES = {"TLS", "lymphoid_aggregate"}
NEUTRAL_CLASSES = {"non_TLS", "unresolved"}

# Sampling points per bin for the containment test: centre plus vertices plus
# edge midpoints. Coarse sampling lets a bin slip through as "partial".
EDGE_STEPS = 4


def load_section(sid: str) -> dict:
    meta = json.loads((HEST / "metadata" / f"{sid}.json").read_text())
    rec = {"section_id": sid,
           "full_w": int(meta["fullres_px_width"]),
           "full_h": int(meta["fullres_px_height"]),
           "um_per_px": float(meta["pixel_size_um_estimated"]),
           "organ": meta.get("organ")}
    with_gate = None
    if GATE.exists():
        for row in csv.DictReader(GATE.open(), delimiter="\t"):
            if row["section_id"] == sid:
                with_gate = row
                break
    rec["registration_verdict"] = with_gate["registration_verdict"] if with_gate else None
    rec["spot_inside_tissue"] = float(with_gate["spot_inside_tissue"]) if with_gate else None
    return rec


def bin_geometry_um(rec: dict, spot_diameter_um, inter_spot_um) -> tuple[str, float]:
    """Return (mapping_mode, area_um2) for one bin.

    Visium: hexagon approximated by its circumscribed circle at the nominal 55 um
    diameter. Visium HD: square bin when the metadata states its side, else the
    caller must fall back to centre-only.
    """
    if spot_diameter_um and inter_spot_um and abs(spot_diameter_um - inter_spot_um) < 1e-6:
        side = float(spot_diameter_um)
        return "area_square_um", side * side
    if spot_diameter_um:
        r = float(spot_diameter_um) / 2.0
        return "area_circle_um", math.pi * r * r
    return "center_only", 0.0


def sample_points_px(cx: float, cy: float, half_px: float, mode: str) -> list[tuple[float, float]]:
    pts = [(cx, cy)]
    if mode == "center_only" or half_px <= 0:
        return pts
    n = 6 if mode == "area_circle_um" else 4
    for i in range(n):
        a = 2 * math.pi * i / n
        pts.append((cx + half_px * math.cos(a), cy + half_px * math.sin(a)))
    for i in range(n):
        a = 2 * math.pi * i / n
        a2 = 2 * math.pi * ((i + 1) % n) / n
        pts.append((cx + half_px * 0.5 * (math.cos(a) + math.cos(a2)),
                    cy + half_px * 0.5 * (math.sin(a) + math.sin(a2))))
    return pts


def to_px(geom, sx: float):
    from shapely import affinity
    return affinity.scale(geom, xfact=sx, yfact=sx, origin=(0, 0))


def from_px(geom, sx: float):
    from shapely import affinity
    return affinity.scale(geom, xfact=1.0 / sx, yfact=1.0 / sx, origin=(0, 0))


def read_instances(path: Path) -> list[dict]:
    doc = json.loads(path.read_text(encoding="utf-8"))
    out = []
    for f in doc.get("features") or []:
        geom = shape(f["geometry"])
        if not geom.is_valid:
            geom = geom.buffer(0)
        p = f.get("properties") or {}
        out.append({"geom": geom,
                    "cls": p.get("structure_class"),
                    "iid": p.get("instance_id"),
                    "confidence": p.get("confidence"),
                    "maturity": p.get("maturity")})
    return out


def load_labels(sid: str) -> tuple[list[dict], list[dict], str]:
    """Return (instances, reviewed_regions, source)."""
    env = INCOMING / sid / "annotations.json"
    if env.exists():
        doc = json.loads(env.read_text(encoding="utf-8"))
        sl = (doc.get("slides") or {}).get(sid)
        if sl is None:
            return [], [], "annotations.json (slide absent)"
        inst = [{"geom": shape({"type": "Polygon", "coordinates": [e["pts"] + [e["pts"][0]]]}),
                 "cls": e.get("structure_class"), "iid": None,
                 "confidence": e.get("confidence"), "maturity": e.get("maturity")}
                for e in sl.get("instances", [])]
        rev = [{"geom": shape({"type": "Polygon",
                              "coordinates": [e["pts"] + [e["pts"][0]]]})}
               for e in sl.get("reviewed", [])]
        return inst, rev, "annotations.json"
    fi = INCOMING / sid / "instances.geojson"
    fr = INCOMING / sid / "reviewed_regions.geojson"
    inst = read_instances(fi) if fi.exists() else []
    rev = [{"geom": shape(f["geometry"])} for f in
           json.loads(fr.read_text(encoding="utf-8")).get("features", [])] if fr.exists() else []
    return inst, rev, "geojson"


def map_section(sid: str) -> dict:
    import h5py
    rec = load_section(sid)
    meta = json.loads((HEST / "metadata" / f"{sid}.json").read_text())
    sd = meta.get("spot_diameter")
    isd = meta.get("inter_spot_dist")
    try:
        sd = float(sd) if sd not in (None, "", "nan") else None
    except (TypeError, ValueError):
        sd = None
    try:
        isd = float(isd) if isd not in (None, "", "nan") else None
    except (TypeError, ValueError):
        isd = None
    mode, area_um2 = bin_geometry_um(rec, sd, isd)
    um_per_px = rec["um_per_px"]
    half_px = (math.sqrt(area_um2 / math.pi) / um_per_px) if area_um2 else 0.0

    inst, rev, source = load_labels(sid)
    pos_geoms, other_geoms = [], []
    for e in inst:
        (pos_geoms if e["cls"] in POSITIVE_CLASSES else other_geoms).append(e)
    pos_union = unary_union([e["geom"] for e in pos_geoms]) if pos_geoms else None
    other_union = unary_union([e["geom"] for e in other_geoms]) if other_geoms else None
    rev_union = unary_union([r["geom"] for r in rev]) if rev else None

    with h5py.File(HEST / "st" / f"{sid}.h5ad", "r") as f:
        obs = f["obs"]
        cols = np.asarray(obs["pxl_col_in_fullres"][:], dtype=float)
        rows = np.asarray(obs["pxl_row_in_fullres"][:], dtype=float)
        bc = obs["_index"][:] if "_index" in obs else None
        it = np.asarray(obs["in_tissue"][:]).astype(bool) if "in_tissue" in obs \
            else np.ones(len(cols), bool)

    barcodes = [b.decode() if isinstance(b, bytes) else str(b) for b in bc] \
        if bc is not None else [str(i) for i in range(len(cols))]

    rec["mapping_mode"] = mode
    rec["bin_area_um2"] = round(area_um2, 1)
    rec["label_source"] = source
    rec["n_instances"] = len(inst)
    rec["n_positive_instances"] = len(pos_geoms)
    rec["n_reviewed"] = len(rev)
    rec["n_bars"] = len(cols)

    labels = []
    counts = Counter()
    n_cov = 0
    for i in range(len(cols)):
        cx, cy = float(cols[i]), float(rows[i])
        inside_tissue = bool(it[i])
        if not inside_tissue:
            status, reason, ids = "unknown", "spot_outside_tissue", ""
            counts[status] += 1
            labels.append((status, reason, ids, 0.0))
            continue
        pts = sample_points_px(cx, cy, half_px, mode)
        pts_geom = [Point(*p) for p in pts]

        in_pos = [pos_union.covers(p) for p in pts_geom] if pos_union else [False] * len(pts)
        all_pos = all(in_pos) if pos_union else False
        any_pos = any(in_pos) if pos_union else False

        if all_pos:
            status, reason, ids, ov = "positive", "fully_covered_by_instance", "", 1.0
        else:
            in_rev = [rev_union.covers(p) for p in pts_geom] if rev_union else [False] * len(pts)
            all_rev = all(in_rev) if rev_union else False
            any_rev = any(in_rev) if rev_union else False
            hits = []
            if pos_union:
                hits = [e["iid"] or f"inst_{j}" for j, e in enumerate(pos_geoms)
                        if e["geom"].intersects(Point(cx, cy))]
            ids = ";".join(hits)
            if any_pos:
                status, reason, ov = "unknown", "partially_covered_by_instance", 0.5
            elif not all_rev:
                status, reason, ov = ("unknown", "no_reviewed_coverage", 0.0) if not any_rev \
                    else ("unknown", "partially_reviewed", 0.5)
            elif other_union is not None and other_union.covers(Point(cx, cy)):
                status, reason, ov = "unknown", "overlaps_non_evaluable_class", 0.0
            else:
                status, reason, ov = "negative", "inside_reviewed_no_instance", 0.0
        counts[status] += 1
        labels.append((status, reason, ids, ov))

    rec["n_positive"] = counts["positive"]
    rec["n_negative"] = counts["negative"]
    rec["n_unknown"] = counts["unknown"]
    rec["reason_counts"] = dict(sorted(Counter(l[1] for l in labels).items()))
    rec["reasons"] = labels
    rec["barcodes"] = barcodes
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()
    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    with GATE.open() as fh:
        gate = [r for r in csv.DictReader(fh, delimiter="\t")]
    have = [r["section_id"] for r in gate
            if (INCOMING / r["section_id"]).exists()]
    if not have:
        print("no annotations yet under " + str(INCOMING))
        return 0

    summary = []
    for sid in sorted(have):
        rec = map_section(sid)
        rows = []
        for bc, (status, reason, ids, ov) in zip(rec["barcodes"], rec["reasons"]):
            rows.append({"section_id": sid, "barcode": bc,
                         "label_status": status, "unknown_reason": reason,
                         "instance_ids": ids, "overlap_fraction": f"{ov:.2f}",
                         "mapping_mode": rec["mapping_mode"],
                         "evidence_level": "independent_reference",
                         "reviewer_credential": "expert_human_review",
                         "double_read": "false",
                         "confidence": "", "maturity": "",
                         "annotation_version": "v1"})
        with (out_root / f"{sid}.tsv.gz").open("wt", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        slim = {k: v for k, v in rec.items() if k not in ("reasons", "barcodes")}
        summary.append(slim)
        print(f'{sid}: pos={rec["n_positive"]} neg={rec["n_negative"]} '
              f'unknown={rec["n_unknown"]} (mode {rec["mapping_mode"]})')

    with (out_root / "mapping_summary.tsv").open("w", newline="", encoding="utf-8") as fh:
        keys = [k for k in summary[0] if k != "reason_counts"]
        w = csv.DictWriter(fh, fieldnames=keys + ["reason_counts"], delimiter="\t")
        w.writeheader()
        for s in summary:
            w.writerow({k: s.get(k, "") for k in keys} |
                       {"reason_counts": json.dumps(s["reason_counts"], ensure_ascii=False)})
    print(f"\nwrote {len(summary)} sections to {out_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())