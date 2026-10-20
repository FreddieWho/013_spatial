#!/usr/bin/env python3
"""Four-cancer HEST identity and image-eligibility inventory.

Reads local metadata only. Does not score genes, copy WSIs, or create labels.
Patient strings are study-local. Blank identity stays unknown.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

import h5py

ROOT = Path("/home/huyudi/013_spatial")
HEST = ROOT / "data/other_sources/hest1k"
TENX = ROOT / "data/other_sources/10x_genomics/manifest_3yr.tsv"
SCORES = ROOT / "infra/tls_fourcancer_20261006/manifest.json"
OUT = ROOT / "infra/tls_four_cancer_labels_20261008/eligibility"
ORGANS = {"Bowel", "Lung", "Breast", "Liver"}
BLANK = {"", "nan", "none", "na", "n/a", "unknown", "null"}


def clean(value) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in BLANK else text


def as_float(value):
    text = clean(value)
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def slug(text: str) -> str:
    base = re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")
    return base[:80] or "untitled"


def url_key(url: str) -> str:
    text = clean(url).rstrip("/").lower()
    if not text or text in BLANK:
        return ""
    return text.split("?")[0].rstrip("/").split("/")[-1]


def sample_tokens(subseries: str) -> list[str]:
    """Sample-like tokens in a subseries label. Not a patient ID."""
    # Underscore is a word character, so do not require a trailing boundary.
    # A1/B1/C1/D1 are capture-area suffixes, not specimen IDs.
    return re.findall(r"(?:P\d+|A\d{5,}|CHC\d+|V\d+T\d+-\d+)", subseries)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def spatial_bounds(h5ad: Path, width: float, height: float) -> dict:
    with h5py.File(h5ad, "r") as handle:
        has_spatial = "obsm" in handle and "spatial" in handle["obsm"]
        obs = handle["obs"]
        if "pxl_col_in_fullres" not in obs or "pxl_row_in_fullres" not in obs:
            return {
                "spatial_coords": "missing",
                "n_coord_spots": 0,
                "frac_inside_image": "",
                "coord_bounds": "not_checked",
            }
        cols = obs["pxl_col_in_fullres"][:]
        rows = obs["pxl_row_in_fullres"][:]
    n = int(len(cols))
    if n == 0 or width is None or height is None:
        return {
            "spatial_coords": "present" if has_spatial else "missing",
            "n_coord_spots": n,
            "frac_inside_image": "",
            "coord_bounds": "not_checked",
        }
    inside = (
        (cols >= 0) & (rows >= 0) & (cols < width) & (rows < height)
    ).sum()
    frac = float(inside) / n
    return {
        "spatial_coords": "present" if has_spatial else "array_only",
        "n_coord_spots": n,
        "frac_inside_image": f"{frac:.6f}",
        "coord_bounds": "inside" if frac >= 0.99 else "out_of_image",
    }


def load_score_flags() -> dict[str, str]:
    if not SCORES.exists():
        return {}
    flags = {}
    for row in json.loads(SCORES.read_text()):
        sid = row["section_id"].removeprefix("hest-")
        usable = row.get("panel_usable")
        flags[sid] = "panel_usable" if usable else "scored_panel_unusable"
    return flags


def load_10x() -> list[dict]:
    with TENX.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with (HEST / "selected_samples.tsv").open() as handle:
        selected = list(csv.DictReader(handle, delimiter="\t"))
    with (HEST / "HEST_v1_3_0.csv").open() as handle:
        catalog = {row["id"]: row for row in csv.DictReader(handle)}
    scores = load_score_flags()
    tenx = load_10x()
    tenx_by_url = defaultdict(list)
    tenx_by_title = defaultdict(list)
    for row in tenx:
        key = url_key(row.get("page_url", ""))
        if key:
            tenx_by_url[key].append(row["dataset_dir"])
        title = clean(row.get("title", "")).lower()
        if title:
            tenx_by_title[title].append(row["dataset_dir"])

    rows = [
        row for row in selected
        if row["organ"] in ORGANS and "Visium" in row["st_technology"]
    ]
    rows.sort(key=lambda row: row["id"])
    sections = []
    overlaps = []
    for row in rows:
        sid = row["id"]
        meta_path = HEST / "metadata" / f"{sid}.json"
        meta = json.loads(meta_path.read_text())
        catalog_row = catalog.get(sid, {})
        patient = clean(meta.get("patient"))
        subseries = clean(meta.get("subseries"))
        tokens = sample_tokens(subseries)
        title = clean(meta.get("dataset_title")) or clean(catalog_row.get("dataset_title"))
        study = clean(meta.get("study_link"))
        source_label = study or title or "unknown_source"
        source_id = slug(source_label)
        # A bare "Patient N" is not identity when the same study uses distinct sample labels.
        patient_status = "known_within_study" if patient else "unknown"
        patient_key = f"{source_id}::{patient}" if patient else "unknown"
        identity_note = "metadata_patient_string" if patient else "blank_patient"
        emb = as_float(meta.get("pixel_size_um_embedded"))
        est = as_float(meta.get("pixel_size_um_estimated"))
        width = as_float(meta.get("fullres_px_width"))
        height = as_float(meta.get("fullres_px_height"))
        wsi = HEST / "wsis" / f"{sid}.tif"
        h5ad = HEST / "st" / f"{sid}.h5ad"
        bounds = spatial_bounds(h5ad, width or 0, height or 0) if h5ad.exists() else {
            "spatial_coords": "missing_file",
            "n_coord_spots": 0,
            "frac_inside_image": "",
            "coord_bounds": "not_checked",
        }
        embedded_unusable = emb is None or emb == 0 or abs(emb - 1.0) < 1e-6
        min_side = min(width, height) if width and height else None
        if not wsi.exists() or not h5ad.exists():
            image_grade = "files_missing"
        elif bounds["coord_bounds"] == "out_of_image":
            image_grade = "registration_uncertain"
        elif est is None or not (0.05 < est <= 0.6) or min_side is None or min_side < 8000:
            image_grade = "coarse_or_uncertain"
        else:
            image_grade = "fullres_candidate"
        disease = clean(row["disease_state"]) or clean(meta.get("disease_state"))
        if disease == "Cancer" and image_grade == "fullres_candidate":
            pilot_pool = "cancer_fullres"
        elif disease == "Cancer":
            pilot_pool = "cancer_not_fullres"
        elif disease == "Healthy":
            pilot_pool = "healthy_background_not_negative"
        elif disease == "Tumor":
            pilot_pool = "tumor_semantics_unchecked"
        else:
            pilot_pool = "other"
        tech_meta = clean(meta.get("st_technology"))
        download = clean(meta.get("download_page_link1"))
        matched = []
        ukey = url_key(download)
        if ukey and ukey in tenx_by_url:
            matched.extend(tenx_by_url[ukey])
        title_key = title.lower()
        if title_key in tenx_by_title:
            matched.extend(tenx_by_title[title_key])
        matched = sorted(set(matched))
        if matched:
            overlaps.append({
                "section_id": sid,
                "hest_download_page": download,
                "hest_title": title,
                "tenx_dataset_dirs": ";".join(matched),
                "match": "possible_same_release",
                "use": "do_not_add_to_hest_count",
            })
        sections.append({
            "section_id": sid,
            "source_id": source_id,
            "source_label": source_label,
            "dataset_title": title,
            "subseries": subseries,
            "sample_tokens": ";".join(tokens),
            "organ": row["organ"],
            "tissue": clean(meta.get("tissue")),
            "disease_state": disease,
            "oncotree_code": clean(row.get("oncotree_code")) or clean(meta.get("oncotree_code")),
            "patient_raw": patient,
            "patient_status": patient_status,
            "patient_key": patient_key,
            "patient_scope": "study_local_not_global",
            "identity_note": identity_note,
            "platform": row["st_technology"],
            "platform_metadata": tech_meta,
            "platform_conflict": "yes" if tech_meta and tech_meta != row["st_technology"] else "no",
            "preservation_method": clean(row.get("preservation_method")) or clean(meta.get("preservation_method")),
            "study_link": study,
            "download_page": download,
            "wsi_path": str(wsi.relative_to(ROOT)),
            "wsi_bytes": wsi.stat().st_size if wsi.exists() else 0,
            "h5ad_path": str(h5ad.relative_to(ROOT)),
            "metadata_sha256": sha256_file(meta_path),
            "wsi_sha256": "deferred_not_read",
            "pixel_size_um_embedded": "" if emb is None else f"{emb:.6g}",
            "pixel_size_um_estimated": "" if est is None else f"{est:.6g}",
            "embedded_pixel_size_unusable": "yes" if embedded_unusable else "no",
            "fullres_px_width": "" if width is None else str(int(width)),
            "fullres_px_height": "" if height is None else str(int(height)),
            "image_grade": image_grade,
            "visual_registration_preview": "NOT_RUN",
            "annotatable": "not_yet",
            "pilot_pool": pilot_pool,
            "prior_molecular_score": scores.get(sid, "absent"),
            "score_used_for_selection": "no",
            "tenx_possible_same_release": "yes" if matched else "no",
            **bounds,
        })

    by_source = defaultdict(list)
    for row in sections:
        by_source[row["source_id"]].append(row)
    for members in by_source.values():
        known = [m for m in members if m["patient_status"] == "known_within_study"]
        by_patient = defaultdict(list)
        for member in known:
            by_patient[member["patient_raw"]].append(member)
        for patient_raw, group in by_patient.items():
            token_sets = {tuple(t for t in m["sample_tokens"].split(";") if t) for m in group}
            token_sets.discard(tuple())
            if len(token_sets) > 1:
                for member in group:
                    member["patient_status"] = "unknown"
                    member["patient_key"] = "unknown"
                    member["identity_note"] = "patient_string_conflicts_with_distinct_sample_labels"
            elif len(token_sets) == 1 and len(group) > 1:
                for member in group:
                    member["identity_note"] = "patient_string_agrees_with_shared_sample_token"

    fields = list(sections[0])
    with (OUT / "sections.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sections)

    groups = defaultdict(list)
    for row in sections:
        groups[row["source_id"]].append(row)
    group_rows = []
    for source_id, members in sorted(groups.items(), key=lambda item: item[0]):
        known = {m["patient_key"] for m in members if m["patient_status"] == "known_within_study"}
        group_rows.append({
            "source_id": source_id,
            "source_label": members[0]["source_label"],
            "n_sections": len(members),
            "organs": ",".join(sorted({m["organ"] for m in members})),
            "platforms": ",".join(sorted({m["platform"] for m in members})),
            "disease_states": ",".join(sorted({m["disease_state"] for m in members})),
            "n_known_patients_within_study": len(known),
            "n_unknown_identity_sections": sum(m["patient_status"] == "unknown" for m in members),
            "known_patient_keys": ";".join(sorted(known)),
            "external_test_candidate": "no",
        })
    with (OUT / "source_groups.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(group_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(group_rows)

    token_hits = defaultdict(list)
    for row in sections:
        for token in [t for t in row["sample_tokens"].split(";") if t]:
            token_hits[(row["source_id"], token)].append(row["section_id"])
    shared = []
    for (source_id, token), sids in sorted(token_hits.items()):
        if len(sids) < 2:
            continue
        members = [r for r in sections if r["section_id"] in sids]
        shared.append({
            "source_id": source_id,
            "sample_token": token,
            "n_sections": len(sids),
            "section_ids": ";".join(sids),
            "platforms": ",".join(sorted({m["platform"] for m in members})),
            "patient_status_after_rule": ",".join(sorted({m["patient_status"] for m in members})),
            "interpretation": "possible_same_material_not_confirmed_patient",
        })
    with (OUT / "shared_sample_tokens.tsv").open("w", newline="") as handle:
        if shared:
            writer = csv.DictWriter(handle, fieldnames=list(shared[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(shared)
        else:
            handle.write("source_id\tsample_token\n")

    with (OUT / "overlap_10x.tsv").open("w", newline="") as handle:
        if overlaps:
            writer = csv.DictWriter(handle, fieldnames=list(overlaps[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(overlaps)
        else:
            handle.write("section_id\tmatch\n")

    summary = {
        "n_sections": len(sections),
        "n_sources": len(group_rows),
        "patient_known": sum(r["patient_status"] == "known_within_study" for r in sections),
        "patient_unknown": sum(r["patient_status"] == "unknown" for r in sections),
        "image_grade": dict(sorted(
            (grade, sum(r["image_grade"] == grade for r in sections))
            for grade in sorted({r["image_grade"] for r in sections})
        )),
        "tenx_overlap_sections": len(overlaps),
        "tenx_manifest_rows": len(tenx),
        "tenx_unique_urls": len({url_key(r.get("page_url", "")) for r in tenx if url_key(r.get("page_url", ""))}),
        "visual_registration_preview": "NOT_RUN",
        "pilot_locked": False,
        "labels_created": False,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
