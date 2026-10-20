#!/usr/bin/env python3
"""Re-audit the four-cancer HEST material itself.

Read-only. Records every defect found instead of filtering silently.
No labels, no scoring, no image copy.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np
import tifffile

ROOT = Path("/home/huyudi/013_spatial")
HEST = ROOT / "data/other_sources/hest1k"
OUT = ROOT / "infra/tls_four_cancer_labels_20261008/eligibility"
ORGANS = {"Bowel", "Lung", "Breast", "Liver"}
TENX = ROOT / "data/other_sources/10x_genomics/manifest_3yr.tsv"
PANEL_GOBP = ROOT / "infra/tls_multiscale_gobp_20261002/panel_gobp.json"
PANEL_CTRL = ROOT / "infra/tls_multiscale_gobp_20261002/panel_controls.json"
CACHE = ROOT / "data/other_sources/hest1k/pixel_size_vis"

ASSETS = {
    "wsi": ("wsis", "{sid}.tif", "essential"),
    "st": ("st", "{sid}.h5ad", "essential"),
    "tissue_contours": ("tissue_seg", "{sid}_contours.geojson", "essential"),
    "thumbnail": ("thumbnails", "{sid}_downscaled_fullres.jpeg", "auxiliary"),
    "tissue_vis": ("tissue_seg", "{sid}_vis.jpg", "auxiliary"),
    "spatial_plot": ("spatial_plots", "{sid}_spatial_plots.png", "auxiliary"),
    "cellvit": ("cellvit_seg", "{sid}_cellvit_seg.parquet", "auxiliary"),
}


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def load_panel() -> dict:
    panel = json.loads(PANEL_GOBP.read_text())
    panel |= json.loads(PANEL_CTRL.read_text())
    return panel


def read_h5ad_structure(path: Path) -> dict:
    out: dict = {}
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        obs_index = None
        for key in ("_index", "index"):
            if key in obs:
                obs_index = obs[key][:]
                break
        var_index = None
        for key in ("_index", "index"):
            if key in f["var"]:
                var_index = f["var"][key][:]
                break

        def decode(arr):
            return [
                v.decode() if isinstance(v, bytes) else str(v) for v in arr
            ] if arr is not None else None

        barcodes = decode(obs_index)
        genes = decode(var_index)
        out["n_obs"] = len(barcodes) if barcodes else 0
        out["n_var"] = len(genes) if genes else 0
        out["dup_barcodes"] = len(barcodes) - len(set(barcodes)) if barcodes else -1
        out["dup_genes"] = len(genes) - len(set(genes)) if genes else -1

        if "obsm" in f and "spatial" in f["obsm"]:
            sp = f["obsm"]["spatial"]
            out["spatial_rows"] = sp.shape[0]
            coords = np.asarray(sp)
            out["coord_nan"] = int(np.isnan(coords).sum())
            out["coord_finite"] = int(np.isfinite(coords).all(axis=1).sum())
            out["coord_min"] = [float(np.nanmin(coords[:, i])) for i in range(coords.shape[1])]
            out["coord_max"] = [float(np.nanmax(coords[:, i])) for i in range(coords.shape[1])]
        else:
            out["spatial_rows"] = None
            out["coord_nan"] = -1
            out["coord_finite"] = -1

        if "in_tissue" in obs:
            it = obs["in_tissue"][:]
            out["in_tissue_sum"] = int(np.asarray(it).sum())
        else:
            out["in_tissue_sum"] = None

        if "pxl_col_in_fullres" in obs and "pxl_row_in_fullres" in obs:
            cols = np.asarray(obs["pxl_col_in_fullres"][:], dtype=float)
            rows = np.asarray(obs["pxl_row_in_fullres"][:], dtype=float)
            out["pxl_rows"] = int(len(cols))
            out["pxl_finite"] = int(np.isfinite(cols).all() and np.isfinite(rows).all())
            out["pxl_min"] = [float(np.nanmin(cols)), float(np.nanmin(rows))]
            out["pxl_max"] = [float(np.nanmax(cols)), float(np.nanmax(rows))]
        else:
            out["pxl_rows"] = 0
            out["pxl_finite"] = -1
            out["pxl_min"] = None
            out["pxl_max"] = None

        xobj = f["X"] if "X" in f else None
        if isinstance(xobj, h5py.Group):
            out["x_kind"] = "sparse_group"
            out["x_dtype"] = str(xobj.attrs.get("h5sparse_format", ""))
            out["x_encoding"] = str(xobj.attrs.get("encoding-type", ""))
            out["x_shape"] = None
        elif xobj is not None:
            out["x_kind"] = "dense_dataset"
            out["x_dtype"] = str(xobj.dtype)
            out["x_encoding"] = ""
            out["x_shape"] = list(xobj.shape)
        else:
            out["x_kind"] = None
            out["x_dtype"] = None
            out["x_encoding"] = ""
            out["x_shape"] = None
        out["uns_spatial_keys"] = list(f["uns"]["spatial"].keys()) if "uns" in f and "spatial" in f["uns"] else None
        sf_path = None
        if "uns" in f and "spatial" in f["uns"]:
            for base in ("", "ST"):
                grp = f["uns"]["spatial"][base] if base else f["uns"]["spatial"]
                if isinstance(grp, h5py.Group) and "scalefactors" in grp:
                    sf_path = grp["scalefactors"]
                    out["scalefactors_under"] = base or "uns/spatial"
                    break
        if sf_path is not None:
            out["scalefactors"] = {}
            for k in sf_path.keys():
                arr = np.asarray(sf_path[k]).ravel()
                # tissue_hires_image_scalef is [scale, offset]; take the scale.
                out["scalefactors"][k] = float(arr[0]) if arr.size else None
            sd = out["scalefactors"].get("spot_diameter_fullres")
            scalef_key = None
            for cand in ("tissue_hires_image_scalef", "tissue_downscaled_fullres_scalef"):
                if isinstance(out["scalefactors"].get(cand), float):
                    scalef_key = cand
                    break
            # Visium spot is 55 um across; this is the authoritative fullres scale.
            out["image_scalef_key"] = scalef_key
            if isinstance(sd, float) and sd > 0:
                out["px_per_um_scalefactor"] = sd / 55.0
            else:
                out["px_per_um_scalefactor"] = None
        else:
            out["scalefactors"] = None
            out["scalefactors_under"] = None
            out["px_per_um_scalefactor"] = None
            out["px_per_um_estimated_field"] = None
    return out, (genes or [])


def read_tif_header(path: Path) -> dict:
    import warnings
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    out = {"tif_reader": "none"}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            with tifffile.TiffFile(path) as tf:
                page = tf.pages[0]
                series = tf.series[0]
                levels = 1
                try:
                    levels = len(series.levels)
                except (AttributeError, ValueError):
                    pass
                out.update({
                    "tif_shape": [int(page.imagewidth), int(page.imagelength)],
                    "tif_levels": levels,
                    "tif_pyramid_pages": len(tf.pages),
                    "tif_compression": str(page.compression),
                    "tif_bits": int(page.bitspersample),
                    "tif_reader": "tifffile",
                })
        except Exception as exc:
            out["tifffile_error"] = f"{type(exc).__name__}: {exc}"
    # tifffile fails on some Deflate pyramids; PIL is the independent fallback.
    try:
        with Image.open(path) as im:
            out["pil_size"] = list(im.size)
            out["pil_levels"] = getattr(im, "n_frames", 1)
            im.seek(min(1, getattr(im, "n_frames", 1) - 1))
            # Header/level navigation only. Decoding a 1.7 GB level is not needed
            # to prove the pyramid is readable.
            out["pil_level1_size"] = list(im.size)
            if out.get("tif_reader") == "none":
                out["tif_reader"] = "PIL_only"
    except Exception as exc:
        out["pil_error"] = f"{type(exc).__name__}: {exc}"
        out["tif_reader"] = "UNREADABLE"
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    panel = load_panel()
    with (HEST / "selected_samples.tsv").open() as handle:
        selected = list(csv.DictReader(handle, delimiter="\t"))
    with (HEST / "HEST_v1_3_0.csv").open() as handle:
        catalog = list(csv.DictReader(handle))

    in_scope = [
        r for r in selected
        if r["organ"] in ORGANS and "Visium" in r["st_technology"]
    ]
    in_scope.sort(key=lambda r: r["id"])
    sel_ids = {r["id"] for r in selected}

    findings: list[dict] = []

    def flag(section_id, check, severity, detail):
        findings.append({
            "section_id": section_id, "check": check,
            "severity": severity, "detail": detail,
        })

    # --- scope: catalog rows in scope of the plan but never downloaded ---
    scope_rows = []
    for row in catalog:
        if row["organ"] not in ORGANS or "Visium" not in row["st_technology"]:
            continue
        sid = row["id"]
        downloaded = sid in sel_ids
        scope_rows.append({
            "section_id": sid,
            "organ": row["organ"],
            "platform": row["st_technology"],
            "disease_state": row["disease_state"],
            "species": row["species"],
            "downloaded": "yes" if downloaded else "no",
            "excluded_because": "" if downloaded else (
                "non_human_species" if row["species"] != "Homo sapiens"
                else "disease_state_not_cancer_or_healthy"),
            "dataset_title": row["dataset_title"],
        })
    for row in scope_rows:
        if row["downloaded"] == "no":
            flag(row["section_id"], "catalog_not_downloaded", "scope",
                 f'{row["organ"]} {row["platform"]} {row["disease_state"]} '
                 f'{row["species"]} {row["excluded_because"]}')

    # --- per-section audit ---
    sections = []
    for row in in_scope:
        sid = row["id"]
        meta = json.loads((HEST / "metadata" / f"{sid}.json").read_text())
        if meta.get("species") != "Homo sapiens":
            flag(sid, "species_not_human", "blocking", str(meta.get("species")))

        record = {"section_id": sid}
        # selected_samples.tsv is the selection of record; metadata JSON must agree.
        for field, key in (("organ", "organ"), ("disease_state", "disease_state"),
                           ("st_technology", "st_technology")):
            mv = str(meta.get(key) or "").strip()
            if mv and mv != row[field]:
                flag(sid, "selection_metadata_disagree", "blocking",
                     f"{field}: selected={row[field]} metadata={mv}")
            record[f"meta_{field}"] = mv or row[field]
        pres = str(row.get("preservation_method") or "").strip()
        if not pres or pres.lower() in {"nan", "none", "na"}:
            flag(sid, "preservation_method_missing", "warning",
                 "neither selection nor metadata records how the tissue was preserved")
        record["preservation_method"] = pres
        missing_assets = []
        for name, (folder, pattern, kind) in ASSETS.items():
            path = HEST / folder / pattern.format(sid=sid)
            exists = path.exists()
            record[f"has_{name}"] = exists
            if not exists:
                missing_assets.append(f"{name}:{kind}")
        record["missing_assets"] = ";".join(missing_assets) if missing_assets else ""
        if missing_assets:
            essential = [a for a in missing_assets if a.endswith("essential")]
            flag(sid, "missing_essential_asset" if essential else "missing_auxiliary_asset",
                 "blocking" if essential else "warning", ";".join(missing_assets))

        width = as_float(meta.get("fullres_px_width"))
        height = as_float(meta.get("fullres_px_height"))
        est = as_float(meta.get("pixel_size_um_estimated"))
        try:
            record.update(read_tif_header(HEST / "wsis" / f"{sid}.tif"))
            if record.get("tif_reader") == "UNREADABLE":
                flag(sid, "wsi_undecodable_by_any_local_reader", "blocking",
                     f"{record.get('tifffile_error', '')} {record.get('pil_error', '')}".strip())
            elif record.get("pil_size") and width and height:
                record["wsi_dims_from_pil"] = str(record["pil_size"])
                if list(map(int, record["pil_size"])) != [int(width), int(height)]:
                    flag(sid, "tif_dims_disagree_with_metadata", "blocking",
                         f'wsi {record["pil_size"]} vs metadata {[width, height]}')
            if record.get("tifffile_error") and record.get("pil_size"):
                flag(sid, "wsi_needs_pil_fallback", "warning",
                     f'tifffile: {record["tifffile_error"]}; PIL decoded OK')
        except Exception as exc:
            flag(sid, "wsi_unreadable", "blocking", f"{type(exc).__name__}: {exc}")

        try:
            struct, genes = read_h5ad_structure(HEST / "st" / f"{sid}.h5ad")
            record.update(struct)
            if struct["spatial_rows"] is not None and struct["spatial_rows"] != struct["n_obs"]:
                flag(sid, "spatial_rows_mismatch", "blocking",
                     f'obsm/spatial {struct["spatial_rows"]} vs obs {struct["n_obs"]}')
            if struct["pxl_rows"] and struct["pxl_rows"] != struct["n_obs"]:
                flag(sid, "pxl_rows_mismatch", "blocking",
                     f'pxl rows {struct["pxl_rows"]} vs obs {struct["n_obs"]}')
            if struct["dup_barcodes"] > 0:
                flag(sid, "duplicate_barcodes", "blocking", str(struct["dup_barcodes"]))
            if struct["coord_nan"] > 0:
                flag(sid, "coord_nan", "blocking", str(struct["coord_nan"]))
            if struct["pxl_finite"] == 0 and struct["pxl_rows"]:
                flag(sid, "pxl_not_finite", "blocking", "non-finite fullres coords")
            nb = as_float(meta.get("nb_genes"))
            if nb and struct["n_var"] and abs(nb - struct["n_var"]) > 1:
                flag(sid, "nb_genes_disagree", "warning",
                     f'metadata {int(nb)} vs h5ad var {struct["n_var"]}')
            under = as_float(meta.get("spots_under_tissue"))
            if under and struct["in_tissue_sum"] and abs(under - struct["in_tissue_sum"]) > 1:
                flag(sid, "spots_under_tissue_disagree", "warning",
                     f'metadata {int(under)} vs obs in_tissue sum {struct["in_tissue_sum"]}')
            if not struct["uns_spatial_keys"]:
                flag(sid, "uns_spatial_missing", "warning", "no uns/spatial keys")
            if struct["scalefactors"] is None:
                flag(sid, "scalefactors_missing", "warning",
                     "no spot_diameter_fullres; pixel scale only from metadata estimate")

            # coordinate containment against declared image size
            if struct["pxl_min"] and width and height:
                over_x = max(0.0, -struct["pxl_min"][0], struct["pxl_max"][0] - width)
                over_y = max(0.0, -struct["pxl_min"][1], struct["pxl_max"][1] - height)
                record["coord_overflow_px"] = round(max(over_x, over_y), 1)
                record["coord_overflow_axis"] = (
                    "x" if over_x >= over_y else "y") if max(over_x, over_y) > 0 else ""
                record["coord_overflow_um"] = (
                    round(max(over_x, over_y) * est, 1) if est else "")
                if max(over_x, over_y) > 1:
                    # One Visium spot diameter (55 um) of slack is a boundary convention.
                    # Beyond that the spots sit in a different frame than the WSI.
                    tol_um = 55.0
                    severity = "warning" if record["coord_overflow_um"] <= tol_um else "blocking"
                    flag(sid, "coords_outside_declared_image", severity,
                         f'overflow {max(over_x, over_y):.0f}px '
                         f'({record["coord_overflow_um"]}um) on {record["coord_overflow_axis"]}, '
                         f'tolerance {tol_um:.0f}um')

            record["px_um_scalefactor"] = ""
            record["px_um_estimated_field"] = f'{est:.5f}' if est else ""
            record["px_um_authoritative"] = f'{est:.5f}' if est else ""
            sd_px = struct["scalefactors"].get("spot_diameter_fullres")
            record["spot_diameter_fullres_px"] = f"{sd_px:.2f}" if isinstance(sd_px, float) else ""
            if isinstance(sd_px, float) and est:
                # For Visium the spot is 55 um across, so 55/px is the fullres scale and
                # must agree with the metadata estimate. Visium HD bins are not 55 um,
                # so the implied constant differs and the 55 um conversion does not apply.
                implied_um = est * sd_px
                record["spot_um_constant_implied"] = f"{implied_um:.1f}"
                if row["st_technology"] == "Visium":
                    record["scale_constant_is_55um"] = "yes" if 52.0 <= implied_um <= 58.0 else "no"
                    record["px_um_from_scalefactor"] = f"{55.0 / sd_px:.5f}"
                    derived = 55.0 / sd_px
                    if abs(derived - est) / est > 0.05:
                        flag(sid, "pixel_scale_fields_disagree", "blocking",
                             f"55/spot_diameter gives {derived:.4f} um/px vs "
                             f"estimated {est:.4f}")
                else:
                    record["scale_constant_is_55um"] = "no"
                    flag(sid, "spot_diameter_not_55um_visium_hd", "warning",
                         f"HD bin: 55 um conversion invalid, implied constant "
                         f"{implied_um:.0f} um; use metadata estimate only")

            # metadata disease_state vs its own free-text comment
            comment = str(meta.get("disease_comment") or "")
            state = row["disease_state"]
            bad_terms = [t for t in ("hepatitis", "cirrhosis", "carcinoma", "cancer", "tumor",
                                     "adenoma", "fibrosis", "failure", "lesion", "metasta")
                         if t in comment.lower()]
            record["disease_comment"] = comment[:120]
            if state == "Healthy" and bad_terms:
                flag(sid, "healthy_label_contradicts_comment", "blocking",
                     f'disease_state=Healthy but comment says: {comment[:80]}')

            # frozen panel coverage, same rule as the scoring script
            index = set(genes)
            fractions = []
            for members in panel.values():
                present = sum(1 for g in members if g in index)
                fractions.append(present / len(members))
            frac = np.asarray(fractions)
            counts = np.asarray([
                sum(1 for g in members if g in index) for members in panel.values()
            ])
            valid = (frac >= 0.8) & (counts >= 10)
            record["panel_sets"] = len(panel)
            record["panel_valid"] = int(valid.sum())
            record["panel_valid_frac"] = round(float(valid.mean()), 4)
            if int(valid.sum()) < 100:
                flag(sid, "panel_coverage_low", "blocking",
                     f'only {int(valid.sum())}/{len(panel)} sets usable')
        except Exception as exc:
            flag(sid, "h5ad_unreadable", "blocking", f"{type(exc).__name__}: {exc}")

        sections.append(record)

    # --- unit structure: capture areas and replicates are not extra patients ---
    rep_groups = defaultdict(list)
    cap_groups = defaultdict(list)
    for row in in_scope:
        meta = json.loads((HEST / "metadata" / f"{row['id']}.json").read_text())
        sub = str(meta.get("subseries") or "")
        study = str(meta.get("study_link") or meta.get("dataset_title") or "unknown")
        import re
        for token in re.findall(r"(?:P\d+|A\d{5,}|CHC\d+|V\d+T\d+-\d+)", sub):
            rep_groups[(study, token)].append(row["id"])
        specimen = sub.split("_")[0]
        for token in re.findall(r"(?:^|_)([A-D]\d)(?=$|_)", sub):
            cap_groups[(study, specimen)].append((row["id"], token))
    for (study, token), sids in sorted(rep_groups.items()):
        if len(sids) > 1:
            flag(";".join(sids), "same_sample_token_multiple_sections", "unit",
                 f'{token} appears in {len(sids)} sections of one study')
    for key, members in sorted(cap_groups.items()):
        if len(members) > 1:
            flag(";".join(s for s, _ in members), "capture_areas_not_independent_sections", "unit",
                 f'{key[1]} appears as capture areas '
                 f'{sorted(t for _, t in members)} in one study')

    # --- 10x manifest integrity ---
    with TENX.open() as handle:
        tenx = list(csv.DictReader(handle, delimiter="\t"))
    url_dirs = defaultdict(list)
    title_dirs = defaultdict(list)
    for row in tenx:
        url_dirs[row["page_url"]].append(row["dataset_dir"])
        title_dirs[row["title"]].append(row["dataset_dir"])
    tenx_issues = {
        "rows": len(tenx),
        "unique_urls": len(url_dirs),
        "urls_with_multiple_dirs": sum(len(v) > 1 for v in url_dirs.values()),
        "duplicate_titles": sum(len(v) > 1 for v in title_dirs.values()),
        "non_ok_status": [r["dataset_dir"] for r in tenx if r["status"] != "ok"],
        "platforms": dict(Counter(r["platform"] for r in tenx)),
        "dirs_present_locally": sum(
            (ROOT / "data/other_sources/10x_genomics" / r["dataset_dir"]).exists()
            for r in tenx),
    }

    # --- write outputs ---
    with (OUT / "catalog_scope_2026.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(scope_rows[0]))
        writer.writeheader()
        writer.writerows(scope_rows)

    sec_fields = sorted({k for r in sections for k in r})
    with (OUT / "data_audit_sections.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sec_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in sections:
            writer.writerow({k: row.get(k, "") for k in sec_fields})

    with (OUT / "data_audit_findings.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["section_id", "check", "severity", "detail"], delimiter="\t")
        writer.writeheader()
        writer.writerows(findings)

    summary = {
        "n_sections_audited": len(sections),
        "catalog_rows_in_plan_scope": len(scope_rows),
        "catalog_rows_never_downloaded": sum(r["downloaded"] == "no" for r in scope_rows),
        "never_downloaded_non_human": sum(
            r["downloaded"] == "no" and r["excluded_because"] == "non_human_species"
            for r in scope_rows),
        "never_downloaded_human_other_state": sum(
            r["downloaded"] == "no" and r["excluded_because"] == "disease_state_not_cancer_or_healthy"
            for r in scope_rows),
        "findings_by_check": dict(sorted(Counter(f["check"] for f in findings).items())),
        "findings_by_severity": dict(sorted(Counter(f["severity"] for f in findings).items())),
        "tenx_manifest": tenx_issues,
        "labels_created": False,
    }
    (OUT / "data_audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()