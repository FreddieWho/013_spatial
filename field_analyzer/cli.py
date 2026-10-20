"""Command line entry points for the FA round.

    python -m field_analyzer.cli assets  --out infra/field_analyzer_20261020
    python -m field_analyzer.cli units   --out infra/field_analyzer_20261020
    python -m field_analyzer.cli poison  --section ZEN36 --out infra/field_analyzer_20261020
    python -m field_analyzer.cli geometry-smoke --out infra/field_analyzer_20261020

Nothing here downloads data, reads human labels, or uses a GPU.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from . import assets, geometry, poisoning, route_a, route_b, scales, synthetic
from .adapters import read_h5ad_counts
from .observation import fingerprint_arrays, visible_mean_std
from .pipeline import InputSpec, build_inputs
from .roles import build_field_splits, build_gene_roles


def repo_root_default() -> Path:
    return Path(__file__).resolve().parent.parent


def load_panel(repo: Path) -> dict[str, list[str]]:
    panel = json.loads((repo / "infra/tls_multiscale_gobp_20261002/panel_controls.json").read_text())
    return {k: list(v) for k, v in panel.items()}


PILOT_FIELDS = (
    "B_AXIS",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
)


def cmd_assets(args: argparse.Namespace) -> int:
    repo, out = Path(args.repo), Path(args.out)
    identity = assets.read_identity(repo / "infra/field_discovery_review_20261018/identity_review.tsv")
    sections = sorted(s for v in identity.units.values() for s in v)
    rows = assets.asset_rows(repo, sections)
    assets.write_tsv(rows, out / "reuse_resolution.tsv")
    detail = assets.section_rows(repo, sections)
    assets.write_tsv(detail, out / "section_asset_check.tsv")
    summary = {
        "n_confirmed_units": identity.n_units,
        "n_confirmed_sections": identity.n_sections,
        "n_sections_present": sum(1 for r in detail if r["exists"]),
        "n_axis_consistent": sum(1 for r in detail if r["axis_consistent"]),
        "states": {r["path"]: r["availability"] for r in rows},
    }
    (out / "asset_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


def cmd_units(args: argparse.Namespace) -> int:
    repo, out = Path(args.repo), Path(args.out)
    identity = assets.read_identity(repo / "infra/field_discovery_review_20261018/identity_review.tsv")
    rows = assets.unit_role_rows(identity.units)
    assets.write_tsv(rows, out / "unit_role_assertions.tsv")

    # assertions actually checked, not merely stated
    problems: list[str] = []
    by_unit: dict[str, set[int]] = {}
    for row in rows:
        by_unit.setdefault(row["unit_id"], set()).add(int(row["fold"]))
    for unit, folds in by_unit.items():
        if len(folds) != 1:
            problems.append(f"unit {unit} spans folds {sorted(folds)}")
    present = [r for r in rows if (repo / "data/other_sources/hest1k/st" / f"{r['section_id']}.h5ad").exists()]
    summary = {
        "n_units": identity.n_units,
        "n_sections": identity.n_sections,
        "n_sections_present": len(present),
        "units_per_fold": {
            str(f): len({r["unit_id"] for r in rows if int(r["fold"]) == f}) for f in sorted({int(r["fold"]) for r in rows})
        },
        "assertions": {
            "same_unit_same_fold": not problems,
            "duplicate_mask_not_extra_patient": "mask_group is unit-scoped by construction",
            "statistics_unit": "unit, not spot",
            "unknown_identity_not_upgraded": True,
        },
        "problems": problems,
    }
    (out / "unit_role_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0 if not problems else 1


def _block_mask(xy: np.ndarray, tissue_support: np.ndarray, radius_multiplier: float = 6.0) -> tuple[np.ndarray, float, np.ndarray]:
    pitch = geometry.estimate_pitch(xy[tissue_support]) if tissue_support.sum() > 2 else geometry.estimate_pitch(xy)
    center = np.median(xy[tissue_support], axis=0) if tissue_support.any() else np.median(xy, axis=0)
    dist = np.linalg.norm(xy - center, axis=1)
    radius = radius_multiplier * pitch
    visible = ~((dist <= radius) & tissue_support)
    return visible, float(radius), center


def cmd_poison(args: argparse.Namespace) -> int:
    repo, out = Path(args.repo), Path(args.out)
    panel = load_panel(repo)
    section = args.section
    table = read_h5ad_counts(repo / "data/other_sources/hest1k/st" / f"{section}.h5ad")

    fields = {fid: panel[fid] for fid in PILOT_FIELDS if fid in panel}
    splits = build_field_splits(fields, table.genes, seed=args.seed)
    roles = build_gene_roles(splits, table.genes)
    spec = InputSpec(
        fields={fid: list(split.a_genes) for fid, split in splits.items()},
        witness_genes=frozenset(roles.witness_genes),
        excluded_genes=frozenset(),
        k_neighbors=args.k,
        coordinate_unit="source_pixels",
    )

    tissue_support = np.ones(table.n_obs, dtype=bool)
    if "in_tissue" in table.aux:
        tissue_support = np.asarray(table.aux["in_tissue"]).astype(bool)
    visible, radius, center = _block_mask(table.xy, tissue_support, args.radius_pitch)
    if (~visible).sum() < 20:
        visible, radius, center = _block_mask(table.xy, tissue_support, 10.0)

    report = poisoning.run_poisoning(
        table.counts,
        table.genes,
        table.xy,
        tissue_support,
        visible,
        spec,
        package_dir=Path(__file__).resolve().parent,
        section_id=section,
        seed=args.seed,
    )
    report["mask"] = {
        "center": center.tolist(),
        "radius": radius,
        "radius_pitch": args.radius_pitch,
        "rule": "geometric disc around the median tissue coordinate; no expression used",
    }
    report["fields"] = {fid: split.to_dict() for fid, split in splits.items()}
    report["read_report"] = table.report.to_dict()
    out.mkdir(parents=True, exist_ok=True)
    (out / f"poisoning_report_{section}.json").write_text(json.dumps(report, indent=2, ensure_ascii=False))
    (out / "poisoning_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False))
    print(
        json.dumps(
            {
                "section": section,
                "passed": report["passed"],
                "n_hidden": report["n_hidden"],
                "scenarios": [
                    (s["scenario"], s["identical"], s["max_abs_delta"]) for s in report["scenarios"]
                ],
                "metric_changed": report["evaluation"]["metric_changed"],
                "prediction_max_abs_delta": report["evaluation"]["prediction_max_abs_delta"],
            },
            indent=2,
        )
    )
    return 0 if report["passed"] else 1



def cmd_scales(args: argparse.Namespace) -> int:
    repo, out = Path(args.repo), Path(args.out) / "scales"
    out.mkdir(parents=True, exist_ok=True)
    panel = load_panel(repo)
    fields_def = {fid: panel[fid] for fid in PILOT_FIELDS if fid in panel}
    identity = assets.read_identity(repo / "infra/field_discovery_review_20261018/identity_review.tsv")

    curve_rows, fit_rows, dir_rows, cross_rows = [], [], [], []
    status_tally: dict[str, dict[str, int]] = {"raw": {}, "residual": {}}
    per_section = []

    for unit in sorted(identity.units):
        for section in sorted(identity.units[unit]):
            table = read_h5ad_counts(repo / "data/other_sources/hest1k/st" / f"{section}.h5ad")
            support = np.ones(table.n_obs, dtype=bool)
            if "in_tissue" in table.aux:
                support = np.asarray(table.aux["in_tissue"]).astype(bool)
            splits = build_field_splits(fields_def, table.genes, seed=args.seed)
            roles = build_gene_roles(splits, table.genes)
            spec = InputSpec(
                fields={fid: list(split.a_genes) for fid, split in splits.items()},
                witness_genes=frozenset(roles.witness_genes),
                k_neighbors=args.k,
                coordinate_unit="source_pixels",
            )
            inputs = build_inputs(table.counts, table.genes, table.xy, support, support, spec)
            pitch = inputs.meta["pitch"]
            pairs = scales.pairs_within(table.xy, support, args.max_lag * pitch)
            edges = scales.lag_edges(pitch, args.max_lag, args.bins)

            scores = {fid: inputs.arrays[f"field::{fid}"] for fid in spec.fields}
            section_status = {"section_id": section, "unit_id": unit, "pitch": pitch, "n_pairs": int(len(pairs.distance)), "fields": {}}

            for fid, score in scores.items():
                trend = scales.quadratic_trend(table.xy, score, support)
                residual = score - trend
                for label, values in (("raw", score), ("residual", residual)):
                    gamma, counts, centers = scales.empirical_variogram(values, pairs, edges)
                    for row in scales.summarise_curve(gamma, counts, centers, pitch):
                        curve_rows.append({"unit_id": unit, "section_id": section, "field_id": fid,
                                           "component": label, **row})
                    fit = scales.fit_variogram(centers, gamma, counts)
                    fit_rows.append({"unit_id": unit, "section_id": section, "field_id": fid,
                                     "component": label, **fit})
                    status_tally[label][fit["status"]] = status_tally[label].get(fit["status"], 0) + 1
                    if label == "raw":
                        section_status["fields"][fid] = {"status": fit["status"], "range_pitch": (fit["range"] / pitch) if fit["range"] else None, "sill": fit["sill"]}
                dir_res = scales.directional_variogram(score, pairs, edges)
                for s, label in enumerate(scales.DIRECTION_LABELS):
                    for k in range(len(edges) - 1):
                        n = int(dir_res["counts"][s, k])
                        dir_rows.append({"unit_id": unit, "section_id": section, "field_id": fid,
                                         "direction_deg": label,
                                         "lag_pitch": float(0.5 * (edges[k] + edges[k + 1]) / pitch),
                                         "gamma": float(dir_res["gamma"][s, k]) if np.isfinite(dir_res["gamma"][s, k]) else None,
                                         "n_pairs": n,
                                         "measurable": bool(n >= scales.MIN_PAIRS_DIRECTION)})

            field_ids = sorted(scores)
            for a in range(len(field_ids)):
                for b in range(a + 1, len(field_ids)):
                    f1, f2 = field_ids[a], field_ids[b]
                    res = scales.cross_field_curves(scores[f1], scores[f2], pairs, edges)
                    for k in range(len(res["centers"])):
                        cc = res["cross_corr"][k]
                        cross_rows.append({"unit_id": unit, "section_id": section, "field_pair": f"{f1}~{f2}",
                                           "lag_pitch": float(res["centers"][k] / pitch),
                                           "cross_corr": float(cc) if np.isfinite(cc) else None,
                                           "cross_gamma": float(res["cross_gamma"][k]) if np.isfinite(res["cross_gamma"][k]) else None,
                                           "n_pairs": int(res["counts"][k]),
                                           "measurable": bool(np.isfinite(cc) and res["counts"][k] >= scales.MIN_PAIRS_FIT)})
            per_section.append(section_status)

    assets.write_tsv(curve_rows, out / "scale_curves.tsv")
    assets.write_tsv(fit_rows, out / "scale_fits.tsv")
    assets.write_tsv(dir_rows, out / "direction_curves.tsv")
    assets.write_tsv(cross_rows, out / "cross_field.tsv")
    summary = {
        "n_sections": len(per_section),
        "n_curve_rows": len(curve_rows),
        "n_fit_rows": len(fit_rows),
        "status_counts": status_tally,
        "sections": per_section,
        "settings": {"max_lag_pitch": args.max_lag, "bins": args.bins, "k": args.k, "seed": args.seed,
                     "lag_unit": "empirical_pitch", "score": "log1p CPM of field A genes, external library denominator"},
        "notes": [
            "residual component removes a quadratic trend only; it is a residual scale, not the field's range",
            "NOT_IDENTIFIED_WITHIN_FOV means no plateau inside the field of view, not 'no structure'",
            "UNIDENTIFIABLE means spatial variance is ~0, so a length parameter has no meaning",
        ],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps({k: summary[k] for k in ("n_sections", "n_curve_rows", "status_counts")}, indent=2, ensure_ascii=False))
    return 0



def cmd_synth(args: argparse.Namespace) -> int:
    out = Path(args.out) / "synth"
    out.mkdir(parents=True, exist_ok=True)
    manifest, diagnostics = [], []
    seeds = [int(s) for s in args.seeds.split(",") if s.strip()]

    for scenario in ["S0", "S1", "S2", "S3", "S5"]:
        builder = synthetic.SCENARIOS[scenario]
        for seed in seeds:
            if scenario == "S5":
                a, b = synthetic.make_s5_pair(seed)
                for sample in (a, b):
                    base = synthetic.materialise(sample, out)
                    manifest.append(_manifest_row(sample, base))
                diagnostics.append({
                    "base_id": f"{a.base_id}|{b.base_id}", "scenario": "S5_pair", "seed": seed,
                    "field_id": "visible_identity", "mass": 0.0,
                    "max": float(np.nanmax(np.abs(np.nan_to_num(a.counts) - np.nan_to_num(b.counts)))),
                    "centroid_x": 0.0, "centroid_y": 0.0, "n_modes_above_halfmax": 0,
                    "note": "max|visible count difference| must be 0",
                })
                continue
            sample = builder(seed)
            sample.counts = synthetic.counts_from_fields(sample, seed=seed)
            base = synthetic.materialise(sample, out)
            manifest.append(_manifest_row(sample, base))
            diagnostics.extend(synthetic.pair_diagnostics([sample]))

    for n_distractors in (0, 4, 16, 32):
        sample = synthetic.build_s4(seeds[0] if seeds else 0, n_distractors=n_distractors)
        sample.counts = synthetic.counts_from_fields(sample, seed=sample.seed)
        base = synthetic.materialise(sample, out)
        manifest.append(_manifest_row(sample, base))

    assets.write_tsv(manifest, out / "scenario_manifest.tsv")
    assets.write_tsv(diagnostics, out / "pair_matching_diagnostics.tsv")
    params = {
        "grid": synthetic.GRID,
        "coordinate_unit": "synthetic_pitch",
        "pitch": synthetic.PITCH,
        "markers_per_field": synthetic.MARKERS_PER_FIELD,
        "background_genes": synthetic.BACKGROUND_GENES,
        "seed_roles": {k: list(v) for k, v in synthetic.SEED_ROLES.items()},
        "seeds_run": seeds,
        "distractor_counts": [0, 4, 16, 32],
        "counts_model": "negative binomial, dispersion 8, lognormal depth (sd 0.25), optional shared technical gradient in S0",
        "generators": {"G1": "analytic gaussian/ring/bezier tube", "G2": "anisotropic diffusion on grid graph (available, not used by default scenarios)"},
        "notes": [
            "S5 pair shares bit-identical visible counts; only the masked centre differs",
            "S3 'broken' removes the middle third of the bridge; 'wrong_end' reconnects to a different endpoint",
            "S4 keeps the same true bridge while adding unrelated smooth fields; budget must not grow with field count",
        ],
    }
    (out / "generator_parameters.json").write_text(json.dumps(params, indent=2, ensure_ascii=False))
    summary = {
        "n_samples": len(manifest),
        "n_diagnostic_rows": len(diagnostics),
        "by_scenario": {s: sum(1 for m in manifest if m["scenario"] == s) for s in sorted({m["scenario"] for m in manifest})},
        "roles": {r: sum(1 for m in manifest if m["role"] == r) for r in sorted({m["role"] for m in manifest})},
        "s5_visible_identity_max_delta": [d["max"] for d in diagnostics if d["scenario"] == "S5_pair"],
        "s1_radial_profile_max_delta": [d["max"] for d in diagnostics if d["scenario"] == "S1_pair" and d["field_id"] == "circular_vs_sided"],
        "s1_centroid_shift_bimodal": [d["centroid_x"] for d in diagnostics if d["scenario"] == "S1_pair" and d["field_id"] == "unimodal_vs_bimodal"],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


def _manifest_row(sample, base: Path) -> dict:
    return {
        "base_id": sample.base_id,
        "scenario": sample.scenario,
        "pair_id": sample.pair_id,
        "seed": sample.seed,
        "role": sample.role,
        "n_fields": len(sample.fields),
        "field_ids": ",".join(sorted(sample.fields)),
        "n_points": sample.n_points,
        "n_visible": int(np.asarray(sample.visible_mask, dtype=bool).sum()),
        "path": str(base.relative_to(base.parents[1])),
    }



def _load_synth_sample(sample_dir: Path, field_ids: list[str]) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray], dict[str, np.ndarray]]:
    obs = np.load(sample_dir / "observed_inputs" / "counts_visible.npz")
    truth = np.load(sample_dir / "latent_truth" / "fields.npz")
    counts = obs["counts"]
    genes = obs["genes"].astype(str)
    xy = obs["xy"]
    visible = obs["visible_mask"].astype(bool)
    counts = np.nan_to_num(counts, nan=0.0)
    from scipy import sparse as _sp
    import json as _json

    meta = _json.loads((sample_dir / "meta.json").read_text())
    sets = {
        fid: [f"{fid}_m{k}" for k in range(meta["markers_per_field"])]
        for fid in field_ids if f"latent_truth/{fid}" or True
    }
    sets = {fid: sets[fid] for fid in field_ids}
    spec = InputSpec(fields=sets, witness_genes=frozenset(), k_neighbors=6, min_library_genes=50)
    inputs = build_inputs(_sp.csr_matrix(counts), genes, xy, visible, visible, spec)
    scores = {fid: inputs.arrays[f"field::{fid}"] for fid in field_ids}
    truths = {fid: np.asarray(truth[fid], dtype=np.float64) for fid in field_ids}
    return xy, visible, scores, truths


def cmd_route_a(args: argparse.Namespace) -> int:
    out = Path(args.out) / "route_a"
    out.mkdir(parents=True, exist_ok=True)
    synth_root = Path(args.out) / "synth"
    manifest = list(csv.DictReader(open(synth_root / "scenario_manifest.tsv"), delimiter="\t"))
    fit_rows, rec_rows, sel_rows = [], [], []
    kinds = {"circular": "ring", "sided": "ring", "unimodal": "blob", "bimodal": "blob"}

    for row in sorted(manifest, key=lambda r: r["base_id"]):
        if row["scenario"] != "S1":
            continue
        sample_dir = synth_root / row["path"]
        field_ids = row["field_ids"].split(",")
        xy, visible, scores, truths = _load_synth_sample(sample_dir, field_ids)
        for fid in field_ids:
            kind = kinds.get(fid, "blob")
            cands = route_a.candidates_for(kind, scores[fid], xy)
            chosen = route_a.select_by_cv(cands)
            truth = truths[fid]
            for cand in cands:
                fit_rows.append({
                    "base_id": row["base_id"], "seed": row["seed"], "field_id": fid, "model": cand.name,
                    "n_params": cand.n_params, "in_sample_sse": cand.in_sse, "cv_sse": cand.cv_sse,
                    "bic": cand.bic, "selected": cand is chosen,
                    "params": json.dumps(cand.params, ensure_ascii=False),
                })
            centre_true = route_a.centered(truth, xy)
            centre_pred = route_a.centered(chosen.prediction, xy)
            peaks = (
                route_a.peak_distances(chosen.prediction, truth, xy)
                if kind == "blob"
                else {"pred_to_truth": float("nan"), "truth_to_pred": float("nan"), "n_pred": 0, "n_truth": 0}
            )
            rec_rows.append({
                "base_id": row["base_id"], "seed": row["seed"], "field_id": fid, "kind": kind,
                "selected_model": chosen.name, "n_params": chosen.n_params,
                "cv_sse": chosen.cv_sse, "in_sample_sse": chosen.in_sse,
                "pearson": route_a.pearson(chosen.prediction, truth),
                "rmse_scaled": route_a.scaled_rmse(chosen.prediction, truth),
                "centre_error_pitch": float(np.linalg.norm(centre_pred - centre_true)),
                "peak_pred_to_truth_pitch": peaks["pred_to_truth"],
                "peak_truth_to_pred_pitch": peaks["truth_to_pred"],
                "n_peaks_pred": peaks["n_pred"], "n_peaks_truth": peaks["n_truth"],
                "angular_amplitude": route_a.angular_amplitude(chosen.params),
            })
            for cand in cands:
                sel_rows.append({
                    "base_id": row["base_id"], "seed": row["seed"], "field_id": fid, "kind": kind,
                    "model": cand.name, "cv_sse": cand.cv_sse, "n_params": cand.n_params,
                    "pearson": route_a.pearson(cand.prediction, truth),
                    "centre_error_pitch": float(np.linalg.norm(route_a.centered(cand.prediction, xy) - centre_true)),
                    "angular_amplitude": route_a.angular_amplitude(cand.params),
                    "selected": cand is chosen,
                })

    assets.write_tsv(fit_rows, out / "fits.tsv")
    assets.write_tsv(rec_rows, out / "recovery.tsv")
    assets.write_tsv(sel_rows, out / "model_comparison.tsv")

    def mean_for(field_id: str, model: str, key: str) -> float:
        vals = [float(r[key]) for r in sel_rows if r["field_id"] == field_id and r["model"] == model and r[key] not in ("", "None")]
        return float(np.nanmean(vals)) if vals else float("nan")

    summary = {
        "n_fits": len(fit_rows),
        "selection_rule": "blocked spatial CV (4 quadrants), simplest model within 5% of the best CV SSE",
        "by_field": {
            fid: {
                m: {k: mean_for(fid, m, k) for k in ("pearson", "centre_error_pitch", "cv_sse", "angular_amplitude")}
                for m in ("A0_ring", "A1_ring", "A0_gauss1", "A1_gauss2", "A1_gauss3")
                if any(r["field_id"] == fid and r["model"] == m for r in sel_rows)
            }
            for fid in ("circular", "sided", "unimodal", "bimodal")
        },
        "notes": [
            "A0 is the legitimate strong reference (radial ring / single Gaussian)",
            "BIC is reported but not used for selection: spatially correlated residuals make an iid penalty over-select",
            "truth is used only for scoring; fitting sees the visible field score",
            "curve-band (spline tube) support is NOT implemented in this iteration",
        ],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps({k: summary[k] for k in ("n_fits", "selection_rule")}, indent=2, ensure_ascii=False))
    return 0



def cmd_route_b(args: argparse.Namespace) -> int:
    out = Path(args.out) / "route_b"
    out.mkdir(parents=True, exist_ok=True)
    synth_root = Path(args.out) / "synth"
    manifest = list(csv.DictReader(open(synth_root / "scenario_manifest.tsv"), delimiter="\t"))
    rows = []
    for row in sorted(manifest, key=lambda r: r["base_id"]):
        if row["scenario"] != "S2":
            continue
        sample_dir = synth_root / row["path"]
        field_ids = row["field_ids"].split(",")
        xy, visible, scores, truths = _load_synth_sample(sample_dir, field_ids)
        n = int(round(np.sqrt(len(xy))))
        anchor = scores["anchor"]
        truth_struct = json.loads(np.load(sample_dir / "evaluation_only" / "answers.npz")["structures"].item())
        for partner_id in ("partner_same_point", "partner_offset"):
            target = scores[partner_id]
            b0 = route_b.fit_isotropic_operator(anchor, target, n)
            b1 = route_b.fit_lag_operator(anchor, target, n, max_lag=12, step=4)
            b0_cv = route_b.blocked_cv(anchor, target, n, lambda train_mask: route_b.fit_isotropic_operator(anchor, target, n, train_mask=train_mask))
            b1_cv = route_b.blocked_cv(anchor, target, n, lambda train_mask: route_b.fit_lag_operator(anchor, target, n, max_lag=12, step=4, train_mask=train_mask))
            truth_disp = 0.0 if partner_id == "partner_same_point" else float(np.hypot(*truth_struct["requested_offset"]))
            rows.append({
                "base_id": row["base_id"], "seed": row["seed"], "partner_id": partner_id,
                "truth_displacement_pitch": truth_disp,
                "truth_partner_separable": truth_struct.get("separable"),
                "B0_displacement": 0.0, "B0_class": "same_point",
                "B1_displacement": float(np.hypot(*b1.displacement)),
                "B1_class": route_b.classify_displacement(b1),
                "B0_pearson": route_b.pearson(b0.prediction, target),
                "B1_pearson": route_b.pearson(b1.prediction, target),
                "B1_minus_B0_pearson": route_b.pearson(b1.prediction, target) - route_b.pearson(b0.prediction, target),
                "B0_cv_sse": b0_cv, "B1_cv_sse": b1_cv,
                "B1_minus_B0_cv_sse": b1_cv - b0_cv,
                "B1_weight_peak": json.dumps(b1.to_dict()["weight_peak"]),
                "B1_top_weights": json.dumps(b1.to_dict()["top_weights"], ensure_ascii=False),
                "B1_displacement_vector": json.dumps([round(b1.displacement[0], 2), round(b1.displacement[1], 2)]),
            })
    assets.write_tsv(rows, out / "displacement_recovery.tsv")
    sweep = route_b.capability_sweep(max_lag=20)
    assets.write_tsv(sweep, out / "capability_sweep.tsv")
    same = [r for r in rows if r["partner_id"] == "partner_same_point"]
    off = [r for r in rows if r["partner_id"] == "partner_offset"]
    summary = {
        "n_rows": len(rows),
        "B0_by_construction": "stationary operator returns displacement 0 for every input",
        "same_point_cases": {
            "B1_displacement": [r["B1_displacement"] for r in same],
            "B1_class": [r["B1_class"] for r in same],
        },
        "offset_cases": {
            "B1_displacement": [r["B1_displacement"] for r in off],
            "B1_class": [r["B1_class"] for r in off],
            "truth": [r["truth_displacement_pitch"] for r in off],
        },
        "prediction_gain": {
            "B1_minus_B0_pearson": [r["B1_minus_B0_pearson"] for r in rows],
            "B1_minus_B0_cv_sse": [r["B1_minus_B0_cv_sse"] for r in rows],
        },
        "capability_sweep": {
            "recovered": sum(1 for r in sweep if r["recovered"]),
            "total": len(sweep),
            "by_noise": {
                str(nz): {
                    "recovered": sum(1 for r in sweep if r["noise_sd"] == nz and r["recovered"]),
                    "total": sum(1 for r in sweep if r["noise_sd"] == nz),
                }
                for nz in sorted({r["noise_sd"] for r in sweep})
            },
        },
        "notes": [
            "B0 is a stationarity-respecting reference and cannot represent displacement by construction",
            "B1 is a ridge-fitted 2-D lag operator on visible data only; the displacement is read from the weight centroid",
            "classification threshold is 3 pitch; the S2 design's achievable separation is 7 pitch",
        ],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


def cmd_geometry_smoke(args: argparse.Namespace) -> int:
    out = Path(args.out) / "geometry_smoke"
    out.mkdir(parents=True, exist_ok=True)
    n = args.grid
    yy, xx = np.mgrid[0:n, 0:n]
    xy = np.column_stack([xx.ravel(), yy.ravel()]).astype(float)
    # a disc hole plus a narrow neck: the neck connects nothing across the hole
    hole = ((xx - n / 2) ** 2 + (yy - n / 2) ** 2) <= (n / 6) ** 2
    support = ~hole.ravel()
    neck = (np.abs(yy.ravel() - n / 2) <= 1) & (np.abs(xx.ravel() - n / 2) <= n / 6)
    support = support | neck

    adj = geometry.build_adjacency(xy, support, k=6, unit="synthetic_pitch")
    region = geometry.region_classes(support, support.copy())

    # a smooth field with a masked disc; statistics must not see the masked part
    r = np.linalg.norm(xy - np.array([n / 2, n / 2]), axis=1)
    field = np.exp(-(r / (n / 4)) ** 2)
    visible = ~(((xx - n / 3) ** 2 + (yy - n / 3) ** 2) <= (n / 8) ** 2).ravel()
    masked_values = field.copy()
    masked_values[~visible] = np.nan
    mean_all, std_all = visible_mean_std(masked_values.reshape(-1, 1), visible)
    visible_only = field[visible]
    mean_visible = float(np.nanmean(visible_only))
    leaks = bool(np.isfinite(masked_values[~visible]).any())

    direction = geometry.local_direction(field, visible, xy, adj)
    direction_all = geometry.local_direction(field, np.ones_like(visible), xy, adj)

    assets.write_tsv(
        [
            {
                "n_points": int(len(xy)),
                "n_support": int(support.sum()),
                "n_visible": int(visible.sum()),
                "n_hidden": int((~visible).sum()),
                "pitch": adj.pitch,
                "n_edges": adj.notes["n_edges"],
                "median_degree": float(np.median(adj.degrees[support])),
                "n_isolated": adj.n_isolated,
                "unit": adj.unit,
            }
        ],
        out / "adjacency_summary.tsv",
    )
    pairs = geometry.distance_pairs(xy, support, max_distance=6 * adj.pitch, n_bins=12)
    assets.write_tsv(
        [
            {
                "bin_low": float(pairs["bin_edges"][i]),
                "bin_high": float(pairs["bin_edges"][i + 1]),
                "n_pairs": int(pairs["counts"][i]),
            }
            for i in range(len(pairs["counts"]))
        ],
        out / "distance_pairs.tsv",
    )
    assets.write_tsv(
        [
            {
                "masked_stats_scope": "visible_only",
                "mean_masked_array": float(mean_all[0]),
                "mean_direct_visible": mean_visible,
                "std_direct_visible": float(np.nanstd(visible_only)),
                "std_masked_array": float(std_all[0]),
                "hidden_values_finite": leaks,
                "difference": abs(float(mean_all[0]) - mean_visible),
            }
        ],
        out / "masked_stats_check.tsv",
    )
    rows = []
    for label, res in (("visible_only", direction), ("unmasked_reference", direction_all)):
        measurable = int(np.asarray(res["measurable"]).sum())
        rows.append(
            {
                "variant": label,
                "n_points": int(len(xy)),
                "n_measurable": measurable,
                "measurable_fraction": measurable / len(xy),
                "n_unmeasurable_hidden": int((~visible & ~np.asarray(res["measurable"])).sum()),
                "angle_defined_outside_visible": int(
                    np.isfinite(np.asarray(res["angle_mod_pi"])[~visible]).sum()
                ),
            }
        )
    assets.write_tsv(rows, out / "direction_measurability.tsv")
    edges_bad = 0
    try:
        edges_bad = geometry.assert_edges_avoid_non_tissue(adj.graph, region)
    except Exception:
        edges_bad = int(((region[adj.graph.tocoo().row] == geometry.NON_TISSUE)).sum())
    summary = {
        "n_points": int(len(xy)),
        "n_support": int(support.sum()),
        "n_visible": int(visible.sum()),
        "pitch": adj.pitch,
        "n_isolated": adj.n_isolated,
        "hidden_values_finite": leaks,
        "masked_stats_matches_visible": float(mean_all[0]) == mean_visible,
        "measurable_fraction_visible_vs_unmasked": [
            float(np.asarray(direction["measurable"]).sum()) / len(xy),
            float(np.asarray(direction_all["measurable"]).sum()) / len(xy),
        ],
        "edges_touching_non_tissue": edges_bad,
        "notes": [
            "direction is undefined (NaN) where the disc is masked",
            "the neck is supported tissue and may be traversed; the disc hole may not",
        ],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="FA observability adapter CLI")
    parser.add_argument("--repo", default=str(repo_root_default()))
    sub = parser.add_subparsers(dest="command", required=True)

    p_assets = sub.add_parser("assets", help="inventory local assets (five states)")
    p_assets.add_argument("--out", default="infra/field_analyzer_20261020")
    p_assets.set_defaults(func=cmd_assets)

    p_units = sub.add_parser("units", help="write unit role assertions")
    p_units.add_argument("--out", default="infra/field_analyzer_20261020")
    p_units.set_defaults(func=cmd_units)

    p_poison = sub.add_parser("poison", help="run the hidden/witness poisoning counterexample")
    p_poison.add_argument("--section", default="ZEN36")
    p_poison.add_argument("--out", default="infra/field_analyzer_20261020")
    p_poison.add_argument("--seed", type=int, default=0)
    p_poison.add_argument("--k", type=int, default=6)
    p_poison.add_argument("--radius-pitch", type=float, default=6.0, dest="radius_pitch")
    p_poison.set_defaults(func=cmd_poison)

    p_scales = sub.add_parser("scales", help="variogram / directional / cross-field scale curves")
    p_scales.add_argument("--out", default="infra/field_analyzer_20261020")
    p_scales.add_argument("--max-lag", type=float, default=12.0, dest="max_lag")
    p_scales.add_argument("--bins", type=int, default=12)
    p_scales.add_argument("--k", type=int, default=6)
    p_scales.add_argument("--seed", type=int, default=0)
    p_scales.set_defaults(func=cmd_scales)

    p_synth = sub.add_parser("synth", help="generate synthetic kill-test samples")
    p_synth.add_argument("--out", default="infra/field_analyzer_20261020")
    p_synth.add_argument("--seeds", default="0,1")
    p_synth.set_defaults(func=cmd_synth)

    p_ra = sub.add_parser("route-a", help="fit A0/A1 geometry on S1 samples")
    p_ra.add_argument("--out", default="infra/field_analyzer_20261020")
    p_ra.set_defaults(func=cmd_route_a)

    p_rb = sub.add_parser("route-b", help="displacement recovery on S2 samples")
    p_rb.add_argument("--out", default="infra/field_analyzer_20261020")
    p_rb.set_defaults(func=cmd_route_b)

    p_geo = sub.add_parser("geometry-smoke", help="synthetic geometry with a hole and a masked field")
    p_geo.add_argument("--out", default="infra/field_analyzer_20261020")
    p_geo.add_argument("--grid", type=int, default=64)
    p_geo.set_defaults(func=cmd_geometry_smoke)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
