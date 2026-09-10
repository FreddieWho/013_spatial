#!/usr/bin/env python3
"""Exploratory outer validation: training-fit shared score applied to held-out sections.

Frozen quantities (gene-space projection, ridge coefficients, 1-D readout
mappings) come only from training-role spots.  Held-out field effects are
reconstructed per gene split from the full-export inference exports and
concatenated to full-panel gene order; the training projection is applied
unchanged.  No refit touches held-out data.  Held-out GT is read for scoring
only and the result is exploratory grade.

Self-consistency gate: the same reconstruction applied to the training
export must reproduce the training coordinates (fail closed otherwise).
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from r04.field_exports import read_heldout_field_export
from r04.runtime import atomic_json
from r04.structure_gt import (
    _fragment_values,
    _read_annotation_csv,
    labels_for_spots,
    load_training_gt_catalog,
)
from r04.structure_readout import (
    _auc,
    _center_by_section,
    _one_dimensional_readout,
    _ridge_fit_predict,
    _standardize_features,
    _weighted_mean,
    _weights,
    effect_coordinates,
)
from scripts.r04_run_structure_readout import _balanced_patient_weights

STRUCTURES = ("TLS", "TUMOR_STROMA_BOUNDARY")
RIDGE_ALPHA = 1.0
RATE_CLIP = 30.0


def _read_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def reconstruct_split_effect(field_mean, evaluation_loading, factor_amplitude,
                             section_ids) -> np.ndarray:
    """Rate-path gene effect for one split (mirrors the readout rate branch)."""
    field = _center_by_section(
        np.exp(np.clip(np.asarray(field_mean, dtype=float), -RATE_CLIP, RATE_CLIP)),
        np.asarray(section_ids),
    )
    loading = np.asarray(evaluation_loading, dtype=float)
    amplitude = np.asarray(factor_amplitude, dtype=float)
    return field @ (loading * amplitude[None, :]).T


def assemble_full_effect(split_exports: list[dict]) -> tuple[np.ndarray, list[str],
                                                             np.ndarray, np.ndarray]:
    """Concatenate per-split eval-gene effects into full-panel gene order."""
    first_ids = np.asarray(split_exports[0]["section_ids"]).astype(str)
    first_bc = np.asarray(split_exports[0]["barcodes"]).astype(str)
    gene_ids = [str(v) for v in list(split_exports[0]["frozen_gene_ids"])]
    full = np.full((len(first_ids), len(gene_ids)), np.nan, dtype=float)
    for export in split_exports:
        if (np.asarray(export["section_ids"]).astype(str) != first_ids).any():
            raise ValueError("split exports have different spot order")
        with_section = reconstruct_split_effect(
            export["field_mean"], export["evaluation_loading"],
            export["factor_amplitude"], export["section_ids"])
        indices = np.asarray(export["evaluation_indices"], dtype=int)
        full[:, indices] = with_section
    if not np.isfinite(full).all():
        raise ValueError("assembled full effect has gaps or non-finite values")
    return full, gene_ids, first_ids, first_bc


def fit_shared_on_training(train_z, train_labels, weights) -> dict:
    """Shared ridge fit identical to the panel readout, plus frozen 1-D mappings."""
    w = _weights(weights, len(train_z), name="training_weights")
    z_std, _ = _standardize_features(np.asarray(train_z, dtype=float),
                                     np.asarray(train_z, dtype=float), w)
    standardized = []
    for name in STRUCTURES:
        y = np.asarray(train_labels[name], dtype=float)
        mean = float(np.average(y, weights=w))
        scale = float(np.sqrt(np.average((y - mean) ** 2, weights=w)))
        standardized.append((y - mean) / (scale if scale > 1e-12 else 1.0))
    shared_target = np.mean(np.vstack(standardized), axis=0)
    _, beta, _ = _ridge_fit_predict(z_std, z_std, shared_target,
                                    alpha=RIDGE_ALPHA, weights=w)
    beta = np.asarray(beta, dtype=float)
    mean_x = _weighted_mean(z_std, w)
    mean_y = float(_weighted_mean(shared_target, weights=w))
    intercept = float(mean_y - mean_x @ beta)
    train_score = z_std @ beta + intercept
    mappings = {}
    for name in STRUCTURES:
        y = np.asarray(train_labels[name], dtype=float)
        coef, *_ = np.linalg.lstsq(
            np.column_stack([np.ones(len(train_score)), train_score])
            * np.sqrt(w)[:, None],
            y * np.sqrt(w), rcond=None)
        mappings[name] = {"intercept": float(coef[0]), "slope": float(coef[1])}
    feature_mean = _weighted_mean(np.asarray(train_z, dtype=float), w)
    feature_var = _weighted_mean((np.asarray(train_z, dtype=float) - feature_mean) ** 2, w)
    feature_scale = np.sqrt(feature_var)
    feature_scale = np.where(feature_scale > 1e-12, feature_scale, 1.0)
    return {"beta": beta, "feature_mean": np.asarray(feature_mean, dtype=float),
            "feature_scale": np.asarray(feature_scale, dtype=float),
            "intercept": intercept, "mappings": mappings,
            "train_weights": np.asarray(w, dtype=float)}


def frozen_test_scores(test_z, frozen) -> np.ndarray:
    """Shared score on new spots with fully frozen parameters."""
    test_std = ((np.asarray(test_z, dtype=float) - frozen["feature_mean"])
                / frozen["feature_scale"])
    return test_std @ frozen["beta"] + frozen["intercept"]


def apply_frozen_shared(test_z, frozen, test_labels, test_weights) -> dict:
    """Apply frozen score + frozen mappings; report per-structure AUCs."""
    test_score = frozen_test_scores(test_z, frozen)
    out = {}
    for name in STRUCTURES:
        mapping = frozen["mappings"][name]
        pred = mapping["intercept"] + mapping["slope"] * test_score
        labels = np.asarray(test_labels[name], dtype=float)
        finite = np.isfinite(labels)
        w = _weights(test_weights, len(pred), name="heldout_weights")[finite]
        out[name] = _auc(labels[finite], pred[finite], w) if int(finite.sum()) else None
    return {"test_score": test_score, "auc_by_structure": out}


def load_heldout_gt_catalog(registry_path: Path, section_ids: list[str]) -> dict:
    """CONFIRMATORY registry rows for explicit held-out sections (exploratory use)."""
    wanted = {str(v) for v in section_ids}
    with registry_path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sections: dict[str, dict[str, dict]] = {"TLS": {}, "TUMOR_STROMA_BOUNDARY": {}}
    allowed: dict[str, list[str]] = {"TLS": [], "TUMOR_STROMA_BOUNDARY": []}
    skipped = 0
    for row in rows:
        structure = str(row.get("structure_id", ""))
        section_id = str(row.get("physical_unit_id", ""))
        if structure not in sections or section_id not in wanted:
            continue
        if row.get("confirmation_status") != "CONFIRMATORY":
            skipped += 1
            continue
        locator = str(row.get("gt_geometry_locator", ""))
        relative_path, positive_values = _fragment_values(locator)
        path = (relative_path if relative_path.is_absolute()
                else registry_path.parent.parent.parent / relative_path)
        if not path.is_file():
            path = registry_path.parent.parent.parent / relative_path
        annotation = _read_annotation_csv(path)
        existing = sections[structure].get(section_id)
        if existing is None:
            sections[structure][section_id] = {
                "path": str(path),
                "positive_values": sorted(set(positive_values)),
                "annotation": annotation,
            }
        else:
            if existing["path"] != str(path):
                raise ValueError(f"multiple GT files disagree for {structure}/{section_id}")
            existing["positive_values"] = sorted(
                set(existing["positive_values"]) | set(positive_values))
        allowed[structure].append(str(row.get("allowed_use", "")))
    return {"schema": "r04.structure_gt_catalog.v1",
            "structures": {n: {"sections": sections[n]} for n in sections},
            "allowed_use_seen": allowed, "skipped_nonconfirmatory": skipped,
            "role_boundary": "heldout_exploratory_scoring_only_never_an_input"}


def permutation_p(test_z, frozen, labels, weights, *, seed: int, draws: int,
                  group_ids) -> dict:
    """Within-section label-permutation null for the frozen outer AUCs.

    The annotated spot set per structure is held fixed; labels are permuted
    among annotated spots within each section."""
    rng = np.random.default_rng(seed)
    test_score = frozen_test_scores(test_z, frozen)
    groups = np.asarray(group_ids).astype(str)
    null: dict[str, list[float]] = {n: [] for n in STRUCTURES}
    vectors = {n: np.asarray(labels[n], dtype=float) for n in STRUCTURES}
    for _ in range(draws):
        for n in STRUCTURES:
            labeled = vectors[n]
            finite = np.flatnonzero(np.isfinite(labeled))
            if len(finite) == 0:
                continue
            order = np.arange(len(finite))
            sections = groups[finite]
            for value in np.unique(sections):
                block = np.flatnonzero(sections == value)
                order[block] = rng.permutation(block)
            permuted = labeled[finite][order]
            mapping = frozen["mappings"][n]
            pred = mapping["intercept"] + mapping["slope"] * test_score[finite]
            value = _auc(permuted, pred,
                         _weights(weights, len(test_score), name="heldout_weights")[finite])
            if value is not None:
                null[n].append(float(value))
    return null


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--effect-panel", type=Path, default=None)
    parser.add_argument("--training-export", type=Path, default=None,
                        help="direct mode: training_full.npz path (D-103 cross-device replay)")
    parser.add_argument("--heldout-export", type=Path, default=None,
                        help="direct mode: heldout_full.npz path (D-103 cross-device replay)")
    parser.add_argument("--split-exports", nargs=2, default=None, metavar=("SPLIT0", "SPLIT1"),
                        help="split-direct mode: two per-split inference npz paths")
    parser.add_argument("--training-manifest", type=Path, default=None)
    parser.add_argument("--expect-input-hash", type=str, default=None)
    parser.add_argument("--expect-config-hash", type=str, default=None)
    parser.add_argument("--fold", type=int, default=0)
    parser.add_argument("--registry", type=Path,
                        default=Path("infra/structure-registry/structure_instances.tsv"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ranks", default="1,2,3")
    parser.add_argument("--seed", type=int, default=20260905)
    parser.add_argument("--draws", type=int, default=200)
    args = parser.parse_args()
    root = args.project_root.resolve()
    split_direct = args.split_exports is not None
    direct = args.heldout_export is not None
    if direct and split_direct:
        raise SystemExit("--heldout-export and --split-exports are mutually exclusive")
    if split_direct and not (args.training_export and args.training_manifest):
        raise SystemExit("split-direct mode needs --training-export and --training-manifest")
    if split_direct and not (args.expect_input_hash and args.expect_config_hash):
        raise SystemExit("split-direct mode needs --expect-input-hash and --expect-config-hash")
    if split_direct:
        training_export = read_heldout_field_export(root / args.training_export)
        split_exports = []
        for sp in args.split_exports:
            split_exports.append(read_heldout_field_export(root / sp))
            meta = _read_json((root / sp).with_suffix(".json"))
            prov = meta.get("provenance", {})
            if prov.get("fit_input_hash") != args.expect_input_hash:
                raise SystemExit("split export input hash differs from expected source")
            if prov.get("fit_config_hash") != args.expect_config_hash:
                raise SystemExit("split export config hash differs from expected source")
            if prov.get("fit_updates") != 0:
                raise SystemExit("split export performed fit updates")
            if prov.get("inference_platform", {}).get("converged") is not True:
                raise SystemExit("split export inference platform not converged")
            if prov.get("validation_gt_read") != "NOT_READ":
                raise SystemExit("split export GT seal broken")
        heldout_export = None
        manifest_arg = root / args.training_manifest
        effect_panel_ref = {"mode": "split_direct_D103",
                            "training": str(args.training_export),
                            "splits": [str(sp) for sp in args.split_exports]}
    elif direct:
        if not (args.training_export and args.heldout_export and args.training_manifest):
            raise SystemExit("direct mode needs --training-export, --heldout-export, --training-manifest")
        if not (args.expect_input_hash and args.expect_config_hash):
            raise SystemExit("direct mode needs --expect-input-hash and --expect-config-hash")
        training_export = read_heldout_field_export(root / args.training_export)
        heldout_export = read_heldout_field_export(root / args.heldout_export)
        training_meta = _read_json((root / args.training_export).with_suffix(".json"))
        heldout_meta = _read_json((root / args.heldout_export).with_suffix(".json"))
        for label, meta in (("training", training_meta), ("heldout", heldout_meta)):
            prov = meta.get("provenance", {})
            if prov.get("fit_input_hash") != args.expect_input_hash:
                raise SystemExit(f"{label} export input hash differs from expected source")
            if prov.get("fit_config_hash") != args.expect_config_hash:
                raise SystemExit(f"{label} export config hash differs from expected source")
            if prov.get("fit_updates") != 0:
                raise SystemExit(f"{label} export performed fit updates")
            if prov.get("validation_gt_read") != "NOT_READ":
                raise SystemExit(f"{label} export GT seal broken")
        ho_plat = heldout_meta.get("provenance", {}).get("inference_platform", {})
        if ho_plat.get("converged") is not True:
            raise SystemExit("held-out export inference platform not converged")
        split_exports = None
        manifest_arg = root / args.training_manifest
        effect_panel_ref = {"mode": "direct_D103", "training": str(args.training_export),
                            "heldout": str(args.heldout_export)}
    else:
        if args.effect_panel is None:
            raise SystemExit("need --effect-panel or direct-mode export paths")
        panel = _read_json(root / args.effect_panel)
        if panel.get("schema") != "r04.frozen_effect_export_panel.v1":
            raise SystemExit("effect panel schema mismatch")
        records = [e for e in panel.get("entries", []) if e.get("fold") == args.fold]
        if len(records) != 1:
            raise SystemExit("effect panel must contain exactly one entry for the fold")
        record = records[0]
        if record.get("status") != "EXPORTED":
            raise SystemExit("entry is not a full export (need held-out inference exports)")
        struct_exports = record.get("structure_exports", [])
        infer_exports = record.get("exports", [])
        if len(struct_exports) != 1 or not infer_exports:
            raise SystemExit("entry export lists are incomplete")
        training_export = read_heldout_field_export(root / str(struct_exports[0]["npz_path"]))
        split_exports = [read_heldout_field_export(root / str(m["npz_path"])) for m in infer_exports]
        heldout_export = None
        manifest_arg = (root / str(panel["training_manifest"]["path"])
                        if isinstance(panel.get("training_manifest"), dict)
                        else root / str(panel["training_manifest"]))
        effect_panel_ref = {"path": str(args.effect_panel)}
    catalog = load_training_gt_catalog(root / args.registry, manifest_arg)
    ranks = [int(t.strip()) for t in str(args.ranks).split(",")]
    output: dict[str, object] = {
        "schema": "r04.explore_outer_validation.v1",
        "status": "EXPLORATORY_OUTER_VALIDATION",
        "evidence_grade": "exploratory_post_hoc_no_confirmatory_claim",
        "fold": args.fold,
        "seed": args.seed,
        "draws": args.draws,
        "effect_panel": effect_panel_ref,
        "cross_device_replay": "D-103" if (direct or split_direct) else None,
        "entries": [],
    }
    output_path = root / args.output
    atomic_json(output_path, output)
    for rank in ranks:
        train_z, _, representation = effect_coordinates(
            training_export, training_export, max_rank=rank)
        # Self-consistency gate: rebuild training coordinates from raw pieces.
        check_full, _, _, _ = assemble_full_effect([training_export])
        xc_check = _center_by_section(
            check_full, np.asarray(training_export["section_ids"]))
        proj_check, *_ = np.linalg.lstsq(xc_check, train_z, rcond=None)
        rel = float(np.linalg.norm(xc_check @ proj_check - train_z)
                    / max(float(np.linalg.norm(train_z)), 1e-12))
        if rel > 1e-6:
            raise SystemExit(f"training reconstruction gate failed: rel={rel:.2e}")
        # Held-out assembly on identical gene order.
        if direct:
            ho_full = np.asarray(heldout_export["spatial_rate_effect_centered"], dtype=float)
            ho_genes = [str(v) for v in list(heldout_export["evaluation_gene_ids"])]
            ho_section_ids = np.asarray(heldout_export["section_ids"]).astype(str)
            ho_barcodes = np.asarray(heldout_export["barcodes"]).astype(str)
            ho_patient_ids_direct = np.asarray(heldout_export["patient_ids"]).astype(str)
            if len(ho_genes) != len(set(ho_genes)):
                raise SystemExit("held-out evaluation genes are not unique")
        else:
            ho_full, ho_genes, ho_section_ids, ho_barcodes = assemble_full_effect(split_exports)
            ho_patient_ids_direct = np.asarray(split_exports[0]["patient_ids"]).astype(str)
        tr_genes = [str(v) for v in list(training_export["frozen_gene_ids"])]
        if sorted(ho_genes) != sorted(tr_genes):
            raise SystemExit("held-out and training gene universes differ")
        if direct:
            order = [ho_genes.index(g) for g in tr_genes]
            ho_full = ho_full[:, order]
        ho_xc = _center_by_section(ho_full, ho_section_ids)
        ho_z = ho_xc @ proj_check
        # Frozen shared fit on training paired spots.
        train_labels = {n: labels_for_spots(
            catalog, n, training_export["section_ids"], training_export["barcodes"])
            for n in STRUCTURES}
        mask = np.isfinite(np.vstack([train_labels[n] for n in STRUCTURES])).all(axis=0)
        paired_patients = np.asarray(training_export["patient_ids"]).astype(str)[mask]
        frozen = fit_shared_on_training(
            train_z[mask],
            {n: train_labels[n][mask] for n in STRUCTURES},
            _balanced_patient_weights(paired_patients))
        # Held-out GT (CONFIRMATORY rows only, scoring only).
        ho_catalog = load_heldout_gt_catalog(
            root / args.registry, [str(v) for v in ho_section_ids])
        ho_labels = {n: labels_for_spots(
            ho_catalog, n, ho_section_ids, ho_barcodes) for n in STRUCTURES}
        ho_patient_ids = ho_patient_ids_direct
        if len(ho_patient_ids) != len(ho_section_ids):
            raise SystemExit("held-out patient IDs are not spot-aligned")
        ho_weights = _balanced_patient_weights(ho_patient_ids)
        result = apply_frozen_shared(ho_z, frozen, ho_labels, ho_weights)
        by_section = {}
        for section in sorted(set(ho_section_ids.tolist())):
            block = np.flatnonzero(ho_section_ids == section)
            patient = str(ho_patient_ids[block[0]]) if len(block) else ""
            entry_sec = {"patient": patient, "n_spots": int(len(block))}
            for n in STRUCTURES:
                lab = np.asarray(ho_labels[n], dtype=float)[block]
                finite = np.isfinite(lab)
                mapping = frozen["mappings"][n]
                pred = (mapping["intercept"] + mapping["slope"]
                        * result["test_score"][block])
                entry_sec[n] = (float(_auc(lab[finite], pred[finite], ho_weights[block][finite]))
                                if int(finite.sum()) > 0 and len(np.unique(lab[finite])) == 2
                                else None)
                entry_sec[n + "_n_pos"] = int(np.sum(lab[finite] == 1)) if int(finite.sum()) else 0
            by_section[section] = entry_sec
        coverage = {}
        for n in STRUCTURES:
            finite = np.isfinite(ho_labels[n])
            coverage[n] = {"annotated_spots": int(finite.sum()),
                           "positive_spots": int(np.sum(ho_labels[n][finite] == 1)),
                           "sections": sorted(set(ho_section_ids[finite].tolist()))}
        null = permutation_p(ho_z, frozen, ho_labels, ho_weights,
                             seed=args.seed + 31 * rank, draws=args.draws,
                             group_ids=ho_section_ids)
        summary = {
            "rank": rank,
            "singular_values": representation["singular_values"],
            "n_training_paired_spots": int(mask.sum()),
            "n_heldout_spots": int(len(ho_section_ids)),
            "heldout_weighting": "patient_equal",
            "heldout_patients": sorted(set(ho_patient_ids.tolist())),
            "coverage_heldout": coverage,
            "by_section": by_section,
            "auc_by_structure": result["auc_by_structure"],
            "null_quantiles_auc": {
                n: ({q: float(np.quantile(v, float(q))) for q in ("0.025", "0.5", "0.975")}
                    if v else None)
                for n, v in null.items()},
            "empirical_upper_p_auc": {
                n: (float(np.mean(np.asarray(v) >= result["auc_by_structure"][n]))
                    if v and result["auc_by_structure"][n] is not None else None)
                for n, v in null.items()},
            "heldout_gt_role_boundary": ho_catalog["role_boundary"],
            "heldout_allowed_use_seen": ho_catalog["allowed_use_seen"],
        }
        output["entries"].append(summary)
        atomic_json(output_path, output)
    output["status"] = "EXPLORATORY_COMPLETE"
    atomic_json(output_path, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
