#!/usr/bin/env python3
"""Composition final audit with nested patient-level cross-fitting (D-108).

Unlike the exploratory smoke (pooled residualizer fit, transductive leak),
the residualizer ``effect ~ ilr(proxy)`` is fitted on inner-train patients
only and the same mapping transforms train and test rows. Test patients never
contribute to any residualizer parameter estimate.

Runs rank1/rank2 (primary)/full-K3 sensitivity, unadjusted vs nested-adjusted,
on existing training-role field exports. Old exploratory artifacts are kept;
this writes a new final-audit artifact. Uncertainty: leave-one-patient-out
empirical spread across inner folds plus spot-level bootstrap CI of the
test-side specific-increment deltas (fixed seed).

Proxies remain proxies (never cell fractions); ILC-class genes are refused.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from r04.composition import apply_residualizer, fit_residualizer
from r04.field_exports import read_heldout_field_export
from r04.runtime import atomic_json
from r04.structure_gt import labels_for_spots, load_training_gt_catalog
from r04.structure_readout import effect_coordinates, evaluate_shared_specific_readout
from scripts.r04_explore_marker_smoke import (
    align_rows,
    load_section_counts,
    module_scores,
)
from scripts.r04_run_structure_readout import (
    _balanced_patient_weights,
    _inner_patient_folds,
)

STRUCTURES = ("TLS", "TUMOR_STROMA_BOUNDARY")
DEFAULT_BOOTSTRAP_DRAWS = 200
BOOTSTRAP_SEED = 20260911
MIN_TRAIN_SPOTS = 10
N_BOOTSTRAP_DRAWS = DEFAULT_BOOTSTRAP_DRAWS


def _sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _readout(z_train, z_test, y_train, y_test, patients_train, patients_test,
             representation):
    return evaluate_shared_specific_readout(
        z_train, z_test, y_train, y_test,
        representation=representation,
        training_weights=_balanced_patient_weights(patients_train),
        heldout_weights=_balanced_patient_weights(patients_test))


def _summarize(readout):
    inc = readout["specific_increment"]
    return {
        "auc_delta": dict(inc["auc_delta_by_target"]),
        "mse_delta": float(inc["mean_mse_delta_specific_minus_shared"]),
    }


def _bootstrap_ci(z_train, z_test, y_train, y_test, patients_train, patients_test,
                  representation, seed):
    rng = np.random.default_rng(seed)
    n = len(z_test)
    draws = {"auc_delta": {n_: [] for n_ in STRUCTURES}, "mse_delta": []}
    for _ in range(N_BOOTSTRAP_DRAWS):
        idx = rng.integers(0, n, n)
        readout = _readout(z_train, np.asarray(z_test)[idx],
                           y_train, {k: np.asarray(v)[idx] for k, v in y_test.items()},
                           patients_train, np.asarray(patients_test)[idx],
                           representation)
        s = _summarize(readout)
        draws["mse_delta"].append(s["mse_delta"])
        for name in STRUCTURES:
            value = s["auc_delta"][name]
            if value is not None:
                draws["auc_delta"][name].append(value)
    ci = {}
    for name in STRUCTURES:
        values = np.asarray(draws["auc_delta"][name], dtype=float)
        ci[name] = ([float(np.quantile(values, 0.025)),
                     float(np.quantile(values, 0.975))] if len(values) else [None, None])
    mse = np.asarray(draws["mse_delta"], dtype=float)
    ci["mse_delta"] = [float(np.quantile(mse, 0.025)), float(np.quantile(mse, 0.975))]
    return {"draws": N_BOOTSTRAP_DRAWS, "seed": seed, "q025_q975": ci}


def _build_pseudo(root, args, export, export_sections, export_barcodes, universe,
                  class_genes):
    n_spots = len(export_barcodes)
    scores = np.full((n_spots, len(class_genes)), np.nan)
    for section_id in sorted(set(export_sections)):
        mat = None
        for cache_dir in [root / d for d in args.cache_dirs]:
            try:
                mat, meta_barcodes, meta_genes = load_section_counts(cache_dir, section_id)
                break
            except (ValueError, OSError, EOFError):
                continue
        if mat is None:
            raise SystemExit(f"section cache missing: {section_id}")
        if meta_genes != universe:
            raise SystemExit(f"gene order differs from universe: {section_id}")
        rows = align_rows(export_sections, export_barcodes, section_id, meta_barcodes)
        scores[rows] = module_scores(mat, args.gene_index, class_genes)
    if not np.isfinite(scores).all():
        raise SystemExit("module score matrix has gaps")
    pseudo = np.maximum(scores, 0.0) + 1e-6
    return pseudo / pseudo.sum(axis=1, keepdims=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--training-exports", nargs="+", type=Path, required=True,
                        help="one training_full.npz per fold, in fold order")
    parser.add_argument("--folds", default="0,4")
    parser.add_argument("--marker-proxy", type=Path, required=True)
    parser.add_argument("--registry", type=Path,
                        default=Path("infra/structure-registry/structure_instances.tsv"))
    parser.add_argument("--training-manifest", type=Path, required=True)
    parser.add_argument("--cache-dirs", nargs="+", type=Path, required=True)
    parser.add_argument("--gene-universe", type=Path, required=True)
    parser.add_argument("--symbol-map", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ranks", default="1,2,full")
    parser.add_argument("--classes", default="T,B,Mye,Epi,Stromal")
    parser.add_argument("--bootstrap-draws", type=int, default=DEFAULT_BOOTSTRAP_DRAWS)
    args = parser.parse_args()
    global N_BOOTSTRAP_DRAWS
    N_BOOTSTRAP_DRAWS = int(args.bootstrap_draws)
    root = args.project_root.resolve()
    folds = [int(t) for t in str(args.folds).split(",") if t.strip() != ""]
    if len(folds) != len(args.training_exports):
        raise SystemExit("folds and training-exports must align")
    wanted_classes = [c.strip() for c in args.classes.split(",") if c.strip()]
    if "ILC" in wanted_classes:
        raise SystemExit("ILC proxy is not viable and must not be gap-filled")

    proxy = _read_json(root / args.marker_proxy)
    if proxy.get("schema") != "r04.marker_proxy.v1":
        raise SystemExit("marker proxy schema mismatch")
    universe = [line.strip() for line in
                (root / args.gene_universe).read_text().splitlines() if line.strip()]
    ensg_to_sym: dict[str, str] = {}
    with open(root / args.symbol_map, newline="") as handle:
        for row in csv.reader(handle):
            if len(row) >= 2:
                ensg_to_sym.setdefault(row[0], row[1].upper())
    gene_index: dict[str, int] = {}
    for i, gene in enumerate(universe):
        sym = ensg_to_sym.get(gene.replace("DEPRECATED_", "").split(".")[0], "")
        if sym:
            gene_index.setdefault(sym, i)
    args.gene_index = gene_index
    class_genes = {}
    for major in wanted_classes:
        entry = proxy["classes"].get(major)
        if entry is None or not entry.get("viable"):
            raise SystemExit(f"class {major} is not viable in the proxy")
        class_genes[major] = entry["in_panel"]
    catalog = load_training_gt_catalog(root / args.registry, root / args.training_manifest)

    output: dict[str, object] = {
        "schema": "r04.composition_final_audit.v1",
        "status": "NESTED_COMPOSITION_AUDIT",
        "evidence_grade": "exploratory_post_hoc_no_confirmatory_claim",
        "decisions": ["D-105", "D-108"],
        "nesting": ("residualizer fitted on inner-train patients only; "
                    "same mapping transforms train and test"),
        "classes": wanted_classes,
        "inputs": {
            str(rel): _sha256(root / rel)
            for rel in [args.marker_proxy, args.training_manifest, args.gene_universe]
        },
        "symbol_map": {"path": str(args.symbol_map),
                         "sha256": _sha256(root / args.symbol_map)},
        "folds": [],
    }
    output_path = root / args.output
    atomic_json(output_path, output)
    for fold, export_rel in zip(folds, args.training_exports):
        export = read_heldout_field_export(root / export_rel)
        export_sections = [str(v) for v in list(export["section_ids"])]
        export_barcodes = [str(v) for v in list(export["barcodes"])]
        pseudo = _build_pseudo(root, args, export, export_sections, export_barcodes,
                               universe, class_genes)
        train_labels = {n: labels_for_spots(
            catalog, n, export["section_ids"], export["barcodes"]) for n in STRUCTURES}
        patient_ids = np.asarray(export["patient_ids"]).astype(str)
        fold_rec: dict[str, object] = {"fold": fold, "export": str(export_rel),
                                       "ranks": []}
        for token in [t.strip() for t in str(args.ranks).split(",")]:
            max_rank = None if token == "full" else int(token)
            rank_name = "full" if token == "full" else f"rank{token}"
            train_z, _, representation = effect_coordinates(export, export,
                                                            max_rank=max_rank)
            paired = np.isfinite(np.vstack([train_labels[n] for n in STRUCTURES])).all(axis=0)
            z = np.asarray(train_z, dtype=float)[paired]
            labels = {n: np.asarray(train_labels[n], dtype=float)[paired]
                      for n in STRUCTURES}
            patients = patient_ids[paired]
            comp = pseudo[paired]
            rank_rec: dict[str, object] = {
                "rank": rank_name, "n_paired_spots": int(paired.sum()),
                "paired_patients": sorted(set(patients.tolist())), "inner_folds": []}
            for inner_train, inner_test, test_patients in _inner_patient_folds(patients):
                fold_item: dict[str, object] = {
                    "test_patients": test_patients,
                    "n_train_spots": int(inner_train.sum()),
                    "n_test_spots": int(inner_test.sum())}
                if (any(len(np.unique(labels[n][inner_train])) < 2
                        or len(np.unique(labels[n][inner_test])) < 2
                        for n in STRUCTURES)
                        or int(inner_train.sum()) < MIN_TRAIN_SPOTS):
                    fold_item["status"] = "NOT_TESTABLE_NO_TARGET_VARIATION"
                    rank_rec["inner_folds"].append(fold_item)
                    continue
                y_train = {n: labels[n][inner_train] for n in STRUCTURES}
                y_test = {n: labels[n][inner_test] for n in STRUCTURES}
                base = _readout(z[inner_train], z[inner_test], y_train, y_test,
                                patients[inner_train], patients[inner_test],
                                representation)
                adj_train_dims, adj_test_dims = [], []
                for dim in range(z.shape[1]):
                    coef = fit_residualizer(z[inner_train][:, dim], comp[inner_train])
                    adj_tr, _ = apply_residualizer(z[inner_train][:, dim],
                                                   comp[inner_train], coef)
                    adj_te, _ = apply_residualizer(z[inner_test][:, dim],
                                                   comp[inner_test], coef)
                    adj_train_dims.append(adj_tr)
                    adj_test_dims.append(adj_te)
                adj = _readout(np.column_stack(adj_train_dims),
                               np.column_stack(adj_test_dims), y_train, y_test,
                               patients[inner_train], patients[inner_test],
                               representation)
                fold_item["status"] = "COMPUTED"
                fold_item["unadjusted"] = _summarize(base)
                fold_item["adjusted"] = _summarize(adj)
                fold_item["adjusted_bootstrap_ci"] = _bootstrap_ci(
                    np.column_stack(adj_train_dims), np.column_stack(adj_test_dims),
                    y_train, y_test, patients[inner_train], patients[inner_test],
                    representation, BOOTSTRAP_SEED + fold * 100 + len(rank_rec["inner_folds"]))
                rank_rec["inner_folds"].append(fold_item)
            fold_rec["ranks"].append(rank_rec)
        output["folds"].append(fold_rec)
        atomic_json(output_path, output)
    output["status"] = "NESTED_AUDIT_COMPLETE"
    atomic_json(output_path, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
