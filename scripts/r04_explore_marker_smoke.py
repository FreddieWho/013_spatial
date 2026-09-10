#!/usr/bin/env python3
"""Marker-proxy composition smoke: adjust training coordinates, rerun readout.

Builds per-spot module scores from voted marker sets (log1p mean over
in-panel genes per class, row-normalized to pseudo-proportions), residualizes
each effect-coordinate dimension with patient-grouped cross-fitting
(r04.composition.crossfit_composition_adjustment), and reruns the identical
inner-crossfit shared/specific readout on adjusted coordinates.

Exploratory grade (D-105): the proxy is a crude stand-in, not a composition
measurement. Barcode alignment between panel caches and field exports is
exact-match fail-closed; gene order is asserted against gene_universe.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from r04.composition import crossfit_composition_adjustment
from r04.field_exports import read_heldout_field_export
from r04.runtime import atomic_json
from r04.structure_gt import labels_for_spots, load_training_gt_catalog
from r04.structure_readout import effect_coordinates, evaluate_shared_specific_readout
from scripts.r04_run_structure_readout import (
    _balanced_patient_weights,
    _inner_patient_folds,
)

STRUCTURES = ("TLS", "TUMOR_STROMA_BOUNDARY")


def _read_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_section_counts(cache_dir: Path, section_id: str) -> tuple:
    """Return (csr counts, metadata barcodes, metadata gene ids) for a section."""
    matches = sorted(cache_dir.glob(f"sections/*_{section_id}.npz"))
    if len(matches) != 1:
        raise ValueError(f"cache lookup ambiguous for {section_id}: {matches}")
    archive = np.load(matches[0], allow_pickle=False)
    mat = sparse.csr_matrix(
        (archive["data"], archive["indices"], archive["indptr"]),
        shape=tuple(archive["shape"]))
    meta = _read_json(matches[0].with_suffix(".json"))
    barcodes = [str(v) for v in meta["barcodes"]]
    if len(barcodes) != mat.shape[0]:
        raise ValueError(f"barcode/count row mismatch for {section_id}")
    return mat, barcodes, [str(v) for v in meta["gene_id"]]


def align_rows(export_sections: list[str], export_barcodes: list[str],
                section_id: str, meta_barcodes: list[str]) -> np.ndarray:
    """Map cache row order onto export positions via (section, barcode) keys."""
    order_index = {(s, b): i for i, (s, b) in enumerate(zip(export_sections, export_barcodes))}
    rows = []
    for barcode in meta_barcodes:
        key = (section_id, barcode)
        if key not in order_index:
            raise SystemExit(f"export/cache barcode mismatch: {section_id}::{barcode}")
        rows.append(order_index[key])
    return np.asarray(rows)


def module_scores(counts: sparse.csr_matrix, gene_index: dict[str, int],
                  class_genes: dict[str, list[str]]) -> np.ndarray:
    """Mean log1p expression over in-panel marker genes per class (spots x classes)."""
    log_counts = counts.copy()
    log_counts.data = np.log1p(log_counts.data)
    cols = []
    for major in sorted(class_genes):
        idx = sorted({gene_index[g] for g in class_genes[major] if g in gene_index})
        if not idx:
            raise ValueError(f"class {major} has no in-panel marker genes")
        cols.append(np.asarray(log_counts[:, idx].mean(axis=1)).ravel())
    return np.column_stack(cols)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--training-export", type=Path, required=True)
    parser.add_argument("--marker-proxy", type=Path, required=True)
    parser.add_argument("--registry", type=Path,
                        default=Path("infra/structure-registry/structure_instances.tsv"))
    parser.add_argument("--training-manifest", type=Path, required=True)
    parser.add_argument("--cache-dirs", nargs="+", type=Path, required=True)
    parser.add_argument("--gene-universe", type=Path, required=True)
    parser.add_argument("--symbol-map", type=Path, required=True,
                        help="ENSG->symbol CSV (same table used for proxy building)")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ranks", default="full,2")
    parser.add_argument("--classes", default="T,B,Mye,Epi,Stromal")
    args = parser.parse_args()
    root = args.project_root.resolve()
    wanted_classes = [c.strip() for c in args.classes.split(",") if c.strip()]

    export = read_heldout_field_export(root / args.training_export)
    proxy = _read_json(root / args.marker_proxy)
    if proxy.get("schema") != "r04.marker_proxy.v1":
        raise SystemExit("marker proxy schema mismatch")
    universe = [l.strip() for l in (root / args.gene_universe).read_text().splitlines()
                if l.strip()]
    import csv
    ensg_to_sym: dict[str, str] = {}
    with open(root / args.symbol_map, newline="") as f:
        for row in csv.reader(f):
            if len(row) >= 2:
                ensg_to_sym.setdefault(row[0], row[1].upper())
    gene_index: dict[str, int] = {}
    for i, g in enumerate(universe):
        sym = ensg_to_sym.get(g.replace("DEPRECATED_", "").split(".")[0], "")
        if sym:
            gene_index.setdefault(sym, i)
    class_genes = {}
    for major in wanted_classes:
        entry = proxy["classes"].get(major)
        if entry is None or not entry.get("viable"):
            raise SystemExit(f"class {major} is not viable in the proxy")
        class_genes[major] = entry["in_panel"]
    catalog = load_training_gt_catalog(root / args.registry, root / args.training_manifest)

    export_sections = [str(v) for v in list(export["section_ids"])]
    export_barcodes = [str(v) for v in list(export["barcodes"])]
    n_spots = len(export_barcodes)
    scores = np.full((n_spots, len(wanted_classes)), np.nan)
    seen_sections = set()
    for section_id in sorted(set(export_sections)):
        mat = meta_barcodes = meta_genes = None
        for cache_dir in [root / d for d in args.cache_dirs]:
            try:
                mat, meta_barcodes, meta_genes = load_section_counts(cache_dir, section_id)
                break
            except (ValueError, OSError, EOFError):
                # Corrupt/truncated shards (e.g. 0-byte writes; EOFError is NOT
                # an OSError subclass) are skipped;
                # a section with no readable shard anywhere fails below.
                continue
        if mat is None:
            raise SystemExit(f"section cache missing: {section_id}")
        if meta_genes != universe:
            raise SystemExit(f"gene order differs from universe: {section_id}")
        rows = align_rows(export_sections, export_barcodes, section_id, meta_barcodes)
        scores[rows] = module_scores(mat, gene_index, class_genes)
        seen_sections.add(section_id)
    if not np.isfinite(scores).all():
        raise SystemExit("module score matrix has gaps")
    pseudo = np.maximum(scores, 0.0) + 1e-6
    pseudo = pseudo / pseudo.sum(axis=1, keepdims=True)

    train_labels = {n: labels_for_spots(
        catalog, n, export["section_ids"], export["barcodes"]) for n in STRUCTURES}
    patient_ids = np.asarray(export["patient_ids"]).astype(str)
    output: dict[str, object] = {
        "schema": "r04.explore_marker_smoke.v1",
        "status": "EXPLORATORY_MARKER_SMOKE",
        "evidence_grade": "exploratory_post_hoc_no_confirmatory_claim",
        "decision": "D-105",
        "classes": wanted_classes,
        "module_score": "mean log1p over in-panel voted genes, row-normalized pseudo-proportions",
        "entries": [],
    }
    output_path = root / args.output
    atomic_json(output_path, output)
    for token in [t.strip() for t in str(args.ranks).split(",")]:
        max_rank = None if token == "full" else int(token)
        rank_name = "full" if token == "full" else f"rank{token}"
        train_z, _, representation = effect_coordinates(export, export, max_rank=max_rank)
        paired = np.isfinite(np.vstack([train_labels[n] for n in STRUCTURES])).all(axis=0)
        paired_z = np.asarray(train_z, dtype=float)[paired]
        paired_labels = {n: np.asarray(train_labels[n], dtype=float)[paired] for n in STRUCTURES}
        paired_patients = patient_ids[paired]
        paired_pseudo = pseudo[paired]
        adj_dims = []
        for dim in range(paired_z.shape[1]):
            adj = crossfit_composition_adjustment(
                paired_z[:, dim], paired_pseudo, groups=paired_patients)
            adj_dims.append(adj.adjusted)
        adj_z = np.column_stack(adj_dims)
        entry = {"rank": rank_name, "folds": []}
        for inner_train, inner_test, test_patients in _inner_patient_folds(paired_patients):
            fold_rec: dict[str, object] = {"test_patients": test_patients}
            for tag, z in (("unadjusted", paired_z), ("adjusted", adj_z)):
                if any(len(np.unique(paired_labels[n][inner_train])) < 2
                       or len(np.unique(paired_labels[n][inner_test])) < 2
                       for n in STRUCTURES):
                    fold_rec[tag] = {"status": "NOT_TESTABLE_NO_TARGET_VARIATION"}
                    continue
                readout = evaluate_shared_specific_readout(
                    z[inner_train], z[inner_test],
                    {n: paired_labels[n][inner_train] for n in STRUCTURES},
                    {n: paired_labels[n][inner_test] for n in STRUCTURES},
                    representation=representation,
                    training_weights=_balanced_patient_weights(paired_patients[inner_train]),
                    heldout_weights=_balanced_patient_weights(paired_patients[inner_test]))
                fold_rec[tag] = {
                    "status": "COMPUTED",
                    "auc_delta": readout["specific_increment"]["auc_delta_by_target"],
                    "mse_delta": readout["specific_increment"][
                        "mean_mse_delta_specific_minus_shared"],
                }
            entry["folds"].append(fold_rec)
        output["entries"].append(entry)
        atomic_json(output_path, output)
    output["status"] = "EXPLORATORY_COMPLETE"
    atomic_json(output_path, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
