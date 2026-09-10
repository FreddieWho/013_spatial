#!/usr/bin/env python3
"""Exploratory robustness deep dive for the fold-0 TLS shared-ecology signal.

For each outer fold and representation rank, refits the shared-ecology ridge
score per inner patient fold (same math as the panel readout) and reports:
per-direction TLS shared AUCs, cosine similarity of the fitted coefficient
vectors in Z-space, gene-space weight profiles recovered exactly from the
centered effect matrix, and cross-fold gene-profile comparison.

Exploratory grade: descriptive only; N=2 paired patients per fold cannot
support confirmatory claims.  Recomputed AUCs are cross-checked against the
frozen deepdive panel and fail closed on mismatch.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from r04.field_exports import read_heldout_field_export
from r04.runtime import atomic_json
from r04.structure_gt import labels_for_spots, load_training_gt_catalog
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
from scripts.r04_run_structure_readout import (
    _balanced_patient_weights,
    _inner_patient_folds,
)

STRUCTURES = ("TLS", "TUMOR_STROMA_BOUNDARY")
RIDGE_ALPHA = 1.0


def _read_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _paired_subset(train_z, labels, patient_ids):
    mask = np.isfinite(np.vstack([labels[name] for name in STRUCTURES])).all(axis=0)
    if not np.any(mask):
        raise ValueError("no paired annotated spots")
    return (
        np.asarray(train_z, dtype=float)[mask],
        {name: np.asarray(labels[name], dtype=float)[mask] for name in STRUCTURES},
        np.asarray(patient_ids).astype(str)[mask],
    )


def _has_target_variation(paired_labels, train, test) -> bool:
    return all(
        len(np.unique(paired_labels[name][train])) >= 2
        and len(np.unique(paired_labels[name][test])) >= 2
        for name in STRUCTURES
    )


def shared_fit(train_z, train_labels, weights):
    """Replicate the panel shared-score fit; return beta and intermediates."""
    train_z = np.asarray(train_z, dtype=float)
    w = _weights(weights, len(train_z), name="training_weights")
    z_std, _ = _standardize_features(train_z, train_z, w)
    standardized = []
    for name in STRUCTURES:
        y = np.asarray(train_labels[name], dtype=float)
        mean = float(np.average(y, weights=w))
        scale = float(np.sqrt(np.average((y - mean) ** 2, weights=w)))
        standardized.append((y - mean) / (scale if scale > 1e-12 else 1.0))
    shared_target = np.mean(np.vstack(standardized), axis=0)
    _, beta, _ = _ridge_fit_predict(
        z_std, z_std, shared_target, alpha=RIDGE_ALPHA, weights=w
    )
    mean_x = _weighted_mean(z_std, w)
    mean_y = float(_weighted_mean(shared_target, weights=w))
    feature_mean = _weighted_mean(train_z, w)
    feature_var = _weighted_mean((train_z - feature_mean) ** 2, w)
    feature_scale = np.sqrt(feature_var)
    feature_scale = np.where(feature_scale > 1e-12, feature_scale, 1.0)
    return np.asarray(beta, dtype=float), z_std, shared_target, mean_x, mean_y, feature_scale


def shared_tls_auc(train_z, test_z, train_labels, test_labels, w_train, w_test):
    """TLS AUC of the shared score through the exact panel readout chain."""
    beta, z_std, shared_target, mean_x, mean_y, _ = shared_fit(train_z, train_labels, w_train)
    _, test_std = _standardize_features(
        np.asarray(train_z, dtype=float), np.asarray(test_z, dtype=float), w_train)
    intercept = mean_y - mean_x @ beta
    train_score = z_std @ beta + intercept
    test_score = test_std @ beta + intercept
    y_train = np.asarray(train_labels["TLS"], dtype=float)
    y_test = np.asarray(test_labels["TLS"], dtype=float)
    pred, _ = _one_dimensional_readout(train_score, test_score, y_train, w_train)
    _, _, _, _, _, feature_scale = shared_fit(train_z, train_labels, w_train)
    return (float(_auc(y_test, pred, w_test)), np.asarray(beta, dtype=float),
            np.asarray(feature_scale, dtype=float))


def recover_gene_projection(centered_effect, z_matrix):
    """Exact gene-space projection P with Z = Xc @ P via least squares."""
    p_hat, residuals, rank, _ = np.linalg.lstsq(centered_effect, z_matrix, rcond=None)
    if rank < z_matrix.shape[1]:
        raise ValueError("effect matrix does not span the representation")
    rel = float(np.linalg.norm(centered_effect @ p_hat - z_matrix)
                / max(float(np.linalg.norm(z_matrix)), 1e-12))
    if rel > 1e-6:
        raise ValueError(f"gene projection recovery inexact: rel={rel:.2e}")
    return np.asarray(p_hat, dtype=float)


def cosine(a: np.ndarray, b: np.ndarray) -> float | None:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    if na <= 1e-12 or nb <= 1e-12:
        return None
    return float(a @ b / (na * nb))


def top_overlap(w1: np.ndarray, w2: np.ndarray, k: int) -> float:
    s1 = set(np.argsort(-np.abs(np.asarray(w1, dtype=float)))[:k].tolist())
    s2 = set(np.argsort(-np.abs(np.asarray(w2, dtype=float)))[:k].tolist())
    if not s1 and not s2:
        return 1.0
    return len(s1 & s2) / len(s1 | s2)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--effect-panel", type=Path, required=True)
    parser.add_argument("--deepdive-panel", type=Path, required=True)
    parser.add_argument("--registry", type=Path,
                        default=Path("infra/structure-registry/structure_instances.tsv"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ranks", default="1,2,3")
    parser.add_argument("--top-k", type=int, default=50)
    args = parser.parse_args()
    root = args.project_root.resolve()
    panel = _read_json(root / args.effect_panel)
    if panel.get("schema") != "r04.frozen_effect_export_panel.v1":
        raise SystemExit("effect panel schema mismatch")
    deepdive = _read_json(root / args.deepdive_panel)
    dd_by_fold = {}
    for entry in deepdive.get("entries", []):
        dd_by_fold[entry["fold"]] = entry
    catalog = load_training_gt_catalog(
        root / args.registry, root / str(panel["training_manifest"]["path"])
        if isinstance(panel.get("training_manifest"), dict)
        else root / str(panel["training_manifest"]),
    )
    ranks = [int(t.strip()) for t in str(args.ranks).split(",")]
    output: dict[str, object] = {
        "schema": "r04.explore_shared_robustness.v1",
        "status": "EXPLORATORY_SHARED_ROBUSTNESS",
        "evidence_grade": "exploratory_post_hoc_no_confirmatory_claim",
        "effect_panel": {"path": str(args.effect_panel)},
        "deepdive_panel": {"path": str(args.deepdive_panel)},
        "ridge_alpha": RIDGE_ALPHA,
        "top_k": args.top_k,
        "note": ("N=2 paired patients per fold: per-direction AUCs and "
                 "coefficient comparisons are descriptive only."),
        "entries": [],
    }
    output_path = root / args.output
    atomic_json(output_path, output)
    gene_weights_by_fold_rank: dict[tuple[int, int], np.ndarray] = {}
    for record in panel.get("entries", []):
        exports = record.get("structure_exports", [])
        if record.get("status") != "TRAINING_EFFECT_EXPORTED" or len(exports) != 1:
            raise SystemExit("effect panel entry is not a single training export")
        export = read_heldout_field_export(root / str(exports[0]["npz_path"]))
        gene_ids = [str(v) for v in list(export["evaluation_gene_ids"])]
        train_labels = {
            name: labels_for_spots(catalog, name, export["section_ids"], export["barcodes"])
            for name in STRUCTURES
        }
        dd_entry = dd_by_fold.get(record["fold"])
        if dd_entry is None:
            raise SystemExit(f"deepdive panel has no entry for fold {record['fold']}")
        entry: dict[str, object] = {
            "restart_index": record["restart_index"],
            "fold": record["fold"],
            "ranks": {},
        }
        for rank in ranks:
            train_z, _, representation = effect_coordinates(export, export, max_rank=rank)
            paired_z, paired_labels, paired_patients = _paired_subset(
                train_z, train_labels, np.asarray(export["patient_ids"])
            )
            paired_mask = np.isfinite(
                np.vstack([train_labels[n] for n in STRUCTURES])).all(axis=0)
            xc = _center_by_section(
                np.asarray(export["spatial_rate_effect_centered"], dtype=float),
                np.asarray(export["section_ids"]),
            )[paired_mask]
            if xc.shape[0] != len(paired_z):
                raise ValueError("centered effect and representation are not spot-aligned")
            proj = recover_gene_projection(xc, paired_z)
            folds = _inner_patient_folds(paired_patients)
            per_direction = []
            betas = []
            dd_ranks = dd_entry["ranks"][f"rank{rank}"]["fold_records"]
            for (inner_train, inner_test, test_patients), dd_fr in zip(folds, dd_ranks):
                if not _has_target_variation(paired_labels, inner_train, inner_test):
                    per_direction.append({"test_patients": test_patients,
                                          "status": "NOT_TESTABLE_NO_TARGET_VARIATION"})
                    betas.append(None)
                    continue
                w_train = _balanced_patient_weights(paired_patients[inner_train])
                w_test = _balanced_patient_weights(paired_patients[inner_test])
                tls_auc, beta, train_scale = shared_tls_auc(
                    paired_z[inner_train], paired_z[inner_test],
                    {n: paired_labels[n][inner_train] for n in STRUCTURES},
                    {n: paired_labels[n][inner_test] for n in STRUCTURES},
                    w_train, w_test,
                )
                want = dd_fr["readout"]["shared_only"]["metrics"]["TLS"]["auc"]
                if tls_auc is None or want is None or not np.isclose(
                        tls_auc, want, rtol=1e-9, atol=1e-12):
                    raise SystemExit(
                        f"AUC cross-check failed fold {record['fold']} rank{rank} "
                        f"{test_patients}: {tls_auc} != deepdive {want}")
                per_direction.append({"test_patients": test_patients, "status": "COMPUTED",
                                      "tls_shared_auc": float(tls_auc),
                                      "n_test": int(np.sum(inner_test)),
                                      "_train_scale": train_scale.tolist()})
                betas.append(np.asarray(beta, dtype=float))
            beta_cos = (cosine(betas[0], betas[1])
                        if len(betas) == 2 and betas[0] is not None and betas[1] is not None
                        else None)
            gene_w = []
            for beta, direction in zip(betas, per_direction):
                if beta is None:
                    gene_w.append(None)
                    continue
                # beta lives in standardized-Z space of its own training fold
                gene_w.append(proj @ (beta / direction["_train_scale"]))
            gene_cos = (cosine(gene_w[0], gene_w[1])
                        if len(gene_w) == 2 and gene_w[0] is not None and gene_w[1] is not None
                        else None)
            overlap = (top_overlap(gene_w[0], gene_w[1], args.top_k)
                       if len(gene_w) == 2 and gene_w[0] is not None and gene_w[1] is not None
                       else None)
            top_genes = []
            for beta_w, direction in zip(gene_w, per_direction):
                if beta_w is None:
                    top_genes.append(None)
                    continue
                order = np.argsort(-np.abs(beta_w))[:20]
                top_genes.append([{"gene": gene_ids[i], "weight": float(beta_w[i])}
                                  for i in order.tolist()])
            if all(g is not None for g in gene_w):
                gene_weights_by_fold_rank[(record["fold"], rank)] = np.mean(
                    np.vstack([g for g in gene_w]), axis=0)
            entry["ranks"][f"rank{rank}"] = {
                "singular_values": representation["singular_values"],
                "per_direction": per_direction,
                "beta_cosine_inner_folds": beta_cos,
                "gene_cosine_inner_folds": gene_cos,
                f"top{args.top_k}_gene_jaccard_inner_folds": overlap,
                "top20_genes_by_abs_weight_per_direction": top_genes,
                "gene_ids_sha256": __import__("hashlib").sha256(
                    "\n".join(gene_ids).encode()).hexdigest(),
            }
        output["entries"].append(entry)
        atomic_json(output_path, output)
    cross_fold = {}
    for rank in ranks:
        a = gene_weights_by_fold_rank.get((0, rank))
        b = gene_weights_by_fold_rank.get((4, rank))
        cross_fold[f"rank{rank}"] = cosine(a, b) if a is not None and b is not None else None
    output["cross_fold_gene_cosine_fold0_vs_fold4"] = cross_fold
    output["status"] = "EXPLORATORY_COMPLETE"
    atomic_json(output_path, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
