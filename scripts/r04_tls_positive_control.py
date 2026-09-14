#!/usr/bin/env python3
"""TLS positive control: can in-panel immune markers predict TLS labels?

Known TLS biology (CXCL13 axis, mature B cells, T cells) is largely ABSENT
from the frozen 4000-gene panel (CXCL13/MS4A1/CD3D/CD19/CCL19/CCL21/LTB all
missing). What the panel DOES contain is a plasma module (IGHM/IGKC/JCHAIN/
MZB1/SDC1/IGLC3/IGHA2) plus PTPRC/CD68. This script tests, on training-role
spots with TLS annotations, whether that available signal predicts TLS.

Interpretation (pre-registered before running):
- plasma/PTPRC AUC clearly > 0.5  => TLS biology IS in the data; the R-04
  null is about cross-patient reproducible FIELDS, not absent biology.
- AUC ~ 0.5 for everything   => biology absent at this resolution/panel, or
  annotation problem; escalate to annotation audit, not to model tuning.

Exploratory grade. No model training, no held-out inference, no GT reuse
beyond scoring. Follows scripts/r04_explore_marker_smoke.py loading pattern.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from sklearn.linear_model import LogisticRegression

from r04.runtime import atomic_json
from r04.structure_gt import labels_for_spots, load_training_gt_catalog
from r04.structure_readout import _auc
from scripts.r04_explore_marker_smoke import (
    align_rows,
    load_section_counts,
)

PLASMA_GENES = ["IGHM", "IGKC", "JCHAIN", "MZB1", "SDC1", "IGLC3", "IGHA2"]
SINGLE_GENES = ["PTPRC", "CD68", "EPCAM", "COL1A1", "ACTA2", "LYZ"]
N_BOOT = 500
SEED = 20260912


def _uauc(target: np.ndarray, prediction: np.ndarray) -> float | None:
    target = np.asarray(target)
    return _auc(target, np.asarray(prediction), np.ones(len(target)))


def _read_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--training-manifest", type=Path, required=True)
    parser.add_argument("--registry", type=Path,
                        default=Path("infra/structure-registry/structure_instances.tsv"))
    parser.add_argument("--cache-dirs", nargs="+", type=Path, required=True)
    parser.add_argument("--gene-universe", type=Path, required=True)
    parser.add_argument("--symbol-list", type=Path, required=True,
                        help="positional SYMBOL list aligned to gene_universe")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()

    universe = [l.strip() for l in (root / args.gene_universe).read_text().splitlines()
                if l.strip()]
    symbols = [l.strip().upper() for l in (root / args.symbol_list).read_text().splitlines()]
    if len(symbols) != len(universe):
        raise SystemExit(f"symbol list length {len(symbols)} != universe {len(universe)}")
    gene_index = {}
    for i, s in enumerate(symbols):
        if s:
            gene_index.setdefault(s, i)
    missing = [g for g in PLASMA_GENES + SINGLE_GENES if g not in gene_index]
    if missing:
        raise SystemExit(f"expected genes absent from panel: {missing}")

    catalog = load_training_gt_catalog(root / args.registry, root / args.training_manifest)
    with open(root / args.training_manifest) as f:
        manifest = json.load(f)
    man_rows = manifest.get("rows", [])
    sec2pat = {e.get("section_id", ""): e.get("patient_id", "") for e in man_rows
               if e.get("section_id", "")}
    sections = sorted(sec2pat)
    if not sections:
        raise SystemExit("no sections found in training manifest")

    feats, labels, patients, sec_ids = [], [], [], []
    n_annot = n_total = 0
    for section_id in sections:
        mat = meta_barcodes = meta_genes = None
        for cache_dir in [root / d for d in args.cache_dirs]:
            try:
                mat, meta_barcodes, meta_genes = load_section_counts(cache_dir, section_id)
                break
            except (ValueError, OSError, EOFError):
                continue
        if mat is None or meta_genes != universe:
            continue
        log = mat.copy().tocsr()
        log.data = np.log1p(log.data)
        lab = np.asarray(labels_for_spots(catalog, "TLS", [section_id] * mat.shape[0],
                                          meta_barcodes),
                         dtype=float)
        finite = np.isfinite(lab)
        n_total += mat.shape[0]
        n_annot += int(finite.sum())
        if not finite.any():
            continue
        cols = []
        plasma_idx = sorted({gene_index[g] for g in PLASMA_GENES})
        cols.append(np.asarray(log[:, plasma_idx].mean(axis=1)).ravel())
        for g in SINGLE_GENES:
            cols.append(np.asarray(log[:, gene_index[g]].toarray()).ravel()
                        if sparse.issparse(log[:, gene_index[g]])
                        else np.asarray(log[:, gene_index[g]]).ravel())
        feats.append(np.column_stack(cols)[finite])
        labels.append(lab[finite])
        patients.extend([sec2pat.get(section_id, section_id)] * int(finite.sum()))
        sec_ids.extend([section_id] * int(finite.sum()))
    if not feats:
        raise SystemExit("no annotated training spots found")
    X = np.vstack(feats)
    y = np.concatenate(labels).astype(int)
    names = ["plasma_module"] + SINGLE_GENES
    if len(np.unique(y)) < 2:
        raise SystemExit("TLS labels lack both classes in training pool")

    rng = np.random.default_rng(SEED)
    results = {}
    for j, name in enumerate(names):
        auc = float(_uauc(y, X[:, j]))
        boot = []
        for _ in range(N_BOOT):
            idx = rng.integers(0, len(y), len(y))
            if len(np.unique(y[idx])) < 2:
                continue
            boot.append(float(_uauc(y[idx], X[idx, j])))
        lo, hi = (float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))) if boot else (None, None)
        results[name] = {"auc": auc, "ci95": [lo, hi], "n_boot": len(boot)}
    Xs = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-12)
    clf = LogisticRegression(max_iter=2000).fit(Xs, y)
    auc_all = float(_uauc(y, clf.predict_proba(Xs)[:, 1]))
    boot = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        try:
            c = LogisticRegression(max_iter=2000).fit(Xs[idx], y[idx])
            boot.append(float(_uauc(y[idx], c.predict_proba(Xs[idx])[:, 1])))
        except Exception:
            continue
    lo, hi = (float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))) if boot else (None, None)
    results["combined_logistic"] = {"auc": auc_all, "ci95": [lo, hi], "n_boot": len(boot)}

    output = {
        "schema": "r04.tls_positive_control.v1",
        "status": "EXPLORATORY_COMPLETE",
        "evidence_grade": "exploratory_post_hoc_no_confirmatory_claim",
        "n_annotated_spots": int(len(y)),
        "n_total_spots_scanned": int(n_total),
        "n_positive": int(y.sum()),
        "n_sections_with_labels": len(set(sec_ids)),
        "n_patients_with_labels": len(set(patients)),
        "features": names,
        "plasma_genes": PLASMA_GENES,
        "results": results,
    }
    atomic_json(root / args.output, output)
    print(f"annotated={len(y)} pos={int(y.sum())} sections={len(set(sec_ids))} patients={len(set(patients))}")
    for name in names + ["combined_logistic"]:
        r = results[name]
        print(f"  {name:18s} AUC={r['auc']:.3f} CI=[{r['ci95'][0]:.3f},{r['ci95'][1]:.3f}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
