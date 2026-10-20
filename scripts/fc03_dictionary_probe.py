#!/usr/bin/env python3
"""Limited K=8 dictionary probe. Does not search K or reopen R-04."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from sklearn.decomposition import NMF

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.fc01_field_recovery import (  # noqa: E402
    HGB,
    PANEL,
    class_weight,
    fit_once,
    labeled_tls,
    ring_features,
)

OUT = ROOT / "infra/fc03_dictionary_probe_20261006"
OX = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
K = 8


def load_peripheral():
    spec = importlib.util.spec_from_file_location("tls_peripheral", ROOT / "scripts/tls_peripheral.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def panel_genes():
    ctrl = json.loads((ROOT / "infra/tls_multiscale_gobp_20261002/panel_controls.json").read_text())
    universe = sorted({gene for genes in ctrl.values() for gene in genes})
    return ctrl, universe, [gene for gene in universe if gene not in set(ctrl[OX])]


def load_gse(peripheral, universe):
    identity = (ROOT / "infra/bioinf-data-index/raw/user_gse175540_20261003/current_26_source_identity_confirmed.tsv").read_text().splitlines()
    header = identity[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in identity[1:] if line.strip() and line.split("\t")[1] == "GSE175540"]
    index = {gene: i for i, gene in enumerate(universe)}
    sections = []
    for row in rows:
        meta = peripheral.source_metadata(row["cohort"], row["section_id"])
        matrix, genes = peripheral.counts(row["cohort"], row["section_id"], meta)
        library = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float32)
        columns = np.array([index.get(str(gene), -1) for gene in genes])
        keep = columns >= 0
        dense = np.zeros((matrix.shape[0], len(universe)), dtype=np.float32)
        if keep.any():
            dense[:, columns[keep]] = matrix[:, keep].toarray()
        positive = library > 0
        dense[positive] = np.log1p(dense[positive] / library[positive, None] * 10000.0)
        geo = np.load(ROOT / "infra/tls_multiscale_gobp_20261002/geometry" / f"{row['section_id']}.npz")
        labels = np.array([str(x) for x in geo["labels"]])
        if len(labels) != len(dense):
            raise SystemExit(f"{row['section_id']} count/label mismatch")
        sections.append({
            "section_id": row["section_id"],
            "unit": row["patient_id"],
            "expr": dense,
            "labels": labels,
            "xy": np.asarray(geo["xy"], dtype=np.float64),
            "ox": np.load(ROOT / "infra/tls_multiscale_gobp_20261002/scores" / f"{row['section_id']}.npz")["scores"][:, PANEL.index(OX)].astype(np.float32),
        })
        print("counts", row["section_id"], dense.shape, flush=True)
    return sections


def fit_nmf(x):
    model = NMF(n_components=K, init="random", random_state=0, max_iter=100, solver="cd")
    weights = model.fit_transform(np.clip(x, 0, None))
    return model, weights.astype(np.float32)


def match_stability(a, b):
    used = set()
    scores = []
    for row in a:
        best, name = 0.0, None
        for i, other in enumerate(b):
            if i in used or np.std(row) <= 1e-8 or np.std(other) <= 1e-8:
                continue
            value = abs(float(np.corrcoef(row, other)[0, 1]))
            if np.isfinite(value) and value > best:
                best, name = value, i
        if name is not None:
            used.add(name)
            scores.append(best)
    return None if not scores else float(np.median(scores))


def stability(train, columns):
    units = sorted({sec["unit"] for sec in train})
    left = {unit for i, unit in enumerate(units) if i % 2 == 0}
    a = np.vstack([sec["expr"][:, columns] for sec in train if sec["unit"] in left])
    b = np.vstack([sec["expr"][:, columns] for sec in train if sec["unit"] not in left])
    if len(a) < 20 or len(b) < 20:
        return None
    model_a = NMF(n_components=K, init="random", random_state=0, max_iter=100, solver="cd").fit(np.clip(a, 0, None))
    model_b = NMF(n_components=K, init="random", random_state=0, max_iter=100, solver="cd").fit(np.clip(b, 0, None))
    return match_stability(model_a.components_, model_b.components_)


def predict_fold(train, test, columns, target):
    model, _ = fit_nmf(np.vstack([sec["expr"][:, columns] for sec in train]))
    prepared = []
    for sec in train + test:
        factors = model.transform(np.clip(sec["expr"][:, columns], 0, None)).astype(np.float32)
        means, empty = ring_features(sec["xy"], factors)
        prepared.append({"factors": factors, "rings": (means, empty), "labels": sec["labels"], "ox": sec["ox"], "section_id": sec["section_id"]})
    n_train = len(train)
    train_p, test_p = prepared[:n_train], prepared[n_train:]
    rows = []
    for model_name, spatial in (("hgb_spot", False), ("hgb_spatial", True)):
        def pack(items):
            blocks = []
            for item in items:
                spot = item["factors"]
                if spatial:
                    means, empty = item["rings"]
                    blocks.append(np.hstack([spot, *means, *[flag[:, None] for flag in empty]]))
                else:
                    blocks.append(spot)
            return np.vstack(blocks)
        x_train = pack(train_p)
        if target == "TLS":
            labels = np.concatenate([item["labels"] for item in train_p])
            mask = labeled_tls(labels)
            y = (labels[mask] == "TLS").astype(int)
            weight = class_weight(y)
            if weight is None:
                return [{"target": target, "model": model_name, "status": "NOT_FITTABLE"}]
            x_fit = x_train[mask]
            kind = "hgb_cls"
        else:
            y = np.concatenate([item["ox"] for item in train_p]).astype(np.float32)
            weight = None
            x_fit = x_train
            kind = "hgb_reg"
        preds, _ = fit_once(kind, x_fit, y, weight, [pack([item]) for item in test_p])
        values = []
        for item, pred in zip(test_p, preds):
            if target == "TLS":
                mask = labeled_tls(item["labels"])
                if mask.sum() == 0 or len(np.unique(item["labels"][mask])) < 2:
                    continue
                from sklearn.metrics import roc_auc_score
                values.append(float(roc_auc_score((item["labels"][mask] == "TLS").astype(int), pred[mask])))
            else:
                values.append(float(np.sqrt(np.mean((pred - item["ox"]) ** 2))))
        rows.append({"target": target, "model": model_name, "status": "OK" if values else "NOT_TESTABLE",
                     "metric": None if not values else float(sum(values) / len(values)), "n_eval_sections": len(values)})
    return rows


def main():
    if __import__("os").environ.get("CUDA_VISIBLE_DEVICES", "") not in ("", "-1"):
        raise SystemExit("GPU disabled")
    HGB["max_iter"] = 120
    OUT.mkdir(parents=True, exist_ok=True)
    peripheral = load_peripheral()
    ctrl, universe, ox_genes = panel_genes()
    gene_index = {gene: i for i, gene in enumerate(universe)}
    ox_columns = np.array([gene_index[gene] for gene in ox_genes])
    all_columns = np.arange(len(universe))
    sections = load_gse(peripheral, universe)
    units = sorted({sec["unit"] for sec in sections})
    results, stabilities = [], []
    for unit in units:
        test = [sec for sec in sections if sec["unit"] == unit]
        train = [sec for sec in sections if sec["unit"] != unit]
        if unit == "GSE175540:c_2" and len(test) != 2:
            raise SystemExit("c_2 split")
        for target, columns in ((OX, ox_columns), ("TLS", all_columns)):
            stable = stability(train, columns)
            stabilities.append({"unit": unit, "target": target, "median_matched_loading_correlation": stable,
                                "unstable": stable is None or stable < 0.5})
            for row in predict_fold(train, test, columns, target):
                row.update({"unit": unit, "cohort": "GSE175540"})
                results.append(row)
            print("fold", unit, target, "stability", stable, flush=True)
    (OUT / "metrics.json").write_text(json.dumps(results, indent=2) + "\n")
    (OUT / "stability.json").write_text(json.dumps(stabilities, indent=2) + "\n")
    (OUT / "STATUS").write_text("FC03_PROBE_COMPLETE\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
