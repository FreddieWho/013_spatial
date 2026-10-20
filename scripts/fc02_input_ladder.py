#!/usr/bin/env python3
"""FC-02 input ladder under the frozen 2026-10-06 contract.

Feature selection uses only the training fold. Test targets are not used to
choose inputs.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.fc01_field_recovery import (
    HGB,
    PANEL,
    WORK as FC01,
    all_folds,
    class_weight,
    fit_once,
    labeled_tls,
    load_sections,
    stack,
)

OUT = Path(__file__).resolve().parents[1] / "infra/fc02_input_ladder_20261006"
OX = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"


def abs_corr(x, y):
    y = np.asarray(y, dtype=np.float64)
    out = np.zeros(x.shape[1], dtype=np.float64)
    if len(y) < 3 or np.std(y) <= 1e-8:
        return out
    for i in range(x.shape[1]):
        col = np.asarray(x[:, i], dtype=np.float64)
        if np.std(col) <= 1e-8:
            continue
        value = np.corrcoef(col, y)[0, 1]
        out[i] = 0.0 if not np.isfinite(value) else abs(float(value))
    return out


def disjoint_names():
    dropped = set(json.loads((FC01 / "contract.json").read_text())["disjoint_drop"][OX])
    return [name for name in PANEL if name != OX and name not in dropped]


def choose_inputs(train, target):
    if target == "TLS":
        labels = np.concatenate([sec["labels"] for sec in train])
        mask = labeled_tls(labels)
        y = (labels[mask] == "TLS").astype(float)
        names = [name for name in PANEL if name != "B_AXIS"]
        x, _ = stack(train, names, False)
        score = abs_corr(x[mask], y)
        ranked = [names[i] for i in np.argsort(-score) if score[i] > 0]
        return {
            "L1": ["B_AXIS"],
            "L5": ["B_AXIS"] + [name for name in ranked if name != "B_AXIS"][:4],
            "Lfull": list(PANEL),
        }
    names = disjoint_names()
    y = np.concatenate([sec["scores"][:, PANEL.index(OX)] for sec in train]).astype(np.float64)
    x, _ = stack(train, names, False)
    score = abs_corr(x, y)
    ranked = [names[i] for i in np.argsort(-score)]
    return {"L1": ranked[:1], "L5": ranked[:5], "Ldisjoint": names}


def rmse(y, pred):
    return float(np.sqrt(np.mean((pred - y) ** 2)))


def auc(y, pred):
    if len(np.unique(y)) < 2:
        return None
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, pred))


def run():
    HGB["max_iter"] = 120
    OUT.mkdir(parents=True, exist_ok=True)
    sections, _genes = load_sections()
    rows = []
    selections = []
    for unit, train, test in all_folds(sections):
        for target in (OX, "TLS"):
            chosen = choose_inputs(train, target)
            selections.append({"unit": unit, "cohort": train[0]["cohort"], "target": target,
                               "inputs": {k: v for k, v in chosen.items()}})
            for level, names in chosen.items():
                if target == "TLS" and any(name == target for name in names):
                    raise SystemExit("TLS label entered inputs")
                if target == OX and OX in names:
                    raise SystemExit("OXPHOS channel entered inputs")
                for model, spatial in (("hgb_spot", False), ("hgb_spatial", True)):
                    kind = "hgb_cls" if target == "TLS" else "hgb_reg"
                    x_train, _ = stack(train, names, spatial)
                    if target == "TLS":
                        labels = np.concatenate([sec["labels"] for sec in train])
                        mask = labeled_tls(labels)
                        y = (labels[mask] == "TLS").astype(int)
                        weight = class_weight(y)
                        if weight is None:
                            rows.append({"unit": unit, "target": target, "level": level, "model": model,
                                         "status": "NOT_FITTABLE", "cohort": train[0]["cohort"]})
                            continue
                        x_fit = x_train[mask]
                    else:
                        y = np.concatenate([sec["scores"][:, PANEL.index(OX)] for sec in train]).astype(np.float32)
                        weight = None
                        x_fit = x_train
                    tests = [stack([sec], names, spatial)[0] for sec in test]
                    preds, _train_pred = fit_once(kind, x_fit, y, weight, tests)
                    values = []
                    for sec, pred in zip(test, preds):
                        if target == "TLS":
                            mask = labeled_tls(sec["labels"])
                            if mask.sum() == 0 or len(np.unique(sec["labels"][mask])) < 2:
                                values.append(None)
                            else:
                                values.append(auc((sec["labels"][mask] == "TLS").astype(int), pred[mask]))
                        else:
                            values.append(rmse(sec["scores"][:, PANEL.index(OX)], pred))
                    finite = [v for v in values if v is not None]
                    rows.append({
                        "unit": unit,
                        "cohort": test[0]["cohort"],
                        "unit_kind": "patient" if test[0]["cohort"] == "GSE175540" else "specimen",
                        "target": target,
                        "level": level,
                        "model": model,
                        "n_inputs": len(names),
                        "n_sections": len(test),
                        "n_evaluable_sections": len(finite),
                        "status": "OK" if finite else "NOT_TESTABLE",
                        "metric": None if not finite else float(sum(finite) / len(finite)),
                    })
                    print("fit", unit, target, level, model, rows[-1]["status"], flush=True)
    (OUT / "metrics.json").write_text(json.dumps(rows, indent=2) + "\n")
    (OUT / "selections.json").write_text(json.dumps(selections, indent=2) + "\n")
    (OUT / "rule_check.json").write_text(json.dumps(combination_rule(rows), indent=2) + "\n")
    (OUT / "STATUS").write_text("FC02_LADDER_COMPLETE\n")
    print(json.dumps(combination_rule(rows), indent=2))


def combination_rule(rows):
    gse = [r for r in rows if r["cohort"] == "GSE175540" and r["status"] == "OK" and r["model"] == "hgb_spatial"]
    index = {(r["unit"], r["target"], r["level"]): r for r in gse}
    checks = []
    for target, large, small, higher in (
        (OX, "Ldisjoint", "L1", False),
        (OX, "Ldisjoint", "L5", False),
        (OX, "L5", "L1", False),
        ("TLS", "Lfull", "L1", True),
        ("TLS", "Lfull", "L5", True),
        ("TLS", "L5", "L1", True),
    ):
        units = sorted({r["unit"] for r in gse if r["target"] == target})
        deltas, directions = [], []
        for unit in units:
            a = index.get((unit, target, large))
            b = index.get((unit, target, small))
            if not a or not b or a["metric"] is None or b["metric"] is None:
                continue
            if higher:
                delta = a["metric"] - b["metric"]
            else:
                if b["metric"] <= 1e-8:
                    continue
                delta = (b["metric"] - a["metric"]) / b["metric"]
            deltas.append(delta)
            directions.append(delta > 0)
        mean = None if not deltas else float(sum(deltas) / len(deltas))
        fraction = None if not directions else float(sum(directions) / len(directions))
        floor = 0.02
        meets = len(deltas) >= 8 and fraction is not None and fraction >= 0.5 and mean is not None and mean >= floor
        checks.append({"target": target, "large": large, "small": small, "n": len(deltas),
                       "mean_delta": mean, "fraction": fraction, "combination_gain": bool(meets)})
    return {"label": "OPERATIONAL_NOT_MECHANISM", "checks": checks}


if __name__ == "__main__":
    run()
