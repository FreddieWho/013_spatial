#!/usr/bin/env python3
"""FC-01 held-out field recovery under the frozen 2026-10-06 contract.

Reads only cached per-spot scores, coordinates and labels. Does not read raw
counts, label-derived distances, endpoints, or full-graph latents.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy import sparse
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "infra/fc01_recovery_20261006"
CONTRACT = json.loads((WORK / "contract.json").read_text())
PANEL = CONTRACT["panel"]
RINGS = [(float(a), float(b)) for a, b in CONTRACT["rings_um"]]
HGB = CONTRACT["models"]["hgb_spot"]


def load_sections():
    identity = (ROOT / CONTRACT["identity_path"]).read_text().splitlines()
    header = identity[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in identity[1:] if line.strip()]
    genes = json.loads((ROOT / CONTRACT["panel_path"]).read_text())
    sections = []
    for row in rows:
        sid = row["section_id"]
        score = np.load(ROOT / CONTRACT["score_dir"] / f"{sid}.npz")
        names = [str(x) for x in score["set_ids"]]
        index = {name: i for i, name in enumerate(names)}
        missing = [name for name in PANEL if name not in index]
        if missing:
            raise SystemExit(f"{sid} missing panel columns: {missing[:3]}")
        values = np.asarray(score["scores"][:, [index[name] for name in PANEL]], dtype=np.float32)
        geo = np.load(ROOT / CONTRACT["geometry_dir"] / f"{sid}.npz")
        xy = np.asarray(geo["xy"], dtype=np.float64)
        labels = np.array([str(x) for x in geo["labels"]])
        if not (len(xy) == len(labels) == len(values)):
            raise SystemExit(f"{sid} row alignment failed")
        unit = row["patient_id"] if row["cohort"] == "GSE175540" else sid
        if row["cohort"] == "GSE175540" and not unit:
            raise SystemExit(f"{sid} missing confirmed patient id")
        if row["cohort"] == "USZ" and row["patient_id"]:
            raise SystemExit(f"{sid} unexpectedly has a patient id")
        sections.append({
            "section_id": sid,
            "cohort": row["cohort"],
            "unit": unit,
            "xy": xy,
            "labels": labels,
            "scores": values,
            "rings": ring_features(xy, values),
        })
        print("loaded", sid, values.shape[0], flush=True)
    return sections, genes


def ring_features(xy, scores):
    tree = cKDTree(xy)
    coo = tree.sparse_distance_matrix(tree, max_distance=RINGS[-1][1], output_type="coo_matrix")
    keep = coo.row != coo.col
    row, col, dist = coo.row[keep], coo.col[keep], coo.data[keep]
    means, empty = [], []
    n = len(xy)
    for lo, hi in RINGS:
        mask = (dist > lo) & (dist <= hi)
        adj = sparse.csr_matrix((np.ones(mask.sum(), dtype=np.float32), (row[mask], col[mask])), shape=(n, n))
        counts = np.asarray(adj.sum(axis=1)).ravel()
        mu = adj @ scores
        positive = counts > 0
        mu[positive] /= counts[positive, None]
        means.append(mu.astype(np.float32))
        empty.append((~positive).astype(np.float32))
    return means, empty


def visible_columns(target, arm):
    drop = {target} if target != "TLS" else set()
    if arm == "disjoint":
        drop.update(CONTRACT["disjoint_drop"][target])
    names = [name for name in PANEL if name not in drop]
    if target in names:
        raise SystemExit("target channel remained visible")
    return names


def assert_gene_disjoint(names, target, genes):
    if target == "TLS":
        return 0
    banned = set(genes[target])
    overlap = 0
    for name in names:
        shared = banned & set(genes[name])
        overlap += len(shared)
        if shared:
            raise SystemExit(f"disjoint arm leaked {name}: {sorted(shared)[:5]}")
    return overlap


def stack(sections, columns, spatial):
    col = [PANEL.index(name) for name in columns]
    blocks, rings = [], []
    for sec in sections:
        spot = sec["scores"][:, col]
        if spatial:
            means, empty = sec["rings"]
            extra = [mu[:, col] for mu in means] + [flag[:, None] for flag in empty]
            blocks.append(np.hstack([spot, *extra]))
        else:
            blocks.append(spot)
        rings.append(len(sec["labels"]))
    return np.vstack(blocks), rings


def labeled_tls(labels):
    return np.isin(labels, ["TLS", "NO_TLS"])


def fit_predict(kind, x_train, y_train, x_test, weight):
    if kind == "prior":
        value = float(np.average(y_train, weights=weight)) if weight is not None else float(y_train.mean())
        return np.full(len(x_test), value, dtype=np.float32), np.full(len(y_train), value, dtype=np.float32)
    if kind == "ridge":
        scaler = StandardScaler()
        model = Ridge(alpha=1.0)
        model.fit(scaler.fit_transform(x_train), y_train)
        return model.predict(scaler.transform(x_test)).astype(np.float32), model.predict(scaler.transform(x_train)).astype(np.float32)
    if kind == "logistic":
        scaler = StandardScaler()
        model = LogisticRegression(C=1.0, solver="lbfgs", max_iter=500)
        model.fit(scaler.fit_transform(x_train), y_train, sample_weight=weight)
        train_x = scaler.transform(x_train)
        test_x = scaler.transform(x_test)
        return model.predict_proba(test_x)[:, 1].astype(np.float32), model.predict_proba(train_x)[:, 1].astype(np.float32)
    params = dict(max_depth=HGB["max_depth"], max_iter=HGB["max_iter"], learning_rate=HGB["learning_rate"],
                  min_samples_leaf=HGB["min_samples_leaf"], random_state=HGB["random_state"], early_stopping=False)
    if kind == "hgb_cls":
        model = HistGradientBoostingClassifier(**params)
        model.fit(x_train, y_train, sample_weight=weight)
        return model.predict_proba(x_test)[:, 1].astype(np.float32), model.predict_proba(x_train)[:, 1].astype(np.float32)
    model = HistGradientBoostingRegressor(**params)
    model.fit(x_train, y_train)
    return model.predict(x_test).astype(np.float32), model.predict(x_train).astype(np.float32)


def class_weight(y):
    y = y.astype(int)
    out = np.ones(len(y), dtype=np.float64)
    for value in (0, 1):
        count = int((y == value).sum())
        if count == 0:
            return None
        out[y == value] = 0.5 * len(y) / count
    return out


def safe_auc(y, pred):
    if len(np.unique(y)) < 2:
        return None
    return float(roc_auc_score(y, pred))


def safe_pearson(y, pred):
    y = np.asarray(y, dtype=np.float64)
    pred = np.asarray(pred, dtype=np.float64)
    if len(y) < 3 or not np.isfinite(y).all() or not np.isfinite(pred).all():
        return None
    if np.std(y) <= 1e-8 or np.std(pred) <= 1e-8:
        return None
    with np.errstate(divide="ignore", invalid="ignore"):
        value = np.corrcoef(y, pred)[0, 1]
    return None if not np.isfinite(value) else float(value)


def components(xy, labels):
    points = np.flatnonzero(labels == "TLS")
    if len(points) == 0:
        return []
    parent = {int(i): int(i) for i in points}

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    if len(points) > 1:
        pairs = cKDTree(xy[points]).query_pairs(130.01)
        for a, b in pairs:
            ia, ib = int(points[a]), int(points[b])
            ra, rb = find(ia), find(ib)
            if ra != rb:
                parent[rb] = ra
    groups = {}
    for i in points:
        groups.setdefault(find(int(i)), []).append(int(i))
    return [np.array(v) for v in groups.values() if len(v) >= 2]


def structure_summary(xy, labels, pred, threshold):
    found = components(xy, labels)
    if not found:
        return {"n_components": 0, "median_component_pred": None, "median_hotspot_distance": None, "hotspot_status": "NO_COMPONENT"}
    if threshold is None or not np.isfinite(threshold):
        status = "DEGENERATE_HOTSPOT"
        distances = [None] * len(found)
    else:
        hot = np.flatnonzero(pred >= threshold)
        status = "OK" if len(hot) else "NO_HOTSPOT"
        distances = []
        for comp in found:
            centroid = xy[comp].mean(axis=0)
            distances.append(None if len(hot) == 0 else float(np.linalg.norm(xy[hot] - centroid, axis=1).min()))
    finite = [d for d in distances if d is not None]
    return {
        "n_components": len(found),
        "median_component_pred": float(np.median([pred[comp].mean() for comp in found])),
        "median_hotspot_distance": None if not finite else float(np.median(finite)),
        "hotspot_status": status,
    }


def evaluate(target, labels, pred, xy, train_pred, train_labels):
    row = {"status": "OK", "n_eval": int(len(pred)), "auc": None, "ap": None, "brier": None,
           "rmse": None, "pearson": None, "n_components": None, "median_component_pred": None,
           "median_hotspot_distance": None, "hotspot_status": None}
    if target == "TLS":
        mask = labeled_tls(labels)
        row["n_eval"] = int(mask.sum())
        if mask.sum() == 0 or len(np.unique(labels[mask])) < 2:
            row["status"] = "NOT_TESTABLE"
            return row
        y = (labels[mask] == "TLS").astype(int)
        p = pred[mask]
        row["auc"] = safe_auc(y, p)
        row["ap"] = float(average_precision_score(y, p))
        row["brier"] = float(brier_score_loss(y, np.clip(p, 0, 1)))
        train_mask = labeled_tls(train_labels)
        threshold = None if train_mask.sum() == 0 else float(np.quantile(train_pred[train_mask], 0.9))
        if threshold is not None and np.nanstd(train_pred[train_mask]) <= 1e-8:
            threshold = None
        row.update(structure_summary(xy, labels, pred, threshold))
        return row
    y = pred * 0 + labels.astype(float) if False else None
    return row


def molecular_metrics(y, pred):
    row = {"status": "OK", "n_eval": int(len(y)), "auc": None, "ap": None, "brier": None,
           "rmse": float(np.sqrt(np.mean((pred - y) ** 2))), "pearson": safe_pearson(y, pred),
           "n_components": None, "median_component_pred": None, "median_hotspot_distance": None,
           "hotspot_status": None}
    if not np.isfinite(row["rmse"]):
        row["status"] = "FAILED_NONFINITE"
    return row


def models_for(target):
    out = [("prior", "prior", False), ("ridge_spot", "logistic" if target == "TLS" else "ridge", False),
           ("hgb_spot", "hgb_cls" if target == "TLS" else "hgb_reg", False),
           ("hgb_spatial", "hgb_cls" if target == "TLS" else "hgb_reg", True)]
    if target == "TLS":
        out.append(("b_axis", "logistic", False))
    return out


def run_fold(train, test, genes, out_dir):
    assert {sec["unit"] for sec in train}.isdisjoint({sec["unit"] for sec in test})
    if test[0]["unit"] == "GSE175540:c_2":
        assert len(test) == 2
    receipts = []
    for target_info in CONTRACT["targets"]:
        target = target_info["id"]
        arms = ["full"] + (["disjoint"] if target != "TLS" else [])
        for arm in arms:
            columns = visible_columns(target, arm)
            if arm == "disjoint":
                assert_gene_disjoint(columns, target, genes)
            if target == "TLS" and "B_AXIS" not in columns:
                raise SystemExit("TLS full input lost B-axis")
            for model_name, kind, spatial in models_for(target):
                used = ["B_AXIS"] if model_name == "b_axis" else columns
                if model_name == "b_axis" and target != "TLS":
                    continue
                t0 = time.time()
                x_train, _ = stack(train, used, spatial)
                y_labels = np.concatenate([sec["labels"] for sec in train])
                y_scores = np.concatenate([sec["scores"][:, PANEL.index(target)] for sec in train]) if target != "TLS" else None
                if target == "TLS":
                    mask = labeled_tls(y_labels)
                    y = (y_labels[mask] == "TLS").astype(int)
                    x_fit = x_train[mask]
                    weight = class_weight(y)
                    if weight is None:
                        receipts.append({"target": target, "arm": arm, "model": model_name, "status": "NOT_FITTABLE"})
                        continue
                else:
                    y = y_scores.astype(np.float32)
                    x_fit = x_train
                    weight = None
                preds = []
                train_pred = None
                for sec in test:
                    x_test, _ = stack([sec], used, spatial)
                    pred, train_pred = fit_predict(kind, x_fit, y, x_test, weight)
                    preds.append(pred)
                # The loop refits conceptually once; fit_predict fits every test section.
                # Refit is deterministic, but avoid repeating it.
                raise RuntimeError("replaced below")
    return receipts


def fit_once(kind, x_fit, y, weight, tests):
    if kind == "prior":
        value = float(np.average(y, weights=weight)) if weight is not None else float(np.mean(y))
        train_pred = np.full(len(y), value, dtype=np.float32)
        return [np.full(len(x), value, dtype=np.float32) for x in tests], train_pred
    if kind == "ridge":
        scaler = StandardScaler()
        model = Ridge(alpha=1.0)
        model.fit(scaler.fit_transform(x_fit), y)
        train_pred = model.predict(scaler.transform(x_fit)).astype(np.float32)
        return [model.predict(scaler.transform(x)).astype(np.float32) for x in tests], train_pred
    if kind == "logistic":
        scaler = StandardScaler()
        model = LogisticRegression(C=1.0, solver="lbfgs", max_iter=500)
        model.fit(scaler.fit_transform(x_fit), y, sample_weight=weight)
        train_x = scaler.transform(x_fit)
        return [model.predict_proba(scaler.transform(x))[:, 1].astype(np.float32) for x in tests], model.predict_proba(train_x)[:, 1].astype(np.float32)
    params = dict(max_depth=HGB["max_depth"], max_iter=HGB["max_iter"], learning_rate=HGB["learning_rate"],
                  min_samples_leaf=HGB["min_samples_leaf"], random_state=HGB["random_state"], early_stopping=False)
    model = HistGradientBoostingClassifier(**params) if kind == "hgb_cls" else HistGradientBoostingRegressor(**params)
    model.fit(x_fit, y, sample_weight=weight) if kind == "hgb_cls" else model.fit(x_fit, y)
    train_pred = model.predict_proba(x_fit)[:, 1].astype(np.float32) if kind == "hgb_cls" else model.predict(x_fit).astype(np.float32)
    test_pred = [model.predict_proba(x)[:, 1].astype(np.float32) if kind == "hgb_cls" else model.predict(x).astype(np.float32) for x in tests]
    return test_pred, train_pred


def execute(train, test, genes):
    rows = []
    arrays = {}
    for target_info in CONTRACT["targets"]:
        target = target_info["id"]
        for arm in (["full"] if target == "TLS" else ["full", "disjoint"]):
            columns = visible_columns(target, arm)
            overlap = 0 if arm == "disjoint" else None
            if arm == "disjoint":
                assert_gene_disjoint(columns, target, genes)
            for model_name, kind, spatial in models_for(target):
                used = ["B_AXIS"] if model_name == "b_axis" else columns
                if target != "TLS" and model_name == "b_axis":
                    continue
                if target in used:
                    raise SystemExit("target visible")
                t0 = time.time()
                x_train, _ = stack(train, used, spatial)
                labels = np.concatenate([sec["labels"] for sec in train])
                if target == "TLS":
                    mask = labeled_tls(labels)
                    y = (labels[mask] == "TLS").astype(int)
                    x_fit, weight = x_train[mask], class_weight(y)
                    train_labels_for_threshold = labels[mask]
                    if weight is None:
                        rows.append(base_row(test, target, arm, model_name, "NOT_FITTABLE", len(used)))
                        continue
                else:
                    y = np.concatenate([sec["scores"][:, PANEL.index(target)] for sec in train]).astype(np.float32)
                    x_fit, weight = x_train, None
                    train_labels_for_threshold = None
                x_tests = [stack([sec], used, spatial)[0] for sec in test]
                # Class weights balance classifiers. The prior must stay the unweighted training mean.
                fit_weight = None if kind == "prior" else weight
                pred_tests, train_pred = fit_once(kind, x_fit, y, fit_weight, x_tests)
                for sec, pred in zip(test, pred_tests):
                    key = f"{target}__{arm}__{model_name}"
                    arrays.setdefault(sec["section_id"], {})[key] = pred
                    if target == "TLS":
                        metric = evaluate(target, sec["labels"], pred, sec["xy"], train_pred, train_labels_for_threshold)
                    else:
                        metric = molecular_metrics(sec["scores"][:, PANEL.index(target)], pred)
                    metric.update(base_row([sec], target, arm, model_name, metric["status"], len(used)))
                    metric["seconds"] = round(time.time() - t0, 3)
                    metric["gene_overlap_pairs"] = overlap
                    rows.append(metric)
                print("fit", test[0]["unit"], target, arm, model_name, round(time.time() - t0, 1), flush=True)
    for sec in test:
        path = None
        payload = arrays[sec["section_id"]]
        payload["labels"] = sec["labels"]
        payload["xy"] = sec["xy"]
        sec["_pred"] = payload
    return rows


def base_row(test, target, arm, model, status, n_features):
    sec = test[0]
    return {"section_id": sec["section_id"], "cohort": sec["cohort"], "unit": sec["unit"],
            "unit_kind": "patient" if sec["cohort"] == "GSE175540" else "specimen",
            "target": target, "arm": arm, "model": model, "status": status, "n_features": n_features}


def write_fold(out_dir, unit, test, rows):
    fold = out_dir / "folds" / unit.replace(":", "_")
    fold.mkdir(parents=True, exist_ok=True)
    (fold / "metrics.json").write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n")
    for sec in test:
        np.savez_compressed(fold / f"{sec['section_id']}.npz", **sec["_pred"])
    (fold / "DONE").write_text(unit + "\n")


def all_folds(sections):
    folds = []
    for cohort in ("GSE175540", "USZ"):
        units = sorted({sec["unit"] for sec in sections if sec["cohort"] == cohort})
        for unit in units:
            test = [sec for sec in sections if sec["cohort"] == cohort and sec["unit"] == unit]
            train = [sec for sec in sections if sec["cohort"] == cohort and sec["unit"] != unit]
            folds.append((unit, train, test))
    return folds


def execute_selected(sections, genes, out_dir, only_unit=None):
    out_dir.mkdir(parents=True, exist_ok=True)
    free = __import__("shutil").disk_usage(ROOT).free
    if free < 1_200_000_000_000:
        raise SystemExit(f"storage floor breached: {free}")
    done = []
    for unit, train, test in all_folds(sections):
        if only_unit and unit != only_unit:
            continue
        marker = out_dir / "folds" / unit.replace(":", "_") / "DONE"
        if marker.exists():
            print("skip", unit, flush=True)
            continue
        print("fold", unit, "train", len(train), "test", len(test), flush=True)
        rows = execute(train, test, genes)
        write_fold(out_dir, unit, test, rows)
        done.append(unit)
    return done


def mean_defined(values):
    kept = [float(v) for v in values if v is not None]
    return None if not kept else float(sum(kept) / len(kept))


def write_tsv(path, rows):
    if not rows:
        path.write_text("")
        return
    keys = list(rows[0])
    lines = ["\t".join(keys)]
    for row in rows:
        lines.append("\t".join("" if row.get(k) is None else str(row.get(k)) for k in keys))
    path.write_text("\n".join(lines) + "\n")


def summarize(out_dir):
    """Aggregate frozen fold metrics. Does not refit or change the contract."""
    paths = sorted((out_dir / "folds").glob("*/metrics.json"))
    if not paths:
        raise SystemExit(f"no fold metrics in {out_dir}")
    section_rows = []
    for path in paths:
        section_rows.extend(json.loads(path.read_text()))
    section_rows = [row for row in section_rows if "section_id" in row]
    grouped = {}
    for row in section_rows:
        key = (row["cohort"], row["unit"], row["unit_kind"], row["target"], row["arm"], row["model"])
        grouped.setdefault(key, []).append(row)
    units = []
    for key, items in grouped.items():
        cohort, unit, kind, target, arm, model = key
        ok = [item for item in items if item.get("status") == "OK"]
        units.append({
            "cohort": cohort, "unit": unit, "unit_kind": kind, "target": target, "arm": arm, "model": model,
            "n_sections": len(items), "n_evaluable_sections": len(ok),
            "status": "OK" if ok else items[0].get("status", "NOT_TESTABLE"),
            "auc": mean_defined(item.get("auc") for item in ok),
            "ap": mean_defined(item.get("ap") for item in ok),
            "brier": mean_defined(item.get("brier") for item in ok),
            "rmse": mean_defined(item.get("rmse") for item in ok),
            "pearson": mean_defined(item.get("pearson") for item in ok),
        })
    write_tsv(out_dir / "metrics_section.tsv", section_rows)
    write_tsv(out_dir / "metrics_unit.tsv", units)
    rule = operational_rule(units)
    (out_dir / "rule_check.json").write_text(json.dumps(rule, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(rule, indent=2, ensure_ascii=False), flush=True)
    return rule


def operational_rule(units):
    """Apply the pre-frozen continuation arithmetic. Not a scientific proof."""
    gse = [row for row in units if row["cohort"] == "GSE175540" and row["status"] == "OK"]
    index = {(row["unit"], row["target"], row["arm"], row["model"]): row for row in gse}
    checks = []
    for target, arm, metric, floor, higher_better in (
        ("TLS", "full", "auc", CONTRACT["gain_rule"]["tls_min_mean_auc_gain"], True),
        ("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "disjoint", "rmse", CONTRACT["gain_rule"]["molecular_min_relative_rmse_reduction"], False),
    ):
        units_present = sorted({key[0] for key in index if key[1] == target and key[2] == arm})
        deltas, directions = [], []
        for unit in units_present:
            spatial = index.get((unit, target, arm, "hgb_spatial"))
            spot = index.get((unit, target, arm, "hgb_spot"))
            if not spatial or not spot or spatial[metric] is None or spot[metric] is None:
                continue
            if higher_better:
                delta = spatial[metric] - spot[metric]
                directions.append(delta > 0)
            else:
                base = spot[metric]
                delta = None if base is None or base <= 1e-8 else (base - spatial[metric]) / base
                if delta is None:
                    continue
                directions.append(delta > 0)
            deltas.append(delta)
        mean_delta = mean_defined(deltas)
        fraction = None if not directions else float(sum(directions) / len(directions))
        meets = (
            len(deltas) >= CONTRACT["gain_rule"]["min_evaluable_patients"]
            and fraction is not None and fraction >= CONTRACT["gain_rule"]["min_patient_fraction_same_direction"]
            and mean_delta is not None and mean_delta >= floor
            and (mean_delta > 0)
        )
        checks.append({"target": target, "arm": arm, "n_evaluable_patients": len(deltas), "mean_delta": mean_delta,
                       "fraction_spatial_better": fraction, "meets_operational_rule": bool(meets),
                       "opens_fc02": bool(meets)})
    return {"label": "OPERATIONAL_NOT_SCIENTIFIC_PROOF", "checks": checks,
            "open_fc02": any(item["opens_fc02"] for item in checks)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["probe", "run", "summarize", "correct"])
    parser.add_argument("--holdout", default="GSE175540:c_5")
    parser.add_argument("--input", default="")
    args = parser.parse_args()
    if args.command == "summarize":
        summarize(Path(args.input) if args.input else WORK / "run")
        return
    if args.command == "correct":
        # Pre-specified single correction: capacity only. Does not rewrite the frozen contract.
        HGB["max_iter"] = 120
        out = WORK / "correction_iter120"
        (out).mkdir(parents=True, exist_ok=True)
        (out / "CORRECTION.txt").write_text(
            "Single pre-specified correction: HistGradientBoosting max_iter 50 -> 120. "
            "Targets, splits, rings, inputs and gain rule unchanged. Not a new search.\n")
        sections, genes = load_sections()
        execute_selected(sections, genes, out)
        (out / "STATUS").write_text("CORRECTION_ITER120_COMPLETE\n")
        summarize(out)
        return
    if __import__("os").environ.get("CUDA_VISIBLE_DEVICES", "") not in ("", "-1"):
        raise SystemExit("GPU is disabled; unset CUDA_VISIBLE_DEVICES or set it to -1")
    sections, genes = load_sections()
    if args.command == "probe":
        out = WORK / "probe_not_evaluation"
        if out.exists():
            import shutil
            shutil.rmtree(out)
        execute_selected(sections, genes, out, args.holdout)
        (out / "STATUS").write_text("TIMING_PROBE_NOT_EVALUATION\n")
    else:
        execute_selected(sections, genes, WORK / "run")
        (WORK / "run" / "STATUS").write_text("OFFICIAL_FOLDS_COMPLETE\n")
        summarize(WORK / "run")


if __name__ == "__main__":
    main()
