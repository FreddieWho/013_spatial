#!/usr/bin/env python3
"""Independent artifact recomputation, not a retraining or human TLS evaluation."""
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "infra/field_discovery_20261018"
OUT = ROOT / "infra/field_discovery_review_20261018"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_table(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def main():
    config = json.loads((WORK / "contract.json").read_text())
    freeze = json.loads((WORK / "freeze_manifest.json").read_text())
    manifest = json.loads((WORK / "input_manifest.json").read_text())
    selected = {r["section_id"]: r for r in manifest if r["status"] == "QUALIFIED"}
    identity_path = ROOT / "infra/tls_four_cancer_labels_20261008/eligibility/sections.tsv"
    identity = {r["section_id"]: {k: r[k] for k in ["patient_status", "patient_key", "patient_scope", "identity_note"]}
                for r in csv.DictReader(identity_path.open(), delimiter="\t")}
    confirmed = {sid: identity[sid]["patient_status"] == "known_within_study"
                 and identity[sid]["patient_key"] == row["unit_id"] for sid, row in selected.items()}
    panels = json.loads((ROOT / config["panel_file"]).read_text())
    definitions = json.loads((WORK / "visible_definitions.json").read_text())
    OUT.mkdir(exist_ok=True)
    checks = Counter()
    assert digest(ROOT / config["split_file"]) == config["split_sha256"]
    assert digest(ROOT / "scripts/field_discovery_20261018.py") == freeze["code_sha256"]
    for name, record in freeze["files"].items():
        p = WORK / name
        assert p.stat().st_size == record["bytes"] and digest(p) == record["sha256"], name
        checks["frozen_file_hashes"] += 1
    for target, views in definitions.items():
        for view, modules in views.items():
            for name, genes in modules.items():
                assert not (set(genes) & set(panels[target])), (target, name)
                checks["target_disjoint_modules"] += 1
    fold_by_unit = defaultdict(set)
    for r in selected.values():
        fold_by_unit[r["unit_id"]].add(r["fold"])
    assert all(len(v) == 1 for v in fold_by_unit.values())
    records = list(csv.DictReader((WORK / "section_metrics.tsv").open(), delimiter="\t"))
    lookup = {(r["section_id"], r["target"], r["view"], r["model"]): r for r in records}
    assert len(lookup) == len(records)
    corrected_records = []
    for record in records:
        corrected_records.append({**record, "patient_known": confirmed[record["section_id"]],
                                  "patient_field_present_v1": record["patient_known"]})
    write_table(OUT / "section_metrics_identity_corrected.tsv", corrected_records)
    write_table(OUT / "identity_review.tsv", [{"section_id": sid, "unit_id": row["unit_id"],
                "patient_field_present_v1": row["patient_known"], "patient_confirmed_within_study": confirmed[sid],
                **identity[sid]} for sid, row in selected.items()])
    unit_squared_errors = defaultdict(list)
    max_metric_difference = 0.0
    for fold in config["training"]["folds"]:
        folder = WORK / f"folds/fold{fold}"
        receipt = json.loads((folder / "group_receipt.json").read_text())
        expected_test = {sid for sid, r in selected.items() if r["fold"] == fold}
        assert set(receipt["test"]) == expected_test
        assert set(receipt["train"]) == set(selected) - expected_test
        train_units = {selected[sid]["unit_id"] for sid in receipt["train"]}
        assert not train_units & {selected[sid]["unit_id"] for sid in receipt["test"]}
        checks["grouped_fold_receipts"] += 1
        ledger = json.loads((folder / "input_receipt.json").read_text())
        assert ledger["human_tls_files_opened"] == []
        assert all(not any(token in p for token in ["annotation/", "incoming/", "tls_labels"]) for p in ledger["inputs"])
        for sid in receipt["test"]:
            row = selected[sid]
            with np.load(WORK / f"cache/{sid}.npz", allow_pickle=False) as z:
                cache_y, barcodes, xy, source_row = [z[k] for k in ["y", "barcodes", "xy_source_pixels", "source_row"]]
            for ti, target in enumerate(config["targets"]):
                for view in config["input_sizes"]:
                    with np.load(folder / f"{sid}_t{ti}_{view}.npz", allow_pickle=False) as z:
                        assert np.array_equal(z["barcodes"], barcodes)
                        assert np.array_equal(z["xy_source_pixels"], xy)
                        assert np.array_equal(z["source_row"], source_row)
                        assert np.array_equal(z["y"], cache_y[:, ti])
                        assert len(z["y"]) == row["n_spots"]
                        checks["complete_aligned_prediction_maps"] += 1
                        y = z["y"].astype(np.float64)
                        for model in z.files:
                            if model in ("y", "barcodes", "xy_source_pixels", "source_row"):
                                continue
                            p = z[model].astype(np.float64)
                            assert p.shape == y.shape and np.all(np.isfinite(p)), (sid, target, model)
                            mse = float(np.mean((y - p) ** 2))
                            centered = float(np.mean(((y - y.mean()) - (p - p.mean())) ** 2))
                            stored = lookup[(sid, target, view, model)]
                            delta = max(abs(mse - float(stored["mse"])), abs(centered - float(stored["centered_mse"])))
                            max_metric_difference = max(max_metric_difference, delta)
                            assert delta < 2e-5, (sid, target, view, model, delta)
                            unit_squared_errors[(target, view, row["unit_id"], model)].append(mse)
                            checks["independent_metric_recomputations"] += 1
    unit_rmse = {k: math.sqrt(np.mean(v)) for k, v in unit_squared_errors.items()}
    unit_confirmed = {unit: all(confirmed[sid] for sid, row in selected.items() if row["unit_id"] == unit)
                      for unit in fold_by_unit}
    identity_strata = []
    for target in config["targets"]:
        for view in config["input_sizes"]:
            for known in [False, True]:
                values = []
                for unit in fold_by_unit:
                    if unit_confirmed[unit] != known:
                        continue
                    baseline = min(unit_rmse[(target, view, unit, model)] for model in config["baselines"])
                    values.append(1 - unit_rmse[(target, view, unit, "space_hgb")] / max(baseline, 1e-12))
                identity_strata.append({"target": target, "view": view, "patient_confirmed_within_study": known,
                                        "n_units": len(values), "macro_relative_rmse_gain": float(np.mean(values)),
                                        "positive_unit_fraction": float(np.mean(np.array(values) > 0))})
    write_table(OUT / "identity_stratified_results.tsv", identity_strata)
    models = sorted({k[3] for k in unit_rmse})
    model_comparison = []
    for target in config["targets"]:
        for view in config["input_sizes"]:
            for model in models:
                values = [v for (t, w, _, m), v in unit_rmse.items() if (t, w, m) == (target, view, model)]
                if values:
                    model_comparison.append({"target": target, "view": view, "model": model,
                                             "n_units": len(values), "macro_unit_rmse": float(np.mean(values))})
    write_table(OUT / "model_comparison.tsv", model_comparison)
    ladder = []
    for target in config["targets"]:
        for model in ["point_hgb", "bag_hgb", "space_hgb"]:
            changes = []
            for unit in fold_by_unit:
                small = unit_rmse[(target, "small5", unit, model)]
                full = unit_rmse[(target, "full15", unit, model)]
                changes.append(1 - full / max(small, 1e-12))
            ladder.append({"target": target, "model": model, "n_units": len(changes),
                           "relative_rmse_gain_full_vs_small": float(np.mean(changes)),
                           "fraction_units_full_better": float(np.mean(np.array(changes) > 0))})
    write_table(OUT / "input_ladder.tsv", ladder)
    result = {"status": "PASS", "checks": dict(checks), "max_metric_absolute_difference": max_metric_difference,
              "n_sections": len(selected), "n_units": len(fold_by_unit),
              "n_known_patient_units": len({r["unit_id"] for sid, r in selected.items() if confirmed[sid]}),
              "raw_patient_string_present_units": len({r["unit_id"] for r in selected.values() if r["patient_known"]}),
              "identity_corrections": [sid for sid, r in selected.items() if r["patient_known"] != confirmed[sid]],
              "identity_metadata_source_sha256": digest(identity_path),
              "identity_note": "v1 patient_known meant a nonempty raw patient field; corrected status uses existing scoped eligibility, never human TLS labels; no training/split/model input changed",
              "n_heldout_spots": sum(r["n_spots"] for r in selected.values()),
              "organs": dict(Counter(r["organ"] for r in selected.values())),
              "audit_code_sha256": digest(Path(__file__)),
              "frozen_run_manifest_sha256": digest(WORK / "freeze_manifest.json"),
              "scope": "stored-artifact recomputation and hash/group/isolation checks, not model retraining or human TLS validation",
              "human_tls_labels_read": False}
    (OUT / "audit.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
