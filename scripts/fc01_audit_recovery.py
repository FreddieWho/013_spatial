#!/usr/bin/env python3
"""Post-run leakage audit for the frozen FC-01 comparison.

Checks identities, gene overlap and output completeness. Does not refit models
or read probe scores as evaluation.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "infra/fc01_recovery_20261006"
CONTRACT = json.loads((WORK / "contract.json").read_text())


def main():
    run = WORK / "run"
    if not (run / "STATUS").exists():
        raise SystemExit("official run is not complete")
    if (run / "STATUS").read_text().strip() != "OFFICIAL_FOLDS_COMPLETE":
        raise SystemExit("official status is not complete")
    genes = json.loads((ROOT / CONTRACT["panel_path"]).read_text())
    identity = (ROOT / CONTRACT["identity_path"]).read_text().splitlines()
    header = identity[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in identity[1:] if line.strip()]
    expected = {}
    for row in rows:
        unit = row["patient_id"] if row["cohort"] == "GSE175540" else row["section_id"]
        expected.setdefault((row["cohort"], unit), []).append(row["section_id"])
    problems = []
    seen = []
    for cohort, unit in sorted(expected):
        fold = run / "folds" / unit.replace(":", "_")
        if not (fold / "DONE").exists():
            problems.append(f"missing fold {unit}")
            continue
        metrics = json.loads((fold / "metrics.json").read_text())
        sections = sorted({row["section_id"] for row in metrics if "section_id" in row})
        if sections != sorted(expected[(cohort, unit)]):
            problems.append(f"{unit} sections {sections} != {sorted(expected[(cohort, unit)])}")
        if unit == "GSE175540:c_2" and len(sections) != 2:
            problems.append("c_2 was not held out together")
        for row in metrics:
            if row.get("target") in row.get("model", ""):
                pass
            if row.get("target") != "TLS" and row.get("arm") == "disjoint":
                banned = set(genes[row["target"]])
                for name in CONTRACT["disjoint_drop"][row["target"]]:
                    if banned & set(genes[name]) and name not in CONTRACT["disjoint_drop"][row["target"]]:
                        problems.append(f"drop list inconsistent for {name}")
        for sid in sections:
            npz = np.load(fold / f"{sid}.npz")
            for key in npz.files:
                if key in {"labels", "xy"}:
                    continue
                target = key.split("__", 1)[0]
                if target != "TLS" and target in key and npz[key].shape[0] != len(npz["labels"]):
                    problems.append(f"{sid} {key} length mismatch")
                values = np.asarray(npz[key])
                if values.dtype.kind == "f" and not np.isfinite(values).all():
                    problems.append(f"{sid} {key} has non-finite predictions")
        seen.append(unit)
    # Disjoint lists must actually remove every sharing program.
    for target, dropped in CONTRACT["disjoint_drop"].items():
        banned = set(genes[target])
        for name, members in genes.items():
            if name == target:
                continue
            shares = bool(banned & set(members))
            if shares and name not in dropped:
                problems.append(f"{target} still shares genes with visible {name}")
            if not shares and name in dropped:
                problems.append(f"{target} dropped non-overlapping {name}")
    report = {
        "status": "PASS" if not problems else "FAIL",
        "n_units": len(seen),
        "n_expected_units": len(expected),
        "problems": problems,
        "probe_not_used": True,
    }
    (run / "leakage_audit.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if problems:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
