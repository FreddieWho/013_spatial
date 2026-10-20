#!/usr/bin/env python3
"""Compare the finished K=8 probe with the fixed-panel iter-120 results.

Does not refit models. Learned factors win only if they beat the fixed panel
on the same predeclared floor.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROBE = ROOT / "infra/fc03_dictionary_probe_20261006"
FIXED = ROOT / "infra/fc01_recovery_20261006/correction_iter120/metrics_unit.tsv"


def mean(values):
    kept = [v for v in values if v is not None]
    return None if not kept else float(sum(kept) / len(kept))


def load_fixed():
    lines = FIXED.read_text().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    out = {}
    for row in rows:
        if row["cohort"] != "GSE175540" or row["status"] != "OK":
            continue
        key = (row["unit"], row["target"], row["arm"], row["model"])
        if row["target"] == "TLS":
            out[key] = float(row["auc"]) if row["auc"] else None
        elif row["arm"] == "disjoint":
            out[key] = float(row["rmse"]) if row["rmse"] else None
    return out


def main():
    if not (PROBE / "STATUS").exists():
        raise SystemExit("probe not complete")
    learned = json.loads((PROBE / "metrics.json").read_text())
    stable = json.loads((PROBE / "stability.json").read_text())
    fixed = load_fixed()
    checks = []
    for target, model, higher, fixed_arm in (
        ("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "hgb_spatial", False, "disjoint"),
        ("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "hgb_spot", False, "disjoint"),
        ("TLS", "hgb_spatial", True, "full"),
        ("TLS", "hgb_spot", True, "full"),
    ):
        deltas, directions = [], []
        for row in learned:
            if row.get("target") != target or row.get("model") != model or row.get("status") != "OK":
                continue
            base = fixed.get((row["unit"], target, fixed_arm, model))
            if base is None or row.get("metric") is None:
                continue
            if higher:
                delta = row["metric"] - base
            else:
                if base <= 1e-8:
                    continue
                delta = (base - row["metric"]) / base
            deltas.append(delta)
            directions.append(delta > 0)
        mean_delta = mean(deltas)
        fraction = None if not directions else float(sum(directions) / len(directions))
        meets = len(deltas) >= 8 and fraction is not None and fraction >= 0.5 and mean_delta is not None and mean_delta >= 0.02
        checks.append({"target": target, "model": model, "n": len(deltas), "mean_delta": mean_delta,
                       "fraction_learned_better": fraction, "beats_fixed_panel": bool(meets)})
    unstable = [row for row in stable if row.get("unstable")]
    report = {
        "label": "LEARNED_VS_FIXED_NOT_DISCOVERY",
        "checks": checks,
        "n_unstable_fold_target": len(unstable),
        "n_stability_rows": len(stable),
        "claim": "no learned-field discovery" if not any(c["beats_fixed_panel"] and c["model"] == "hgb_spatial" for c in checks) else "operational_gain_only_if_stable",
    }
    if any(c["beats_fixed_panel"] and c["model"] == "hgb_spatial" for c in checks) and len(unstable) > len(stable) / 2:
        report["claim"] = "partial_candidate_only"
    (PROBE / "fixed_comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
