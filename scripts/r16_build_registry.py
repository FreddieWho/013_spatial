#!/usr/bin/env python3
"""Merge Tier-1 (axis rulers) and Tier-2 (census groups) rows into the R-16
field registry TSV (r16.field_registry.v1). Read-only over both artifacts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

HEADER = ("pattern_id\ttier\tkind\taxis_name\tlineage\tn_sections\tn_patients\t"
          "resolution_support\ttop_markers\tgeometry\teffect_size\tnull_type\t"
          "null_draws\tnull_p\tfdr_q\tcross_patient_shape_corr\tresidual_dimension\t"
          "uncertainty_ci95\tevidence_grade\tnotes")
COLS = HEADER.split("\t")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier1", type=Path, default=ROOT / "infra/r16/registry_rows_tier1_20260914.json")
    ap.add_argument("--tier2", type=Path, default=ROOT / "infra/r16/registry_rows_tier2_20260914.json")
    ap.add_argument("--out", type=Path, default=ROOT / "infra/r16/field_registry.tsv")
    args = ap.parse_args()

    rows = []
    for path in (args.tier1, args.tier2):
        for r in json.load(open(path)):
            rows.append([str(r.get(c, "NA")) for c in COLS])
    n_t1 = sum(1 for r in rows if r[1] == "1")
    n_t2 = len(rows) - n_t1
    args.out.write_text(
        "# r16.field_registry.v1 (2026-09-14, D-115 two-tier; grades capped at "
        "EXPLORATORY_REPRODUCED; CLAIM_CANDIDATE requires separate approval)\n"
        + HEADER + "\n"
        + "\n".join("\t".join(r) for r in rows) + "\n")
    print(f"registry rows: {len(rows)} (tier1={n_t1}, tier2={n_t2}) -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
