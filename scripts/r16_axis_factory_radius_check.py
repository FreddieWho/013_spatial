#!/usr/bin/env python3
"""Radius check for promoted Tier-1 v2 groups (design v2.1, D-118).

Recomputes candidate loadings at k=30 for the sections feeding promoted
groups and reports per-group agreement (median best-match |cos| between
k=15 and k=30 loadings of the same section). Escalation rule: median < 0.9
-> report only, do not silently proceed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from r16 import axis_factory as F
from r16 import census as C

ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=ROOT / "infra/r16/census_cache_hvg10k/sections")
    ap.add_argument("--artifact", type=Path, default=ROOT / "infra/r16/axis_factory_20260915.json")
    ap.add_argument("--loadings", type=Path,
                    default=ROOT / "infra/r16/axis_factory_loadings_20260915.npz")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "infra/r16/axis_factory_radius_check_20260915.json")
    args = ap.parse_args()

    d = json.load(open(args.artifact))
    z = np.load(args.loadings)
    keys = [str(k) for k in z["keys"]]
    L15 = z["loadings"].astype(np.float64)
    key_index = {k: i for i, k in enumerate(keys)}
    promoted = [g for g in d["groups"] if g["evidence_grade"] == "EXPLORATORY_REPRODUCED"]
    stems = sorted({c["stem"] for c in d["candidates"]})

    # recompute k30 loadings per section (PCA only, no LISA)
    k30 = {}
    for stem in stems:
        mat, barcodes, genes, coords, meta = C.load_section(args.cache, stem)
        x_log = np.log1p(mat.toarray().astype(np.float32))
        x_z = C.log1p_zscore(mat)
        m30 = F.zscore_columns(np.asarray((F.spatial_mean_weights(coords, 30) @ x_log)))
        per_lam = {}
        for lam in F.LAMBDAS:
            _, loadings, _ = F.candidate_pcs(
                F.augment(x_z, m30.astype(np.float32), lam), F.N_PC)
            per_lam[lam] = loadings
        k30[stem] = per_lam

    report = {"groups": {}, "overall": {}}
    all_best = []
    for g in promoted:
        members = [c for c in d["candidates"]
                   if c["patient_id"] in set(g["patients"])]
        # restrict to stems actually contributing: use group member stems
        bests = []
        seen_stems = sorted({c["stem"] for c in members})
        for stem in seen_stems:
            li = [key_index[k] for k in keys if k.startswith(stem + "::")]
            for i in li:
                lam = d["candidates"][i]["lambda"]
                pc = d["candidates"][i]["pc"]
                a = L15[i] / np.linalg.norm(L15[i])
                B = k30[stem][lam]
                bn = B / np.linalg.norm(B, axis=0, keepdims=True)
                bests.append(float(np.abs(a @ bn).max()))
        bests = np.array(bests)
        all_best.extend(bests.tolist())
        report["groups"][g["pattern_id"]] = {
            "n_compared": int(len(bests)),
            "median_best_match_cos": float(np.median(bests)),
            "frac_above_0.9": float((bests >= 0.9).mean()),
        }
    all_best = np.array(all_best)
    report["overall"] = {
        "median_best_match_cos": float(np.median(all_best)),
        "frac_above_0.9": float((all_best >= 0.9).mean()),
        "n_compared": int(len(all_best)),
        "escalation": bool(float(np.median(all_best)) < 0.9),
    }
    args.out.write_text(json.dumps(report, indent=2, ensure_ascii=False))
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
