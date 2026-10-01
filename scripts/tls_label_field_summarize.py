#!/usr/bin/env python3
"""Summarise the v8.1 label-anchor field test, including the pre-declared
sensitivities (drop weak-label USZ sections; GSE-only)."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/tls_label_field_20261001"
WEAK = {"usz-LC4", "usz-LC5"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    null_rows = list(csv.DictReader(open(OUT / f"label_field_null{args.tag}.tsv"), delimiter="\t"))
    sec_rows = list(csv.DictReader(open(OUT / f"label_field_sections{args.tag}.tsv"), delimiter="\t"))

    def run(rows, label):
        ps = np.array([float(r["p_one_sided"]) for r in rows if r["p_one_sided"] != ""])
        obs = np.array([float(r["obs"]) for r in rows if r["obs"] != ""])
        nl = np.array([float(r["null_median"]) for r in rows if r["null_median"] != ""])
        print(f"[{label}] section-tests with p: {len(ps)}  frac p<=0.05 = "
              f"{(ps <= 0.05).mean():.3f}  median p {np.median(ps):.3f}")
        per_set = {}
        for set_id in sorted({r["set_id"] for r in rows}):
            sub = [r for r in rows if r["set_id"] == set_id]
            p = np.array([float(r["p_one_sided"]) for r in sub if r["p_one_sided"] != ""])
            o = np.array([float(r["obs"]) for r in sub if r["obs"] != ""])
            n = np.array([float(r["null_median"]) for r in sub if r["null_median"] != ""])
            frac = float((p <= 0.05).mean()) if len(p) else float("nan")
            claim = int(len(p) >= 8 and frac >= 0.5 and len(o) and len(n)
                        and np.median(o) > np.median(n))
            per_set[set_id] = {"n_eval": len(p), "frac_p05": frac, "claim": claim,
                               "med_obs": float(np.median(o)) if len(o) else float("nan")}
        return per_set

    all_ps = run(null_rows, "primary: all 26 sections")
    ctrl = {k: v for k, v in all_ps.items() if k.startswith("CTRL_") or k.startswith("GENE_")}
    print("\ncontrols / single genes:")
    for k in sorted(ctrl):
        v = ctrl[k]
        print(f"  {k:28} n_eval={v['n_eval']:2} frac_p05={v['frac_p05']:.2f} "
              f"med_obs={v['med_obs']:+.5f} claim={v['claim']}")
    top = sorted(all_ps.items(), key=lambda kv: -kv[1]["frac_p05"])[:12]
    print("\ntop by frac(p<=0.05):")
    for k, v in top:
        print(f"  {k[:56]:56} n_eval={v['n_eval']:2} frac_p05={v['frac_p05']:.2f} "
              f"med_obs={v['med_obs']:+.5f} claim={v['claim']}")
    clam = [k for k, v in all_ps.items() if v["claim"]]
    print(f"\nreadouts meeting the pre-declared claim rule: {clam if clam else 'NONE'}")

    # per-cohort
    coh = {}
    for r in null_rows:
        if r["p_one_sided"] == "":
            continue
        coh.setdefault(r["cohort"], []).append(float(r["p_one_sided"]))
    print("\nper cohort (all readouts):")
    for c, v in sorted(coh.items()):
        v = np.array(v)
        print(f"  {c:12} n={len(v):4} frac p<=0.05 = {(v <= 0.05).mean():.3f}")

    # pre-declared sensitivities
    print()
    run([r for r in null_rows if r["section_id"] not in WEAK], "sens A: drop weak-label LC4/LC5")
    run([r for r in null_rows if r["cohort"] == "GSE175540"], "sens B: GSE175540 only")
    run([r for r in null_rows if r["cohort"] == "USZ"], "sens C: USZ only")

    # per-set own_contrast sign consistency from the sections table
    print("\nsign consistency of the conditioned contrast (own_contrast):")
    for set_id in sorted({r["set_id"] for r in sec_rows}):
        v = np.array([float(r["own_contrast"]) for r in sec_rows if r["set_id"] == set_id])
        if len(v):
            print(f"  {set_id[:56]:56} pos_frac={float((v > 0).mean()):.2f} med={np.median(v):+.5f}")


if __name__ == "__main__":
    main()
