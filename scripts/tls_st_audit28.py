#!/usr/bin/env python3
"""Full 28-section audit of the ST-detection probe (D-218).

GSE (20 sections): sensitivity/specificity vs TLS/NO_TLS labels.
USZ (8 sections): TLS-spot coverage + false-positive rate on explicit
non-TLS (TUM+NOR); UNASSIGNED/blank never scored.
Pass bar (frozen): median GSE sens >= 0.5 AND median GSE spec >= 0.9.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.tls_st_detect_probe import run as probe_run, tls_score, smooth, components  # noqa: E402
from scripts.tls_peripheral import source_metadata, counts  # noqa: E402
from sklearn.cluster import KMeans

OUT = ROOT / "infra/tls_st_audit28_20261006"


def usz_audit(sid):
    m = source_metadata("USZ", sid)
    X, genes = counts("USZ", sid, m)
    score = tls_score(X, np.array(genes))
    sm = smooth(score, m["xy"])
    km = KMeans(n_clusters=5, n_init=10, random_state=0).fit(sm.reshape(-1, 1))
    order = np.argsort([sm[km.labels_ == k].mean() for k in range(5)])
    mask = km.labels_ == int(order[-1])
    lab = m["labels"]
    tls = lab == "TLS"
    neg = np.isin(lab, ["TUM", "NOR"])
    comps = components(m["xy"], mask)
    return {"section": sid, "cohort": "USZ",
            "tls_spots": int(tls.sum()), "coverage": float(mask[tls].mean()) if tls.any() else None,
            "neg_spots": int(neg.sum()),
            "fpr_on_tumor_normal": float(mask[neg].mean()) if neg.any() else None,
            "n_comp": len(comps),
            "comp_sizes": sorted([len(c) for c in comps], reverse=True)[:6]}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    idrows = list(csv.DictReader(open(
        ROOT / "infra/bioinf-data-index/raw/user_gse175540_20261003/current_26_source_identity_confirmed.tsv"),
        delimiter="\t"))
    gse = [r["section_id"].split("-", 1)[1]
           for r in idrows if r["cohort"] == "GSE175540"]
    gse += ["GSM5924036_ffpe_c_21", "GSM5924048_frozen_b_13"]
    usz = sorted(r["section_id"] for r in idrows if r["cohort"] == "USZ")
    print(f"GSE={len(gse)} USZ={len(usz)}", flush=True)
    rows = []
    for g in gse:
        try:
            r = probe_run(g)
            r["cohort"] = "GSE175540"
            rows.append(r)
        except Exception as e:
            rows.append({"gsm": g, "cohort": "GSE175540", "status": f"ERROR {e}"})
    for s in usz:
        try:
            rows.append(usz_audit(s))
        except Exception as e:
            rows.append({"section": s, "cohort": "USZ", "status": f"ERROR {e}"})
    (OUT / "audit.json").write_text(json.dumps(rows, indent=1) + "\n")
    g = [r for r in rows if r.get("cohort") == "GSE175540" and "sens" in r]
    med_sens = float(np.median([r["sens"] for r in g])) if g else None
    med_spec = float(np.median([r["spec"] for r in g])) if g else None
    verdict = (med_sens is not None and med_sens >= 0.5
               and med_spec is not None and med_spec >= 0.9)
    print(f"GSE n={len(g)} med_sens={med_sens} med_spec={med_spec} PASS={verdict}", flush=True)
    for r in rows:
        print(json.dumps(r)[:220], flush=True)
    (OUT / "STATUS").write_text(f"AUDIT n={len(rows)} PASS={verdict}\n")


if __name__ == "__main__":
    main()
