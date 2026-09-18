#!/usr/bin/env python3
"""T3-1: candidate union — 191 基因 + 联合臂/Tier-2 marker，带去重表。

输入：laneB_residual JSON（191 全记录）+ joint A/B/G/H + Tier-2 census。
输出：infra/r16/recovery_20260918/programs/candidate_union.tsv
列：gene | in_laneB191 | laneB_I_med | laneB_npat_coh | sources | n_sources
    | npat_max | axis_overlap | axis_name
不过滤 191（全部保留）；只组织 + 标注来源。
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/r16/recovery_20260918/programs"

AXIS_OF = {}
proxy = json.load(open(ROOT / "infra/r04/marker_proxy_combined.json"))
for k, v in proxy["classes"].items():
    if k == "ILC":
        continue
    for g in v.get("voted_genes", []):
        AXIS_OF.setdefault(g, []).append(k)
for g in ["IGHM", "IGKC", "JCHAIN", "MZB1", "SDC1", "IGLC3", "IGHA2"]:
    AXIS_OF.setdefault(g, []).append("Plasma")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    lb = json.load(open(ROOT / "infra/r16/laneB_residual_20260917.json"))
    lb191 = {x["gene"]: x for x in lb["results"]
             if x.get("grade") == "EXPLORATORY_REPRODUCED"}
    t2 = json.load(open(ROOT / "infra/r16/composition_pattern_census_20260914.json"))
    t2rep = [g for g in t2["pattern_groups"]
             if g.get("evidence_grade") == "EXPLORATORY_REPRODUCED"]

    src: dict[str, dict] = {}
    for g in lb191:
        src.setdefault(g, {"sources": set(), "npat_max": 0})
        src[g]["sources"].add("laneB191")
        src[g]["npat_max"] = max(src[g]["npat_max"], lb191[g]["n_patients_coherent"])
    for g in t2rep:
        for m in g["top_markers"]:
            s = src.setdefault(m, {"sources": set(), "npat_max": 0})
            s["sources"].add(f"tier2:{g['pattern_id']}(npat={g['n_patients']})")
            s["npat_max"] = max(s["npat_max"], g["n_patients"])
    for f, tag in [("infra/r16/joint_embed_20260915.json", "AB"),
                   ("infra/r16/joint_embed_G_20260916.json", "G"),
                   ("infra/r16/joint_embed_H_20260917.json", "H")]:
        j = json.load(open(ROOT / f))
        for c in j["candidates"]:
            if c.get("n_patients", 0) < 2:
                continue
            cid = f"{tag}:{c.get('resolution', c.get('col'))}/{c['cluster_id']}"
            for m in c["markers"][:10]:
                s = src.setdefault(m, {"sources": set(), "npat_max": 0})
                s["sources"].add(f"{cid}(npat={c['n_patients']})")
                s["npat_max"] = max(s["npat_max"], c["n_patients"])

    with open(OUT / "candidate_union.tsv", "w") as fh:
        fh.write("gene\tin_laneB191\tlaneB_I_med\tlaneB_npat_coh\t"
                 "n_sources\tnpat_max\taxis_overlap\tsources\n")
        for g in sorted(src):
            in_lb = g in lb191
            fh.write(f"{g}\t{int(in_lb)}\t"
                     f"{lb191[g]['resid_moran_median'] if in_lb else 'NA'}\t"
                     f"{lb191[g]['n_patients_coherent'] if in_lb else 'NA'}\t"
                     f"{len(src[g]['sources'])}\t{src[g]['npat_max']}\t"
                     f"{';'.join(AXIS_OF.get(g, [])) or '-'}\t"
                     f"{'|'.join(sorted(src[g]['sources']))}\n")
    n_lb = sum(1 for g in src if g in lb191)
    print(f"union={len(src)} genes ({n_lb} in laneB191, "
          f"{len(src)-n_lb} joint/tier2-only) -> {OUT/'candidate_union.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
