#!/usr/bin/env python3
"""GO-BUILD: GOBP 全量通路 -> 先验程序定义冻结（T3-GO）。

- GMT: c5.go.v2025.1.Hs.symbols.gmt（GOBP_ 前缀，7583 条）
- coverage: 通路基因 ∩ 三队列交集（Vanderbilt cache10k ∩ STCRC全片并集 ∩ USZ全片并集）
  <3 基因 -> NOT_TESTABLE（用户拍板阈值）
- half-split: seed 20260920 shuffle，input = ceil(n/2)（保证 input≥2），readout = 其余
  input/readout 零重叠（防泄漏铁律，机器检查）
- 去重：covered 基因集完全相同的 GO ID 只跑一份（id map 保留，上报时展开）；
  原始结果全留，家族归并只做 report 层（用户拍板）
输出：infra/r16/recovery_20260918/programs_go/{cohort_genes.json, definitions.json,
  coverage_report.tsv, not_testable.tsv, id_map.json, build_provenance.json}
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/r16/recovery_20260918/programs_go"
GMT = Path("/home/huyudi/006/data/pathway/c5.go.v2025.1.Hs.symbols.gmt")
SPLIT_SEED = 20260920


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    # 1. GMT 解析（只取 GOBP_）
    gobp: dict[str, list[str]] = {}
    with open(GMT) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0].startswith("GOBP_"):
                gobp[p[0]] = sorted(set(g for g in p[2:] if g))
    print(f"GOBP pathways: {len(gobp)}")

    # 2. 三队列基因集（Vanderbilt cache10k 由 T0 契约保证；外部两队列全片并集）
    tc = json.load(open("/tmp/three_cohort_genes.json"))
    cache = set(tc["cache10k"])
    st = set(tc["STCRC"])
    us = set(tc["USZ"])
    three = cache & st & us
    json.dump({"cache10k": sorted(cache), "STCRC_union": sorted(st),
               "USZ_union": sorted(us), "three_intersection": sorted(three)},
              open(OUT / "cohort_genes.json", "w"))
    print(f"three-cohort intersection: {len(three)}")

    # 3. coverage + half-split
    rng = np.random.default_rng(SPLIT_SEED)
    defs, cover_rows, not_rows = {}, [], []
    for gid in sorted(gobp):
        raw = gobp[gid]
        cov = sorted(set(raw) & three)
        if len(cov) < 3:
            not_rows.append((gid, len(raw), len(cov)))
            continue
        order = rng.permutation(len(cov))
        n_in = (len(cov) + 1) // 2  # ceil：input 取大半，保证 input>=2
        inp = sorted(cov[i] for i in order[:n_in])
        ro = sorted(cov[i] for i in order[n_in:])
        assert not (set(inp) & set(ro)), gid
        assert len(inp) >= 2 and len(ro) >= 1, gid
        defs[gid] = {"program_id": gid, "source": "GOBP_MSIGDB_v2025.1",
                     "input_genes": inp, "readout_genes": ro,
                     "covered_size": len(cov), "raw_size": len(raw),
                     "split_seed": SPLIT_SEED, "split_rule": "shuffle_then_ceil_input",
                     "frozen_for_external_at": "", "definition_hash": ""}
        cover_rows.append((gid, len(raw), len(cov), len(inp), len(ro)))
    print(f"testable: {len(defs)}, NOT_TESTABLE(<3): {len(not_rows)}")

    # 4. 去重：covered 集相同的只跑一份
    by_set: dict[frozenset, list[str]] = {}
    for gid, d in defs.items():
        key = frozenset(d["input_genes"]) | frozenset(d["readout_genes"])
        by_set.setdefault(key, []).append(gid)
    uniq = sorted(by_set.values(), key=lambda v: v[0])
    id_map = {v[0]: v for v in uniq}  # representative -> all GO IDs sharing the set
    n_dup_extra = sum(len(v) - 1 for v in uniq)
    json.dump(id_map, open(OUT / "id_map.json", "w"))
    print(f"unique covered sets: {len(uniq)} (+{n_dup_extra} duplicate IDs mapped)")

    json.dump(defs, open(OUT / "definitions.json", "w"))
    with open(OUT / "coverage_report.tsv", "w") as f:
        f.write("go_id\traw_size\tcovered_size\tn_input\tn_readout\trepresentative\n")
        rep_of = {}
        for rep, ids in id_map.items():
            for i in ids:
                rep_of[i] = rep
        for gid, raw, cov, ni, nr in sorted(cover_rows):
            f.write(f"{gid}\t{raw}\t{cov}\t{ni}\t{nr}\t{rep_of[gid]}\n")
    with open(OUT / "not_testable.tsv", "w") as f:
        f.write("go_id\traw_size\tcovered_size\treason\n")
        for gid, raw, cov in sorted(not_rows):
            f.write(f"{gid}\t{raw}\t{cov}\tcovered<3\n")
    prov = {"gmt": str(GMT), "gmt_sha256": sha256(GMT),
            "n_gobp": len(gobp), "n_testable": len(defs),
            "n_unique_sets": len(uniq), "n_not_testable": len(not_rows),
            "three_intersection_size": len(three),
            "split_seed": SPLIT_SEED, "coverage_rule": "covered>=3",
            "cohort_manifests": ["infra/r04/role_manifests/internal_validation_manifest.json",
                                 "infra/r04/role_manifests/external_validation_manifest.json",
                                 "infra/r16/census_cache_hvg10k (Vanderbilt 47/30)"]}
    json.dump(prov, open(OUT / "build_provenance.json", "w"), indent=2,
              ensure_ascii=False)
    print(f"GO-BUILD done -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
