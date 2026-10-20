"""生成患者级分组 fold 分配（SPLIT_FREEZE_20261014.md 的执行件）。

确定性：seed 固定；输入为 sections.tsv + 声明合并表。任何修订必须改版本号。
"""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/tls_hooknet_rna_triage_20261014"
SECTIONS = ROOT / "infra/tls_four_cancer_labels_20261008/eligibility/sections.tsv"
IDX = ROOT / "dist/tls_annotation_pack/data/index.json"
VERSION = "v1"
SEED = 20261014
K = 5

FORCE_SAME_UNIT = [("TENX23", "TENX24"), ("TENX13", "TENX14")]


def main():
    slides = [s["section_id"] for s in json.loads(IDX.read_text())["slides"]]
    with open(SECTIONS) as f:
        sec = {r["section_id"]: r for r in csv.DictReader(f, delimiter="\t")}

    unit_of = {}
    for sid in slides:
        r = sec[sid]
        unit_of[sid] = (r["patient_key"] if r["patient_status"].startswith("known")
                        else f"singleton::{sid}")
    for a, b in FORCE_SAME_UNIT:
        if a in unit_of and b in unit_of:
            unit_of[b] = unit_of[a]

    units = {}
    for sid, u in unit_of.items():
        units.setdefault(u, {"organ": sec[sid]["organ"], "slides": []})["slides"].append(sid)

    # 分层：按 organ 分桶后各桶内按哈希序轮转发牌，保证 fold 器官结构均衡
    bucket = {}
    for u, d in units.items():
        bucket.setdefault(d["organ"], []).append(u)
    fold_of = {}
    for organ, us in bucket.items():
        us.sort(key=lambda u: hashlib.sha256(f"{SEED}|{u}".encode()).hexdigest())
        for i, u in enumerate(us):
            fold_of[u] = i % K

    rows = []
    for sid in sorted(slides):
        u = unit_of[sid]
        rows.append({"section_id": sid, "unit_id": u,
                     "organ": sec[sid]["organ"], "fold": fold_of[u]})
    out = OUT / f"split_folds_{VERSION}.tsv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["section_id", "unit_id", "organ", "fold"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    sha = hashlib.sha256(out.read_bytes()).hexdigest()
    n_units = len({r["unit_id"] for r in rows})
    from collections import Counter
    print(f"units: {n_units}, slides: {len(rows)}, folds: {K}")
    print("fold sizes:", dict(Counter(r["fold"] for r in rows)))
    print("per-organ units:", {o: len(v) for o, v in bucket.items()})
    print(f"written {out} sha256={sha}")


if __name__ == "__main__":
    main()
