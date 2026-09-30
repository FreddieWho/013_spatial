#!/usr/bin/env python3
"""T4: 冻结程序的独立队列复现 + M0/M1/M2 增量预测（02 T4，03 §6）。

设计（防泄漏为硬约束）：
- 程序 input 基因权重冻结（Vanderbilt 均值权重）；readout 基因与 input 零重叠
  （机器检查）；与轴基因的重叠走 leave-axis-out 敏感性分支。
- 三层终点分开：分子共表达再现（readout内相关）/ 空间形状再现（program轮廓
  Moran + 与组成锚点的关系）/ 预测增量（M0=Q+C vs M1=Q+C+程序非空间活性
  vs M2=+邻域空间项，留出患者，读出=readout均值）。
- 患者先汇总；两队列全程序全披露；exposure 记录（已知轴已在两队列打分过）。
- P-plasma：只验分子再现（NOT_SEPARABLE，T3 已降级），不争增量。
输出：external_replication.tsv, incremental_prediction.tsv,
      coverage_shift_report.md
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy import sparse

from r16 import axes as A
from r16 import census as C
from r16.section_io import load_section_symbols

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/r16/recovery_20260918"

READOUT = {
    "P-epi": ["MUC2", "MUC4", "TFF1", "CLDN3", "KRT19", "CDH1"],
    "P-stromal": ["COL6A2", "COL5A1", "FBN1", "BGN", "LUM", "FN1"],
    "P-plasma": ["IGHM", "IGHA2", "MZB1", "SDC1", "TNFRSF17", "CD38"],
    "P-smmhc": ["TAGLN", "CNN1", "DES", "MYH11", "TPM2"],
    "P-mhc2": ["HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "CD86"],
    "P-stress": ["JUN", "JUNB", "EGR1", "IER2", "DUSP1", "ATF3"],
}


def prog_score(counts, gidx: dict, genes: list[str]) -> np.ndarray:
    cols = [gidx[g] for g in genes if g in gidx]
    if not cols:
        return np.full(counts.shape[0], np.nan)
    X = counts[:, cols]
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    return np.log1p(X.astype(float)).mean(axis=1)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=C.SEED)
    args = ap.parse_args()
    t0 = time.time()
    defs = json.load(open(OUT / "programs/program_definitions.json"))
    # 冻结：填时间戳 + hash + readout（T4 前最后一次写入定义）
    frozen_at = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    for pn, gl in READOUT.items():
        allinp = {g for d in defs.values() for g in d["input_genes"]}
        assert not (set(gl) & allinp), f"{pn} readout overlaps input!"
        assert not (set(gl) & set(defs[pn]["input_genes"])), f"{pn} self overlap!"
        defs[pn]["readout_genes"] = gl
        defs[pn]["frozen_for_external_at"] = frozen_at
    h = hashlib.sha256(json.dumps(defs, sort_keys=True).encode()).hexdigest()[:16]
    for pn in defs:
        defs[pn]["definition_hash"] = h
    json.dump(defs, open(OUT / "programs/program_definitions.json", "w"),
              ensure_ascii=False, indent=2)
    print(f"froze definitions at {frozen_at} hash={h}")

    proxy = json.load(open(ROOT / "infra/r04/marker_proxy_combined.json"))
    top = None  # 外部队列不用 Vanderbilt HVG 过滤：用各片实际基因（D-116 口径）
    axis_full = {k: list(v["voted_genes"]) for k, v in proxy["classes"].items()
                 if k not in A.RETIRED_AXES}
    axis_full["Plasma"] = list(A.PLASMA_GENES)

    rep_rows, inc_rows, cov_rows = [], [], []
    for mf, cohort in [("infra/r04/role_manifests/internal_validation_manifest.json", "ST-CRC"),
                       ("infra/r04/role_manifests/external_validation_manifest.json", "USZ")]:
        man = json.load(open(ROOT / mf))
        # 每片：input score / readout score / Q / C(六轴，用片内基因)
        sec_data = []
        for row in man["rows"]:
            sec = load_section_symbols(row)
            genes = [g if isinstance(g, str) else g.decode() for g in sec.gene_id]
            gidx = {}
            for i, g in enumerate(genes):
                if g not in gidx:
                    gidx[g] = i
            counts = sec.counts.tocsr() if sparse.issparse(sec.counts) \
                else sparse.csr_matrix(sec.counts)
            lib = np.asarray(counts.sum(axis=1)).ravel()
            ng = np.asarray((counts > 0).sum(axis=1)).ravel()
            coords = np.asarray(sec.coords, dtype=float)
            # Q
            Q = np.column_stack([np.ones(len(coords)), np.log10(lib + 1),
                                 (ng - ng.mean()) / (ng.std() + 1e-9)])
            # C：六轴（log1p_zscore 片内，与发现集同定义；leave-axis-out 另分支记）
            Xd = counts.toarray().astype(np.float32)
            np.log1p(Xd, out=Xd)
            mu, sd = Xd.mean(0), Xd.std(0)
            sd[sd < 1e-8] = 1.0
            Xz = (Xd - mu) / sd
            C6 = np.column_stack([A.axis_scores(Xz, gidx, gs)[0][:, None]
                                  for gs in axis_full.values()])
            scores = {}
            for pn in defs:
                scores[pn + ":in"] = prog_score(counts, gidx, defs[pn]["input_genes"])
                scores[pn + ":ro"] = prog_score(counts, gidx, defs[pn]["readout_genes"])
            sec_data.append({"patient": sec.patient_id, "section": sec.section_id,
                             "n": len(coords), "coords": coords, "Q": Q, "C6": C6,
                             "scores": scores,
                             "cov": {pn: f"{sum(1 for g in defs[pn]['input_genes'] if g in gidx)}/{len(defs[pn]['input_genes'])}"
                                     for pn in defs}})
            cov_rows.append({"cohort": cohort, "section": sec.section_id,
                             "patient": sec.patient_id, "n_spots": len(coords),
                             ** {f"cov_{pn}": sec_data[-1]["cov"][pn] for pn in defs}})
        print(f"[{cohort}] loaded {len(sec_data)} sections", flush=True)
        w_cache = {}
        # 分子再现：患者内 input-readout 相关（中位）+ readout 内相关
        for pn in defs:
            for sd in sec_data:
                si, ro = sd["scores"][pn + ":in"], sd["scores"][pn + ":ro"]
                ok = np.isfinite(si) & np.isfinite(ro)
                sd.setdefault("mol", {})[pn] = float(np.corrcoef(si[ok], ro[ok])[0, 1]) \
                    if ok.sum() > 10 else np.nan
            by_pat = {}
            for sd in sec_data:
                by_pat.setdefault(sd["patient"], []).append(sd["mol"][pn])
            for p, vs in sorted(by_pat.items()):
                rep_rows.append({"program": pn, "cohort": cohort, "patient": p,
                                 "n_sections": len(vs),
                                 "mol_in_ro_corr_med": round(float(np.nanmedian(vs)), 4),
                                 "level": "molecular"})
        # 空间再现：readout 轮廓 Moran（患者中位）+ input 轮廓 Moran
        for pn in defs:
            for sd in sec_data:
                key = (sd["section"],)
                if key not in w_cache:
                    w_cache[key] = C.spatial_weights(sd["coords"], C.KNN_SPATIAL)
                w = w_cache[key]
                sd.setdefault("spa", {})[pn + ":in"] = C.moran_i(sd["scores"][pn + ":in"], w)
                sd["spa"][pn + ":ro"] = C.moran_i(sd["scores"][pn + ":ro"], w)
            by_pat_in, by_pat_ro = {}, {}
            for sd in sec_data:
                by_pat_in.setdefault(sd["patient"], []).append(sd["spa"][pn + ":in"])
                by_pat_ro.setdefault(sd["patient"], []).append(sd["spa"][pn + ":ro"])
            for p in sorted(by_pat_in):
                rep_rows.append({"program": pn, "cohort": cohort, "patient": p,
                                 "n_sections": len(by_pat_in[p]),
                                 "spat_in_moran_med": round(float(np.nanmedian(by_pat_in[p])), 4),
                                 "spat_ro_moran_med": round(float(np.nanmedian(by_pat_ro[p])), 4),
                                 "level": "spatial"})
        # 增量预测：留一患者；读出=readout均值；M0=Q+C6, M1=+input分, M2=+邻域均值
        pats = sorted({sd["patient"] for sd in sec_data})
        for pn in defs:
            if pn == "P-plasma":
                inc_rows.append({"program": pn, "cohort": cohort, "patient": "ALL",
                                 "n_train_patients": "NA", "mse_M0": "NA",
                                 "mse_M1": "NA", "mse_M2": "NA",
                                 "note": "NOT_SEPARABLE降级：只验分子再现，不争增量"})
                continue
            for p_held in pats:
                tr = [sd for sd in sec_data if sd["patient"] != p_held]
                te = [sd for sd in sec_data if sd["patient"] == p_held]
                Xtr0 = np.vstack([np.column_stack([sd["Q"], sd["C6"]]) for sd in tr])
                ytr = np.concatenate([sd["scores"][pn + ":ro"] for sd in tr])
                Xte0 = np.vstack([np.column_stack([sd["Q"], sd["C6"]]) for sd in te])
                yte = np.concatenate([sd["scores"][pn + ":ro"] for sd in te])
                Ptr_in = np.concatenate([sd["scores"][pn + ":in"] for sd in tr])
                Pte_in = np.concatenate([sd["scores"][pn + ":in"] for sd in te])
                # M2 邻域项：同片 input 分的 kNN 均值（可见点分子值允许，03 M2 定义）
                def neigh_mean(sds, key):
                    from sklearn.neighbors import NearestNeighbors
                    outs = []
                    for sd in sds:
                        nn = NearestNeighbors(n_neighbors=9).fit(sd["coords"])
                        _, idx = nn.kneighbors(sd["coords"])
                        v = sd["scores"][key]
                        outs.append(v[idx[:, 1:]].mean(axis=1))
                    return np.concatenate(outs)
                Ntr = neigh_mean(tr, pn + ":in")
                Nte = neigh_mean(te, pn + ":in")
                Xtr1 = np.column_stack([Xtr0, Ptr_in])
                Xte1 = np.column_stack([Xte0, Pte_in])
                Xtr2 = np.column_stack([Xtr1, Ntr])
                Xte2 = np.column_stack([Xte1, Nte])
                mses = []
                for Xtr, Xte in [(Xtr0, Xte0), (Xtr1, Xte1), (Xtr2, Xte2)]:
                    beta, _, _, _ = np.linalg.lstsq(Xtr, ytr, rcond=None)
                    mses.append(float(((Xte @ beta - yte) ** 2).mean()))
                inc_rows.append({"program": pn, "cohort": cohort, "patient": p_held,
                                 "n_train_patients": len(pats) - 1,
                                 "mse_M0": round(mses[0], 5), "mse_M1": round(mses[1], 5),
                                 "mse_M2": round(mses[2], 5),
                                 "dM1_pct": round(100 * (mses[0] - mses[1]) / mses[0], 2),
                                 "dM2_pct": round(100 * (mses[1] - mses[2]) / mses[1], 2)
                                 if mses[1] > 0 else "NA",
                                 "note": ""})
    with open(OUT / "external_replication.tsv", "w") as f:
        cols = ["program", "cohort", "patient", "n_sections", "level",
                "mol_in_ro_corr_med", "spat_in_moran_med", "spat_ro_moran_med"]
        f.write("\t".join(cols) + "\n")
        for r in rep_rows:
            f.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")
    with open(OUT / "incremental_prediction.tsv", "w") as f:
        cols = ["program", "cohort", "patient", "n_train_patients",
                "mse_M0", "mse_M1", "mse_M2", "dM1_pct", "dM2_pct", "note"]
        f.write("\t".join(cols) + "\n")
        for r in inc_rows:
            f.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")
    with open(OUT / "coverage_shift_report.md", "w") as f:
        f.write("# T4 coverage（input基因在外部队列的覆盖；readout三队列全覆盖已验）\n\n")
        f.write("cohort|section|patient|n_spots|" + "|".join(f"cov_{p}" for p in defs) + "\n")
        f.write("---|---|---|---" + "|---" * len(defs) + "\n")
        for r in cov_rows:
            f.write(f"{r['cohort']}|{r['section'].split('::')[-1][:28]}|"
                    f"{r['patient'].split('::')[-1][:8]}|{r['n_spots']}|"
                    + "|".join(str(r[f"cov_{p}"]) for p in defs) + "\n")
        f.write("\n## exposure 记录\n已知六轴曾在两队列打分（D-116）；新程序定义Vanderbilt冻结、外部首次应用，\n"
                "但M0/M1/M2的回归系数在留出框架内估计（非零样本外推，03 §6.1 局部校准口径）。\n"
                "P-plasma 不争增量（NOT_SEPARABLE）。USZ跨癌种：CRC特有上皮程序不要求在肾/肺成立。\n")
    print(f"T4 done ({time.time()-t0:.0f}s): rep={len(rep_rows)} inc={len(inc_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
