#!/usr/bin/env python3
"""T3-3: 可分性 + 非线性基线敏感性 + 全患者效应表 + 程序冻结定义。

- leave-module-out：程序基因从重叠轴定义中剔除，轴分相关>0.9 才算可分，
  否则记 NOT_SEPARABLE_WITH_THIS_PROXY（03 §2.2 / 02 T3-5）。
- 非线性基线敏感性：D6 + 平方项（低自由度联合非线性），残余 Moran 若崩塌，
  则线性残余≠发现（02 T3-4）。
- program_patient_effects.tsv：47 片逐片三层 Moran → 患者中位（A18：全表，不截断）。
- program_definitions.json：冻结定义（frozen_for_external_at 为空，T4 前冻结）。
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from r16 import axes as A
from r16 import census as C

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/r16/recovery_20260918/programs"

PROGS = {
    "P-epi": (["EPCAM", "CEACAM5", "CEACAM6", "CLDN4", "CLDN7", "KRT8", "KRT18",
               "MUC13", "TFF3", "ELF3", "FXYD3", "LGALS4", "ANXA2", "PHGR1",
               "LGALS3BP", "GPX2"],
              "identity-compatible: enterocyte/goblet differentiation state; "
              "state-compatible: epithelial activation"),
    "P-stromal": (["COL1A1", "COL1A2", "COL3A1", "COL6A3", "DCN", "SPARC",
                   "IGFBP7", "VIM", "C1R"],
                  "identity-compatible: fibroblast/matrix abundance; "
                  "state-compatible: matrix remodeling activity"),
    "P-plasma": (["IGHA1", "IGKC", "IGLC2", "IGLC3", "JCHAIN", "IGHV3-30",
                  "IGHV3-35", "IGKV1D-12", "IGKV3OR2-268",
                  "DEPRECATED_ENSG00000211685", "DEPRECATED_ENSG00000211890"],
                 "identity-compatible: plasma cell abundance"),
    "P-smmhc": (["ACTA2", "FLNA", "MYL9"],
                "identity-compatible: smooth muscle/myofibroblast; "
                "state-compatible: contractile state"),
    "P-mhc2": (["CD74", "HLA-DRA"],
               "state-compatible: antigen-presentation activity; "
               "identity-compatible: APC (incl. B/myeloid) abundance"),
    "P-stress": (["FOS", "DEPRECATED_ENSG00000170345"],
                 "state-compatible: immediate-early stress response; "
                 "technical-compatible: procedure/permeabilization gradient "
                 "not excluded"),
}


def main() -> int:
    cache = ROOT / "infra/r16/census_cache_hvg10k/sections"
    proxy = json.load(open(ROOT / "infra/r04/marker_proxy_combined.json"))
    union = json.load(open(ROOT / "infra/r16/hvg_union_genes.json"))
    rank = np.load(ROOT / "infra/r16/hvg_final_rank_20260914.npy")
    top10k = {union[i] for i in np.flatnonzero(rank < 10000)}
    axis_defs = {k: [g for g in v["voted_genes"] if g in top10k]
                 for k, v in proxy["classes"].items() if k not in A.RETIRED_AXES}
    axis_defs["Plasma"] = [g for g in A.PLASMA_GENES if g in top10k]
    stems = C.list_stems(cache)
    _, _, genes0, _, _ = C.load_section(cache, stems[0])

    # ---- 1. leave-module-out separability (one mid-size section suffices
    # for the correlation diagnostic; axis defs are global)
    stem0 = sorted(stems, key=lambda s: s)[len(stems) // 2]
    mat, _, genes, coords, meta = C.load_section(cache, stem0)
    gx = {g: i for i, g in enumerate(genes)}
    x_z = C.log1p_zscore(mat)
    print(f"separability section: {stem0[:30]} n={len(coords)}")
    sep = {}
    for pn, (gl, interp) in PROGS.items():
        ov = {a: sorted(set(gs) & set(gl)) for a, gs in axis_defs.items()}
        ov = {a: v for a, v in ov.items() if v}
        rec = {"overlapping_axes": {a: {"n_overlap": len(v), "genes": v}
                                    for a, v in ov.items()},
               "axes_separable": True, "notes": []}
        for a, v in ov.items():
            rest = [g for g in axis_defs[a] if g not in set(gl)]
            s_full, n_full = A.axis_scores(x_z, gx, axis_defs[a])
            s_loo, n_loo = A.axis_scores(x_z, gx, rest)
            if n_loo == 0:
                rec["axes_separable"] = False
                rec["notes"].append(f"{a}: no genes left -> NOT_SEPARABLE_WITH_THIS_PROXY")
                continue
            r = float(np.corrcoef(s_full, s_loo)[0, 1])
            rec["overlapping_axes"][a]["loo_corr"] = round(r, 4)
            rec["overlapping_axes"][a]["n_left"] = n_loo
            if r < 0.9:
                rec["axes_separable"] = False
                rec["notes"].append(f"{a}: loo_corr={r:.3f} <0.9 -> NOT_SEPARABLE_WITH_THIS_PROXY")
        sep[pn] = rec
        flag = "SEPARABLE" if rec["axes_separable"] else "NOT_SEPARABLE"
        print(f"  {pn}: {flag} overlap={ {a: v['n_overlap'] for a, v in rec['overlapping_axes'].items()} }")
    json.dump(sep, open(OUT / "separability.json", "w"), ensure_ascii=False, indent=2)

    # ---- 2+3. 全片循环：三层 Moran + 非线性敏感性（4 代表片） + 患者效应
    nl_stems = {"0005_HTAN_8270_AS_10_filtered_trimmed.h5ad",
                "0021_HTAN_6723_KL_3_filtered_trimmed.h5ad",
                "0009_HTAN_8578_AS_1_filtered_trimmed_WD33473.h5ad",
                "0046_HTAN_7003_AS_8_filtered_trimmed.h5ad"}
    eff_rows = []
    nl_rows = []
    for si, stem in enumerate(stems):
        mat, bcs, genes, coords, meta = C.load_section(cache, stem)
        gx = {g: i for i, g in enumerate(genes)}
        Xl = np.log1p(mat.toarray().astype(np.float32))
        lib = np.asarray(mat.sum(axis=1)).ravel()
        ng = np.asarray((mat > 0).sum(axis=1)).ravel()
        x_z = C.log1p_zscore(mat)
        D6 = np.column_stack([A.axis_scores(x_z, gx, gs)[0][:, None]
                              for gs in axis_defs.values()])
        Q = np.column_stack([np.ones(len(coords)), np.log10(lib + 1),
                             (ng - ng.mean()) / (ng.std() + 1e-9)])
        D = np.column_stack([Q, D6])
        # 非线性基线：D + 六轴平方项（低自由度联合非线性，02 T3-4）
        Z6 = (D6 - D6.mean(axis=0)) / (D6.std(axis=0) + 1e-9)
        Dnl = np.column_stack([D, Z6 ** 2])
        w = C.spatial_weights(coords, C.KNN_SPATIAL)
        for pn, (gl, interp) in PROGS.items():
            cols = [gx[g] for g in gl if g in gx]
            cov = f"{len(cols)}/{len(gl)}"
            if not cols:
                eff_rows.append((pn, meta, len(coords), cov, *[None] * 6, "NO_COVERAGE"))
                continue
            raw = Xl[:, cols].mean(axis=1)
            rQ = raw - Q @ np.linalg.lstsq(Q, raw, rcond=None)[0]
            rQC = raw - D @ np.linalg.lstsq(D, raw, rcond=None)[0]
            eff_rows.append((pn, meta, len(coords), cov,
                             round(float(C.moran_i(raw, w)), 4),
                             round(float(C.moran_i(rQ, w)), 4),
                             round(float(C.moran_i(rQC, w)), 4),
                             None, None, "OK"))
            if stem in nl_stems:
                rNL = raw - Dnl @ np.linalg.lstsq(Dnl, raw, rcond=None)[0]
                nl_rows.append({"section": stem[:30], "program": pn,
                                "moran_QC_lin": round(float(C.moran_i(rQC, w)), 4),
                                "moran_QC_quad": round(float(C.moran_i(rNL, w)), 4)})
        if (si + 1) % 10 == 0:
            print(f"  [{si+1}/{len(stems)}]", flush=True)
    with open(OUT / "program_section_effects.tsv", "w") as f:
        f.write("program\tpatient\tsection\tn_spots\tgene_coverage\t"
                "moran_raw\tmoran_resid_Q\tmoran_resid_QC\tstatus\n")
        for pn, meta, n, cov, m0, m1, m2, _, _, st in eff_rows:
            f.write(f"{pn}\t{meta['patient_id']}\t{meta['section_id']}\t{n}\t{cov}\t"
                    f"{m0}\t{m1}\t{m2}\t{st}\n")
    # 患者中位表（患者等权）
    pat: dict = defaultdict(list)
    for pn, meta, n, cov, m0, m1, m2, _, _, st in eff_rows:
        if st == "OK":
            pat[(pn, meta["patient_id"])].append((m0, m1, m2, n))
    with open(OUT / "program_patient_effects.tsv", "w") as f:
        f.write("program\tpatient\tn_sections\tmed_spots\t"
                "med_moran_raw\tmed_moran_Q\tmed_moran_QC\tstatus\n")
        for (pn, p), vs in sorted(pat.items()):
            a = np.median([v[0] for v in vs])
            b = np.median([v[1] for v in vs])
            c = np.median([v[2] for v in vs])
            f.write(f"{pn}\t{p}\t{len(vs)}\t{int(np.median([v[3] for v in vs]))}\t"
                    f"{a:.4f}\t{b:.4f}\t{c:.4f}\tOK\n")
    print("nonlinear sensitivity (QC_lin -> QC_quad):")
    for r in nl_rows:
        print(f"  {r['section'][:22]} {r['program']}: {r['moran_QC_lin']} -> {r['moran_QC_quad']}")
    json.dump(nl_rows, open(OUT / "nonlinear_sensitivity.json", "w"), indent=2)

    # ---- 4. program_definitions.json（冻结结构，时间戳空，T4 前冻结）
    defs = {}
    for pn, (gl, interp) in PROGS.items():
        cols = [g for g in gl if g in set(genes0)]
        defs[pn] = {
            "program_id": pn, "source_candidates": "laneB191+Tier2-13+jointABGH",
            "discovery_patients": "HTAN_VANDERBILT_CRC 30 patients",
            "gene_weights": {g: 1.0 / len(cols) for g in cols},
            "input_genes": cols, "readout_genes": "TBD_T4_DISJOINT",
            "direction_rule": "mean of log1p, higher = stronger",
            "normalization_rule": "log1p(counts) per spot; Q/QC residual variants defined",
            "nuisance_rule": "Q=[1,log10(libsize),ngenes_z]; QC=Q+6axes(linear); QCnlin sensitivity",
            "geometry_definition": "signed hops to program-own q0.7 contour (T4); NA until then",
            "coverage_policy": f"all {len(cols)} genes in cache10k; section min genes {len(cols)}/{len(gl)}",
            "interpretation": interp,
            "separability": sep[pn],
            "frozen_for_external_at": "", "definition_hash": "",
        }
    json.dump(defs, open(OUT / "program_definitions.json", "w"),
              ensure_ascii=False, indent=2)
    print(f"effects + definitions done -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
