#!/usr/bin/env python3
"""GO-RUN: 6870 unique GOBP 程序 × T4 同款三层验证（发现集 QC + ST-CRC/USZ M0/M1/M2）。

复用 T4 机器（逐字同口径）：
- Vanderbilt 发现集：program_patient_effects 同款（raw/Q/QC 三层 Moran，患者中位）
- ST-CRC/USZ：分子再现（input-readout corr）+ 空间（Moran）+ 留一患者 M0/M1/M2
- 流式：三队列矩阵各读一次；Vanderbilt 用 cache npz（counts），外部用 section_io
- 输出：programs_go/{discovery_effects.tsv, external_replication.tsv,
  incremental_prediction.tsv} + run log
性能：6870 × (47 Vanderbilt + 22 外部) 片；Vanderbilt 逐片 OLS 三层是主要成本。
"""

from __future__ import annotations

import argparse
import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import sparse

from r16 import axes as A
from r16 import census as C
from r16.section_io import load_section_symbols

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/r16/recovery_20260918/programs_go"


def prog_score_mat(X: np.ndarray, cols: list[int]) -> np.ndarray:
    return np.log1p(X[:, cols].astype(float)).mean(axis=1)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="smoke: only first N programs")
    ap.add_argument("--tag", type=str, default="full")
    args = ap.parse_args()
    t0 = time.time()

    defs = json.load(open(OUT / "definitions.json"))
    id_map = json.load(open(OUT / "id_map.json"))
    reps = sorted(id_map)  # 只跑 representative
    if args.limit:
        reps = reps[:args.limit]
    print(f"programs to run: {len(reps)}", flush=True)

    proxy = json.load(open(ROOT / "infra/r04/marker_proxy_combined.json"))
    union = json.load(open(ROOT / "infra/r16/hvg_union_genes.json"))
    rank = np.load(ROOT / "infra/r16/hvg_final_rank_20260914.npy")
    top10k = {union[i] for i in np.flatnonzero(rank < 10000)}
    axis_defs = {k: [g for g in v["voted_genes"] if g in top10k]
                 for k, v in proxy["classes"].items() if k not in A.RETIRED_AXES}
    axis_defs["Plasma"] = [g for g in A.PLASMA_GENES if g in top10k]
    axis_full = {k: list(v["voted_genes"]) for k, v in proxy["classes"].items()
                 if k not in A.RETIRED_AXES}
    axis_full["Plasma"] = list(A.PLASMA_GENES)

    # ---- Vanderbilt 发现集：cache npz 流式 ----
    cache = ROOT / "infra/r16/census_cache_hvg10k/sections"
    stems = C.list_stems(cache)
    _, _, genes0, _, _ = C.load_section(cache, stems[0])
    genes0 = list(genes0)
    g0 = {g: i for i, g in enumerate(genes0)}
    # 预查 coverage（定义保证三队列交集内，Vanderbilt 全覆盖）
    prog_cols = {}
    for gid in reps:
        d = defs[gid]
        ic = [g0[g] for g in d["input_genes"]]
        rc = [g0[g] for g in d["readout_genes"]]
        prog_cols[gid] = (ic, rc)

    disc_rows = []  # (gid, patient, med_raw, med_Q, med_QC, nsec)
    pat_acc: dict = defaultdict(lambda: defaultdict(list))
    for si, stem in enumerate(stems):
        mat, _, _, coords, meta = C.load_section(cache, stem)
        gx = None  # cache 基因序全局一致，用 g0
        X = mat.toarray().astype(np.float32)
        Xl = np.log1p(X)
        lib = np.asarray(mat.sum(axis=1)).ravel()
        ng = np.asarray((mat > 0).sum(axis=1)).ravel()
        x_z = C.log1p_zscore(mat)
        D6 = np.column_stack([A.axis_scores(x_z, g0, gs)[0][:, None]
                              for gs in axis_defs.values()])
        Q = np.column_stack([np.ones(len(coords)), np.log10(lib + 1),
                             (ng - ng.mean()) / (ng.std() + 1e-9)])
        D = np.column_stack([Q, D6])
        w = C.spatial_weights(coords, C.KNN_SPATIAL)
        bQ = {}
        for gid in reps:
            ic, rc = prog_cols[gid]
            raw = Xl[:, ic].mean(axis=1)
            rQ = raw - Q @ np.linalg.lstsq(Q, raw, rcond=None)[0]
            rQC = raw - D @ np.linalg.lstsq(D, raw, rcond=None)[0]
            pat_acc[gid][meta["patient_id"]].append(
                (C.moran_i(raw, w), C.moran_i(rQ, w), C.moran_i(rQC, w)))
        if (si + 1) % 10 == 0:
            print(f"  disc [{si+1}/{len(stems)}] ({time.time()-t0:.0f}s)", flush=True)
    with open(OUT / f"discovery_effects_{args.tag}.tsv", "w") as f:
        f.write("program\tn_input\tn_readout\tpatient\tn_sections\t"
                "med_moran_raw\tmed_moran_Q\tmed_moran_QC\n")
        for gid in reps:
            d = defs[gid]
            for p, vs in sorted(pat_acc[gid].items()):
                a = float(np.median([v[0] for v in vs]))
                b = float(np.median([v[1] for v in vs]))
                c = float(np.median([v[2] for v in vs]))
                f.write(f"{gid}\t{len(d['input_genes'])}\t{len(d['readout_genes'])}\t"
                        f"{p}\t{len(vs)}\t{a:.4f}\t{b:.4f}\t{c:.4f}\n")
    print(f"discovery done ({time.time()-t0:.0f}s)", flush=True)

    # ---- 外部两队列 ----
    rep_rows, inc_rows = [], []
    for mf, cohort in [("infra/r04/role_manifests/internal_validation_manifest.json", "ST-CRC"),
                       ("infra/r04/role_manifests/external_validation_manifest.json", "USZ")]:
        man = json.load(open(ROOT / mf))
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
            Q = np.column_stack([np.ones(len(coords)), np.log10(lib + 1),
                                 (ng - ng.mean()) / (ng.std() + 1e-9)])
            Xd = counts.toarray().astype(np.float32)
            np.log1p(Xd, out=Xd)
            mu, sd = Xd.mean(0), Xd.std(0)
            sd[sd < 1e-8] = 1.0
            Xz = (Xd - mu) / sd
            C6 = np.column_stack([A.axis_scores(Xz, gidx, gs)[0][:, None]
                                  for gs in axis_full.values()])
            Xc = counts.toarray().astype(float)
            scores = {}
            for gid in reps:
                d = defs[gid]
                ic = [gidx[g] for g in d["input_genes"] if g in gidx]
                rc = [gidx[g] for g in d["readout_genes"] if g in gidx]
                scores[gid] = (np.log1p(Xc[:, ic]).mean(axis=1) if ic else np.full(len(coords), np.nan),
                               np.log1p(Xc[:, rc]).mean(axis=1) if rc else np.full(len(coords), np.nan),
                               len(ic), len(rc))
            sec_data.append({"patient": sec.patient_id, "section": sec.section_id,
                             "n": len(coords), "coords": coords, "Q": Q, "C6": C6,
                             "scores": scores})
        print(f"[{cohort}] loaded {len(sec_data)} sections", flush=True)
        w_cache = {}
        for gid in reps:
            # 分子 + 空间
            for sd in sec_data:
                si, ro = sd["scores"][gid][:2]
                ok = np.isfinite(si) & np.isfinite(ro)
                sd.setdefault("mol", {})[gid] = float(np.corrcoef(si[ok], ro[ok])[0, 1]) \
                    if ok.sum() > 10 else np.nan
                key = (sd["section"],)
                if key not in w_cache:
                    w_cache[key] = C.spatial_weights(sd["coords"], C.KNN_SPATIAL)
                w = w_cache[key]
                sd.setdefault("spa", {})[gid + ":in"] = C.moran_i(si, w)
                sd["spa"][gid + ":ro"] = C.moran_i(ro, w)
            by = defaultdict(list)
            byi, byr = defaultdict(list), defaultdict(list)
            for sd in sec_data:
                by[sd["patient"]].append(sd["mol"][gid])
                byi[sd["patient"]].append(sd["spa"][gid + ":in"])
                byr[sd["patient"]].append(sd["spa"][gid + ":ro"])
            for p in sorted(by):
                rep_rows.append({"program": gid, "cohort": cohort, "patient": p,
                                 "n_sections": len(by[p]),
                                 "mol": round(float(np.nanmedian(by[p])), 4),
                                 "sp_in": round(float(np.nanmedian(byi[p])), 4),
                                 "sp_ro": round(float(np.nanmedian(byr[p])), 4)})
            # 增量：留一患者
            pats = sorted({sd["patient"] for sd in sec_data})
            for p_held in pats:
                tr = [sd for sd in sec_data if sd["patient"] != p_held]
                te = [sd for sd in sec_data if sd["patient"] == p_held]
                Xtr0 = np.vstack([np.column_stack([sd["Q"], sd["C6"]]) for sd in tr])
                ytr = np.concatenate([sd["scores"][gid][1] for sd in tr])
                Xte0 = np.vstack([np.column_stack([sd["Q"], sd["C6"]]) for sd in te])
                yte = np.concatenate([sd["scores"][gid][1] for sd in te])
                Ptr = np.concatenate([sd["scores"][gid][0] for sd in tr])
                Pte = np.concatenate([sd["scores"][gid][0] for sd in te])

                def neigh(sds, k_):
                    from sklearn.neighbors import NearestNeighbors
                    outs = []
                    for sd in sds:
                        nn = NearestNeighbors(n_neighbors=9).fit(sd["coords"])
                        _, idx = nn.kneighbors(sd["coords"])
                        outs.append(sd["scores"][gid][k_] [idx[:, 1:]].mean(axis=1))
                    return np.concatenate(outs)
                Ntr, Nte = neigh(tr, 0), neigh(te, 0)
                Xtr1 = np.column_stack([Xtr0, Ptr])
                Xte1 = np.column_stack([Xte0, Pte])
                Xtr2 = np.column_stack([Xtr1, Ntr])
                Xte2 = np.column_stack([Xte1, Nte])
                ms = []
                for Xtr, Xte in [(Xtr0, Xte0), (Xtr1, Xte1), (Xtr2, Xte2)]:
                    beta, _, _, _ = np.linalg.lstsq(Xtr, ytr, rcond=None)
                    ms.append(float(((Xte @ beta - yte) ** 2).mean()))
                inc_rows.append({"program": gid, "cohort": cohort, "patient": p_held,
                                 "dM1": round(100 * (ms[0] - ms[1]) / ms[0], 2) if ms[0] > 0 else "NA",
                                 "dM2": round(100 * (ms[1] - ms[2]) / ms[1], 2) if ms[1] > 0 else "NA"})
        print(f"[{cohort}] done ({time.time()-t0:.0f}s)", flush=True)

    with open(OUT / f"external_replication_{args.tag}.tsv", "w") as f:
        f.write("program\tcohort\tpatient\tn_sections\tmol_in_ro\tspat_in\tspat_ro\n")
        for r in rep_rows:
            f.write(f"{r['program']}\t{r['cohort']}\t{r['patient']}\t{r['n_sections']}\t"
                    f"{r['mol']}\t{r['sp_in']}\t{r['sp_ro']}\n")
    with open(OUT / f"incremental_prediction_{args.tag}.tsv", "w") as f:
        f.write("program\tcohort\tpatient\tdM1_pct\tdM2_pct\n")
        for r in inc_rows:
            f.write(f"{r['program']}\t{r['cohort']}\t{r['patient']}\t{r['dM1']}\t{r['dM2']}\n")
    print(f"GO-RUN {args.tag} done ({time.time()-t0:.0f}s): "
          f"rep={len(rep_rows)} inc={len(inc_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
