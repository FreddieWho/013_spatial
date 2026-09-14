#!/usr/bin/env python3
"""R-16 pilot: Seurat-style integration HVG (seurat_v3 + SelectIntegrationFeatures).

Per section: seurat_v3 HVG ranking on raw counts (loess log10(var)~log10(mean),
standardized variance). Combine: recurrence-first ranking — primary key = number
of PATIENTS where the gene is top-K in >=1 section (section votes collapse to
patient votes so multi-section patients don't dominate); tie-break = median
within-section rank. Fully deterministic.

Exploratory pilot to size the R-16 census universe (top 8000 / 10000) and to
check where canonical TLS markers land. No claim grade.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/r16/hvg_recurrence_pilot_20260914.json"

PER_SECTION_TOP_K = 5000
TLS_CANONICAL = [
    "CXCL13", "MS4A1", "CD3D", "CD3E", "CD19", "CCL19", "CCL21", "LTB",
    "CXCR5", "CCR7", "CD79A", "BANK1", "SELL", "CR2",
]


def seurat_v3_ranks(mean: np.ndarray, var: np.ndarray) -> np.ndarray:
    """Rank genes by seurat_v3 standardized variance. Returns rank array (0=best)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.log10(mean)
        y = np.log10(var)
    ok = np.isfinite(x) & np.isfinite(y) & (mean > 0) & (var > 0)
    xs, ys = x[ok], y[ok]
    fitted = lowess(ys, xs, frac=0.3, return_sorted=False)
    # spread around the fit
    resid = ys - fitted
    sd = np.sqrt(lowess(resid ** 2, xs, frac=0.3, return_sorted=False).clip(min=1e-12))
    std_var = np.full(len(mean), -np.inf)
    std_var[ok] = np.clip(resid / sd, -1e6, np.sqrt(len(mean)) * 10)
    order = np.argsort(-std_var)
    rank = np.empty(len(mean), dtype=np.int64)
    rank[order] = np.arange(len(mean))
    rank[~ok] = len(mean)  # unranked genes go last
    return rank


def main() -> int:
    man = json.load(open(ROOT / "infra/r04/role_manifests/training_manifest.json"))
    rows = man["rows"]
    union = json.load(open("/tmp/det_union_genes.json"))
    uidx = {s: i for i, s in enumerate(union)}
    U = len(union)

    patient_of = [r["patient_id"] for r in rows]
    patients = sorted(set(patient_of))
    pat_topk = {p: np.zeros(U, dtype=bool) for p in patients}  # gene top-K in >=1 section of patient
    rank_sum = np.zeros(U)
    rank_cnt = np.zeros(U)

    for r in rows:
        with h5py.File(r["matrix_locator"]) as f:
            X = f["X"]
            if isinstance(X, h5py.Group):
                from scipy.sparse import csr_matrix
                M = csr_matrix((X["data"][:], X["indices"][:], X["indptr"][:]),
                               shape=tuple(X.attrs["shape"]))
                mean = np.asarray(M.mean(axis=0)).ravel()
                var = np.asarray(M.multiply(M).mean(axis=0)).ravel() - mean ** 2
            else:
                A = X[:].astype(np.float64)
                mean = A.mean(axis=0)
                var = A.var(axis=0)
            syms = [s.decode() if isinstance(s, bytes) else s for s in f["var"]["_index"][:]]
        rk = seurat_v3_ranks(mean, var)
        idx = np.array([uidx[s] for s in syms])
        top_mask = rk < PER_SECTION_TOP_K
        pat_topk[r["patient_id"]][idx[top_mask]] = True
        rank_sum[idx] += rk
        rank_cnt[idx] += 1

    recurrence = np.zeros(U, dtype=int)
    for p in patients:
        recurrence += pat_topk[p].astype(int)
    median_rank = np.where(rank_cnt > 0, rank_sum / np.maximum(rank_cnt, 1), np.inf)

    order = sorted(range(U), key=lambda i: (-recurrence[i], median_rank[i]))
    final_rank = np.empty(U, dtype=int)
    for pos, i in enumerate(order):
        final_rank[i] = pos

    in_panel = set(l.strip() for l in open(ROOT / "infra/r04/composition_ref/cache_genes_symbols.txt") if l.strip())
    tls_report = {}
    for g in TLS_CANONICAL:
        i = uidx[g]
        tls_report[g] = {
            "patient_recurrence": int(recurrence[i]),
            "final_rank": int(final_rank[i]),
            "in_top8000": bool(final_rank[i] < 8000),
            "in_top10000": bool(final_rank[i] < 10000),
            "in_frozen4000": g in in_panel,
        }

    rec_dist = {str(k): int((recurrence == k).sum()) for k in range(0, max(recurrence) + 1)}
    genes_at = {n: int((final_rank < n).sum()) for n in [4000, 8000, 10000, 12000]}
    top_genes_8000 = sorted([union[i] for i in range(U) if final_rank[i] < 8000])
    overlap_4000 = len(set(top_genes_8000[:4000]) & in_panel)

    artifact = {
        "schema": "r16.hvg_recurrence_pilot.v1",
        "status": "EXPLORATORY_PILOT_NOT_CLAIM",
        "method": "seurat_v3 per-section top-K + patient-recurrence combine (SelectIntegrationFeatures analogue)",
        "per_section_top_k": PER_SECTION_TOP_K,
        "n_sections": len(rows),
        "n_patients": len(patients),
        "n_union_genes": U,
        "recurrence_distribution_patients": rec_dist,
        "genes_within_rank": genes_at,
        "overlap_top4000_vs_frozen4000": overlap_4000,
        "tls_canonical": tls_report,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=False))
    np.save(ROOT / "infra/r16/hvg_final_rank_20260914.npy", final_rank)

    print(f"patients={len(patients)} sections={len(rows)} union={U}")
    print(f"overlap top4000 vs frozen4000: {overlap_4000}")
    print(f"{'gene':8s} {'rec':>4s} {'rank':>6s} {'<=8k':>5s} {'<=10k':>6s}")
    for g in TLS_CANONICAL:
        r = tls_report[g]
        print(f"{g:8s} {r['patient_recurrence']:4d} {r['final_rank']:6d} "
              f"{str(r['in_top8000']):>5s} {str(r['in_top10000']):>6s}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
