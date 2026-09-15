#!/usr/bin/env python3
"""Apply frozen Tier-1 axes independently to a validation lineage (D-114 C).

No pooled clustering, no Vanderbilt HVG-10k filter. Same marker modules, genes
present in each section. Exploratory replication — does not write discovery
registry rows.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from r04.structure_readout import _auc
from r16 import axes as A
from r16 import census as C
from r16.section_io import load_section_symbols, usz_tls_labels

ROOT = Path(__file__).resolve().parent.parent
PC_BINS = np.arange(-4, 9)


def _axis_defs(proxy: dict) -> dict[str, list[str]]:
    defs = {k: list(v["voted_genes"]) for k, v in proxy["classes"].items()}
    defs["Plasma"] = list(A.PLASMA_GENES)
    return defs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--marker-proxy", type=Path, default=ROOT / "infra/r04/marker_proxy_combined.json")
    ap.add_argument("--null-draws", type=int, default=C.NULL_DRAWS)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    man = json.load(open(args.manifest))
    rows = man["rows"]
    axis_defs = _axis_defs(json.load(open(args.marker_proxy)))
    rng = np.random.default_rng(C.SEED)

    coh = {a: [] for a in axis_defs}
    stab = {a: [] for a in axis_defs}
    n_genes_used = {a: [] for a in axis_defs}
    pat_profiles: dict[str, np.ndarray] = {}
    tls_aucs: list[dict] = []
    section_fail: list[dict] = []

    for row in rows:
        try:
            sec = load_section_symbols(row)
        except Exception as exc:  # disclosed per-section fail-closed
            section_fail.append({"section_id": row["section_id"], "error": f"{type(exc).__name__}: {exc}"})
            continue
        needed = {g for gs in axis_defs.values() for g in gs} | {"CXCL13"}
        first: dict[str, int] = {}
        keep: list[int] = []
        genes: list[str] = []
        for i, g in enumerate(sec.gene_id):
            if g not in needed or g in first:
                continue
            first[g] = i
            keep.append(i)
            genes.append(g)
        counts = sparse_or_csr(sec.counts)[:, keep]
        gidx = {g: i for i, g in enumerate(genes)}
        x_z = C.log1p_zscore(counts)
        w = C.spatial_weights(np.asarray(sec.coords), C.KNN_SPATIAL)
        scores = {}
        for a, gs in axis_defs.items():
            s, n_used = A.axis_scores(x_z, gidx, gs)
            scores[a] = s
            n_genes_used[a].append(n_used)
            I, p = C.moran_permutation_p(s, w, args.null_draws, rng)
            coh[a].append({"patient_id": sec.patient_id, "section_id": sec.section_id,
                           "moran_i": I, "moran_p": p, "n_genes": n_used})
            stab[a].append(A.contour_stability(s)[0])
        if "CXCL13" in gidx:
            bmask = A.contour_mask(scores["B"], 0.7)
            if 0 < int(bmask.sum()) < len(bmask):
                signed = A.signed_hop_distance(w, bmask)
                prof = A.binned_profile(x_z[:, gidx["CXCL13"]], signed, PC_BINS)
                acc = pat_profiles.setdefault(sec.patient_id, np.zeros((len(PC_BINS), 2)))
                for i, (v, n) in enumerate(prof):
                    if n and np.isfinite(v):
                        acc[i, 0] += v * n
                        acc[i, 1] += n
        if sec.lineage == "TLS_VISIUM_USZ":
            y = usz_tls_labels(sec.section_id, sec.barcode, ROOT)
            if y is not None and np.isfinite(y).sum() >= 20 and len(np.unique(y[np.isfinite(y)])) == 2:
                fin = np.isfinite(y)
                yw = y[fin]
                wt = np.ones(int(fin.sum()))
                rec = {"section_id": sec.section_id, "patient_id": sec.patient_id,
                       "n_labeled": int(fin.sum()), "n_tls": int(yw.sum())}
                for a in ("B", "T", "Plasma"):
                    rec[f"auc_{a}"] = _auc(yw, scores[a][fin], wt)
                tls_aucs.append(rec)

    axes_summary = {}
    for a in axis_defs:
        per_pat: dict[str, list[float]] = {}
        for r in coh[a]:
            per_pat.setdefault(r["patient_id"], []).append(r["moran_p"])
        n_coh = sum(1 for v in per_pat.values() if float(np.median(v)) <= 0.01)
        axes_summary[a] = {
            "n_genes_median": int(np.median(n_genes_used[a])) if n_genes_used[a] else 0,
            "patients_coherent": n_coh,
            "n_patients": len(per_pat),
            "moran_i_median": float(np.median([r["moran_i"] for r in coh[a]])) if coh[a] else None,
            "moran_p_median": float(np.median([r["moran_p"] for r in coh[a]])) if coh[a] else None,
            "contour_stability_mean": float(np.mean(stab[a])) if stab[a] else None,
        }

    pc_peaks = {}
    for pat, acc in pat_profiles.items():
        prof = np.where(acc[:, 1] > 0, acc[:, 0] / np.maximum(acc[:, 1], 1), np.nan)
        if np.isfinite(prof).any():
            pc_peaks[pat] = int(PC_BINS[int(np.nanargmax(prof))])
    peaks = np.array(list(pc_peaks.values()), dtype=int) if pc_peaks else np.array([], dtype=int)
    positive_control = {
        "gene": "CXCL13", "axis": "B", "contour_quantile": 0.7,
        "n_patients": int(len(peaks)),
        "peak_hop_distribution": {str(int(k)): int(v) for k, v in
                                  zip(*np.unique(peaks, return_counts=True))} if len(peaks) else {},
        "n_peak_inside_or_edge_le1": int((peaks <= 1).sum()) if len(peaks) else 0,
    }

    artifact = {
        "schema": "r16.axis_replication.v1",
        "status": "EXPLORATORY_REPLICATION_NOT_CLAIM",
        "lineage": rows[0]["lineage"] if rows else None,
        "n_sections": len(rows),
        "n_sections_ok": len(rows) - len(section_fail),
        "section_errors": section_fail,
        "parameters": {
            "marker_filter": "section_present_voted_genes_no_vanderbilt_hvg_cap",
            "null_type": "value_permutation", "null_draws": args.null_draws,
            "pooled_clustering": False,
        },
        "axes": axes_summary,
        "positive_control": positive_control,
        "usz_tls_auc": tls_aucs,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, ensure_ascii=False))
    print(f"lineage={artifact['lineage']} ok={artifact['n_sections_ok']}/{artifact['n_sections']} fail={len(section_fail)}")
    for a, s in axes_summary.items():
        print(f"  {a:8s} coherent={s['patients_coherent']}/{s['n_patients']} "
              f"I={s['moran_i_median']:.2f} genes={s['n_genes_median']} stab={s['contour_stability_mean']:.2f}")
    pc = positive_control
    print(f"  PC CXCL13~B: n={pc['n_patients']} peak<=+1: {pc['n_peak_inside_or_edge_le1']}")
    if tls_aucs:
        b = [r["auc_B"] for r in tls_aucs if r.get("auc_B") is not None]
        print(f"  USZ TLS AUC B median={float(np.median(b)):.3f} n={len(b)}")
    return 0


def sparse_or_csr(mat):
    from scipy import sparse
    return mat.tocsr() if sparse.issparse(mat) else sparse.csr_matrix(mat)


if __name__ == "__main__":
    raise SystemExit(main())
