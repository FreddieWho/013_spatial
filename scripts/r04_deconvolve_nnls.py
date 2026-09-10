#!/usr/bin/env python3
"""NNLS deconvolution against the frozen CRC reference (exploratory).

For each spot, solves non-negative least squares of the log1p count vector
onto the six Major-class mean profiles from the reference bundle, then
row-normalizes to pseudo-proportions. Transparent limitations (recorded in
output, not hidden): linear mixing assumption, no platform-effect modeling
between scRNA reference and Visium spots, mean-profile collapse of
within-class heterogeneity, no uncertainty propagation. This is a baseline
against which fancier methods can be compared, not a measurement.

Inputs are hash-pinned (reference provenance.json) and outputs carry full
provenance. Exploratory grade only.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import optimize, sparse

MAJORS = ("T", "B", "Mye", "ILC", "Epi", "Stromal")


def _read_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def load_reference(ref_dir: Path) -> dict:
    d = np.load(ref_dir / "reference_counts.npz", allow_pickle=True)
    prov = _read_json(ref_dir / "provenance.json")
    if prov.get("schema") != "r04.composition_reference.v1":
        raise ValueError("reference provenance schema mismatch")
    mat = sparse.csr_matrix(
        (d["counts_data"], d["counts_indices"], d["counts_indptr"]),
        shape=tuple(d["counts_shape"]))
    return {
        "matrix": mat,  # genes x cells
        "genes": [str(v) for v in list(d["panel_symbols"])],
        "major": np.asarray(d["major"].tolist()),
        "provenance": prov,
    }


def differential_markers(ref: dict, top_n: int) -> list[str]:
    """Top-N Welch t-stat genes per Major class on log1p counts."""
    from scipy import stats
    mat = ref["matrix"].tocsr()
    log_data = mat.copy()
    log_data.data = np.log1p(log_data.data)
    majors = ref["major"]
    genes = ref["genes"]
    picked: dict[str, None] = {}
    for mj in MAJORS:
        idx = np.flatnonzero(majors == mj)
        rest = np.flatnonzero(majors != mj)
        if len(idx) < 10 or len(rest) < 10:
            continue  # too few cells to estimate a contrast; never pick blind
        take = min(20000, len(idx), len(rest))
        rng = np.random.default_rng(20260910)
        a = log_data[:, rng.choice(idx, take, replace=False)].toarray()
        b = log_data[:, rng.choice(rest, take, replace=False)].toarray()
        with np.errstate(all="ignore"):
            t_stat, _ = stats.ttest_ind(a, b, axis=1, equal_var=False)
        t_stat = np.nan_to_num(t_stat, nan=-np.inf)
        mean_diff = a.mean(axis=1) - b.mean(axis=1)
        order = np.lexsort((-mean_diff, -t_stat))
        kept = 0
        for j in order:
            if mean_diff[j] <= 0:
                continue
            picked.setdefault(genes[int(j)], None)
            kept += 1
            if kept >= top_n:
                break
    return list(picked)


def mean_profiles(ref: dict, space: str = "proportion") -> tuple[np.ndarray, list[str]]:
    """Mean profile per Major class (genes x classes).

    proportion space (default): mean counts L1-normalized per class and per
    spot before NNLS. Removes total-RNA-content bias that otherwise lets
    high-RNA classes (epithelial) absorb immune signal. logmean space keeps
    the legacy behavior (mean of log1p) for comparison.
    """
    mat = ref["matrix"].tocsr()
    majors = ref["major"]
    if space == "proportion":
        col_sums = np.asarray(mat.sum(axis=0)).ravel()
        col_sums[col_sums == 0] = 1.0
        normed = mat.multiply(1.0 / col_sums).tocsr()  # multiply may return coo
        profiles = []
        for mj in MAJORS:
            rows = np.flatnonzero(majors == mj)
            if len(rows) == 0:
                raise ValueError(f"reference has no cells for {mj}")
            prof = np.asarray(normed[:, rows].mean(axis=1)).ravel()
            profiles.append(prof / prof.sum())
        return np.column_stack(profiles), list(MAJORS)
    if space != "logmean":
        raise ValueError(f"unknown profile space: {space}")
    log_data = mat.copy()
    log_data.data = np.log1p(log_data.data)
    profiles = []
    for mj in MAJORS:
        rows = np.flatnonzero(majors == mj)
        if len(rows) == 0:
            raise ValueError(f"reference has no cells for {mj}")
        profiles.append(np.asarray(log_data[:, rows].mean(axis=1)).ravel())
    return np.column_stack(profiles), list(MAJORS)


def deconvolve_nnls(spot_matrix, profiles: np.ndarray, space: str = "proportion") -> np.ndarray:
    """Row-normalized NNLS weights (spots x classes)."""
    if sparse.issparse(spot_matrix):
        spot_matrix = spot_matrix.toarray()
    spot_matrix = np.asarray(spot_matrix, dtype=float)
    if space == "proportion":
        row_sums = spot_matrix.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        targets = spot_matrix / row_sums
    else:
        targets = np.log1p(np.maximum(spot_matrix, 0.0))
    out = np.zeros((targets.shape[0], profiles.shape[1]))
    for i in range(targets.shape[0]):
        weights, _ = optimize.nnls(profiles, targets[i])
        total = weights.sum()
        out[i] = weights / total if total > 0 else np.full_like(weights, 1.0 / len(weights))
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref-dir", type=Path, required=True)
    parser.add_argument("--counts", type=Path, required=True,
                        help="spot x gene count matrix (.npz with data/indices/indptr/shape)")
    parser.add_argument("--genes", type=Path, required=True,
                        help="gene symbols in count-matrix column order, one per line")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-shared-genes", type=int, default=100)
    parser.add_argument("--profile-space", choices=["proportion", "logmean"],
                        default="proportion")
    parser.add_argument("--de-genes-per-class", type=int, default=0,
                        help="restrict to top-N Welch t-stat genes per class (0 = all shared)")
    args = parser.parse_args()
    ref = load_reference(args.ref_dir)
    arc = np.load(args.counts, allow_pickle=False)
    spots = sparse.csr_matrix(
        (arc["data"], arc["indices"], arc["indptr"]), shape=tuple(arc["shape"]))
    raw_lines = [l.strip() for l in Path(args.genes).read_text().splitlines()]
    if len(raw_lines) != spots.shape[1]:
        raise SystemExit("gene list length differs from count matrix width")
    # Positional mapping: empty entries (unmapped panel genes) never match;
    # duplicated symbols keep the first column (recorded below).
    first_pos: dict[str, int] = {}
    dup_symbols = 0
    for i, g in enumerate(raw_lines):
        if not g:
            continue
        if g in first_pos:
            dup_symbols += 1
            continue
        first_pos[g] = i
    genes = first_pos  # symbol -> column; placeholder replaced below
    ref_genes = ref["genes"]
    ref_set = set(ref_genes)
    common = [g for g in genes if g in ref_set]
    if len(common) < args.min_shared_genes:
        raise SystemExit(f"only {len(common)} shared genes; refusing")
    spot_idx = [genes[g] for g in common]
    ref_pos = {g: i for i, g in enumerate(ref_genes)}
    ref_idx = [ref_pos[g] for g in common]
    print(f"shared genes: {len(common)} (dup symbols collapsed: {dup_symbols})", flush=True)
    if args.de_genes_per_class > 0:
        keep = differential_markers(ref, args.de_genes_per_class)
        keep_idx = [i for i, g in enumerate(ref["genes"]) if g in set(keep)]
        if len(keep_idx) < 100:
            raise SystemExit("DE marker selection returned too few genes")
        ref = dict(ref, matrix=ref["matrix"][keep_idx, :].tocsr(),
                   genes=[ref["genes"][i] for i in keep_idx])
        # Re-derive intersections: pre-restriction indices are stale.
        ref_set = set(ref["genes"])
        common = [g for g in genes if g in ref_set]
        ref_pos = {g: i for i, g in enumerate(ref["genes"])}
        ref_idx = [ref_pos[g] for g in common]
        spot_idx = [genes[g] for g in common]
        print(f"DE restriction: {len(keep_idx)} genes, {len(common)} in spots", flush=True)
    profiles_full, _ = mean_profiles(ref, space=args.profile_space)
    profiles = profiles_full[ref_idx]
    props = deconvolve_nnls(spots[:, spot_idx], profiles, space=args.profile_space)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, proportions=props,
                        classes=np.asarray(list(MAJORS)),
                        genes_used=np.asarray(common))
    (args.output.with_suffix(".json")).write_text(json.dumps({
        "schema": "r04.deconv_nnls.v1",
        "evidence_grade": "exploratory_post_hoc_no_confirmatory_claim",
        "method": "per-spot NNLS onto Major mean profiles, row-normalized",
        "profile_space": args.profile_space,
        "limitations": [
            "linear mixing assumption",
            "no platform-effect modeling between scRNA reference and Visium spots",
            "mean-profile collapse of within-class heterogeneity",
            "no uncertainty propagation",
        ],
        "n_spots": int(props.shape[0]),
        "n_genes_shared": len(common),
        "reference_provenance": ref["provenance"],
    }, indent=1))
    means = props.mean(axis=0)
    print("mean composition:", {m: round(float(v), 3) for m, v in zip(MAJORS, means)})
    print("DECONV_OK", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
