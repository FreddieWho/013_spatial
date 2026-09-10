#!/usr/bin/env python3
"""Unnamed-field reproducibility audit (D-109): rotation-invariant, label-free.

Compares the K=3 loading subspaces of two fitted training exports (fold 0 vs
fold 4) without TLS/TSB labels, without factor-index matching, and without
cross-patient coordinate alignment. Per-fold rank-r subspaces are defined
intrinsically by SVD of the gene loading matrix; similarity is measured by
canonical correlations (invariant to within-subspace rotation, sign and
factor permutation). A gene-row permutation null (fixed seed) calibrates
chance-level alignment.

Pre-committed verdict rule per candidate subspace (rank-1, rank-2):
  SURVIVES          min observed canonical corr > max(null q975, 0.5)
  DOES_NOT_SURVIVE  min observed canonical corr <= null q975
  NOT_IDENTIFIABLE  otherwise (above chance but weak alignment)
Rank-3 values are reported descriptive-only (third direction already closed
by D-106/D-107; no verdict is issued here).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from r04.runtime import atomic_json

ALIGN_FLOOR = 0.5
N_PERM = 200
PERM_SEED = 20260911


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_loading(path: Path) -> tuple[np.ndarray, list[str]]:
    with np.load(path, allow_pickle=False) as archive:
        loading = np.asarray(archive["evaluation_loading"], dtype=float)
        genes = [str(v) for v in list(archive["frozen_gene_ids"])]
    if loading.ndim != 2 or loading.shape[1] != 3 or loading.shape[0] != len(genes):
        raise ValueError(f"unexpected loading layout: {path}")
    if not np.isfinite(loading).all():
        raise ValueError(f"non-finite loadings: {path}")
    return loading, genes


def _subspace_basis(loading: np.ndarray, rank: int) -> np.ndarray:
    u, s, _ = np.linalg.svd(loading, full_matrices=False)
    if s[rank - 1] <= 0:
        raise ValueError("degenerate loading subspace")
    return u[:, :rank]


def _canonical_corrs(basis_a: np.ndarray, basis_b: np.ndarray) -> list[float]:
    if basis_a.shape[0] != basis_b.shape[0] or basis_a.shape[1] != basis_b.shape[1]:
        raise ValueError("subspace bases are incompatible")
    s = np.linalg.svd(basis_a.T @ basis_b, compute_uv=False)
    return [float(min(1.0, max(0.0, v))) for v in s]


def _verdict(min_obs: float, null_q975: float) -> str:
    if min_obs <= null_q975:
        return "DOES_NOT_SURVIVE"
    if min_obs > max(null_q975, ALIGN_FLOOR):
        return "SURVIVES"
    return "NOT_IDENTIFIABLE"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--export-a", type=Path, required=True)
    parser.add_argument("--export-b", type=Path, required=True)
    parser.add_argument("--labels", nargs=2, default=["fold0", "fold4"])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--n-perm", type=int, default=N_PERM)
    parser.add_argument("--seed", type=int, default=PERM_SEED)
    args = parser.parse_args()
    loading_a, genes_a = _load_loading(args.export_a)
    loading_b, genes_b = _load_loading(args.export_b)
    if genes_a != genes_b:
        raise SystemExit("frozen gene universes differ; subspaces not comparable")
    rng = np.random.default_rng(int(args.seed))
    result: dict[str, object] = {
        "schema": "r04.unnamed_field_audit.v1",
        "status": "UNNAMED_REPRODUCIBILITY_AUDIT",
        "evidence_grade": "exploratory_post_hoc_no_confirmatory_claim",
        "decisions": ["D-107", "D-109"],
        "exports": {args.labels[0]: str(args.export_a), args.labels[1]: str(args.export_b)},
        "export_sha256": {args.labels[0]: _sha256(args.export_a),
                          args.labels[1]: _sha256(args.export_b)},
        "n_genes": len(genes_a),
        "method": ("SVD-defined rank-r loading subspaces; canonical correlations; "
                   "gene-row permutation null"),
        "candidates": [],
    }
    for rank in (1, 2, 3):
        basis_a = _subspace_basis(loading_a, rank)
        basis_b = _subspace_basis(loading_b, rank)
        observed = _canonical_corrs(basis_a, basis_b)
        null_min = []
        for _ in range(int(args.n_perm)):
            perm = rng.permutation(len(genes_b))
            null_basis = _subspace_basis(loading_b[perm], rank)
            null_min.append(min(_canonical_corrs(basis_a, null_basis)))
        null_min = np.asarray(null_min)
        q975 = float(np.quantile(null_min, 0.975))
        entry: dict[str, object] = {
            "candidate": f"rank{min(rank, 2)}_subspace" if rank <= 2 else "rank3_descriptive",
            "rank": rank,
            "observed_canonical_corrs": observed,
            "min_observed": min(observed),
            "null_min_q975": q975,
            "null_min_max": float(null_min.max()),
            "n_perm": int(args.n_perm),
            "seed": int(args.seed),
        }
        entry["verdict"] = _verdict(min(observed), q975) if rank <= 2 else "DESCRIPTIVE_ONLY"
        result["candidates"].append(entry)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(output, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
