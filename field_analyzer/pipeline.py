"""The single information entry shared by routes A/B/C.

Everything a model may see is produced here, from the *raw* matrix plus the
visible mask.  The function re-derives every statistic after masking, so any
leak introduced by "compute first, mask later" shows up as a poisoning failure.

Five entry points that must stay independent of hidden values
(FA_CONTRACT_ADDENDUM_20261020.md §1):
  1. library size / normalization denominator (feature genes only)
  2. per-spot depth and detection rate
  3. whole-slide mean/std of features
  4. expression-weighted adjacency
  5. query eligibility (geometry only)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np
from scipy import sparse

from . import geometry
from .observation import (
    ContractError,
    MaskViolation,
    visible_detection_rate,
    visible_library_size,
    visible_mean_std,
)


@dataclass(frozen=True)
class InputSpec:
    """Declared field definitions and the witness exclusion set."""

    fields: Mapping[str, Sequence[str]]
    witness_genes: frozenset[str] = frozenset()
    excluded_genes: frozenset[str] = frozenset()
    k_neighbors: int = 6
    max_edge_pitch: float = 3.0
    normalization_scale: float = 1e4
    mmd_scale_quantile: float = 0.5
    coordinate_unit: str = "source_pixels"
    external_library: bool = True
    min_library_genes: int = 100

    def feature_mask(self, genes: np.ndarray) -> np.ndarray:
        genes = np.asarray(genes).astype(str)
        blocked = self.witness_genes | self.excluded_genes
        mask = np.array([g not in blocked for g in genes], dtype=bool)
        leaked = sorted(self.witness_genes.intersection(genes))
        if leaked:
            # witness genes may exist in the matrix, they just must not be features
            pass
        return mask


@dataclass
class InputBundle:
    arrays: dict[str, np.ndarray]
    meta: dict[str, Any] = field(default_factory=dict)


def _field_gene_columns(genes: np.ndarray, members: Sequence[str], allowed: np.ndarray) -> np.ndarray:
    members = set(map(str, members))
    return np.array([(g in members) and ok for g, ok in zip(genes, allowed)], dtype=bool)


def build_inputs(
    counts: Any,
    genes: np.ndarray,
    xy: np.ndarray,
    tissue_support: np.ndarray,
    visible_mask: np.ndarray,
    spec: InputSpec,
) -> InputBundle:
    """Assemble every model-visible array. Hidden rows are NaN throughout."""
    counts = sparse.csr_matrix(counts)
    genes = np.asarray(genes).astype(str)
    xy = np.asarray(xy, dtype=np.float64)
    tissue_support = np.asarray(tissue_support, dtype=bool)
    visible_mask = np.asarray(visible_mask, dtype=bool)
    n_obs, n_genes = counts.shape
    if genes.shape != (n_genes,):
        raise ContractError("gene axis mismatch")
    if xy.shape != (n_obs, 2):
        raise ContractError("coordinate shape mismatch")

    # feature axis for the *library denominator*: everything except witnesses
    feature_mask = spec.feature_mask(genes)
    if not feature_mask.any():
        raise ContractError("no feature genes remain after witness exclusion")

    # scored genes: union of the declared fields (the only columns we densify)
    score_members = set().union(*(set(map(str, v)) for v in spec.fields.values())) if spec.fields else set()
    score_mask = np.array([(g in score_members) and ok for g, ok in zip(genes, feature_mask)], dtype=bool)
    if not score_mask.any():
        raise ContractError("no scored genes survive witness exclusion")
    score_genes = genes[score_mask]

    # (1) library size, visible spots only. In witness mode the denominator must
    # exclude every analysed field gene (A and B), not just the witness half.
    if spec.external_library:
        field_union = set(map(str, spec.witness_genes))
        for members in spec.fields.values():
            field_union.update(map(str, members))
        library_mask = np.array([g not in field_union for g in genes], dtype=bool)
        library_scope = "all_genes_excluding_all_analyzed_field_genes"
    else:
        library_mask = feature_mask
        library_scope = "all_genes_except_witness"
    if library_mask.sum() < spec.min_library_genes:
        raise ContractError(
            f"external library has only {int(library_mask.sum())} genes "
            f"(< min_library_genes={spec.min_library_genes})"
        )
    library = visible_library_size(counts, visible_mask, library_mask)
    # (2) depth and detection, visible spots only
    depth = library.copy()
    # detection uses the same target-independent denominator as the library:
    # using the feature axis would keep the field's own A genes in the ratio.
    detection = visible_detection_rate(counts, visible_mask, library_mask)

    sub = counts[:, score_mask].tocsr()
    dense = np.asarray(sub.todense(), dtype=np.float64)
    dense[~visible_mask] = np.nan
    # (3) whole-slide statistics of the scored genes, visible spots only
    slide_mean, slide_std = visible_mean_std(dense, visible_mask)

    with np.errstate(invalid="ignore", divide="ignore"):
        factor = np.where(
            np.isfinite(library), spec.normalization_scale / np.maximum(library, 1.0), np.nan
        )
        normalized = np.log1p(dense * factor[:, None])
    normalized[~visible_mask] = np.nan

    arrays: dict[str, np.ndarray] = {
        "library": library,
        "depth": depth,
        "detection": detection,
        "score_genes": score_genes,
        "slide_mean": slide_mean,
        "slide_std": slide_std,
    }
    scores: dict[str, np.ndarray] = {}
    for field_id in sorted(spec.fields):
        columns = _field_gene_columns(score_genes, spec.fields[field_id], score_mask[score_mask])
        score = np.full(n_obs, np.nan)
        if columns.any() and visible_mask.any():
            with np.errstate(invalid="ignore"):
                score[visible_mask] = np.nanmean(normalized[visible_mask][:, columns], axis=1)
        scores[field_id] = score
        arrays[f"field::{field_id}"] = scores[field_id]

    # (4) structural adjacency (geometry only) and expression-weighted adjacency
    support = tissue_support | visible_mask
    adjacency = geometry.build_adjacency(
        xy, support, k=spec.k_neighbors, max_edge_pitch=spec.max_edge_pitch, unit=spec.coordinate_unit
    )
    arrays["degree"] = adjacency.degrees
    arrays["isolated"] = adjacency.isolated
    arrays["pitch"] = np.array([adjacency.pitch])

    for field_id in sorted(scores):
        score = scores[field_id]
        weighted = _expression_weighted_graph(adjacency.graph, score, spec.mmd_scale_quantile)
        arrays[f"neighbor::{field_id}"] = geometry.neighbor_sum(
            weighted, np.nan_to_num(score, nan=0.0), valid=visible_mask
        )
        arrays[f"neighbor_struct::{field_id}"] = geometry.neighbor_sum(
            adjacency.graph, np.nan_to_num(score, nan=0.0), valid=visible_mask
        )

    # (5) query eligibility: geometric rule only, never value-driven
    eligible = tissue_support & np.isfinite(xy).all(axis=1)
    arrays["query_eligible"] = eligible

    meta = {
        "n_obs": int(n_obs),
        "n_visible": int(visible_mask.sum()),
        "n_genes": int(n_genes),
        "n_scored_genes": int(score_mask.sum()),
        "field_ids": sorted(spec.fields),
        "pitch": float(adjacency.pitch),
        "n_isolated": int(adjacency.n_isolated),
        "coordinate_unit": spec.coordinate_unit,
        "normalization_scale": float(spec.normalization_scale),
        "witness_excluded": len(spec.witness_genes),
        "k_neighbors": int(spec.k_neighbors),
        "max_edge_pitch": float(spec.max_edge_pitch),
        "n_edges": int(adjacency.notes["n_edges"]),
        "slide_stats_scope": "scored_genes",
        "library_scope": library_scope,
        "detection_scope": library_scope,
        "n_library_genes": int(library_mask.sum()),
    }
    return InputBundle(arrays=arrays, meta=meta)


def _expression_weighted_graph(graph: sparse.csr_matrix, score: np.ndarray, quantile: float) -> sparse.csr_matrix:
    """Gaussian similarity on the structural graph, visible values only."""
    s = np.asarray(score, dtype=np.float64)
    if not np.isfinite(s).any():
        return graph.copy()
    finite = s[np.isfinite(s)]
    scale = float(np.quantile(np.abs(finite - np.median(finite)), quantile)) if len(finite) > 1 else 0.0
    if not np.isfinite(scale) or scale <= 0:
        return graph.copy()
    coo = graph.tocoo()
    ok = np.isfinite(s[coo.row]) & np.isfinite(s[coo.col])
    diff = np.abs(s[coo.row] - s[coo.col])
    w = np.where(ok, np.exp(-diff / scale), 0.0)
    weighted = sparse.coo_matrix((w, (coo.row, coo.col)), shape=graph.shape).tocsr()
    degree = np.asarray(weighted.sum(axis=1)).ravel()
    return (sparse.diags(1.0 / np.maximum(degree, 1e-12)) @ weighted).tocsr()


def build_truth(
    counts: Any,
    genes: np.ndarray,
    hidden_rows: np.ndarray,
    spec: InputSpec,
) -> dict[str, np.ndarray]:
    """Evaluation-only truth. Never call this from a model path.

    Field truth uses the library computed over all spots *and the same score
    operator*, so the only difference from the model's input is information the
    model does not have.
    """
    counts = sparse.csr_matrix(counts)
    genes = np.asarray(genes).astype(str)
    feature_mask = spec.feature_mask(genes)
    score_members = set().union(*(set(map(str, v)) for v in spec.fields.values())) if spec.fields else set()
    score_mask = np.array([(g in score_members) and ok for g, ok in zip(genes, feature_mask)], dtype=bool)
    score_genes = genes[score_mask]
    all_rows = np.ones(counts.shape[0], dtype=bool)
    if spec.external_library:
        field_union = set(map(str, spec.witness_genes))
        for members in spec.fields.values():
            field_union.update(map(str, members))
        library_mask = np.array([g not in field_union for g in genes], dtype=bool)
    else:
        library_mask = feature_mask
    library = visible_library_size(counts, all_rows, library_mask)
    dense = np.asarray(counts[:, score_mask].todense(), dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        factor = np.where(np.isfinite(library), spec.normalization_scale / np.maximum(library, 1.0), np.nan)
        normalized = np.log1p(dense * factor[:, None])
    out: dict[str, np.ndarray] = {}
    for field_id in sorted(spec.fields):
        columns = _field_gene_columns(score_genes, spec.fields[field_id], score_mask[score_mask])
        if not columns.any():
            out[field_id] = np.full(len(hidden_rows), np.nan)
            continue
        with np.errstate(invalid="ignore"):
            out[field_id] = np.nanmean(normalized[hidden_rows][:, columns], axis=1)
    return out


def ridge_predict(inputs: InputBundle, field_id: str, alpha: float = 1.0) -> np.ndarray:
    """Deterministic reference model used to test prediction invariance.

    Visible-only closed-form ridge on (neighbour mean, depth, detection), then
    applied to every point.  No randomness, no hidden data.
    """
    y = inputs.arrays[f"field::{field_id}"]
    neighbor = inputs.arrays[f"neighbor_struct::{field_id}"]
    depth = inputs.arrays["depth"]
    detection = inputs.arrays["detection"]
    n = len(y)
    usable = np.isfinite(y) & np.isfinite(neighbor) & np.isfinite(depth) & np.isfinite(detection)
    if usable.sum() < 5:
        return np.full(n, np.nan)
    X = np.column_stack(
        [np.ones(usable.sum()), neighbor[usable], np.log1p(depth[usable]), detection[usable]]
    )
    target = y[usable]
    A = X.T @ X + alpha * np.eye(X.shape[1])
    A[0, 0] -= alpha  # do not penalise the intercept
    beta = np.linalg.solve(A, X.T @ target)

    X_all = np.column_stack(
        [
            np.ones(n),
            np.nan_to_num(neighbor, nan=np.nanmedian(neighbor[usable])),
            np.log1p(np.nan_to_num(depth, nan=np.nanmedian(depth[usable]))),
            np.nan_to_num(detection, nan=np.nanmedian(detection[usable])),
        ]
    )
    return X_all @ beta
