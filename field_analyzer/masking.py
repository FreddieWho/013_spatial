"""Masked model-input pipeline and its five entry points (FA-00 §1).

The contract names five places where hidden expression can leak into a model even when the
expression matrix itself is masked:

1. library size / normalisation denominator (including the external reference denominator that
   must exclude the whole witness set B);
2. per-point depth and detection;
3. whole-section mean / variance / quantiles (the old ``feature_views`` broadcast these);
4. expression-weighted adjacency or graph weights;
5. query eligibility / ``ok_rows``.

Every entry point here is a separate function so the poisoning counterexample can fingerprint it
individually.  The rule is always *mask first, then aggregate*: rows outside
``visible_mask & tissue_support`` never enter a sum, a quantile, a graph weight, or a model fit.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Sequence

import numpy as np
from scipy import sparse

from .adapters import reference_intersection
from .geometry import TissueGeometry, build_adjacency
from .observation import (GeneRole, ObservationContractError, array_digest, json_digest,
                          text_digest)


@dataclass(frozen=True)
class InputConfig:
    """Everything the model is allowed to know, fixed before any expression is read."""

    feature_genes: tuple[str, ...]
    witness_genes: tuple[str, ...]
    reference_genes: tuple[str, ...]
    radius_um: float = 120.0
    k_neighbors: int = 6
    min_detection: float = 0.0
    quantiles: tuple[float, ...] = (0.05, 0.5, 0.95)
    target_field: int = 0
    ridge_lambda: float = 1.0
    coordinate_unit: str = 'um'
    scale_provenance: str = 'synthetic_um'
    tls_reference_path: str | None = None
    label_file_reference: str | None = None

    def provenance(self) -> dict:
        return dict(feature_genes=list(self.feature_genes),
                    witness_genes=list(self.witness_genes),
                    reference_genes=list(self.reference_genes),
                    radius_um=float(self.radius_um), k_neighbors=int(self.k_neighbors),
                    min_detection=float(self.min_detection),
                    quantiles=list(map(float, self.quantiles)),
                    target_field=int(self.target_field), ridge_lambda=float(self.ridge_lambda),
                    coordinate_unit=self.coordinate_unit,
                    scale_provenance=self.scale_provenance,
                    tls_reference_path=self.tls_reference_path,
                    label_file_reference=self.label_file_reference)


def _validate_counts(counts, genes, n_obs: int) -> sparse.csr_matrix:
    matrix = counts.tocsr() if sparse.issparse(counts) else sparse.csr_matrix(counts)
    if matrix.shape[0] != n_obs:
        raise ObservationContractError(
            f'counts have {matrix.shape[0]} rows but {n_obs} points were declared; '
            'refusing to guess the orientation')
    if matrix.shape[1] != len(genes):
        raise ObservationContractError(
            f'counts have {matrix.shape[1]} columns but {len(genes)} genes were declared')
    if matrix.nnz and (not np.all(np.isfinite(matrix.data)) or np.any(matrix.data < 0)):
        raise ObservationContractError('counts must be finite and non-negative')
    return matrix


def _visible_rows(visible_mask, tissue_support) -> np.ndarray:
    visible = np.asarray(visible_mask, dtype=bool)
    support = np.asarray(tissue_support, dtype=bool)
    if visible.shape != support.shape:
        raise ObservationContractError('visible_mask and tissue_support must have equal shape')
    if np.any(visible & ~support):
        raise ObservationContractError('visible points must lie inside tissue support')
    rows = np.flatnonzero(visible)
    if rows.size == 0:
        raise ObservationContractError('no visible points')
    return rows


def _gene_columns(genes: Sequence[str], wanted: Sequence[str], where: str) -> np.ndarray:
    index = {str(g): i for i, g in enumerate(genes)}
    missing = [g for g in wanted if str(g) not in index]
    if missing:
        raise ObservationContractError(f'{where}: {len(missing)} genes absent, e.g. {missing[:5]}')
    return np.array([index[str(g)] for g in wanted], dtype=np.int64)


def entry_library(counts, genes, rows, reference_genes) -> dict:
    """Entry 1: reference-only library size (never counts a witness gene)."""
    columns = _gene_columns(genes, reference_genes, 'reference_genes')
    ref = counts[rows][:, columns].tocsr()
    library = np.asarray(ref.sum(axis=1)).ravel().astype(np.float64)
    return dict(library=library, reference_genes=tuple(reference_genes),
                n_reference=len(reference_genes), n_zero_library=int((library <= 0).sum()))


def entry_depth_detection(counts, rows, reference_columns) -> dict:
    """Entry 2: per-point depth and detection from reference genes on visible rows only."""
    ref = counts[rows][:, reference_columns].tocsr()
    total = np.asarray(ref.sum(axis=1)).ravel().astype(np.float64)
    detected = np.asarray((ref > 0).sum(axis=1)).ravel().astype(np.float64)
    return dict(depth=np.log1p(total), total=total,
                detection=detected / max(len(reference_columns), 1))


def entry_section_stats(features_visible: np.ndarray, quantiles: Sequence[float]) -> dict:
    """Entry 3: whole-section statistics computed on visible rows only."""
    if features_visible.shape[0] == 0:
        raise ObservationContractError('no visible rows for section statistics')
    if not np.all(np.isfinite(features_visible)):
        raise ObservationContractError('section statistics may only see finite visible values')
    mean = features_visible.mean(axis=0)
    var = features_visible.var(axis=0)
    quant = np.quantile(features_visible, list(quantiles), axis=0)
    return dict(mean=mean, var=var, quantiles=quant, quantile_levels=tuple(map(float, quantiles)))


def entry_graph(geometry: TissueGeometry, values_visible: np.ndarray, cfg: InputConfig) -> dict:
    """Entry 4: expression-weighted adjacency built only from visible values."""
    graph = build_adjacency(geometry, radius=cfg.radius_um, k=cfg.k_neighbors,
                            values=values_visible, weight_mode='expression')
    return graph


def entry_query_eligibility(rows, library: np.ndarray, detection: np.ndarray,
                            min_detection: float) -> dict:
    """Entry 5: query list and eligibility mask; never reads hidden expression."""
    eligible = (library > 0) & (detection >= float(min_detection))
    return dict(query_rows=rows.astype(np.int64), eligible=eligible,
                ok_rows=rows[eligible].astype(np.int64), n_ok=int(eligible.sum()),
                n_dropped=int((~eligible).sum()))


def _standardise(train: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = train.mean(axis=0)
    scale = train.std(axis=0)
    scale = np.where(scale > 0, scale, 1.0)
    return (train - mean) / scale, mean, scale


def _ridge_ols(train: np.ndarray, target: np.ndarray, lam: float) -> tuple[np.ndarray, float]:
    design = np.column_stack([train, np.ones(train.shape[0])])
    gram = design.T @ design + lam * np.eye(design.shape[1])
    gram[-1, -1] -= lam  # do not penalise the intercept
    rhs = design.T @ target
    beta = np.linalg.solve(gram, rhs)
    return beta[:-1], float(beta[-1])


def build_model_inputs(counts, genes: Sequence[str], barcodes: Sequence[str], xy,
                       tissue_support, visible_mask, role: GeneRole, cfg: InputConfig) -> dict:
    """Full masked pipeline; returns arrays and per-entry fingerprints.

    ``counts`` is the *raw* matrix (visible plus hidden).  Rows outside ``visible_mask`` are
    dropped before any statistic is formed, so poisoning a hidden entry cannot change an
    observable quantity.
    """
    n_obs = len(np.asarray(barcodes))
    matrix = _validate_counts(counts, genes, n_obs)
    rows = _visible_rows(visible_mask, tissue_support)
    # --- role enforcement (never a boolean switch)
    role.assert_features(cfg.feature_genes, 'feature_genes')
    role.assert_no_witness(cfg.reference_genes, 'reference_genes')
    if set(cfg.feature_genes) != set(role.a_feature) & set(cfg.feature_genes):
        raise ObservationContractError('feature_genes must be A_FEATURE genes')
    overlap = reference_intersection(cfg.reference_genes, role)
    if overlap != 0:
        raise ObservationContractError(
            f'reference genes intersect A or B by {overlap}; the reference set must be disjoint')
    if set(cfg.witness_genes) != set(role.b_witness) & set(cfg.witness_genes):
        raise ObservationContractError('witness_genes must come from the B_WITNESS role')
    # --- entry 1 and 2
    library_info = entry_library(matrix, genes, rows, cfg.reference_genes)
    reference_columns = _gene_columns(genes, cfg.reference_genes, 'reference_genes')
    depth_info = entry_depth_detection(matrix, rows, reference_columns)
    # --- visible-only feature matrix
    feature_columns = _gene_columns(genes, cfg.feature_genes, 'feature_genes')
    counts_rows = matrix[rows][:, feature_columns].tocsr()
    library = library_info['library']
    scale = 1e4 / np.maximum(library, 1.0)
    normalized = np.log1p(np.asarray(counts_rows.todense(), dtype=np.float64) * scale[:, None])
    normalized[library <= 0] = np.nan
    features_visible = normalized
    # --- entry 3
    finite_rows = np.all(np.isfinite(features_visible), axis=1)
    stats_info = entry_section_stats(features_visible[finite_rows], cfg.quantiles)
    # --- entry 4
    field_value = np.full(n_obs, np.nan)
    field_value[rows[finite_rows]] = features_visible[finite_rows].mean(axis=1)
    geometry = TissueGeometry(xy=np.asarray(xy, dtype=np.float64),
                              coordinate_unit=cfg.coordinate_unit,
                              scale_provenance=cfg.scale_provenance,
                              tissue_support=np.asarray(tissue_support, dtype=bool),
                              visible_mask=np.asarray(visible_mask, dtype=bool))
    values_for_graph = np.array(field_value, copy=True)
    if not np.all(np.isfinite(values_for_graph[np.asarray(visible_mask, dtype=bool)])):
        raise ObservationContractError('graph values must be finite on visible points')
    graph = entry_graph(geometry, values_for_graph, cfg)
    # --- entry 5
    eligible_info = entry_query_eligibility(rows, np.where(np.isfinite(features_visible[:, 0]),
                                                           library, 0.0),
                                            depth_info['detection'], cfg.min_detection)
    # --- deterministic model on eligible visible rows (no hidden input, no test leakage)
    ok_local = np.flatnonzero(eligible_info['eligible'])
    if ok_local.size < 3:
        raise ObservationContractError('fewer than three eligible query points')
    design = features_visible[ok_local]
    target = design[:, cfg.target_field]
    standardised, feat_mean, feat_scale = _standardise(design)
    beta, intercept = _ridge_ols(standardised, target, cfg.ridge_lambda)
    predictions = standardised @ beta + intercept
    # --- visible-only matrix in store form (NaN at hidden rows)
    values_visible = np.full((n_obs, len(cfg.feature_genes)), np.nan)
    values_visible[rows] = features_visible
    inputs = dict(
        counts_visible=sparse.csr_matrix(np.where(np.isfinite(normalized), normalized, 0.0)),
        normalized_visible=values_visible, visible_mask=np.asarray(visible_mask, dtype=bool),
        tissue_support=np.asarray(tissue_support, dtype=bool), barcodes=list(map(str, barcodes)),
        xy=np.asarray(xy, dtype=np.float64), library=library,
        depth=depth_info['depth'], dets=depth_info['detection'],
        section_mean=stats_info['mean'], section_var=stats_info['var'],
        section_quantiles=stats_info['quantiles'],
        adjacency=graph['adjacency'], adjacency_degree=graph['degree'],
        graph_no_neighbor=graph['no_neighbor'], ok_rows=eligible_info['ok_rows'],
        query_rows=eligible_info['query_rows'], model_beta=beta, intercept=intercept,
        predictions=predictions, feature_mean=feat_mean, feature_scale=feat_scale,
        eligible=eligible_info['eligible'],
    )
    inputs['fingerprints'] = fingerprint_inputs(inputs, role, cfg)
    inputs['role_counts'] = role.counts()
    inputs['reference_intersection'] = overlap
    inputs['n_visible_rows'] = int(rows.size)
    inputs['n_hidden_rows'] = int(n_obs - rows.size)
    inputs['graph_missing_marks'] = int(graph['n_missing'])
    return inputs


def fingerprint_inputs(inputs: Mapping, role: GeneRole, cfg: InputConfig) -> dict:
    """The explicit fingerprint list from the addendum §5."""
    adjacency = inputs['adjacency']
    return dict(
        values_visible=array_digest(np.asarray(inputs['normalized_visible'], dtype=np.float64)),
        visible_mask=array_digest(np.asarray(inputs['visible_mask'])),
        tissue_support=array_digest(np.asarray(inputs['tissue_support'])),
        barcodes=text_digest(list(inputs['barcodes'])),
        xy=array_digest(np.asarray(inputs['xy'], dtype=np.float64)),
        gene_role=role.digest(),
        input_provenance=json_digest(cfg.provenance()),
        feature_transform=json_digest(dict(
            feature_mean=array_digest(np.asarray(inputs['feature_mean'], dtype=np.float64)),
            feature_scale=array_digest(np.asarray(inputs['feature_scale'], dtype=np.float64)),
            library=array_digest(np.asarray(inputs['library'], dtype=np.float64)),
            depth=array_digest(np.asarray(inputs['depth'], dtype=np.float64)),
            detection=array_digest(np.asarray(inputs['dets'], dtype=np.float64)),
            section_mean=array_digest(np.asarray(inputs['section_mean'], dtype=np.float64)),
            section_var=array_digest(np.asarray(inputs['section_var'], dtype=np.float64)),
            section_quantiles=array_digest(np.asarray(inputs['section_quantiles'],
                                                      dtype=np.float64)))),
        adjacency=json_digest(dict(
            shape=list(adjacency.shape),
            data=array_digest(np.asarray(adjacency.data, dtype=np.float64)),
            indices=array_digest(np.asarray(adjacency.indices)),
            indptr=array_digest(np.asarray(adjacency.indptr)),
            degree=array_digest(np.asarray(inputs['adjacency_degree'])),
            no_neighbor=array_digest(np.asarray(inputs['graph_no_neighbor'])))),
        query=json_digest(dict(
            query_rows=array_digest(np.asarray(inputs['query_rows'])),
            eligible=array_digest(np.asarray(inputs['eligible'])),
            ok_rows=array_digest(np.asarray(inputs['ok_rows'])))),
        model=json_digest(dict(
            beta=array_digest(np.asarray(inputs['model_beta'], dtype=np.float64)),
            intercept=float(inputs['intercept']))),
        predictions=array_digest(np.asarray(inputs['predictions'], dtype=np.float64)),
    )


def fingerprint_digests(fingerprints: Mapping) -> dict:
    """Collapse the fingerprint tree into one hash per top-level object."""
    out = {}
    for key, value in fingerprints.items():
        if isinstance(value, dict) and 'sha256' in value:
            out[key] = value['sha256']
        elif isinstance(value, dict):
            out[key] = json_digest(value)['sha256']
        else:
            out[key] = json_digest(value)['sha256']
    return out


def compare_fingerprints(left: Mapping, right: Mapping) -> dict:
    """Element-wise comparison; returns per-object equality and the reason on a mismatch."""
    result = {}
    keys = sorted(set(left) | set(right))
    for key in keys:
        a, b = left.get(key), right.get(key)
        if a is None or b is None:
            result[key] = dict(equal=False, reason='missing')
            continue
        if isinstance(a, dict) and isinstance(b, dict):
            sub = compare_fingerprints(a, b)
            equal = all(v['equal'] for v in sub.values())
            result[key] = dict(equal=equal, detail=sub if not equal else 'all_subobjects_equal')
        else:
            result[key] = dict(equal=(a == b), left=a, right=b) if a != b else dict(equal=True)
    return result
