"""FA-13: conditional cross-field increment outside the target's own neighbourhood.

Question (docs/FA13_CONDITIONAL_RELATION_20261023.md): with the target field's own visible
neighbourhood already in the model (M0), does adding another field's visible neighbourhood (M1)
improve prediction at fixed query points Q?

Leak-safety rules implemented here (they are the reason this module exists instead of reusing
``section_scores``):

  * everything starts from raw counts and the observation mask; Q is ``in_tissue`` minus the
    fixed FA-12 hold-out mask (same section id + seed 0 + 10 %), so Q is identical to FA-12;
  * a field score is the mean over its gene half of ``log1p(1e4 * count / background library)``
    where the background library excludes **all analysed field genes**;
  * feature-side genes of a direction drop the witness union of that direction (a gene that is a
    witness for any analysed field never enters a feature axis);
  * depth coefficients, standardisation, ridge coefficients and every neighbour aggregation are
    fitted/applied on O = in_tissue minus Q only.  Q counts, depth and witness values never enter
    a feature, a fit or a prediction; they only appear in the offline truth used for evaluation;
  * neighbour aggregation is k nearest *observed* points, self always excluded, and a query point
    without a valid neighbour is reported rather than filled with a zero;
  * neither model searches alpha, k, direction, scales or field subsets.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

from . import roles
from .spatial_lag import holdout_mask

PANEL_PATH = 'infra/tls_multiscale_gobp_20261002/panel_controls.json'
FIELDS = {
    'EMT': 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
    'IFNG': 'HALLMARK_INTERFERON_GAMMA_RESPONSE',
    'OXPHOS': 'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
}
PRIMARY_RULE = ('EMT', 'IFNG')
DEPTH_CONTROL_RULE = ('IFNG', 'OXPHOS')
RULES = (PRIMARY_RULE, DEPTH_CONTROL_RULE)
K_NEIGHBOURS = 8
RIDGE_ALPHA = 1.0
Q_FRACTION = 0.1
Q_SEED = 0
#: pre-fixed descriptive stability rule (no magnitude threshold, no p-value claim)
STABLE_UNITS = 6
MIN_FIT_ROWS = 50
MIN_EVAL_ROWS = 10


def load_panel(path: str | Path = PANEL_PATH) -> dict[str, list[str]]:
    panel = json.loads(Path(path).read_text())
    return {short: list(panel[key]) for short, key in FIELDS.items()}


def _mean_log_normalised(counts: sparse.csr_matrix, columns: np.ndarray,
                         library: np.ndarray) -> tuple[np.ndarray, int]:
    """Mean of ``log1p(1e4 * count / library)`` over ``columns``; NaN where library <= 0."""
    out = np.full(counts.shape[0], np.nan)
    if columns.size == 0:
        return out, 0
    sub = counts[:, columns].astype(np.float64).tocsr()
    ok = library > 0
    if ok.any():
        denom = library[ok][:, None]
        dense = np.asarray(sub[ok].todense(), dtype=np.float64)
        out[ok] = np.log1p(1e4 * dense / denom).mean(axis=1)
    return out, int(columns.size)


@dataclass
class SectionDesign:
    """Per-section raw scores, the fixed query mask and the fitted-on-O calibration."""

    section_id: str
    unit_id: str
    xy: np.ndarray
    in_tissue: np.ndarray
    q_mask: np.ndarray
    log_lib: np.ndarray
    scores: dict[tuple[str, str], np.ndarray]
    gene_sets: dict[tuple[str, str], tuple[str, ...]]
    barcodes: np.ndarray | None = None
    meta: dict = field(default_factory=dict)

    @property
    def n_points(self) -> int:
        return int(self.xy.shape[0])

    @property
    def observed(self) -> np.ndarray:
        return np.asarray(self.in_tissue, dtype=bool) & ~np.asarray(self.q_mask, dtype=bool)


def build_section(section_id: str, unit_id: str, counts, genes: Sequence[str], xy,
                  in_tissue, panel: Mapping[str, Sequence[str]], *,
                  seed: int = Q_SEED, q_fraction: float = Q_FRACTION,
                  barcodes: Sequence[str] | None = None) -> SectionDesign:
    """Build a section design from raw counts under the FA-13 masking rules."""
    counts = sparse.csr_matrix(counts)
    genes = np.asarray(genes, dtype=object)
    if counts.shape[1] != len(genes):
        raise ValueError(f'counts have {counts.shape[1]} columns but {len(genes)} genes')
    xy = np.asarray(xy, dtype=np.float64)
    in_tissue = np.asarray(in_tissue, dtype=bool)
    if xy.shape != (counts.shape[0], 2):
        raise ValueError(f'xy shape {xy.shape} does not match {counts.shape[0]} rows')
    if in_tissue.shape != (counts.shape[0],):
        raise ValueError('in_tissue must be per row')
    splits = roles.build_field_splits(panel, genes, seed=seed)
    field_names = sorted(panel)
    a_union: set[str] = set()
    b_union: set[str] = set()
    for name in field_names:
        a_union.update(splits[name].a_genes)
        b_union.update(splits[name].b_genes)
    panel_union = set().union(*(set(map(str, panel[name])) for name in field_names))
    background = [g for g in genes if str(g) not in panel_union]
    bg_index = {str(g): i for i, g in enumerate(genes)}
    bg_columns = np.array([bg_index[g] for g in background], dtype=np.int64)
    library = np.asarray(counts[:, bg_columns].sum(axis=1)).ravel().astype(np.float64) \
        if bg_columns.size else np.zeros(counts.shape[0])
    log_lib = np.log1p(library)
    scores: dict[tuple[str, str], np.ndarray] = {}
    gene_sets: dict[tuple[str, str], tuple[str, ...]] = {}
    for name in field_names:
        split = splits[name]
        half_sets = {
            ('A', 'feature'): [g for g in split.a_genes if g not in b_union],
            ('A', 'truth'): list(split.a_genes),
            ('B', 'feature'): [g for g in split.b_genes if g not in a_union],
            ('B', 'truth'): list(split.b_genes),
        }
        for (half, role), gene_list in half_sets.items():
            columns = np.array([bg_index[g] for g in gene_list if g in bg_index], dtype=np.int64)
            values, n_used = _mean_log_normalised(counts, columns, library)
            scores[(name, half, role)] = values
            gene_sets[(name, half, role)] = tuple(g for g in gene_list if g in bg_index)
    q_mask = holdout_mask(counts.shape[0], q_fraction, seed, section_id) & in_tissue
    return SectionDesign(section_id=section_id, unit_id=unit_id, xy=xy, in_tissue=in_tissue,
                         q_mask=q_mask, log_lib=log_lib, scores=scores, gene_sets=gene_sets,
                         barcodes=None if barcodes is None else np.asarray(barcodes, dtype=object),
                         meta=dict(n_background=len(background), n_panel=counts.shape[1] - len(background),
                                   splits={name: dict(n_a=len(splits[name].a_genes),
                                                      n_b=len(splits[name].b_genes))
                                           for name in field_names},
                                   n_feature_genes={f'{k[0]}:{k[1]}': len(v)
                                                    for k, v in gene_sets.items()},
                                   q_points=int(q_mask.sum()), observed_points=int(
                                       (in_tissue & ~q_mask).sum())))


def calibrate(values: np.ndarray, log_lib: np.ndarray, fit_mask: np.ndarray, *,
              conditioned: bool = True) -> tuple[np.ndarray, dict]:
    """Residualise against depth and z-score using **fit rows only** (O rows)."""
    v = np.asarray(values, dtype=np.float64).copy()
    fit = np.asarray(fit_mask, dtype=bool) & np.isfinite(v)
    params = dict(n_fit=int(fit.sum()), conditioned=bool(conditioned))
    if fit.sum() < 10:
        return np.full_like(v, np.nan), params
    if conditioned:
        x = np.asarray(log_lib, dtype=np.float64)[fit]
        y = v[fit]
        xm, ym = float(x.mean()), float(y.mean())
        var = float(np.var(x))
        slope = float(np.cov(x, y, bias=True)[0, 1] / var) if var > 0 else 0.0
        intercept = ym - slope * xm
        params.update(depth_intercept=intercept, depth_slope=slope)
        v = v - (intercept + slope * np.asarray(log_lib, dtype=np.float64))
    fit = np.asarray(fit_mask, dtype=bool) & np.isfinite(v)
    mu, sd = float(v[fit].mean()), float(v[fit].std())
    if not np.isfinite(sd) or sd <= 0:
        return np.full_like(v, np.nan), params
    params.update(z_mean=mu, z_std=sd)
    out = (v - mu) / sd
    return out, params


def neighbour_feature(xy: np.ndarray, in_tissue: np.ndarray, observed: np.ndarray,
                      values: np.ndarray, *, k: int = K_NEIGHBOURS) -> tuple[np.ndarray, np.ndarray]:
    """Mean of the ``k`` nearest *observed* neighbours; the point itself is always excluded."""
    xy = np.asarray(xy, dtype=np.float64)
    inside = np.flatnonzero(np.asarray(in_tissue, dtype=bool))
    obs = np.asarray(observed, dtype=bool) & np.asarray(in_tissue, dtype=bool)
    values = np.asarray(values, dtype=np.float64)
    mean = np.full(xy.shape[0], np.nan)
    count = np.zeros(xy.shape[0], dtype=np.int64)
    if inside.size == 0 or obs.sum() < k:
        return mean, count
    tree = cKDTree(xy[inside])
    kk = min(int(k) + 1, inside.size)
    dist, idx = tree.query(xy, k=kk)
    if dist.ndim == 1:
        dist, idx = dist[:, None], idx[:, None]
    for i in range(xy.shape[0]):
        neigh = inside[idx[i]]
        keep = (neigh != i) & obs[neigh] & np.isfinite(values[neigh])
        neigh = neigh[keep]
        if neigh.size == 0:
            continue
        neigh = neigh[:k]
        mean[i] = float(np.mean(values[neigh]))
        count[i] = int(neigh.size)
    return mean, count


def fold_design(section: SectionDesign, target: str, predictor: str, *,
                feature_half: str = 'A', truth_half: str = 'B', conditioned: bool = True,
                k: int = K_NEIGHBOURS) -> dict:
    """Design matrices for one section under one direction (feature half -> truth half)."""
    obs = section.observed
    same_raw = section.scores[(target, feature_half, 'feature')]
    cross_raw = section.scores[(predictor, feature_half, 'feature')]
    truth_raw = section.scores[(target, truth_half, 'truth')]
    same_cal, same_par = calibrate(same_raw, section.log_lib, obs, conditioned=conditioned)
    cross_cal, cross_par = calibrate(cross_raw, section.log_lib, obs, conditioned=conditioned)
    truth_cal, truth_par = calibrate(truth_raw, section.log_lib, obs, conditioned=conditioned)
    same_nb, same_n = neighbour_feature(section.xy, section.in_tissue, obs, same_cal, k=k)
    cross_nb, cross_n = neighbour_feature(section.xy, section.in_tissue, obs, cross_cal, k=k)
    return dict(section_id=section.section_id, unit_id=section.unit_id,
                y=truth_cal, same=same_nb, cross=cross_nb, obs=obs, q=section.q_mask,
                n_same=same_n, n_cross=cross_n,
                params=dict(same=same_par, cross=cross_par, truth=truth_par))


def _ridge(beta_design: np.ndarray, y: np.ndarray, alpha: float) -> np.ndarray:
    design = np.column_stack([np.ones(beta_design.shape[0]), beta_design])
    gram = design.T @ design + alpha * np.eye(design.shape[1])
    gram[0, 0] -= alpha
    return np.linalg.solve(gram, design.T @ y)


def _r2(pred: np.ndarray, truth: np.ndarray) -> tuple[float, int]:
    ok = np.isfinite(pred) & np.isfinite(truth)
    if ok.sum() < MIN_EVAL_ROWS:
        return float('nan'), int(ok.sum())
    ss_res = float(np.sum((pred[ok] - truth[ok]) ** 2))
    ss_tot = float(np.sum((truth[ok] - np.mean(truth[ok])) ** 2))
    return (1.0 - ss_res / ss_tot) if ss_tot > 0 else float('nan'), int(ok.sum())


def _rmse(pred: np.ndarray, truth: np.ndarray) -> float:
    ok = np.isfinite(pred) & np.isfinite(truth)
    if ok.sum() == 0:
        return float('nan')
    return float(np.sqrt(np.mean((pred[ok] - truth[ok]) ** 2)))


def run_rule(sections: Sequence[SectionDesign], target: str, predictor: str, *,
             feature_half: str = 'A', truth_half: str = 'B', conditioned: bool = True,
             k: int = K_NEIGHBOURS, alpha: float = RIDGE_ALPHA,
             protocol: str = 'conditioned') -> tuple[list[dict], list[dict]]:
    """Leave-one-unit-out M0 vs M1 on the same query points.

    Returns ``(unit_rows, section_rows)``.  M0 uses only the target's own neighbour summary;
    M1 adds the predictor field's neighbour summary.  Both are fitted on the same rows and
    scored on the same Q points.
    """
    designs = [fold_design(s, target, predictor, feature_half=feature_half,
                           truth_half=truth_half, conditioned=conditioned, k=k)
               for s in sections]
    units = sorted({d['unit_id'] for d in designs})
    unit_rows: list[dict] = []
    section_rows: list[dict] = []
    for held in units:
        train = [d for d in designs if d['unit_id'] != held]
        test = [d for d in designs if d['unit_id'] == held]
        X0, X1, y = [], [], []
        for d in train:
            sel = d['obs'] & np.isfinite(d['y']) & np.isfinite(d['same']) & np.isfinite(d['cross'])
            if sel.any():
                X0.append(d['same'][sel])
                X1.append(np.column_stack([d['same'][sel], d['cross'][sel]]))
                y.append(d['y'][sel])
        if not y or sum(len(v) for v in y) < MIN_FIT_ROWS:
            unit_rows.append(dict(target=target, predictor=predictor, feature_half=feature_half,
                                  truth_half=truth_half, protocol=protocol, held_unit=held,
                                  n_fit=0, n_test=0, r2_m0=float('nan'), r2_m1=float('nan'),
                                  delta_r2=float('nan'), rmse_m0=float('nan'), rmse_m1=float('nan'),
                                  status='NOT_TESTABLE_FIT_ROWS'))
            continue
        X0 = np.concatenate(X0)
        X1 = np.concatenate(X1)
        y = np.concatenate(y)
        beta0 = _ridge(X0[:, None], y, alpha)
        beta1 = _ridge(X1, y, alpha)
        pred0, pred1, truth = [], [], []
        for d in test:
            sel = d['q'] & np.isfinite(d['same']) & np.isfinite(d['cross']) & np.isfinite(d['y'])
            if not sel.any():
                continue
            p0 = beta0[0] + beta0[1] * d['same'][sel]
            p1 = beta1[0] + beta1[1] * d['same'][sel] + beta1[2] * d['cross'][sel]
            pred0.append(p0)
            pred1.append(p1)
            truth.append(d['y'][sel])
            r2_0, n_sec = _r2(p0, d['y'][sel])
            r2_1, _ = _r2(p1, d['y'][sel])
            section_rows.append(dict(
                target=target, predictor=predictor, feature_half=feature_half,
                truth_half=truth_half, protocol=protocol, unit_id=held,
                section_id=d['section_id'], n_test=int(sel.sum()),
                r2_m0=r2_0, r2_m1=r2_1, delta_r2=(r2_1 - r2_0) if np.isfinite(r2_0)
                and np.isfinite(r2_1) else float('nan'),
                rmse_m0=_rmse(p0, d['y'][sel]), rmse_m1=_rmse(p1, d['y'][sel]),
                n_no_neighbour=int((d['q'] & ~np.isfinite(d['same'])).sum())))
        if not truth:
            unit_rows.append(dict(target=target, predictor=predictor, feature_half=feature_half,
                                  truth_half=truth_half, protocol=protocol, held_unit=held,
                                  n_fit=int(y.size), n_test=0, r2_m0=float('nan'),
                                  r2_m1=float('nan'), delta_r2=float('nan'),
                                  rmse_m0=float('nan'), rmse_m1=float('nan'),
                                  status='NOT_TESTABLE_NO_QUERY'))
            continue
        pred0 = np.concatenate(pred0)
        pred1 = np.concatenate(pred1)
        truth = np.concatenate(truth)
        r2_0, n_test = _r2(pred0, truth)
        r2_1, _ = _r2(pred1, truth)
        unit_rows.append(dict(
            target=target, predictor=predictor, feature_half=feature_half,
            truth_half=truth_half, protocol=protocol, held_unit=held,
            n_fit=int(y.size), n_test=int(n_test),
            r2_m0=r2_0, r2_m1=r2_1,
            delta_r2=(r2_1 - r2_0) if np.isfinite(r2_0) and np.isfinite(r2_1) else float('nan'),
            rmse_m0=_rmse(pred0, truth), rmse_m1=_rmse(pred1, truth),
            delta_rmse=(_rmse(pred1, truth) - _rmse(pred0, truth)),
            beta0=float(beta0[1]), beta1_same=float(beta1[1]), beta1_cross=float(beta1[2]),
            status='OK'))
    return unit_rows, section_rows


def run_single_feature(sections: Sequence[SectionDesign], target: str, predictor: str, *,
                       feature: str, feature_half: str = 'A', truth_half: str = 'B',
                       conditioned: bool = True, k: int = K_NEIGHBOURS,
                       alpha: float = RIDGE_ALPHA) -> list[dict]:
    """Reference model with a single neighbour feature (``same`` or ``cross`` only).

    Used only to show what each feature achieves alone; the decision comparison is M0 (same)
    versus M1 (same + cross) in :func:`run_rule`.
    """
    if feature not in ('same', 'cross'):
        raise ValueError("feature must be 'same' or 'cross'")
    designs = [fold_design(s, target, predictor, feature_half=feature_half,
                           truth_half=truth_half, conditioned=conditioned, k=k)
               for s in sections]
    rows: list[dict] = []
    for held in sorted({d['unit_id'] for d in designs}):
        train = [d for d in designs if d['unit_id'] != held]
        test = [d for d in designs if d['unit_id'] == held]
        xs, ys = [], []
        for d in train:
            sel = d['obs'] & np.isfinite(d['y']) & np.isfinite(d[feature])
            if sel.any():
                xs.append(d[feature][sel])
                ys.append(d['y'][sel])
        if not ys or sum(len(v) for v in ys) < MIN_FIT_ROWS:
            rows.append(dict(held_unit=held, feature=feature, n_fit=0, n_test=0, r2=float('nan'),
                             status='NOT_TESTABLE_FIT_ROWS'))
            continue
        x = np.concatenate(xs)
        y = np.concatenate(ys)
        beta = _ridge(x[:, None], y, alpha)
        pred, truth = [], []
        for d in test:
            sel = d['q'] & np.isfinite(d[feature]) & np.isfinite(d['y'])
            if sel.any():
                pred.append(beta[0] + beta[1] * d[feature][sel])
                truth.append(d['y'][sel])
        if not truth:
            rows.append(dict(held_unit=held, feature=feature, n_fit=int(y.size), n_test=0,
                             r2=float('nan'), status='NOT_TESTABLE_NO_QUERY'))
            continue
        r2, n = _r2(np.concatenate(pred), np.concatenate(truth))
        rows.append(dict(held_unit=held, feature=feature, n_fit=int(y.size), n_test=int(n),
                         r2=r2, status='OK'))
    return rows


def summarise(unit_rows: Sequence[Mapping], *, stable_units: int = STABLE_UNITS) -> dict:
    """Descriptive summary per rule/direction/protocol; no magnitude threshold."""
    out: dict[str, dict] = {}
    keys = sorted({(r['target'], r['predictor'], r['feature_half'], r['truth_half'],
                    r['protocol']) for r in unit_rows})
    for key in keys:
        rows = [r for r in unit_rows if (r['target'], r['predictor'], r['feature_half'],
                                         r['truth_half'], r['protocol']) == key and r['status'] == 'OK']
        n = len(rows)
        deltas = [float(r['delta_r2']) for r in rows if np.isfinite(r['delta_r2'])]
        r2m0 = [float(r['r2_m0']) for r in rows if np.isfinite(r['r2_m0'])]
        r2m1 = [float(r['r2_m1']) for r in rows if np.isfinite(r['r2_m1'])]
        positive = int(sum(d > 0 for d in deltas))
        label = f'{key[0]}<-{key[1]}|{key[2]}/{key[3]}|{key[4]}'
        out[label] = dict(
            target=key[0], predictor=key[1], feature_half=key[2], truth_half=key[3],
            protocol=key[4], n_units=n, n_positive=positive,
            median_delta_r2=float(np.median(deltas)) if deltas else float('nan'),
            mean_delta_r2=float(np.mean(deltas)) if deltas else float('nan'),
            median_r2_m0=float(np.median(r2m0)) if r2m0 else float('nan'),
            median_r2_m1=float(np.median(r2m1)) if r2m1 else float('nan'),
            unit_delta_r2={r['held_unit']: float(r['delta_r2']) for r in rows
                           if np.isfinite(r['delta_r2'])},
            stable=(n >= stable_units and positive >= stable_units),
            decision_rule=f'stable = >= {stable_units} of >= {stable_units} units with delta_r2 > 0 '
                          '(descriptive consistency, no magnitude threshold)',
            n_not_testable=int(sum(1 for r in unit_rows
                                   if (r['target'], r['predictor'], r['feature_half'],
                                       r['truth_half'], r['protocol']) == key
                                   and r['status'] != 'OK')))
    return out
