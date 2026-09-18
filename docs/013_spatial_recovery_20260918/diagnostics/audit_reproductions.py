#!/usr/bin/env python3
"""Small mathematical/software counterexamples, NOT a reanalysis of patient data.

These reproduce patterns observed in the reviewed repository. They do not
import repository modules and are not proposed production inference code.
Usage: python diagnostics/audit_reproductions.py --out evidence/toy_audit.json
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial import cKDTree
from scipy.spatial.distance import squareform


def groups(loadings: np.ndarray, *, sign_invariant: bool) -> np.ndarray:
    """Rows are axes, NOT directed biological cluster contrasts."""
    a = np.asarray(loadings, dtype=float)
    norms = np.linalg.norm(a, axis=1, keepdims=True)
    if a.ndim != 2 or len(a) < 2 or np.any(norms == 0):
        raise ValueError('Need at least two nonzero loading rows')
    a = a / norms
    sim = a @ a.T
    if sign_invariant:
        sim = np.abs(sim)
    dist = 1 - np.clip(sim, -1, 1)
    np.fill_diagonal(dist, 0)
    return fcluster(linkage(squareform(dist, checks=False), method='average'),
                    t=0.25, criterion='distance')


def contour_jaccard(x: np.ndarray) -> list[float]:
    masks = [x >= np.quantile(x, q) for q in (0.6, 0.7, 0.8)]
    return [float(np.logical_and(a, b).sum()/np.logical_or(a, b).sum())
            for a, b in zip(masks[:-1], masks[1:])]


def cohort_stat(profiles: np.ndarray) -> float:
    """Input: patient x bin; every patient has one vote."""
    return float(np.max(np.median(profiles, axis=0)))


def matched_null_stats(surrogates: np.ndarray) -> np.ndarray:
    """Input: joint draw x patient x bin, output: joint draw."""
    return np.max(np.median(surrogates, axis=1), axis=1)


def monte_carlo_p(observed: float, null_stats: np.ndarray) -> float:
    z = np.asarray(null_stats, dtype=float)
    if z.ndim != 1 or len(z) < 1 or not np.isfinite(z).all():
        raise ValueError('null_stats must be a finite nonempty vector')
    return float((1 + np.count_nonzero(z >= observed))/(len(z) + 1))


def run(seed: int = 20260918) -> dict:
    rng = np.random.default_rng(seed)
    v = np.array([1., -2., 3., 4.]); v /= np.linalg.norm(v)
    L = np.stack([v, -v])
    signed_groups = groups(L, sign_invariant=False)
    abs_groups = groups(L, sign_invariant=True)
    orient = np.where((L @ v) >= 0, 1., -1.)
    sign_case = dict(legacy_groups=int(len(set(signed_groups))),
                     sign_invariant_groups=int(len(set(abs_groups))),
                     legacy_mean_norm=float(np.linalg.norm(L.mean(axis=0))),
                     oriented_mean_norm=float(np.linalg.norm((L*orient[:,None]).mean(axis=0))))

    basis = np.eye(8)[:, :2]
    angle = np.pi/4
    rotation = np.array([[np.cos(angle), -np.sin(angle)],
                         [np.sin(angle), np.cos(angle)]])
    rotated = basis @ rotation
    rotation_case = dict(best_individual_axis_cosine=float(np.abs(basis.T@rotated).max()),
                         canonical_correlations=np.linalg.svd(basis.T@rotated, compute_uv=False).tolist())

    q, _ = np.linalg.qr(rng.normal(size=(600, 6)) - 0)  # then center/re-orthogonalize
    q, _ = np.linalg.qr(q-q.mean(axis=0))
    target = q.sum(axis=1)/np.sqrt(6)
    corrs = np.corrcoef(np.column_stack([q, target]), rowvar=False)[-1, :-1]
    residual = target-q@np.linalg.lstsq(q, target, rcond=None)[0]
    mixture_case = dict(max_individual_axis_correlation=float(np.abs(corrs).max()),
                        joint_axis_r_squared=float(1-(residual@residual)/(target@target)))

    ordered = np.arange(1200, dtype=float)
    shuffled = rng.permutation(ordered)
    j1, j2 = contour_jaccard(ordered), contour_jaccard(shuffled)
    contour_case = dict(ordered_jaccards=j1, shuffled_jaccards=j2,
                        mean_jaccard=float(np.mean(j1)), expected_by_algebra=17/24)

    # Entirely synthetic Gaussian patient profiles; not a tissue null model.
    null_profiles = rng.normal(size=(999, 30, 13))
    legacy = np.max(null_profiles, axis=2).ravel()  # WRONG cohort statistic
    correct = matched_null_stats(null_profiles)
    # A low, consistent positive shape distinguishes the distributions.
    shape = 0.9*np.exp(-0.5*((np.arange(13)-6)/2)**2)
    observed_profiles = shape[None, :] + 0.1*rng.normal(size=(30,13))
    stat = cohort_stat(observed_profiles)
    aggregation_case = dict(
        n_patients=30, n_joint_draws=999,
        observed_cohort_statistic=stat,
        legacy_null_count=int(len(legacy)), correct_null_count=int(len(correct)),
        legacy_null_median=float(np.median(legacy)), correct_null_median=float(np.median(correct)),
        legacy_p=monte_carlo_p(stat, legacy), matched_statistic_p=monte_carlo_p(stat, correct),
        warning='Illustration of statistic mismatch, NOT a calibrated spatial test or effect estimate.')

    # Repeat one patient's identical section: section-weighted output changes.
    pat_values = np.array([1., 0., -1.])
    section_repeated = np.array([1.]*8+[0., -1.])
    patient_case = dict(original_section_median=float(np.median(pat_values)),
                        repeated_section_median=float(np.median(section_repeated)),
                        patient_first_median=float(np.median(pat_values)))

    # Nearest-neighbor remapping can duplicate indices; no matching constraints.
    coords = np.array([[0.,0.],[1.,0.],[2.,0.],[3.,0.], [0.,1.],[1.,1.], [0.,2.]])
    moved = -(coords-coords.mean(axis=0))+coords.mean(axis=0)+np.array([0.1,-0.1])
    dist, idx = cKDTree(coords).query(moved, k=1)
    mapping_case = dict(n_spots=len(coords), n_unique_source_indices=int(len(np.unique(idx))),
                        indices=idx.tolist(), is_bijection=bool(len(np.unique(idx))==len(coords)),
                        note='Demonstrates missing bijection guarantee; not a measurement of real remap acceptance.')
    return dict(scope='SYNTHETIC_MATHEMATICAL_AUDIT_ONLY', seed=seed,
                patient_data_used=False, repository_training_run=False,
                pc_sign=sign_case, pc_rotation=rotation_case,
                known_axis_mixture=mixture_case, contour_stability=contour_case,
                laneA_statistic_mismatch=aggregation_case,
                patient_weighting=patient_case, nearest_neighbor_remap=mapping_case)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--seed', type=int, default=20260918)
    args = p.parse_args()
    result = run(args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
