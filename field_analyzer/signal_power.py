"""FA-07B: how much observation is needed for the displacement to be identifiable?

FA-03B showed that B1 recovers a 7-pitch displacement from the *latent* field but not
from a 12-marker negative-binomial observation.  This module answers the practical
follow-up: how many markers (and how much depth) are required before the observed
score carries the displacement?

Sweep (same S2 latent fields, same operator, only the observation changes):
  markers per field in {12, 36, 96}
  depth multiplier   in {1, 3}
Success criterion: the B1 weight peak lands within 2 pitch of the true shift (7).
"""

from __future__ import annotations

from typing import Any, Iterable

import numpy as np
from scipy import sparse

from . import route_b, synthetic
from .depth_conditioning import residualise
from .pipeline import InputSpec, build_inputs


def score_fields(sample: synthetic.SyntheticSample, *, depth_mu: float, seed: int) -> dict[str, np.ndarray]:
    """Build counts at the requested observation power and score every field."""
    counts = synthetic.counts_from_fields(sample, depth_mu=depth_mu, seed=seed)
    genes = np.array(sample.gene_names())
    xy = synthetic.grid_xy(synthetic.GRID)
    support = np.ones(len(xy), dtype=bool)
    spec = InputSpec(
        fields=sample.field_gene_sets(),
        witness_genes=frozenset(),
        k_neighbors=6,
        min_library_genes=50,
    )
    inputs = build_inputs(sparse.csr_matrix(counts), genes, xy, support, support, spec)
    return {fid: inputs.arrays[f"field::{fid}"] for fid in spec.fields}


def sweep(
    seed: int = 0,
    *,
    markers: Iterable[int] = (12, 36, 96),
    depths: Iterable[float] = (1.0, 3.0),
    sigmas: Iterable[float] = (9.0,),
    max_lag: int = 12,
    step: int = 4,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for sigma in sigmas:
      for m in markers:
        sample = synthetic.build_s2(seed, sigma=sigma)
        sample.markers_per_field = int(m)
        truth_shift = float(np.hypot(*sample.structures["requested_offset"]))
        for depth_mu in depths:
            scores = score_fields(sample, depth_mu=depth_mu, seed=seed)
            n = synthetic.GRID
            anchor = scores["anchor"]
            for partner_id, expected in (("partner_same_point", 0.0), ("partner_offset", truth_shift)):
                fit = route_b.fit_lag_operator(
                    anchor, scores[partner_id], n, max_lag=max_lag, step=step
                )
                cv_b0 = route_b.blocked_cv(
                    anchor, scores[partner_id], n,
                    lambda train_mask: route_b.fit_isotropic_operator(
                        anchor, scores[partner_id], n, max_lag=max_lag, train_mask=train_mask
                    ),
                )
                cv_b1 = route_b.blocked_cv(
                    anchor, scores[partner_id], n,
                    lambda train_mask: route_b.fit_lag_operator(
                        anchor, scores[partner_id], n, max_lag=max_lag, step=step, train_mask=train_mask
                    ),
                )
                peak = fit.weight_peak
                peak_err = float(np.hypot(peak[0] - (0 if expected == 0 else 0), peak[1] - expected))
                cc = route_b.cross_correlation_peak(anchor, scores[partner_id], n, max_lag=max_lag)
                cc_err = float(np.hypot(0.0, cc["magnitude"] - expected))
                rows.append(
                    {
                        "seed": seed,
                        "sigma_pitch": float(sigma),
                        "markers_per_field": int(m),
                        "depth_mu": float(depth_mu),
                        "partner": partner_id,
                        "expected_shift_pitch": expected,
                        "peak_lag": list(peak),
                        "peak_error_pitch": peak_err,
                        "recovered": bool(peak_err <= 2.0),
                        "B1_disp": float(np.hypot(*fit.displacement)),
                        "B0_cv_sse": cv_b0,
                        "B1_cv_sse": cv_b1,
                        "cv_ratio": (cv_b1 / cv_b0) if cv_b0 > 0 else float("nan"),
                        "cc_corr": cc["corr"],
                        "cc_dy": cc["dy"],
                        "cc_dx": cc["dx"],
                        "cc_magnitude": cc["magnitude"],
                        "cc_error_pitch": cc_err,
                        "recovered_by_cc": bool(cc_err <= 2.0),
                    }
                )
    return rows
