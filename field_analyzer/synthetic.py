"""Synthetic kill tests (FA-02).

Purpose: decide whether a route has a *capability* before touching real data.
Two independent generators:

  G1  analytic / spline construction (geometric truth is explicit)
  G2  a different mechanism (anisotropic diffusion on the grid graph, so a route
      cannot win by being the same model that generated the data)

Scenarios follow 06_SYNTHETIC_KILL_TESTS.md:
  S0 null structure, S1 non-uniform shape, S2 non-trivial coupling,
  S3 connectivity, S4 distractor inflation, S5 non-identifiable gap.

Every sample keeps three things physically apart:
  observed_inputs/    what a model may read (counts on visible spots only)
  latent_truth/       field values used for scoring (never read by a model)
  evaluation_only/    masks, structure parameters, pair answers

Nothing here reads real data, human labels, or the network.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

PITCH = 1.0
GRID = 64
SEED_ROLES = {
    "train": tuple(range(0, 8)),
    "development": (8, 9, 10, 11),
    "evaluation": tuple(range(12, 28)),
}
MARKERS_PER_FIELD = 12
BACKGROUND_GENES = 128


def role_of_seed(seed: int) -> str:
    for role, seeds in SEED_ROLES.items():
        if seed in seeds:
            return role
    raise ValueError(f"seed {seed} does not belong to any role; roles are disjoint by design")


# --------------------------------------------------------------------------- #
# G1: analytic fields
# --------------------------------------------------------------------------- #
def grid_xy(n: int = GRID) -> np.ndarray:
    yy, xx = np.mgrid[0:n, 0:n]
    return np.column_stack([xx.ravel(), yy.ravel()]).astype(float)


def gaussian(xy: np.ndarray, center: Sequence[float], sigma: float, amplitude: float = 1.0) -> np.ndarray:
    d2 = np.sum((xy - np.asarray(center, dtype=float)) ** 2, axis=1)
    return amplitude * np.exp(-d2 / (2.0 * sigma ** 2))


def ring(xy: np.ndarray, center: Sequence[float], radius: float, thickness: float, amplitude: float = 1.0) -> np.ndarray:
    d = np.linalg.norm(xy - np.asarray(center, dtype=float), axis=1)
    return amplitude * np.exp(-0.5 * ((d - radius) / thickness) ** 2)


def sided_ring(
    xy: np.ndarray,
    center: Sequence[float],
    radius: float,
    thickness: float,
    *,
    angular_order: int = 2,
    modulation: float = 0.5,
    phase: float = 0.0,
) -> np.ndarray:
    """Ring whose thickness/amplitude vary with angle, radial mean preserved.

    The angular term has zero mean over each ring, so the radial profile of the
    pair member matches the circular ring while the shape differs.
    """
    delta = xy - np.asarray(center, dtype=float)
    d = np.linalg.norm(delta, axis=1)
    theta = np.arctan2(delta[:, 1], delta[:, 0])
    mod = 1.0 + modulation * np.cos(angular_order * theta + phase)
    mod = np.clip(mod, 0.05, None)
    return np.exp(-0.5 * ((d - radius) / (thickness * mod)) ** 2)


def match_radial_profile(
    target: np.ndarray, reference: np.ndarray, xy: np.ndarray, *, n_bins: int = 24
) -> np.ndarray:
    """Rescale ``target`` so its radial profile equals ``reference``'s, bin by bin.

    Used to build the S1 control: same radial mean, different angular shape.
    The residual (after binning) is reported in the pair diagnostics, never assumed
    to be exactly zero.
    """
    center = np.sum(xy * reference[:, None], axis=0) / max(np.sum(reference), 1e-12)
    d = np.linalg.norm(xy - center, axis=1)
    edges = np.linspace(0, d.max() + 1e-9, n_bins + 1)
    idx = np.clip(np.digitize(d, edges) - 1, 0, n_bins - 1)
    prof_target = np.bincount(idx, weights=target, minlength=n_bins) / np.maximum(
        np.bincount(idx, minlength=n_bins), 1
    )
    prof_ref = np.bincount(idx, weights=reference, minlength=n_bins) / np.maximum(
        np.bincount(idx, minlength=n_bins), 1
    )
    ratio = np.where(prof_target > 1e-12, prof_ref / np.maximum(prof_target, 1e-12), 1.0)
    return target * ratio[idx]


def tube(xy: np.ndarray, path: np.ndarray, width: float) -> np.ndarray:
    """Distance to a polyline, mapped through a Gaussian cross-section."""
    out = np.full(xy.shape[0], np.inf)
    for a, b in zip(path[:-1], path[1:]):
        ab = b - a
        denom = float(ab @ ab)
        if denom == 0:
            d = np.linalg.norm(xy - a, axis=1)
        else:
            t = np.clip(((xy - a) @ ab) / denom, 0.0, 1.0)
            proj = a + t[:, None] * ab
            d = np.linalg.norm(xy - proj, axis=1)
        out = np.minimum(out, d)
    return np.exp(-0.5 * (out / width) ** 2)


def bezier(points: np.ndarray, n: int = 60) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    t = np.linspace(0.0, 1.0, n)[:, None]
    degree = len(points) - 1
    out = np.zeros((n, 2))
    from math import comb

    for i, p in enumerate(points):
        basis = comb(degree, i) * (t ** i) * ((1 - t) ** (degree - i))
        out += basis * p
    return out


def bimodal(xy: np.ndarray, center: Sequence[float], sigma: float, separation: float, angle: float = 0.0) -> np.ndarray:
    c = np.asarray(center, dtype=float)
    direction = np.array([np.cos(angle), np.sin(angle)])
    a = c + direction * separation / 2.0
    b = c - direction * separation / 2.0
    return 0.5 * (gaussian(xy, a, sigma) + gaussian(xy, b, sigma))


# --------------------------------------------------------------------------- #
# G2: a different mechanism (anisotropic diffusion on the grid graph)
# --------------------------------------------------------------------------- #
def diffuse_field(
    n: int = GRID,
    *,
    sources: Sequence[Sequence[float]] = ((20, 20),),
    steps: int = 120,
    conductance_x: float = 1.0,
    conductance_y: float = 0.35,
    decay: float = 0.02,
) -> np.ndarray:
    """Steady/transient diffusion with anisotropic conductance.

    This is deliberately *not* the geometric model used by route A: mass spreads
    along the grid with direction-dependent conductance, so a generator-identical
    shortcut cannot produce it.
    """
    f = np.zeros((n, n), dtype=np.float64)
    for (x, y) in sources:
        i, j = int(round(y)), int(round(x))
        f[max(i - 1, 0):i + 2, max(j - 1, 0):j + 2] += 1.0
    cx = conductance_x / 4.0
    cy = conductance_y / 4.0
    for _ in range(steps):
        lap = np.zeros_like(f)
        lap[:, 1:] += cx * (f[:, :-1] - f[:, 1:])
        lap[:, :-1] += cx * (f[:, 1:] - f[:, :-1])
        lap[1:, :] += cy * (f[:-1, :] - f[1:, :])
        lap[:-1, :] += cy * (f[1:, :] - f[:-1, :])
        f += lap - decay * f
    f -= f.min()
    if f.max() > 0:
        f /= f.max()
    return f.ravel()


# --------------------------------------------------------------------------- #
# observation model
# --------------------------------------------------------------------------- #
@dataclass
class SyntheticSample:
    base_id: str
    scenario: str
    seed: int
    role: str
    pair_id: str
    fields: dict[str, np.ndarray]
    counts: np.ndarray  # (n_points, n_genes)
    visible_mask: np.ndarray
    structures: dict[str, Any] = field(default_factory=dict)
    meta: dict[str, Any] = field(default_factory=dict)

    @property
    def n_points(self) -> int:
        return int(self.counts.shape[0])

    def gene_names(self) -> list[str]:
        names = [f"{fid}_m{k}" for fid in sorted(self.fields) for k in range(MARKERS_PER_FIELD)]
        names += [f"bg_{k}" for k in range(BACKGROUND_GENES)]
        return names

    def field_gene_sets(self) -> dict[str, list[str]]:
        return {fid: [f"{fid}_m{k}" for k in range(MARKERS_PER_FIELD)] for fid in sorted(self.fields)}


def counts_from_fields(
    sample: SyntheticSample,
    *,
    depth_mu: float = 1.0,
    depth_sd: float = 0.25,
    dispersion: float = 8.0,
    dropout: float = 0.0,
    seed: int = 0,
    technical_gradient: np.ndarray | None = None,
) -> np.ndarray:
    """Negative-binomial counts whose means are built from the latent fields."""
    rng = np.random.default_rng(seed + 991)
    n = sample.n_points
    n_fields = len(sample.fields)
    n_markers = n_fields * MARKERS_PER_FIELD
    mean = np.zeros((n, n_markers + BACKGROUND_GENES), dtype=np.float64)

    for k, fid in enumerate(sorted(sample.fields)):
        f = np.clip(np.asarray(sample.fields[fid], dtype=float), 0.0, None)
        loadings = 1.0 + 0.6 * np.linspace(-1, 1, MARKERS_PER_FIELD)
        mean[:, k * MARKERS_PER_FIELD:(k + 1) * MARKERS_PER_FIELD] = (
            np.outer(f, loadings) * 6.0 + 0.5
        )
    bg = rng.gamma(2.0, 0.25, size=BACKGROUND_GENES)
    mean[:, n_markers:] = bg[None, :]

    depth = rng.lognormal(mean=0.0, sigma=depth_sd, size=n) * depth_mu
    if technical_gradient is not None:
        depth = depth * np.asarray(technical_gradient, dtype=float)
    mean = mean * depth[:, None]

    theta = dispersion
    lam = rng.gamma(shape=theta, scale=mean / theta)
    counts = rng.poisson(lam).astype(np.float64)
    if dropout > 0:
        keep = rng.random(counts.shape) > dropout
        counts = np.where(keep, counts, 0.0)
    return counts


# --------------------------------------------------------------------------- #
# scenario builders
# --------------------------------------------------------------------------- #
def _empty_sample(scenario: str, seed: int, pair_id: str, suffix: str = "") -> SyntheticSample:
    base_id = hashlib.sha1(f"{scenario}|{seed}|{pair_id}|{suffix}".encode()).hexdigest()[:12]
    return SyntheticSample(
        base_id=base_id,
        scenario=scenario,
        seed=seed,
        role=role_of_seed(seed),
        pair_id=pair_id,
        fields={},
        counts=np.zeros((GRID * GRID, 1)),
        visible_mask=np.ones(GRID * GRID, dtype=bool),
    )


def build_s0(seed: int, pair_id: str = "s0_null") -> SyntheticSample:
    """Independent smooth fields sharing a technical depth gradient; no relation."""
    rng = np.random.default_rng(seed)
    xy = grid_xy()
    sample = _empty_sample("S0", seed, pair_id)
    for k in range(4):
        centers = rng.uniform(10, GRID - 10, size=(3, 2))
        f = sum(gaussian(xy, c, rng.uniform(6, 14), 1.0 / 3) for c in centers)
        sample.fields[f"f{k}"] = f / max(f.max(), 1e-9)
    yy, xx = np.mgrid[0:GRID, 0:GRID]
    sample.meta["technical_gradient"] = (1.0 + 0.8 * (xx.ravel() / GRID)).tolist()
    sample.structures = {"relation": "none"}
    return sample


def build_s1(seed: int, pair_id: str = "s1_shape") -> SyntheticSample:
    """Matched radial profile: circular ring vs angularly thickened ring."""
    rng = np.random.default_rng(seed + 7)
    xy = grid_xy()
    center = np.array([GRID / 2 + rng.uniform(-4, 4), GRID / 2 + rng.uniform(-4, 4)])
    radius = rng.uniform(12, 18)
    thickness = rng.uniform(2.5, 4.0)
    sample = _empty_sample("S1", seed, pair_id)
    base_ring = ring(xy, center, radius, thickness)
    shaped = sided_ring(
        xy, center, radius, thickness, angular_order=2, modulation=0.6, phase=rng.uniform(0, np.pi)
    )
    sample.fields["circular"] = base_ring
    sample.fields["sided"] = match_radial_profile(shaped, base_ring, xy)
    # second pair: unimodal vs bimodal with the same centroid and mass
    sigma = rng.uniform(4, 6)
    sep = 2.2 * sigma
    uni = gaussian(xy, center, sigma * 1.25)
    bi = bimodal(xy, center, sigma, sep, angle=rng.uniform(0, np.pi))
    bi = bi * (uni.sum() / max(bi.sum(), 1e-9))
    sample.fields["unimodal"] = uni
    sample.fields["bimodal"] = bi
    sample.structures = {"center": center.tolist(), "radius": radius, "thickness": thickness}
    return sample


def build_s2(seed: int, pair_id: str = "s2_coupling") -> SyntheticSample:
    """Same-point correlation matched by *noise level*, with a pure spatial offset.

    Matching the correlation by mixing in the anchor would create a spurious lag-0
    component and make the offset pair indistinguishable from a same-point pair by
    construction.  Here the offset partner is a pure shifted copy plus independent
    noise, and the noise level is solved so that its same-point correlation equals
    the same-point partner's.
    """
    rng = np.random.default_rng(seed + 13)
    xy = grid_xy()
    c = np.array([GRID / 2, GRID / 2]) + rng.uniform(-5, 5, 2)
    shift = np.array([7.0, 0.0])
    base = gaussian(xy, c, 9.0)
    shifted = gaussian(xy, c + shift, 9.0)
    noise = rng.normal(size=base.shape)
    noise /= np.std(noise)

    base_noise = rng.normal(size=base.shape)
    base_noise /= np.std(base_noise)
    target_noise = 0.3
    partner_same = base + target_noise * base_noise
    target = float(np.corrcoef(base, partner_same)[0, 1])

    lo, hi = 0.0, 20.0
    best = (lo, shifted, float(np.corrcoef(base, shifted)[0, 1]))
    for _ in range(80):
        s = 0.5 * (lo + hi)
        cand = shifted + s * noise
        corr = float(np.corrcoef(base, cand)[0, 1])
        best = (s, cand, corr)
        if corr > target:
            lo = s
        else:
            hi = s
    s, candidate, achieved = best

    sample = _empty_sample("S2", seed, pair_id)
    sample.fields["anchor"] = base
    sample.fields["partner_same_point"] = partner_same
    sample.fields["partner_offset"] = candidate
    sample.structures = {
        "requested_offset": shift.tolist(),
        "same_point_corr_target": target,
        "same_point_corr_matched": achieved,
        "noise_sd_on_offset_partner": float(s),
        "construction": "pure shifted copy + independent noise; correlation matched by noise level only",
        "achieved_peak_shift_pitch": float(shift[0]),
        "separable": True,
    }
    return sample


def _bridge_paths(rng: np.random.Generator) -> dict[str, np.ndarray]:
    left = np.array([12.0, GRID / 2])
    right = np.array([GRID - 12.0, GRID / 2])
    wrong = np.array([GRID - 12.0, GRID / 2 - 16.0])
    mid = (left + right) / 2
    return {
        "continuous": bezier([left, mid + np.array([0.0, 4.0]), right]),
        "broken": None,  # handled by masking the middle
        "wrong_end": bezier([left, mid + np.array([0.0, -6.0]), wrong]),
        "y_branch": None,
        "ring": None,
        "left": left,
        "right": right,
        "wrong": wrong,
    }


def build_s3(seed: int, pair_id: str = "s3_connectivity") -> SyntheticSample:
    """True bridge vs broken vs wrong-end vs Y-branch vs redundant ring."""
    rng = np.random.default_rng(seed + 29)
    xy = grid_xy()
    paths = _bridge_paths(rng)
    width = 1.6
    continuous = tube(xy, paths["continuous"], width)
    broken = continuous.copy()
    # remove the middle third along the bridge's axis
    mid_x = 0.5 * (paths["continuous"][:, 0].min() + paths["continuous"][:, 0].max())
    broken[np.abs(xy[:, 0] - mid_x) < 4.0] = 0.0
    wrong = tube(xy, paths["wrong_end"], width)
    ring_path = bezier(
        [paths["left"], [GRID / 2, GRID / 2 - 20], paths["right"], [GRID / 2, GRID / 2 + 18], paths["left"]]
    )
    redundant = np.maximum(continuous, 0.7 * tube(xy, ring_path, width))
    y_branch_path = bezier([paths["left"], [GRID / 2, GRID / 2], [GRID / 2 + 6, GRID / 2 + 16]])
    y_branch = np.maximum(continuous, 0.9 * tube(xy, y_branch_path, width))

    sample = _empty_sample("S3", seed, pair_id)
    sample.fields["endpoint_left"] = gaussian(xy, paths["left"], 3.5)
    sample.fields["endpoint_right"] = gaussian(xy, paths["right"], 3.5)
    sample.fields["endpoint_wrong"] = gaussian(xy, paths["wrong"], 3.5)
    sample.fields["bridge_continuous"] = continuous
    sample.fields["bridge_broken"] = broken
    sample.fields["bridge_wrong_end"] = wrong
    sample.fields["bridge_redundant"] = redundant
    sample.fields["bridge_y_branch"] = y_branch
    sample.structures = {
        "left": paths["left"].tolist(),
        "right": paths["right"].tolist(),
        "wrong": paths["wrong"].tolist(),
        "width_pitch": width,
        "gap_half_width": 4.0,
    }
    return sample


def build_s4(seed: int, pair_id: str = "s4_distractors", n_distractors: int = 0) -> SyntheticSample:
    """A true bridge plus N unrelated smooth fields (field-count inflation)."""
    sample = build_s3(seed, pair_id=f"{pair_id}_n{n_distractors}")
    sample.scenario = "S4"
    rng = np.random.default_rng(seed + 101 + n_distractors)
    xy = grid_xy()
    for k in range(n_distractors):
        c = rng.uniform(8, GRID - 8, size=2)
        f = gaussian(xy, c, rng.uniform(5, 12))
        f = f / max(f.max(), 1e-9)
        sample.fields[f"distractor_{k:02d}"] = f
    sample.structures["n_distractors"] = n_distractors
    return sample


def build_s5(seed: int, pair_id: str = "s5_nonidentifiable") -> SyntheticSample:
    """Identical visible data, different hidden structure (bridge vs no bridge)."""
    sample = build_s3(seed, pair_id=pair_id)
    xy = grid_xy()
    mask = np.ones(GRID * GRID, dtype=bool)
    hole = np.linalg.norm(xy - np.array([GRID / 2, GRID / 2]), axis=1) <= 8.0
    mask[hole] = False
    sample.fields = {"bridge_continuous": sample.fields["bridge_continuous"]}
    sample.fields["bridge_absent"] = np.zeros(GRID * GRID)
    sample.visible_mask = mask
    sample.structures["mask"] = "disc centre radius 8"
    sample.structures["same_visible_values"] = True
    return sample


def build_s3_parallel(seed: int, pair_id: str = "s3_parallel") -> SyntheticSample:
    """Bridge supported by a *parallel* companion field on one half only.

    Neither field alone lies on the corridor for the whole route: field A covers the
    left half on the corridor, field B covers the right half but sits 3 pitch to the
    side.  A single-field or max-aggregate router follows B off the corridor; a
    relation-aware common support can use the known parallel offset to stay on it.

    This is our own extension (not one of the package's S0-S5).
    """
    rng = np.random.default_rng(seed + 71)
    xy = grid_xy()
    left = np.array([12.0, GRID / 2])
    right = np.array([GRID - 12.0, GRID / 2])
    path = bezier([left, [GRID / 2, GRID / 2 + 4.0], right])
    corridor = tube(xy, path, 1.6)

    xs = xy[:, 0]
    on_left = (xs <= GRID / 2 + 2).astype(float)
    on_right = (xs >= GRID / 2 - 2).astype(float)
    field_a = corridor * on_left
    companion = np.zeros_like(corridor)
    shifted = tube(xy - np.array([0.0, 6.0]), path, 1.6)
    companion = shifted * on_right
    noise = 0.02 * rng.normal(size=corridor.shape)
    sample = _empty_sample("S3p", seed, pair_id)
    sample.fields["corridor"] = corridor
    sample.fields["support_on_corridor"] = np.clip(field_a + noise, 0, None)
    sample.fields["support_parallel_right"] = np.clip(companion + noise, 0, None)
    sample.structures = {
        "left": left.tolist(),
        "right": right.tolist(),
        "parallel_offset_pitch": [0.0, 6.0],
        "coverage": "field A left half, field B right half offset +6 pitch in y",
        "extension": "our own S3 variant, not part of the package's S0-S5",
    }
    return sample


def make_s5_pair(seed: int) -> tuple[SyntheticSample, SyntheticSample]:
    """Two samples with bit-identical visible counts and different hidden structure."""
    a = build_s3(seed, pair_id="s5_nonidentifiable")
    a.scenario = "S5"
    xy = grid_xy()
    mask = np.ones(GRID * GRID, dtype=bool)
    mask[np.linalg.norm(xy - np.array([GRID / 2, GRID / 2]), axis=1) <= 8.0] = False
    a.fields = {"bridge": a.fields["bridge_continuous"]}
    a.visible_mask = mask
    counts = counts_from_fields(a, seed=seed + 5)
    a.counts = counts
    b = SyntheticSample(
        base_id=hashlib.sha1(f"S5|{seed}|absent".encode()).hexdigest()[:12],
        scenario="S5",
        seed=seed,
        role=role_of_seed(seed),
        pair_id="s5_nonidentifiable",
        fields={"bridge": np.zeros(GRID * GRID)},
        counts=counts.copy(),
        visible_mask=mask.copy(),
        structures={"variant": "bridge absent in hidden disc", "same_visible_values": True},
    )
    b.meta["pair_with"] = a.base_id
    a.meta["pair_with"] = b.base_id
    a.structures["variant"] = "bridge present in hidden disc"
    a.structures["same_visible_values"] = True
    return a, b


SCENARIOS = {
    "S0": build_s0,
    "S1": build_s1,
    "S2": build_s2,
    "S3": build_s3,
    "S4": build_s4,
    "S5": build_s5,
}


# --------------------------------------------------------------------------- #
# materialisation
# --------------------------------------------------------------------------- #
def materialise(sample: SyntheticSample, root: Path, *, counts_seed: int | None = None) -> Path:
    """Write observed_inputs / latent_truth / evaluation_only for one sample."""
    base = root / "samples" / f"{sample.scenario}_{sample.base_id}"
    (base / "observed_inputs").mkdir(parents=True, exist_ok=True)
    (base / "latent_truth").mkdir(parents=True, exist_ok=True)
    (base / "evaluation_only").mkdir(parents=True, exist_ok=True)

    genes = np.array(sample.gene_names())
    counts = sample.counts
    if counts.shape[1] != len(genes):
        gradient = sample.meta.get("technical_gradient")
        counts = counts_from_fields(
            sample,
            seed=counts_seed if counts_seed is not None else sample.seed,
            technical_gradient=np.asarray(gradient) if gradient is not None else None,
        )
        sample.counts = counts
    if counts.shape[1] != len(genes):
        raise RuntimeError("count matrix width does not match the gene list")

    visible = np.asarray(sample.visible_mask, dtype=bool)
    observed = counts.copy()
    observed[~visible] = np.nan
    np.savez_compressed(
        base / "observed_inputs" / "counts_visible.npz",
        counts=observed,
        genes=genes,
        xy=grid_xy(GRID),
        visible_mask=visible,
        coordinate_unit="synthetic_pitch",
    )
    np.savez_compressed(
        base / "latent_truth" / "fields.npz",
        **{fid: np.asarray(v, dtype=np.float64) for fid, v in sample.fields.items()},
    )
    np.savez_compressed(
        base / "evaluation_only" / "answers.npz",
        visible_mask=visible,
        structures=json.dumps(sample.structures),
        meta=json.dumps({k: v for k, v in sample.meta.items() if k != "technical_gradient"}),
    )
    (base / "meta.json").write_text(
        json.dumps(
            {
                "base_id": sample.base_id,
                "scenario": sample.scenario,
                "pair_id": sample.pair_id,
                "seed": sample.seed,
                "role": sample.role,
                "n_fields": len(sample.fields),
                "field_ids": sorted(sample.fields),
                "n_points": sample.n_points,
                "n_visible": int(visible.sum()),
                "n_genes": len(genes),
                "markers_per_field": MARKERS_PER_FIELD,
                "background_genes": BACKGROUND_GENES,
            },
            indent=2,
        )
    )
    return base


def pair_diagnostics(samples: Sequence[SyntheticSample]) -> list[dict[str, Any]]:
    """Are the matched pairs actually matched? Report the residuals, do not assume."""
    xy = grid_xy()
    rows: list[dict[str, Any]] = []
    for s in samples:
        for name, f in sorted(s.fields.items()):
            rows.append(
                {
                    "base_id": s.base_id,
                    "scenario": s.scenario,
                    "seed": s.seed,
                    "field_id": name,
                    "mass": float(np.sum(f)),
                    "max": float(np.max(f)),
                    "centroid_x": float(np.sum(xy[:, 0] * f) / max(np.sum(f), 1e-12)),
                    "centroid_y": float(np.sum(xy[:, 1] * f) / max(np.sum(f), 1e-12)),
                    "n_modes_above_halfmax": int(_count_modes(f)),
                }
            )
    for s in samples:
        if s.scenario == "S2":
            rows.append(
                {
                    "base_id": s.base_id,
                    "scenario": "S2_pair",
                    "seed": s.seed,
                    "field_id": "anchor_vs_partners",
                    "mass": float(np.corrcoef(s.fields["anchor"], s.fields["partner_same_point"])[0, 1]),
                    "max": float(np.corrcoef(s.fields["anchor"], s.fields["partner_offset"])[0, 1]),
                    "centroid_x": float(s.structures.get("achieved_peak_shift_pitch", float("nan"))),
                    "centroid_y": float(s.structures.get("requested_offset", [0, 0])[1]),
                    "n_modes_above_halfmax": 0,
                    "note": "'mass'=same-point corr of same-point partner, 'max'=same-point corr of offset partner, "
                            "'centroid_x'=achieved peak shift in pitch, 'centroid_y'=requested offset",
                }
            )
            continue
        if s.scenario != "S1":
            continue
        a, b = s.fields["circular"], s.fields["sided"]
        radial_a = _radial_profile(a, xy)
        radial_b = _radial_profile(b, xy)
        rows.append(
            {
                "base_id": s.base_id,
                "scenario": "S1_pair",
                "seed": s.seed,
                "field_id": "circular_vs_sided",
                "mass": float(np.sum(a) - np.sum(b)),
                "max": float(np.max(np.abs(radial_a - radial_b))),
                "centroid_x": float(np.sum(a) - np.sum(b)),
                "centroid_y": 0.0,
                "n_modes_above_halfmax": 0,
                "note": "max|radial profile difference| in column 'max'; mass difference in 'mass'",
            }
        )
        u, v = s.fields["unimodal"], s.fields["bimodal"]
        rows.append(
            {
                "base_id": s.base_id,
                "scenario": "S1_pair",
                "seed": s.seed,
                "field_id": "unimodal_vs_bimodal",
                "mass": float(np.sum(u) - np.sum(v)),
                "max": float(np.max(v)),
                "centroid_x": float(
                    np.sum(xy[:, 0] * u) / max(np.sum(u), 1e-12) - np.sum(xy[:, 0] * v) / max(np.sum(v), 1e-12)
                ),
                "centroid_y": 0.0,
                "n_modes_above_halfmax": int(_count_modes(v)),
                "note": "centroid shift in column centroid_x",
            }
        )
    return rows


def _radial_profile(field: np.ndarray, xy: np.ndarray, n_bins: int = 24) -> np.ndarray:
    center = np.sum(xy * field[:, None], axis=0) / max(np.sum(field), 1e-12)
    d = np.linalg.norm(xy - center, axis=1)
    bins = np.linspace(0, d.max() + 1e-9, n_bins + 1)
    idx = np.clip(np.digitize(d, bins) - 1, 0, n_bins - 1)
    sums = np.bincount(idx, weights=field, minlength=n_bins)
    counts = np.bincount(idx, minlength=n_bins).astype(float)
    return np.where(counts > 0, sums / np.maximum(counts, 1), 0.0)


def _count_modes(field: np.ndarray) -> int:
    """Count local maxima above half maximum on the 64x64 grid."""
    img = field.reshape(GRID, GRID)
    peak = img.max()
    if peak <= 0:
        return 0
    above = img > 0.5 * peak
    count = 0
    for i in range(1, GRID - 1):
        for j in range(1, GRID - 1):
            if not above[i, j]:
                continue
            window = img[i - 1:i + 2, j - 1:j + 2]
            if img[i, j] >= window.max() - 1e-12:
                count += 1
    return count
