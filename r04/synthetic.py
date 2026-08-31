"""Deterministic fixtures for R-04 implementation tests."""

from __future__ import annotations

import numpy as np
from scipy import sparse

from .types import SectionData


def make_overlapping_sections(
    *,
    sections: int = 4,
    spots_per_section: int = 64,
    genes: int = 40,
    seed: int = 20260807,
    composition_only: bool = False,
    antagonistic: bool = False,
) -> tuple[list[SectionData], np.ndarray]:
    rng = np.random.default_rng(seed)
    result: list[SectionData] = []
    truth: list[np.ndarray] = []
    grid_size = int(np.ceil(np.sqrt(spots_per_section)))
    base = np.array([(i, j) for i in range(grid_size) for j in range(grid_size)], dtype=float)[:spots_per_section]
    gene_program_a = np.linspace(-1.0, 1.0, genes)
    gene_program_b = np.cos(np.linspace(0, np.pi * 2, genes))
    for index in range(sections):
        coords = base + rng.normal(0, 0.03, base.shape)
        x = coords[:, 0] / max(1, grid_size - 1)
        y = coords[:, 1] / max(1, grid_size - 1)
        field_a = np.sin(2 * np.pi * x) + 0.5 * y
        field_b = np.cos(2 * np.pi * y) - 0.25 * x
        if composition_only:
            field_a = np.zeros_like(field_a)
            field_b = np.zeros_like(field_b)
        truth.append(np.column_stack([field_a, field_b]))
        log_rate = 1.2 + 0.2 * rng.normal(size=(spots_per_section, genes))
        log_rate += np.outer(field_a, gene_program_a) * 0.12
        if antagonistic:
            log_rate += np.outer(field_b, gene_program_b) * 0.12
        else:
            log_rate += np.outer(field_b, np.maximum(gene_program_b, 0)) * 0.08
        library = rng.lognormal(mean=0.0, sigma=0.15, size=spots_per_section)
        rate = np.exp(log_rate) * library[:, None]
        counts = rng.poisson(rate).astype(np.int32)
        result.append(SectionData(
            section_id=f"S{index+1:02d}", patient_id=f"P{index+1:02d}", block_id=f"B{index+1:02d}",
            lineage="SYNTHETIC", barcode=tuple(f"BC{index:02d}_{i:04d}" for i in range(spots_per_section)),
            coords=coords, counts=sparse.csr_matrix(counts), gene_id=tuple(f"G{i:04d}" for i in range(genes)),
        ))
    return result, np.concatenate(truth, axis=0)


def _irregular_coordinates(
    rng: np.random.Generator,
    domain: str,
    spots: int,
) -> np.ndarray:
    accepted: list[np.ndarray] = []
    while sum(len(chunk) for chunk in accepted) < spots:
        proposal = rng.uniform(-1.0, 1.0, size=(spots * 4, 2))
        accepted.append(proposal[irregular_domain_mask(proposal, domain)])
    return np.concatenate(accepted, axis=0)[:spots]


def _irregular_field_basis(coords: np.ndarray, domain: str) -> tuple[np.ndarray, np.ndarray]:
    """Create smooth effects whose support is not a radial or grid label."""
    x, y = np.asarray(coords, dtype=float).T
    if domain == "crescent":
        first = np.exp(-((x + 0.48) ** 2 + (y - 0.25) ** 2) / 0.12)
        second = np.exp(-((x + 0.38) ** 2 + (y + 0.35) ** 2) / 0.12)
    elif domain == "branch":
        first = np.exp(-((y - 0.85 * x - 0.02) ** 2) / 0.025)
        second = np.exp(-((y + 0.78 * x + 0.08) ** 2) / 0.025)
    elif domain == "disconnected":
        first = np.exp(-(((x + 0.48) / 0.25) ** 2 + ((y + 0.05) / 0.30) ** 2))
        second = np.exp(-(((x - 0.46) / 0.22) ** 2 + ((y - 0.12) / 0.25) ** 2))
    else:
        raise ValueError(f"unsupported irregular calibration domain: {domain}")
    return first - first.mean(), second - second.mean()


def make_identifiability_calibration_sections(
    *,
    mode: str,
    domain: str = "disconnected",
    generator: str = "mnsf_correct",
    sections: int = 3,
    spots_per_section: int = 80,
    genes: int = 32,
    seed: int = 20260807,
) -> tuple[list[SectionData], np.ndarray, dict[str, object]]:
    """Generate explicit continuous-field controls for K and recovery gates.

    ``mode`` is one of ``no_field``, ``one_disconnected``,
    ``two_independent`` or ``two_collinear``.  The returned truth contains
    continuous function values, never spot labels.  ``generator`` records
    whether gene effects are non-negative, signed, or deliberately
    misspecified for a robustness control.
    """
    modes = {"no_field", "one_disconnected", "two_independent", "two_collinear"}
    generators = {"mnsf_correct", "signed_correct", "misspecified"}
    if mode not in modes or generator not in generators:
        raise ValueError("unsupported calibration mode or generator")
    if sections < 1 or spots_per_section < 4 or genes < 4:
        raise ValueError("identifiability calibration dimensions are too small")
    rng = np.random.default_rng(seed)
    gene_id = tuple(f"G{i:04d}" for i in range(genes))
    result: list[SectionData] = []
    truth_parts: list[np.ndarray] = []
    gene_position = np.linspace(0.0, 1.0, genes)
    program_a = np.linspace(0.15, 0.85, genes)
    program_b = np.sin(np.linspace(0.0, 2.0 * np.pi, genes))
    positive_loading = np.column_stack([
        0.03 + np.exp(-0.5 * ((gene_position - 0.25) / 0.12) ** 2),
        0.03 + np.exp(-0.5 * ((gene_position - 0.72) / 0.14) ** 2),
    ])
    positive_loading /= positive_loading.sum(axis=0, keepdims=True)
    positive_amplitude = np.asarray([60.0, 45.0])
    if generator == "misspecified":
        program_a = np.linspace(-1.0, 1.0, genes)
    for index in range(sections):
        coords = _irregular_coordinates(rng, domain, spots_per_section)
        coords += rng.normal(0.0, 0.008, size=coords.shape)
        first, second = _irregular_field_basis(coords, domain)
        if mode == "no_field":
            fields = np.empty((spots_per_section, 0), dtype=float)
        elif mode == "one_disconnected":
            fields = first[:, None]
        elif mode == "two_collinear":
            fields = np.column_stack([first, first])
        else:
            fields = np.column_stack([first, second])
        truth_parts.append(fields)
        if generator == "mnsf_correct":
            active = fields.shape[1]
            spatial_rate = np.zeros((spots_per_section, genes), dtype=float)
            for factor in range(active):
                spatial_rate += np.exp(fields[:, factor, None]) * (
                    positive_loading[:, factor] * positive_amplitude[factor]
                )[None, :]
            rate = 2.5 + spatial_rate
        elif fields.shape[1] == 0:
            log_rate = np.full((spots_per_section, genes), 1.2)
            rate = np.exp(log_rate)
        elif fields.shape[1] == 1:
            log_rate = 1.2 + np.outer(fields[:, 0], program_a) * 0.55
            rate = np.exp(log_rate)
        else:
            log_rate = 1.2 + np.outer(fields[:, 0], program_a) * 0.55
            log_rate += np.outer(fields[:, 1], program_b) * 0.45
            rate = np.exp(log_rate)
        if generator == "misspecified":
            rate = np.exp(log_rate + 0.35 * np.sin(3.0 * coords[:, 0:1] + coords[:, 1:2]))
        library = rng.lognormal(mean=0.0, sigma=0.12, size=spots_per_section)
        counts = rng.poisson(rate * library[:, None]).astype(np.int32)
        result.append(SectionData(
            section_id=f"CAL_{index + 1:02d}",
            patient_id=f"CAL_P{index + 1:02d}",
            block_id=f"CAL_B{index + 1:02d}",
            lineage=f"CALIBRATION_{domain}",
            barcode=tuple(f"CALBC_{index:02d}_{spot:04d}" for spot in range(spots_per_section)),
            coords=coords,
            counts=sparse.csr_matrix(counts),
            gene_id=gene_id,
            library_size=library,
            metadata={"calibration_mode": mode, "generator": generator, "domain": domain},
        ))
    truth = np.concatenate(truth_parts, axis=0)
    metadata: dict[str, object] = {
        "mode": mode,
        "domain": domain,
        "generator": generator,
        "generator_version": (
            "crescent_supported_fields_v3"
            if domain == "crescent" else "independent_gene_effects_v2"
        ),
        "k_true": int(truth.shape[1]),
        "field_standard_deviation": truth.std(axis=0).tolist() if truth.size else [],
        "field_definition": "smooth continuous coordinate functions on an irregular support",
    }
    if generator == "mnsf_correct":
        active = truth.shape[1]
        true_loading = positive_loading[:, :active]
        true_amplitude = positive_amplitude[:active]
        true_effect = (
            np.exp(truth) @ (true_loading * true_amplitude).T
            if active else np.empty((len(truth), genes), dtype=float)
        )
        metadata.update({
            "true_loading": true_loading.tolist(),
            "true_amplitude": true_amplitude.tolist(),
            "true_baseline": 2.5,
            "k_eff_true": int(np.linalg.matrix_rank(
                true_effect - true_effect.mean(axis=0, keepdims=True)
            )) if active else 0,
        })
    else:
        metadata["k_eff_true"] = int(np.linalg.matrix_rank(truth)) if truth.size else 0
    return result, truth, metadata


def irregular_domain_mask(coords: np.ndarray, domain: str) -> np.ndarray:
    """Return a non-grid domain mask for continuous-field calibration.

    The domains are unions/differences of smooth inequalities, not labels or
    clusters.  They are used only to test whether a field model follows an
    irregular support and remains continuous inside it.
    """
    coords = np.asarray(coords, dtype=float)
    x, y = coords[:, 0], coords[:, 1]
    if domain == "crescent":
        outer = (x / 0.92) ** 2 + (y / 0.72) ** 2 <= 1.0
        inner = ((x - 0.22) / 0.70) ** 2 + ((y + 0.02) / 0.54) ** 2 <= 1.0
        return outer & ~inner
    if domain == "branch":
        trunk = (x < 0.25) & (np.abs(y - 0.10 * x) < 0.12)
        upper = (x >= 0.0) & (np.abs(y - (0.85 * x + 0.02)) < 0.11)
        lower = (x >= 0.0) & (np.abs(y + (0.78 * x) + 0.08) < 0.11)
        return trunk | upper | lower
    if domain == "disconnected":
        left = ((x + 0.48) / 0.32) ** 2 + ((y + 0.05) / 0.45) ** 2 <= 1.0
        right = ((x - 0.46) / 0.28) ** 2 + ((y - 0.12) / 0.36) ** 2 <= 1.0
        return left | right
    raise ValueError(f"unsupported irregular domain: {domain}")


def make_irregular_continuous_field_sections(
    *,
    domain: str = "crescent",
    sections: int = 3,
    spots_per_section: int = 80,
    genes: int = 24,
    seed: int = 20260807,
) -> tuple[list[SectionData], np.ndarray]:
    """Create counts driven by smooth fields on an irregular 2-D support.

    The two planted fields are explicit continuous functions of coordinates:
    ``f1(x,y)=sin(2x+y)+0.35 exp(-d1^2/0.18)`` and
    ``f2(x,y)=cos(x-1.4y)-0.25 tanh(3(x+y))``.  The support can be a crescent,
    a branching tube, or disconnected components, so no radial/rectangular
    assumption is needed by the calibration.
    """
    if spots_per_section < 2 or genes < 4:
        raise ValueError("irregular calibration needs at least 2 spots and 4 genes")
    rng = np.random.default_rng(seed)
    gene_id = tuple(f"G{i:04d}" for i in range(genes))
    program_a = np.linspace(-1.0, 1.0, genes)
    program_b = np.sin(np.linspace(0.0, 2.0 * np.pi, genes))
    result: list[SectionData] = []
    truth: list[np.ndarray] = []
    for index in range(sections):
        accepted: list[np.ndarray] = []
        while sum(len(chunk) for chunk in accepted) < spots_per_section:
            proposal = rng.uniform(-1.0, 1.0, size=(spots_per_section * 4, 2))
            accepted.append(proposal[irregular_domain_mask(proposal, domain)])
        coords = np.concatenate(accepted, axis=0)[:spots_per_section]
        coords += rng.normal(0.0, 0.012, size=coords.shape)
        x, y = coords[:, 0], coords[:, 1]
        f1 = np.sin(2.0 * x + y) + 0.35 * np.exp(-((x + 0.35) ** 2 + (y - 0.2) ** 2) / 0.18)
        f2 = np.cos(x - 1.4 * y) - 0.25 * np.tanh(3.0 * (x + y))
        fields = np.column_stack([f1, f2])
        truth.append(fields)
        log_rate = 1.1 + np.outer(f1, program_a) * 0.22 + np.outer(f2, program_b) * 0.18
        library = rng.lognormal(mean=0.0, sigma=0.18, size=spots_per_section)
        counts = rng.poisson(np.exp(log_rate) * library[:, None]).astype(np.int32)
        result.append(SectionData(
            section_id=f"IRREG_{index + 1:02d}",
            patient_id=f"IRREG_P{index + 1:02d}",
            block_id=f"IRREG_B{index + 1:02d}",
            lineage=f"IRREG_{domain}",
            barcode=tuple(f"IRBC_{index:02d}_{spot:04d}" for spot in range(spots_per_section)),
            coords=coords,
            counts=sparse.csr_matrix(counts),
            gene_id=gene_id,
        ))
    return result, np.concatenate(truth, axis=0)
