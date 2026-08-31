"""Outcome-blind training-only gene universe selection for R-04."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse

from .types import R04ContractError, SectionData


MITOCHONDRIAL_GENE_IDS = frozenset({
    "ENSG00000198888", "ENSG00000198763", "ENSG00000198804", "ENSG00000198712",
    "ENSG00000228253", "ENSG00000198899", "ENSG00000198938", "ENSG00000198840",
    "ENSG00000212907", "ENSG00000198886", "ENSG00000198786", "ENSG00000198695",
    "ENSG00000198727",
})


@dataclass(frozen=True)
class GeneSelectionReport:
    gene_id: tuple[str, ...]
    common_gene_count: int
    eligible_gene_count: int
    requested_gene_count: int
    patient_count: int
    patient_prevalence_threshold: float
    spot_detection_threshold: float
    excluded_technical_count: int


def _technical_gene(gene: str) -> bool:
    upper = gene.upper().split(".", 1)[0]
    return upper in MITOCHONDRIAL_GENE_IDS or upper.startswith(("MT-", "MT_", "MT")) or upper.startswith(("RNA18S", "RNA28S", "RN45S", "RNA5S"))


def select_training_genes(
    sections: list[SectionData],
    *,
    n_genes: int = 4000,
    min_patient_fraction: float = 0.70,
    min_spot_fraction: float = 0.05,
) -> GeneSelectionReport:
    if not sections or n_genes <= 0:
        raise R04ContractError("training gene selection requires sections and a positive gene count")
    patients = sorted({section.patient_id for section in sections})
    if len(patients) < 2:
        raise R04ContractError("training gene selection requires at least two patients")
    common = set(sections[0].gene_id)
    for section in sections:
        if len(set(section.gene_id)) != len(section.gene_id):
            raise R04ContractError(f"duplicate gene IDs in {section.section_id}")
        common &= set(section.gene_id)
    genes = tuple(sorted(common))
    if not genes:
        raise R04ContractError("training sections have no common genes")
    position = {gene: index for index, gene in enumerate(genes)}
    n = len(genes)
    total_spots = 0
    total_sum = np.zeros(n, dtype=float)
    total_sum_sq = np.zeros(n, dtype=float)
    patient_detected = {patient: np.zeros(n, dtype=float) for patient in patients}
    patient_spots = {patient: 0 for patient in patients}

    for section in sections:
        counts = section.counts.tocsr() if sparse.issparse(section.counts) else sparse.csr_matrix(section.counts)
        libraries = np.asarray(counts.sum(axis=1)).ravel().astype(float)
        if np.any(libraries <= 0):
            raise R04ContractError(f"non-positive library size in {section.section_id}")
        section_columns = np.asarray([index for index, gene in enumerate(section.gene_id) if gene in position], dtype=int)
        columns = np.asarray([position[section.gene_id[index]] for index in section_columns], dtype=int)
        counts = counts[:, section_columns]
        normalized = counts.multiply((10000.0 / libraries)[:, None]).tocsr()
        normalized.data = np.log1p(normalized.data)
        section_sum = np.asarray(normalized.sum(axis=0)).ravel()
        section_sum_sq = np.asarray(normalized.multiply(normalized).sum(axis=0)).ravel()
        detected = np.asarray((counts > 0).sum(axis=0)).ravel().astype(float)
        total_sum[columns] += section_sum
        total_sum_sq[columns] += section_sum_sq
        patient_detected[section.patient_id][columns] += detected
        patient_spots[section.patient_id] += counts.shape[0]
        total_spots += counts.shape[0]

    prevalence = np.zeros(n, dtype=float)
    for patient in patients:
        prevalence += (patient_detected[patient] / patient_spots[patient]) >= min_spot_fraction
    prevalence /= len(patients)
    variance = (total_sum_sq - (total_sum * total_sum) / total_spots) / max(total_spots - 1, 1)
    variance = np.maximum(variance, 0.0)
    technical = np.asarray([_technical_gene(gene) for gene in genes], dtype=bool)
    eligible = (prevalence >= min_patient_fraction) & ~technical
    eligible_indices = np.flatnonzero(eligible)
    if len(eligible_indices) < n_genes:
        raise R04ContractError(
            f"training gene universe has only {len(eligible_indices)} eligible genes; {n_genes} requested"
        )
    ranked = sorted(eligible_indices.tolist(), key=lambda index: (-float(variance[index]), genes[index]))
    selected = tuple(genes[index] for index in ranked[:n_genes])
    return GeneSelectionReport(
        gene_id=selected,
        common_gene_count=len(genes),
        eligible_gene_count=len(eligible_indices),
        requested_gene_count=n_genes,
        patient_count=len(patients),
        patient_prevalence_threshold=min_patient_fraction,
        spot_detection_threshold=min_spot_fraction,
        excluded_technical_count=int(technical.sum()),
    )
