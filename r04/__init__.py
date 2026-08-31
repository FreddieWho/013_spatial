"""R-04 latent spatial-field discovery control plane.

The package deliberately exposes continuous fields and their uncertainty.  It
does not expose a clustering-based candidate API.
"""

from .gene_selection import GeneSelectionReport, select_training_genes
from .types import CandidateField, FieldFit, SectionData

__all__ = [
    "CandidateField", "FieldFit", "SectionData",
    "GeneSelectionReport", "select_training_genes",
]
