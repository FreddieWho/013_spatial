"""Field analyzer package (FA round, 2026-10-20).

Scope of this package: the *single information entry* used by all routes
(A: deformable geometry, B: cross-field relation operator, C: multi-field
connectivity). Nothing in this package may read hidden-region expression.

Contract authority:
  docs/FA_CONTRACT_ADDENDUM_20261020.md   (project-side, wins on conflict)
  docs/field_analyzer_killtest_20261020/05_COMMON_PROTOCOL.md
"""

from .observation import (
    ROLE_A_FEATURE,
    ROLE_B_WITNESS,
    ROLE_BACKGROUND,
    ROLE_EXCLUDED,
    ROLES,
    GeneRoles,
    HiddenTruth,
    MaskViolation,
    ObservationBundle,
    fingerprint_arrays,
    fingerprint_np,
    load_observation,
    save_observation,
    split_for_evaluation,
    visible_library_size,
    visible_mean_std,
    visible_detection_rate,
)
from .adapters import AxisSpec, ReadReport, read_h5ad_counts
from . import geometry

__all__ = [
    "ROLE_A_FEATURE",
    "ROLE_B_WITNESS",
    "ROLE_BACKGROUND",
    "ROLE_EXCLUDED",
    "ROLES",
    "GeneRoles",
    "HiddenTruth",
    "MaskViolation",
    "ObservationBundle",
    "fingerprint_arrays",
    "fingerprint_np",
    "load_observation",
    "save_observation",
    "split_for_evaluation",
    "visible_library_size",
    "visible_mean_std",
    "visible_detection_rate",
    "AxisSpec",
    "ReadReport",
    "read_h5ad_counts",
    "geometry",
]
