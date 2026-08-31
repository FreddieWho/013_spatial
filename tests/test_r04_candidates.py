from __future__ import annotations

import numpy as np

from r04.candidates import match_model_factors
from r04.types import FieldFit


def _fit(model: str, factor: str, section: str, patient: str, loading: np.ndarray, field: np.ndarray) -> FieldFit:
    return FieldFit(
        model,
        factor,
        section,
        np.asarray(loading, dtype=float),
        np.asarray(field, dtype=float),
        np.ones(len(field), dtype=float) * 0.1,
        1.0,
        float(np.var(field)),
        input_hash="same-input",
        diagnostics={"patient_id": patient},
    )


def test_candidate_matching_is_one_to_one_and_joint_sign_invariant() -> None:
    sections = [("S1", "P1"), ("S2", "P2")]
    m1 = [_fit("m", "MNSF_01", sid, pid, [1, 0, 0], [1, 2, 3, 4]) for sid, pid in sections]
    m2 = [_fit("m", "MNSF_02", sid, pid, [0, 1, 0], [4, 3, 2, 1]) for sid, pid in sections]
    s1 = [_fit("s", "SIGNED_01", sid, pid, [0, -1, 0], -np.asarray([4, 3, 2, 1])) for sid, pid in sections]
    s2 = [_fit("s", "SIGNED_02", sid, pid, [-1, 0, 0], -np.asarray([1, 2, 3, 4])) for sid, pid in sections]
    candidates = match_model_factors(m1 + m2, s1 + s2, loading_threshold=0.8, field_threshold=0.8)
    assert [candidate.tier for candidate in candidates] == ["BOTH", "BOTH"]
    assert {candidate.member_factor_ids for candidate in candidates} == {
        ("MNSF_01", "SIGNED_02"),
        ("MNSF_02", "SIGNED_01"),
    }
    assert all(candidate.diagnostics["sign_flip"] == -1 for candidate in candidates)
    assert all(candidate.diagnostics["assignment_margin"] >= 0.9 for candidate in candidates)
