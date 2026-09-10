import json
from pathlib import Path

from scripts.r04_finalize_phase import main as gate_main


def _write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))
    return path


def _closure(root):
    return _write(root / "k.json", {
        "schema": "r04.k_closure_canonical.v1",
        "status": "K_CLOSED_WORKING_K3_RANK2_READOUT"})


def _unnamed(root):
    return _write(root / "u.json", {
        "schema": "r04.unnamed_field_audit.v1",
        "candidates": [
            {"candidate": "rank1_subspace", "verdict": "SURVIVES"},
            {"candidate": "rank2_subspace", "verdict": "NOT_IDENTIFIABLE"},
            {"candidate": "rank3_descriptive", "verdict": "DESCRIPTIVE_ONLY"}]})


def _comp(root, hits, patients=None):
    # hits: {(fold, rank, idx): (delta, q025)} for TSB; TLS always null.
    folds = []
    for fold in (0, 4):
        pair = patients[fold] if patients else (f"P{fold}a", f"P{fold}b")
        ranks = []
        for rank in ("rank1", "rank2"):
            inner = []
            for i, patient in enumerate(pair):
                delta, q025 = hits.get((fold, rank, i), (0.0, -0.1))
                inner.append({
                    "status": "COMPUTED",
                    "test_patients": [patient],
                    "adjusted": {"auc_delta": {"TLS": 0.0,
                                               "TUMOR_STROMA_BOUNDARY": delta}},
                    "adjusted_bootstrap_ci": {"q025_q975": {
                        "TLS": [-0.1, 0.1],
                        "TUMOR_STROMA_BOUNDARY": [q025, q025 + 0.2]}}})
            ranks.append({"rank": rank, "paired_patients": list(pair),
                          "inner_folds": inner})
        folds.append({"fold": fold, "ranks": ranks})
    return _write(root / "c.json", {
        "schema": "r04.composition_final_audit.v1",
        "status": "NESTED_AUDIT_COMPLETE", "folds": folds})


def _outer(root, name):
    return _write(root / name, {
        "schema": "r04.explore_outer_validation.v1", "status": "EXPLORATORY_COMPLETE",
        "entries": [{"rank": 1}, {"rank": 2}, {"rank": 3}]})


def _run_gate(tmp_path, monkeypatch, comp_path, label):
    import sys
    k = _closure(tmp_path)
    u = _unnamed(tmp_path)
    o0 = _outer(tmp_path, "o0.json")
    o4 = _outer(tmp_path, "o4.json")
    out = tmp_path / "gate.json"
    argv = ["r04_finalize_phase",
            "--project-root", str(tmp_path),
            "--k-closure", str(k), "--composition-audit", str(comp_path),
            "--unnamed-audit", str(u), "--output", str(out)]
    monkeypatch.setattr(sys, "argv", argv)
    import scripts.r04_finalize_phase as gate_mod
    monkeypatch.setattr(gate_mod, "OUTER_FILES",
                        {0: "o0.json", 4: "o4.json"})
    assert gate_main() == 0
    return json.loads(out.read_text())


def test_negative_outcome_closes_r04(tmp_path, monkeypatch):
    gate = _run_gate(tmp_path, monkeypatch, _comp(tmp_path, {}), "neg")
    assert gate["status"] == "R04_COMPLETE_NO_REPRODUCIBLE_RESIDUAL_FIELD"
    assert gate["surviving_candidates"] == []
    assert gate["next_nodes_authorized"] == []
    assert "R-05:NOT_TRIGGERED_NO_SURVIVING_R04_FIELD" in gate["blocked_or_not_triggered_nodes"]
    assert gate["coordinate_null"] == "NOT_TRIGGERED"
    assert gate["independent_lineage"] == "NOT_TRIGGERED_NO_SURVIVING_CANDIDATE"


def test_surviving_candidate_requires_both_folds(tmp_path, monkeypatch):
    # Stable residual in fold 0 only -> tested but unreplicated.
    comp = _comp(tmp_path, {(0, "rank2", 0): (0.15, 0.05)})
    gate = _run_gate(tmp_path, monkeypatch, comp, "one")
    assert gate["status"] == "R04_COMPLETE_NO_REPRODUCIBLE_RESIDUAL_FIELD"
    assert gate["surviving_candidates"] == []
    comp2 = _comp(tmp_path, {(0, "rank2", 0): (0.15, 0.05),
                             (4, "rank2", 1): (0.12, 0.02)})
    gate2 = _run_gate(tmp_path, monkeypatch, comp2, "two")
    assert gate2["status"] == "R04_COMPLETE_WITH_REPRODUCIBLE_RESIDUAL_FIELD"
    assert gate2["surviving_candidates"] == ["TUMOR_STROMA_BOUNDARY_rank2_residual"]
    assert gate2["next_nodes_authorized"] == ["R-05"]


def test_same_patient_cannot_replicate_with_itself(tmp_path, monkeypatch):
    comp = _comp(tmp_path, {(0, "rank2", 1): (0.15, 0.05),
                            (4, "rank2", 1): (0.12, 0.02)},
                 patients={0: ("P0a", "SHARED"), 4: ("P4a", "SHARED")})
    gate = _run_gate(tmp_path, monkeypatch, comp, "self")
    assert gate["surviving_candidates"] == []
    assert gate["status"] == "R04_COMPLETE_NO_REPRODUCIBLE_RESIDUAL_FIELD"


def test_trivial_effect_below_floor_does_not_survive(tmp_path, monkeypatch):
    comp = _comp(tmp_path, {(0, "rank2", 0): (0.008, 0.002),
                            (4, "rank2", 0): (0.006, 0.001)})
    gate = _run_gate(tmp_path, monkeypatch, comp, "tiny")
    assert gate["surviving_candidates"] == []
    assert gate["status"] == "R04_COMPLETE_NO_REPRODUCIBLE_RESIDUAL_FIELD"


def test_untested_candidate_is_not_identifiable(tmp_path, monkeypatch):
    import json
    comp_path = _comp(tmp_path, {})
    comp = json.loads(comp_path.read_text())
    for fold in comp["folds"]:
        for rank in fold["ranks"]:
            for inner in rank["inner_folds"]:
                inner["status"] = "NOT_TESTABLE_NO_TARGET_VARIATION"
                inner.pop("adjusted", None)
                inner.pop("adjusted_bootstrap_ci", None)
    comp_path.write_text(json.dumps(comp))
    gate = _run_gate(tmp_path, monkeypatch, comp_path, "untested")
    assert gate["status"] == "R04_BLOCKED_IDENTIFIABILITY" 


def test_missing_artifact_fails_closed(tmp_path, monkeypatch):
    import sys
    argv = ["r04_finalize_phase", "--project-root", str(tmp_path),
            "--k-closure", str(tmp_path / "absent.json"),
            "--composition-audit", str(tmp_path / "c.json"),
            "--unnamed-audit", str(tmp_path / "u.json"),
            "--output", str(tmp_path / "gate.json")]
    monkeypatch.setattr(sys, "argv", argv)
    try:
        gate_main()
    except ValueError as exc:
        assert "missing evidence artifact" in str(exc)
    else:
        raise AssertionError("missing artifact must fail closed")
