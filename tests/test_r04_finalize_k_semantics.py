from scripts.r04_finalize_k_semantics import build_final_artifact


def test_final_handoff_links_scores_subspaces_and_current_bridge_block() -> None:
    aggregate = {
        "schema": "r04.k_semantics_and_downstream_robustness.v1",
        "k_model": {"structure_mapping": "many_to_many"},
        "k_eff": {"value": None},
        "selected_k": None,
        "available_evidence": {
            "score_summary": {"cell_count": 20},
            "subspace_summary_by_fold": {"0": {}, "4": {}},
        },
    }
    entries = []
    for fold in (0, 4):
        for restart in range(5):
            entries.append({
                "fold": fold,
                "restart_index": restart,
                "coverage_training": {"paired": {"patients": ["p1", "p2"]}},
                "full_k3": {"inner_crossfit": {"status": "COMPUTED"}},
                "stable_rank2_sensitivity": {"inner_crossfit": {"status": "COMPUTED"}},
            })
    structure = {
        "schema": "r04.structure_readout_panel.v1",
        "status": "COMPLETE_DESCRIPTIVE_NOT_FORMAL_STRUCTURE_CLAIM",
        "scope": "training_role_inner_crossfit_only",
        "downstream_k_robust": "not_tested",
        "k2_bridge_required": True,
        "k2_bridge_reason": "insufficient paired patients",
        "entries": entries,
        "internal_external_validation_gt": "SEALED_NOT_READ",
    }
    stage = {
        "schema": "r04.k2_bridge_stage_a.v1",
        "status": "BLOCKED_CPU_RUNTIME",
        "blocked_reason": "stopped before checkpoint",
        "attempted_cell": {"status": "STOPPED_BEFORE_CHECKPOINT"},
        "gpu_required_for_bounded_completion": True,
    }
    manifest = {
        "schema": "r04.k2_bridge_stage_a_manifest.v1",
        "status": "FROZEN_PRECOMPUTE_MANIFEST",
        "entries": [{"fold": 0}],
    }
    result = build_final_artifact(
        aggregate,
        structure,
        stage,
        manifest,
        source_artifacts=[{"path": "x.json", "sha256": "a" * 64}],
    )
    assert result["schema"] == "r04.k_semantics_and_downstream_robustness.v2"
    assert result["selected_k"] is None
    assert result["score_summary"]["cell_count"] == 20
    assert set(result["subspace_stability_by_fold"]) == {"0", "4"}
    assert result["structure_readout"]["downstream_k_robust"] == "not_tested"
    assert result["k2_bridge"]["stage_status"] == "BLOCKED_CPU_RUNTIME"
    assert result["k2_bridge"]["scientific_result"].startswith("NOT_AVAILABLE")


def test_final_handoff_rejects_completed_k2_as_current_block() -> None:
    base = {
        "schema": "r04.k_semantics_and_downstream_robustness.v1",
        "k_model": {},
        "k_eff": {},
        "available_evidence": {"score_summary": {}, "subspace_summary_by_fold": {}},
    }
    structure = {
        "schema": "r04.structure_readout_panel.v1",
        "entries": [{}] * 10,
        "internal_external_validation_gt": "SEALED_NOT_READ",
    }
    stage = {"schema": "r04.k2_bridge_stage_a.v1", "status": "COMPLETE"}
    manifest = {"schema": "r04.k2_bridge_stage_a_manifest.v1", "status": "FROZEN_PRECOMPUTE_MANIFEST", "entries": []}
    try:
        build_final_artifact(base, structure, stage, manifest, source_artifacts=[])
    except ValueError as exc:
        assert "CPU runtime block" in str(exc)
    else:
        raise AssertionError("completed K2 stage must not be relabeled as blocked")


def _closure_cell(fold, k, scores, in_hash="044d28c5" + "0" * 56,
                  obj_hash="b5a3fb68" + "0" * 56):
    return {
        "status": "FIT_AND_SCORED",
        "k_model": k,
        "fold": fold,
        "input_hash": in_hash,
        "inference_platform": {"converged": True, "split_converged": [True, True]},
        "fit_diagnostics": {"steps": 8400, "objective_input_hash": obj_hash},
        "patient_scores": dict(scores),
    }


def _closure_runs(base0, base4):
    # Deltas reproduce the recorded D-106 means up to the binding tolerance.
    k02_f0 = {"p1": 0.0, "p2": 0.0, "p3": 0.0, "p4": -99.63, "p5": 0.0, "p6": 0.0}
    k02_f4 = {"p1": 7.27, "p2": 186.82, "p3": 51.16, "p4": 94.93, "p5": 16.64, "p6": 74.94}
    k3_f0 = {"p1": 0.0, "p2": 0.0, "p3": 0.0, "p4": -173.78, "p5": 0.0, "p6": 0.0}
    k3_f4 = {"p1": 11.89, "p2": 186.82, "p3": 51.16, "p4": 94.93, "p5": 16.64, "p6": 74.94}
    f4_in, f4_obj = "0595b0ae" + "1" * 56, "593d99cd" + "1" * 56
    runs = {}
    for tag, f0d, f4d in (("k02", k02_f0, k02_f4), ("k3", k3_f0, k3_f4)):
        hi = 2 if tag == "k02" else 3
        k0_f0 = {p: 0.0 for p in f0d}
        k0_f4 = {p: 0.0 for p in f4d}
        runs[tag] = {
            0: {"top": {"status": "WIRING_ONLY_NOT_SCIENTIFIC",
                        "manifest_hash": "5a0dc1de" + "2" * 56, "selected_k": None},
                "cells": {0: _closure_cell(0, 0, k0_f0),
                          hi: _closure_cell(0, hi, {p: k0_f0[p] + f0d[p] for p in f0d})}},
            4: {"top": {"status": "WIRING_ONLY_NOT_SCIENTIFIC",
                        "manifest_hash": "5a0dc1de" + "2" * 56, "selected_k": None},
                "cells": {0: _closure_cell(4, 0, k0_f4, in_hash=f4_in, obj_hash=f4_obj),
                          hi: _closure_cell(4, hi, {p: k0_f4[p] + f4d[p] for p in f4d},
                                            in_hash=f4_in, obj_hash=f4_obj)}},
        }
    return runs, f4_in, f4_obj


def _closure_audits(f4_in, f4_obj):
    def audit(fold, factors, in_hash, obj_hash):
        return {"fold": fold, "factors": factors, "checkpoints": [{
            "global_step": 8400, "objective_input_hash_status": "VERIFIED",
            "input_hash": in_hash, "objective_input_hash": obj_hash}]}
    return [audit(0, 0, "044d28c5" + "0" * 56, "b5a3fb68" + "0" * 56),
            audit(0, 2, "044d28c5" + "0" * 56, "b5a3fb68" + "0" * 56),
            audit(4, 0, f4_in, f4_obj), audit(4, 2, f4_in, f4_obj)]


def test_k_closure_recomputes_recorded_means() -> None:
    from scripts.r04_finalize_k_semantics import build_k_closure
    runs, f4_in, f4_obj = _closure_runs(None, None)
    result = build_k_closure(runs, _closure_audits(f4_in, f4_obj), source_artifacts=[])
    assert result["schema"] == "r04.k_closure_canonical.v1"
    assert result["status"] == "K_CLOSED_WORKING_K3_RANK2_READOUT"
    assert result["working_k_model"] == 3
    assert result["primary_readout_rank"] == 2
    assert result["global_k_eff"] == "NOT_IDENTIFIABLE"
    assert result["selected_k"] is None
    assert result["k_search_closed_for_r04"] is True
    assert abs(result["paired_means"]["k02_minus_k0_fold0"] + 16.61) < 0.005
    assert abs(result["paired_means"]["k02_minus_k0_fold4"] - 71.96) < 0.005


def test_k_closure_rejects_score_drift() -> None:
    from scripts.r04_finalize_k_semantics import build_k_closure
    runs, f4_in, f4_obj = _closure_runs(None, None)
    runs["k02"][4]["cells"][2]["patient_scores"] = dict(
        runs["k02"][4]["cells"][2]["patient_scores"])
    runs["k02"][4]["cells"][2]["patient_scores"]["p2"] = -500.0
    try:
        build_k_closure(runs, _closure_audits(f4_in, f4_obj), source_artifacts=[])
    except ValueError as exc:
        assert "drifts from D-106" in str(exc)
    else:
        raise AssertionError("tampered scores must fail closed")


def test_k_closure_rejects_fold4_hash_split() -> None:
    from scripts.r04_finalize_k_semantics import build_k_closure
    runs, f4_in, f4_obj = _closure_runs(None, None)
    runs["k3"][4]["cells"][3]["input_hash"] = "deadbeef" + "3" * 56
    try:
        build_k_closure(runs, _closure_audits(f4_in, f4_obj), source_artifacts=[])
    except ValueError as exc:
        assert "fold-4 cells disagree" in str(exc)
    else:
        raise AssertionError("split fold-4 hashes must fail closed")
