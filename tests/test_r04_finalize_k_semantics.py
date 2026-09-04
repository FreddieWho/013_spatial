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
