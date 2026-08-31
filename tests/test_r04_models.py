from __future__ import annotations

import pytest
import numpy as np

from r04.models import MNSFConfig, MNSFEstimator, SignedResidualGPConfig, SignedResidualGPEstimator
from r04.models.mnsf import predict_log_mean_from_fields
from r04.models.mnsf import _diagonal_gp_kl, _rate_initial_components
from r04.serialization import read_frozen_model, write_frozen_model
from r04.synthetic import make_overlapping_sections


tf = pytest.importorskip("tensorflow")
pytest.importorskip("tensorflow_probability")


def test_diagonal_gp_kl_matches_identity_prior() -> None:
    prior_inverse = tf.eye(3)
    loc = tf.zeros((3, 2))
    unit = _diagonal_gp_kl(tf, loc, tf.ones((3, 2)), prior_inverse, tf.constant(0.0))
    wider = _diagonal_gp_kl(tf, loc, tf.fill((3, 2), 2.0), prior_inverse, tf.constant(0.0))
    assert float(unit.numpy()) == pytest.approx(0.0, abs=1e-6)
    assert float(wider.numpy()) > 0.0


def test_dual_models_emit_continuous_uncertain_fields() -> None:
    sections, _ = make_overlapping_sections(sections=3, spots_per_section=25, genes=18, antagonistic=True)
    mnsf = MNSFEstimator(MNSFConfig(factors=2, inducing_points=8, steps=4, posterior_draws=8, seed=17)).fit(sections)
    signed = SignedResidualGPEstimator(SignedResidualGPConfig(factors=2, inducing_points=8, steps=4, posterior_draws=8, seed=19)).fit(sections)
    assert mnsf.loading_ is not None and mnsf.loading_.shape == (18, 2)
    assert signed.loading_ is not None and signed.loading_.shape == (18, 2)
    assert len(mnsf.fields_) == len(sections)
    assert len(signed.fields_) == len(sections)
    for fields in [*mnsf.fields_, *signed.fields_]:
        for field in fields:
            assert field.field_mean.shape == field.field_sd.shape
            assert field.field_mean.ndim == 1
            assert field.diagnostics["uncertainty_kind"] == "variational_inducing_weights"


def test_frozen_models_infer_new_sections_without_refitting(tmp_path) -> None:
    sections, _ = make_overlapping_sections(
        sections=3, spots_per_section=20, genes=12, antagonistic=True
    )
    estimators = [
        MNSFEstimator(MNSFConfig(factors=2, inducing_points=6, steps=3, posterior_draws=4, seed=31)),
        SignedResidualGPEstimator(SignedResidualGPConfig(factors=2, inducing_points=6, steps=3, posterior_draws=4, seed=37)),
    ]
    for index, estimator in enumerate(estimators):
        estimator.fit(sections[:2])
        path = tmp_path / f"model_{index}.json"
        write_frozen_model(path, estimator)
        restored = read_frozen_model(path)
        before = estimator.loading_.copy()
        inferred = restored.infer(sections[2:], steps=2, posterior_draws=3)
        assert len(inferred) == 1
        assert len(inferred[0]) == 2
        assert restored.fit_input_hash_ != ""
        assert restored.loading_.shape == before.shape
        np.testing.assert_allclose(restored.loading_, before)
        for field in inferred[0]:
            assert field.field_mean.shape == field.field_sd.shape == (20,)
            assert field.diagnostics["fit_input_hash"] == restored.fit_input_hash_
            assert field.diagnostics["inference_steps"] == 2


def test_legacy_frozen_models_fail_closed_after_gp_parameterization_change() -> None:
    sections, _ = make_overlapping_sections(
        sections=1, spots_per_section=12, genes=8, antagonistic=True
    )
    for estimator in [
        MNSFEstimator(MNSFConfig(factors=2, inducing_points=5, steps=2, posterior_draws=2, seed=47)),
        SignedResidualGPEstimator(SignedResidualGPConfig(factors=2, inducing_points=5, steps=2, posterior_draws=2, seed=53)),
    ]:
        estimator.fit(sections)
        state = estimator.frozen_state()
        state["config"].pop("gp_parameterization")
        with pytest.raises(ValueError, match="D-050 GP parameterization"):
            type(estimator).from_frozen_state(state)


def test_mnsf_field_prediction_is_continuous_and_gene_aligned() -> None:
    sections, _ = make_overlapping_sections(sections=2, spots_per_section=12, genes=8, antagonistic=True)
    estimator = MNSFEstimator(MNSFConfig(factors=2, inducing_points=5, steps=3, posterior_draws=3, seed=73)).fit(sections[:1])
    inferred = estimator.infer(sections[1:], steps=2, posterior_draws=2)
    log_mu = predict_log_mean_from_fields(
        sections[1:], inferred, estimator.loading_, estimator.factor_amplitude_, estimator.gene_baseline_
    )
    assert log_mu.shape == (12, 8)
    assert np.isfinite(log_mu).all()


def test_mnsf_compiled_execution_matches_eager_reference_on_small_fixture() -> None:
    sections, _ = make_overlapping_sections(
        sections=2, spots_per_section=12, genes=8, antagonistic=True
    )
    common = dict(
        factors=2,
        inducing_points=5,
        steps=3,
        posterior_draws=2,
        gene_batch_size=4,
        diagnostic_interval=1,
        seed=211,
    )
    eager = MNSFEstimator(MNSFConfig(**common, execution_mode="eager")).fit(sections)
    compiled = MNSFEstimator(
        MNSFConfig(**common, execution_mode="compiled")
    ).fit(sections)
    np.testing.assert_allclose(
        compiled.loading_, eager.loading_, rtol=1e-5, atol=1e-6
    )
    np.testing.assert_allclose(
        compiled.factor_amplitude_, eager.factor_amplitude_, rtol=1e-5, atol=1e-6
    )
    np.testing.assert_allclose(
        compiled.gene_baseline_, eager.gene_baseline_, rtol=1e-5, atol=1e-6
    )
    np.testing.assert_allclose(
        compiled.diagnostics_["loss_trace"],
        eager.diagnostics_["loss_trace"],
        rtol=1e-5,
        atol=1e-6,
    )
    assert compiled.diagnostics_["execution_mode"] == "compiled"

    inferred_eager = eager.infer(sections, steps=2, posterior_draws=2)
    inferred_compiled = compiled.infer(sections, steps=2, posterior_draws=2)
    np.testing.assert_allclose(
        inferred_compiled[0][0].field_mean,
        inferred_eager[0][0].field_mean,
        rtol=1e-5,
        atol=1e-6,
    )
    np.testing.assert_allclose(
        inferred_compiled[1][1].field_mean,
        inferred_eager[1][1].field_mean,
        rtol=1e-5,
        atol=1e-6,
    )
    assert compiled.inference_diagnostics_["execution_mode"] == "compiled"


def test_mnsf_compiled_execution_preserves_k0_and_staged_schedule() -> None:
    sections, _ = make_overlapping_sections(
        sections=2, spots_per_section=10, genes=7, antagonistic=True
    )
    estimator = MNSFEstimator(MNSFConfig(
        factors=0,
        nonspatial_rank=1,
        steps=3,
        posterior_draws=1,
        diagnostic_interval=1,
        optimization_schedule="staged_shared_first",
        shared_steps=2,
        execution_mode="compiled",
        seed=223,
    )).fit(sections)
    assert estimator.fields_ == [[], []]
    assert estimator.diagnostics_["execution_mode"] == "compiled"
    assert len(estimator.diagnostics_["loss_trace"]) == 3
    inferred = estimator.infer(sections, steps=2, posterior_draws=1)
    assert inferred == [[], []]
    assert estimator.inference_diagnostics_["execution_mode"] == "compiled"


def test_mnsf_k0_uses_same_nuisance_channel_and_emits_no_spatial_fields() -> None:
    sections, _ = make_overlapping_sections(
        sections=2, spots_per_section=14, genes=9, antagonistic=True
    )
    estimator = MNSFEstimator(MNSFConfig(
        factors=0,
        nonspatial_rank=1,
        inducing_points=5,
        steps=2,
        posterior_draws=2,
        seed=79,
    )).fit(sections[:1])
    assert estimator.loading_ is not None
    assert estimator.loading_.shape == (9, 0)
    assert estimator.factor_amplitude_.shape == (0,)
    assert estimator.nonspatial_loading_ is not None
    assert all(group == [] for group in estimator.fields_)

    inferred = estimator.infer(sections[1:], steps=2, posterior_draws=2)
    assert inferred == [[]]
    assert estimator.inference_nonspatial_component_ is not None
    log_mu = predict_log_mean_from_fields(
        sections[1:],
        inferred,
        estimator.loading_,
        estimator.factor_amplitude_,
        estimator.gene_baseline_,
        nonspatial_component=estimator.inference_nonspatial_component_,
        nonspatial_loading=estimator.nonspatial_loading_,
    )
    assert log_mu.shape == (14, 9)
    assert np.isfinite(log_mu).all()


def test_frozen_models_support_explicit_common_observed_panel() -> None:
    sections, _ = make_overlapping_sections(sections=2, spots_per_section=16, genes=10, antagonistic=True)
    for estimator in [
        MNSFEstimator(MNSFConfig(factors=2, inducing_points=5, steps=2, posterior_draws=2, seed=61)),
        SignedResidualGPEstimator(SignedResidualGPConfig(factors=2, inducing_points=5, steps=2, posterior_draws=2, seed=67)),
    ]:
        estimator.fit(sections)
        observed_genes = tuple(estimator.gene_id_[:-2])
        projected = estimator.subset_to_genes(observed_genes)
        reduced = [
            type(section)(
                section.section_id,
                section.patient_id,
                section.block_id,
                section.lineage,
                section.barcode,
                section.coords,
                section.counts[:, :-2],
                observed_genes,
                section.library_size,
                section.metadata,
            )
            for section in sections
        ]
        inferred = projected.infer(reduced, steps=1, posterior_draws=2)
        assert projected.diagnostics_["conditional_panel"] is True
        assert len(inferred) == len(reduced)
        assert inferred[0][0].loading.shape == (8,)


def test_checkpoint_resume_restores_effective_step(tmp_path) -> None:
    sections, _ = make_overlapping_sections(sections=2, spots_per_section=12, genes=8, antagonistic=True)
    configs = [
        MNSFConfig(factors=2, inducing_points=4, steps=2, posterior_draws=1, seed=71, checkpoint_dir=str(tmp_path / "mnsf")),
        SignedResidualGPConfig(factors=2, inducing_points=4, steps=2, posterior_draws=1, seed=73, checkpoint_dir=str(tmp_path / "signed")),
    ]
    for config in configs:
        estimator_type = MNSFEstimator if isinstance(config, MNSFConfig) else SignedResidualGPEstimator
        first = estimator_type(config).fit(sections)
        resumed = estimator_type(config).fit(sections)
        assert first.loading_ is not None and resumed.loading_ is not None
        assert resumed.fields_
        assert resumed.diagnostics_["steps"] == 2


def test_mnsf_compiled_checkpoint_resume_preserves_optimizer_state(tmp_path, monkeypatch) -> None:
    sections, _ = make_overlapping_sections(
        sections=2, spots_per_section=10, genes=7, antagonistic=True
    )
    common = dict(
        factors=2,
        inducing_points=4,
        steps=4,
        posterior_draws=1,
        checkpoint_steps=2,
        execution_mode="compiled",
        seed=227,
    )
    uninterrupted = MNSFEstimator(MNSFConfig(**common)).fit(sections)
    checkpointed = MNSFConfig(
        **common, checkpoint_dir=str(tmp_path / "compiled_mnsf")
    )
    original_apply = tf.keras.optimizers.Adam.apply_gradients
    calls = 0

    def interrupt_on_third_update(optimizer, *args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 3:
            raise RuntimeError("simulated interruption")
        return original_apply(optimizer, *args, **kwargs)

    monkeypatch.setattr(
        tf.keras.optimizers.Adam, "apply_gradients", interrupt_on_third_update
    )
    with pytest.raises(RuntimeError, match="simulated interruption"):
        MNSFEstimator(checkpointed).fit(sections)
    monkeypatch.setattr(tf.keras.optimizers.Adam, "apply_gradients", original_apply)
    resumed = MNSFEstimator(checkpointed).fit(sections)
    np.testing.assert_allclose(
        resumed.loading_, uninterrupted.loading_, rtol=1e-5, atol=1e-6
    )
    np.testing.assert_allclose(
        resumed.gene_baseline_, uninterrupted.gene_baseline_, rtol=1e-5, atol=1e-6
    )
    assert resumed.diagnostics_["best_step"] == uninterrupted.diagnostics_["best_step"]


def test_mnsf_checkpoint_resume_preserves_precrash_best_state(tmp_path, monkeypatch) -> None:
    sections, _ = make_overlapping_sections(
        sections=2, spots_per_section=12, genes=8, antagonistic=True
    )
    common = dict(
        factors=2,
        inducing_points=4,
        steps=4,
        posterior_draws=1,
        seed=79,
        checkpoint_steps=2,
    )
    uninterrupted = MNSFEstimator(MNSFConfig(**common)).fit(sections)
    checkpointed = MNSFConfig(
        **common, checkpoint_dir=str(tmp_path / "interrupted_mnsf")
    )

    original_apply = tf.keras.optimizers.Adam.apply_gradients
    calls = 0

    def interrupt_on_third_update(optimizer, *args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 3:
            raise RuntimeError("simulated interruption")
        return original_apply(optimizer, *args, **kwargs)

    monkeypatch.setattr(
        tf.keras.optimizers.Adam, "apply_gradients", interrupt_on_third_update
    )
    with pytest.raises(RuntimeError, match="simulated interruption"):
        MNSFEstimator(checkpointed).fit(sections)
    monkeypatch.setattr(tf.keras.optimizers.Adam, "apply_gradients", original_apply)

    resumed = MNSFEstimator(checkpointed).fit(sections)
    np.testing.assert_allclose(resumed.loading_, uninterrupted.loading_, atol=1e-6)
    np.testing.assert_allclose(
        resumed.factor_amplitude_, uninterrupted.factor_amplitude_, atol=1e-6
    )
    np.testing.assert_allclose(
        resumed.gene_baseline_, uninterrupted.gene_baseline_, atol=1e-6
    )
    assert resumed.diagnostics_["loss_best"] == pytest.approx(
        uninterrupted.diagnostics_["loss_best"], abs=1e-6
    )
    assert resumed.diagnostics_["best_step"] == uninterrupted.diagnostics_["best_step"]


def test_mnsf_one_step_best_state_matches_the_evaluated_preupdate_parameters() -> None:
    sections, _ = make_overlapping_sections(
        sections=2, spots_per_section=12, genes=8, antagonistic=True
    )
    config = MNSFConfig(
        factors=2,
        inducing_points=4,
        steps=1,
        posterior_draws=1,
        seed=83,
    )
    counts = np.concatenate([
        section.counts.toarray()
        if hasattr(section.counts, "toarray")
        else np.asarray(section.counts)
        for section in sections
    ], axis=0)
    library = np.concatenate([
        np.asarray(section.library_size, dtype=float)
        if section.library_size is not None
        else np.asarray(section.counts.sum(axis=1)).ravel().astype(float)
        for section in sections
    ])
    initial_loading, _, _, _ = _rate_initial_components(
        counts, library, config.factors, seed=config.seed
    )
    fitted = MNSFEstimator(config).fit(sections)
    np.testing.assert_allclose(fitted.loading_, initial_loading, atol=1e-6)
    assert fitted.diagnostics_["best_step"] == 1


def test_mnsf_continuation_preserves_optimizer_state_and_starts_at_checkpoint(
    tmp_path,
) -> None:
    sections, _ = make_overlapping_sections(
        sections=2, spots_per_section=12, genes=8, antagonistic=True
    )
    source = tmp_path / "source"
    target = tmp_path / "target"
    MNSFEstimator(
        MNSFConfig(
            factors=0,
            inducing_points=5,
            nonspatial_rank=1,
            steps=3,
            posterior_draws=2,
            checkpoint_dir=str(source),
            checkpoint_steps=1,
            seed=101,
        )
    ).fit(sections)
    continued = MNSFEstimator(
        MNSFConfig(
            factors=0,
            inducing_points=5,
            nonspatial_rank=1,
            steps=6,
            learning_rate=0.01,
            learning_rate_final=0.0025,
            learning_rate_decay_steps=3,
            posterior_draws=2,
            checkpoint_dir=str(target),
            checkpoint_steps=1,
            resume_checkpoint_dir=str(source),
            seed=101,
        )
    ).fit(sections)
    assert continued.diagnostics_["start_step"] == 3
    assert continued.diagnostics_["optimizer_steps_this_call"] == 3
    assert continued.diagnostics_["resume_checkpoint_dir"] == str(source)
    assert continued.diagnostics_["learning_rate_trace"][0]["learning_rate"] == pytest.approx(
        0.0025, abs=1e-6
    )


def test_staged_shared_first_freezes_spatial_blocks_and_preserves_k0_baseline() -> None:
    sections, _ = make_overlapping_sections(
        sections=2, spots_per_section=10, genes=7, antagonistic=True
    )
    common = dict(
        inducing_points=4,
        steps=4,
        posterior_draws=1,
        diagnostic_interval=1,
        seed=109,
    )
    joint_k0 = MNSFEstimator(MNSFConfig(factors=0, **common)).fit(sections)
    staged_k0 = MNSFEstimator(MNSFConfig(
        factors=0,
        optimization_schedule="staged_shared_first",
        shared_steps=2,
        **common,
    )).fit(sections)
    np.testing.assert_allclose(
        staged_k0.diagnostics_["loss_trace"],
        joint_k0.diagnostics_["loss_trace"],
        atol=1e-7,
    )
    assert all(
        item["optimization_stage"] == "joint"
        for item in staged_k0.diagnostics_["optimization_stage_trace"]
    )

    staged_k2 = MNSFEstimator(MNSFConfig(
        factors=2,
        optimization_schedule="staged_shared_first",
        shared_steps=2,
        **common,
    )).fit(sections)
    assert [
        item["optimization_stage"]
        for item in staged_k2.diagnostics_["optimization_stage_trace"]
    ] == ["shared_first", "shared_first", "joint", "joint"]
    first, second = staged_k2.diagnostics_["parameter_norm_trace"][:2]
    assert second["loading"] == pytest.approx(first["loading"], abs=1e-7)
    assert second["gp_loc"] == pytest.approx(first["gp_loc"], abs=1e-7)
    assert second["gp_scale"] == pytest.approx(first["gp_scale"], abs=1e-7)


def test_mnsf_emits_dense_fixed_objective_trace_without_changing_fit_contract() -> None:
    sections, _ = make_overlapping_sections(
        sections=2, spots_per_section=8, genes=6, antagonistic=True
    )
    estimator = MNSFEstimator(MNSFConfig(
        factors=0,
        inducing_points=3,
        steps=2,
        posterior_draws=1,
        evaluation_interval=1,
        evaluation_mc_draws=2,
        seed=113,
    )).fit(sections)
    trace = estimator.diagnostics_["evaluation_trace"]
    assert len(trace) == 2
    assert all(item["objective_kind"] == "dense_full_panel_fixed_stateless_mc" for item in trace)
    assert all(item["objective_sd"] == pytest.approx(0.0, abs=1e-7) for item in trace)
    assert estimator.diagnostics_["deterministic_objective_platform"]["converged"] is False
