import copy
import hashlib
from dataclasses import replace

import numpy as np
import pytest
import torch

from r04.checkpoint_audit import evaluate_checkpoint, nuisance_gauge_summary
from r04.models import MNSFConfig, MNSFEstimator
from scripts.r04_checkpoint_audit import _audit_payload
from r04.synthetic import make_overlapping_sections


def test_nuisance_gauge_shift_preserves_additive_rate() -> None:
    baseline = np.asarray([2.0, 3.0, 5.0])
    loading = np.asarray([[0.2], [0.3], [0.5]])
    nuisance = np.asarray([[1.5], [2.0], [3.0]])
    shift = 1.0
    before = baseline[None, :] + nuisance @ loading.T
    after = (
        (baseline + shift * loading[:, 0])[None, :]
        + (nuisance - shift) @ loading.T
    )
    np.testing.assert_allclose(before, after)

    summary = nuisance_gauge_summary(
        np.log(np.expm1(baseline)),
        np.log(loading),
        np.log(np.expm1(nuisance)),
    )
    assert summary["rank"] == 1
    assert summary["h_min"][0] == pytest.approx(1.50001, abs=1e-6)


def test_checkpoint_audit_payload_records_ridge() -> None:
    payload = _audit_payload(
        manifest_hash="manifest",
        gene_count=4000,
        factors=3,
        fold=4,
        folds=5,
        group_seed=20260807,
        inducing_points=16,
        lengthscale=3.0,
        ridge=1e-5,
        mc_draws=4,
        checkpoints=[],
    )

    assert payload["schema"] == "r04.checkpoint_audit.v1"
    assert payload["ridge"] == 1e-5


def test_torch_checkpoint_audit_reads_v3_parameter_layout_and_p_only_fallback(tmp_path) -> None:
    sections, _ = make_overlapping_sections(
        sections=2, spots_per_section=8, genes=6, antagonistic=True
    )
    checkpoint_dir = tmp_path / "checkpoint"
    MNSFEstimator(
        MNSFConfig(
            factors=2,
            inducing_points=3,
            steps=2,
            posterior_draws=1,
            checkpoint_dir=str(checkpoint_dir),
            checkpoint_steps=1,
            seed=127,
        )
    ).fit(sections)

    checkpoint_path = checkpoint_dir / "checkpoint.pt"
    before_sha = hashlib.sha256(checkpoint_path.read_bytes()).hexdigest()
    audited = evaluate_checkpoint(
        sections,
        checkpoint_dir,
        factors=2,
        inducing_points=3,
        lengthscale=2.0,
        ridge=1e-5,
        mc_draws=2,
        seed=129,
    )
    assert audited["global_step"] == 2
    assert audited["n_spots"] == 16
    assert audited["n_genes"] == 6
    assert np.isfinite(audited["objective_mean"])
    assert audited["parameter_state_source"] == "params"
    assert audited["objective_input_hash_status"] == "VERIFIED"
    assert hashlib.sha256(checkpoint_path.read_bytes()).hexdigest() == before_sha

    payload = torch.load(
        str(checkpoint_path), map_location="cpu", weights_only=True
    )
    payload.pop("parameter_layout")
    fallback_path = tmp_path / "checkpoint_without_layout.pt"
    torch.save(payload, str(fallback_path))
    fallback = evaluate_checkpoint(
        sections,
        fallback_path,
        factors=2,
        inducing_points=3,
        lengthscale=2.0,
        ridge=1e-5,
        mc_draws=2,
        seed=129,
    )
    assert fallback["global_step"] == audited["global_step"]
    assert fallback["objective_mean"] == pytest.approx(
        audited["objective_mean"], abs=1e-6
    )

    legacy = copy.deepcopy(payload)
    legacy.pop("objective_input_hash")
    legacy_path = tmp_path / "legacy_v3_without_objective_hash.pt"
    torch.save(legacy, str(legacy_path))
    legacy_audit = evaluate_checkpoint(
        sections,
        legacy_path,
        factors=2,
        inducing_points=3,
        lengthscale=2.0,
        ridge=1e-5,
        mc_draws=1,
        seed=129,
    )
    assert legacy_audit["objective_input_hash_status"] == "NOT_VERIFIED_LEGACY_V3"


def test_torch_checkpoint_audit_rejects_legacy_tensorflow_state(tmp_path) -> None:
    legacy_path = tmp_path / "legacy_tf.pt"
    torch.save(
        {"checkpoint_schema": "r04.mnsf_checkpoint.v2", "backend": "tensorflow"},
        str(legacy_path),
    )
    sections, _ = make_overlapping_sections(
        sections=1, spots_per_section=8, genes=6, antagonistic=True
    )
    with pytest.raises(ValueError, match="legacy TensorFlow"):
        evaluate_checkpoint(
            sections,
            legacy_path,
            factors=2,
            inducing_points=3,
            lengthscale=2.0,
            ridge=1e-5,
            mc_draws=1,
        )


def test_torch_checkpoint_audit_never_fills_params_from_best_state(tmp_path) -> None:
    sections, _ = make_overlapping_sections(
        sections=2, spots_per_section=8, genes=6, antagonistic=True
    )
    checkpoint_dir = tmp_path / "checkpoint"
    MNSFEstimator(
        MNSFConfig(
            factors=2,
            inducing_points=3,
            steps=2,
            posterior_draws=1,
            checkpoint_dir=str(checkpoint_dir),
            checkpoint_steps=1,
            seed=131,
        )
    ).fit(sections)
    payload = torch.load(
        str(checkpoint_dir / "checkpoint.pt"), map_location="cpu", weights_only=True
    )
    payload["params"].pop("p_0")
    broken = tmp_path / "missing_params.pt"
    torch.save(payload, str(broken))
    with pytest.raises(ValueError, match="missing (?:parameter )?tensor p_0"):
        evaluate_checkpoint(
            sections,
            broken,
            factors=2,
            inducing_points=3,
            lengthscale=2.0,
            ridge=1e-5,
            mc_draws=1,
        )


def test_torch_checkpoint_audit_rejects_malformed_layout_and_one_sided_nuisance(tmp_path) -> None:
    sections, _ = make_overlapping_sections(
        sections=2, spots_per_section=8, genes=6, antagonistic=True
    )
    checkpoint_dir = tmp_path / "checkpoint"
    MNSFEstimator(
        MNSFConfig(
            factors=2,
            inducing_points=3,
            steps=2,
            posterior_draws=1,
            checkpoint_dir=str(checkpoint_dir),
            checkpoint_steps=1,
            seed=137,
        )
    ).fit(sections)
    original = torch.load(
        str(checkpoint_dir / "checkpoint.pt"), map_location="cpu", weights_only=True
    )
    malformed = copy.deepcopy(original)
    malformed["parameter_layout"]["q_loc"] = ["p_4"]
    malformed_path = tmp_path / "malformed_layout.pt"
    torch.save(malformed, str(malformed_path))
    with pytest.raises(ValueError, match="one key per section"):
        evaluate_checkpoint(
            sections,
            malformed_path,
            factors=2,
            inducing_points=3,
            lengthscale=2.0,
            ridge=1e-5,
            mc_draws=1,
        )

    one_sided = copy.deepcopy(original)
    one_sided["parameter_layout"]["raw_h"] = None
    one_sided_path = tmp_path / "one_sided_nuisance.pt"
    torch.save(one_sided, str(one_sided_path))
    with pytest.raises(ValueError, match="raw_v and raw_h together"):
        evaluate_checkpoint(
            sections,
            one_sided_path,
            factors=2,
            inducing_points=3,
            lengthscale=2.0,
            ridge=1e-5,
            mc_draws=1,
        )


def test_torch_checkpoint_audit_hashes_explicit_library_size(tmp_path) -> None:
    sections, _ = make_overlapping_sections(
        sections=2, spots_per_section=8, genes=6, antagonistic=True
    )
    checkpoint_dir = tmp_path / "checkpoint"
    MNSFEstimator(
        MNSFConfig(
            factors=2,
            inducing_points=3,
            steps=2,
            posterior_draws=1,
            checkpoint_dir=str(checkpoint_dir),
            checkpoint_steps=1,
            seed=139,
        )
    ).fit(sections)
    changed = [
        replace(section, library_size=np.full(len(section.coords), 2.0))
        for section in sections
    ]
    with pytest.raises(ValueError, match="objective input hash"):
        evaluate_checkpoint(
            changed,
            checkpoint_dir,
            factors=2,
            inducing_points=3,
            lengthscale=2.0,
            ridge=1e-5,
            mc_draws=1,
        )
