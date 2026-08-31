from __future__ import annotations

from scripts.r04_late_lr_continuation import _diagnostic_payload


def test_continuation_diagnostic_preserves_requested_factor_count() -> None:
    payload = _diagnostic_payload(
        factors=3,
        fold=0,
        steps=7200,
        expected_start_step=4800,
        learning_rate=0.01,
        learning_rate_final=None,
        learning_rate_decay_steps=None,
        gene_batch_size=512,
        evaluation_interval=300,
        evaluation_mc_draws=4,
        fit_converged=False,
        source_checkpoint_dir="source/k3/fold0",
        restart_index=2,
        optimization_seed=12345,
    )

    assert payload["k_model"] == 3
    assert payload["source_checkpoint_dir"] == "source/k3/fold0"
    assert payload["wiring_only"] is True
    assert payload["restart_index"] == 2
    assert payload["optimization_seed"] == 12345
