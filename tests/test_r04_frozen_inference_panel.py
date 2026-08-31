from __future__ import annotations

import json
from pathlib import Path

from scripts.r04_frozen_inference_panel import (
    _cell_specs,
    _child_environment,
    _completed_output,
    _resume_command,
)


def test_frozen_panel_worker_has_four_paired_cells_and_no_fit_extension() -> None:
    specs = _cell_specs(3, Path("out"), Path("source"))
    assert [(item["factors"], item["fold"]) for item in specs] == [
        (0, 0), (3, 0), (0, 4), (3, 4)
    ]
    command = _resume_command(
        python="python",
        repo=Path("/repo"),
        restart_index=3,
        spec=specs[1],
        panel_audit=Path("panel.json"),
        acceptance=Path("acceptance.json"),
    )
    assert command[command.index("--fit-steps") + 1] == "8400"
    assert command[command.index("--expected-start-step") + 1] == "8400"
    assert command[command.index("--inference-steps") + 1] == "2400"
    assert command[command.index("--manifest-json") + 1].endswith(
        "infra/r04/role_manifests/training_manifest.json"
    )
    assert "--practical-acceptance" in command
    assert "--allow-source-environment-hash" in command


def test_frozen_panel_cpu_default_can_explicitly_leave_gpu_visible(monkeypatch) -> None:
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    cpu_env = _child_environment(Path("/repo"), use_gpu=False)
    gpu_env = _child_environment(Path("/repo"), use_gpu=True)
    assert cpu_env["CUDA_VISIBLE_DEVICES"] == "-1"
    assert "CUDA_VISIBLE_DEVICES" not in gpu_env
    assert gpu_env["PYTHONPATH"] == "/repo"


def test_completed_output_rejects_stale_provenance(tmp_path: Path) -> None:
    expected = {
        "identity": {"restart_index": 3, "k_model": 3, "fold": 0},
        "source_sha256": "source",
        "checkpoint_bundle_sha256": "checkpoint",
        "strict_audit_sha256": "audit",
        "practical_acceptance_sha256": "acceptance",
        "protocol": {"inference_steps": 2400, "learning_rate": 0.01},
    }
    (tmp_path / "r04_resume_infer.json").write_text(json.dumps({
        "status": "FROZEN_HELDOUT_DIAGNOSTIC_NOT_K_SELECTION",
        "inference_authority": "PRACTICAL_POSTHOC_ACCEPTANCE",
        "strict_audit_status": "STOPPED_OFF_PLATFORM",
        "selected_k": None,
        "frozen_inference_provenance": expected,
    }), encoding="utf-8")
    (tmp_path / "cell.json").write_text(json.dumps({
        "inference_platform": {"split_count": 2},
    }), encoding="utf-8")
    assert _completed_output(tmp_path, expected)
    stale = dict(expected)
    stale["strict_audit_sha256"] = "changed"
    assert not _completed_output(tmp_path, stale)
