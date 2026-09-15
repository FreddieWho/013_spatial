"""Tests for the GraphST bridge (D-121). Locally runnable: anndata only,
no torch/scanpy/GraphST needed (run-script torch path is stubbed)."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent


def test_export_one_section_roundtrip(tmp_path):
    from scripts.r16_export_graphst_bridge import main as _  # noqa
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "expmod", ROOT / "scripts/r16_export_graphst_bridge.py")
    # module executes main() only under __main__; import just checks syntax
    assert spec is not None
    # run export on first real section via CLI (repo root on PYTHONPATH, as usual)
    import os
    env = dict(os.environ, PYTHONPATH=str(ROOT))
    r = subprocess.run(
        [sys.executable, str(ROOT / "scripts/r16_export_graphst_bridge.py"),
         "--out", str(tmp_path / "gbridge"), "--sections",
         "0000_HTAN_8270_AS_2_filtered_trimmed.h5ad"],
        capture_output=True, text=True, cwd=ROOT, timeout=600, env=env)
    assert r.returncode == 0, r.stderr[-2000:]
    import anndata as ad
    a = ad.read_h5ad(tmp_path / "gbridge" /
                      "0000_HTAN_8270_AS_2_filtered_trimmed.h5ad.h5ad")
    assert a.n_obs == 665 and a.n_vars == 10000
    assert np.all(np.asarray(a.X.todense()) == np.floor(np.asarray(a.X.todense())))
    assert "spatial" in a.obsm and a.obsm["spatial"].shape == (665, 2)
    assert list(a.obs.columns) == ["section_id", "patient_id", "barcode"]
    assert len(set(a.var_names)) == 10000
    man = json.loads((tmp_path / "gbridge" / "manifest.json").read_text())
    assert len(man) == 1 and man[0]["n_spots"] == 665


def test_run_script_with_stubbed_graphst(tmp_path):
    """Full wrapper path with a fake GraphST module (no torch/scanpy needed)."""
    import types
    fake_mod = types.ModuleType("GraphST")

    class FakeGS:
        def __init__(self, adata, **kw):
            self.adata = adata.copy()
            self.kw = kw

        def train(self):
            rng = np.random.default_rng(0)
            n = self.adata.n_obs
            # two-blobs embedding keyed to sections half/half -> 2 clusters
            emb = np.zeros((n, 8))
            emb[:n // 2, 0] = 5.0
            emb[n // 2:, 1] = 5.0
            emb += rng.normal(scale=0.1, size=emb.shape)
            self.adata.obsm["emb"] = emb
            return self.adata

    # self-contained: export one section first (seconds)
    import os
    env = dict(os.environ, PYTHONPATH=str(ROOT))
    r = subprocess.run(
        [sys.executable, str(ROOT / "scripts/r16_export_graphst_bridge.py"),
         "--out", str(tmp_path / "gbridge"), "--sections",
         "0000_HTAN_8270_AS_2_filtered_trimmed.h5ad"],
        capture_output=True, text=True, cwd=ROOT, timeout=600, env=env)
    assert r.returncode == 0, r.stderr[-2000:]
    fake_mod.GraphST = FakeGS
    sys.modules["GraphST"] = fake_mod
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "grun", ROOT / "scripts/r16_graphst_run.py")
        grun = importlib.util.module_from_spec(spec)
        sys.argv = ["grun", "--data", str(tmp_path / "gbridge"),
                    "--out", str(tmp_path / "gout"), "--epochs", "2",
                    "--device", "cpu", "--dry-run"]
        # torch import inside main would fail locally? torch CPU exists; ok either way
        spec.loader.exec_module(grun)
        grun.main()
    finally:
        sys.modules.pop("GraphST", None)
    import pandas as pd
    df = pd.read_csv(tmp_path / "gout" / "labels_armG.csv")
    assert list(df.columns) == ["section", "orig_barcode",
                                "GRAPHST_025", "GRAPHST_05", "GRAPHST_10"]
    assert len(df) == 665
    params = json.loads((tmp_path / "gout" / "params.json").read_text())
    assert params["epochs"] == 3 and params["seed"] == 20260914


def test_leiden_on_embedding_helper():
    sys.path.insert(0, str(ROOT / "scripts"))
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "grun2", ROOT / "scripts/r16_graphst_run.py")
    src = (ROOT / "scripts/r16_graphst_run.py").read_text()
    assert "def leiden_on_embedding" in src and "def resolve_device" in src
    # exercise the real helper without torch: load module with stubbed torch/GraphST
    import types
    sys.modules.setdefault("torch", types.ModuleType("torch"))
    grun = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(grun)
    except Exception:
        pytest.skip("run-script import needs torch/scanpy stack")
        return
    rng = np.random.default_rng(1)
    emb = np.vstack([rng.normal(0, 0.3, size=(60, 8)),
                     rng.normal(5, 0.3, size=(60, 8))]).astype(float)
    out = grun.leiden_on_embedding(emb)
    assert set(out) == {"GRAPHST_025", "GRAPHST_05", "GRAPHST_10"}
    assert len(set(out["GRAPHST_10"].tolist()) - {-1}) >= 2
