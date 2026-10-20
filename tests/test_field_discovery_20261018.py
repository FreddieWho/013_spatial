"""Focused tests for target isolation, sample grouping, and source count layouts."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/field_discovery_20261018.py"
spec = importlib.util.spec_from_file_location("field_discovery", SCRIPT)
fd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fd)


class DiscoveryTests(unittest.TestCase):
    def test_target_perturbation_does_not_change_inputs_or_qc(self):
        rng = np.random.default_rng(7)
        x = sparse.csr_matrix(rng.poisson(3, (30, 12)))
        genes = np.array([f"G{i}" for i in range(12)])
        panel = {"T": list(genes[:3]), "A": list(genes[2:8]), "B": list(genes[7:])}
        a = fd.visible_features(x, genes, panel, "T", ["A", "B"])
        changed = x.toarray(); changed[:, :3] += 10000
        b = fd.visible_features(sparse.csr_matrix(changed), genes, panel, "T", ["A", "B"])
        np.testing.assert_array_equal(a, b)
        self.assertFalse(np.array_equal(fd.target_readout(x, genes, panel["T"]),
                                       fd.target_readout(sparse.csr_matrix(changed), genes, panel["T"])))

    def test_neighborhoods_exclude_self_and_do_not_cross_sections(self):
        xy = np.array([[0, 0], [1, 0], [2, 0], [3, 0]], float)
        near, outer, pitch = fd.neighbor_graphs(xy)
        self.assertEqual(pitch, 1)
        np.testing.assert_array_equal(near.diagonal(), 0)
        np.testing.assert_array_equal(outer.diagonal(), 0)
        self.assertEqual(near.multiply(outer).nnz, 0)
        np.testing.assert_allclose(np.asarray(near.sum(axis=1)).ravel(), 1)
        changed = np.array([[1], [2], [3], [100]], dtype=float)
        self.assertEqual((near @ changed)[0, 0], 2)

    def test_unit_cannot_cross_folds(self):
        rows = [{"section_id": "a", "unit_id": "p", "fold": 0},
                {"section_id": "b", "unit_id": "p", "fold": 1}]
        with self.assertRaisesRegex(ValueError, "crosses"):
            fd.assert_grouped(rows)

    def test_square_csc_is_not_silently_treated_as_csr(self):
        arr = np.array([[0, 1, 2], [3, 0, 4], [5, 6, 0]], dtype=np.float32)
        with tempfile.TemporaryDirectory() as tmp:
            for fmt in ["csc", "csr"]:
                path = Path(tmp) / f"{fmt}.h5"
                mat = getattr(sparse, f"{fmt}_matrix")(arr)
                with h5py.File(path, "w") as f:
                    g = f.create_group("X"); g.attrs["encoding-type"] = f"{fmt}_matrix"
                    g.attrs["shape"] = mat.shape
                    for key in ["data", "indices", "indptr"]:
                        g.create_dataset(key, data=getattr(mat, key))
                with h5py.File(path) as f:
                    np.testing.assert_array_equal(fd.read_matrix(f["X"], 3, 3).toarray(), arr)

    def test_allowlist_rejects_annotation_and_path_escape(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            guard = fd.InputGuard(root, {root / "allowed.json"}, root / "out")
            with self.assertRaisesRegex(PermissionError, "allowlist"):
                guard.check(root / "annotation/incoming/manual.json")
            with self.assertRaises(PermissionError):
                guard.check(root / "../outside.json")


if __name__ == "__main__":
    unittest.main()
