#!/usr/bin/env python3
"""Four-cancer ST scoring (D-216): frozen 1742 panel on HEST Visium
bowel/lung/breast/liver. Outputs are candidate weak-label inputs only.
No labels are produced; all spots stay unknown.
"""
import csv
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.tls_peripheral import module  # noqa: E402
from scripts.tls_multiscale_gobp import build_endpoints  # noqa: E402
from scripts.gobp_halo_screen import score_sets  # noqa: E402

HEST = ROOT / "data/other_sources/hest1k"
OUT = ROOT / "infra/tls_fourcancer_20261006"
ORGANS = ("Bowel", "Lung", "Breast", "Liver")


def section_list():
    rows = list(csv.DictReader(open(HEST.parent / "hest1k" / "selected_samples.tsv"), delimiter="\t"))
    out = [(r["id"], r["organ"]) for r in rows
           if r["organ"] in ORGANS and r["st_technology"] == "Visium"]
    return out


def load_hest(sid):
    with h5py.File(HEST / "st" / f"{sid}.h5ad", "r") as f:
        X = f["X"]
        n_obs = len(f["obs"]["_index"])
        if isinstance(X, h5py.Group):
            data, indices, indptr = X["data"][:], X["indices"][:], X["indptr"][:]
            n_var = len(f["var"]["_index"])
            if len(indptr) - 1 == n_obs:
                X = sparse.csr_matrix((data, indices, indptr), shape=(n_obs, n_var))
            elif len(indptr) - 1 == n_var:
                X = sparse.csc_matrix((data, indices, indptr), shape=(n_obs, n_var)).tocsr()
            else:
                raise ValueError(f"indptr {len(indptr)} matches neither n_obs {n_obs} nor n_var {n_var}")
            assert X.shape == (n_obs, n_var), X.shape
        else:
            X = sparse.csr_matrix(X[:])
            if X.shape[0] != n_obs:
                X = X.T.tocsr()
        genes = np.array([g.decode() if isinstance(g, bytes) else str(g)
                          for g in f["var"]["_index"][:]])
        obs = f["obs"]
        bc = np.array([b.decode() if isinstance(b, bytes) else str(b) for b in obs["_index"][:]])
        ar = np.asarray(obs["array_row"][:]).astype(int)
        ac = np.asarray(obs["array_col"][:]).astype(int)
        it = np.asarray(obs["in_tissue"][:]).astype(int) == 1
    order = np.flatnonzero(it)
    coords = np.stack([ar[order], ac[order]], axis=1)
    xy = module("tls_pool_expand").hex_xy(coords[:, 0], coords[:, 1], 100)
    return X[order], genes, xy, bc[order]


def main():
    (OUT / "scores").mkdir(parents=True, exist_ok=True)
    (OUT / "geometry").mkdir(parents=True, exist_ok=True)
    panel = json.loads((ROOT / "infra/tls_multiscale_gobp_20261002/panel_gobp.json").read_text())
    control = json.loads((ROOT / "infra/tls_multiscale_gobp_20261002/panel_controls.json").read_text())
    allsets = panel | control
    assert len(allsets) == 1742
    secs = section_list()
    print(f"sections: {len(secs)}", flush=True)
    manifest = []
    prev = {}
    mf = OUT / "manifest.json"
    if mf.exists():
        try:
            for r in json.loads(mf.read_text()):
                prev[r.get("section_id")] = r
        except Exception:
            pass
    t0 = time.time()
    for sid, organ in secs:
        out = OUT / "scores" / f"hest-{sid}.npz"
        if out.exists() and (OUT / "geometry" / f"hest-{sid}.npz").exists():
            print("skip", sid, flush=True)
            if f"hest-{sid}" in prev:
                manifest.append(prev[f"hest-{sid}"])
            continue
        t1 = time.time()
        try:
            X, genes, xy, bc = load_hest(sid)
        except Exception as e:
            print(f"FAIL load {sid}: {e}", flush=True)
            manifest.append({"section_id": f"hest-{sid}", "organ": organ, "status": "LOAD_FAIL"})
            continue
        index = {g: i for i, g in enumerate(genes)}
        members = [np.array([index[g] for g in gl if g in index], dtype=int)
                   for gl in allsets.values()]
        S, k = score_sets(X, members)
        fraction = k / np.array([len(v) for v in allsets.values()])
        valid = (fraction >= .8) & (k >= 10) & (np.nanstd(S, axis=0) > 1e-12)
        lib = np.asarray(X.sum(axis=1)).ravel().astype(np.float32)
        np.savez_compressed(out, scores=S.astype(np.float32), set_ids=np.array(list(allsets)),
                            coverage=fraction, n_genes=k, valid=valid, libsize=lib,
                            detected=np.diff(X.indptr))
        labels = np.array([""] * len(xy))
        graph, ends, members_c, edge = build_endpoints(xy, labels)
        sparse.save_npz(OUT / "geometry" / f"hest-{sid}_graph.npz", graph)
        np.savez_compressed(OUT / "geometry" / f"hest-{sid}.npz", xy=xy, labels=labels,
                            order=np.arange(len(xy)),
                            endpoint_members=np.concatenate(members_c) if members_c else np.array([], int),
                            endpoint_indptr=np.r_[0, np.cumsum([len(v) for v in members_c])])
        manifest.append({"section_id": f"hest-{sid}", "organ": organ, "status": "OK",
                         "n_spots": int(len(xy)), "n_valid_sets": int(valid.sum()),
                         "seconds": round(time.time() - t1, 1)})
        print(f"done {sid} {S.shape} valid={int(valid.sum())} {time.time()-t1:.0f}s", flush=True)
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    ok = sum(1 for m in manifest if m.get("status") == "OK")
    (OUT / "STATUS").write_text(f"SCORED {ok}/{len(manifest)}\n")
    print(f"total {time.time()-t0:.0f}s ok={ok}/{len(manifest)}", flush=True)


if __name__ == "__main__":
    main()
