import h5py
import numpy as np
for sid in ["TENX13", "TENX39", "NCBI776"]:
    p = f"/home/huyudi/013_spatial/data/other_sources/hest1k/st/{sid}.h5ad"
    f = h5py.File(p, "r")
    X = f["X"]
    if isinstance(X, h5py.Group):
        import scipy.sparse as sp
        m = sp.csr_matrix((X["data"][:], X["indices"][:], X["indptr"][:]),
                          shape=tuple(X.attrs["shape"]))
    else:
        m = X[:]
    lib = np.asarray(m.sum(axis=1)).ravel()
    det = np.asarray((m > 0).sum(axis=1)).ravel()
    print(f"{sid}: shape={m.shape} lib_med={np.median(lib):.0f} "
          f"det_med={np.median(det):.0f} frac>500={(lib > 500).mean():.2f} "
          f"obs={list(f['obs'].keys())[:8]}", flush=True)
