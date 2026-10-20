import h5py
import numpy as np

p = "/home/huyudi/013_spatial/data/external_tatls2026/adata_tls_full.h5ad"
f = h5py.File(p, "r")
print("TOP:", list(f.keys()), flush=True)
X = f["X"]
if isinstance(X, h5py.Group):
    print("X sparse:", X.attrs.get("shape"), X.attrs.get("h5py_version"), list(X.keys()), flush=True)
else:
    print("X dense:", X.shape, X.dtype, flush=True)
o = f["obs"]
keys = list(o.keys())
print("n_obs_keys:", len(keys), "n_obs:", len(o["__index__"]) if "__index__" in o else "?", flush=True)
for k in keys:
    if k.startswith("__"):
        continue
    d = o[k]
    if isinstance(d, h5py.Group):
        print(" obs:", k, "<categorical>", list(d.keys()), flush=True)
    else:
        print(" obs:", k, d.shape, str(d.dtype), flush=True)
print("obsm:", list(f["obsm"].keys()), flush=True)
print("uns:", list(f["uns"].keys()), flush=True)


def show_cats(group, name, maxn=12):
    if name not in group:
        return
    d = group[name]
    if "categories" in d:
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in d["categories"][:]]
        codes = d["codes"][:]
        vals, counts = np.unique(codes, return_counts=True)
        print(f" {name}: ncat={len(cats)}", flush=True)
        for v, c in zip(vals[:maxn], counts[:maxn]):
            label = cats[v] if 0 <= v < len(cats) else f"NA({v})"
            print(f"   {label}: {c}", flush=True)
        if len(vals) > maxn:
            print(f"   ... +{len(vals)-maxn} more", flush=True)
    else:
        try:
            v = d[:]
            print(f" {name}: dtype={d.dtype}", flush=True)
        except Exception as e:
            print(f" {name}: unreadable {e}", flush=True)


for k in keys:
    if not k.startswith("__"):
        show_cats(o, k)
