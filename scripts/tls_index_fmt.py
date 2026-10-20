import h5py
f=h5py.File("/home/huyudi/013_spatial/data/external_tatls2026/adata_tls_full.h5ad","r")
idx=f["obs"]["_index"][:]
def s(x): return x.decode() if isinstance(x,bytes) else str(x)
print("n:",len(idx),flush=True)
print("first8:",[s(x) for x in idx[:8]],flush=True)
print("last3:",[s(x) for x in idx[-3:]],flush=True)
