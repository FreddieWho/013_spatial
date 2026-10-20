import numpy as np

p = np.load("infra/tls_weaklabel_probe_20261006/hooknet_crop/tls_pred.npy")
n = np.load("infra/tls_weaklabel_probe_20261006/hooknet_crop/neg_pred.npy")
print("shape", p.shape, flush=True)
s = int(round(p.shape[0] ** 0.5))
print("sqrt", s, s * s == p.shape[0], flush=True)
for name, a in (("tls", p), ("neg", n)):
    m = a.reshape(s, s, 3)
    c = m[s // 2 - 100:s // 2 + 100, s // 2 - 100:s // 2 + 100]
    print(name, "full-mean", m.reshape(-1, 3).mean(0).round(4),
          "center-mean", c.reshape(-1, 3).mean(0).round(4),
          "centermax", c.reshape(-1, 3).max(0).round(4), flush=True)
