"""R-16 joint-embedding discovery track (D-119).

Pooled single-cell-style integration over the 47 HTAN Vanderbilt sections
ONLY (discovery role; validation lineages never touched):
  Arm A: pool log1p -> ComBat(batch=section) -> pooled z-score -> PCA30
         -> kNN15 -> Leiden x3 resolutions (min cluster 20).
  Arm B: same, but on BANKSY-style neighbor-augmented features
         [own, nb-mean] (spatial kNN, self excluded, row-normalized).
Internal judge: split-half patient replication (15/15, seeded) run
independently per arm; candidates match across halves by best cosine.
External judge (ST-CRC/USZ) stays a separate, later step.

Deterministic: all stochastic steps take explicit seeds.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

from r16 import census as C

SEED = C.SEED
COMBAT_KEY = "section"
N_PC = 30
AUGMENT_LAMBDA = 0.5
KNN_NB = 6


def pool_log1p(mats: list) -> np.ndarray:
    """Stack per-section log1p dense float32. Memory: ~4.5GB for 47 sections."""
    parts = []
    for m in mats:
        x = m.toarray().astype(np.float32)
        np.log1p(x, out=x)
        parts.append(x)
    return np.vstack(parts)


def combat_correct(x_log: np.ndarray, batches: np.ndarray, max_iter: int = 30):
    """ComBat batch correction, parametric EB (Johnson et al. 2007).

    In-house numpy (no scanpy dep: scanpy import is broken in this env via
    matplotlib/GLIBCXX mismatch; recorded in D-119). Fully deterministic.
    """
    x = np.asarray(x_log, dtype=np.float64)
    batches = np.asarray(batches)
    uniq, bidx = np.unique(batches, return_inverse=True)
    n, g = x.shape
    alpha = x.mean(axis=0)
    var_pool = x.var(axis=0)
    var_pool[var_pool < 1e-8] = 1e-8
    sig = np.sqrt(var_pool)
    z = (x - alpha) / sig
    nb = len(uniq)
    gamma_hat = np.vstack([z[bidx == i].mean(axis=0) for i in range(nb)])
    delta2_hat = np.vstack([z[bidx == i].var(axis=0) for i in range(nb)])
    delta2_hat = np.clip(delta2_hat, 1e-8, None)
    gamma_bar = gamma_hat.mean(axis=1)
    tau2_bar = np.clip(gamma_hat.var(axis=1), 1e-8, None)
    m = delta2_hat.mean(axis=1)
    v = np.clip(delta2_hat.var(axis=1), 1e-8, None)
    lam_bar = m * m / v + 2.0
    theta_bar = (m * m * m + m * v) / v
    gamma_star = gamma_hat.copy()
    delta2_star = delta2_hat.copy()
    counts = np.array([(bidx == i).sum() for i in range(nb)], dtype=float)
    for _ in range(max_iter):
        g_old = gamma_star.copy()
        for i in range(nb):
            zi = z[bidx == i]
            gamma_star[i] = (counts[i] * tau2_bar[i] * gamma_hat[i] + delta2_star[i] * gamma_bar[i]) / (counts[i] * tau2_bar[i] + delta2_star[i])
            sse = ((zi - gamma_star[i]) ** 2).sum(axis=0)
            delta2_star[i] = (theta_bar[i] + 0.5 * sse) / (counts[i] / 2.0 + lam_bar[i] - 1.0)
        if np.max(np.abs(gamma_star - g_old)) < 1e-4:
            break
    out = np.empty_like(z)
    for i in range(nb):
        out[bidx == i] = (z[bidx == i] - gamma_star[i]) / np.sqrt(delta2_star[i])
    return (out * sig + alpha).astype(np.float32)


def neighbor_mean_log1p(x_log: np.ndarray, coords: np.ndarray, k: int = KNN_NB):
    """Row-normalized spatial kNN mean of log1p, self EXCLUDED (factory convention)."""
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, idx = nn.kneighbors(coords)
    rows = np.repeat(np.arange(len(coords)), k)
    cols = idx[:, 1:].ravel()
    w = sparse.csr_matrix((np.full(len(rows), 1.0 / k), (rows, cols)),
                          shape=(len(coords), len(coords)))
    return np.asarray((w @ x_log))


def augment_blocks(x_z: np.ndarray, n_z: np.ndarray, lam: float = AUGMENT_LAMBDA):
    """BANKSY-style concat [sqrt(1-l)*own, sqrt(l)*neighbor]."""
    a = float(np.sqrt(1.0 - lam))
    b = float(np.sqrt(lam))
    return np.hstack([a * x_z, b * n_z])


def zscore_pooled(x: np.ndarray) -> np.ndarray:
    mu = x.mean(axis=0)
    sd = x.std(axis=0)
    sd[sd < 1e-8] = 1.0
    return (x - mu) / sd


def patient_halves(patients: list[str], seed: int = SEED):
    """Deterministic 15/15 split of the 30 patients."""
    uniq = sorted(set(patients))
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(uniq))
    half = len(uniq) // 2
    a = {uniq[i] for i in order[:half]}
    return a, set(uniq) - a


def best_match_table(sigs_a: np.ndarray, sigs_b: np.ndarray):
    """Best cosine of each row of A against all rows of B (zero rows -> 0)."""
    na = np.linalg.norm(sigs_a, axis=1, keepdims=True)
    nb = np.linalg.norm(sigs_b, axis=1, keepdims=True)
    na[na == 0] = 1.0
    nb[nb == 0] = 1.0
    sim = (sigs_a / na) @ (sigs_b / nb).T
    return sim.max(axis=1)
