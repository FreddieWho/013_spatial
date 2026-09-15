"""R-16 Tier-1 v2 ruler factory core (design v2.1, D-118).

Data-driven axis discovery on neighbor-augmented representation:
per-section per-lambda PCA (top-10 PCs = candidate axes) + LISA tile
coherence (value-permutation null + BH-FDR, valid for autocorrelation
tests per I-021 refined) + cross-patient loading matching.

Conventions (all disclosed, frozen):
- Neighborhood mean EXCLUDES self (own term is separate, BANKSY-style).
- Mixing: Z = sqrt(1-lam) * X_z + sqrt(lam) * M_z, lam in {0, 0.2, 0.5, 0.8}.
- LISA null permutes values (exchangeable under no-autocorrelation H0).
- Tile criterion: >=1 hotspot (connected significant component) of >=20 spots.
- Deterministic: all stochastic steps take explicit seeds.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree  # noqa: F401 (re-exported for tests)
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

from r16 import census as C

SEED = C.SEED
LAMBDAS = (0.0, 0.2, 0.5, 0.8)
K_PRIMARY = 15
K_RADIUS_CHECK = 30
N_PC = 10
MIN_HOTSPOT_SIZE = 20
NULL_DRAWS = C.NULL_DRAWS
# LISA tests one hypothesis PER SPOT, so the permutation floor 1/(draws+1)
# must sit below the BH line: 200 draws (floor 0.005) provably recovers
# nothing on a 256-spot test case; 999 draws (floor 0.001) recovers it fully.
LISA_DRAWS = 999
FDR_ALPHA = 0.05


# ---------------------------------------------------------------- neighborhoods


def spatial_mean_weights(coords: np.ndarray, k: int) -> sparse.csr_matrix:
    """Row-normalized spatial kNN mean weights, self EXCLUDED."""
    n = len(coords)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, idx = nn.kneighbors(coords)
    rows = np.repeat(np.arange(n), k)
    cols = idx[:, 1:].ravel()
    w = sparse.csr_matrix((np.full(len(rows), 1.0 / k), (rows, cols)), shape=(n, n))
    return w


def augment(x_z: np.ndarray, m_z: np.ndarray, lam: float) -> np.ndarray:
    """BANKSY-style mixing of own field and neighborhood field."""
    return float(np.sqrt(1.0 - lam)) * x_z + float(np.sqrt(lam)) * m_z


def zscore_columns(x: np.ndarray) -> np.ndarray:
    mu = x.mean(axis=0)
    sd = x.std(axis=0)
    sd[sd < 1e-8] = 1.0
    return (x - mu) / sd


# ---------------------------------------------------------------- candidate axes


def candidate_pcs(z_aug: np.ndarray, n_pc: int = N_PC,
                  seed: int = SEED) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """PCA on the augmented matrix. Returns (scores n×k, loadings p×k, var_ratio k)."""
    k = int(min(n_pc, min(z_aug.shape) - 1))
    model = PCA(n_components=k, random_state=seed)
    scores = model.fit_transform(z_aug)
    return scores, model.components_.T.copy(), model.explained_variance_ratio_.copy()


# ---------------------------------------------------------------- LISA tile coherence


def bh_fdr(p: np.ndarray, alpha: float = FDR_ALPHA) -> np.ndarray:
    """Benjamini-Hochberg reject mask."""
    p = np.asarray(p, dtype=float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order]
    thresh = (np.arange(1, m + 1) / m) * alpha
    below = ranked <= thresh
    if not below.any():
        return np.zeros(m, dtype=bool)
    kmax = np.flatnonzero(below).max()
    out = np.zeros(m, dtype=bool)
    out[order[:kmax + 1]] = True
    return out


def lisa_tile_stats(s: np.ndarray, w: sparse.csr_matrix, draws: int,
                    rng: np.random.Generator,
                    min_hotspot: int = MIN_HOTSPOT_SIZE,
                    alpha: float = FDR_ALPHA) -> dict:
    """Local Moran (LISA) on one field + hotspot summary.

    Hotspot = connected component (>= min_hotspot spots) of FDR-significant
    spots on the SYMMETRIZED spatial graph. Purely descriptive localization;
    promotion evidence comes from cross-patient loading recurrence, never
    from these p-values alone.
    """
    s = np.asarray(s, dtype=float)
    n = len(s)
    sd = s.std()
    z = (s - s.mean()) / (sd if sd > 0 else 1.0)
    lag = np.asarray((w @ z)).ravel()
    obs = z * lag
    null = np.empty((draws, n))
    for d in range(draws):
        zp = rng.permutation(z)
        null[d] = zp * np.asarray((w @ zp)).ravel()
    p_two = (np.sum(np.abs(null) >= np.abs(obs), axis=0) + 1) / (draws + 1)
    sig = bh_fdr(p_two, alpha)  # draws must satisfy 1/(draws+1) << alpha/n
    n_sig = int(sig.sum())
    sym = w.maximum(w.T)
    n_hot, max_hot = 0, 0
    if n_sig:
        sub = sym[sig][:, sig]
        n_comp, lab = connected_components(sub, directed=False)
        sizes = np.bincount(lab)
        big = sizes[sizes >= min_hotspot]
        n_hot, max_hot = int(len(big)), int(big.max()) if len(big) else 0
    return {
        "n_sig": n_sig,
        "frac_sig": float(n_sig / n),
        "n_hotspots": n_hot,
        "max_hotspot": max_hot,
        "tile_ok": bool(n_hot >= 1),
    }


# ---------------------------------------------------------------- matching helpers


def abs_cosine(a: np.ndarray, b: np.ndarray) -> float:
    """Sign-aligned cosine (handles PCA sign ambiguity)."""
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    return float(abs(a @ b / (na * nb)))


def top_loading_genes(loading: np.ndarray, genes: list[str], topn: int = 10) -> list[str]:
    order = np.argsort(-np.abs(loading))[:topn]
    return [genes[i] for i in order]
