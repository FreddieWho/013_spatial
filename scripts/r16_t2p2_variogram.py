#!/usr/bin/env python3
"""T2P2: variogram-matched Gaussian surrogate on real tissue geometry.

One method only (recovery contract): stationary exponential+nugget covariance
fit to the empirical semivariogram, one Cholesky per field, surrogates L@z.
Old rotation/remap is not the null. Thresholds live in the contract written by
`freeze` and are not rewritten by `run`.

Usable only if ALL pre-declared gates pass. Otherwise association p is
PERMANENTLY_NOT_CALIBRATED. Do not relax thresholds after seeing results.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

from r16 import axes as A
from r16 import census as C

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/r16/t2p2_20260928"
CACHE = ROOT / "infra/r16/census_cache_hvg10k/sections"
DIAG_PATIENTS = [
    "HTAN_VANDERBILT_CRC::PATIENT::4aa571d63ab1aa20",
    "HTAN_VANDERBILT_CRC::PATIENT::4f81822a6f2ecd63",
    "HTAN_VANDERBILT_CRC::PATIENT::0bd3490036cbd4de",
    "HTAN_VANDERBILT_CRC::PATIENT::286667cab598bb4a",
    "HTAN_VANDERBILT_CRC::PATIENT::0def2e0b01043545",
    "HTAN_VANDERBILT_CRC::PATIENT::07dc676d92f9d153",
]
REAL_GENES = ("CXCL13", "COL1A1")
N_BINS = 8
N_PAIRS = 80000
IMPL_REPS = 100
MISSPEC_REPS = 80
POWER_REPS = 40
SURROGATES = 99
JITTER = 1e-5


def contract_body():
    return {
        "version": 1,
        "status": "FROZEN_BEFORE_RESULTS",
        "method": "stationary exponential+nugget; empirical semivariogram; Cholesky surrogates L@z",
        "geometry": "6 T2 diagnostic patients, first manifest section, effect-blind",
        "coords": "Visium array_row/col to hex pitch when source obs has both and barcodes align; else cache pixel coords, flagged, not silently treated as pitch",
        "statistic": "mean(score | B-contour signed hop<=0) - mean(score | hop>=4); contour q=0.7",
        "real_genes": list(REAL_GENES),
        "gates": {
            "impl_fpr_reps": IMPL_REPS,
            "misspec_reps": MISSPEC_REPS,
            "power_reps": POWER_REPS,
            "surrogates": SURROGATES,
            "fpr_window": [0.02, 0.10],
            "power_min": 0.40,
            "bump_sd": 1.0,
            "fit_rrmse_max": 0.35,
            "min_sections_fit": 5,
            "surrogate_variogram_rrmse_max": 0.35,
            "surrogate_abs_corr_max": 0.20,
        },
        "families": ["implementation_grf", "graph_smooth", "zero_inflated_smooth"],
        "usable_iff": "implementation FPR in window AND both misspecification families in window AND power>=0.40 AND at least one real gene fit-passes on >=5/6 sections AND surrogate match gates pass on that gene",
        "failure_action": "else PERMANENTLY_NOT_CALIBRATED for Lane A association p; do not relax thresholds; known-model calibration does not substitute",
        "not_a_discovery_run": "no GO expansion, no p-values reported as biological findings",
    }


def freeze():
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "run_contract.json"
    if path.exists():
        print("contract already frozen; not rewritten")
        return
    body = contract_body()
    body["frozen_at_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    path.write_text(json.dumps(body, indent=2))
    print("frozen", path)


def semivariogram(xy, z, seed):
    n = len(z)
    rng = np.random.default_rng(seed)
    i = rng.integers(0, n, N_PAIRS)
    j = rng.integers(0, n, N_PAIRS)
    ok = i != j
    i, j = i[ok], j[ok]
    h = np.linalg.norm(xy[i] - xy[j], axis=1)
    g = 0.5 * (z[i] - z[j]) ** 2
    edges = np.quantile(h, np.linspace(0, 1, N_BINS + 1))
    for k in range(1, len(edges)):
        if edges[k] <= edges[k - 1]:
            edges[k] = edges[k - 1] + 1e-6
    lags, gamma, counts = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (h >= a) & (h < b if b < edges[-1] else h <= b)
        if m.sum() < 30:
            continue
        lags.append(float(np.median(h[m])))
        gamma.append(float(np.mean(g[m])))
        counts.append(int(m.sum()))
    return np.asarray(lags), np.asarray(gamma), counts


def fit_exp(lags, gamma):
    if len(lags) < 4 or not np.isfinite(gamma).all() or gamma.max() <= 0:
        return None
    best = None
    nug_grid = np.linspace(0, max(float(gamma[0]), 1e-8), 5)
    sill0 = max(float(gamma.max() - gamma[0]), 1e-6)
    for nug in nug_grid:
        for sill in np.linspace(0.5 * sill0, 1.8 * sill0, 7):
            for ell in np.geomspace(max(float(lags[1]), 1e-3), max(float(lags[-1]), 1e-2), 8):
                pred = nug + sill * (1 - np.exp(-lags / ell))
                rrmse = float(np.sqrt(np.mean((pred - gamma) ** 2)) / (np.mean(gamma) + 1e-12))
                if best is None or rrmse < best[0]:
                    best = (rrmse, float(nug), float(sill), float(ell))
    return best


def cov_matrix(xy, nugget, sill, ell):
    # Chunked Euclidean distances; full matrix is required for Cholesky.
    n = len(xy)
    C = np.empty((n, n), dtype=np.float64)
    step = 512
    for a in range(0, n, step):
        b = min(n, a + step)
        d = np.linalg.norm(xy[a:b, None, :] - xy[None, :, :], axis=2)
        C[a:b] = sill * np.exp(-d / max(ell, 1e-6))
    diag = nugget + JITTER * (sill + nugget + 1.0)
    np.fill_diagonal(C, np.diag(C) + diag)
    return C


def factor(xy, nugget, sill, ell):
    C = cov_matrix(xy, nugget, sill, ell)
    try:
        return np.linalg.cholesky(C)
    except np.linalg.LinAlgError:
        np.fill_diagonal(C, np.diag(C) + 1e-3 * (sill + 1.0))
        return np.linalg.cholesky(C)


def contrast(values, near, far):
    if near.sum() < 20 or far.sum() < 20:
        return None
    return float(values[near].mean() - values[far].mean())


def two_sided_p(t, nulls):
    nulls = np.asarray(nulls, dtype=float)
    nulls = nulls[np.isfinite(nulls)]
    if not len(nulls) or not np.isfinite(t):
        return None
    return float((1 + np.sum(np.abs(nulls) >= abs(t))) / (1 + len(nulls)))


def fpr(ps):
    ps = np.asarray([p for p in ps if p is not None], dtype=float)
    if not len(ps):
        return None, 0
    return float(np.mean(ps <= 0.05)), int(len(ps))


def in_window(rate, window):
    return rate is not None and window[0] <= rate <= window[1]


def hex_or_pixel(meta, barcodes, pixel):
    path = meta.get("source_matrix_locator")
    if not path or not Path(path).exists():
        return pixel, "pixel_fallback_no_source"
    try:
        with h5py.File(path, "r") as h:
            obs = h["obs"]
            if "array_row" not in obs or "array_col" not in obs or "_index" not in obs:
                return pixel, "pixel_fallback_no_array"
            idx = [b.decode() if isinstance(b, bytes) else str(b) for b in obs["_index"][:]]
            rows = np.asarray(obs["array_row"][:])
            cols = np.asarray(obs["array_col"][:])
    except Exception:
        return pixel, "pixel_fallback_read_error"
    lut = {b: (int(rows[i]), int(cols[i])) for i, b in enumerate(idx)}
    xy = np.empty_like(pixel)
    for i, b in enumerate(barcodes):
        key = b.decode() if isinstance(b, bytes) else str(b)
        if key not in lut:
            return pixel, "pixel_fallback_barcode_mismatch"
        r, c = lut[key]
        xy[i, 0] = c + (0.5 if r % 2 else 0.0)
        xy[i, 1] = r * (3 ** 0.5) / 2
    return xy, "hex_pitch"


def knn_idx(xy, k=6):
    tree = cKDTree(xy)
    _, idx = tree.query(xy, k=k + 1)
    return idx[:, 1:]


def smooth(eps, idx, steps=3, w=0.7):
    z = eps.astype(float).copy()
    for _ in range(steps):
        z = (1 - w) * z + w * z[idx].mean(axis=1)
    return z


def load_diag():
    man = json.loads((CACHE.parent / "cache_manifest.json").read_text())
    first = {}
    for row in man["rows"]:
        first.setdefault(row["patient_id"], row["stem"])
    proxy = json.loads((ROOT / "infra/r04/marker_proxy_combined.json").read_text())
    bgenes = list(proxy["classes"]["B"]["voted_genes"])
    sections = []
    for pid in DIAG_PATIENTS:
        stem = first[pid]
        mat, bcs, genes, pixel, meta = C.load_section(CACHE, stem)
        xy, coord_kind = hex_or_pixel(meta, bcs, pixel)
        gidx = {g: i for i, g in enumerate(genes)}
        xz = C.log1p_zscore(mat)
        sB, nB = A.axis_scores(xz, gidx, bgenes)
        w = C.spatial_weights(xy, C.KNN_SPATIAL)
        signed = A.signed_hop_distance(w, A.contour_mask(sB, 0.7))
        near = signed <= 0
        far = signed >= 4
        real = {}
        for g in REAL_GENES:
            if g in gidx:
                real[g] = xz[:, gidx[g]].astype(float)
        sections.append(dict(
            stem=stem, patient=pid, xy=xy, coord_kind=coord_kind,
            n=len(bcs), nB=nB, near=near, far=far, real=real,
            idx=knn_idx(xy),
        ))
        print(f"loaded {stem[:24]} n={len(bcs)} coord={coord_kind} Bgenes={nB} "
              f"near={int(near.sum())} far={int(far.sum())}", flush=True)
    return sections


def p_from_field(field, sec, rng, refit):
    t = contrast(field, sec["near"], sec["far"])
    if t is None:
        return None
    if refit:
        lags, gamma, _ = semivariogram(sec["xy"], field, int(rng.integers(1, 10**9)))
        fit = fit_exp(lags, gamma)
        if fit is None:
            return None
        _, nug, sill, ell = fit
        L = factor(sec["xy"], nug, sill, ell)
    else:
        L = sec["impl_L"]
    nulls = []
    for _ in range(SURROGATES):
        z = L @ rng.normal(size=sec["n"])
        c = contrast(z, sec["near"], sec["far"])
        if c is not None:
            nulls.append(c)
    return two_sided_p(t, nulls)


def run():
    contract = json.loads((OUT / "run_contract.json").read_text())
    gates = contract["gates"]
    assert gates["fpr_window"] == [0.02, 0.10]
    assert gates["fit_rrmse_max"] == 0.35
    t0 = time.time()
    sections = load_diag()
    rng = np.random.default_rng(20260928)

    fit_rows = []
    for sec in sections:
        for g, vals in sec["real"].items():
            lags, gamma, counts = semivariogram(sec["xy"], vals, 20260928)
            fit = fit_exp(lags, gamma)
            if fit is None:
                fit_rows.append(dict(stem=sec["stem"], gene=g, status="FIT_FAILED"))
                continue
            rrmse, nug, sill, ell = fit
            # surrogate match: 20 draws from fitted model
            L = factor(sec["xy"], nug, sill, ell)
            rrmses, corrs = [], []
            for d in range(20):
                z = L @ rng.normal(size=sec["n"])
                lg, gg, _ = semivariogram(sec["xy"], z, 1000 + d)
                # compare on shared lag count via refit prediction at observed lags
                pred = nug + sill * (1 - np.exp(-lags / ell))
                # surrogate variogram interpolated by its own fit
                sf = fit_exp(lg, gg)
                if sf is None:
                    continue
                spred = sf[1] + sf[2] * (1 - np.exp(-lags / sf[3]))
                rrmses.append(float(np.sqrt(np.mean((spred - pred) ** 2)) / (np.mean(pred) + 1e-12)))
                corrs.append(float(np.corrcoef(z, vals)[0, 1]))
            fit_rows.append(dict(
                stem=sec["stem"], gene=g, status="FIT", rrmse=rrmse,
                nugget=nug, sill=sill, length=ell, n_bins=len(lags),
                pair_counts=counts, coord_kind=sec["coord_kind"],
                surrogate_variogram_rrmse_med=float(np.median(rrmses)) if rrmses else None,
                surrogate_abs_corr_med=float(np.median(np.abs(corrs))) if corrs else None,
            ))
            print(f"  fit {g} {sec['stem'][:16]} rrmse={rrmse:.3f}", flush=True)
        # implementation covariance: median NN * 4, unit sill
        nn = np.median(np.linalg.norm(sec["xy"][sec["idx"][:, 0]] - sec["xy"], axis=1))
        sec["impl_params"] = (0.05, 1.0, float(max(nn * 4, 1e-3)))
        sec["impl_L"] = factor(sec["xy"], *sec["impl_params"])

    def collect(kind, maker, refit, reps):
        ps = []
        for sec in sections:
            for r in range(reps):
                field = maker(sec, rng)
                ps.append(p_from_field(field, sec, rng, refit=refit))
            print(f"  {kind} {sec['stem'][:16]} done reps={reps}", flush=True)
        rate, n = fpr(ps)
        return {"family": kind, "fpr": rate, "n": n, "p_values": ps}

    print("implementation GRF", flush=True)
    impl = collect(
        "implementation_grf",
        lambda sec, rng: sec["impl_L"] @ rng.normal(size=sec["n"]),
        refit=False,
        reps=IMPL_REPS,
    )
    print("graph smooth misspec", flush=True)
    graph = collect(
        "graph_smooth",
        lambda sec, rng: smooth(rng.normal(size=sec["n"]), sec["idx"]),
        refit=True,
        reps=MISSPEC_REPS,
    )
    print("zero-inflated smooth", flush=True)

    def zi(sec, rng):
        eps = rng.normal(size=sec["n"])
        eps[rng.random(sec["n"]) < 0.4] = 0.0
        return smooth(eps, sec["idx"])

    zinf = collect("zero_inflated_smooth", zi, refit=True, reps=MISSPEC_REPS)
    print("power", flush=True)

    def bumped(sec, rng):
        base = smooth(rng.normal(size=sec["n"]), sec["idx"])
        sd = float(np.std(base) + 1e-8)
        out = base.copy()
        out[sec["near"]] += gates["bump_sd"] * sd
        return out

    power = collect("power_bump", bumped, refit=True, reps=POWER_REPS)
    detect, n_pow = fpr(power["p_values"])

    # gene fit gate
    gene_pass = {}
    for g in REAL_GENES:
        rows = [r for r in fit_rows if r.get("gene") == g and r.get("status") == "FIT"]
        n_ok = sum(1 for r in rows if r["rrmse"] <= gates["fit_rrmse_max"]
                   and r.get("surrogate_variogram_rrmse_med") is not None
                   and r["surrogate_variogram_rrmse_med"] <= gates["surrogate_variogram_rrmse_max"]
                   and r["surrogate_abs_corr_med"] <= gates["surrogate_abs_corr_max"])
        gene_pass[g] = n_ok
    fit_ok = any(v >= gates["min_sections_fit"] for v in gene_pass.values())
    window = gates["fpr_window"]
    impl_ok = in_window(impl["fpr"], window)
    graph_ok = in_window(graph["fpr"], window)
    zi_ok = in_window(zinf["fpr"], window)
    pow_ok = detect is not None and detect >= gates["power_min"]
    usable = bool(impl_ok and graph_ok and zi_ok and pow_ok and fit_ok)
    if not impl_ok:
        verdict = "BLOCKED_IMPLEMENTATION"
    elif usable:
        verdict = "USABLE"
    else:
        verdict = "PERMANENTLY_NOT_CALIBRATED"

    def slim(block):
        return {k: block[k] for k in ("family", "fpr", "n")}

    summary = {
        "verdict": verdict,
        "usable": usable,
        "seconds": round(time.time() - t0, 1),
        "implementation": slim(impl),
        "graph_smooth": slim(graph),
        "zero_inflated_smooth": slim(zinf),
        "power_detection_at_0.05": detect,
        "power_n": n_pow,
        "gene_sections_passing_fit": gene_pass,
        "gates_passed": {
            "implementation_fpr": impl_ok,
            "graph_smooth_fpr": graph_ok,
            "zero_inflated_fpr": zi_ok,
            "power": pow_ok,
            "real_gene_fit": fit_ok,
        },
        "coord_kinds": sorted({s["coord_kind"] for s in sections}),
    }
    (OUT / "fit_rows.json").write_text(json.dumps(fit_rows, indent=2))
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    # reliability-style histogram, not a biological p
    with open(OUT / "calibration_fpr.tsv", "w") as f:
        f.write("family\tfpr\tn\n")
        for block in (impl, graph, zinf):
            f.write(f"{block['family']}\t{block['fpr']}\t{block['n']}\n")
        f.write(f"power_bump\t{detect}\t{n_pow}\n")
    print(json.dumps(summary, indent=2))
    return 0 if verdict != "BLOCKED_IMPLEMENTATION" else 2


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "freeze":
        freeze()
    elif len(sys.argv) > 1 and sys.argv[1] == "run":
        raise SystemExit(run())
    else:
        raise SystemExit("usage: freeze | run")
