#!/usr/bin/env python3
"""v9: composition control for the label-anchor field test.

Adds per-spot NNLS composition (6 classes, frozen GSE236581 reference) to the
conditioning set, then repeats the v8.2 section-level contrast and null.
Contract: infra/tls_label_field_20261001/run_contract_v9.json.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy import sparse, optimize
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/tls_label_field_20261001"
COMP = OUT / "composition"
CONTRACT = OUT / "run_contract_v9.json"
REF_DIR = ROOT / "infra/r04/composition_ref"
V7_TARGETS = ROOT / "infra/gobp_halo_20261001/target_sets_v7.tsv"
HALLMARK = ROOT / "data/geneset/h.all.v2026.1.Hs.symbols.gmt"
SINGLE_GENES = ["CXCL13", "CCL19", "CCL21", "LTB", "CD74", "HLA-DRA", "IFNG",
                "CXCL9", "GZMB", "CD8A", "MKI67", "COL1A1"]
BIN_EDGES = np.arange(0, 1001, 100)
NBINS = len(BIN_EDGES) - 1
NEAR_I, FAR_I = 2, 7
MIN_BIN_SPOTS = 30
MIN_SECTIONS = 8
RING_FACTOR = 1.3
CLASSES = ("T", "B", "Mye", "ILC", "Epi", "Stromal")


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_gmt(path):
    out = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3:
            out[p[0]] = [g for g in p[2:] if g]
    return out


def load_reference():
    z = np.load(REF_DIR / "reference_counts.npz", allow_pickle=False)
    M = sparse.csr_matrix((z["counts_data"], z["counts_indices"], z["counts_indptr"]),
                          shape=tuple(z["counts_shape"]))
    symbols = [str(s) for s in z["panel_symbols"]]
    major = np.array([str(s) for s in z["major"]])
    row_sums = np.asarray(M.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0
    N = sparse.diags(1.0 / row_sums) @ M
    profiles = []
    for cls in CLASSES:
        rows = np.flatnonzero(major == cls)
        prof = np.asarray(N[:, rows].mean(axis=1)).ravel()
        profiles.append(prof / prof.sum())
    return np.column_stack(profiles), symbols, dict(zip(CLASSES, range(len(CLASSES))))


def spot_proportions(profiles, Msub):
    dense = Msub.toarray() if sparse.issparse(Msub) else np.asarray(Msub, float)
    dense = np.asarray(dense, dtype=float)
    rs = dense.sum(axis=1, keepdims=True)
    rs[rs == 0] = 1.0
    targets = dense / rs
    out = np.zeros((targets.shape[0], profiles.shape[1]))
    for i in range(targets.shape[0]):
        w, _ = optimize.nnls(profiles, targets[i])
        t = w.sum()
        out[i] = w / t if t > 0 else np.full_like(w, 1.0 / len(w))
    return out


def decile_design(values, n_bins=10):
    q = np.quantile(values, np.linspace(0, 1, n_bins + 1)[1:-1])
    d = np.digitize(values, q)
    D = np.zeros((len(values), n_bins))
    D[np.arange(len(values)), d] = 1.0
    return D


def trend_design(xy):
    u = (xy[:, 0] - xy[:, 0].mean()) / (np.ptp(xy[:, 0]) / 2 + 1e-9)
    v = (xy[:, 1] - xy[:, 1].mean()) / (np.ptp(xy[:, 1]) / 2 + 1e-9)
    return np.column_stack([np.ones(len(xy)), u, v, u * u, u * v, v * v])


def main():
    global OUT, COMP
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=200)
    ap.add_argument("--sections", default="")
    ap.add_argument("--tag", default="")
    ap.add_argument("--skip-deconv", action="store_true")
    ap.add_argument("--mode", choices=["marker_proxy", "nnls"], default="marker_proxy",
                    help="marker_proxy: tissue-agnostic module scores as composition covariates (v9.1); nnls: frozen CRC reference (invalid out of domain, kept for the record)")
    ap.add_argument("--de-genes-per-class", type=int, default=100,
                    help="marker-panel variant for the nnls mode (0 = full shared panel)")
    ap.add_argument("--profile-space", choices=["proportion", "logmean"], default="proportion")
    ap.add_argument("--output-dir", type=Path, required=True,
                    help="New audit directory; historical outputs must not be overwritten")
    args = ap.parse_args()
    OUT = args.output_dir
    OUT.mkdir(parents=True, exist_ok=False)
    COMP = OUT / "composition"
    if args.draws < 1:
        raise ValueError("draws must be positive")
    rng = np.random.default_rng(20261001)

    v9 = json.loads(CONTRACT.read_text())
    assert v9["status"] == "FROZEN_BEFORE_SCORING"
    screen = load_module("screen", ROOT / "scripts/gobp_halo_screen.py")
    nullmod = load_module("nullmod", ROOT / "scripts/gobp_halo_null.py")
    labmod = load_module("labmod", ROOT / "scripts/tls_label_field.py")
    pool = load_module("pool", ROOT / "scripts/tls_pool_expand.py")
    hallmark = load_module("hallmark", ROOT / "scripts/tls_field_hallmark.py")

    # composition covariates
    if args.mode == "marker_proxy":
        mp = json.loads((ROOT / "infra/r04/marker_proxy_combined.json").read_text())
        axmod = load_module("axmod", ROOT / "r16/axes.py")
        cov_defs = [("COV_T", mp["classes"]["T"]["voted_genes"]),
                    ("COV_Mye", mp["classes"]["Mye"]["voted_genes"]),
                    ("COV_Epi", mp["classes"]["Epi"]["voted_genes"]),
                    ("COV_Stromal", mp["classes"]["Stromal"]["voted_genes"]),
                    ("COV_Plasma", list(axmod.PLASMA_GENES))]
        profiles = None
        _classes = None
        panel_symbols = None
        print(f"composition mode=marker_proxy covariates={[c for c, _ in cov_defs]}", flush=True)
    else:
        deconv = load_module("deconv", ROOT / "scripts/r04_deconvolve_nnls.py")
        ref = deconv.load_reference(REF_DIR)
        profiles_full, _classes = deconv.mean_profiles(ref, space=args.profile_space)
        if args.de_genes_per_class > 0:
            keep = set(deconv.differential_markers(ref, args.de_genes_per_class))
            kept_idx = [i for i, g in enumerate(ref["genes"]) if g in keep]
            profiles = profiles_full[kept_idx]
            panel_symbols = [ref["genes"][i] for i in kept_idx]
        else:
            profiles = profiles_full
            panel_symbols = list(ref["genes"])
        cov_defs = None
        print(f"reference: profiles {profiles.shape} classes={list(_classes)} "
              f"space={args.profile_space} de_genes={args.de_genes_per_class}", flush=True)

    targets = list(csv.DictReader(open(V7_TARGETS), delimiter="\t"))
    readouts = [(t["set_id"], t["genes"].split(",")) for t in targets]
    hm = read_gmt(HALLMARK)
    for key in ("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                "HALLMARK_G2M_CHECKPOINT", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"):
        readouts.append((key, hm[key]))
    for g in SINGLE_GENES:
        readouts.append((f"GENE_{g.replace('-', '_')}", [g]))

    diag = list(csv.DictReader(open(ROOT / "infra/tls_pool_expand_20260930/pool_diagnostics.tsv"),
                               delimiter="\t"))
    ok = [r for r in diag if r["status"] == "ok" and r["cohort"] in ("USZ", "GSE175540")]
    if args.sections:
        keep = set(args.sections.split(","))
        ok = [r for r in ok if r["section_id"] in keep]

    COMP.mkdir(parents=True, exist_ok=True)
    rows, null_rows = [], []
    exclusions = []
    t_start = time.time()
    for row in ok:
        sid, cohort = row["section_id"], row["cohort"]
        t0 = time.time()
        loaded = hallmark.load_section_by_id(pool, cohort, sid)
        if loaded is None or loaded[0] is None:
            exclusions.append({"section_id": sid, "reason": "loader_failed"})
            continue
        _, X, genes, xy, pitch, labels, index = loaded
        y = np.asarray(labels.get("y"), dtype=object)
        tls_mask = (y == "TLS")
        if tls_mask.sum() < 20:
            exclusions.append({"section_id": sid, "reason": "fewer_than_20_TLS_spots"})
            continue
        gidx = dict(index)
        if args.mode == "nnls":
            comp_path = COMP / f"{sid}.npz"
            if args.skip_deconv and comp_path.exists():
                props = np.load(comp_path)["proportions"]
            else:
                sel_rows, sel_cols = [], []
                for i, sym in enumerate(panel_symbols):
                    v = gidx.get(sym)
                    if v is None:
                        continue
                    for r in (v if isinstance(v, list) else [v]):
                        sel_rows.append(r)
                        sel_cols.append(i)
                Ssel = sparse.csr_matrix((np.ones(len(sel_rows), dtype=np.float32),
                                          (sel_rows, sel_cols)),
                                         shape=(X.shape[1], len(panel_symbols)))
                Msub = (X @ Ssel).tocsr()
                props = deconv.deconvolve_nnls(Msub, profiles, space=args.profile_space)
                np.savez_compressed(comp_path, proportions=props,
                                    classes=np.array(_classes),
                                    n_genes_used=len(set(sel_cols)))
            cov_labels = list(_classes)
        else:
            props = None
            cov_labels = [c for c, _ in cov_defs]
        def member(gl):
            idxs = []
            for g in gl:
                v = gidx.get(g)
                if v is None:
                    continue
                idxs.extend(v if isinstance(v, list) else [v])
            return np.array(sorted(idxs), dtype=np.int32)

        members = [member(gl) for _, gl in readouts]
        if args.mode == "marker_proxy":
            members = members + [member(gl) for _, gl in cov_defs]
        members = members + [member(pool.signature_symbols())]
        S, k = screen.score_sets(X, members)
        s = S[:, -1].astype(np.float64)
        if args.mode == "marker_proxy":
            props = S[:, len(readouts):len(readouts) + len(cov_defs)]
        if not np.isfinite(s).all():
            continue
        B, cnt, dist, keepv = labmod.distance_bins(xy, tls_mask, pitch)
        if min(cnt[NEAR_I], cnt[FAR_I]) < MIN_BIN_SPOTS:
            exclusions.append({"section_id": sid, "reason": "insufficient_distance_bins", "n_near": int(cnt[NEAR_I]), "n_far": int(cnt[FAR_I])})
            print(f"{sid}: bins too small, skipped", flush=True)
            continue
        Bn = B @ sparse.diags(1.0 / np.maximum(cnt, 1))
        trend = trend_design(xy)
        dcomp = np.column_stack([decile_design(props[:, c]) for c in range(props.shape[1])])
        X_own = np.column_stack([decile_design(s), trend])
        X_comp = np.column_stack([decile_design(s), dcomp, trend])
        Q_own = screen.residual_operator(X_own)
        Q_comp = screen.residual_operator(X_comp)
        pairs, pd = nullmod.pair_list(xy, rng)
        cache = {}
        for j, (set_id, _) in enumerate(readouts):
            if len(members[j]) == 0:
                continue
            r = S[:, j].astype(np.float64)
            if not np.isfinite(r).all():
                continue
            R_own = screen.residualize(r, Q_own)
            R_comp = screen.residualize(r, Q_comp)

            def contrast(vec):
                m = np.asarray(Bn.T @ vec)
                return float(m[NEAR_I] - m[FAR_I])

            rows.append({"section_id": sid, "cohort": cohort, "set_id": set_id,
                         "own_contrast": round(contrast(R_own), 6),
                         "comp_contrast": round(contrast(R_comp), 6),
                         "raw_contrast": round(contrast(r), 6),
                         "n_near": int(cnt[NEAR_I]), "n_far": int(cnt[FAR_I]),
                         "k_used": int(k[j])})
            for arm, R in (("own", R_own), ("comp", R_comp)):
                obs = contrast(R)
                corr = nullmod.band_corr(R, pairs, pd) if pairs is not None else {}
                sigma = nullmod.sigma_from_bands(corr)
                note = "NOT_CALIBRATED"
                if np.isfinite(sigma):
                    key = (arm, round(sigma, 1))
                    if key not in cache:
                        cache[key] = nullmod.smoothed_noise(xy, sigma, args.draws, rng)
                    S_sim = nullmod.quantile_map(cache[key], R)
                else:
                    null_rows.append({"section_id": sid, "cohort": cohort, "set_id": set_id,
                                      "arm": arm, "obs": round(obs, 6), "null_median": "",
                                      "p_one_sided": "", "sigma": "",
                                      "note": "NOT_TESTABLE_SIGMA_UNFIT"})
                    continue
                S_sim = screen.residualize(S_sim, Q_own if arm == "own" else Q_comp)
                means = np.asarray(Bn.T @ S_sim)
                null_con = means[NEAR_I, :] - means[FAR_I, :]
                p = float((1 + np.sum(null_con >= obs)) / (args.draws + 1))
                null_rows.append({"section_id": sid, "cohort": cohort, "set_id": set_id,
                                  "arm": arm, "obs": round(obs, 6),
                                  "null_median": round(float(np.median(null_con)), 6),
                                  "p_one_sided": round(p, 4),
                                  "sigma": round(sigma, 1) if np.isfinite(sigma) else "",
                                  "note": note})
        print(f"{sid}: near={int(cnt[NEAR_I])} far={int(cnt[FAR_I])} "
              + (",".join(f"{c}={props[:, i].mean():.2f}" for i, c in enumerate(cov_labels))
                 if props is not None else "marker_proxy")
              + f" elapsed={time.time()-t0:.1f}s", flush=True)

    seen = {(r["section_id"], r["set_id"]) for r in rows}
    missing_readouts = [{"section_id": sid, "set_id": set_id,
                         "status": "NOT_TESTABLE_NO_FINITE_SCORE"}
                        for sid in sorted({r["section_id"] for r in rows})
                        for set_id, _ in readouts if (sid, set_id) not in seen]
    (OUT / "readout_exclusions.json").write_text(json.dumps(missing_readouts, indent=2) + "\n")
    (OUT / "exclusions.json").write_text(json.dumps(exclusions, indent=2) + "\n")
    (OUT / "audit_status.json").write_text(json.dumps({
        "status": "NOT_CALIBRATED", "decision": "D-167", "draws": args.draws,
        "n_sections_requested": len(ok), "excluded": len(exclusions),
        "changes": ["Visium doubled-column geometry", "rank-aware projection",
                    "sigma factor corrected", "project surrogate through same design",
                    "no white-noise fallback"],
        "interpretation": "Exploratory sensitivity only; no formal positive or negative claim."
    }, indent=2) + "\n")
    with (OUT / f"label_field_comp_sections{args.tag}.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    with (OUT / f"label_field_comp_null{args.tag}.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(null_rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(null_rows)
    out = []
    for set_id, _ in readouts:
        rec = {"set_id": set_id}
        for arm in ("own", "comp"):
            sub = [r for r in null_rows if r["set_id"] == set_id and r["arm"] == arm]
            sub = [r for r in sub if r["p_one_sided"] != ""]
            ps = np.array([float(r["p_one_sided"]) for r in sub])
            obs = np.array([float(r["obs"]) for r in sub])
            nl = np.array([float(r["null_median"]) for r in sub])
            frac = float((ps <= 0.05).mean()) if len(ps) else float("nan")
            med_obs = float(np.median(obs)) if len(obs) else float("nan")
            med_null = float(np.median(nl)) if len(nl) else float("nan")
            rec[f"{arm}_n"] = len(ps)
            rec[f"{arm}_frac_p05"] = round(frac, 3) if frac == frac else ""
            rec[f"{arm}_median"] = round(med_obs, 6) if med_obs == med_obs else ""
            rec[f"{arm}_claim"] = ""
            rec[f"{arm}_status"] = "NOT_CALIBRATED"
            rec[f"{arm}_screen_rule_pass"] = int(len(ps) >= MIN_SECTIONS and frac >= 0.5 and med_obs > med_null)
        out.append(rec)
    out.sort(key=lambda r: (-(r["comp_frac_p05"] if isinstance(r["comp_frac_p05"], float) else -1),
                            -(r["comp_median"] if isinstance(r["comp_median"], float) else -9)))
    with (OUT / f"label_field_comp_summary{args.tag}.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader(); w.writerows(out)
    print(f"wrote sections={len(rows)} null={len(null_rows)} elapsed={time.time()-t_start:.0f}s")
    for r in out[:10]:
        print(f"  {r['set_id'][:46]:46} own n={r['own_n']:2} f={r['own_frac_p05']} | "
              f"comp n={r['comp_n']:2} f={r['comp_frac_p05']} med={r['comp_median']} "
              f"claim={r['comp_claim']}")


if __name__ == "__main__":
    main()
