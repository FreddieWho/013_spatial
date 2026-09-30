#!/usr/bin/env python3
"""Composition-residual TLS contrast on the frozen GSE175540 strict set.

Question, fixed before scores: after removing marker-composition and depth,
does the adaptive-immune TLS advantage remain, and do raw and depth agree?

Not a new discovery screen. Not localization. Not cross-cancer.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/gobp_composition_residual_20260928"
PRIOR = ROOT / "infra/gobp_external_score_20260927/run_contract.json"
AXES = ("B", "T", "Mye", "Epi", "Stromal", "Plasma")
PRIMARY = "GOBP_ADAPTIVE_IMMUNE_RESPONSE"
MIN_AXIS_GENES = 2
MIN_AXES = 3


def load_ext():
    path = ROOT / "scripts/external_score_gse175540.py"
    spec = importlib.util.spec_from_file_location("extscore", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def axis_gene_map():
    from r16 import axes as A
    proxy = json.loads((ROOT / "infra/r04/marker_proxy_combined.json").read_text())
    out = {k: list(proxy["classes"][k]["voted_genes"]) for k in ("B", "T", "Mye", "Epi", "Stromal")}
    out["Plasma"] = list(A.PLASMA_GENES)
    return out


def freeze():
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "run_contract.json"
    if path.exists():
        print("contract already frozen")
        return
    prior = json.loads(PRIOR.read_text())
    assert PRIMARY in prior["strict_programs"]
    axes = axis_gene_map()
    body = {
        "version": 1,
        "status": "FROZEN_BEFORE_RESIDUAL_SCORES",
        "frozen_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "question": "Does GOBP_ADAPTIVE_IMMUNE_RESPONSE TLS-vs-NO_TLS patient SMD remain after within-patient linear removal of Q+C6, in both raw and depth modes?",
        "scope": "Same 16 strict patients and 31 strict programs as D-148. Blanks unknown. No block-level or cross-cancer claim. Aggregation arms not repeated.",
        "patients": prior["strict_patients"],
        "programs": prior["strict_programs"],
        "primary_program": PRIMARY,
        "axes": AXES,
        "axis_rule": "Drop ILC. For each program, remove that program's input+readout symbols from each axis. Axis with <2 genes left is dropped, not zero-filled. Fewer than 3 remaining axes => NOT_SEPARABLE.",
        "q_rule": "Intercept, log10(lib_excl+1), z-scored detected-gene count. lib_excl and detected-gene count exclude that program's input+readout symbols (duplicate feature rows summed).",
        "score_rule": "raw = mean log1p(count). depth = mean log1p(count*1e4/lib_excl). Axis scores use the same mode and the same lib_excl as the program score.",
        "fit_rule": "OLS of program score on Q+C within each patient, fit only on explicit TLS and NO_TLS spots matched to matrix and in-tissue coordinates. Residual = score - fitted. Then same SMD as D-148 on residuals.",
        "remain_rule": {
            "per_mode": "REMAINS if residual median SMD>0.2, retained fraction vs same-mode unadjusted median SMD>=0.5, and >=14/16 patients have residual SMD>0. COLLAPSED if retained fraction<0.5 or positive patients<12. Else ATTENUATED.",
            "overall": "SURVIVES_BOTH_MODES only if both modes are REMAINS. Otherwise do not claim a composition-orthogonal TLS advantage.",
            "family": "31-program residual sign tests are descriptive context with Holm within mode. They do not change the primary remain rule.",
        },
        "input_hashes": {
            "prior_contract": sha256(PRIOR),
            "go_representatives.json": sha256(ROOT / "infra/repair_20260921/go_representatives.json"),
            "marker_proxy_combined.json": sha256(ROOT / "infra/r04/marker_proxy_combined.json"),
        },
        "axis_sizes": {k: len(v) for k, v in axes.items()},
    }
    path.write_text(json.dumps(body, indent=2))
    print("frozen", path)


def gene_vector(csc, rows):
    if not rows:
        return np.zeros(csc.shape[1])
    return np.asarray(csc[rows, :].sum(axis=0)).ravel()


def design_and_scores(ext, csc, symbol_rows, defs, program, axes, mode):
    inputs = defs[program]["input_genes"]
    ban = set(inputs) | set(defs[program]["readout_genes"])
    dup = [s for s in inputs if len(symbol_rows.get(s, [])) > 1]
    if dup:
        return None, f"strict_duplicate:{dup[0]}"
    in_rows = []
    for s in inputs:
        rows = symbol_rows.get(s, [])
        if len(rows) != 1:
            return None, f"missing_input:{s}"
        in_rows.append(rows[0])
    axis_cols = []
    axis_used = {}
    for name in AXES:
        genes = [g for g in axes[name] if g not in ban and len(symbol_rows.get(g, [])) == 1]
        axis_used[name] = len(genes)
        if len(genes) < MIN_AXIS_GENES:
            continue
        axis_cols.append([symbol_rows[g][0] for g in genes])
    if len(axis_cols) < MIN_AXES:
        return None, "NOT_SEPARABLE"
    ban_rows = []
    for s in ban:
        ban_rows.extend(symbol_rows.get(s, []))
    total = np.asarray(csc.sum(axis=0)).ravel()
    lib_excl = np.maximum(total - gene_vector(csc, ban_rows), 0.0)
    detected = np.asarray((csc > 0).sum(axis=0)).ravel().astype(float)
    if ban_rows:
        detected = detected - np.asarray((csc[ban_rows, :] > 0).sum(axis=0)).ravel()
    detected = np.maximum(detected, 0.0)

    def mean_score(rows):
        mat = csc[rows, :].toarray()
        if mode == "depth":
            mat = mat * (1e4 / np.maximum(lib_excl, 1.0))[None, :]
        return np.log1p(mat).mean(axis=0)

    score = mean_score(in_rows)
    C = np.column_stack([mean_score(rows) for rows in axis_cols])
    ng = (detected - detected.mean()) / (detected.std() + 1e-9)
    Q = np.column_stack([np.ones(csc.shape[1]), np.log10(lib_excl + 1.0), ng])
    X = np.column_stack([Q, C])
    return dict(score=score, X=X, n_axes=len(axis_cols), axis_used=axis_used), None


def classify(resid_smd, unadj_smd, n_pos, n):
    if unadj_smd is None or abs(unadj_smd) < 1e-8:
        retained = None
    else:
        retained = resid_smd / unadj_smd
    if resid_smd > 0.2 and retained is not None and retained >= 0.5 and n_pos >= 14 and n == 16:
        return "REMAINS", retained
    if retained is not None and (retained < 0.5 or n_pos < 12):
        return "COLLAPSED", retained
    return "ATTENUATED", retained


def run():
    contract = json.loads((OUT / "run_contract.json").read_text())
    assert contract["status"] == "FROZEN_BEFORE_RESIDUAL_SCORES"
    ext = load_ext()
    defs = json.loads((ROOT / "infra/repair_20260921/go_representatives.json").read_text())
    axes = axis_gene_map()
    sources = [s for s in ext.qual_sources() if s["patient_id"] in set(contract["patients"])]
    assert len(sources) == 16
    qc = {r["patient_id"]: r for r in csv.DictReader(
        open(ROOT / "infra/gobp_external_qualification_20260925/sample_qc.tsv"), delimiter="\t")}
    rows = []
    for src in sources:
        pid = src["patient_id"]
        csc, barcodes, col, symbol_rows, keep = ext.load_sample(src)
        tls = [b for b, v in keep.items() if v == "TLS"]
        no = [b for b, v in keep.items() if v == "NO_TLS"]
        assert len(tls) == int(qc[pid]["matched_tls_spots"])
        assert len(no) == int(qc[pid]["matched_no_tls_spots"])
        ix = np.array([col[b] for b in tls + no])
        tls_ix = np.array([col[b] for b in tls])
        no_ix = np.array([col[b] for b in no])
        for program in contract["programs"]:
            for mode in ("raw", "depth"):
                packed, err = design_and_scores(ext, csc, symbol_rows, defs, program, axes, mode)
                if packed is None:
                    rows.append(dict(program=program, patient=pid, mode=mode, status=err))
                    continue
                s = packed["score"]
                X = packed["X"]
                beta, _, _, _ = np.linalg.lstsq(X[ix], s[ix], rcond=None)
                fitted = X[ix] @ beta
                ss_tot = float(((s[ix] - s[ix].mean()) ** 2).sum())
                r2 = 1 - float(((s[ix] - fitted) ** 2).sum()) / ss_tot if ss_tot > 0 else None
                resid = s - X @ beta
                unadj = ext.endpoints(s[tls_ix], s[no_ix])
                adj = ext.endpoints(resid[tls_ix], resid[no_ix])
                rows.append(dict(
                    program=program, patient=pid, mode=mode, status="OK",
                    n_axes=packed["n_axes"],
                    unadj_smd=None if unadj is None else unadj["smd"],
                    resid_smd=None if adj is None else adj["smd"],
                    resid_auc=None if adj is None else adj["auc"],
                    qc_r2=r2,
                ))
        print("done", pid, flush=True)
    with open(OUT / "patient_residuals.tsv", "w", newline="") as f:
        cols = ["program", "patient", "mode", "status", "n_axes", "unadj_smd", "resid_smd", "resid_auc", "qc_r2"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})
    summary = []
    by = {}
    for r in rows:
        by.setdefault((r["program"], r["mode"]), []).append(r)
    for (program, mode), rs in sorted(by.items()):
        ok = [r for r in rs if r["status"] == "OK" and r["resid_smd"] is not None and r["unadj_smd"] is not None]
        if len(ok) < 16:
            summary.append(dict(program=program, mode=mode, status="INCOMPLETE", n=len(ok)))
            continue
        resid = np.array([r["resid_smd"] for r in ok])
        unadj = np.array([r["unadj_smd"] for r in ok])
        n_pos = int((resid > 0).sum())
        label, retained = classify(float(np.median(resid)), float(np.median(unadj)), n_pos, len(ok))
        p = float(binomtest(n_pos, len(ok), 0.5).pvalue)
        summary.append(dict(
            program=program, mode=mode, status=label, n=len(ok), n_pos=n_pos,
            median_unadj_smd=float(np.median(unadj)),
            median_resid_smd=float(np.median(resid)),
            retained_fraction=retained, sign_p=p,
            median_qc_r2=float(np.median([r["qc_r2"] for r in ok])),
        ))
    # Holm within mode, descriptive
    for mode in ("raw", "depth"):
        items = [s for s in summary if s["mode"] == mode and "sign_p" in s]
        if not items:
            continue
        adj = ext.holm(np.array([s["sign_p"] for s in items]))
        for s, a in zip(items, adj):
            s["sign_p_holm"] = float(a)
    primary = {s["mode"]: s for s in summary if s["program"] == PRIMARY}
    both = (
        primary.get("raw", {}).get("status") == "REMAINS"
        and primary.get("depth", {}).get("status") == "REMAINS"
    )
    verdict = "SURVIVES_BOTH_MODES" if both else "DOES_NOT_SURVIVE"
    with open(OUT / "program_summary.tsv", "w", newline="") as f:
        cols = ["program", "mode", "status", "n", "n_pos", "median_unadj_smd",
                "median_resid_smd", "retained_fraction", "sign_p", "sign_p_holm", "median_qc_r2"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for s in summary:
            w.writerow({c: s.get(c, "") for c in cols})
    receipt = {
        "finished_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "verdict": verdict,
        "primary": primary,
        "n_rows": len(rows),
    }
    (OUT / "receipt.json").write_text(json.dumps(receipt, indent=2))
    print(json.dumps({"verdict": verdict, "primary": primary}, indent=2))


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2 or sys.argv[1] not in ("freeze", "run"):
        raise SystemExit("usage: freeze | run")
    freeze() if sys.argv[1] == "freeze" else run()
