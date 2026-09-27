"""GSE175540 patient-level external TLS-label association for the 45 GO BP candidates.

Arms (frozen before any expression score; see run_contract.json):
  Arm 1 STRICT: 16 strict two-class patients x 31 strictly-unambiguous programs.
  Arm 2 expanded, same 16 patients, all 45 programs, two aggregation routes:
    AGG-SUM : duplicate symbol -> sum raw counts across its GRCh38
              Gene-Expression feature rows, then one gene in the score.
    AGG-DROP: duplicate symbols excluded from the input gene list
              (program scored on remaining genes; NOT_TESTABLE if <2 remain).

Scoring modes (GO-stage convention, D-140): raw = mean log1p(count) over
effective input genes; depth_normalized = mean log1p(count*1e4/lib_excl) with
lib_excl = total UMI minus all frozen input+readout symbols (duplicates summed),
identical in both aggregation arms so only the numerator differs.

Endpoint (frozen): per patient, per program, per arm, per mode, on explicit
TLS vs explicit NO_TLS spots matched to matrix + in-tissue coordinates
(blanks stay unknown):
  primary     SMD = (mean_TLS - mean_NO_TLS)/pooled spot SD
  sensitivity AUC of the input score separating TLS from NO_TLS.
Cross-patient: median SMD, sign counts, two-sided sign p (descriptive) + Holm
within each arm+mode SMD family. Patient is the unit; no block-level claims;
renal overlap with USZ noted so no cross-cancer claim.

Usage: freeze (write contract) | run (compute or verify) .
Env: LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse
from scipy.stats import binomtest, rankdata

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "infra/gobp_external_score_20260927"
QUAL = ROOT / "infra/gobp_external_qualification_20260925"
DEFS = ROOT / "infra/repair_20260921/go_representatives.json"
CANDCTX = ROOT / "infra/gobp_stage_20260922/candidate_context.tsv"
MIN_CLASS_N = 3
MIN_GENES = 2


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def decode(values):
    return [v.decode() if isinstance(v, bytes) else v for v in values]


def load_labels(path):
    with gzip.open(path, "rt", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    return {r["Barcode"].strip(): (r["TLS_2_cat"] or "").strip() for r in rows}


def load_positions(path):
    result = {}
    with gzip.open(path, "rt", newline="", encoding="utf-8-sig") as handle:
        for row in csv.reader(handle):
            result[row[0].strip()] = row[1]
    return result


def cand45():
    with open(CANDCTX, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return sorted(r["program"] for r in rows if r["molecular_robust_gt20"] == "True")


def qual_sources():
    with open(QUAL / "run_contract.json") as f:
        return json.load(f)["selected_sources"]


def strict_patients():
    with open(QUAL / "sample_qc.tsv", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return sorted(r["patient_id"] for r in rows
                  if r["strict_two_class_patient_eligible"] == "true"), rows


def strict_programs(cands):
    with open(QUAL / "program_qualification_summary.tsv", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    elig = {r["program_id"] for r in rows
            if r["strictly_eligible_across_all_18_sources"] == "true"}
    return sorted(set(cands) & elig)


def load_sample(src):
    with h5py.File(ROOT / src["matrix_path"], "r") as handle:
        matrix = handle["matrix"]
        barcodes = decode(matrix["barcodes"][:])
        fg = matrix["features"]
        names = decode(fg["name"][:])
        ftypes = decode(fg["feature_type"][:])
        genomes = decode(fg["genome"][:])
        shape = [int(v) for v in matrix["shape"][:]]
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
    assert shape == [len(names), len(barcodes)]
    assert len(indptr) == len(barcodes) + 1
    csc = sparse.csc_matrix((data.astype(float), indices, indptr), shape=shape)
    symbol_rows = defaultdict(list)
    for i, (nm, ft, gn) in enumerate(zip(names, ftypes, genomes)):
        if ft == "Gene Expression" and gn == "GRCh38":
            symbol_rows[nm].append(i)
    labels = load_labels(ROOT / src["source_path"])
    positions = load_positions(ROOT / src["positions_path"])
    keep = {}
    for bc, val in labels.items():
        if val in ("TLS", "NO_TLS") and bc in set(barcodes) \
                and positions.get(bc) == "1":
            keep[bc] = val
    col = {bc: i for i, bc in enumerate(barcodes)}
    return csc, barcodes, col, symbol_rows, keep


def gene_counts(csc, rows):
    if not rows:
        return np.zeros(csc.shape[1])
    return np.asarray(csc[rows, :].sum(axis=0)).ravel()


def holm(pvals):
    order = np.argsort(pvals)
    m = len(pvals)
    adj = np.empty(m)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (m - rank) * pvals[idx])
        adj[idx] = min(running, 1.0)
    return adj


def endpoints(s_tls, s_no):
    n1, n2 = len(s_tls), len(s_no)
    if n1 < MIN_CLASS_N or n2 < MIN_CLASS_N:
        return None
    m1, m2 = float(np.mean(s_tls)), float(np.mean(s_no))
    v1, v2 = float(np.var(s_tls, ddof=1)), float(np.var(s_no, ddof=1))
    pooled = ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)
    if pooled <= 0:
        return None
    smd = (m1 - m2) / np.sqrt(pooled)
    ranks = rankdata(np.concatenate([s_tls, s_no]))
    auc = (ranks[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n2)
    return dict(n_tls=n1, n_no=n2, mean_tls=m1, mean_no=m2,
                smd=smd, auc=float(auc))


def freeze():
    OUT.mkdir(parents=True, exist_ok=True)
    cands = cand45()
    assert len(cands) == 45
    defs = json.load(open(DEFS))
    assert all(p in defs for p in cands)
    sources = qual_sources()
    assert len(sources) == 18
    strict_pats, _ = strict_patients()
    assert len(strict_pats) == 16
    strict_progs = strict_programs(cands)
    assert len(strict_progs) == 31, strict_progs
    contract = {
        "version": 1,
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "FROZEN_BEFORE_SCORING",
        "scope": "Patient-level TLS-label association of 45-candidate input scores "
                 "on GSE175540. No model fitting on GSE175540. Blanks unknown. "
                 "No block-level inference. Renal overlap with USZ: independent "
                 "patients only, not cross-cancer validation.",
        "arms": {
            "strict": "16 strict patients x 31 strictly-unambiguous programs.",
            "agg_sum": "Same 16 patients x all 45 programs; duplicate symbol "
                       "counts summed across GRCh38 Gene-Expression rows.",
            "agg_drop": "Same 16 patients x all 45 programs; duplicate symbols "
                        "excluded from input list; patient-program NOT_TESTABLE "
                        "if <2 effective input genes remain.",
        },
        "scoring_modes": {
            "raw": "mean log1p(count) over effective input genes.",
            "depth_normalized": "mean log1p(count*1e4/lib_excl); lib_excl = total "
                                "UMI minus all frozen input+readout symbols "
                                "(duplicates summed), identical in both agg arms.",
        },
        "endpoint": {
            "primary": "per-patient SMD of input score, TLS minus NO_TLS.",
            "sensitivity": "per-patient AUC.",
            "cross_patient": "median SMD, sign counts, two-sided sign p "
                             "descriptive + Holm within each arm+mode SMD family.",
            "min_class_n": MIN_CLASS_N, "min_genes": MIN_GENES,
        },
        "strict_patients": strict_pats,
        "strict_programs": strict_progs,
        "candidates": cands,
        "input_hashes": {
            "go_representatives.json": sha256(DEFS),
            "candidate_context.tsv": sha256(CANDCTX),
            "qualification_run_contract.json": sha256(QUAL / "run_contract.json"),
            "sample_qc.tsv": sha256(QUAL / "sample_qc.tsv"),
            "program_qualification_summary.tsv":
                sha256(QUAL / "program_qualification_summary.tsv"),
            "scoring_script": sha256(Path(__file__)),
        },
    }
    with open(OUT / "run_contract.json", "w") as f:
        json.dump(contract, f, indent=2, ensure_ascii=False)
    print(f"contract frozen: {len(cands)} candidates, "
          f"{len(strict_pats)} patients, {len(strict_progs)} strict programs")


def run():
    contract = json.load(open(OUT / "run_contract.json"))
    for name, digest in contract["input_hashes"].items():
        if name == "scoring_script":
            continue
        print(f"input pinned: {name} {digest[:12]}")
    defs = json.load(open(DEFS))
    sources = [s for s in qual_sources()
               if s["patient_id"] in set(contract["strict_patients"])]
    assert len(sources) == 16
    qc = {r["patient_id"]: r for r in
          list(csv.DictReader(open(QUAL / "sample_qc.tsv"), delimiter="\t"))}
    pat_rows, sum_rows = [], []
    for src in sources:
        pid = src["patient_id"]
        csc, barcodes, col, symbol_rows, keep = load_sample(src)
        tls_bc = sorted(b for b, v in keep.items() if v == "TLS")
        no_bc = sorted(b for b, v in keep.items() if v == "NO_TLS")
        assert len(tls_bc) == int(qc[pid]["matched_tls_spots"]), pid
        assert len(no_bc) == int(qc[pid]["matched_no_tls_spots"]), pid
        tls_ix = np.array([col[b] for b in tls_bc])
        no_ix = np.array([col[b] for b in no_bc])
        total = np.asarray(csc.sum(axis=0)).ravel()
        frozen_syms = set()
        for p in contract["candidates"]:
            frozen_syms.update(defs[p]["input_genes"] + defs[p]["readout_genes"])
        lib_excl = total.copy()
        for sym in frozen_syms:
            rows = symbol_rows.get(sym, [])
            if rows:
                lib_excl -= gene_counts(csc, rows)
        lib_excl = np.maximum(lib_excl, 0.0)
        for prog in contract["candidates"]:
            inputs = defs[prog]["input_genes"]
            dup = sorted({s for s in inputs if len(symbol_rows.get(s, [])) > 1})
            routes = {}
            if prog in set(contract["strict_programs"]):
                assert not dup, (prog, dup)
                rows = [symbol_rows[s][0] for s in inputs]
                routes["strict"] = (rows, len(inputs))
            sum_vecs, drop_rows = [], []
            for s in inputs:
                rws = symbol_rows.get(s, [])
                if len(rws) > 1:
                    sum_vecs.append(gene_counts(csc, rws))
                elif len(rws) == 1:
                    sum_vecs.append(gene_counts(csc, rws))
                    drop_rows.append(rws[0])
            routes["agg_sum"] = (None, len(inputs))
            routes["agg_drop"] = (drop_rows, len(drop_rows))
            for arm, spec in routes.items():
                if arm == "agg_sum":
                    ng = spec[1]
                    if ng < MIN_GENES or not sum_vecs:
                        continue
                    mat = np.vstack(sum_vecs)
                else:
                    rows, ng = spec[0], spec[1]
                    if ng < MIN_GENES or not rows:
                        continue
                    mat = csc[rows, :].toarray()
                for mode in ("raw", "depth_normalized"):
                    if mode == "raw":
                        sm = np.log1p(mat)
                    else:
                        scale = 1e4 / np.maximum(lib_excl, 1.0)
                        sm = np.log1p(mat * scale[None, :])
                    score = sm.mean(axis=0)
                    ep = endpoints(score[tls_ix], score[no_ix])
                    if ep is None:
                        continue
                    pat_rows.append(dict(program=prog, patient=pid, arm=arm,
                                         mode=mode, n_input_genes=ng,
                                         n_dup_symbols=len(dup), **ep))
    with open(OUT / "patient_endpoints.tsv", "w", newline="") as f:
        cols = ["program", "patient", "arm", "mode", "n_input_genes",
                "n_dup_symbols", "n_tls", "n_no", "mean_tls", "mean_no",
                "smd", "auc"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in pat_rows:
            w.writerow({c: r[c] for c in cols})
    by = defaultdict(list)
    for r in pat_rows:
        by[(r["program"], r["arm"], r["mode"])].append(r)
    fams = defaultdict(list)
    for key, rs in by.items():
        smds = np.array([r["smd"] for r in rs])
        n_pos = int((smds > 0).sum())
        n = len(smds)
        p = binomtest(n_pos, n, 0.5).pvalue
        fams[(rs[0]["arm"], rs[0]["mode"])].append((key, p))
        sum_rows.append(dict(program=key[0], arm=key[1], mode=key[2],
                             n_patients=n, n_pos=n_pos,
                             median_smd=float(np.median(smds)),
                             median_auc=float(np.median([r["auc"] for r in rs])),
                             sign_p=p))
    holm_map = {}
    for fam, items in fams.items():
        adj = holm(np.array([p for _, p in items]))
        for (key, _), a in zip(items, adj):
            holm_map[key] = float(a)
    with open(OUT / "arm_summary.tsv", "w", newline="") as f:
        cols = ["program", "arm", "mode", "n_patients", "n_pos",
                "median_smd", "median_auc", "sign_p", "sign_p_holm"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in sorted(sum_rows, key=lambda r: (r["arm"], r["mode"], r["program"])):
            r["sign_p_holm"] = holm_map[(r["program"], r["arm"], r["mode"])]
            w.writerow(r)
    outputs = ["patient_endpoints.tsv", "arm_summary.tsv", "run_contract.json"]
    receipt = {"finished_at_utc": datetime.now(timezone.utc).isoformat(),
               "n_patient_rows": len(pat_rows), "n_summary_rows": len(sum_rows),
               "sha256": {n: sha256(OUT / n) for n in outputs}}
    with open(OUT / "receipt.json", "w") as f:
        json.dump(receipt, f, indent=2)
    print(f"done: {len(pat_rows)} patient rows, {len(sum_rows)} summary rows")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=("freeze", "run"))
    args = ap.parse_args()
    if args.action == "freeze":
        freeze()
    else:
        if (OUT / "receipt.json").exists():
            receipt = json.load(open(OUT / "receipt.json"))
            for name, digest in receipt["sha256"].items():
                assert sha256(OUT / name) == digest, name
            print("outputs verified against receipt; no recompute")
        else:
            run()


if __name__ == "__main__":
    main()
