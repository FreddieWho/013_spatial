#!/usr/bin/env python3
"""Nested patient-LOPO single-input-gene baseline for the frozen GO BP 45."""

import csv
import hashlib
import json
import shutil
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from r16.recovery.corrected import fit_from_moments
from r16.recovery.gobp_stage import FeatureEngine, moments
from r16.recovery.repair_pipeline import axes, sections


OUT = ROOT / "infra/gobp_single_marker_20260924"
STAGE = ROOT / "infra/gobp_stage_20260922"
SOURCE = ROOT / "infra/repair_20260921"
DEFINITION_PATH = SOURCE / "go_representatives.json"
MAX_FREE_BYTES = 1_200_000_000_000
N_BASE = 9
N_FULL = 10


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    rows = list(rows)
    if not rows:
        raise ValueError(f"refusing to write an empty table: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def atomic_json(path, value):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)


def json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def section_paths(row):
    return [SOURCE / "cache" / f"{row['key']}_{suffix}" for suffix in ("counts.npy", "source_genes.json")]


def make_manifest(candidate_ids, definitions):
    return [{
        "program_id": pid,
        "definition_hash": definitions[pid]["definition_hash"],
        "covered_size": definitions[pid]["covered_size"],
        "input_genes": sorted(definitions[pid]["input_genes"]),
        "readout_genes": definitions[pid]["readout_genes"],
    } for pid in candidate_ids]


def collect_inputs(sec_rows):
    fixed = [
        SOURCE / "sections.json", SOURCE / "genes.json", DEFINITION_PATH,
        ROOT / "infra/r04/marker_proxy_combined.json",
        STAGE / "run_contract.json", STAGE / "candidate_context.tsv",
        STAGE / "molecular_patients.tsv", STAGE / "stage_receipt.json",
        ROOT / "docs/plan.md", ROOT / "docs/decisions.md",
    ]
    paths = fixed[:]
    for row in sec_rows:
        paths.extend(section_paths(row))
    missing = [p for p in paths if not p.is_file()]
    if missing:
        raise FileNotFoundError(missing[0])
    return {p.relative_to(ROOT).as_posix(): digest(p) for p in sorted(set(paths))}


def frozen_contract(manifest, input_hashes):
    helper_paths = [
        Path(__file__).resolve(),
        ROOT / "r16/recovery/gobp_stage.py",
        ROOT / "r16/recovery/corrected.py",
        ROOT / "r16/recovery/repair_pipeline.py",
    ]
    return {
        "version": 1,
        "audit_date": "2026-09-24",
        "status": "FROZEN_BEFORE_COMPUTATION",
        "scope": "Conditional exploratory single-gene comparison for the 45 programs selected by the completed GO BP stage; raw and depth-normalized; ST-CRC and USZ patient-level LOPO.",
        "candidate_selection": "Reuse molecular_robust_gt20=true from candidate_context.tsv. This 45-program list was selected using all patients and is post-selected; the comparison is not independent confirmation.",
        "models": {
            "M0": "intercept + log1p background library + log1p detected genes + the six frozen clean composition axes",
            "M1": "M0 + the existing full-program mean input-gene score",
            "single_gene": "M0 + one gene selected only from that program's frozen input_genes; input and readout genes remain excluded from M0 QC/composition features",
        },
        "split_and_selection": "Outer leave-one-patient-out within each cohort, matching the completed GO stage. Within each outer training set, inner leave-one-patient-out chooses the input gene with the lowest mean held-out MSE; fit that gene on all outer-training patients and score the untouched outer patient. Equal section moments within patient, then equal patient weight. Ties use ascending gene symbol.",
        "scoring": "Reuse FeatureEngine transforms for raw and depth_normalized; compare patient-level test MSE. Improvement over selected single gene = 100 * (single_gene_MSE - M1_MSE) / single_gene_MSE; positive means the full program has lower MSE.",
        "scope_limits": [
            "No program-selection re-run inside outer folds; the fixed 45 are post-selected and results remain exploratory.",
            "No p-values or confirmatory family-wise claims across overlapping GO programs.",
            "CPU/local only; no external data, GPU, retraining, or edits to canonical prior-stage outputs.",
        ],
        "candidate_programs_sha256": json_hash(manifest),
        "input_sha256": input_hashes,
        "code_sha256": {p.relative_to(ROOT).as_posix(): digest(p) for p in helper_paths},
    }


def save_npz_atomic(path, **arrays):
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    temp.replace(path)


def write_section_mode(engine, row, mode, candidate_ids, definitions, input_digest):
    section_dir = OUT / "section_moments"
    section_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{row['key']}_{mode}"
    npz_path = section_dir / f"{tag}.npz"
    receipt_path = section_dir / f"{tag}.json"
    if npz_path.exists() or receipt_path.exists():
        if not (npz_path.is_file() and receipt_path.is_file()):
            raise RuntimeError(f"partial checkpoint pair exists: {tag}")
        prior = json.loads(receipt_path.read_text(encoding="utf-8"))
        if prior.get("input_sha256") != input_digest or prior.get("output_sha256") != digest(npz_path):
            raise RuntimeError(f"checkpoint hash mismatch: {tag}")
        return prior

    n_programs = len(candidate_ids)
    max_genes = max(len(definitions[pid]["input_genes"]) for pid in candidate_ids)
    base_gram = np.full((n_programs, N_BASE, N_BASE), np.nan)
    base_rhs = np.full((n_programs, N_BASE), np.nan)
    target_second = np.full(n_programs, np.nan)
    full_gram = np.full((n_programs, N_FULL, N_FULL), np.nan)
    full_rhs = np.full((n_programs, N_FULL), np.nan)
    marker_cross = np.full((n_programs, max_genes, N_BASE), np.nan)
    marker_second = np.full((n_programs, max_genes), np.nan)
    marker_rhs = np.full((n_programs, max_genes), np.nan)
    marker_count = np.zeros(n_programs, dtype=np.int16)
    statuses = np.full(n_programs, "NOT_COMPUTED", dtype="U40")
    index = {pid: i for i, pid in enumerate(candidate_ids)}
    started = time.monotonic()

    for valid, x, y, bad in engine.blocks(mode, batch_size=n_programs):
        for pid, reason in bad:
            statuses[index[pid]] = reason
        if not valid:
            continue
        g0, b0, yy = moments(x[:, :, :N_BASE], y)
        gm1, bm1, _ = moments(x[:, :, :N_FULL], y)
        for local, pid in enumerate(valid):
            i = index[pid]
            gene_names = sorted(definitions[pid]["input_genes"])
            if len(gene_names) > max_genes:
                raise ValueError("marker array width is smaller than a program input set")
            gene_columns = [engine.gx[g] for g in gene_names]
            expression = engine.c[:, gene_columns].toarray()
            if mode == "raw":
                marker_x = np.log1p(expression)
            else:
                background_library = np.expm1(x[:, local, 1])
                scale = 1e4 / np.maximum(background_library, 1.0)
                marker_x = np.log1p(expression * scale[:, None])
            x0 = x[:, local, :N_BASE]
            y0 = y[:, local]
            n = len(y0)
            base_gram[i], base_rhs[i], target_second[i] = g0[local], b0[local], yy[local]
            full_gram[i], full_rhs[i] = gm1[local], bm1[local]
            marker_cross[i, :len(gene_names)] = (x0.T @ marker_x / n).T
            marker_second[i, :len(gene_names)] = np.mean(marker_x * marker_x, axis=0)
            marker_rhs[i, :len(gene_names)] = marker_x.T @ y0 / n
            marker_count[i] = len(gene_names)
            statuses[i] = "COMPUTED"

    save_npz_atomic(
        npz_path,
        program_ids=np.asarray(candidate_ids),
        status=statuses,
        marker_count=marker_count,
        base_gram=base_gram,
        base_rhs=base_rhs,
        target_second=target_second,
        full_gram=full_gram,
        full_rhs=full_rhs,
        marker_cross=marker_cross,
        marker_second=marker_second,
        marker_rhs=marker_rhs,
    )
    receipt = {
        "status": "COMPLETED",
        "section": row["key"],
        "cohort": row["cohort"],
        "patient": row["patient"],
        "mode": mode,
        "programs_computed": int(np.sum(statuses == "COMPUTED")),
        "program_status_counts": dict(Counter(statuses.tolist())),
        "input_sha256": input_digest,
        "output_sha256": digest(npz_path),
        "elapsed_seconds": time.monotonic() - started,
    }
    atomic_json(receipt_path, receipt)
    return receipt


def augment(stats):
    """Build M0 + one-marker moments for every input gene in one program."""
    n = len(stats["marker_second"])
    gram = np.zeros((n, N_FULL, N_FULL))
    rhs = np.zeros((n, N_FULL))
    gram[:, :N_BASE, :N_BASE] = stats["base_gram"]
    gram[:, :N_BASE, N_BASE] = stats["marker_cross"]
    gram[:, N_BASE, :N_BASE] = stats["marker_cross"]
    gram[:, N_BASE, N_BASE] = stats["marker_second"]
    rhs[:, :N_BASE] = stats["base_rhs"]
    rhs[:, N_BASE] = stats["marker_rhs"]
    return gram, rhs


def model_mse(train_stats, test_stats, gene_index=None, full_program=False):
    if full_program:
        train_gram = train_stats["full_gram"]
        train_rhs = train_stats["full_rhs"]
        test_gram = test_stats["full_gram"]
        test_rhs = test_stats["full_rhs"]
    elif gene_index is None:
        train_gram = train_stats["base_gram"]
        train_rhs = train_stats["base_rhs"]
        test_gram = test_stats["base_gram"]
        test_rhs = test_stats["base_rhs"]
    else:
        tg, tr = augment(train_stats)
        vg, vr = augment(test_stats)
        beta = fit_from_moments(tg[gene_index], tr[gene_index])
        err = test_stats["target_second"] - 2 * beta @ vr[gene_index] + beta @ vg[gene_index] @ beta
        return max(0.0, float(err))
    beta = fit_from_moments(train_gram, train_rhs)
    err = test_stats["target_second"] - 2 * beta @ test_rhs + beta @ test_gram @ beta
    return max(0.0, float(err))


def mean_stats(stats_by_patient, patient_ids):
    return {key: np.mean([stats_by_patient[p][key] for p in patient_ids], axis=0)
            for key in ("base_gram", "base_rhs", "target_second", "full_gram", "full_rhs",
                        "marker_cross", "marker_second", "marker_rhs")}


def read_checkpoint_stats(sec_rows, candidate_ids):
    section_stats = defaultdict(list)
    checkpoint_rows = []
    for row in sec_rows:
        for mode in ("raw", "depth_normalized"):
            tag = f"{row['key']}_{mode}"
            npz_path = OUT / "section_moments" / f"{tag}.npz"
            receipt_path = OUT / "section_moments" / f"{tag}.json"
            rec = json.loads(receipt_path.read_text(encoding="utf-8"))
            if rec["status"] != "COMPLETED" or digest(npz_path) != rec["output_sha256"]:
                raise RuntimeError(f"invalid section checkpoint: {tag}")
            checkpoint_rows.append(rec)
            with np.load(npz_path, allow_pickle=False) as archive:
                if archive["program_ids"].tolist() != candidate_ids:
                    raise RuntimeError(f"program-order mismatch in {tag}")
                for i, pid in enumerate(candidate_ids):
                    status = str(archive["status"][i])
                    item = {"status": status}
                    for key in ("base_gram", "base_rhs", "target_second", "full_gram", "full_rhs",
                                "marker_cross", "marker_second", "marker_rhs"):
                        item[key] = archive[key][i].copy()
                    count = int(archive["marker_count"][i])
                    item["marker_cross"] = item["marker_cross"][:count]
                    item["marker_second"] = item["marker_second"][:count]
                    item["marker_rhs"] = item["marker_rhs"][:count]
                    section_stats[(mode, row["cohort"], row["patient"], pid)].append(item)
    patient_stats = {}
    for key, sections_for_patient in section_stats.items():
        statuses = {x["status"] for x in sections_for_patient}
        if statuses != {"COMPUTED"}:
            patient_stats[key] = {"status": ";".join(sorted(statuses))}
            continue
        patient_stats[key] = {"status": "COMPUTED"}
        for field in ("base_gram", "base_rhs", "target_second", "full_gram", "full_rhs",
                      "marker_cross", "marker_second", "marker_rhs"):
            patient_stats[key][field] = np.mean([x[field] for x in sections_for_patient], axis=0)
    return patient_stats, checkpoint_rows


def run_nested_comparison(candidate_ids, definitions, sec_rows, patient_stats, frozen_rows):
    genes_by_program = {pid: sorted(definitions[pid]["input_genes"]) for pid in candidate_ids}
    frozen = {(r["mode"], r["cohort"], r["patient"], r["program"]): r for r in frozen_rows}
    patients_by_cohort = defaultdict(list)
    for row in sec_rows:
        patients_by_cohort[row["cohort"]].append(row["patient"])
    patients_by_cohort = {c: sorted(set(ps)) for c, ps in patients_by_cohort.items()}
    patient_rows, trace_rows, check_rows = [], [], []

    for mode in ("raw", "depth_normalized"):
        for cohort in ("ST-CRC", "USZ"):
            patients = patients_by_cohort[cohort]
            for pid in candidate_ids:
                program_genes = genes_by_program[pid]
                for test_patient in patients:
                    key = (mode, cohort, test_patient, pid)
                    test = patient_stats.get(key, {"status": "MISSING_PATIENT_STATS"})
                    reference = frozen.get(key)
                    if test["status"] != "COMPUTED" or reference is None or reference["status"] != "COMPUTED":
                        patient_rows.append({"mode": mode, "cohort": cohort, "program": pid,
                                             "patient": test_patient, "status": test["status"],
                                             "selected_marker": "", "n_input_genes": len(program_genes)})
                        continue
                    train_patients = [p for p in patients if p != test_patient]
                    if len(train_patients) < 3:
                        raise ValueError(f"too few outer-training patients for {cohort}")

                    # Full-program M1 and M0 are recomputed from the same section moments.
                    train_stats = mean_stats(patient_stats, [(mode, cohort, p, pid) for p in train_patients])
                    m0 = model_mse(train_stats, test)
                    m1 = model_mse(train_stats, test, full_program=True)
                    frozen_m0, frozen_m1 = float(reference["mse_M0"]), float(reference["mse_M1"])
                    check_rows.append({
                        "mode": mode, "cohort": cohort, "program": pid, "patient": test_patient,
                        "recomputed_mse_M0": m0, "frozen_mse_M0": frozen_m0,
                        "abs_diff_M0": abs(m0 - frozen_m0),
                        "recomputed_mse_M1": m1, "frozen_mse_M1": frozen_m1,
                        "abs_diff_M1": abs(m1 - frozen_m1),
                        "pass": str(np.isclose(m0, frozen_m0, rtol=1e-7, atol=1e-8) and
                                     np.isclose(m1, frozen_m1, rtol=1e-7, atol=1e-8)).lower(),
                    })

                    inner_errors = []
                    for inner_patient in train_patients:
                        inner_train = [p for p in train_patients if p != inner_patient]
                        inner_train_stats = mean_stats(patient_stats,
                            [(mode, cohort, p, pid) for p in inner_train])
                        inner_test = patient_stats[(mode, cohort, inner_patient, pid)]
                        gene_gram, gene_rhs = augment(inner_train_stats)
                        beta = fit_from_moments(gene_gram, gene_rhs)
                        test_gram, test_rhs = augment(inner_test)
                        errors = (inner_test["target_second"] - 2 * np.einsum("pi,pi->p", beta, test_rhs) +
                                  np.einsum("pi,pij,pj->p", beta, test_gram, beta))
                        inner_errors.append(np.maximum(errors, 0.0))
                    inner_mean = np.mean(inner_errors, axis=0)
                    order = np.argsort(inner_mean, kind="stable")
                    best_idx = int(order[0])
                    selected_gene = program_genes[best_idx]
                    outer_train_stats = mean_stats(patient_stats,
                        [(mode, cohort, p, pid) for p in train_patients])
                    single_mse = model_mse(outer_train_stats, test, gene_index=best_idx)
                    gain_over_marker = (100 * (single_mse - m1) / single_mse
                                        if single_mse > 1e-12 else float("nan"))
                    marker_gain_over_m0 = (100 * (m0 - single_mse) / m0 if m0 > 1e-12 else float("nan"))
                    program_gain_over_m0 = (100 * (m0 - m1) / m0 if m0 > 1e-12 else float("nan"))
                    patient_rows.append({
                        "mode": mode, "cohort": cohort, "program": pid, "patient": test_patient,
                        "status": "COMPUTED", "selected_marker": selected_gene,
                        "selection_inner_mean_mse": float(inner_mean[best_idx]),
                        "selection_runner_up": program_genes[int(order[1])] if len(order) > 1 else "",
                        "n_inner_patients": len(inner_errors), "n_input_genes": len(program_genes),
                        "mse_M0": m0, "mse_M1": m1, "mse_single_marker": single_mse,
                        "delta_M1_over_M0_pct": program_gain_over_m0,
                        "delta_single_over_M0_pct": marker_gain_over_m0,
                        "M1_improvement_over_selected_marker_pct": gain_over_marker,
                    })
                    trace_rows.append({
                        "mode": mode, "cohort": cohort, "program": pid, "outer_patient": test_patient,
                        "selected_marker": selected_gene,
                        "selected_inner_mean_mse": float(inner_mean[best_idx]),
                        "runner_up": program_genes[int(order[1])] if len(order) > 1 else "",
                        "runner_up_inner_mean_mse": float(inner_mean[int(order[1])]) if len(order) > 1 else "",
                        "n_inner_patients": len(inner_errors), "selection_scope": "outer_training_patients_only",
                    })
    return patient_rows, trace_rows, check_rows


def summarize(patient_rows):
    computed = [r for r in patient_rows if r.get("status") == "COMPUTED"]
    by_program = defaultdict(list)
    for row in computed:
        by_program[(row["mode"], row["cohort"], row["program"])].append(row)
    program_summary = []
    for (mode, cohort, pid), values in sorted(by_program.items()):
        gains = np.array([float(r["M1_improvement_over_selected_marker_pct"]) for r in values])
        program_summary.append({
            "mode": mode, "cohort": cohort, "program": pid, "n_patients": len(values),
            "median_M1_improvement_over_marker_pct": float(np.nanmedian(gains)),
            "n_patients_M1_lower_MSE": int(np.sum(gains > 0)),
            "n_patients_marker_lower_MSE": int(np.sum(gains < 0)),
            "min_patient_gain_pct": float(np.nanmin(gains)),
            "max_patient_gain_pct": float(np.nanmax(gains)),
            "selected_markers": ";".join(sorted({r["selected_marker"] for r in values})),
        })
    by_group = defaultdict(list)
    for row in program_summary:
        by_group[(row["mode"], row["cohort"])].append(row)
    cohort_summary = []
    for (mode, cohort), values in sorted(by_group.items()):
        medians = np.array([float(r["median_M1_improvement_over_marker_pct"]) for r in values])
        cohort_summary.append({
            "mode": mode, "cohort": cohort, "n_programs": len(values),
            "median_program_gain_over_marker_pct": float(np.median(medians)),
            "programs_with_positive_median_gain": int(np.sum(medians > 0)),
            "programs_with_negative_median_gain": int(np.sum(medians < 0)),
            "programs_with_zero_median_gain": int(np.sum(medians == 0)),
            "interpretation": "descriptive only; programs overlap and were selected on these cohorts",
        })
    selected = defaultdict(list)
    for row in computed:
        selected[(row["mode"], row["cohort"], row["selected_marker"])].append(row)
    gene_frequency = []
    for (mode, cohort, gene), values in sorted(selected.items()):
        gene_frequency.append({
            "mode": mode, "cohort": cohort, "gene_symbol": gene,
            "n_program_patient_folds_selected": len(values),
            "n_distinct_programs": len({r["program"] for r in values}),
            "n_distinct_outer_patients": len({r["patient"] for r in values}),
        })
    return program_summary, cohort_summary, gene_frequency


def main():
    if shutil.disk_usage(ROOT).free < MAX_FREE_BYTES:
        raise RuntimeError("Storage reserve below 1.2 TB")
    sec_rows = sections(["ST-CRC", "USZ"])
    stage_contract = json.loads((STAGE / "run_contract.json").read_text(encoding="utf-8"))
    if stage_contract["programs"] != 6870 or stage_contract["modes"] != ["raw", "depth_normalized"]:
        raise ValueError("unexpected completed GO stage contract")
    definitions = json.loads(DEFINITION_PATH.read_text(encoding="utf-8"))
    candidate_rows = read_tsv(STAGE / "candidate_context.tsv")
    candidate_ids = sorted(r["program"] for r in candidate_rows if r["molecular_robust_gt20"].lower() == "true")
    if len(candidate_ids) != 45 or not set(candidate_ids) <= set(definitions):
        raise ValueError(f"expected the frozen 45 candidates, found {len(candidate_ids)}")
    manifest = make_manifest(candidate_ids, definitions)
    input_hashes = collect_inputs(sec_rows)
    code_paths = [Path(__file__), ROOT / "r16/recovery/gobp_stage.py",
                  ROOT / "r16/recovery/corrected.py", ROOT / "r16/recovery/repair_pipeline.py"]
    contract = frozen_contract(manifest, input_hashes)

    if not OUT.exists():
        OUT.mkdir(parents=True)
        atomic_json(OUT / "program_manifest.json", manifest)
        atomic_json(OUT / "run_contract.json", contract)
    else:
        prior_contract = OUT / "run_contract.json"
        manifest_path = OUT / "program_manifest.json"
        if not prior_contract.is_file() or not manifest_path.is_file():
            raise RuntimeError("output directory exists without a frozen run contract")
        if json.loads(prior_contract.read_text(encoding="utf-8")) != contract:
            raise RuntimeError("frozen run contract changed; use a new output directory")
        if json.loads(manifest_path.read_text(encoding="utf-8")) != manifest:
            raise RuntimeError("program manifest mismatch")
        closed_receipt = OUT / "receipt.json"
        if closed_receipt.is_file():
            receipt = json.loads(closed_receipt.read_text(encoding="utf-8"))
            if receipt.get("run_contract_sha256") != digest(prior_contract):
                raise RuntimeError("closed receipt does not match the frozen contract")
            for relative, expected in receipt.get("outputs_sha256", {}).items():
                path = OUT / relative
                if not path.is_file() or digest(path) != expected:
                    raise RuntimeError(f"closed output hash mismatch: {relative}")
            if receipt.get("status") != "COMPLETED_EXPLORATORY_POSTSELECTED":
                raise RuntimeError("unexpected closed-run status")
            print(json.dumps({"status": receipt["status"], "programs": receipt["n_programs"],
                              "patient_rows": receipt["n_single_marker_patient_rows"],
                              "outputs": len(receipt["outputs_sha256"]), "verified_existing_run": True}, ensure_ascii=False))
            return

    genes = json.loads((SOURCE / "genes.json").read_text(encoding="utf-8"))
    gene_axes = axes()
    tasks_done = 0
    for row in sec_rows:
        if shutil.disk_usage(ROOT).free < MAX_FREE_BYTES:
            raise RuntimeError("Storage reserve fell below 1.2 TB; stopping")
        counts_path, source_gene_path = section_paths(row)
        if digest(counts_path) != row["counts_sha256"]:
            raise RuntimeError(f"source count hash changed: {row['key']}")
        counts = np.load(counts_path, mmap_mode="r")
        source_genes = set(json.loads(source_gene_path.read_text(encoding="utf-8")))
        engine = FeatureEngine(counts, genes, gene_axes, {pid: definitions[pid] for pid in candidate_ids},
                               unavailable=set(genes) - source_genes)
        input_digest = json_hash({
            "counts": input_hashes[counts_path.relative_to(ROOT).as_posix()],
            "source_genes": input_hashes[source_gene_path.relative_to(ROOT).as_posix()],
            "program_manifest": json_hash(manifest),
            "section_row": row,
        })
        for mode in ("raw", "depth_normalized"):
            rec = write_section_mode(engine, row, mode, candidate_ids, definitions, input_digest)
            tasks_done += 1
            print(f"checkpoint {tasks_done}/44 {row['key']} {mode} {rec['programs_computed']}/45", flush=True)
    frozen_rows = read_tsv(STAGE / "molecular_patients.tsv")
    frozen_rows = [r for r in frozen_rows if r["program"] in set(candidate_ids)]
    patient_stats, checkpoints = read_checkpoint_stats(sec_rows, candidate_ids)
    patient_rows, trace_rows, check_rows = run_nested_comparison(
        candidate_ids, definitions, sec_rows, patient_stats, frozen_rows)
    if not check_rows or not all(r["pass"] == "true" for r in check_rows):
        raise RuntimeError("recomputed M0/M1 scores did not match the frozen GO stage")
    if len([r for r in patient_rows if r.get("status") == "COMPUTED"]) != 1350:
        raise RuntimeError("expected 45 programs x 15 patients x 2 modes")

    program_summary, cohort_summary, gene_frequency = summarize(patient_rows)
    write_tsv(OUT / "single_marker_patients.tsv", patient_rows)
    write_tsv(OUT / "selection_trace.tsv", trace_rows)
    write_tsv(OUT / "molecular_reconstruction_check.tsv", check_rows)
    write_tsv(OUT / "single_marker_program_summary.tsv", program_summary)
    write_tsv(OUT / "single_marker_cohort_summary.tsv", cohort_summary)
    write_tsv(OUT / "selected_gene_frequency.tsv", gene_frequency)

    outputs = sorted(p for p in OUT.rglob("*") if p.is_file() and p.name != "receipt.json")
    receipt = {
        "status": "COMPLETED_EXPLORATORY_POSTSELECTED",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "audit_date": "2026-09-24",
        "run_contract_sha256": digest(OUT / "run_contract.json"),
        "input_sha256": input_hashes,
        "n_programs": len(candidate_ids),
        "n_sections": len(sec_rows),
        "n_patients": {c: len({r["patient"] for r in sec_rows if r["cohort"] == c}) for c in ("ST-CRC", "USZ")},
        "n_section_mode_checkpoints": len(checkpoints),
        "n_single_marker_patient_rows": len([r for r in patient_rows if r.get("status") == "COMPUTED"]),
        "n_reconstruction_rows": len(check_rows),
        "molecular_reconstruction_max_abs_diff_M0": max(float(r["abs_diff_M0"]) for r in check_rows),
        "molecular_reconstruction_max_abs_diff_M1": max(float(r["abs_diff_M1"]) for r in check_rows),
        "candidate_selection_note": "The 45 were selected using all patients in the completed GO BP stage; this run is a conditional, exploratory comparison, not an independent confirmatory test.",
        "compute_limits": {"GPU_used": False, "external_data_acquired": False, "old_models_retrained": False,
                           "free_bytes_after": shutil.disk_usage(ROOT).free},
        "outputs_sha256": {p.relative_to(OUT).as_posix(): digest(p) for p in outputs},
    }
    if receipt["compute_limits"]["free_bytes_after"] < MAX_FREE_BYTES:
        raise RuntimeError("Storage reserve below 1.2 TB at closeout")
    atomic_json(OUT / "receipt.json", receipt)
    print(json.dumps({"status": receipt["status"], "programs": receipt["n_programs"],
                      "patient_rows": receipt["n_single_marker_patient_rows"],
                      "max_abs_diff_M0": receipt["molecular_reconstruction_max_abs_diff_M0"],
                      "max_abs_diff_M1": receipt["molecular_reconstruction_max_abs_diff_M1"],
                      "outputs": len(outputs)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
