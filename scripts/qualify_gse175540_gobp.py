#!/usr/bin/env python3
"""Qualify local GSE175540 TLS labels and GO BP feature coverage, without scoring."""

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "infra/gobp_external_qualification_20260925"
RAW = ROOT / "data/GEO/GSE175540/raw"
GT_SOURCES = ROOT / "infra/structure-registry/gt_source_audit.tsv"
OUTER_SPLITS = ROOT / "infra/structure-registry/outer_splits.tsv"
GO_DEFS = ROOT / "infra/repair_20260921/go_representatives.json"
GO_CONTEXT = ROOT / "infra/gobp_stage_20260922/candidate_context.tsv"
CONTRACT = RUN / "run_contract.json"
RECEIPT = RUN / "receipt.json"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_hash(value):
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    rows = list(rows)
    if not rows:
        raise ValueError(f"refusing to write an empty table: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def decode(values):
    return [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values]


def one_match(pattern):
    matches = sorted(RAW.glob(pattern))
    if len(matches) != 1:
        raise ValueError(f"expected one local file for {pattern}, found {len(matches)}")
    return matches[0]


def candidates():
    definitions = json.loads(GO_DEFS.read_text(encoding="utf-8"))
    context = [
        row for row in read_tsv(GO_CONTEXT)
        if row["molecular_robust_gt20"].lower() == "true"
    ]
    ids = sorted(row["program"] for row in context)
    if len(ids) != 45 or not set(ids) <= set(definitions):
        raise ValueError(f"expected 45 frozen programs, found {len(ids)}")
    manifest = []
    for pid in ids:
        definition = definitions[pid]
        inputs = sorted(set(definition["input_genes"]))
        readouts = sorted(set(definition["readout_genes"]))
        if set(inputs) & set(readouts):
            raise ValueError(f"input/readout overlap in frozen definition: {pid}")
        manifest.append({
            "program_id": pid,
            "definition_hash": definition["definition_hash"],
            "input_genes": inputs,
            "readout_genes": readouts,
        })
    return definitions, manifest


def source_manifest():
    sources = read_tsv(GT_SOURCES)
    splits = read_tsv(OUTER_SPLITS)
    splits_by_id = defaultdict(list)
    for split in splits:
        splits_by_id[split["physical_unit_id"]].append(split)

    selected = []
    for source in sources:
        physical = source["physical_unit_id"]
        if not (
            source["audit_status"] == "AUDITABLE_GT"
            and source["structure_id"] == "TLS"
            and physical.startswith("external_geo::GSM")
            and source["path"].startswith("data/GEO/GSE175540/")
        ):
            continue
        matches = splits_by_id.get(physical, [])
        if len(matches) != 1:
            raise ValueError(f"expected one frozen R-02 row for {physical}, found {len(matches)}")
        split = matches[0]
        if split["logical_unit_id"] != "GEO::GSE175540":
            continue
        gsm = physical.split("::", 1)[1]
        annotation = ROOT / source["path"]
        if not annotation.name.startswith(gsm + "_"):
            raise ValueError(f"annotation accession/path mismatch: {physical} {annotation}")
        selected.append({
            "gt_source_id": source["gt_source_id"],
            "source_path": source["path"],
            "source_sha256": source["sha256"],
            "physical_unit_id": physical,
            "patient_id": split["patient_id"],
            "logical_unit_id": split["logical_unit_id"],
            "primary_role": split["primary_role"],
            "block_id": split["block_id"],
            "block_level_eligible": split["block_level_eligible"],
            "gsm": gsm,
            "matrix_path": one_match(gsm + "*_filtered_feature_bc_matrix.h5").relative_to(ROOT).as_posix(),
            "positions_path": one_match(gsm + "*_tissue_positions_list.csv.gz").relative_to(ROOT).as_posix(),
        })
    selected.sort(key=lambda row: row["gsm"])
    if len(selected) != 18:
        raise ValueError(f"expected 18 auditable GSE175540 TLS sources, found {len(selected)}")
    if len({row["patient_id"] for row in selected}) != 18:
        raise ValueError("GSE175540 TLS source rows do not map to 18 distinct patients")
    if any(row["primary_role"] != "external_validation" for row in selected):
        raise ValueError("one or more GSE175540 labels lack the frozen external_validation role")
    if any(row["block_level_eligible"] != "no" or row["block_id"] for row in selected):
        raise ValueError("unexpected GSE175540 block-level eligibility")
    return selected


def current_inputs():
    _, manifest = candidates()
    selected = source_manifest()
    paths = [GT_SOURCES, OUTER_SPLITS, GO_DEFS, GO_CONTEXT]
    for row in selected:
        paths.extend([
            ROOT / row["source_path"],
            ROOT / row["matrix_path"],
            ROOT / row["positions_path"],
        ])
    paths.append(Path(__file__).resolve())
    hashes = {
        path.relative_to(ROOT).as_posix(): sha256(path)
        for path in sorted(set(paths))
    }
    return selected, manifest, hashes


def freeze():
    if CONTRACT.exists() or RECEIPT.exists():
        raise FileExistsError(f"qualification run already exists: {RUN}")
    selected, manifest, hashes = current_inputs()
    RUN.mkdir(parents=True, exist_ok=True)
    contract = {
        "version": 1,
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "FROZEN_BEFORE_QUALIFICATION_AUDIT",
        "scope": "Local source, label-key, matrix-integrity, and exact GO-symbol eligibility audit for the 18 R-02 GSE175540 TLS sources. No expression score, model fit, hypothesis test, download, or GPU use.",
        "frozen_candidate_count": len(manifest),
        "candidate_manifest_sha256": json_hash(manifest),
        "source_manifest_sha256": json_hash(selected),
        "rules": {
            "labels": "Only explicit TLS and NO_TLS values are classes. Blank annotation values remain unknown and are never treated as negative.",
            "sample_key": "Accession must match the frozen R-02 physical unit and the annotation, filtered matrix, and tissue-position filenames. Explicitly labeled barcodes must map to both matrix and in-tissue coordinates; missing labels are reported and not imputed.",
            "gene_symbols": "Use exact frozen program symbols against the matrix features/name entries restricted to Gene Expression and GRCh38. No synonym mapping. A symbol with multiple feature IDs is ambiguous; do not silently select or aggregate it.",
            "strict_two_class_sample": "A patient is eligible for a two-class comparison only when its explicit TLS and NO_TLS labels both exist and every explicit label barcode maps to the filtered matrix and in-tissue coordinates.",
            "strict_program": "A program is eligible only when all frozen input and readout symbols are present and each occurs as exactly one matrix feature in every R-02 patient sample.",
            "independence": "Patient is the unit. GSE175540 has no block-level eligibility; no block-level inference is allowed.",
        },
        "selected_sources": selected,
        "input_sha256": hashes,
        "computation": "Read-only sequential inspection of local gzip CSV and 10x HDF5 metadata/sparse counts; no derived matrices are saved.",
    }
    CONTRACT.write_text(json.dumps(contract, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": contract["status"], "patients": len(selected), "programs": len(manifest), "inputs": len(hashes)}, ensure_ascii=False))


def verify_contract():
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    selected, manifest, hashes = current_inputs()
    if hashes != contract["input_sha256"]:
        raise RuntimeError("frozen qualification input/code hashes changed")
    if json_hash(selected) != contract["source_manifest_sha256"]:
        raise RuntimeError("frozen GSE175540 source manifest changed")
    if json_hash(manifest) != contract["candidate_manifest_sha256"]:
        raise RuntimeError("frozen GO candidate manifest changed")
    return contract, selected, manifest


def load_labels(path):
    with gzip.open(path, "rt", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not {"Barcode", "TLS_2_cat"} <= set(reader.fieldnames or []):
            raise ValueError(f"unexpected TLS annotation schema: {path}")
        rows = list(reader)
    barcodes = [row["Barcode"].strip() for row in rows]
    if any(not barcode for barcode in barcodes) or len(barcodes) != len(set(barcodes)):
        raise ValueError(f"empty or duplicated annotation barcode: {path}")
    values = [(row["TLS_2_cat"] or "").strip() for row in rows]
    return dict(zip(barcodes, values))


def load_positions(path):
    result = {}
    with gzip.open(path, "rt", newline="", encoding="utf-8-sig") as handle:
        for row in csv.reader(handle):
            if len(row) != 6:
                raise ValueError(f"expected six tissue-position columns in {path}")
            barcode = row[0].strip()
            if not barcode or barcode in result or row[1] not in {"0", "1"}:
                raise ValueError(f"invalid/duplicated tissue-position key in {path}: {barcode}")
            float(row[4])
            float(row[5])
            result[barcode] = {"in_tissue": row[1], "array_row": row[2], "array_col": row[3]}
    return result


def audit_sample(row, definitions, manifest):
    annotation_path = ROOT / row["source_path"]
    matrix_path = ROOT / row["matrix_path"]
    positions_path = ROOT / row["positions_path"]
    labels = load_labels(annotation_path)
    positions = load_positions(positions_path)
    label_counts = Counter(value if value else "UNKNOWN_BLANK" for value in labels.values())
    unexpected = sorted(set(labels.values()) - {"TLS", "NO_TLS", ""})
    if unexpected:
        raise ValueError(f"unexpected TLS label categories in {annotation_path}: {unexpected}")

    with h5py.File(matrix_path, "r") as handle:
        matrix = handle["matrix"]
        barcodes = decode(matrix["barcodes"][:])
        feature_group = matrix["features"]
        names = decode(feature_group["name"][:])
        feature_ids = decode(feature_group["id"][:])
        feature_types = decode(feature_group["feature_type"][:])
        genomes = decode(feature_group["genome"][:])
        shape = [int(value) for value in matrix["shape"][:]]
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]

    if len(barcodes) != len(set(barcodes)):
        raise ValueError(f"duplicated matrix barcodes: {matrix_path}")
    if len(names) != len(feature_ids) or len(names) != len(feature_types) or len(names) != len(genomes):
        raise ValueError(f"HDF5 feature metadata length mismatch: {matrix_path}")
    if shape != [len(names), len(barcodes)]:
        raise ValueError(f"HDF5 matrix shape disagrees with metadata: {matrix_path}")
    if len(indptr) != len(barcodes) + 1 or int(indptr[0]) != 0 or int(indptr[-1]) != len(data):
        raise ValueError(f"HDF5 sparse pointer length/end mismatch: {matrix_path}")
    if len(indices) != len(data) or np.any(np.diff(indptr) < 0):
        raise ValueError(f"invalid HDF5 sparse index arrays: {matrix_path}")
    if len(indices) and (int(indices.min()) < 0 or int(indices.max()) >= len(names)):
        raise ValueError(f"HDF5 feature index outside bounds: {matrix_path}")
    counts_valid = bool(np.issubdtype(data.dtype, np.integer) and np.all(data > 0))
    if not counts_valid:
        raise ValueError(f"matrix contains non-integer or non-positive stored UMI values: {matrix_path}")
    if any(feature_type != "Gene Expression" for feature_type in feature_types):
        raise ValueError(f"non-gene-expression feature row in filtered matrix: {matrix_path}")

    barcode_to_col = {barcode: index for index, barcode in enumerate(barcodes)}
    symbol_rows = defaultdict(list)
    eligible_feature_rows = 0
    for index, (name, feature_type, genome) in enumerate(zip(names, feature_types, genomes)):
        if feature_type == "Gene Expression" and genome == "GRCh38":
            symbol_rows[name].append(index)
            eligible_feature_rows += 1
    symbol_counts = {name: len(ix) for name, ix in symbol_rows.items()}
    all_candidate_symbols = {
        gene for program in manifest
        for gene in program["input_genes"] + program["readout_genes"]
    }
    duplicate_candidate_symbols = sorted(
        gene for gene in all_candidate_symbols if symbol_counts.get(gene, 0) > 1
    )

    unmatched_by_class = Counter()
    unmatched_positions_by_class = Counter()
    labels_in_tissue = 0
    matched_explicit = 0
    matched_tls = []
    matched_no_tls = []
    matrix_positions_missing = 0
    matrix_outside_tissue = 0
    cumulative = np.empty(len(data) + 1, dtype=np.int64)
    cumulative[0] = 0
    np.cumsum(data, dtype=np.int64, out=cumulative[1:])
    cell_totals = np.diff(cumulative[indptr])
    detected_per_cell = np.diff(indptr)

    for barcode, label in labels.items():
        klass = label or "UNKNOWN_BLANK"
        pos = positions.get(barcode)
        if pos is None:
            unmatched_positions_by_class[klass] += 1
        elif pos["in_tissue"] == "1":
            labels_in_tissue += 1
        if barcode not in barcode_to_col:
            unmatched_by_class[klass] += 1
            continue
        if pos is None:
            matrix_positions_missing += 1
        elif pos["in_tissue"] != "1":
            matrix_outside_tissue += 1
        col = barcode_to_col[barcode]
        if label == "TLS":
            matched_tls.append(col)
            matched_explicit += 1
        elif label == "NO_TLS":
            matched_no_tls.append(col)
            matched_explicit += 1

    feature_rows = {gene: symbol_rows.get(gene, []) for gene in all_candidate_symbols}
    program_rows = []
    collisions = []
    for program in manifest:
        pid = program["program_id"]
        input_genes = program["input_genes"]
        readout_genes = program["readout_genes"]
        input_missing = [gene for gene in input_genes if not feature_rows[gene]]
        readout_missing = [gene for gene in readout_genes if not feature_rows[gene]]
        input_ambiguous = [gene for gene in input_genes if len(feature_rows[gene]) > 1]
        readout_ambiguous = [gene for gene in readout_genes if len(feature_rows[gene]) > 1]
        program_rows.append({
            "gsm": row["gsm"],
            "patient_id": row["patient_id"],
            "program_id": pid,
            "input_gene_count": len(input_genes),
            "input_symbols_missing": len(input_missing),
            "input_missing_symbols": ";".join(input_missing),
            "input_ambiguous_symbols": ";".join(input_ambiguous),
            "readout_gene_count": len(readout_genes),
            "readout_symbols_missing": len(readout_missing),
            "readout_missing_symbols": ";".join(readout_missing),
            "readout_ambiguous_symbols": ";".join(readout_ambiguous),
            "exact_symbol_coverage": str(not input_missing and not readout_missing).lower(),
            "one_feature_per_symbol": str(not input_ambiguous and not readout_ambiguous).lower(),
        })
        for gene in sorted(set(input_ambiguous + readout_ambiguous)):
            ids = [feature_ids[index] for index in feature_rows[gene]]
            collisions.append({
                "gsm": row["gsm"],
                "patient_id": row["patient_id"],
                "gene_symbol_exact": gene,
                "feature_ids": ";".join(ids),
                "program_id": pid,
                "roles": ";".join(
                    role for role, genes in (("input", input_genes), ("readout", readout_genes))
                    if gene in genes
                ),
            })

    mapped_explicit = len(matched_tls) + len(matched_no_tls)
    matrix_label_fraction = mapped_explicit / max(label_counts["TLS"] + label_counts["NO_TLS"], 1)
    exact_file_key = (
        row["gsm"] in annotation_path.name
        and row["gsm"] in matrix_path.name
        and row["gsm"] in positions_path.name
    )
    both_classes = bool(matched_tls and matched_no_tls)
    all_explicit_labels_map = unmatched_by_class["TLS"] == 0 and unmatched_by_class["NO_TLS"] == 0
    all_explicit_positions_map = (
        unmatched_positions_by_class["TLS"] == 0
        and unmatched_positions_by_class["NO_TLS"] == 0
        and matrix_positions_missing == 0
        and matrix_outside_tissue == 0
    )
    strict_pair_eligible = bool(
        exact_file_key and both_classes and all_explicit_labels_map and all_explicit_positions_map
    )

    def median_or_blank(values):
        return float(np.median(values)) if len(values) else ""

    qc = {
        "gsm": row["gsm"],
        "patient_id": row["patient_id"],
        "frozen_role": row["primary_role"],
        "block_id": row["block_id"],
        "label_source_sha256_matches_registry": str(sha256(annotation_path) == row["source_sha256"]).lower(),
        "sample_accession_matches_all_file_keys": str(exact_file_key).lower(),
        "matrix_features": shape[0],
        "matrix_barcodes": shape[1],
        "gene_expression_grch38_feature_rows": eligible_feature_rows,
        "unique_exact_gene_symbols": len(symbol_rows),
        "candidate_duplicate_symbols": ";".join(duplicate_candidate_symbols),
        "matrix_nnz": len(data),
        "matrix_total_umi": int(data.sum(dtype=np.int64)),
        "matrix_counts_integrity": str(counts_valid).lower(),
        "annotation_rows": len(labels),
        "annotation_tls": label_counts["TLS"],
        "annotation_no_tls": label_counts["NO_TLS"],
        "annotation_blank_unknown": label_counts["UNKNOWN_BLANK"],
        "annotation_rows_unmatched_matrix": sum(unmatched_by_class.values()),
        "unmatched_tls": unmatched_by_class["TLS"],
        "unmatched_no_tls": unmatched_by_class["NO_TLS"],
        "explicit_label_matrix_coverage": f"{mapped_explicit}/{label_counts['TLS'] + label_counts['NO_TLS']}",
        "explicit_label_matrix_coverage_fraction": f"{matrix_label_fraction:.8f}",
        "annotation_rows_unmatched_positions": sum(unmatched_positions_by_class.values()),
        "explicit_labels_not_in_tissue": matrix_outside_tissue,
        "labels_in_tissue_positions": labels_in_tissue,
        "matched_tls_spots": len(matched_tls),
        "matched_no_tls_spots": len(matched_no_tls),
        "matched_tls_median_umi": median_or_blank(cell_totals[matched_tls]),
        "matched_no_tls_median_umi": median_or_blank(cell_totals[matched_no_tls]),
        "matched_tls_median_detected_features": median_or_blank(detected_per_cell[matched_tls]),
        "matched_no_tls_median_detected_features": median_or_blank(detected_per_cell[matched_no_tls]),
        "both_explicit_classes": str(both_classes).lower(),
        "all_explicit_labels_map": str(all_explicit_labels_map).lower(),
        "all_explicit_labels_have_in_tissue_coordinates": str(all_explicit_positions_map).lower(),
        "strict_two_class_patient_eligible": str(strict_pair_eligible).lower(),
    }
    return qc, program_rows, collisions


def run():
    contract, selected, manifest = verify_contract()
    if RECEIPT.exists():
        receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
        if receipt.get("frozen_contract_sha256") != sha256(CONTRACT):
            raise RuntimeError("qualification receipt does not match the frozen contract")
        for rel, expected in receipt["output_sha256"].items():
            path = ROOT / rel
            if not path.is_file() or sha256(path) != expected:
                raise RuntimeError(f"existing qualification output hash mismatch: {rel}")
        print(json.dumps({
            "status": receipt["status"],
            "verified_existing_run": True,
            "strict_patients": receipt["n_strict_two_class_patients"],
            "strict_programs": receipt["n_programs_with_unambiguous_exact_symbols_all_18"],
        }, ensure_ascii=False))
        return

    sample_qc = []
    program_coverage = []
    collisions = []
    definitions = json.loads(GO_DEFS.read_text(encoding="utf-8"))
    for row in selected:
        qc, rows, dupes = audit_sample(row, definitions, manifest)
        sample_qc.append(qc)
        program_coverage.extend(rows)
        collisions.extend(dupes)

    qc_by_gsm = {row["gsm"]: row for row in sample_qc}
    program_summary = []
    for program in manifest:
        pid = program["program_id"]
        rows = [row for row in program_coverage if row["program_id"] == pid]
        missing_samples = [row["gsm"] for row in rows if row["exact_symbol_coverage"] != "true"]
        ambiguous_samples = [row["gsm"] for row in rows if row["one_feature_per_symbol"] != "true"]
        strict_patients = [
            row["gsm"] for row in rows
            if row["one_feature_per_symbol"] == "true"
            and qc_by_gsm[row["gsm"]]["strict_two_class_patient_eligible"] == "true"
        ]
        ambiguous_symbols = sorted({
            symbol
            for row in rows
            for field in ("input_ambiguous_symbols", "readout_ambiguous_symbols")
            for symbol in row[field].split(";") if symbol
        })
        program_summary.append({
            "program_id": pid,
            "frozen_input_genes": len(program["input_genes"]),
            "frozen_readout_genes": len(program["readout_genes"]),
            "samples_missing_any_exact_symbol": len(missing_samples),
            "samples_with_ambiguous_feature_symbol": len(ambiguous_samples),
            "ambiguous_symbols": ";".join(ambiguous_symbols),
            "strict_two_class_patients_with_unambiguous_mapping": len(strict_patients),
            "strictly_eligible_across_all_18_sources": str(
                not missing_samples and not ambiguous_samples
            ).lower(),
        })

    write_tsv(RUN / "sample_qc.tsv", sample_qc)
    write_tsv(RUN / "program_coverage.tsv", program_coverage)
    write_tsv(RUN / "program_qualification_summary.tsv", program_summary)
    if collisions:
        write_tsv(RUN / "feature_symbol_collisions.tsv", collisions)

    strict_patients = sum(row["strict_two_class_patient_eligible"] == "true" for row in sample_qc)
    strict_programs = sum(row["strictly_eligible_across_all_18_sources"] == "true" for row in program_summary)
    fully_symbol_covered = sum(row["samples_missing_any_exact_symbol"] == 0 for row in program_summary)
    status = "QUALIFICATION_PARTIAL_WITH_EXCLUSIONS"
    output_paths = [
        RUN / "sample_qc.tsv",
        RUN / "program_coverage.tsv",
        RUN / "program_qualification_summary.tsv",
    ]
    if (RUN / "feature_symbol_collisions.tsv").is_file():
        output_paths.append(RUN / "feature_symbol_collisions.tsv")
    receipt = {
        "status": status,
        "frozen_contract_sha256": sha256(CONTRACT),
        "input_sha256": contract["input_sha256"],
        "n_r02_external_tls_sources": len(selected),
        "n_distinct_patients": len({row["patient_id"] for row in selected}),
        "n_strict_two_class_patients": strict_patients,
        "n_programs": len(program_summary),
        "n_programs_with_all_exact_symbols_all_18": fully_symbol_covered,
        "n_programs_with_unambiguous_exact_symbols_all_18": strict_programs,
        "n_feature_symbol_collision_rows": len(collisions),
        "expression_scores_run": False,
        "external_association_scores_run": False,
        "output_sha256": {
            path.relative_to(ROOT).as_posix(): sha256(path) for path in output_paths
        },
    }
    RECEIPT.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": status,
        "r02_patients": receipt["n_distinct_patients"],
        "strict_two_class_patients": strict_patients,
        "programs_with_all_symbols": fully_symbol_covered,
        "programs_unambiguous_all18": strict_programs,
        "feature_symbol_collision_rows": len(collisions),
        "expression_scores_run": False,
    }, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("freeze", "run"))
    args = parser.parse_args()
    if args.command == "freeze":
        freeze()
    else:
        run()


if __name__ == "__main__":
    main()
