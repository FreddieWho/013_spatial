#!/usr/bin/env python3
"""Create a reproducible census of robust GO BP gene overlap and local GT labels."""

import csv
import hashlib
import json
import shutil
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "infra/gobp_followup_20260924"
GO_DEFS = ROOT / "infra/repair_20260921/go_representatives.json"
GO_CONTEXT = ROOT / "infra/gobp_stage_20260922/candidate_context.tsv"
GT_SOURCES = ROOT / "infra/structure-registry/gt_source_audit.tsv"
OUTER_SPLITS = ROOT / "infra/structure-registry/outer_splits.tsv"
GSE_RAW = ROOT / "data/GEO/GSE175540/raw"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def jaccard(left, right):
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def summarize_gt():
    sources = read_tsv(GT_SOURCES)
    splits = read_tsv(OUTER_SPLITS)
    by_physical = defaultdict(list)
    for row in splits:
        by_physical[row["physical_unit_id"]].append(row)

    detail_rows = []
    groups = defaultdict(list)
    for source in sources:
        if source["audit_status"] != "AUDITABLE_GT":
            continue
        physical_id = source["physical_unit_id"]
        matches = by_physical.get(physical_id, [])
        if len(matches) > 1:
            raise ValueError(f"physical unit maps to multiple frozen R02 rows: {physical_id}")
        split = matches[0] if matches else {}
        source_path = ROOT / source["path"]
        local_hash = sha256(source_path) if source_path.is_file() else ""
        hash_matches = bool(local_hash) and local_hash == source["sha256"]
        detail = {
            "gt_source_id": source["gt_source_id"],
            "source_class": source["source_class"],
            "structure_id": source["structure_id"],
            "audit_status": source["audit_status"],
            "source_path": source["path"],
            "registry_sha256": source["sha256"],
            "local_sha256": local_hash,
            "local_sha_matches": str(hash_matches).lower(),
            "physical_unit_id": physical_id,
            "logical_unit_id": split.get("logical_unit_id", ""),
            "patient_id": split.get("patient_id", ""),
            "block_id": split.get("block_id", ""),
            "identity_granularity": split.get("identity_granularity", ""),
            "primary_role": split.get("primary_role", "UNMAPPED"),
            "block_level_eligible": split.get("block_level_eligible", "unknown"),
            "physical_link_status": source["physical_link_status"],
            "provenance_status": source["provenance_status"],
            "same_assay_status": source["same_assay_status"],
        }
        detail_rows.append(detail)
        groups[(detail["logical_unit_id"] or "UNMAPPED", source["source_class"], source["structure_id"])].append(detail)

    fields = list(detail_rows[0]) if detail_rows else []
    write_tsv(OUT / "gt_sources.tsv", fields, sorted(detail_rows, key=lambda r: (r["logical_unit_id"], r["structure_id"], r["gt_source_id"])))

    summary_rows = []
    for (lineage, source_class, structure), items in sorted(groups.items()):
        summary_rows.append({
            "logical_unit_id": lineage,
            "source_class": source_class,
            "structure_id": structure,
            "n_auditable_sources": len(items),
            "n_physical_units": len({r["physical_unit_id"] for r in items if r["physical_unit_id"]}),
            "n_patient_ids": len({r["patient_id"] for r in items if r["patient_id"]}),
            "roles": ";".join(sorted({r["primary_role"] for r in items})),
            "n_block_level_eligible": sum(r["block_level_eligible"] == "yes" for r in items),
            "n_not_block_level_eligible": sum(r["block_level_eligible"] == "no" for r in items),
            "n_unmapped": sum(r["primary_role"] == "UNMAPPED" for r in items),
            "n_local_hash_matches": sum(r["local_sha_matches"] == "true" for r in items),
        })
    summary_fields = list(summary_rows[0]) if summary_rows else []
    write_tsv(OUT / "gt_lineage_summary.tsv", summary_fields, summary_rows)

    cohorts = defaultdict(list)
    for row in detail_rows:
        cohorts[row["logical_unit_id"] or "UNMAPPED"].append(row)
    cohort_rows = []
    for lineage, items in sorted(cohorts.items()):
        cohort_rows.append({
            "logical_unit_id": lineage,
            "n_auditable_sources": len(items),
            "structures": ";".join(sorted({r["structure_id"] for r in items})),
            "n_physical_units": len({r["physical_unit_id"] for r in items if r["physical_unit_id"]}),
            "n_patient_ids": len({r["patient_id"] for r in items if r["patient_id"]}),
            "roles": ";".join(sorted({r["primary_role"] for r in items})),
            "n_block_level_eligible": sum(r["block_level_eligible"] == "yes" for r in items),
            "n_not_block_level_eligible": sum(r["block_level_eligible"] == "no" for r in items),
            "n_local_hash_matches": sum(r["local_sha_matches"] == "true" for r in items),
        })
    write_tsv(OUT / "gt_cohort_summary.tsv", list(cohort_rows[0]), cohort_rows)

    raw_files = sorted(p for p in GSE_RAW.rglob("*") if p.is_file()) if GSE_RAW.exists() else []
    categories = {
        "spot_annotation_csv_gz": [p for p in raw_files if "annot" in p.name.lower() and p.name.lower().endswith(".csv.gz")],
        "filtered_feature_bc_matrix_h5": [p for p in raw_files if p.name.endswith("filtered_feature_bc_matrix.h5")],
        "tissue_positions_csv_gz": [p for p in raw_files if p.name.endswith("tissue_positions_list.csv.gz")],
    }
    inventory_rows = []
    for category, paths in categories.items():
        inventory_rows.append({
            "resource_class": category,
            "n_local_files": len(paths),
            "total_bytes": sum(p.stat().st_size for p in paths),
            "relative_paths": ";".join(p.relative_to(ROOT).as_posix() for p in paths),
        })
    write_tsv(OUT / "gse175540_local_inventory.tsv", list(inventory_rows[0]), inventory_rows)
    return detail_rows, summary_rows, cohort_rows, inventory_rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    definitions = json.loads(GO_DEFS.read_text(encoding="utf-8"))
    context = read_tsv(GO_CONTEXT)
    robust_rows = [row for row in context if row["molecular_robust_gt20"].lower() == "true"]
    robust_ids = sorted(row["program"] for row in robust_rows)
    if len(robust_ids) != 45 or not set(robust_ids) <= set(definitions):
        raise ValueError(f"expected 45 frozen robust programs, found {len(robust_ids)}")
    context_by_id = {row["program"]: row for row in robust_rows}

    program_sets = {}
    program_rows = []
    for program_id in robust_ids:
        definition = definitions[program_id]
        inputs = set(definition["input_genes"])
        readouts = set(definition["readout_genes"])
        program_sets[program_id] = {"input": inputs, "readout": readouts, "combined": inputs | readouts}
        row = context_by_id[program_id]
        program_rows.append({
            "program_id": program_id,
            "covered_size": definition["covered_size"],
            "raw_size": definition["raw_size"],
            "input_gene_count": len(inputs),
            "readout_gene_count": len(readouts),
            "combined_distinct_gene_count": len(inputs | readouts),
            "molecular_minimum_median_improvement": row["M1_minimum_median_improvement"],
            "molecular_all_conditions_complete": row["M1_all_conditions_complete"],
            "definition_hash": definition["definition_hash"],
        })
    write_tsv(OUT / "go45_programs.tsv", list(program_rows[0]), program_rows)

    frequencies = {mode: Counter() for mode in ("input", "readout", "combined")}
    for sets in program_sets.values():
        for mode, genes in sets.items():
            frequencies[mode].update(genes)
    all_genes = sorted(set().union(*(s["combined"] for s in program_sets.values())))
    frequency_rows = [{
        "gene_symbol": gene,
        "input_program_count": frequencies["input"][gene],
        "readout_program_count": frequencies["readout"][gene],
        "combined_program_count": frequencies["combined"][gene],
    } for gene in all_genes]
    write_tsv(OUT / "go45_gene_frequency.tsv", list(frequency_rows[0]), frequency_rows)

    pair_rows = []
    for left, right in combinations(robust_ids, 2):
        result = {"program_a": left, "program_b": right}
        for mode in ("input", "readout", "combined"):
            a, b = program_sets[left][mode], program_sets[right][mode]
            result[f"{mode}_intersection_n"] = len(a & b)
            result[f"{mode}_jaccard"] = f"{jaccard(a, b):.8f}"
        pair_rows.append(result)
    write_tsv(OUT / "go45_pairwise_overlap.tsv", list(pair_rows[0]), pair_rows)

    gt_rows, gt_summary, gt_cohorts, gse_inventory = summarize_gt()
    combined_scores = [float(row["combined_jaccard"]) for row in pair_rows]
    input_scores = [float(row["input_jaccard"]) for row in pair_rows]
    readout_scores = [float(row["readout_jaccard"]) for row in pair_rows]
    covered = [int(r["covered_size"]) for r in program_rows]
    inputs_occurrences = sum(len(s["input"]) for s in program_sets.values())
    readout_occurrences = sum(len(s["readout"]) for s in program_sets.values())
    combined_occurrences = sum(len(s["combined"]) for s in program_sets.values())
    top_input = sorted(frequencies["input"].items(), key=lambda x: (-x[1], x[0]))[:10]

    outputs = sorted(p for p in OUT.iterdir() if p.is_file() and p.name != "receipt.json")
    receipt = {
        "audit_date": "2026-09-24",
        "script": Path(__file__).relative_to(ROOT).as_posix(),
        "script_sha256": sha256(Path(__file__)),
        "inputs_sha256": {p.relative_to(ROOT).as_posix(): sha256(p) for p in (GO_DEFS, GO_CONTEXT, GT_SOURCES, OUTER_SPLITS)},
        "method": {
            "candidate_selection": "candidate_context.tsv molecular_robust_gt20 == true; definitions are the frozen GO BP representatives",
            "gene_overlap": "per-program input, readout, and input union readout sets; pairwise Jaccard = intersection / union",
            "gt_join": "AUDITABLE_GT source rows joined to the frozen R02 outer split on physical_unit_id; source file SHA-256 independently recomputed",
            "scope": "local inventory only; no model fit, external download, or GPU use",
        },
        "go45": {
            "n_programs": len(program_rows),
            "covered_size_min_median_max": [min(covered), median(covered), max(covered)],
            "input_gene_count_min_median_max": [min(len(s['input']) for s in program_sets.values()), median(len(s['input']) for s in program_sets.values()), max(len(s['input']) for s in program_sets.values())],
            "readout_gene_count_min_median_max": [min(len(s['readout']) for s in program_sets.values()), median(len(s['readout']) for s in program_sets.values()), max(len(s['readout']) for s in program_sets.values())],
            "input_occurrences_unique_genes": [inputs_occurrences, len(frequencies["input"])],
            "readout_occurrences_unique_genes": [readout_occurrences, len(frequencies["readout"])],
            "combined_occurrences_unique_genes": [combined_occurrences, len(frequencies["combined"])],
            "pair_count": len(pair_rows),
            "pairwise_jaccard_median": {"input": median(input_scores), "readout": median(readout_scores), "combined": median(combined_scores)},
            "combined_pairs_ge_0_25": sum(score >= 0.25 for score in combined_scores),
            "combined_pairs_ge_0_50": sum(score >= 0.50 for score in combined_scores),
            "top_input_genes": [{"symbol": gene, "program_count": n} for gene, n in top_input],
        },
        "gt": {
            "auditable_source_rows": len(gt_rows),
            "local_sha256_matches": sum(r["local_sha_matches"] == "true" for r in gt_rows),
            "mapped_to_r02": sum(r["primary_role"] != "UNMAPPED" for r in gt_rows),
            "lineage_structure_groups": len(gt_summary),
            "lineages": len(gt_cohorts),
            "gse175540_local_file_counts": {r["resource_class"]: r["n_local_files"] for r in gse_inventory},
        },
        "outputs_sha256": {p.name: sha256(p) for p in outputs},
        "free_bytes_after": shutil.disk_usage(ROOT).free,
    }
    (OUT / "receipt.json").write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"out": OUT.relative_to(ROOT).as_posix(), "go45": receipt["go45"], "gt": receipt["gt"], "free_bytes_after": receipt["free_bytes_after"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
