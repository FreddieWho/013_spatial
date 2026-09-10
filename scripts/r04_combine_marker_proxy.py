#!/usr/bin/env python3
"""Combine CRC marker lists from three public databases by majority vote.

Sources (all versioned, provenance recorded):
  1. CellMarker 2.0 human CRC subset (literature-curated cell_name -> Symbol).
  2. PanglaoDB canonical markers, GI tract + immune/vasculature/connective
     organs (sensitivity/specificity quantified).
  3. CellTypist Human_Colorectal_Cancer model top-50 genes per cell type
     (logistic coefficients).

Each source's native cell labels are mapped onto six Major classes
(T/B/Mye/ILC/Epi/Stromal) via an explicit, auditable mapping table. Genes
voted by >=2 sources form the proxy set per class. Coverage against the
4000-gene spatial panel decides per-class viability (recorded, not forced).

Exploratory grade only. No expression data from this project is used here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

MAJORS = ("T", "B", "Mye", "ILC", "Epi", "Stromal")

# Deliberately unmapped: no matching Major exists in the GSE236581 reference
# (which has exactly these six classes), so mapping them anywhere would
# corrupt cross-route comparability. Recorded, not silently dropped.
DROPPED_BY_DESIGN = ("endothelial", "neuron", "adipocyte", "adipose",
                     "chondrocyte", "cholangiocyte", "cardiomyocyte",
                     "hepatocyte", "follicular cell")

# Native label (lowercased substring rules per source) -> Major.
CELLMARKER_RULES: list[tuple[str, str]] = [
    ("t cell", "T"), ("treg", "T"), ("tex", "T"), ("th17", "T"),
    ("t follicular", "T"), ("gamma delta", "T"), ("nkt", "T"),
    ("natural killer", "ILC"), ("cytokine induced killer", "T"),
    ("effector t", "T"), ("teff", "T"), ("tfr", "T"),
    ("eosinophil", "Mye"), ("basophil", "Mye"), ("granulocyte", "Mye"),
    ("endocrine cell", "Epi"), ("coloncyte", "Epi"),
    ("germinal center", "B"),
    ("b cell", "B"), ("plasma", "B"), ("b-1", "B"),
    ("macrophage", "Mye"), ("monocyte", "Mye"), ("dendritic", "Mye"),
    ("myeloid", "Mye"), ("mast cell", "Mye"), ("neutrophil", "Mye"),
    ("microglia", "Mye"), ("spp1", "Mye"), ("cDC", "Mye"),
    ("nk cell", "ILC"), ("ilc", "ILC"),
    ("epithelial", "Epi"), ("enterocyte", "Epi"), ("goblet", "Epi"),
    ("paneth", "Epi"), ("enteroendocrine", "Epi"), ("stem cell", "Epi"),
    ("progenitor", "Epi"), ("cancer", "Epi"), ("tumor cell", "Epi"),
    ("cms1", "Epi"), ("cms2", "Epi"), ("cms3", "Epi"), ("cms4", "Epi"),
    ("ta cell", "Epi"), ("tuft", "Epi"), ("m cell", "Epi"),
    ("fibroblast", "Stromal"), ("stromal", "Stromal"), ("myofibroblast", "Stromal"),
    ("smooth muscle", "Stromal"), ("pericyte", "Stromal"), ("glial", "Stromal"),
    ("mesenchymal", "Stromal"), ("caf", "Stromal"),
]

PANGLAO_RULES: list[tuple[str, str]] = [
    ("t cells", "T"), ("t cell", "T"), ("regulatory t", "T"), ("t-helper", "T"),
    ("cytotoxic", "T"), ("memory", "T"), ("gamma-delta", "T"),
    ("b cells", "B"), ("b cell", "B"), ("plasma", "B"), ("germinal", "B"),
    ("macrophages", "Mye"), ("monocytes", "Mye"), ("dendritic", "Mye"),
    ("myeloid", "Mye"), ("mast cells", "Mye"), ("neutrophils", "Mye"),
    ("microglia", "Mye"), ("basophils", "Mye"), ("eosinophils", "Mye"),
    ("nk cells", "ILC"), ("nk cell", "ILC"), ("ilc", "ILC"),
    ("enterocytes", "Epi"), ("goblet", "Epi"), ("paneth", "Epi"),
    ("enteroendocrine", "Epi"), ("epithelial", "Epi"), ("crypt", "Epi"),
    ("tuft", "Epi"), ("chief", "Epi"), ("parietal", "Epi"), ("foveolar", "Epi"),
    ("fibroblasts", "Stromal"), ("stromal", "Stromal"), ("myofibroblasts", "Stromal"),
    ("smooth muscle", "Stromal"), ("pericytes", "Stromal"), ("glia", "Stromal"),
    ("stellate", "Stromal"), ("mesenchymal", "Stromal"),
]

CELLTYPIST_RULES: list[tuple[str, str]] = [
    ("t cells", "T"), ("t cell", "T"), ("regulatory t", "T"), ("t helper", "T"),
    ("t follicular", "T"), ("gamma delta", "T"), ("proliferating", "T"),
    ("plasma", "B"), (" b", "B"), ("b,", "B"),
    ("macrophage", "Mye"), ("monocyte", "Mye"), ("cdc", "Mye"),
    ("myeloid", "Mye"), ("mast cells", "Mye"), ("neutrophil", "Mye"),
    ("pro-inflammatory", "Mye"), ("spp1", "Mye"),
    ("nk cells", "ILC"),
    ("enterocyte", "Epi"), ("goblet", "Epi"), ("stem-like", "Epi"),
    ("cms1", "Epi"), ("cms2", "Epi"), ("cms3", "Epi"), ("cms4", "Epi"),
    ("intermediate", "Epi"),
    ("myofibroblast", "Stromal"), ("smooth muscle", "Stromal"),
    ("pericyte", "Stromal"), ("stromal", "Stromal"),
    ("endothelial", "Stromal"), ("ec", "Stromal"),
    ("glial", "Stromal"),
]


def apply_rules(label: str, rules: list[tuple[str, str]]) -> str | None:
    text = f" {label.lower()} "
    for needle, major in rules:
        if needle in text:
            return major
    return None


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db-dir", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--symbol-map", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--min-votes", type=int, default=2)
    parser.add_argument("--min-panel-genes", type=int, default=10)
    args = parser.parse_args()

    db = args.db_dir
    per_source: dict[str, dict[str, set[str]]] = {}
    unmapped: dict[str, list[str]] = {}

    # 1. CellMarker2.0 CRC subset.
    cm = pd.read_parquet(db / "cellmarker_crc_subset.parquet")
    cm_sets: dict[str, set[str]] = {m: set() for m in MAJORS}
    cm_unmapped: list[str] = []
    for _, row in cm.iterrows():
        major = apply_rules(str(row["cell_name"]), CELLMARKER_RULES)
        sym = str(row["Symbol"]).strip().upper()
        if major is None:
            cm_unmapped.append(str(row["cell_name"]))
        elif sym and sym != "NAN":
            cm_sets[major].add(sym)
    per_source["cellmarker2.0"] = cm_sets
    unmapped["cellmarker2.0"] = sorted(set(cm_unmapped))

    # 2. PanglaoDB full table, filtered to CRC-relevant organs here
    # (GI + immune/vasculature/connective/blood/epithelium/liver).
    pg = pd.read_csv(db / "PanglaoDB_markers.tsv.gz", sep="\t")
    pg_sets: dict[str, set[str]] = {m: set() for m in MAJORS}
    pg_unmapped: list[str] = []
    keep_organs = {"GI tract", "Immune system", "Vasculature", "Connective tissue",
                   "Blood", "Epithelium", "Liver"}
    for _, row in pg.iterrows():
        if str(row["organ"]) not in keep_organs:
            continue
        major = apply_rules(str(row["cell type"]), PANGLAO_RULES)
        sym = str(row["official gene symbol"]).strip().upper()
        if major is None:
            pg_unmapped.append(str(row["cell type"]))
        elif sym and sym != "NAN":
            pg_sets[major].add(sym)
    per_source["panglaodb"] = pg_sets
    unmapped["panglaodb"] = sorted(set(pg_unmapped))

    # 3. CellTypist CRC top-50.
    with open(db / "celltypist_crc_top50.json") as f:
        ct = json.load(f)["markers"]
    ct_sets: dict[str, set[str]] = {m: set() for m in MAJORS}
    ct_unmapped: list[str] = []
    for cell_type, genes in ct.items():
        if cell_type.lower() in {"unknown", "intermediate"}:
            ct_unmapped.append(f"{cell_type} (skipped: non-informative)")
            continue
        major = apply_rules(cell_type, CELLTYPIST_RULES)
        if major is None:
            ct_unmapped.append(cell_type)
            continue
        for gene, _ in genes:
            ct_sets[major].add(str(gene).strip().upper())
    per_source["celltypist_crc"] = ct_sets
    unmapped["celltypist_crc"] = sorted(set(ct_unmapped))

    # Panel mapping ENSG -> symbol (upper).
    import csv
    ensg_to_sym: dict[str, str] = {}
    with open(args.symbol_map, newline="") as f:
        for row in csv.reader(f):
            if len(row) >= 2:
                ensg_to_sym.setdefault(row[0], row[1].upper())
    panel_ids = [l.strip() for l in Path(args.panel).read_text().splitlines() if l.strip()]
    panel_syms = {ensg_to_sym[g.replace("DEPRECATED_", "").split(".")[0]]
                  for g in panel_ids
                  if g.replace("DEPRECATED_", "").split(".")[0] in ensg_to_sym}

    # Vote + coverage.
    classes = {}
    for major in MAJORS:
        votes: dict[str, int] = {}
        sources: dict[str, list[str]] = {}
        for source, sets in per_source.items():
            for gene in sets[major]:
                votes[gene] = votes.get(gene, 0) + 1
                sources.setdefault(gene, []).append(source)
        voted = sorted(g for g, v in votes.items() if v >= args.min_votes)
        in_panel = sorted(g for g in voted if g in panel_syms)
        classes[major] = {
            "voted_genes": voted,
            "n_voted": len(voted),
            "in_panel": in_panel,
            "n_in_panel": len(in_panel),
            "viable": len(in_panel) >= args.min_panel_genes,
            "per_source_counts": {s: len(per_source[s][major]) for s in per_source},
            "gene_sources": {g: sources[g] for g in in_panel},
        }
    result = {
        "schema": "r04.marker_proxy.v1",
        "evidence_grade": "exploratory_post_hoc_no_confirmatory_claim",
        "sources": {
            "cellmarker2.0": {"file": "Cell_marker_Human.xlsx",
                              "sha256": sha256_file(db / "Cell_marker_Human.xlsx")},
            "panglaodb": {"file": "PanglaoDB_markers.tsv.gz",
                          "sha256": sha256_file(db / "PanglaoDB_markers.tsv.gz")},
            "celltypist_crc": {"file": "celltypist_crc_top50.json",
                               "sha256": sha256_file(db / "celltypist_crc_top50.json")},
        },
        "min_votes": args.min_votes,
        "min_panel_genes": args.min_panel_genes,
        "dropped_by_design_no_matching_major": list(DROPPED_BY_DESIGN),
        "unmapped_labels": unmapped,
        "classes": classes,
        "overall_viable": all(classes[m]["viable"] for m in MAJORS),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=1))
    for major in MAJORS:
        c = classes[major]
        print(f"{major:8s} voted={c['n_voted']:4d} in_panel={c['n_in_panel']:4d} "
              f"viable={c['viable']} src={c['per_source_counts']}")
    print("OVERALL_VIABLE:", result["overall_viable"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
