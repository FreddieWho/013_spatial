#!/usr/bin/env python3
"""Find auditable cross-source duplicate candidates from bundled metadata."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import tempfile
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit


MAX_METADATA_BYTES = 20 * 1024 * 1024
FIELDS = [
    "candidate_group_id",
    "left_record_id",
    "right_record_id",
    "relation_type",
    "evidence_grade",
    "resolution",
    "leakage_group_id",
    "canonical_source_url",
    "left_raw_patient",
    "left_raw_subseries",
    "notes",
]

# Audited non-singleton mappings. These are explicit, versioned interpretations
# of bundled region/run/subseries metadata, not expression-derived matches.
EXACT_MEMBER = {
    "TENX184": "Visium_HD_Human_Ovarian_Cancer_FF",
    "TENX185": "Visium_HD_Human_Ovarian_Cancer_FF_Min_Depth",
    "TENX168": "Visium_HD_Human_Lung_Cancer_HD_Only_Experiment1",
    "TENX169": "Visium_HD_Human_Lung_Cancer_post_Xenium_v1_Experiment1",
    "TENX170": "Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2",
    "TENX171": "Visium_HD_Human_Lung_Cancer_post_Xenium_Prime_5K_Experiment2",
    "TENX98": "Xenium_V1_FFPE_Human_Breast_IDC_Big_2",
    "TENX99": "Xenium_V1_FFPE_Human_Breast_IDC_Big_1",
    "TENX94": "Xenium_V1_FFPE_Human_Breast_ILC",
    "TENX95": "Xenium_V1_FFPE_Human_Breast_IDC",
    "TENX89": "CytAssist_FFPE_Human_Colon_Rep2",
    "TENX90": "CytAssist_FFPE_Human_Colon_Rep1",
    "TENX91": "CytAssist_FFPE_Human_Colon_Post_Xenium_Rep2",
    "TENX92": "CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1",
    "TENX105": "Xenium_V1_hKidney_cancer_section",
    "TENX189": "Xenium_Prime_Human_Lung_Cancer_FFPE",
    "TENX190": "Xenium_V1_Human_Lung_Cancer_FFPE",
    "TENX120": "Xenium_V1_hLiver_cancer_section_FFPE",
    "TENX111": "Xenium_V1_hColon_Cancer_Add_on_FFPE",
    "TENX182": "Visium_HD_3prime_Human_Ovarian_Cancer_FF",
    "TENX183": "Visium_HD_3prime_Human_Ovarian_Cancer_FF_Min_Depth",
    "TENX96": "Xenium_V1_FFPE_Human_Breast_ILC_With_Addon",
    "TENX97": "Xenium_V1_FFPE_Human_Breast_IDC_With_Addon__1.0.2",
    "TENX134": "Xenium_V1_hBoneMarrow_acute_lymphoid_leukemia_section",
}
POSSIBLE_MEMBER = {
    ("TENX184", "Visium_HD_Human_Ovarian_Cancer_FF_Min_Depth"),
    ("TENX185", "Visium_HD_Human_Ovarian_Cancer_FF"),
    ("TENX182", "Visium_HD_3prime_Human_Ovarian_Cancer_FF_Min_Depth"),
    ("TENX183", "Visium_HD_3prime_Human_Ovarian_Cancer_FF"),
    ("TENX168", "Visium_HD_Human_Lung_Cancer_post_Xenium_v1_Experiment1"),
    ("TENX169", "Visium_HD_Human_Lung_Cancer_HD_Only_Experiment1"),
    ("TENX170", "Visium_HD_Human_Lung_Cancer_post_Xenium_Prime_5K_Experiment2"),
    ("TENX171", "Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2"),
    ("TENX98", "Xenium_V1_FFPE_Human_Breast_IDC_Big_1"),
    ("TENX99", "Xenium_V1_FFPE_Human_Breast_IDC_Big_2"),
    ("TENX89", "CytAssist_FFPE_Human_Colon_Post_Xenium_Rep2"),
    ("TENX91", "CytAssist_FFPE_Human_Colon_Rep2"),
    ("TENX90", "CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1"),
    ("TENX92", "CytAssist_FFPE_Human_Colon_Rep1"),
    ("TENX111", "Xenium_V1_hColon_Cancer_Base_FFPE"),
}


def canonicalize_url(raw_url: Any) -> str:
    if raw_url is None:
        return ""
    value = str(raw_url).strip()
    if not value or value.lower() == "nan":
        return ""
    parsed = urlsplit(value if "://" in value else f"https://{value}")
    host = parsed.netloc.lower()
    if host.startswith("www."):
        host = host[4:]
    path = re.sub(r"/+", "/", parsed.path).rstrip("/")
    return f"{host}{path}"


def _raw(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value)


def _classify(rows: list[dict[str, str]]) -> None:
    group_sizes: dict[str, int] = {}
    for row in rows:
        group = row["candidate_group_id"]
        group_sizes[group] = group_sizes.get(group, 0) + 1

    adjacency: dict[str, set[str]] = {}
    for row in rows:
        hest = row["left_record_id"].removeprefix("HEST::")
        tenx = row["right_record_id"].removeprefix("TENX::")
        if group_sizes[row["candidate_group_id"]] == 1 or EXACT_MEMBER.get(hest) == tenx:
            row["relation_type"] = "same_source_dataset_lineage"
            row["evidence_grade"] = "E2_corroborated"
            row["resolution"] = "CONFIRMED_SOURCE_DATASET"
        elif (hest, tenx) in POSSIBLE_MEMBER:
            row["relation_type"] = "possible_physical_relation"
            row["evidence_grade"] = "E2_corroborated"
            row["resolution"] = "POSSIBLE_PHYSICAL_RELATION"
        else:
            row["relation_type"] = "same_source_page_only"
            row["evidence_grade"] = "E0_unknown"
            row["resolution"] = "NO_EVIDENCE"
        if row["resolution"] != "NO_EVIDENCE":
            left = row["left_record_id"]
            right = row["right_record_id"]
            adjacency.setdefault(left, set()).add(right)
            adjacency.setdefault(right, set()).add(left)

    component_by_node: dict[str, str] = {}
    unseen = set(adjacency)
    while unseen:
        seed = min(unseen)
        stack = [seed]
        component: set[str] = set()
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            stack.extend(adjacency.get(node, ()))
        unseen.difference_update(component)
        digest = hashlib.sha256("|".join(sorted(component)).encode("utf-8")).hexdigest()[:16]
        leakage = f"LEAKAGE::{digest}"
        for node in component:
            component_by_node[node] = leakage
    for row in rows:
        if row["resolution"] == "NO_EVIDENCE":
            row["leakage_group_id"] = ""
        else:
            row["leakage_group_id"] = component_by_node[row["left_record_id"]]


def _validate_small(path: Path, suffix: str) -> None:
    if path.suffix.lower() != suffix:
        raise ValueError(f"unexpected metadata suffix: {path}")
    if path.stat().st_size > MAX_METADATA_BYTES:
        raise ValueError(f"metadata file too large: {path}")


def find_hest_tenx_candidates(
    hest_metadata_dir: Path,
    tenx_manifest: Path,
) -> list[dict[str, str]]:
    _validate_small(tenx_manifest, ".tsv")
    with tenx_manifest.open(newline="", encoding="utf-8-sig") as handle:
        tenx_rows = list(csv.DictReader(handle, delimiter="\t"))

    tenx_by_url: dict[str, list[dict[str, str]]] = {}
    for row in tenx_rows:
        canonical = canonicalize_url(row.get("page_url"))
        dataset_dir = (row.get("dataset_dir") or "").strip()
        if canonical and dataset_dir:
            tenx_by_url.setdefault(canonical, []).append(row)

    candidates: list[dict[str, str]] = []
    for metadata_path in sorted(hest_metadata_dir.glob("*.json")):
        _validate_small(metadata_path, ".json")
        with metadata_path.open(encoding="utf-8") as handle:
            metadata = json.load(handle)
        if not isinstance(metadata, dict):
            raise ValueError(f"expected JSON object: {metadata_path}")
        canonical = canonicalize_url(metadata.get("download_page_link1"))
        if not canonical or canonical not in tenx_by_url:
            continue
        hest_id = str(metadata.get("id") or metadata_path.stem).strip()
        group_hash = hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:16]
        for tenx_row in tenx_by_url[canonical]:
            dataset_dir = tenx_row["dataset_dir"].strip()
            candidates.append(
                {
                    "candidate_group_id": f"URL::{group_hash}",
                    "left_record_id": f"HEST::{hest_id}",
                    "right_record_id": f"TENX::{dataset_dir}",
                    "relation_type": "same_source_page",
                    "evidence_grade": "",
                    "resolution": "",
                    "leakage_group_id": "",
                    "canonical_source_url": canonical,
                    "left_raw_patient": _raw(metadata.get("patient")),
                    "left_raw_subseries": _raw(metadata.get("subseries")),
                    "notes": (
                        "Exact bundled source-page match; physical sample identity "
                        "still requires sample-level metadata"
                    ),
                }
            )

    _classify(candidates)
    candidates.sort(
        key=lambda row: (
            row["candidate_group_id"],
            row["left_record_id"],
            row["right_record_id"],
        )
    )
    return candidates


def write_candidates(rows: list[dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = ""
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            newline="",
            encoding="utf-8",
            dir=output_path.parent,
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = handle.name
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(
                sorted(
                    rows,
                    key=lambda row: (
                        row["candidate_group_id"],
                        row["left_record_id"],
                        row["right_record_id"],
                    ),
                )
            )
        os.replace(temporary, output_path)
    finally:
        if temporary:
            Path(temporary).unlink(missing_ok=True)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--hest-metadata",
        type=Path,
        default=root / "data/other_sources/hest1k/metadata",
    )
    parser.add_argument(
        "--tenx-manifest",
        type=Path,
        default=root / "data/other_sources/10x_genomics/manifest_3yr.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "infra/sample-registry/staging/"
            "hest_tenx_duplicate_candidates.tsv"
        ),
    )
    args = parser.parse_args()
    output = args.output.resolve()
    try:
        output.relative_to(root)
    except ValueError as exc:
        raise ValueError("output escapes project root") from exc

    rows = find_hest_tenx_candidates(args.hest_metadata, args.tenx_manifest)
    write_candidates(rows, output)
    groups = len({row["candidate_group_id"] for row in rows})
    counts: dict[str, int] = {}
    for row in rows:
        resolution = row["resolution"]
        counts[resolution] = counts.get(resolution, 0) + 1
    summary = " ".join(f"{key}={counts[key]}" for key in sorted(counts))
    print(
        f"cross-source candidates={len(rows)} source_page_groups={groups} "
        f"{summary}"
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
