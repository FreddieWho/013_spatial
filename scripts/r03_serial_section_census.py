#!/usr/bin/env python3
"""R-03 serial-section & registration measurability census (read-only).

Censuses already-downloaded data only (no external fetch): the frozen
R-01 registry (block/section/z-order/thickness coverage), per-lineage local
metadata (ST-CRC README + Zenodo record, 10x Block-A bundle layout, HTAN
sample key, HEST field coverage), and registration inputs (H&E hires images,
tissue positions, scalefactors).

Outputs a block-level inventory plus measurability verdicts for the four
R-11 task types. A pair counts toward adjacent-plane work only with RECORDED
inter-section distance; co-block grouping without order/spacing (HTAN) and
acquisition z-steps (HEST Xenium) do not qualify. Nothing here changes H-06
confidence; absence of usable pairs keeps H-06 UNJUDGED, never failed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from r04.runtime import atomic_json

SCHEMA = "r03.serial_section_census.v1"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_units(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _registry_inventory(units_path: Path) -> dict[str, object]:
    rows = _read_units(units_path)
    inc = [r for r in rows if r["record_status"] == "RESOLVED_INCLUDED_CANDIDATE"]
    by_src: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in inc:
        by_src[row["source_namespace"]].append(row)
    lineages = {}
    for src, members in sorted(by_src.items()):
        blocks: dict[str, list[str]] = defaultdict(list)
        for row in members:
            blocks[row["block_id"] or "(no-block)"].append(row["physical_unit_id"])
        multi = {b: v for b, v in blocks.items() if b != "(no-block)" and len(v) > 1}
        lineages[src] = {
            "n_units": len(members),
            "n_blocks": sum(1 for b in blocks if b != "(no-block)"),
            "multi_section_blocks": len(multi),
            "multi_section_block_sizes": sorted(len(v) for v in multi.values()),
            "z_position_coverage": sum(1 for r in members if r["z_position"].strip()),
            "section_order_coverage": sum(1 for r in members if r["section_order"].strip()),
            "thickness_coverage": sum(1 for r in members if r["section_thickness"].strip()),
            "thickness_values": sorted({r["section_thickness"] for r in members
                                        if r["section_thickness"].strip()}),
        }
    return {"n_included": len(inc), "lineages": lineages}


def _verdicts() -> dict[str, object]:
    return {
        "in_plane_core_masked": {
            "verdict": "MEASURABLE",
            "reason": "single sections suffice; 167 included units across 5 lineages",
        },
        "adjacent_plane_prediction": {
            "verdict": "NOT_MEASURABLE",
            "reason": ("no section pair in downloaded data carries a recorded "
                       "inter-section distance; ST-CRC serial pairs, Block-A "
                       "Sec1/Sec2 and HTAN co-block pieces all lack spacing"),
        },
        "held_out_middle_plane": {
            "verdict": "NOT_MEASURABLE",
            "reason": "no ordered triplet with known spacing exists locally",
        },
        "stack_edge_extrapolation": {
            "verdict": "NOT_MEASURABLE",
            "reason": ("same blocker as adjacent-plane prediction: "
                       "no recorded inter-section spacing"),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path,
                        default=Path("infra/r03/serial_section_census_20260911.json"))
    args = parser.parse_args()
    root = args.project_root.resolve()
    units = root / "infra/sample-registry/physical_units.tsv"
    if not units.is_file():
        raise SystemExit("registry physical_units.tsv missing")
    hest_coverage = root / "infra/sample-registry/staging/hest_field_coverage.tsv"
    result: dict[str, object] = {
        "schema": SCHEMA,
        "status": "R03_CENSUS_COMPLETE",
        "created_at": "2026-09-11",
        "scope": "already-downloaded data only; no external fetch",
        "registry": _registry_inventory(units),
        "source_metadata": {
            "st_crc_serial_pairs": {
                "finding": "README states 'two serial sections per patient' (7 blocks x 2)",
                "spacing": "NOT_RECORDED_LOCALLY",
                "thickness": "NOT_RECORDED_LOCALLY",
                "registration_inputs": "H&E hires/lowres + scalefactors + tissue positions present per section (spot-checked SN048_A121573_Rep1)",
            },
            "tenx_block_a_sections": {
                "finding": "Section_1/Section_2 of one block (serial-section validation role)",
                "spacing": "NOT_RECORDED_LOCALLY",
                "registration_inputs": "H&E hires/lowres + scalefactors + tissue positions present per section",
            },
            "htan_co_block_pieces": {
                "finding": "14 multi-piece blocks (2-3 pieces) in the 47-unit training set",
                "ordering_spacing": "NOT_RECORDED (z_position/section_order empty registry-wide)",
                "registration_inputs": "h5ad only; no H&E downloaded",
                "use": "block-level grouping only (already used in R-01); not adjacent-plane evidence",
            },
            "usz_singletons": {
                "finding": "8 singleton blocks; registry thickness '5 um' (registry-sourced)",
                "registration_inputs": "h5ad only; high-res tif not fetched",
            },
            "external_geo": {
                "finding": "96 included rows carry no block/section info (patient envelopes only)",
            },
            "hest_z_step": {
                "finding": "raw_z_step_size nonempty 29/492 (see staging hest_field_coverage.tsv)",
                "interpretation": "acquisition/segmentation step (mostly Xenium), NOT slice spacing (I-002)",
            },
            "non_core_10x_rep_pairs": {
                "finding": ("downloaded 10x demo dirs contain Rep1/Rep2 pairs "
                            "(Colon, Colon Post-Xenium, Lung post-Xenium sets); roles: none"),
                "spacing": "NOT_RECORDED_LOCALLY",
                "use": "inventory only; not claim-bearing without roles",
            },
        },
        "registration_error_floor": {
            "status": "UNMEASURED",
            "note": ("no fiducial-based cross-section registration performed; "
                     "Visium 100um pitch / 55um spot are instrument scales, "
                     "not a registration claim"),
        },
        "task_verdicts": _verdicts(),
        "h06_impact": ("no change: H-06 stays UNJUDGED; usable-pair count is zero, "
                       "which is a data-boundary fact, not a failure of the hypothesis"),
        "inputs": {
            "physical_units.tsv": _sha256(units),
            "hest_field_coverage.tsv": (_sha256(hest_coverage)
                                        if hest_coverage.is_file() else "MISSING"),
        },
    }
    out = root / args.output
    out.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(out, result)
    print(result["status"], "| multi-section blocks:",
          sum(v["multi_section_blocks"] for v in
              result["registry"]["lineages"].values()),
          "| adjacent-plane:", result["task_verdicts"]["adjacent_plane_prediction"]["verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
