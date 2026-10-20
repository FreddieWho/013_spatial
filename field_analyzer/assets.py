"""Asset inventory for FA-00.

Produces the five-state availability table required by
FA_CONTRACT_ADDENDUM_20261020.md §6 and the per-section pilot check.

State vocabulary:
  VERIFIED_LOCAL               present and usable as declared
  PRESENT_BUT_INSUFFICIENT     present, but fewer/partial units than the contract
  NOT_PRESENT                  path absent
  REFERENCED_LOCAL_NOT_VERIFIED not checked here (kept for bookkeeping only)
  TRACKED_SOURCE_ONLY          source code in git, no data payload
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

from .adapters import scan_h5ad_header

CONFIRMED_STATUS = "known_within_study"


@dataclass
class UnitTable:
    units: dict[str, list[str]]
    rows: list[dict[str, str]]

    @property
    def n_units(self) -> int:
        return len(self.units)

    @property
    def n_sections(self) -> int:
        return sum(len(v) for v in self.units.values())


def read_identity(path: Path) -> UnitTable:
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    units: dict[str, list[str]] = {}
    for row in rows:
        if row["patient_confirmed_within_study"] != "True":
            continue
        units.setdefault(row["unit_id"], []).append(row["section_id"])
    for key in units:
        units[key] = sorted(units[key])
    return UnitTable(units=units, rows=rows)


def assign_folds(units: Sequence[str], n_folds: int = 5) -> dict[str, int]:
    """Deterministic cross-unit fold assignment (sorted unit id, round robin)."""
    ordered = sorted(units)
    return {unit: i % n_folds for i, unit in enumerate(ordered)}


def unit_role_rows(units: dict[str, list[str]], n_folds: int = 5) -> list[dict[str, Any]]:
    folds = assign_folds(list(units), n_folds=n_folds)
    rows: list[dict[str, Any]] = []
    for unit in sorted(units):
        for section in sorted(units[unit]):
            rows.append(
                {
                    "unit_id": unit,
                    "section_id": section,
                    "fold": folds[unit],
                    "role": "grouped_holdout",
                    "is_replicate_of": "",
                    "note": "same unit => same fold; mask/rotation variants inherit the unit role",
                }
            )
    return rows


def _stat_row(path: Path) -> tuple[bool, int, float]:
    if not path.exists():
        return False, 0, 0.0
    stat = path.stat()
    return True, int(stat.st_size), float(stat.st_mtime)


def _manifest_version(path: Path) -> str:
    for candidate in ("prepare_version.json", "freeze_manifest.json", "contract.json"):
        f = path / candidate if path.is_dir() else path.parent / candidate
        if f.exists():
            try:
                data = json.loads(f.read_text())
            except Exception:
                continue
            for key in ("version", "run_id", "prepared_at", "freeze_id", "contract_version"):
                if key in data:
                    return f"{candidate}:{key}={data[key]}"
    return ""


ASSETS: Sequence[dict[str, str]] = (
    {
        "path": "data/other_sources/hest1k/st/{section_id}.h5ad",
        "what": "raw Visium counts + spatial coordinates (HEST 1k)",
        "reuse_decision": "ADAPT",
        "provenance": "R06; read via field_analyzer.adapters with explicit axis",
    },
    {
        "path": "infra/field_discovery_review_20261018/identity_review.tsv",
        "what": "corrected within-study patient identity",
        "reuse_decision": "REUSE",
        "provenance": "R08; unknown stays unknown",
    },
    {
        "path": "infra/field_source_holdout_20261018",
        "what": "previous source holdout splits and cached predictions",
        "reuse_decision": "REBUILD",
        "provenance": "R16; full-channel predictions are not block-recovery predictions",
    },
    {
        "path": "infra/tls_multiscale_gobp_20261002/panel_controls.json",
        "what": "51 fixed gene-set definitions (4 proxy fields among them)",
        "reuse_decision": "ADAPT",
        "provenance": "R07; A/B split and witness exclusion applied on top",
    },
    {
        "path": "scripts/field_discovery_20261018.py",
        "what": "legacy reader / feature views / neighbor graphs",
        "reuse_decision": "ADAPT",
        "provenance": "R06; silent transpose and self-fallback removed",
    },
    {
        "path": "scripts/field_candidate_subfield_20261019.py",
        "what": "legacy candidate subfield probe (A/B cache)",
        "reuse_decision": "REBUILD",
        "provenance": "R07; whole-cache load forbidden, role registry replaces it",
    },
    {
        "path": "scripts/field_biological_review_20261019.py",
        "what": "legacy field descriptors (centroid, patch stats, band autocorr)",
        "reuse_decision": "ADAPT",
        "provenance": "R08; residual scale may not be reported as field scale",
    },
    {
        "path": "scripts/tls_multiscale_gobp.py",
        "what": "legacy tissue graph and label-derived endpoints",
        "reuse_decision": "SEALED",
        "provenance": "R09; human-label endpoints barred from discovery",
    },
    {
        "path": "scripts/tls_multiscale_forward.py",
        "what": "legacy corridor path and sideband operators",
        "reuse_decision": "HISTORICAL_ONLY",
        "provenance": "R10; label-defined endpoints, symmetric sidebands",
    },
    {
        "path": "infra/r16/t2p2_20260928",
        "what": "legacy variogram matched null",
        "reuse_decision": "HISTORICAL_ONLY",
        "provenance": "R11; p-values and weak-sill lengths not reused",
    },
    {
        "path": "docs/TLS_MULTISCALE_BRIDGE_CALIBRATION_20261003.md",
        "what": "historical bridge calibration failure",
        "reuse_decision": "HISTORICAL_ONLY",
        "provenance": "R12; negative evidence, not evidence of absence",
    },
    {
        "path": "docs/TLS_PATIENT_VALIDATION_20261003.md",
        "what": "historical GSE patient metadata",
        "reuse_decision": "HISTORICAL_ONLY",
        "provenance": "R13; no labels or scores enter selection",
    },
)


def asset_rows(repo_root: Path, sections: Sequence[str] | None = None) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    known_sections = list(sections or [])
    for spec in ASSETS:
        raw = spec["path"]
        if "{section_id}" in raw:
            present = [s for s in known_sections if (repo_root / raw.format(section_id=s)).exists()]
            n_expected = len(known_sections)
            n_present = len(present)
            if n_expected == 0:
                state = "REFERENCED_LOCAL_NOT_VERIFIED"
            elif n_present == n_expected:
                state = "VERIFIED_LOCAL"
            elif n_present > 0:
                state = "PRESENT_BUT_INSUFFICIENT"
            else:
                state = "NOT_PRESENT"
            total_bytes = sum(
                (repo_root / raw.format(section_id=s)).stat().st_size for s in present
            )
            rows.append(
                {
                    "path": raw,
                    "what": spec["what"],
                    "availability": state,
                    "exists": bool(n_present),
                    "bytes": total_bytes,
                    "mtime": "",
                    "manifest_version": "identity_review.tsv",
                    "reuse_decision": spec["reuse_decision"],
                    "provenance": spec["provenance"],
                    "notes": f"{n_present}/{n_expected} confirmed sections present"
                    + (f"; missing={sorted(set(known_sections) - set(present))}" if n_present < n_expected else ""),
                }
            )
            continue
        path = repo_root / raw
        exists, size, mtime = _stat_row(path)
        state = "VERIFIED_LOCAL" if exists else "NOT_PRESENT"
        rows.append(
            {
                "path": raw,
                "what": spec["what"],
                "availability": state,
                "exists": exists,
                "bytes": size,
                "mtime": mtime if exists else "",
                "manifest_version": _manifest_version(path) if exists else "",
                "reuse_decision": spec["reuse_decision"],
                "provenance": spec["provenance"],
                "notes": "" if exists else "path absent on this machine",
            }
        )
    return rows


def section_rows(repo_root: Path, sections: Sequence[str]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for section in sorted(sections):
        path = repo_root / "data/other_sources/hest1k/st" / f"{section}.h5ad"
        info = scan_h5ad_header(path)
        rows.append(
            {
                "section_id": section,
                "path": str(path.relative_to(repo_root)) if path.exists() else str(path),
                "exists": info.get("exists", False),
                "bytes": info.get("bytes", 0),
                "encoding": info.get("encoding", ""),
                "shape": "x".join(str(v) for v in info.get("shape", ())),
                "n_barcodes": info.get("n_barcodes", 0),
                "n_genes": info.get("n_genes", 0),
                "has_spatial": info.get("has_spatial", False),
                "spatial_shape": "x".join(str(v) for v in info.get("spatial_shape", ())),
                "axis_consistent": info.get("declared_axis_consistent", False),
                "obs_columns": ",".join(info.get("obs_columns", [])[:8]),
            }
        )
    return rows


def write_tsv(rows: Iterable[dict[str, Any]], path: Path) -> Path:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return path
    keys = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in keys})
    return path
