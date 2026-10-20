#!/usr/bin/env python3
"""Combine identity, data-audit and registration checks into one annotatable table.

Read-only. Produces the section-level gate that the pilot choice would sit on.
No pilot is locked here and no label is created.
"""

from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path("/home/huyudi/013_spatial")
HEST = ROOT / "data/other_sources/hest1k"
ELIG = ROOT / "infra/tls_four_cancer_labels_20261008/eligibility"
ORGANS = {"Bowel", "Lung", "Breast", "Liver"}

# Gates, fixed before reading the combined table.
MIN_INSIDE = 0.90          # fraction of in-tissue spots on the tissue polygon
MIN_PANEL_SETS = 100       # usable sets of the frozen 1742-set panel
CELL_LEVEL_UM = 0.6        # at or below this is cell-level texture
# A human can still judge aggregate architecture and density well above 0.6
# um/px; 1.0 um/px was set after rendering 0.9 um/px fields for review.
HUMAN_MAX_UM = 1.0


def specimen_key(source_id: str, subseries: str) -> str:
    sub = subseries or ""
    stem = sub.split("_")[0]
    m = re.match(r"([A-Za-z]*\d+)", stem)
    return f"{source_id}::{m.group(1) if m else stem or 'unnamed'}"


def main() -> None:
    ids = {}
    for name, key in (("sections.tsv", None),
                      ("data_audit_sections.tsv", None),
                      ("registration/registration_check.tsv", None)):
        for row in csv.DictReader((ELIG / name).open(), delimiter="\t"):
            ids.setdefault(row["section_id"], {}).update(row)

    rows = []
    for sid, rec in sorted(ids.items()):
        meta = json.loads((HEST / "metadata" / f"{sid}.json").read_text())
        inside = rec.get("identity_inside")
        panel_ok = int(rec.get("panel_valid") or 0) >= MIN_PANEL_SETS
        px_um = float(rec["px_um_authoritative"]) if rec.get("px_um_authoritative") else None
        reg_ok = inside is not None and float(inside) >= MIN_INSIDE
        coarse = px_um is None or px_um > CELL_LEVEL_UM
        # Healthy tissue is background only. It is never an annotatable TLS target
        # and never a negative.
        disease = rec.get("disease_state", "")
        usable_disease = disease in {"Cancer", "Tumor"}
        if not reg_ok:
            status = "blocked_registration"
        elif not panel_ok:
            status = "blocked_panel"
        elif coarse:
            status = "annotatable_coarse_only"
        elif usable_disease:
            status = "annotatable_candidate"
        else:
            status = "background_only_not_negative"

        # Second gate, for a human reading H&E. The 0.6 um/px bar exists because
        # HookNet needs 40x texture; it does not limit what a pathologist can
        # read, and a targeted expression panel is irrelevant to morphology.
        # The human gate asks a different question: does the image register, is
        # it legible, and is it a plausible TLS context.
        if not reg_ok:
            hstatus = "blocked_registration"
        elif px_um is not None and px_um <= HUMAN_MAX_UM:
            hstatus = ("annotatable_candidate" if usable_disease
                       else "background_only_not_negative")
        else:
            hstatus = "annotatable_coarse_only"
        if hstatus == "annotatable_coarse_only" and px_um is not None:
            hreason = f"resolution {px_um:.3f} um/px above the human bar {HUMAN_MAX_UM}"
        elif not panel_ok:
            hreason = "targeted panel is irrelevant to morphology; kept for human reading"
        elif usable_disease:
            hreason = "registers, legible, tumour context"
        else:
            hreason = f"{disease} context: background only, never a negative"
        row["gate_status_human"] = hstatus
        row["gate_reason_human"] = hreason
        rows.append({
            "section_id": sid,
            "organ": rec.get("organ"),
            "platform": rec.get("platform"),
            "disease_state": disease,
            "specimen_key": specimen_key(rec.get("source_id", ""), str(meta.get("subseries") or "")),
            "patient_status": rec.get("patient_status"),
            "patient_key": rec.get("patient_key"),
            "identity_note": rec.get("identity_note"),
            "px_um": f"{px_um:.4f}" if px_um else "",
            "spot_inside_tissue": inside or "",
            "registration_verdict": rec.get("verdict"),
            "panel_valid_sets": rec.get("panel_valid"),
            "panel_usable": "yes" if panel_ok else "no",
            "tenx_possible_same_release": rec.get("tenx_possible_same_release"),
            "gate_status": status,
            "gate_status_human": hstatus,
            "gate_reason_human": hreason,
            "pilot_selected": "no",
            "labels_present": "no",
        })

    with (ELIG / "annotatable_gate.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    cands = [r for r in rows if r["gate_status"] == "annotatable_candidate"]
    coarse = [r for r in rows if r["gate_status"] == "annotatable_coarse_only"]
    units = defaultdict(list)
    for r in cands:
        units[r["specimen_key"]].append(r["section_id"])
    known_units = {k: v for k, v in units.items()
                   if any(x["patient_status"] == "known_within_study"
                          for x in rows if x["section_id"] in v)}

    summary = {
        "gates": {"min_inside": MIN_INSIDE, "min_panel_sets": MIN_PANEL_SETS,
                  "cell_level_um": CELL_LEVEL_UM,
                  "human_max_um": HUMAN_MAX_UM},
        "status_counts": dict(Counter(r["gate_status"] for r in rows)),
        "status_counts_human": dict(Counter(r["gate_status_human"] for r in rows)),
        "by_organ": {o: dict(Counter(r["gate_status"] for r in rows if r["organ"] == o))
                     for o in sorted(ORGANS)},
        "annotatable_candidates": len(cands),
        "distinct_specimens_among_candidates": len(units),
        "specimens_with_known_patient": len(known_units),
        "known_patient_specimens": {k: v for k, v in sorted(known_units.items())},
        "unknown_identity_candidates": sorted(
            r["section_id"] for r in cands if r["patient_status"] != "known_within_study"),
        "coarse_only": sorted(r["section_id"] for r in coarse),
        "pilot_locked": False,
        "external_test_available": False,
        "labels_created": False,
    }
    (ELIG / "annotatable_gate_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()