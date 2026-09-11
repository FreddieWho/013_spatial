#!/usr/bin/env python3
"""Build provisional WUSTL sample crosswalk + deep-3D inventory (R-03 Phase B).

Inputs (all public, no auth):
  author_mapping/ST_subclone_publication/Data_access/Sample_ID_Lookup_table_v1.xlsx
  author_mapping/3d-analysis/scripts/{visium,cosmx,xenium,he}_reference/*.csv
Outputs (tracked):
  WUSTL_sample_crosswalk.tsv, WUSTL_deep3d_inventory.tsv
Rules: unknowns stay empty; U-number section hints are E1 only; z_position /
section_order are always null here (see §18 of the task brief).
"""
from __future__ import annotations

import csv
import re
import shutil
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOOKUP = (HERE / "author_mapping/ST_subclone_publication/Data_access"
          / "Sample_ID_Lookup_table_v1.xlsx")
TRACK = HERE / "author_mapping/3d-analysis/scripts"
META = HERE / "metadata"
SHEET_NAMES = ["Readme", "Visium", "sc_snRNA_seq", "bulk RNA-seq",
               "CODEX", "WXS", "Xenium"]
PINS = {
    "ST_subclone_publication":
        "7067cd16f06ec4aa2500aa1e5c9a6eb1e6e42cfa",
    "3d-analysis": "14f04c1b787991898380305702c31d52c19fa694",
    "mushroom": "dc3598950fd233f2255a8e53006c915cc425a4a7",
}
DEEP3D = ["HT704B1", "HT891Z1", "HT913Z1", "HT206B1", "HT397B1"]
ASSAY_CANON = {
    "visium": "Visium", "visiumhd": "VisiumHD", "xenium": "Xenium",
    "cosmx": "CosMx", "codex": "CODEX", "he": "H&E",
    "lightsheet": "lightsheet", "3' snrna-seq": "3' snRNA-seq",
    "bulk rna-seq": "bulk RNA-seq", "wxs": "WXS",
}
SPATIAL_ASSAYS = {"Visium", "VisiumHD", "Xenium", "CosMx", "CODEX",
                  "H&E", "lightsheet"}


def canon_assay(raw: str) -> str:
    import html
    key = html.unescape(raw or "").strip().lower()
    return ASSAY_CANON.get(key, html.unescape(raw or "").strip())


PLATFORM_BY_ASSAY = {"Visium": "Visium", "CODEX": "CODEX", "Xenium": "Xenium"}

CROSSWALK_COLS = ["paper_sample_id", "participant_id", "biospecimen_id",
                  "parent_biospecimen_id", "specimen_id", "block_candidate",
                  "section_candidate", "assay", "platform", "source_file",
                  "source_commit", "evidence_grade"]


def read_xlsx_sheets(path: Path) -> dict[str, list[dict[str, str]]]:
    z = zipfile.ZipFile(path)
    strs = re.findall(r"<t[^>]*>(.*?)</t>",
                      z.read("xl/sharedStrings.xml").decode("utf-8"))
    out: dict[str, list[dict[str, str]]] = {}
    for idx, name in enumerate(SHEET_NAMES, start=1):
        try:
            xml = z.read(f"xl/worksheets/sheet{idx}.xml").decode("utf-8")
        except KeyError:
            continue
        rows: dict[int, dict[str, str]] = {}
        for m in re.finditer(r'<c r="([A-Z]+)(\d+)"([^>]*)>(.*?)</c>', xml):
            col, row, attrs, body = m.groups()
            sm = re.search(r't="s"><v>(\d+)</v>', "<c " + attrs + ">" + body)
            if sm:
                val = strs[int(sm.group(1))]
            else:
                vm = re.search(r"<v>([^<]*)</v>", body)
                val = vm.group(1) if vm else ""
            rows.setdefault(int(row), {})[col] = val
        ordered = [rows[r] for r in sorted(rows)]
        out[name] = ordered
    return out


def section_hint(wustl_id: str) -> tuple[str, str]:
    m = re.search(r"U(\d+)", wustl_id or "")
    if m:
        return m.group(1), "E1"
    return "", "E0"


def main() -> int:
    META.mkdir(parents=True, exist_ok=True)
    sheets = read_xlsx_sheets(LOOKUP)
    # copy tracking CSVs into metadata/
    copied = []
    for sub in ["visium_reference", "cosmx_reference", "xenium_reference",
                "he_reference"]:
        for csv_path in sorted((TRACK / sub).glob("*.csv")):
            dest = META / f"3d-analysis_{sub}_{csv_path.name}"
            shutil.copyfile(csv_path, dest)
            copied.append(dest.name)
    rows: list[dict[str, str]] = []

    def add(paper, participant, bio, assay, platform, source, commit, grade,
            section="", block=""):
        if not paper:
            return
        sec, gsec = (section, grade) if section else section_hint(paper)
        rows.append({
            "paper_sample_id": paper,
            "participant_id": participant,
            "biospecimen_id": bio,
            "parent_biospecimen_id": "",
            "specimen_id": "",
            "block_candidate": block,
            "section_candidate": sec,
            "assay": assay,
            "platform": platform,
            "source_file": source,
            "source_commit": commit,
            "evidence_grade": gsec,
        })

    lookup_src = ("ST_subclone_publication/Data_access/"
                  "Sample_ID_Lookup_table_v1.xlsx")
    for sheet, recs in sheets.items():
        if sheet == "Readme" or not recs:
            continue
        last_assay = ""
        for rec in recs[1:]:
            # merged Excel cells leave assay empty: forward-fill within sheet
            if rec.get("A", "").strip():
                last_assay = canon_assay(rec["A"])
            assay = last_assay
            # skip stray header echoes
            if rec.get("C", "") in ("", "WUSTL Biospecimen ID"):
                continue
            add(rec.get("C", ""), rec.get("D", ""), rec.get("E", ""),
                assay, PLATFORM_BY_ASSAY.get(assay, ""), lookup_src,
                PINS["ST_subclone_publication"], "E1")
    track_pin = PINS["3d-analysis"]
    level3 = {
        "visium": "Visium", "cosmx": "CosMx", "xenium": "Xenium",
    }
    for fname in sorted(META.glob("3d-analysis_*Level 3*.csv")):
        assay = canon_assay(
            next((a for k, a in level3.items() if k in fname.name.lower()), ""))
        with open(fname, newline="") as f:
            for rec in csv.DictReader(f):
                wustl = (rec.get("WUSTL Specimen") or "").strip()
                fid = (rec.get("HTAN_DATA_FILE_ID") or "").strip()
                bio = "" if fid.lower().startswith("matt fills") else fid
                plat = canon_assay(rec.get("PLATFORM") or assay)
                add(wustl, (rec.get("WUSTL Participant") or "").strip(), bio,
                    assay, plat if plat in SPATIAL_ASSAYS else "",
                    f"3d-analysis tracking: {fname.name}", track_pin, "E1")
    survey = META / ("3d-analysis_he_reference_HTAN2 3D Prostate Breast "
                     "Data Survey - Datafiles.csv")
    if survey.is_file():
        with open(survey, newline="") as f:
            for rec in csv.DictReader(f):
                exp = canon_assay(rec.get("experiment") or "")
                add((rec.get("Specimen name") or "").strip(),
                    (rec.get("participant name") or "").strip(), "", exp,
                    exp if exp in SPATIAL_ASSAYS else "",
                    f"3d-analysis survey: {survey.name}", track_pin, "E1")
    cosmx2d = META / "3d-analysis_cosmx_reference_2d data tracking - cosmx.csv"
    if cosmx2d.is_file():
        with open(cosmx2d, newline="") as f:
            for rec in csv.DictReader(f):
                sec_id = (rec.get("Section_ID") or "").strip()
                m = re.search(r"U([A-Za-z]*\d+)", sec_id)
                assay = canon_assay(rec.get("Assay_Type") or "CosMx")
                add(sec_id or (rec.get("Section_ID_clean") or "").strip(),
                    (rec.get("Participant_ID") or "").strip(), "",
                    assay, assay if assay in SPATIAL_ASSAYS else "",
                    f"3d-analysis tracking: {cosmx2d.name}", track_pin, "E1",
                    section=(m.group(1) if m else ""),
                    block=(rec.get("Block_ID") or "").strip())
    with open(HERE / "WUSTL_sample_crosswalk.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CROSSWALK_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    deep = [r for r in rows if r["participant_id"] in DEEP3D]
    with open(HERE / "WUSTL_deep3d_inventory.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CROSSWALK_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(deep)
    from collections import Counter
    print(f"crosswalk rows={len(rows)} assays={Counter(r['assay'] for r in rows)}")
    print(f"deep3d rows={len(deep)} participants={Counter(r['participant_id'] for r in deep)}")
    print(f"tracking CSVs copied={len(copied)}")
    prov = HERE / "provenance" / "crosswalk_provenance.json"
    prov.parent.mkdir(parents=True, exist_ok=True)
    import json, hashlib
    prov.write_text(json.dumps({
        "lookup_sha256": hashlib.sha256(LOOKUP.read_bytes()).hexdigest(),
        "pins": PINS,
        "tracked_csvs": copied,
        "section_hint_rule": "U<digits> in WUSTL ID -> section_candidate, E1 only",
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
