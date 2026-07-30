#!/usr/bin/env python3
"""Inventory allowlisted, small, bundled metadata for roadmap node R-01.

The script rejects expression matrices, images, generic archives, HDF5 files,
oversized inputs and paths outside the project root. Configured ``.xlsx``
metadata workbooks may be hashed as OOXML bytes, but are not parsed here.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any, Iterable


MAX_METADATA_BYTES = 20 * 1024 * 1024
MAX_COLLECTION_BYTES = 100 * 1024 * 1024
ALLOWED_SUFFIXES = {
    "table": {".csv", ".tsv"},
    "file": {".md", ".txt"},
    "metadata_workbook": {".xlsx"},
    "json_collection": {".json"},
    "text_collection": {".md", ".txt"},
}

OUTPUT_FIELDS = [
    "source_id",
    "kind",
    "path_pattern",
    "member_count",
    "row_count",
    "bytes_total",
    "fields",
    "identity_fields",
    "section_fields",
    "gt_risk_fields",
    "sha256",
    "status",
    "notes",
]

IDENTITY_TOKENS = {
    "accession",
    "block",
    "case",
    "donor",
    "id",
    "patient",
    "sample",
    "section",
    "source",
    "study",
    "subject",
}
SECTION_TOKENS = {
    "z",
    "z_step_size",
    "section",
    "section_order",
    "slice",
    "z_position",
    "zposition",
    "spacing",
    "thickness",
    "order",
}
GT_RISK_TOKENS = {
    "annotation",
    "class",
    "cluster",
    "count",
    "gt",
    "label",
    "location",
    "mask",
    "maturation",
    "pathology",
    "presence",
    "region",
    "score",
    "segmentation",
    "tls",
}


def _safe_path(root: Path, relative: str) -> Path:
    root = root.resolve()
    path = (root / relative).resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"configured path escapes project root: {relative}") from exc
    return path


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_member(path: Path, kind: str) -> None:
    suffix = path.suffix.lower()
    allowed = ALLOWED_SUFFIXES[kind]
    if suffix not in allowed:
        raise ValueError(
            f"disallowed suffix for {kind}: {path.name} ({suffix or 'none'})"
        )
    size = path.stat().st_size
    if size > MAX_METADATA_BYTES:
        raise ValueError(
            f"metadata file too large: {path} ({size} > {MAX_METADATA_BYTES})"
        )


def _collection_sha256(root: Path, members: Iterable[Path]) -> str:
    digest = hashlib.sha256()
    for member in members:
        digest.update(str(member.relative_to(root)).encode("utf-8"))
        digest.update(b"\0")
        digest.update(_sha256_file(member).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def _matches(fields: set[str], tokens: set[str]) -> str:
    matched = []
    for field in fields:
        normalized = field.lower().replace("-", "_").strip()
        parts = set(normalized.split("_"))
        if normalized in tokens or parts.intersection(tokens):
            matched.append(field)
    return ";".join(sorted(matched))


def _base_record(source: dict[str, Any]) -> dict[str, Any]:
    return {
        "source_id": source["source_id"],
        "kind": source["kind"],
        "path_pattern": source.get("path", source.get("glob", "")),
        "member_count": 0,
        "row_count": 0,
        "bytes_total": 0,
        "fields": "",
        "identity_fields": "",
        "section_fields": "",
        "gt_risk_fields": "",
        "sha256": "",
        "status": "MISSING",
        "notes": source.get("notes", ""),
    }


def _finish_record(
    record: dict[str, Any],
    fields: set[str],
    members: list[Path],
    root: Path,
) -> dict[str, Any]:
    record["member_count"] = len(members)
    record["bytes_total"] = sum(path.stat().st_size for path in members)
    record["fields"] = ";".join(sorted(fields))
    record["identity_fields"] = _matches(fields, IDENTITY_TOKENS)
    record["section_fields"] = _matches(fields, SECTION_TOKENS)
    record["gt_risk_fields"] = _matches(fields, GT_RISK_TOKENS)
    record["sha256"] = _collection_sha256(root, members)
    record["status"] = "PRESENT"
    return record


def inspect_source(root: Path, source: dict[str, Any]) -> dict[str, Any]:
    root = root.resolve()
    record = _base_record(source)
    kind = source["kind"]

    if kind in {"table", "file", "metadata_workbook"}:
        path = _safe_path(root, source["path"])
        if not path.is_file():
            return record
        members = [path]
    elif kind in {"json_collection", "text_collection"}:
        pattern = source["glob"]
        if Path(pattern).is_absolute() or ".." in Path(pattern).parts:
            raise ValueError(f"unsafe glob: {pattern}")
        members = sorted(path.resolve() for path in root.glob(pattern) if path.is_file())
        if not members:
            return record
        for member in members:
            member.relative_to(root)
    else:
        raise ValueError(f"unsupported source kind: {kind}")

    for member in members:
        _validate_member(member, kind)
    total_bytes = sum(member.stat().st_size for member in members)
    if total_bytes > MAX_COLLECTION_BYTES:
        raise ValueError(
            "metadata collection too large: "
            f"{source['source_id']} ({total_bytes} > {MAX_COLLECTION_BYTES})"
        )

    fields: set[str] = set()
    if kind == "table":
        delimiter = "\t" if members[0].suffix.lower() == ".tsv" else ","
        with members[0].open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            fields.update(reader.fieldnames or [])
            record["row_count"] = sum(1 for _ in reader)
    elif kind == "json_collection":
        for member in members:
            with member.open(encoding="utf-8") as handle:
                payload = json.load(handle)
            if not isinstance(payload, dict):
                raise ValueError(f"expected JSON object: {member}")
            fields.update(str(key) for key in payload)
        record["row_count"] = len(members)
    elif kind == "text_collection":
        record["row_count"] = sum(
            sum(1 for line in member.open(encoding="utf-8", errors="replace") if line.strip())
            for member in members
        )
    elif kind in {"file", "metadata_workbook"}:
        record["row_count"] = 1

    return _finish_record(record, fields, members, root)


def run_inventory(
    root: Path,
    config_path: Path,
    output_path: Path,
) -> list[dict[str, Any]]:
    root = root.resolve()
    config_path = config_path.resolve()
    try:
        config_path.relative_to(root)
    except ValueError as exc:
        raise ValueError("catalog escapes project root") from exc
    output_path = output_path.resolve()
    try:
        output_path.relative_to(root)
    except ValueError as exc:
        raise ValueError("output escapes project root") from exc
    if output_path == config_path:
        raise ValueError("output collides with catalog input")

    with config_path.open(encoding="utf-8") as handle:
        sources = json.load(handle)
    if not isinstance(sources, list):
        raise ValueError("catalog must contain a JSON list")

    source_ids = [source.get("source_id") for source in sources]
    if any(not source_id for source_id in source_ids):
        raise ValueError("every source requires source_id")
    if len(source_ids) != len(set(source_ids)):
        raise ValueError("source_id values must be unique")

    for source in sources:
        if "path" in source:
            configured = _safe_path(root, source["path"])
            if configured == output_path:
                raise ValueError("output collides with metadata input")
        elif "glob" in source:
            pattern = source["glob"]
            if Path(pattern).is_absolute() or ".." in Path(pattern).parts:
                raise ValueError(f"unsafe glob: {pattern}")
            if any(path.resolve() == output_path for path in root.glob(pattern)):
                raise ValueError("output collides with metadata input")

    records = [inspect_source(root, source) for source in sources]
    records.sort(key=lambda row: row["source_id"])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_name = ""
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
            temp_name = handle.name
            writer = csv.DictWriter(
                handle,
                fieldnames=OUTPUT_FIELDS,
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(records)
        os.replace(temp_name, output_path)
    finally:
        if temp_name:
            Path(temp_name).unlink(missing_ok=True)
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="project root",
    )
    parser.add_argument(
        "--catalog",
        type=Path,
        default=Path("infra/sample-registry/source_catalog.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("infra/sample-registry/source_metadata_inventory.tsv"),
    )
    args = parser.parse_args()

    root = args.root.resolve()
    catalog = args.catalog if args.catalog.is_absolute() else root / args.catalog
    output = args.output if args.output.is_absolute() else root / args.output
    records = run_inventory(root, catalog, output)
    present = sum(record["status"] == "PRESENT" for record in records)
    print(f"metadata sources: {present}/{len(records)} present")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
