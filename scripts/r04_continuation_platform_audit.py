#!/usr/bin/env python3
"""CLI for the read-only R-04 continuation platform audit."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from r04.continuation_platform_audit import audit_continuation_manifest
from r04.runtime import atomic_json

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--base-dir", type=Path, default=None)
    args = parser.parse_args()
    manifest = json.loads(args.manifest_json.read_text(encoding="utf-8"))
    result = audit_continuation_manifest(manifest, args.base_dir or args.manifest_json.parent)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(args.output_json, result)
    return 0 if result["status"] == "ALL_PASS" else 1
if __name__ == "__main__":
    raise SystemExit(main())
