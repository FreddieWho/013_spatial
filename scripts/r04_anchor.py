#!/usr/bin/env python3
"""Run post-freeze anchor evaluation; never used by discovery."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from r04.anchor_eval import evaluate_anchor_correspondence
from r04.runtime import atomic_json
from r04.serialization import read_json


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-registry", type=Path, required=True)
    parser.add_argument("--expected-hash", required=True)
    parser.add_argument("--field", type=Path, required=True)
    parser.add_argument("--anchor", type=Path, required=True)
    parser.add_argument("--groups", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    registry = read_json(args.candidate_registry)
    field = np.load(args.field)["values"]
    anchor = np.load(args.anchor)["values"]
    groups = np.load(args.groups)["values"]
    result = evaluate_anchor_correspondence(
        candidate_registry_hash=registry["candidate_registry_hash"],
        expected_registry_hash=args.expected_hash,
        field_values=field, anchor_values=anchor, groups=groups,
    )
    atomic_json(args.output, {"schema": "r04.anchor.v1", "status": result.status, "reasons": list(result.reasons), "metrics": dict(result.metrics)})
    return 0 if result.status == "ANCHOR_EVALUATED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
