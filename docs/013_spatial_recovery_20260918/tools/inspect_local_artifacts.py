#!/usr/bin/env python3
"""Read-only, metadata-only inventory for the 013_spatial recovery tasks.

No network, model fitting, expression-matrix loading, deletion, or Git writes.
A missing expected path does not prove that the underlying data do not exist.
"""
from __future__ import annotations
import argparse
import ast
import csv
import hashlib
import json
import shutil
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

REVIEW_COMMIT = '73996528634ffe0362f1a0341771658171f3b5e2'
SOURCE_PATHS = [
    'AGENTS.md', 'STATUS.md', 'TODO.md', 'LEADS.md',
    'docs/plan.md', 'docs/roadmap.md', 'docs/decisions.md',
    'docs/ISSUES.md', 'docs/R16_DISCOVERY_DESIGN_20260914.md',
    'scripts/r16_axis_factory.py', 'scripts/r16_laneA_screen.py',
    'scripts/r16_laneB_screen.py', 'scripts/r16_build_registry.py',
    'r16/census.py', 'r16/axes.py', 'r16/axis_factory.py', 'r16/section_io.py',
]
ARTIFACT_PATTERNS = [
    'infra/r16/field_registry.tsv', 'infra/r16/axis_factory*.npz',
    'infra/r16/axis_factory*.json', 'infra/r16/laneA*.json',
    'infra/r16/laneB*.json', 'infra/r16/joint_embed*.json',
    'infra/r16/axis_replication*.json', 'infra/r16/*manifest*.json',
    'infra/r16/census_cache_hvg10k/*manifest*.json',
    'infra/r16/census_cache_hvg10k/sections/*manifest*.json',
]


def git_read(root: Path, *args: str) -> dict:
    try:
        r = subprocess.run(['git', '-C', str(root), *args], check=False,
                           capture_output=True, text=True, timeout=15)
        return {'returncode': r.returncode, 'stdout': r.stdout.strip(),
                'stderr': r.stderr.strip()}
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {'error': f'{type(exc).__name__}: {exc}'}


def source_info(root: Path, relative: str) -> dict:
    p = root / relative
    if not p.is_file():
        return {'path': relative, 'state': 'NOT_AT_EXPECTED_PATH'}
    b = p.read_bytes()
    out = {'path': relative, 'state': 'PRESENT', 'size_bytes': len(b),
           'sha256': hashlib.sha256(b).hexdigest()}
    if p.suffix == '.py':
        try:
            tree = ast.parse(b.decode('utf-8'))
            out['function_locations'] = [
                {'name': n.name, 'start_line': n.lineno, 'end_line': n.end_lineno}
                for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
        except (SyntaxError, UnicodeError) as exc:
            out['parse_error'] = str(exc)
    return out


def registry_info(path: Path) -> dict:
    if not path.is_file():
        return {'state': 'NOT_AT_EXPECTED_PATH'}
    counts = {'kind': Counter(), 'tier': Counter(), 'evidence_grade': Counter()}
    rows = 0
    try:
        with path.open(encoding='utf-8') as f:
            lines = (line for line in f if line.strip() and not line.startswith('#'))
            for r in csv.DictReader(lines, delimiter='\t'):
                rows += 1
                for col in counts:
                    counts[col][r.get(col) or 'MISSING'] += 1
        return {'state': 'READ', 'row_count': rows,
                'counts': {k: dict(v) for k, v in counts.items()},
                'warning': 'Rows are heterogeneous records, not independent field discoveries.'}
    except (OSError, UnicodeError, csv.Error) as exc:
        return {'state': 'ERROR', 'error': f'{type(exc).__name__}: {exc}'}


def inventory(root: Path) -> dict:
    if not (root / 'docs/plan.md').is_file():
        raise ValueError('Expected repository docs/plan.md is missing. Supply the repository root.')
    head = git_read(root, 'rev-parse', 'HEAD')
    head_sha = head.get('stdout') if head.get('returncode') == 0 else None
    free = shutil.disk_usage(root).free
    artifacts, seen = [], set()
    for pattern in ARTIFACT_PATTERNS:
        paths = sorted(root.glob(pattern))
        if not paths:
            artifacts.append({'pattern': pattern, 'state': 'NO_MATCH_AT_EXPECTED_LOCATION'})
        for path in paths:
            if not path.is_file() or path in seen:
                continue
            seen.add(path)
            artifacts.append({'path': str(path.relative_to(root)), 'state': 'PRESENT',
                              'size_bytes': path.stat().st_size,
                              'contents_loaded': False})
    cache = root / 'infra/r16/census_cache_hvg10k/sections'
    cache_meta_count = (sum(1 for p in cache.glob('*.json') if 'manifest' not in p.name)
                        if cache.is_dir() else None)
    return {
        'schema': '013_spatial.recovery_preflight.v1',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'repo_root': str(root), 'review_commit': REVIEW_COMMIT,
        'current_head': head_sha, 'head_matches_review': head_sha == REVIEW_COMMIT,
        'head_diagnostic': head,
        'tracked_worktree_status': git_read(root, 'status', '--porcelain', '--untracked-files=no'),
        'disk_free_bytes': free, 'analysis_reserve_required_bytes': 1_200_000_000_000,
        'analysis_reserve_met': free >= 1_200_000_000_000,
        'disk_note': 'This small metadata read is allowed; before heavy analysis enforce the repository reserve.',
        'sources': [source_info(root, p) for p in SOURCE_PATHS],
        'artifacts': artifacts, 'cache_section_metadata_count': cache_meta_count,
        'registry': registry_info(root / 'infra/r16/field_registry.tsv'),
        'scope': 'METADATA_ONLY_NO_PATIENT_ANALYSIS',
        'counts_normalization_status': 'REQUIRES_MANIFEST_AND_BUILDER_REVIEW',
        'notes': ['Do not overwrite results if HEAD has advanced; inspect the actual diff.',
                  'Missing default paths are not evidence of absence of data.',
                  'No expression matrices, NPZ contents, or hidden labels were loaded.'],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--overwrite', action='store_true')
    args = ap.parse_args()
    root = args.repo.resolve()
    out = args.out.resolve()
    if out.exists() and not args.overwrite:
        ap.error('Output already exists; use a new path or explicitly pass --overwrite.')
    try:
        data = inventory(root)
    except (ValueError, OSError) as exc:
        ap.error(str(exc))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(f'Wrote metadata inventory: {out}')
    print(f'Current HEAD: {data["current_head"] or "unavailable"}')
    print(f'Registry records: {data["registry"].get("row_count", "unavailable")}')
    print('This is NOT a scientific result and does not certify data normalization.')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
