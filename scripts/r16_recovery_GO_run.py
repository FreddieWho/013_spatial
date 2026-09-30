#!/usr/bin/env python3
"""Corrected D-140 entrypoint. Outputs: infra/repair_20260921/.

The historical implementation and its outputs are retained as audit evidence.
Requires `python -m r16.recovery.repair_pipeline prepare` before analysis.
"""
import argparse
from r16.recovery import repair_pipeline as replay

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    replay.replay(go=True)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
