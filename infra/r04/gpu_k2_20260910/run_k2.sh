#!/usr/bin/env bash
# R-04 K=2 Stage A GPU 主运行（D-104，L-001）：预检 → fold-0 K=0/K=2 配对 → fold-4 K=0/K=2 配对 → 端点审计
# 在实例上执行：nohup bash infra/r04/gpu_k2_20260910/run_k2.sh > /tmp/run_k2.log 2>&1 &
set -euo pipefail
PROJ=/home/huyudi/013_spatial
cd "$PROJ"
export PYTHONPATH="$PROJ"
export R04_STORAGE_RESERVE_BYTES=100000000000
export CUBLAS_WORKSPACE_CONFIG=:4096:8
PY=/root/miniconda3/bin/python
MANIFEST=infra/r04/panel_cache_training/training_panel_manifest.json
GENES=infra/r04/gene_universe.txt
STAMP=$(date +%Y%m%d)

# ---------- 预检：fold-0 2-step wiring，校验已知 hash 锚点 ----------
PRE=/tmp/r04_gpu_preflight_${STAMP}
rm -rf "$PRE"
$PY scripts/r04_real_k_search.py --manifest-json "$MANIFEST" --gene-list "$GENES" \
  --output-dir "$PRE" --k-values 0 --folds 5 --gene-folds 2 --group-seed 20260807 \
  --split-seed 20260817 --restart-index 0 --fold-values 0 --steps 2 --inference-steps 2 \
  --checkpoint-steps 600 --inducing-points 16 --lengthscale 3.0 --learning-rate 0.01 \
  --gene-batch-size 512 --diagnostic-interval 100 --evaluation-interval 0 --evaluation-mc-draws 4 \
  --optimization-schedule joint --shared-steps 0
$PY - "$PRE" <<'PY'
import json, sys
from pathlib import Path
pre = Path(sys.argv[1])
top = json.load(open(pre / "r04_real_k_search.json"))
cell = json.load(open(pre / "cells/k0/fold0/cell.json"))
oih = cell.get("fit_diagnostics", {}).get("objective_input_hash", "")
ok = top["status"] == "WIRING_ONLY_NOT_SCIENTIFIC"
ok &= top["manifest_hash"].startswith("5a0dc1de")
ok &= cell["input_hash"].startswith("044d28c5")
ok &= oih.startswith("b5a3fb68")
print("PREFLIGHT", "PASS" if ok else "FAIL", cell["input_hash"][:16], oih[:16], top["manifest_hash"][:16])
sys.exit(0 if ok else 3)
PY
rm -rf "$PRE"

# ---------- 主运行：两 fold 各做 K=0/K=2 同协议配对 ----------
for FOLD in 0 4; do
  OUT=infra/r04/torch_representative_fold${FOLD}_${STAMP}_k02
  mkdir -p "$OUT"
  (nvidia-smi --query-gpu=memory.used --format=csv,noheader -l 30 > "$OUT/vram.log" 2>/dev/null & echo $! > /tmp/vram.pid) || true
  $PY scripts/r04_real_k_search.py --manifest-json "$MANIFEST" --gene-list "$GENES" \
    --output-dir "$OUT" --k-values 0,2 --folds 5 --gene-folds 2 --group-seed 20260807 \
    --split-seed 20260817 --restart-index 0 --fold-values "$FOLD" --steps 8400 --inference-steps 2400 \
    --checkpoint-steps 600 --inducing-points 16 --lengthscale 3.0 --learning-rate 0.01 \
    --gene-batch-size 512 --diagnostic-interval 100 --evaluation-interval 0 --evaluation-mc-draws 4 \
    --optimization-schedule joint --shared-steps 0
  kill "$(cat /tmp/vram.pid)" 2>/dev/null || true
  cp /tmp/run_k2.log "$OUT/run.log" 2>/dev/null || true
  for K in 0 2; do
    $PY scripts/r04_checkpoint_audit.py --manifest-json "$MANIFEST" --gene-list "$GENES" \
      --output-json "infra/r04/torch_k${K}_fold${FOLD}_checkpoint_audit_${STAMP}.json" \
      --checkpoint "$OUT/checkpoints/k${K}/fold${FOLD}/checkpoint.pt" \
      --factors "$K" --fold "$FOLD" --mc-draws 4
  done
done

echo RUN_ALL_DONE
