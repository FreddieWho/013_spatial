#!/usr/bin/env bash
# R-04 fold-4 GPU 主运行：预检（fail-closed hash 锚点）→ K=0/K=3 配对跑 → 端点审计
# 在实例上执行：nohup bash infra/r04/gpu_fold4_20260905/run_fold4.sh > /tmp/run_fold4.log 2>&1 &
set -euo pipefail
PROJ=/home/huyudi/013_spatial
cd "$PROJ"
export PYTHONPATH="$PROJ"
export R04_STORAGE_RESERVE_BYTES=100000000000   # 实例盘 <1.2TB，D-086 先例
export CUBLAS_WORKSPACE_CONFIG=:4096:8          # seed_everything 启用确定性算法，cublas 必须配
PY=/root/miniconda3/bin/python
MANIFEST=infra/r04/panel_cache_training/training_panel_manifest.json
GENES=infra/r04/gene_universe.txt
OUT=infra/r04/torch_representative_fold4_20260905_k03

# ---------- 预检：fold-0 2-step wiring，校验已知 hash 锚点（验证上传数据完整） ----------
PRE=/tmp/r04_gpu_preflight_$(date +%Y%m%d)
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

# ---------- 主运行：fold 4，K=0/K=3 同协议配对（脚本强制 in-run K=0 baseline） ----------
mkdir -p "$OUT"
(nvidia-smi --query-gpu=memory.used --format=csv,noheader -l 30 > "$OUT/vram.log" 2>/dev/null & echo $! > /tmp/vram.pid) || true
$PY scripts/r04_real_k_search.py --manifest-json "$MANIFEST" --gene-list "$GENES" \
  --output-dir "$OUT" --k-values 0,3 --folds 5 --gene-folds 2 --group-seed 20260807 \
  --split-seed 20260817 --restart-index 0 --fold-values 4 --steps 8400 --inference-steps 2400 \
  --checkpoint-steps 600 --inducing-points 16 --lengthscale 3.0 --learning-rate 0.01 \
  --gene-batch-size 512 --diagnostic-interval 100 --evaluation-interval 0 --evaluation-mc-draws 4 \
  --optimization-schedule joint --shared-steps 0
kill "$(cat /tmp/vram.pid)" 2>/dev/null || true

# ---------- 端点审计：对两个 K 的 fold-4 最终 checkpoint（重复性 + hash 核验） ----------
for K in 0 3; do
  $PY scripts/r04_checkpoint_audit.py --manifest-json "$MANIFEST" --gene-list "$GENES" \
    --output-json "infra/r04/torch_k${K}_fold4_checkpoint_audit_$(date +%Y%m%d).json" \
    --checkpoint "$OUT/checkpoints/k${K}/fold4/checkpoint.pt" \
    --factors "$K" --fold 4 --mc-draws 4
done

cp /tmp/run_fold4.log "$OUT/run.log" 2>/dev/null || true
echo RUN_ALL_DONE
