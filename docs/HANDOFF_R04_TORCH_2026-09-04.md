# R-04 Torch 工作交接

更新时间：2026-09-04（Asia/Shanghai）  
适用范围：PyTorch 迁移后的 R-04 真实数据代表性复跑，以及进入空间场稳健性分析前的计算收尾。

这是一份交接文档，不是新的科学方案。科学目标、判据和路线仍以 `docs/plan.md`、`docs/roadmap.md`、`docs/decisions.md`、`docs/ISSUES.md` 和 `STATUS.md` 为准。

## 先看结论

- 当前本机没有正在运行的 R-04 计算进程。
- fold 0 的 K=0 拟合 checkpoint 已达到 `global_step=8400`；最终 checkpoint 位于：
  `infra/r04/torch_representative_fold0_20260902_cont_k0_7800/checkpoints/k0/fold0/`。
- 随后的 K=0、2400-step 留出推断/评分没有完成。上次进程超过 7200 秒限制后退出，输出目录只有复制的 checkpoint，没有 `cell.json`，不能作为评分结果或聚合输入。
- 当前没有新的 Torch 科学结论，`selected_k` 仍为 `null`。不要把遗留的顶层 `status=RUNNING` 文件当成仍在运行的进程。
- 下一真正决策点是计算资源：CPU 继续运行会非常慢；若批准 GPU，先解决 V100 与当前 CUDA 环境的兼容性，再完成 K=0 留出推断和 K=3 fold 0。

## 1. 科学边界

模型表示的是空间分子效应：

\[
h(\mu_{ig})=o_i+b_g+n_{ig}+\sum_{k=1}^{K_{model}}\lambda_{gk}f_k(s_i).
\]

- `K_model` 是允许的潜在空间效应维度，也就是表示容量；`K_eff` 是由留出预测、restart/fold 稳定性、空间置换 null 和重复性共同支持的有效维度。
- K 不是组织结构、结构实例或“真实生物场”的数量。一个结构可以使用多个潜在方向，多个结构也可以共享一个方向。
- factor 编号、正负号和旋转不唯一。后续主要看效应子空间、重建的空间效应和独立 GT 下的结构特异读出，不给单个 factor 直接命名为 TLS、血管或坏死。
- 两个来源完全共线时，应报告联合场或不可识别；不能通过增大 K 强行拆分。
- 代表性 fold 0 运行即使成功，也只是 Torch 链路和参数协议的诊断，不能单独宣布 K 或证明生物假设。

## 2. 当前状态矩阵

| 项目 | 状态 | 可接受解释 |
| --- | --- | --- |
| R-01 / R-02 | 已完成的上游数据、角色和 GT 约束 | 可用于 R-04；不重新定义结构标签 |
| 旧 TensorFlow 五 seed 面板 | 已封存 | 历史稳定性诊断，不与新 Torch 数值混合 |
| Torch 2-step K=0/K=3 fold 0 预检 | 已完成 | `wiring-only`，只证明接口可通 |
| Torch K=0 fold 0 fit | checkpoint 到 8400 | 拟合诊断，不是留出科学结果 |
| Torch K=0 fold 0 留出 inference/scoring | 未完成 | 上次 CPU 运行超时；禁止聚合 |
| Torch K=3 fold 0 fit/inference/scoring | 未开始 | 完成 K=0 后用完全相同协议运行 |
| fold 4 K=0/K=3 | 未开始 | fold 0 结果和资源可接受后再运行 |
| `selected_k` | `null` | 保留；不因当前缺少结果而填值 |
| 空间置换 null、组成拆分、独立 lineage 复现 | 未开始 | R-04 代表性计算和结构读出门禁之后进行 |
| GPU | 未启用 | `AGENTS.md` 当前禁止未经批准启用 |

注意：当前 `STATUS.md` 的“正在补最后 600 步”描述落后于实际 artifact。接手者先以本文件和 checkpoint 文件为准，并在下一次状态更新时修正 `STATUS.md`，不要修改历史 roadmap 记录。

## 3. 权威输入与结果位置

### 必读顺序

1. `AGENTS.md`
2. `STATUS.md`
3. `docs/plan.md`
4. `docs/roadmap.md`
5. `docs/decisions.md`
6. `docs/ISSUES.md`
7. `infra/r04/README.md`
8. `scripts/r04_real_k_search.py`
9. `r04/models/mnsf.py`

### 数据与协议输入

- training manifest：`infra/r04/panel_cache_training/training_panel_manifest.json`
- gene universe：`infra/r04/gene_universe.txt`
- training panel：47 个切片、30 个患者、4,000 个固定基因；总规模约 113,344 spots。
- Torch 环境记录：`infra/r04/environment.lock.json`
- 虚拟环境：`/home/huyudi/.venvs/013_spatial-r04-torch`

### 当前 Torch 输出

- 2-step wiring-only 预检：`/tmp/r04_torch_real_preflight_20260902/`（临时诊断，不是科学结果）。
- K=0 首段 checkpoint（step 4200）：`infra/r04/torch_representative_fold0_20260902/`。
- K=0 continuation（step 7800）：`infra/r04/torch_representative_fold0_20260902_cont_k0/`。
- K=0 最终 continuation（step 8400）：`infra/r04/torch_representative_fold0_20260902_cont_k0_7800/`。
  - `cell.json`：`FIT_ONLY_DIAGNOSTIC`；只包含最后 600 步的 fit 诊断，未做 inference。
  - checkpoint：Torch v3、`params`/`best_state`/Adam state 和输入 hash 均保留。
- K=0 评分失败/未完成目录：`infra/r04/torch_representative_fold0_20260902_k0_score/`。
  - 顶层 JSON 留有 `status=RUNNING`，但没有 `cell.json`，也没有活动进程；这是超时后的遗留标记。
  - 不得把该目录作为 `FIT_AND_SCORED` 结果、不得加入聚合。

### 历史但仍需保留的证据

- `infra/r04/restart_stability_panel_20260831.json`
- `infra/r04/cross_fold_k0_k3_diagnostic_20260825.json`
- `infra/r04/continuation_platform_audit_20260831_v2/`
- `remote_outputs/restart_cont_r{2,3,4}_frozen_infer2400_20260831_v3/`

历史冻结推断 12 个单元已完成，但严格训练审计仍为 `STOPPED_OFF_PLATFORM`，实际接受范围只是冻结 held-out inference/scoring；它们不能代替新的 Torch 匹配重跑。

## 4. 固定计算协议

除非追加 decision 并说明不可替代性，不要改变以下参数：

| 参数 | 值 |
| --- | --- |
| outer folds / representative fold | 5 / fold 0（必要时 fold 4） |
| gene cross-fit | 2；`split_seed=20260817` |
| group seed | `20260807` |
| K | K=0 公平基线、K=3；K=2 仅条件触发 |
| fit steps | 8400 |
| held-out inference steps | 2400；每个 gene split 各一次 |
| gene batch | 512 |
| inducing points / lengthscale | 16 / 3.0 |
| learning rate | 0.01 |
| objective | NB + GP；`nonspatial_rank=1` |
| schedule | `joint`，`shared_steps=0` |
| diagnostic / evaluation interval | 100 / 0 |
| evaluation MC draws | 4 |
| checkpoint interval | 600 |

代表性运行使用 `--fold-values 0`，因此脚本必须继续输出 `wiring_only=true`；这不是完整 outer-fold 科学估计。

## 5. 下一步执行顺序

### A. 资源与环境预检

1. 未经用户批准不要启用 GPU。
2. 若使用 V100：当前环境记录是 CUDA 13 构建，而 Volta/V100 不应直接使用该构建；准备支持 Volta 的 CUDA 12.6 legacy PyTorch 环境，并重新记录环境 hash。
3. 先跑 2-step 和 50–100-step GPU 预检，确认 `torch.cuda.is_available()`、设备名、显存、输入/objective hash 和 checkpoint schema，再决定是否继续长跑。
4. 若无法使用 GPU，不要把超长 CPU 运行伪装成已启动后台任务；留下命令、资源限制和 `BLOCKED_CPU_RUNTIME` 状态。

### B. 完成当前代表性链路

1. 用 step-8400 K=0 checkpoint 完成两路 2400-step held-out inference/scoring；必须生成完整 `cell.json`，并审计输入、配置、环境和 checkpoint hash。
2. 用同一 manifest、fold、gene split、NB/GP、优化器和训练步数运行 K=3 fold 0 fit，再运行两路 held-out inference/scoring。
3. 比较 K=0/K=3 的 fit objective platform、inference platform、患者级配对分数和方向；仍标记为代表性/诊断结果，不填 `selected_k`。
4. 若 fold 0 结果和资源时间可接受，再考虑 fold 4；不要自动扩展到 K>3。

### C. 科学门禁之后

- 先做 training-role GT 的共享成分、结构特异残余和 K 敏感性读出。
- 只有当主要结构结论依赖不稳定第三方向时才运行 K=2 Stage A。
- 结构结论 K-robust 且空间协议冻结后，才进入空间置换 null、组成—状态拆分和独立 lineage 复现。
- 不能用单 factor 命名代替子空间读出，也不能用 K 未唯一确定作为无限期阻塞理由。

## 6. 资源估计与风险

- 单张 V100 32GB 对当前单 cell、K≤3 和 gene batch=512 的显存容量预计足够；不要并行塞多个 cell。
- 主机内存建议 48GB；32GB 是当前代码路径的实际下限。代码会把 panel counts 转成 dense float32，并且推断阶段使用完整 4,000-gene 矩阵，内存余量不足容易在数据转换或临时张量阶段失败。
- 当前 CPU 观测是：K=0 fit 8400 步约数小时，2400-step 双 gene split inference 超过 7200 秒仍未产出 cell。GPU 没有实测速度，因此不能承诺固定倍数。
- 租用时长建议：只做 fold 0 的 K=0/K=3，预留 24 小时；若要连 fold 4，预留 48 小时。空间 null 等重复实验另行按短预检结果估算。
- 项目必须始终保留至少 1.2TB 存储余量；不删除旧 checkpoint、历史审计或 remote output 来腾空间。

## 7. 可复现入口

激活固定环境：

```bash
source /home/huyudi/.venvs/013_spatial-r04-torch/bin/activate
```

代表性 K=3 fold 0 的入口（输出目录必须新建，不能覆盖历史目录）：

```bash
python scripts/r04_real_k_search.py \
  --manifest-json infra/r04/panel_cache_training/training_panel_manifest.json \
  --gene-list infra/r04/gene_universe.txt \
  --output-dir infra/r04/torch_representative_fold0_20260904_k3 \
  --k-values 3 \
  --folds 5 --gene-folds 2 \
  --group-seed 20260807 --split-seed 20260817 \
  --restart-index 0 --fold-values 0 \
  --steps 8400 --inference-steps 2400 \
  --checkpoint-steps 600 \
  --inducing-points 16 --lengthscale 3.0 \
  --learning-rate 0.01 --gene-batch-size 512 \
  --diagnostic-interval 100 --evaluation-interval 0 \
  --evaluation-mc-draws 4 \
  --optimization-schedule joint --shared-steps 0
```

该命令是 fold-0 wiring-only/代表性诊断入口，不是正式完整 K 搜索。K=0 的已有 checkpoint 不应被该命令从头覆盖；先明确采用可审计的 checkpoint-only inference 入口，再运行评分。

## 8. 不得做的事

- 不把旧 TensorFlow 数值和新 Torch 数值混成一个 panel。
- 不把 `FIT_ONLY_DIAGNOSTIC`、`RUNNING` 遗留标记或未生成 `cell.json` 的目录称为完成。
- 不用 fold 0 单个代表性结果宣布 K、结构数量或生物机制。
- 不改变 gene universe、fold、gene split、objective、阈值或 optimizer 来追求更快通过。
- 不因 GPU 尚未批准而静默租用；也不在 V100 上直接复用未经验证的 CUDA 13 环境。
- 不提前运行空间 null、组成拆分、结构命名或独立 lineage 复现来绕过当前代表性链路。

## 9. 交接验收

接手者完成本阶段前，至少应能指出：

1. K=0 checkpoint 是否达到 8400，且 checkpoint 审计是否通过；
2. K=0/K=3 每个 cell 是否都有完整 fit、inference、scoring 和 provenance；
3. 哪些结果是 wiring-only、哪些是历史冻结推断、哪些可以进入科学聚合；
4. `selected_k` 是否仍为 `null`，以及为什么不能把 K 当结构数；
5. GPU 环境是否通过兼容性和短预检；
6. 在空间 null、组成拆分和独立复现之前，不得扩大 H-01/H-03/H-04/H-05 的科学结论。

