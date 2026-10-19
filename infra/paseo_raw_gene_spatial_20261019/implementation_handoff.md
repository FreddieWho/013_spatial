# 同原始基因输入空间对照：第一阶段实现交接（R-01 修复后）

任务书：`infra/paseo_raw_gene_spatial_20261019/implementation_task.md`。科学合同（主代理）：`docs/FIELD_RAW_GENE_SPATIAL_CONTRACT_20261019.md`（含 D-269、D-270）；合同预审：`infra/paseo_raw_gene_spatial_20261019/review/preflight.md`；首轮正式审核（NOT_ACCEPTED，唯一结果阻断 R-01）：`infra/paseo_raw_gene_spatial_20261019/review/reviewed_implementation_round1.json`。
状态：**R-01 已修复并按 D-270 复核口径重跑定向检查与小样本贯通；全量未运行；本版本尚未被独立审核**。本文件不代表审核通过，也不宣布任何科学结论。

## 0. 本轮勘误（如实记录，含被审核否定的措辞）

1. **等价检查容差确实被放宽过**：原机器界限 `1e-6`，实测复制块残差 6.18e-6 → 旧界限下检查失败；随后改为 `max(2×求解器地板, 1e-5)`。首轮审核承认数值依据（同设计换求解器差 3.40e-6、正则化 Gram 条件数 ≈5.2e5、错缩放对照 2.11e-4 ≈34× 残差），但**否定“这不是放宽”的说法**。现按 D-270 固定为 **1e-5**、不再自适应扩大；求解器地板与负对照只作记录。上版交接与本节的旧措辞以本节为准。
2. **邻域块写成 “WZ/VZ” 漏了 QC**：实际代码（自始）在图乘上作用的是 **X = [Z, QC]**，即 **W[Z,QC]、V[Z,QC]**；旧文字有误，以代码为准，现已改口径并在检查中显式验证（`candidate_neighbour_blocks_are_graph_aggregates_of_x` = (0.0, 0.0)）。
3. **全量模型行数是 2500 而不是 5000**：本轮只有 full15 一个视图 → 125 次切片访问 × 4 目标 × 5 模型 = **2500** 行；旧成本估算里的 5000 是两视图口径，评价时间估计相应减半。

## 1. R-01 修复内容（设计结构）

- **X = 冻结可见基因 + QC**；**E = 既有器官/平台指示列**（复用冻结 point 的 E）。
- 同点设计 = **`[X, E]`**。
- 真实邻域设计 = **`[X, WX, VX]` 各乘 `1/√3` 之后，再拼接 E 一次**；背景/置换对应替换两个邻域块（`[X, BG, BG, E]`、`[X, PWX, PVX, E]`）。
- **E 在块缩放之外、只出现一次、不进入图乘**；E 与同点模型使用完全相同的训练期标准化与惩罚（同一 `StandardScaler` 口径、同一 `alpha`、同截距）。
- 实现口径：`fit_linear(..., scaled_width)` 只对前 `scaled_width` 列乘 `1/√3`（训练行常量列剔除后按保留列重算该宽度）；`prepare_scaled`/`model_predict` 与拟合严格同口径。
- 邻域块仍是 **W[Z,QC]、V[Z,QC]**；背景是全片 X 的均值（非 512 训练行）。

## 2. 版本与 Git

- 本次提交：`scripts/field_raw_gene_spatial_20261019.py`、`scripts/run_field_raw_gene_spatial_20261019.sh`、`infra/field_raw_gene_spatial_20261019/{contract.json,input_checks.json,prepare_receipt.json,splits.json}`、`.../checks/{checks.json,solver_diagnostics.tsv}`、`.../archive/*.log`（含本轮新失败记录）、`.../smoke/`（小样本小产物）、本文件。
- 运行与失败日志（`stages.log`、`prepare.log`、`check.log`、`smoke.log`）默认被 `.gitignore` 的 `*.log` 排除，按“保留失败记录”要求继续用 `git add -f` 强制入库。
- 未提交（工作区，按项目惯例）：`checks/tiny/*.npz`、`smoke/splits/*.npz`、`smoke/evaluation/*`、`evaluation/`、`PREPARED`/`COMPLETE` 标记；逐文件清单见 `smoke/freeze_manifest.json` 与 `input_checks.json`。
- 我未修改主代理的科学合同或根文档；D-270 由主代理提交。

## 3. 参数与输入（不变部分）

`alpha=10`、截距拟合、`Ridge(solver='auto')`（稠密 n>p → cholesky；`tol=1e-6`、`max_iter=5000` 备查）、训练抽样 = 冻结 `sample_indices(..., seed=20261018, 512)`、权重 = 冻结单位等权/片等权、标准化仅训练行、训练行常量列剔除只按训练行决定；平滑 = `0.5*(point + neighbor_mean(near, point))`；置换 = 每片一次固定种子、同时作用于 near/outer、实际算子沿用同一 self fallback、不改 X/E/条码/目标顺序。
输入：每目标只读冻结 full15 可见矩阵 + QC；**原始 HDF5 不在允许表内**；划分用冻结 `splits.json`（不重抽）——内部 5 折恰好划分 72 片、来源 6 折恰好划分 53 片、训练/测试单位互斥、来源臂留出来源互斥；内部 65 单位 / 来源 47 单位；确认 14 片 / 8 单位；来源不足 19 片保留内部；乳腺/肝脏保留器官外推。

## 4. 本轮验证（定向检查 49 项 / 0 失败 + 短贯通）

命令：`bash scripts/run_field_raw_gene_spatial_20261019.sh check`（exit=0）、`... smoke`（exit=0）。

**R-01 直接证据（新增）**
| 检查 | 观测 |
|---|---|
| `block_and_context_column_layout` | (点设计宽, 邻域设计宽, E宽) = (1649, 4935, 6)，缩放宽度 0 与 3×1643 ✓ |
| `context_columns_same_standardisation_and_penalty_as_point` | 同点与邻域模型标准化后的 E 列最大差 **0.0** |
| `x_block_scaled_by_one_over_sqrt3_only` | X 块经 1/√3 后的最大差 **0.0** |
| `context_block_appended_once_unscaled` | 缩放宽度 (4929, 4929)、E 宽 6、E 与 `E` 原始列差 **0.0** |
| `candidate_neighbour_blocks_are_graph_aggregates_of_x` | (W 块−W[Z,QC], V 块−V[Z,QC]) = **(0.0, 0.0)** |
| `triple_copy_block_equivalence_fixed_tolerance` | 残差 **5.77e-6 ≤ 固定 1e-5**；求解器地板 3.46e-6（仅记录）；**错缩放对照 1.49e-4 ≈26×** |
| `permutation_keeps_x_block` / `permutation_keeps_context_block` | 0.0 / 0.0 |

**其余已关闭项（复跑仍在）**：冻结同点预测差 **1.69e-5 ≤ 1e-4**、每片 RMSE 差 **2.04e-6 ≤ 1e-5**；置换保持 nnz/度/权重/对称性与置换图孤立点 self fallback（差 0.0，孤立点数 (2,2)）；篡改留出目标后设计矩阵位级不变；细节分量 Pearson 与手工一致；q75 严格大于分位（85% 零值场 = 非零点 150）；范围绝对偏差真实可评例 (-6.245, 6.245)、相消反例 (mean_abs=0.266667, abs_mean=0.0)、0 计入分母 (3,3)、缺失片 NOT_ESTIMABLE (None, True)；聚合与已审模块一致 (0.2)、单位键不含 source_key、元数据冲突 `mixed`。

**短贯通（截断训练 8 片 / 测试 2 片 / 2 个划分）**：fit 内部折 46.0 s + 来源折 41.6 s；evaluate 0.8 + 0.5 s；summarize 0.6 s；80 行 section 行、80 行 unit 行、1680 行 macro、896 行 paired、40 行 claim；范围读出 **80/80 可评**；诊断显示 point (feats=1649, scaled=0, ctx=6) 与三块模型 (feats=4935, scaled=4929, ctx=6) ✓。

## 5. 本轮新失败与修复（继续保留）

| 轮次 | 现象 | 原因与修复 |
|---|---|---|
| check（R-01 修复中） | `NameError: CONTEXT_DIM` | 常量未定义 → 在常量区定义 `CONTEXT_DIM` |
| check | `TypeError: tuple indices must be integers`（`stack_training(...)['point_ridge']`） | `stack_training` 已改为返回 `(stacks, widths)` → 更新冻结一致性探针调用点 |
| check | `TypeError: unhashable type: 'slice'` | 某处仍传 `section_designs(...)['point_ridge']` 字典而非 `['design']` → 修正 |
| check | `DeprecationWarning: Conversion of an array ... n_iter_` | lsqr 的 `n_iter_` 是数组 → 用 `np.ravel(...)[0]` |

以上失败日志分别保存在 `archive/check_stale_stack_call_r01_fix.log`、`archive/check_design_dict_slice_failure.log` 等；`stages.log` 保留全部返回码（含本轮 3 次 `check exit=1`）。首轮审核文件与既有失败归档均未删除。

## 6. 全量成本（按实测外推，已按 2500 模型行口径）

- fit：真实规模实测 point 2.9 s、三块模型 11.8 s（28 488 训练行，8 线程，cholesky）；按各划分训练行数缩放合计 ≈**25 min CPU**；`workers=2` 时 ≈**12–15 min 墙钟**（最长单划分 ≈2.7 min + 特征/预测开销）。
- evaluate：**2500** 模型行 × ≈0.03 s ≈ **1.4 min CPU**（旧 5000 行口径减半）。
- summarize ≈1 min；磁盘 <0.5 GB（合同预算 4 GB）；内存单 worker 峰值 ≈2 GB（本机可用 691 GB）。

## 7. 疑点 / NOT_RUN（交复核）

1. **本版本尚未被独立审核**；首轮审核文件只为旧版本（`7ba553f`/`1758023`）背书，不替本版本背书。
2. 求解器仍是 `auto`（cholesky）；在该病态设计上 cholesky 与 lsqr 预测差 ≈3.5e-6（相对 ≈3e-5）。若复核要求更高数值可复现性可换 lsqr（约 3× 时间），本轮未改（不得按分数选参）。
3. 分块标准化后空间块与同点块容量并不相等（合同已声明）；1/√3 只抵消复制通道造成的有效正则减弱。
4. 置换为单次机制破坏对照，不作 p 值/零模型。
5. 本轮关闭 shuffle 参考（`do_shuffle=False`）：`*_shuffled_*`、`rec_jaccard_q75_excess` 不产出。
6. 人工 TLS 验证、独立 marker 隔离、未知患者身份、正式空间显著性、物理尺度/视野外推断：NOT_RUN。
7. 小样本贯通为截断训练集（8 片），不是任何折的真实结果；全量未运行。
8. `readout_claim_summary.tsv` 为描述性合取（无阈值、无 p 值），小样本 40 条中 1 条为真，样本截断下不可解释。

## 8. 复核所需文件

- 实现与合同：`scripts/field_raw_gene_spatial_20261019.py`、`scripts/run_field_raw_gene_spatial_20261019.sh`、`infra/field_raw_gene_spatial_20261019/contract.json`。
- 证据：`infra/field_raw_gene_spatial_20261019/checks/checks.json`（49/0）、`.../checks/solver_diagnostics.tsv`、`.../smoke/{smoke_summary.json,REPORT.md,macro_summary.tsv,macro_paired.tsv,unit_metrics.tsv,prediction_metrics.tsv,paired_metrics.tsv,readout_claim_summary.tsv,denominators.tsv}`、`.../input_checks.json`、`.../prepare_receipt.json`。
- 失败与历史：`.../stages.log`、`.../archive/*.log`、`review/reviewed_implementation_round1.json`、`review/preflight.md`。
- 版本范围：本次提交与 `git show --stat`；未提交大文件在 §2 列明。

## 9. 交付后冻结

本文件交付后停止修改实现/机器合同/检查/交接与既有小结果，等待独立复核；复核通过前**不运行全量**，不宣布科学结论。
