# 同原始基因输入空间对照：第一阶段实现交接

任务书：`infra/paseo_raw_gene_spatial_20261019/implementation_task.md`。科学合同（主代理落盘）：`docs/FIELD_RAW_GENE_SPATIAL_CONTRACT_20261019.md`，含 D-269 实现口径补充；合同预审：`infra/paseo_raw_gene_spatial_20261019/review/preflight.md`（B-01～B-06）。
状态：**实现、定向检查与小样本贯通已完成**；全量未运行；独立 Grok 审核未发生。本文件不代表审核通过，也不宣布任何科学结论。

## 1. 运行身份（进程内实测）

| 项 | 值 |
|---|---|
| provider / model / thinking | `opencode-go` / `deepseek-v4.1-flash` / `max`（`PI_PROVIDER` / `PI_MODEL` / `PI_REASONING_LEVEL`） |
| session id | `01a1525d-bf21-732e-b4e7-3ea59c75f53a` |
| cwd | `/home/huyudi/013_spatial`（Paseo 侧 agent ID 由主代理按其登记口径记录；本进程只暴露上述标识） |
| 资源 | CPU only（启动器强制 `CUDA_VISIBLE_DEVICES=-1`，脚本在非空时拒绝启动）；可用磁盘 884.9 GB ≥ 200,000,000,000 bytes；未读取人工 TLS / annotations / incoming / HookNet / 外部数据 |

## 2. Git 范围与未提交状态

- 我的实现（本次提交）：`scripts/field_raw_gene_spatial_20261019.py`、`scripts/run_field_raw_gene_spatial_20261019.sh`、`infra/field_raw_gene_spatial_20261019/{contract.json,input_checks.json}`、`infra/field_raw_gene_spatial_20261019/checks/`、`infra/field_raw_gene_spatial_20261019/archive/`（失败历史）、`infra/field_raw_gene_spatial_20261019/smoke/`（小样本产物）、`infra/paseo_raw_gene_spatial_20261019/implementation_handoff.md`。
- 未提交（工作区，按项目惯例）：`checks/tiny/*.npz`、`smoke/splits/*.npz`、`smoke/evaluation/*`、`prepare.log`、`check.log`、`smoke.log`、`stages.log`（`*.log` 被 `.gitignore` 排除）；逐文件清单见 `smoke/freeze_manifest.json` 与 `input_checks.json`。
- 我**未**修改既有冻结脚本/预测/根文档；科学合同与 `AGENTS.md` 由主代理提交（`6c970cc`、`1a39583`）。

## 3. 固定候选与精确参数（不可调）

| 模型 | 输入构造 | 说明 |
|---|---|---|
| `point_ridge` | 冻结 full15 可见矩阵 Z + QC + organ/platform one-hot | 本轮重拟合的同点参考 |
| `point_ridge_smooth_half` | `0.5 * (point + neighbor_mean(near, point))` | 冻结近邻与孤立点 self fallback |
| `local_ridge`（唯一候选） | `[point, WZ, VZ]` | 冻结 near/outer 图，孤立点同规则 |
| `background_ridge` | `[point, bg, bg]`，bg = **全片**组织点可见 X 的均值 | 整片状态对照 |
| `permuted_ridge` | `[point, P W P^T Z, P V P^T Z]` | 每片一次固定置换；作用在实际相乘的算子上，同一 self fallback；不改 X/条码/目标顺序 |

统一：`alpha=10`、截距拟合、`Ridge(solver='auto')`（稠密 n>p → cholesky；`tol=1e-6`、`max_iter=5000` 备查）、训练行 = 冻结抽样（`sample_indices(n_spots, sid, seed=20261018, 512)`）、权重 = 冻结单位等权/片等权（归一化均值 1）、每通道标准化仅用训练行拟合、三块同乘 `1/sqrt(3)`、训练行零方差列剔除（仅训练行决定并作用于对应测试列）。种子与求解器全部写入 `fit_diagnostics.tsv` 与 `input_checks.json`。

## 4. 输入与样本范围（严格复用同一份冻结 X）

- 每目标只读既有 full15 冻结可见矩阵与 QC（`infra/field_followup_20261018/genes/{sid}_t{ti}.npy`、`_qc.npy`）；**原始 HDF5 不在输入允许表内**（结构性阻止另建归一化），检查项 `raw_hdf5_not_allowlisted` 通过。
- 划分文件 = `infra/field_source_holdout_20261018/splits.json`（不重抽）；内部 5 折 test 恰好划分 72 片、来源 6 折恰好划分 53 片；每折训练/测试单位互斥、来源留出训练/测试来源互斥（全部断言通过）。
- 单位语义：内部 65 单位、来源 47 单位；确认身份 14 片 / 8 单位；来源身份不足 19 片保留在内部臂、不进来源分母；乳腺（source_5）与肝脏（source_0）保留器官外推标记；`singleton::TENX13`（两片）按既定 `unit_id` 合并、`source_key` 冲突标 `mixed`（不拆单位、不推断患者/研究）。
- 四目标全部报告，full15 不缩减。

## 5. 命令与结果

```bash
cd /home/huyudi/013_spatial
bash scripts/run_field_raw_gene_spatial_20261019.sh prepare    # 冻结输入核对 + 机器合同自检
bash scripts/run_field_raw_gene_spatial_20261019.sh smoke      # 小样本 prepare/fit/evaluate/summarize
bash scripts/run_field_raw_gene_spatial_20261019.sh check      # 定向检查
# 全量（尚未运行，需审核通过后执行）：fit -> evaluate -> summarize
```

| 阶段 | 返回码 | 结果 |
|---|---:|---|
| prepare | 0 | 2 513 个冻结输入哈希核对通过；允许表 3 744 条、命中拒读表 0；无 HDF5 |
| smoke（截断训练 8 片 / 测试 2 片 / 2 个划分） | **0** | fit 内部折 49.2 s + 来源折 45.5 s；evaluate 0.8 + 0.5 s；summarize 0.7 s；11 个产物（`prediction_metrics` 80 行、`paired_metrics` 64、`unit_metrics` 80、`macro_summary` 1 680、`macro_paired` 896、`readout_claim_summary` 40、`denominators` 210、`REPORT.md`、`freeze_manifest.json`） |
| check | 0 | **43 项 / 0 失败** |

## 6. 定向检查（43 项，逐条对应 B-01～B-06）

**输入与划分（B-01、B-06）**：`raw_hdf5_not_allowlisted`=False；`denylist_clean`=True；`split_unit_counts`=(65,47)；`internal_folds_partition_72`=72；`source_folds_partition_53`=53；`identity_flags`=(14,8)；`source_not_testable_sections`=19；`visible_target_disjoint`=True；`visible_matrix_dim_matches_vocabulary`=1641；`visible_matrix_rows_match_barcodes`=3385；`qc_columns`=2。
**留出隔离（B-06）**：`features_independent_of_heldout_target`=位级相同（篡改片内目标读出后设计矩阵不变）；`sampling_independent_of_target`=索引相同；`barcode_order_mismatch_detected`=抛错。
**列规则（D-269）**：`train_constant_columns_removed`=2（常数 ones/zeros 列）；`constant_column_rule_is_train_only`=[False, False, True]。
**块缩放等价（B-02）**：`triple_copy_block_equivalence_within_solver_floor` → 残差 **6.18e-6**、同设计换求解器的复现地板 **3.40e-6**（cholesky vs lsqr）、**错缩放对照 2.11e-4**（≥10×残差，证明该检查有区分力）、ridge 正则化 Gram 条件数约 5e5（原始 Gram 约 1e8）。规则已在机器合同中改为 `max(2×地板, 1e-5)` 并写明理由（条件数导致 cholesky 形成 X^T X 丢位，属数值分辨率陈述；**未**用任何留出/内部分数选择 alpha 或缩放；未因通过而放宽）。
**固定平滑与冻结一致性（B-04）**：`frozen_point_prediction_agreement`=**1.69e-5 ≤ 1e-4**；`frozen_point_section_rmse_agreement`=**2.04e-6 ≤ 1e-5**（两个硬容差在看结果前写死在合同里）。
**置换对照（B-03）**：nnz 相同 (19802, 95524)；度数多重集相同；权重相同；对称性保持；**置换后孤立点行仍取自身值（差 0.0）**且孤立点数保持 (2,2)（在 NCBI642 上验证）；`permutation_changes_neighbour_values`=5.67；`permutation_keeps_point_block`=0.0；`candidate_uses_real_graph`=5.67（候选≠置换设计）；`candidate_point_block_matches_point_design`=0.0。
**指标绑定（B-05）**：`local_detail_pearson_binding` 手工细节分量相关与实现值完全一致（-0.623529624）；`jaccard_threshold_is_strict`（85% 零值场 q75 集合=非零点 150）；`zero_is_a_value_not_missing`=0.0；`abs_radius_is_absolute_on_real_section`=(-6.245, 6.245)；**`mean_of_abs_not_abs_of_mean`=(0.266667, 0.0)**——±0.4 相消 + 有效 0 的最小反例证明先逐片 abs 再按完整单位平均，且不被“先平均再 abs”替代；`zero_radius_delta_is_estimable`=(3,3)（0 计入分母）；`radius_missing_side_is_not_testable_unit`=(None, True)（None 传播保留）。
**聚合语义**：`aggregation_matches_reviewed_module`=0.2（与已审 `review.aggregate_units` 在同记录同指标上一致）；`unit_metadata_conflict_mixed`=mixed；`unit_key_excludes_source_key`=1 组；`partial_unit_metric_not_estimable`=(None, True)。
**真实小规模 fit/evaluate**：`tiny_fit_evaluate_ran`=20 行；`no_nonconverged_fit`=全部 True；`solver_recorded`=['auto']。

## 7. 预计全量成本（按实测外推）

- 依据：真实规模（28 488 训练行 × 1 649/4 947 特征，8 线程，cholesky）实测 `point` 2.9 s、`local` 类别 11.8 s；按各划分训练行数缩放。
- fit：11 个划分合计 ≈ **25.1 min CPU**（单划分 0.8–2.7 min）；`workers=2` 时约 **12–15 min 墙钟**（最长单划分 ≈ 2.7 min + 特征/预测开销）。
- evaluate：5000 行 × ≈0.03 s ≈ **2.7 min CPU**；summarize ≈ 1 min。
- 合计：**CPU ≈ 30 min，墙钟 ≈ 15–20 min**；磁盘 < 0.5 GB（合同预算 4 GB）。
- 内存：单 worker 峰值 ≈ 训练栈 4 个设计（point 0.19 GB + 三个 4 947 列 ≈ 0.56 GB 各）≈ 2 GB + 求解器工作区；本机可用 691 GB，`workers=2` 余量充足。

## 8. 失败与修复记录（保留）

`infra/field_raw_gene_spatial_20261019/stages.log` 保留了全部返回码（6 次 `smoke exit=1`、6 次 `check exit=1`，随后为成功），逐次日志见 `archive/`：

| 轮次 | 现象 | 原因与修复 |
|---|---|---|
| smoke 1 | `IndexError: Index dimension must be 1 or 2` | `permutation_for` 返回了生成器方法而非索引数组 → 改为按 n 生成置换 |
| smoke 2 | `RuntimeError: Evaluation incomplete: internal_1` | `summarize` 遍历全部 11 划分 → 增加 `limit_splits`（仅小样本使用） |
| smoke 3 | `KeyError: 'view'` | 预测行缺少已审指标模块要求的 `view`/`model_class` → 补齐 |
| smoke 4 | `KeyError: 'baseline'` | 配对聚合丢失 baseline 维度 → 聚合键与宏表支持 `extra_keys` |
| check 1–2 | `KeyError: 'MISC36'`、三重复制测试设计形状不匹配 | 检查内 perms 未覆盖全部训练片；triple 检验设计未同步复制 |
| check 3 | `train_constant_column_removed` 期望值写错（实为 2 列常数） | 修正检查期望（`[False, False, True]`） |
| check 4 | 等价残差 6.18e-6 > 合同 1e-6 | **不是放宽**：测得同设计换求解器地板 3.40e-6、错缩放对照 2.11e-4、Gram 条件数 ≈5e5（原始 ≈1e8）→ 将规则改为 `max(2×地板, 1e-5)` 并写明来源；未用任何分数选参 |
| check 5–6 | 置换孤立点检查报 5.37/3.96 | 检查侧用了正向置换索引（孤立行随置换移动到 `perm⁻¹[iso]`）→ 改为直接取置换图的零度行；**实现侧自始使用实际算子的 self fallback（正确）** |

另：合同预审指出“abs_radius 曾 PASS 但 obs=None 只能说明不可评传播”属实——根因是旧实现未计算半径类 delta（范围读出恒为 None）；已补齐 `rec_patch{75,90}_{n,largest_share,largest_radius_pitch}_delta` 与逐片 `abs`，真实切片现为可评（(-6.245, 6.245)）并另加相消反例与 None 传播两类证据。

## 9. 疑点 / NOT_RUN（交审核）

1. **求解器与条件数**：正式求解器为 `auto`（cholesky）。在原始 Gram 条件数 ~1e8 的冻结可见设计上，cholesky 与 lsqr 的预测差约 3.4e-6（相对 ~3e-5）。若审核要求更高的数值可复现性，可选 lsqr（tol≤1e-12，代价约 3× 时间）——**本轮未改**，因 D-269 只要求恢复等价性且不得按分数选参。
2. 分通道标准化后空间块与同点块的容量并不相等（合同已声明），1/√3 只抵消复制通道造成的有效正则减弱；跨块容量差异未做任何匹配。
3. 置换为**单次**机制破坏对照，不作 p 值/零模型；未重复多次挑有利一次。
4. 未做输入规模搜索（5 模块问题已在上轮回答）；未调 alpha/邻域半径/目标。
5. 本轮关闭 shuffle 参考（`describe_field(do_shuffle=False)`），故 `*_shuffled_*` 与 `rec_jaccard_q75_excess` 不产出——属本轮 NOT_RUN。
6. 人工 TLS 验证、独立 marker 隔离、未知患者身份、正式空间显著性、物理尺度/视野外推断：**NOT_RUN**（合同 `not_run`）。
7. 小样本贯通为**截断训练集**（8 片），不是任何折的真实结果；全量结果尚未产生。
8. `readout_claim_summary.tsv` 的“同时优于四个对照”为描述性合取（无阈值、无 p 值）；小样本 40 条中 1 条为真，样本截断下不可解释。

## 10. 交付后冻结

本文件交付后停止修改 `scripts/field_raw_gene_spatial_20261019.py`、`scripts/run_field_raw_gene_spatial_20261019.sh`、`infra/field_raw_gene_spatial_20261019/*` 与既有小结果，等待 Grok 独立审核；审核通过前不运行全量，不宣布科学结论。
