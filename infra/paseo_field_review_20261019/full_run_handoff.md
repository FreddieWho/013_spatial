# 全量描述性运行交接（full_run_handoff）

任务书：`infra/paseo_field_review_20261019/full_run_task.md`。授权依据：`review/reviewed_revision.json` = `ACCEPTED_FOR_DESCRIPTIVE_RUN_ONLY`，`allow_full_descriptive_analysis = true`，未关闭结果相关问题为空；其绑定的脚本/协议/合同/烟测驱动/检查记录/烟测清单 6 项在启动前与磁盘逐一核对一致（全部 match）。
版本管理：按 D-265 起用项目 Git 记录范围，**不再维护手工哈希清单**。本文件不替主代理判定任何生物学结论。

## 1. 执行结果（三阶段全部成功）

| 阶段 | 命令 | 返回码 | 墙钟 | CPU 秒 | 日志 |
|---|---|---:|---:|---:|---|
| prepare | `bash scripts/run_field_biological_review_20261019.sh prepare` | **0** | 130 s | 138.8 s | `infra/field_biological_review_20261019/prepare.log` |
| run | `bash scripts/run_field_biological_review_20261019.sh run` | **0** | 115 s | — | `.../run.log`（11 个划分，asks 见 §4） |
| summarize | `bash scripts/run_field_biological_review_20261019.sh summarize` | **0** | 12 s | — | `.../summarize.log` |

串行编排与返回码留痕：`infra/field_biological_review_20261019/full_run_stages.log`（`prepare exit=0 seconds=130` / `run exit=0 seconds=115` / `summarize exit=0 seconds=12`）。
启动器固定 `PYTHONNOUSERSITE=1`、`OMP_NUM_THREADS=2`、`OPENBLAS/MKL_NUM_THREADS=1`、`CUDA_VISIBLE_DEVICES=-1`；全程 CPU、本地数据、未读取人工 TLS / annotations / incoming（`input_checks.json` 记录 `human_tls_files_opened: []`，拒读表与允许表白名单双检通过）。

## 2. Git 范围与未提交状态（D-265）

- 分支/基线：`main`，启动前 HEAD = `5a2d97a`。
- **提交 A（代码·协议·合同·检查证据）**：`15dad63` “场生物学复核：协议、脚本与检查证据入库（D-265 起用 Git 记范围；全量描述性运行的代码/合同冻结版）”——18 个文件；含 `scripts/field_biological_review_20261019.py`、`scripts/smoke_field_biological_review_20261019.py`、`scripts/run_field_biological_review_20261019.sh`、`docs/FIELD_BIOLOGICAL_REVIEW_PROTOCOL_20261019.md`、`infra/field_biological_review_20261019/contract.json`、`checks/`（含历史失败归档）、`smoke/smoke_summary.json`、`full_run_derived/denominators_by_arm_model_metric.tsv`。
- **提交 B（本文件与运行台账）**：`git log -1`（提交信息标为“全量描述性运行交接与台账”），含本文件、`input_checks.json`、`freeze_manifest.json`、`REPORT.md`、`STATUS`、`prepare_receipt.json`、`measured_fields.tsv`、`measured_summary.tsv`、`within_unit_consistency.tsv`、`cross_stage_consistency.json`、`diagnostic_sensitivity.tsv`、`full_run_derived/`。
- **被 `.gitignore` 排除（`*.log`，故不入 Git，仅存留工作区）**：`prepare.log`、`run.log`、`summarize.log`、`full_run_stages.log`、`smoke.log`；三阶段返回码同时写在本文件 §1 与上述日志尾部。
- **未提交（有意不入 Git，保持与既有三轮运行目录同样的处理）**：体积较大的逐行表与逐划分目录——`prediction_metrics.tsv`(11 MB)、`prediction_module_relations.tsv`(27 MB)、`macro_summary.tsv`(4.5 MB)、`macro_strata.tsv`(13 MB)、`unit_metrics.tsv`(4.6 MB)、`measured_module_relations.tsv`(4.2 MB)、`splits/`(38 MB)、`smoke/` 内大表(13 MB)。全部路径与逐文件 sha256 见**已提交的** `infra/field_biological_review_20261019/freeze_manifest.json`（148 项）；这些文件仅存在于本工作区，未提交状态如实注明。
- 其他协作者在同一工作区的未提交改动（`AGENTS.md`、`STATUS.md`、`docs/*`、`TODO.md`、`LEADS.md` 等）**未被本次提交包含、未被修改**。
- 范围查看方式：`git show --stat 15dad63`、`git show --stat <提交 B>`；与 `5a2d97a` 的全量工作区差异 `git status --short`。

## 3. 产物清单（真实行数/范围）

| 产物 | 行数 | 范围 | 大小 |
|---|---:|---|---:|
| `measured_fields.tsv` | 288 | 72 片 × 4 目标 | 509 KB |
| `measured_module_relations.tsv` | 5 760 | 288 × (full15 15 模块 + small5 5 模块) | 4.2 MB |
| `prediction_metrics.tsv` | 5 000 | 1000 冻结预测包 × 5 模型 | 11 MB |
| `prediction_module_relations.tsv` | 50 000 | 5000 × 10 模块均值 | 27 MB |
| `unit_metrics.tsv` | 4 520 | 单位 × 5 模型（内部 2 640 + 来源 1 880） | 4.6 MB |
| `macro_summary.tsv` | 19 240 | 划分级 16 280 + 臂级 2 960 | 4.5 MB |
| `macro_strata.tsv` | 69 560 | 分层行（其中 13 570 行显式标 `no_estimable_units_in_stratum`） | 13 MB |
| `within_unit_consistency.tsv` | 1 096 | 7 个双片单位（三种 `profile_column`） | 308 KB |
| `measured_summary.tsv` | 128 | 目标 × 器官/平台/来源/资格分层 | 77 KB |
| `diagnostic_sensitivity.tsv` | 9 | 3 个诊断标记 × 3 个相邻切点 | 1.1 KB |
| `cross_stage_consistency.json` | — | `compared=30030, max_abs_difference=0.0` | 53 B |
| `REPORT.md` / `STATUS` | — | 描述性小结（含 §6.1 边界）；`DESCRIPTIVE_REVIEW_COMPLETE_NO_HUMAN_TLS_NO_PASS_THRESHOLD` | 17.6 / 0.06 KB |
| `freeze_manifest.json` | — | 148 项；`human_tls_labels_read=false`、`biological_pass_threshold=null` | 23.7 KB |

冻结清单一致性核对：148 项中仅 `summarize.log` 与 `full_run_stages.log` 两个**日志**在哈希写入后仍有输出（自身 stdout 与串行编排追加），属预期；全部数据产物逐文件 sha256 与清单一致。

## 4. 分账：臂、单位、来源、身份

| 层 | 内部臂 | 来源臂 |
|---|---|---|
| 切片 | 72 | 53 |
| 唯一单位（`unit_id`） | 65 | 47 |
| 划分 | internal_0..4（15/13/13/13/11 单位） | source_0..5（2/6/2/2/34/1 单位） |
| 每（目标×视图×模型）单位组 | 66（见下方说明） | 47 |
| 来源身份不足（`NOT_TESTABLE_SOURCE_IDENTITY`） | 19 片保留在实测清单与内部臂 | 不进入来源分母 |
| 器官外推（`cross_organ_extrapolation`） | 不适用 | 乳腺 (source_5 1 片)、肝脏 (source_0 2 片)：`same_organ_train_sections=0`；臂级标签 `false`/`mixed` 各 1 480 行 |
| 平台 | Visium 252 + Visium HD 36（行数，含 4 目标） | 同上（按来源切片） |
| 身份 | 确认 14 片 / 8 单位；其余 57 单位为未知，未升级为独立患者 | 来源臂沿用同一身份表 |

**必须如实说明的单位计权细节（本次全量运行新暴露，未修改已验收代码）**：`singleton::TENX13` 的两片 `TENX13`/`TENX14` 的 `metadata` 标题不同，冻结 manifest 由此得到不同 `source_key`（`title:Human Breast Cancer (Block A Section 1/2)`）。聚合键包含 `source_key`，故该单位在内部臂被拆成 **2 个 1 片单位组**（因而内部臂每组 66 而非 65）。影响量化（内部臂臂级宏平均，40 个目标×视图×模型组合汇总，把两片先合并为单位再平均 vs 现状）：`rec_rmse` Δ=+2.46e-04、`rec_centered_rmse` Δ=+8.8e-05、`rec_pearson` Δ=−2.1e-04、`rec_jaccard_q75` Δ=+4.0e-04、`rec_centroid_displacement_pitch` Δ=+3.4e-03（`gene_ridge`）；`fixed_local_hgb` 同类量级（`rec_jaccard_q75` Δ=+8.2e-04、`rec_centroid_displacement_pitch` Δ=+4.9e-03）；`train_mean` 的 RMSE 类 Δ≈+3.2e-04。即**一个单位被赋予双倍单位权**，相对偏差约 0.01%–0.2%。是否需要在合同层合并同 `unit_id` 跨 `source_key` 的切片，属需审核决定的问题，本轮未改。

## 5. 完整 / 部分 / 不可评 分母

逐（臂 × 模型 × 指标）的完整表：`infra/field_biological_review_20261019/full_run_derived/denominators_by_arm_model_metric.tsv`（也可从已提交的 `unit_metrics.tsv` 的 `<metric>__n`、`partial_unit_metrics`、`not_estimable_metrics` 复算）。

关键事实：

- **无部分可评单位落在恢复类指标上**：`partial_unit_metrics` 全库仅 26 个单位行，且全部是 `geo_decay_half_pitch`（1 个双片单位、一个切片可评另一个不可评）；恢复/关系类指标的单位要么整单位可评、要么整单位不可评。
- `rec_rmse` / `rec_centered_rmse` / `rec_mae` / `rec_bias`：内部 528、来源 376 个单位组**全部可评**（各模型）。
- 相关类（`rec_pearson`、`rec_spearman`、`rec_local_pattern_pearson`、`rec_jaccard_q75/q80/q90`、`prof_relation_pearson`、`rec_autocorr_*_delta`、`rec_centroid_displacement_pitch`）：内部 528 + 来源 376 单位组全部可评（4 个非 `train_mean` 模型）；`train_mean` 因常量化按合同 `sd_floor` 结构性不可评（内部 528、来源 376 全缺失）。
- **`rec_decay_half_pitch_delta` 全库不可评**（内部 528 + 来源 376 全缺失）：因为实测 `geo_decay_half_pitch` 在 288 个实测场中仅 2 个可评（其余 286 个为「首带 ≤0」或「band 上限内未衰减」）。空间尺度在本轮由 `geo_corr_band1..6` 与 `geo_patch75/90_*_radius_pitch` 承担，**单个衰减尺度数不可用**。
- `train_mean` 的完整类指标为 RMSE/MAE/bias/`rec_sd_ratio`（=0）与部分几何列；不进入相关排名。

## 6. 实测场资格与连续证据（72 片 × 4 目标）

| 目标 | 任意目标基因检出点占比 中位（最小） | 有效基因数 中位（最小） | 检出基因数 中位 | 实测 sd < 1e-6 的场 |
|---|---|---|---|---|
| B_AXIS | 0.999（0.612） | 2.90（1.19） | 47 | 0 |
| EMT | 1.000（0.962） | 32.83（7.66） | 199 | 0 |
| IFN-γ | 1.000（0.950） | 24.74（10.83） | 195 | 0 |
| OXPHOS | 1.000（0.919） | 75.41（4.20） | 191 | 0 |

- 实测状态：`DESCRIBED` 2、`DESCRIBED_WITH_NOT_ESTIMABLE_METRICS` 286、`NOT_TESTABLE_*` **0**（没有实测场因常量/无目标基因被判不可评）。
- 诊断标记（非门槛、不删行）：首位基因占比 >0.5 共 21 场；有效基因数 <2 共 18 场（与前者同场出现 18 场）；切点敏感性见 `diagnostic_sensitivity.tsv`。
- 预测状态：5 000 行中 `NOT_ESTIMABLE_CONSTANT_PREDICTION` 1 000（全部 `train_mean`）、`DESCRIBED_WITH_NOT_ESTIMABLE_METRICS` 3 994、`DESCRIBED` 6。
- 身份：确认 14 片 / 8 单位；未知 57 单位，切片数未冒称患者数。

## 7. NOT_RUN / NOT_TESTABLE 全表

| 项 | 状态 | 依据 |
|---|---|---|
| 人工 TLS 结构验证与人工评价 | **NOT_RUN** | 全程未读取（拒读表 + 允许表，`human_tls_files_opened: []`） |
| 独立 marker 隔离与重训 | **NOT_RUN** | 既有模型用完整目标基因集训练，事后拆分不构成独立留出 |
| 任何训练/重训/调参/新模型 | **NOT_RUN** | 只重读 1000 份冻结预测 |
| 正式空间 p 值 / 显著性 | **NOT_RUN** | 打乱为未校准的描述性扰动，不产生 p 值 |
| 复杂零模型、密度/支撑保持校准 | **NOT_RUN** | 合同 `not_run` |
| 校准不确定性、视野外/三维恢复 | **NOT_RUN** | 合同 `not_run` |
| 来源身份不足的 19 片来源留出评价 | **NOT_TESTABLE_SOURCE_IDENTITY** | 保留在实测清单与内部臂 |
| 肺部来源留出 | **NOT_TESTABLE** | 6 片肺全部在来源身份不足集合内 |
| 乳腺/肝脏来源留出 | 器官外推（`cross_organ_extrapolation=true` 的来源） | 留出时无同器官训练切片 |
| `rec_decay_half_pitch_delta` | **NOT_ESTIMABLE 全库** | 实测衰减尺度仅 2/288 可评（§5） |
| `train_mean` 的相关类指标 | **NOT_ESTIMABLE**（结构性） | 常量预测，`sd_floor` 规则 |
| 外层裹 `geo_decay_half_pitch` 的部分单位 | 部分可评（26 个单位行） | 已按 R-02 规则记 `partial_unit_metrics`，不进等权均值 |
| 物理尺度/微米换算/超分辨率 | **NOT_RUN，且不主张** | 仅以本片经验点间距为单位 |

## 8. 失败与新发现

1. **全量三阶段无失败**（返回码全 0）；无分区中止、无空表、无 `NOT_ESTIMABLE_UNIT` 误标。
2. 历史失败（本轮修复前，已归档、未删除）：`checks/archive/` 下 7 次运行记录，其中 `run-13923425088ea8`、`run-08c20cf71466c8` 各 12 项失败（隔离点邻域语义、容差分层）；`checks/archive/smoke_fail_protocol_not_allowlisted.log`（`summarize` 读协议被允许表拒绝，已修）。
3. **本次全量暴露的新问题（未修改已验收代码，交审核）**：`unit_id` 相同但 `source_key` 不同的切片被拆成两个单位组，使 `singleton::TENX13` 获得双倍单位权（§4 末段，影响 0.01%–0.2%）。是否合并属合同层决定。
4. 次要记录：`freeze_manifest.json` 的两条日志项（`summarize.log`、`full_run_stages.log`）在哈希写入后仍继续追加，故与磁盘不一致；数据产物全部一致。

## 9. 资源实际使用

- 磁盘：启动前可用 **888.51 GB**；`summarize` 记录 **888.4 GB**；结束后 **888.22 GB**（净增约 0.29 GB，远低于合同 `estimated_max_new_bytes = 3e9`）。全程高于 200,000,000,000 bytes 下限。
- CPU/时间：prepare 130 s 墙钟（138.8 CPU 秒，3 进程）；run 115 s；summarize 12 s；合计约 **4 分 17 秒**墙钟（与烟测外推的 ~13 min CPU 上限一致）。
- 输入核对：`input_checks.json` 记录 `frozen_hashes 2948`、`expression_hashes 72` 全部通过；面板文件与冻结 discovery 收据一致；白名单 3 738 条、命中拒读表 0 条。
- GPU：无（启动器强制 `CUDA_VISIBLE_DEVICES=-1`，脚本自身在非空时拒绝启动）。

## 10. 边界（不替主代理判定）

本轮产物是**既有开发材料的描述性再分析**：无生物学成功/失败二分门槛、无正式空间 p 值、无独立 marker 验证、无 TLS 真值；预测与可见模块的相关含输入回声，不作为新关系；**数学上可算不等于有足够生物学参照**（低覆盖、单基因主导、边界支撑、未确认物理尺度仍限制解释）。主代理据此分析 HF-01/HF-02/HF-04 的现有边界，本文件不宣布任何顶层假设成立。
