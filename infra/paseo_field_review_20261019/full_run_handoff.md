# 全量描述性运行交接（full_run_handoff）

任务书：`infra/paseo_field_review_20261019/full_run_task.md`。授权依据：`review/reviewed_revision.json` = `ACCEPTED_FOR_DESCRIPTIVE_RUN_ONLY`，`allow_full_descriptive_analysis = true`，未关闭结果相关问题为空；启动前其绑定的脚本/协议/合同/烟测驱动/检查记录/烟测清单 6 项与磁盘逐一核对一致。
版本管理：按 D-265 用项目 Git 记录范围，**不维护手工代码哈希清单**。本文件不替主代理判定任何生物学结论。

**修订 1（聚合键修复）**：全量交接后暴露结果相关问题——`singleton::TENX13` 的两片因切片标题不同（`source_key` 不同）被拆成两个单位组，内部臂 65 个既定单位变成 66 个。已按现行「完整单位等权」合同修复**新分析代码的聚合键**，并复用既有全量表重做汇总（仅重跑 `summarize`，未重跑 prepare/run、未改冻结预测/历史脚本/根文档）。修复内容、验证与前后影响见 §C 与 §D。

## 1. 执行结果（三阶段全部成功）

| 阶段 | 命令 | 返回码 | 墙钟 | CPU 秒 | 日志（`*.log` 被 .gitignore 排除，仅存工作区） |
|---|---|---:|---:|---:|---|
| prepare | `bash scripts/run_field_biological_review_20261019.sh prepare` | **0** | 130 s | 138.8 s | `infra/field_biological_review_20261019/prepare.log` |
| run | `bash scripts/run_field_biological_review_20261019.sh run` | **0** | 115 s | — | `.../run.log` |
| summarize（初版） | `bash scripts/run_field_biological_review_20261019.sh summarize` | **0** | 12 s | — | `.../summarize.log` |
| summarize（聚合键修复后重跑） | 同上 | **0** | 14 s | — | `.../summarize_fix.log` |

编排留痕：`infra/field_biological_review_20261019/full_run_stages.log`（含 `summarize(re-run) exit=0 seconds=14`）。
启动器固定 `PYTHONNOUSERSITE=1`、`OMP_NUM_THREADS=2`、`OPENBLAS/MKL_NUM_THREADS=1`、`CUDA_VISIBLE_DEVICES=-1`；全程 CPU、本地数据、未读取人工 TLS / annotations / incoming（`input_checks.json`：`human_tls_files_opened: []`）。

## 2. 产物清单（修复后真实行数/范围）

| 产物 | 行数 | 范围 | 是否随修复改变 |
|---|---:|---|---|
| `measured_fields.tsv` | 288 | 72 片 × 4 目标 | 否（字节相同） |
| `measured_module_relations.tsv` | 5 760 | 288 ×（full15 15 + small5 5 模块） | 否 |
| `prediction_metrics.tsv` | 5 000 | 1000 冻结包 × 5 模型 | 否 |
| `prediction_module_relations.tsv` | 50 000 | 5000 × 10 模块均值 | 否 |
| `unit_metrics.tsv` | **4 480** | 单位 × 5 模型（内部 2 600 + 来源 1 880） | **是**（修复前 4 520） |
| `macro_summary.tsv` | 19 240 | 划分级 16 280 + 臂级 2 960 | **是**（数值） |
| `macro_strata.tsv` | **78 440** | 分层行（含 `source_split` 新层；15 294 行带空层 reason） | **是**（行数与数值） |
| `within_unit_consistency.tsv` | 1 096 | 7 个双片单位（三种 `profile_column`） | 否（字节相同） |
| `measured_summary.tsv` | 128 | 目标 × 器官/平台/来源/资格分层 | 否（字节相同） |
| `diagnostic_sensitivity.tsv` | 9 | 3 诊断标记 × 3 相邻切点 | 否 |
| `cross_stage_consistency.json` | — | `compared=30030, max_abs_difference=0.0` | 否 |
| `REPORT.md` / `STATUS` | — | 描述性小结（含 §6.1 边界） | 数值更新 |
| `freeze_manifest.json` | — | 150 项；`human_tls_labels_read=false`、`biological_pass_threshold=null` | 更新 |

**修复的字节级影响范围**（以修复前 `freeze_manifest.json` 为对照，逐文件 sha256 复核）：148 项中 **138 项字节不变**，10 项改变——`unit_metrics.tsv`、`macro_summary.tsv`、`macro_strata.tsv`、`REPORT.md`、`checks/{check.log,checks.json,prediction_metrics_sample.tsv,prediction_module_relations_sample.tsv}`、`full_run_stages.log`、`summarize.log`；新增 2 项（`STATUS`、`summarize_fix.log`，均为清单自指时点差异）。**实测表、预测表、`splits/*`、单位内一致性与实测汇总均未被修复触及。**

## 3. 分账：臂、单位、来源、身份（修复后）

| 层 | 内部臂 | 来源臂 |
|---|---|---|
| 切片 | 72 | 53 |
| 既定单位（`unit_id`） | 65 | 47 |
| 每（目标×视图×模型）单位组 | **65（修复前 66）** | 47 |
| 划分 | internal_0..4（15/13/13/13/11 单位） | source_0..5（2/6/2/2/34/1） |
| 来源身份不足 | 19 片保留在实测清单与内部臂 | 不进入来源分母 |
| 器官外推 | 不适用 | 乳腺（source_5）、肝脏（source_0）留出时无同器官训练；臂级标签 `false`/`mixed` 各 1 480 行 |
| 平台 | Visium 252 + Visium HD 36（行数含 4 目标） | 同左 |
| 身份 | 确认 14 片 / 8 单位；其余 57 单位为未知，未升级为独立患者 | 同一身份表 |

- 单位级元数据冲突：`source_key`/`organ`/`platform`/`patient_known`/`confirmed_unit` 作为元数据保留，同一既定单位内不一致时标 `mixed` 并给 `*_n_distinct`；来源均衡与来源分层改用**划分级 `held_source`（`source_split_key`）**，不使用切片级 `source_key`。当前唯一 `mixed` 来自 `singleton::TENX13`（1 480 个分层行，仅内部臂），未被平铺成独立患者或研究身份。

## 4. 完整 / 部分 / 不可评 分母

逐（臂 × 模型 × 指标）表：`full_run_derived/denominators_by_arm_model_metric.tsv`（370 行，已按修复后单位表重算）。

- 恢复/关系类指标：内部 520 + 来源 376 单位组在 4 个非 `train_mean` 模型下全部可评；`train_mean` 按 `sd_floor` 结构性不可评。
- `partial_unit_metrics`：修复后 **28** 个单位行（修复前 26），**全部为 `geo_decay_half_pitch`**；合并后的 `singleton::TENX13` 因两片衰减可用性不同而正确暴露为 partial（不再被拆成两个 1 片单位而隐藏）。
- **`rec_decay_half_pitch_delta` 全库不可评**（内部 520 + 来源 376）：实测 `geo_decay_half_pitch` 仅 2/288 可评（多数为「首带 ≤0」或「band 上限内未衰减」）；空间尺度由 `geo_corr_band1..6` 与 `geo_patch75/90_*_radius_pitch` 承担。

## 5. 实测场资格与连续证据（72 片 × 4 目标，未受修复影响）

| 目标 | 任意目标基因检出点占比 中位（最小） | 有效基因数 中位（最小） | 检出基因数 中位 |
|---|---|---|---|
| B_AXIS | 0.999（0.612） | 2.87（1.19） | 47 |
| EMT | 1.000（0.962） | 32.79（7.66） | 199 |
| IFN-γ | 1.000（0.950） | 24.66（10.83） | 195 |
| OXPHOS | 1.000（0.919） | 75.21（4.20） | 191 |

实测状态：`DESCRIBED` 2、`DESCRIBED_WITH_NOT_ESTIMABLE_METRICS` 286、`NOT_TESTABLE_*` 0。诊断标记（非门槛）：首位基因占比 >0.5 共 21 场；有效基因数 <2 共 18 场。预测状态：`NOT_ESTIMABLE_CONSTANT_PREDICTION` 1 000（全部 `train_mean`）、`DESCRIBED_WITH_NOT_ESTIMABLE_METRICS` 3 994、`DESCRIBED` 6。

## 6. NOT_RUN / NOT_TESTABLE 全表

| 项 | 状态 | 依据 |
|---|---|---|
| 人工 TLS 结构验证、人工评价 | **NOT_RUN** | 全程未读取 |
| 独立 marker 隔离与重训 | **NOT_RUN** | 既有模型用完整目标基因集训练 |
| 训练/重训/调参/新模型 | **NOT_RUN** | 只重读 1000 份冻结预测 |
| 正式空间 p 值/显著性 | **NOT_RUN** | 打乱为未校准描述性扰动 |
| 复杂零模型、密度/支撑校准、校准不确定性、视野外/三维 | **NOT_RUN** | 合同 `not_run` |
| 来源身份不足 19 片的来源留出评价 | **NOT_TESTABLE_SOURCE_IDENTITY** | 保留在实测清单与内部臂 |
| 肺部来源留出（6 片肺） | **NOT_TESTABLE** | 全在来源身份不足集合内 |
| 乳腺/肝脏来源留出 | 器官外推标记 | 无同器官训练切片 |
| `rec_decay_half_pitch_delta` | **NOT_ESTIMABLE 全库** | 实测衰减尺度仅 2/288 可评 |
| `train_mean` 相关类指标 | **NOT_ESTIMABLE（结构性）** | 常量预测 + `sd_floor` |
| `geo_decay_half_pitch` 部分单位 | 部分可评（28 行） | 记 `partial_unit_metrics`，不进等权均值 |
| 物理尺度/微米/超分辨率 | **NOT_RUN，且不主张** | 仅以本片经验点间距为单位 |

## C. 聚合键修复（本轮）

- **代码**：`scripts/field_biological_review_20261019.py`
  - 单位身份键改为 `split, arm, target, view, model, model_class, cross_organ_extrapolation, unit_id`；`source_key`/`organ`/`platform`/`patient_known`/`confirmed_unit` 移出键，作为元数据聚合（一致取值 / `mixed` / `unknown`，并给 `*_n_distinct`）。
  - 新增单位列 `section_ids`、`source_split_key`（来源臂取划分级 `held_source`）；`macro_summary` 的来源均衡与 `macro_strata` 的 `source_split` 层改用 `source_split_key`，并保留 `source_metadata` 层呈现切片级来源元数据。
  - 未改任何指标公式、阈值、冻结预测读取路径或 `sd_floor`。
- **协议**（`docs/FIELD_BIOLOGICAL_REVIEW_PROTOCOL_20261019.md` §5、§7）：写明单位身份只由既定 `unit_id` 与分析键确定、元数据冲突标 `mixed`、来源均衡用 `source_split_key`。
- **定向检查**（`bash scripts/run_field_biological_review_20261019.sh check`，**1 969 项 / 0 失败**）：
  - 聚合契约 20 项，新增：冲突 `source_key` 必须合成单一单位组（`conflicting_source_key_single_unit_group`）、`n_sections=2`、`source_key=mixed`、`source_key_n_distinct=2`、单位值等于两片均值（0.1/0.3→0.2）、一致 `source_key` 保留原值、`source_split_key` 用于均衡、`source_split` 分层不出现 `mixed` 元数据值。
  - 真实样本（7 片，含 TENX13/TENX14 与 TENX171 外环孤立点）：`unit_groups_equal_established_unit_identity`（360 = 360）、`before_fix_key_would_split_units`（修复前键会给 400 组）、`real_two_section_unit_merged_mixed_source_key`（`singleton::TENX13`：2 片、`mixed`）。

## D. 前后影响（真实全量表核算）

完整表：`full_run_derived/unit_key_fix_impact.tsv`（臂 × 模型 × 指标的前后均值与单位数）；核验证据：`full_run_derived/unit_key_fix_checks.tsv`。

| 观察 | 修复前 | 修复后 |
|---|---|---|
| 单位行 | 4 520 | **4 480**（−40 = 内部臂 TENX13 伪行） |
| 内部臂每组合单位数 | 66 | **65**（来源臂 47 不变） |
| `singleton::TENX13` | 2 行 × 1 片，`source_key` 为两个标题 | **1 行 × 2 片**，`source_key=mixed`，`n_distinct=2`，`rec_rmse`=0.092598 = 两片均值 |
| partial 单位行（全为 `geo_decay_half_pitch`） | 26 | 28 |
| 分层行 / 带 reason 行 | 69 560 / 13 570 | 78 440 / 15 294 |
| 臂级宏平均（内部臂，`gene_ridge`） | — | `rec_rmse` 0.084655→0.084901（Δ+2.5e-04）；`rec_pearson` 0.773878→0.773665（Δ−2.1e-04）；`rec_jaccard_q75` 0.515458→0.515852（Δ+4.0e-04）；`rec_centroid_displacement_pitch` 0.757744→0.761098（Δ+3.4e-03） |
| 臂级宏平均（内部臂，`fixed_local_hgb`） | — | `rec_jaccard_q75` Δ+8.2e-04、`rec_centroid_displacement_pitch` Δ+4.9e-03、`rec_rmse` Δ+2.0e-04 |
| 来源臂全部指标 | — | **Δ = 0（完全不变）** |

即修复只改变内部臂中一个受影响单位的权重（相对量级 0.01%–0.2%），来源臂、实测表、预测表与全部逐片/逐模块明细不变。

## 7. 资源实际使用

- 磁盘：启动前 888.51 GB → summarize 时 888.4 GB → 修复后重跑 summarize 结束时 **888.2 GB**（净增约 0.3 GB；合同 `estimated_max_new_bytes = 3e9`）。全程高于 200,000,000,000 bytes。
- 时间：prepare 130 s + run 115 s + summarize 12 s + 修复后 summarize 14 s ≈ **4 分 31 秒**墙钟；检查（7 片版）另耗约 3 分钟。
- 输入核对：`input_checks.json` 记录 2 948 个冻结文件 + 72 个表达文件全部通过；白名单 3 738 条、命中拒读表 0。
- GPU：无（启动器强制 `CUDA_VISIBLE_DEVICES=-1`）。

## 8. Git 范围与未提交状态（D-265）

- `15dad63`：协议、脚本、启动器、合同、检查证据（含历史失败归档）、分母表（18 文件）。
- `6a89be3`：全量运行交接与台账（11 文件）。
- **本次修复提交**（提交信息“聚合键修复：既定 unit_id 不再因标题/source_key 重复赋权”）：更新后的脚本与协议、`checks/`（1 969 项证据）、`full_run_derived/{unit_key_fix_impact.tsv,unit_key_fix_checks.tsv,denominators_by_arm_model_metric.tsv}`、本文件，以及更新后的 `unit_metrics.tsv`、`macro_summary.tsv`、`macro_strata.tsv`、`REPORT.md`、`STATUS`、`freeze_manifest.json`。
- **仍未提交（工作区）**：`prediction_metrics.tsv`、`prediction_module_relations.tsv`、`measured_module_relations.tsv`、`splits/`、`smoke/` 大表、`full_run_derived/before_unit_key_fix/`（修复前聚合副本，供本地对照）与全部 `*.log`（被 `.gitignore` 排除）；路径与 sha256 见已提交的 `freeze_manifest.json`。
- 其他协作者的未提交改动（`AGENTS.md`、`STATUS.md`、`docs/*` 等）**未被本提交包含、未被修改**。
- 查看范围：`git show --stat 15dad63`、`git show --stat 6a89be3`、`git show --stat HEAD`；`git status --short`。

## 9. 失败与新发现

1. 全量三阶段无失败（返回码全 0）；修复后 `summarize` 重跑亦为 0。
2. 历史失败（已归档、未删除）：`checks/archive/` 下 7 次运行记录，其中 `run-13923425088ea8`、`run-08c20cf71466c8` 各 12 项失败（隔离点邻域语义、容差分层）；`checks/archive/smoke_fail_protocol_not_allowlisted.log`（`summarize` 读协议被允许表拒绝）。
3. 本次暴露并已修复：聚合键重复赋权（§C/§D）。修复后不再存在“同一既定单位因标题被拆成两个单位组”的情形；来源信息保留为元数据并在冲突时标 `mixed`，未推断患者或研究身份。
4. 记录项：`freeze_manifest.json` 中 `summarize.log`、`summarize_fix.log`、`full_run_stages.log` 等日志项在哈希写入后仍会继续追加，属自指时点差异；数据产物逐文件一致。

## 10. 边界（不替主代理判定）

本轮产物是既有开发材料的**描述性再分析**：无生物学成功/失败二分门槛、无正式空间 p 值、无独立 marker 验证、无 TLS 真值；预测与可见模块的相关含输入回声，不作为新关系；**数学上可算不等于有足够生物学参照**。本文件（含修复与前后影响）供独立 Grok 复核；复核完成前不宣布任何科学结论，也不替主代理判定 HF-01/HF-02/HF-04 的成立与否。

## 主代理记录勘误与复核结果（2026-10-19）

Grok对Git `5fee80f` 复核为 `ACCEPTED_FOR_DESCRIPTIVE_ANALYSIS_WITH_LIMITS`，无未关闭的结果相关问题，详见 `review/reviewed_unit_aggregation_fix.json`。根据审核更正本文件中内部衰减指标分母528→520，以及四目标有效基因数中位数；仅修正文档旧数，数据和代码未变。科学解释见 `docs/FIELD_BIOLOGICAL_REVIEW_RESULTS_20261019.md`；半衰尺度与自相关是细节分量指标，不代表原场的生物尺度。
