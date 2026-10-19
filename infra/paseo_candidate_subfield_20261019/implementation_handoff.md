# 两分量候选场探针（D-279）实施与结果交接

- 执行代理：Pi `opencode-go/deepseek-v4.1-flash`（thinking `max`），cwd `/home/huyudi/013_spatial`
- 版本：**Git `054a16b`**（代码+产物）；口径写死在 `f46e60d`/`df20880`，权威科学合同为 `4732584`（D-279）
- 放行依据：`infra/paseo_candidate_subfield_20261019/review/preflight.json` = `FIT_RELEASED`（`computation_authorized=true`）；`run` 阶段有释放门强制校验
- 复核状态：**已复核，无必改项、不需重跑**（Grok `high`，`infra/paseo_candidate_subfield_20261019/review/implementation_and_results.json`，范围 `054a16b`）
- 结论上限：**描述性结果可用，但不能称候选场已被核验**；不是独立病理真值、独立生物场或 HF-01；未使用 2026-10-06 FC-03 分数/因子/K=8/0.5 线

## 1 按 D-279 实现的口径

| 环节 | 实现（代码位置） |
|---|---|
| A/B 拆分 | `split_genes`：基因名排序 + `default_rng(0).permutation`，A=ceil(2n/3)、余数归 B；B_AXIS 33/16、三个 Hallmark 134/66；读表达前固定 |
| 分母 | `cache_section`：`L_external = 点总计数 − 本目标整个面板计数`；A、B 图统一 `log1p(1e4·count/L_external)`；`L_external<=0` 的点从该目标两侧全部剔除 |
| 可测 A | `assemble`：训练单位内**每一片** ≥20 非零点且单位方差>0 才算该单位可测；≥半数训练单位才进入 A 名单；A<6 → 四臂 `NOT_TESTABLE_INSUFFICIENT_A` |
| 抽样 | `n_sample=min(512, 各训练单位最小可用点数)`，每单位同数，双切片按片均分，种子 0；测试单位全点、不重抽 |
| 缩放 | 每基因训练抽样**正值均值**，同一缩放施于测试 |
| 模型 | `fit_nmf`：K=1/K=2，`nndsvda/cd/random_state=0/max_iter=500/tol=1e-4`；未收敛/零载荷/非有限 → `NOT_TESTABLE` |
| 对齐 | 折内 5 次 90% 行重抽样 vs 完整训练拟合，按共享 A 载荷两种一一配对取余弦和大者（平局原序）；跨折仅描述（`crossfold_similarity.py`） |
| 参照 | A 均值、K1、单基因（**未缩放**训练图均值最高）、全 A Ridge(alpha=10,StandardScaler 训练期) —— 同折同 A 同抽样 |
| B | 仅作目标值；测试片 ≥20 非零且方差>0 才可评；q75 严格 `>`，空集/常量不可评 |
| 汇总 | 逐基因→逐片→逐单位（任一片不可评则该基因该单位剔除）→完整单位等权；`summary_*` + 分层 `summary_strata_main.tsv` |
| 技术 | `depth=log1p(点总计数−A∪B 面板计数)`、`detection=面板外阳性比例`，单独报告 |

## 2 复现命令

```bash
python scripts/field_candidate_subfield_20261019.py cache --workers 4   # 72 片面板缓存（本机约 3 min，113 MB，不入库）
python scripts/field_candidate_subfield_20261019.py smoke              # 单目标单折贯通（B_AXIS × internal_0）
python scripts/field_candidate_subfield_20261019.py run --workers 3    # 四目标 × 11 划分（放行门校验）
python scripts/field_candidate_subfield_20261019.py summarize
python infra/field_candidate_subfield_20261019/strata_summary.py main
python infra/field_candidate_subfield_20261019/crossfold_similarity.py
python infra/field_candidate_subfield_20261019/build_report.py
python scripts/field_candidate_verify.py main                          # 473 项检查
```

## 3 运行结果

- 44/44 单元 `OK`；K1、K2 全部收敛；可用 A：B_AXIS 26–30、EMT 128–133、IFN 127–129、OXPHOS 128–134；`n_sample` 每单元 512。
- 逐片 133 750 行：`OK` 91.7%，`NOT_TESTABLE_SPARSE_OR_CONSTANT` 7 475、`NOT_TESTABLE_GENE_ABSENT` 2 075、`NOT_ESTIMABLE_TIES_OR_CONSTANT` 980、`NOT_TESTABLE_TRAIN_ABSENT` 625。
- 逐单位 119 840 行，其中不可评 10 580；完整单位汇总 11 488 行。

### 各目标（完整单位等权，B 基因逐单位汇总，Pearson 中位）

| 目标 | A 均值 | K1 | 单基因 | K2 | 全 A Ridge |
|---|---|---|---|---|---|
| B_AXIS | **0.131** | 0.105 | 0.051 | 0.106 | **0.140** |
| EMT | 0.283 | 0.283 | 0.214 | **0.294** | 0.286 |
| IFN-γ | 0.127 | 0.128 | 0.070 | 0.136 | **0.146** |
| OXPHOS | 0.261 | 0.264 | 0.055 | 0.256 | 0.230 |

### K=2 是否带来增益（同目标同划分同 B 基因配对）

| 目标 | K2 − A 均值（中位差 / K2 更高比例） | K2 − 全 A Ridge |
|---|---|---|
| B_AXIS | −0.014 / 30% | −0.018 / 25% |
| EMT | +0.011 / 70% | +0.003 / 53% |
| IFN-γ | +0.002 / 55% | −0.008 / 40% |
| OXPHOS | −0.001 / 45% | +0.022 / 87% |

### 稳定性与技术

- 两分量载荷余弦中位：B_AXIS 0.27、EMT 0.51、IFN 0.64、OXPHOS 0.65（越高越接近冗余；**不设阈值、不据此命名两个场**）。
- 折内重抽样配对余弦 ≈ 0.99998（组成在同一折内可重复）；跨折 220 对全部可评，最佳配对余弦中位 0.971–0.994。
- 分量与技术协变量相关中位（|r|）：B_AXIS 0.15–0.23、EMT 0.19–0.61、IFN 0.35–0.48、OXPHOS 0.55–0.71；OXPHOS 的技术相关最高，其 K2 对全 A Ridge 的表面优势**可能部分是技术共变**，按合同不得据此单独裁决。

## 4 定向验证（`verification_main.tsv`：473 项，0 失败）

覆盖：A∩B=∅ 与 ceil 规则（4 目标）；冻结划分复用、单位不跨训练/测试、来源留出确持出（11 划分）；`L_external` 与 A/B 图按原始 h5ad 重算逐点一致（8 组）；候选特征列只含 A 基因（44 单元）；每片 ≥20 与半数训练单位可测规则重算一致（44 单元）；`n_sample` 同数与测试不重抽（44 单元）；缩放=训练正值均值（132 例）；测试投影用训练冻结载荷（NNLS 复核，44 单元）；q75 严格大于；宏观均值等于完整单位等权（非基因×片直接平均）；`run` 释放门在未放行时必须中止；旧 FC-03 记号检查。

## 5 复核范围建议（Grok）

1. 泄漏：B 是否进入 A 图分母、候选输入、训练标准化、测试投影、分量匹配或选模（`scripts/field_candidate_subfield_20261019.py` 的 `cache_section/assemble/run_split`）。
2. 分组：`infra/field_source_holdout_20261018/*/split.json` 是否原样复用、单位是否不跨训练/测试（`field_candidate_verify.py` S02 组）。
3. 归一化/条码：`L_external` 定义与 `L<=0` 剔除、测试不重抽、q75 严格（S02/S05 组）。
4. 计权：逐基因→逐片→逐单位→完整单位等权是否与 `summary_per_unit.tsv`/`summary_macro.tsv` 一致（S06 组）。
5. 不可评语义：稀疏/缺失/并列是否只记 `NOT_TESTABLE` 而未当阴性（`summary_per_section.tsv` 的 status 列）。
6. 口径边界：报告是否只称候选场、是否避免 FC-03 与 HF-01 表述（`REPORT.md`）。

已知限制：技术协变量可能带走真实生物变化，本轮只并列报告；B 基因在训练片缺失（625 行）与测试片稀疏（7 475 行）使部分基因无法评价；OXPHOS 的分量技术相关最高，其相对优势需谨慎解读。

## 6 复核结论（Grok `high`，2026-10-19，范围 `054a16b`）

状态 `DESCRIPTIVE_RESULTS_USABLE_NO_CANDIDATE_ADVANTAGE`：`allow_descriptive_report=true`、`allow_scientific_claim=false`、`required_fix=false`。已关闭项：A/B 33/16 与 134/66 且交集为空、A/B 图用排除整面板的外部文库、抽查可用 A 不含 B、Ridge 特征只来自 A 且标准化/投影都在训练期冻结、折内配对只用 A 载荷、跨折脚本不读测试图或测试 B、44 单元全部可运行且单位不跨训练测试、任一片不可评即剔除该基因该单位、全量验证按 11 划分×4 目标而非 smoke 单单元期望。

独立重算的中位差与其一致：B_AXIS −0.0141（K2 更高 29.7%）、EMT +0.0111（70.0%）、IFN-γ +0.0017（54.9%）、OXPHOS −0.0007（44.5%）；对全 A Ridge：B_AXIS −0.0183、EMT +0.0029、IFN −0.0082、OXPHOS +0.0222（K2 更高 86.9%）。

使用边界（按复核原文）：

- 固定 K=2 **没有稳定优于同一批 A 基因的简单均值**；全 A Ridge 在 B_AXIS、EMT、IFN-γ 上至少不弱于 K=2。
- OXPHOS 的 +0.0222 与该目标**最高的深度/检出技术相关同时出现**，不能写成"更干净的分子身份"。
- 折内重抽样接近 1 只说明同一次拟合可重复，**不是跨单位发现了场**；两分量载荷余弦仅报告，不命名两个场。
- 内部与来源分开后方向不变；乳腺 source_5、肝脏 source_0 仍是器官外推，不当同器官复现；可用 A 的 26–30 是半数训练单位规则，不是旧可行性表的 16–18。
- 该结论限制的是**这次固定的两分量表示**，不是所有场发现，也不支持 HF-01；旧 FC-03 的氧化磷酸化优势仍不计入。
- 复核不承诺无缺陷（`defect_free_guarantee=false`）。代码与产物保持冻结，不再修改。
