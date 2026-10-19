# 四个代理场体检（第一段）：执行交接（R-01 对照修复后）

**状态：第一段描述性表可用；匹配对照比较仅按目标作描述性使用（R-01 已关闭）。**
复核：`review/implementation_and_results.json`（R-01 阻断）→ 修复后由同一审核代理关闭 R-01，并给出使用边界：
对照对比可作描述性使用，但 **EMT 仅 2/72、OXPHOS 仅 7/72 可评，不能代表全体切片**；**不得恢复跨目标 +0.047，也不得恢复修复前“OXPHOS 低于对照”的说法**；片内同向仍不是生物协同、独立病理真值或 HF-01；**第二段仍不放行**。本文件所载数字未因该关闭而改动。

**声明上限：片内同向、留一相似、两半相关高于对照，都不是生物协同、独立病理真值或 HF-01。** 切片数不冒称患者数。

## R-01 修复（本轮唯一改动）

- **缺陷**：`match_cells` 对两半分别从同一目标外基因池抽取，同一次对照可能把同一基因同时放进两张对照图 → 对照相关被抬高（审核独立指出）。
- **修复**：新增 `match_cells_joint`——每次抽取在**一次抽样**中按“分档单元”取 `n_half_a + n_half_b` 个，前 `n_half_a` 给 A 半、其余给 B 半，**两半合计不重复**；若某档候选不足则整次对照不可评（不改档、不补抽）。对照种子独立记录（`control_seed`），逐次写入 `control_shared_genes`（应恒为 0）。
- **范围控制**：新增命令 `recompute-controls` 只重写 `subset_draws.tsv` 的对照列与 `subset_consistency.tsv` 的对照汇总列；目标两半的掩码用同一批种子重生成并**断言与存档的 `n_half_a`/`n_half_b` 完全一致**（不符即报错），因此目标集合的 10 次拆分未变。
- **未变证据**：与存档 `archive/pre_r01_controls/` 逐行比较，**非对照列 0 处变化**（draws 2 880 行、consistency 288 行）。
- **对照证据**：1 130 次可评对照抽取，`control_shared_genes` 最大值为 **0**；可评覆盖因规则更严而下降（见下），这正是修复生效的表现。
- **被取代的旧结论（不得再用）**：修复前“OXPHOS 低于对照（可评 31 片、中位差 −0.039）”与“跨目标 +0.047 / 64%”已作废，存档保留但不用于解释。

## 1. 实现与版本

- 计算脚本 `scripts/field_proxy_audit_20261019.py`（含 `recompute-controls`）；机器合同 `infra/field_proxy_audit_20261019/contract.json`（v2，逐条写入 D-277）；定向验证 `verify_results.py`。
- 旧代理读数沿用已审定义（基因图 = `log1p(1e4 × 原始计数 / 同点全基因文库)`；集合图为算术平均）；单基因/留一/两半图同计数、同分母、同组织点与条码顺序，仅换基因集合；零文库点不补值。
- 已审支持表复用 `infra/field_biological_review_20261019/measured_fields.tsv`（逐行记录 `approved_support_source`）；其有效基因数/首位份额按**原始计数总额**加权，本表“读出加权”列按 log1p 均值，二者分列、不混用。
- 面板外缺失与面板内零计数分列；未检出不是生物阴性。
- 可测基因 = 面板内 + ≥20 个组织点非零 + 方差>0；互斥子集各需 ≥3（合计 ≥6），否则 `NOT_ESTIMABLE`（本轮 0 例）。
- 技术对照 = **每张可测基因图一次片内回归**（截距 + 排除本目标基因后的点级文库 log1p + 非目标检出比例）；原图与残差分开，残差不替代原图。
- 固定对照 = 每（切片×目标）**10 次固定种子**均分并全部逐次入库；非目标对照按检出比例×平均标准化表达的四分位档匹配，**单次两半不重复**。

## 2. 覆盖与分母

| 项 | 值 |
|---|---|
| 切片 × 目标 | **72 × 4 = 288** 行 |
| 单位 | **65**（每目标 65 行，共 260 单位×目标行） |
| 臂 | 内部 72 片/65 单位；来源 53 片/47 单位（6 个可辨来源）；来源不足 19 片仅内部 |
| 身份 | 确认 14 片 / 8 单位；其余 57 单位未知 |
| 器官外推 | `source_0`（肝）、`source_5`（乳腺）；肺来源不可评 |
| 单位×指标状态 | COMPLETE 4 044、INCOMPLETE_UNIT 4、NOT_ESTIMABLE_UNIT（对照不可评）随修复后覆盖变化 |
| 结构化不可评价 | `not_testable.tsv` **0** 行（无切片×目标低于 ≥6 可测基因） |
| **对照可评覆盖（修复后）** | 288 个 case 中 **175 个无可评对照抽取**；B_AXIS 66/72、IFN-γ 38/72、**EMT 2/72**、**OXPHOS 7/72** |

## 3. 主要描述性结果（详见 `REPORT.md`）

**代理支持（72 片中位）**

| 目标 | 面板内 | 任意检出 | 可测（≥20 点且方差>0） | 已审有效基因数 | 已审首位份额 | 读出加权首位份额 |
|---|---:|---:|---:|---:|---:|---:|
| B_AXIS | 47 | 47 | 38 | 2.87 | 0.535 | 0.236 |
| IFN-γ | 196 | 195 | 188 | 24.66 | 0.120 | 0.044 |
| EMT | 200 | 199 | 190 | 32.79 | 0.100 | 0.033 |
| OXPHOS | 191 | 191 | 191 | 75.21 | 0.055 | 0.022 |

**聚合代理图自身的技术可解释方差（`proxy_technical_r2`，片内一次回归的 R²）**：B_AXIS 中位 **0.11**、EMT **0.17**、IFN-γ **0.43**、OXPHOS **0.76**（0.40–0.93）。逐基因技术 R² 中位仅 0.007–0.051——**两者不可混为“基因都被技术梯度主导”**；聚合图的高 R² 不能外推到逐基因，残差变小也不能单独判定生物或技术来源。

**片内基因对协同（中位相关，原图 → 残差图）**：B_AXIS 0.037→0.032；IFN-γ 0.022→0.009；EMT 0.037→0.030；OXPHOS 0.063→0.011。

**删除最大贡献基因**：B_AXIS 被删份额 0.236，剩余图与可测均值相关 0.985、q75 Jaccard 0.869、质心位移 0.34 pitch；三个 Hallmark 删除份额 0.022–0.044，相关 ≈1.000、Jaccard 0.976–0.983。注意 B_AXIS 的 0.985 是 **log1p 均值图** 的留一稳定性，**不推翻已审的计数加权有效基因数 2.87 / 首位份额 0.535**。

**互斥两半 vs 匹配对照（修复后；仅按目标单独读，不跨目标合并）**

| 目标 | 两半相关中位 | 可评对照切片 | 可比切片里观测更高 | 观测减对照中位 |
|---|---:|---:|---:|---:|
| B_AXIS | 0.599 | 66/72 | 92% | **+0.264** |
| IFN-γ | 0.790 | 38/72 | 82% | +0.050 |
| EMT | 0.916 | **2/72** | 100% | +0.164（仅 2 片） |
| OXPHOS | 0.873 | **7/72** | 57% | +0.007 |

EMT 的对照中位（0.68）只来自 2 片、OXPHOS（0.84）只来自 7 片；**不能读成全体切片**。B_AXIS 是唯一覆盖较高（66/72）且差值较大的目标，但它同时是已审支持最窄的目标（有效基因数 2.87）。

## 4. 解释边界

- 全部为描述性：无 PASS 阈值、无显著性、无独立病理真值；高值集合是操作阈值（严格 > q75）。
- **R-01 已关闭（复核代理确认），不再修改对照代码**：对照差可按目标作描述性使用；EMT（2/72）与 OXPHOS（7/72）的可评覆盖很窄，不得读作全体切片；不恢复跨目标 +0.047，不恢复修复前“OXPHOS 低于对照”的说法。
- B_AXIS 的差值在内部与来源臂一致，但只说明“同集合的基因图在同一切片更一致”，不是独立分子身份或场间关系证据。
- OXPHOS/IFN-γ 的聚合代理图对非目标深度/检出高度敏感（R² 0.76/0.43），其代理图解释必须与技术成分同读；B_AXIS 支持窄。
- **第一段不锁定第二段**：不提供候选基因成员、权重或候选数；第二段未授权。

## 5. 定向验证（`verify_results.py`，22 项 / 0 失败）

覆盖 288 行/72 片/4 目标、每目标 65 单位、单位状态键合法；冻结读数复现 **≤6.0e-08**；已审支持表 288/288 连接且 `approved_genes_detected`/`in_panel` 逐行相等、来源逐行记录；删除表/子集表各 288 行；**每 case 恰 10 次抽取且 10 个不同种子**、每半 ≥3 基因；**`control_shared_genes` 全部为 0（1 130 次可评对照）**；对照可评次数 0–10 全部记录；`not_testable` 原因集合合法（空）；条码与冻结 cache 相等；每张图与 `base.target_readout` 一致；残差对两协变量正交 4.8e-16；技术协变量确实排除目标基因（重建差 0.0）。
另：`recompute-controls` 与存档比较**非对照列 0 处变化**（draws 2 880 行、consistency 288 行）。

## 6. 命令、计时与资源

```bash
cd /home/huyudi/013_spatial
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=-1
PY=/opt/anaconda3/bin/python
$PY scripts/field_proxy_audit_20261019.py prepare              # 冻结输入哈希核对
$PY scripts/field_proxy_audit_20261019.py audit                # 72 片，4 workers
$PY scripts/field_proxy_audit_20261019.py summarize
$PY scripts/field_proxy_audit_20261019.py recompute-controls   # R-01：只重算对照列
$PY scripts/field_proxy_audit_20261019.py summarize            # 单位/分层/报告重算
$PY infra/field_proxy_audit_20261019/verify_results.py
```
| 阶段 | 退出码 | 墙钟 | 进程 User | 进程 Sys | 最大 RSS |
|---|---:|---:|---:|---:|---:|
| prepare | 0 | 80 s | 34.4 s | 8.8 s | 165 MB |
| audit | 0 | 78 s | 237.8 s | 60.0 s | 2.52 GB |
| summarize | 0 | 3 s | 2.6 s | 0.2 s | 206 MB |
| recompute-controls | 0 | 29 s | — | — | — |
| verify | 0 | ~5 s | — | — | — |

CPU、本地数据；人工 TLS / annotations / incoming / HookNet 未读取；产物约 22 MB，磁盘余量 ≥200 GB。

## 7. 交付物

- 结果：`REPORT.md`、`section_metrics.tsv`(288)、`synergy.tsv`(288)、`gene_coherence.tsv`(43 014)、`deletion_sensitivity.tsv`(288)、`subset_consistency.tsv`(288)、`subset_draws.tsv`(2 880)、`not_testable.tsv`(0)、`unit_metrics.tsv`(260)、`unit_macro.tsv`、`unit_strata.tsv`。
- 核对：`verify_results.json`（22 项）、`input_checks.json`、`prepare_receipt.json`、`recompute_controls.log`、`stages.log`、`*_time.txt`。
- R-01 证据：`archive/pre_r01_controls/`（修复前的 draws/consistency/unit/REPORT 副本，仅作对照，不用于解释）。
- 未提交：`sections/`（逐片中间目录）、`smoke/`（2 片贯通检查）。
- 第二段：未执行、未授权。
