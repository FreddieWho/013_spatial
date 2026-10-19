# 四个代理场体检（第一段，描述性）

只描述代理集合的内部一致性。**不把片内同向写成生物协同、独立病理真值或 HF-01**；不命名真实场；高值集合是操作阈值（严格大于 q75），不是生物真值。

覆盖：72 片 × 4 目标；单位 65（跨目标）；内部 72 片/65 单位、来源 53 片/47 单位；来源不足 19 片仅内部；肺不可作来源；乳腺/肝脏为器官外推。

## 1. 代理支持与可测基因

| 目标 | 片数 | 面板内中位 | 任意检出中位（已审口径） | 可测基因中位（≥20 点且方差>0） | 已审有效基因数中位 | 已审首位份额中位 | 读出加权首位份额中位 | 原读数复算最大差 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| B_AXIS | 72 | 47 | 47 | 38 | 2.87 | 0.535 | 0.236 | 5.9e-08 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | 72 | 196 | 195 | 188 | 24.66 | 0.120 | 0.044 | 5.9e-08 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 72 | 200 | 199 | 190 | 32.79 | 0.100 | 0.033 | 6.0e-08 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | 72 | 191 | 191 | 191 | 75.21 | 0.055 | 0.022 | 6.0e-08 |

上述“已审”列直接取自 `infra/field_biological_review_20261019/measured_fields.tsv`（来源逐行记录在 `approved_support_source`），其有效基因数与首位份额按原始计数总额加权；“读出加权”列按本表统一读出（log1p 均值）计算，两者定义不同、分列报告，不混用。

## 2. 片内基因空间协同（基因对相关分布；原图 vs 单次片内回归残差）

| 目标 | 基因对中位 | 原图中位 | 原图 q10–q90 | 残差中位 | 残差 q10–q90 | 技术可解释方差中位 |
|---|---:|---:|---|---:|---|---:|
| B_AXIS | 684 | 0.030 | -0.019–0.146 | 0.024 | -0.020–0.135 | 0.007 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | 17578 | 0.022 | -0.019–0.095 | 0.009 | -0.024–0.061 | 0.017 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 18050 | 0.037 | -0.030–0.205 | 0.030 | -0.027–0.175 | 0.017 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | 18145 | 0.051 | -0.002–0.129 | 0.009 | -0.024–0.053 | 0.051 |

技术对照：每张可测基因图**一次**片内回归（截距 + 排除本目标基因后的点级文库 log1p + 非目标检出比例）。该回归也会去掉与细胞密度耦合的生物信号，残差相似或消失均不能单独判定生物/技术来源；共同的全基因归一化分母仍可能引入依赖。

## 3. 删除最大贡献基因（可测集合内）

| 目标 | 片数 | 被删份额中位 | 与可测均值相关中位 | q75 Jaccard 中位 | 质心位移中位(pitch) | 剩余基因中位 |
|---|---:|---:|---:|---:|---:|---:|
| B_AXIS | 72 | 0.236 | 0.985 | 0.869 | 0.340 | 36 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | 72 | 0.044 | 1.000 | 0.976 | 0.038 | 187 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 72 | 0.033 | 1.000 | 0.983 | 0.045 | 190 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | 72 | 0.022 | 1.000 | 0.980 | 0.021 | 190 |

与“原代理全定义图”的相关另记在 `deletion_sensitivity.tsv`（`corr_vs_proxy_full_definition`）：基因集合不同、公式相同，单独报告。

## 4. 互斥子集一致性（10 次固定种子全部报告 + 四分位档匹配的非目标对照）

| 目标 | 片数 | 两半相关中位 | q10–q90（片间） | q75 Jaccard 中位 | 匹配对照相关中位 | 可评对照次数中位 |
|---|---:|---:|---|---:|---:|---:|
| B_AXIS | 72 | 0.599 | 0.498–0.643 | 0.392 | 0.286 | 10.0 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | 72 | 0.790 | 0.771–0.804 | 0.497 | 0.724 | 10.0 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 72 | 0.916 | 0.894–0.925 | 0.676 | 0.682 | 0.0 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | 72 | 0.873 | 0.870–0.881 | 0.576 | 0.836 | 0.0 |

**对照可评覆盖（必须与上表同读，且不得跨目标合并）**

| 目标 | 切片×目标 | 有可评对照 | 无可评对照 | 可比切片里观测更高 | 观测减对照中位 |
|---|---:|---:|---:|---:|---:|
| B_AXIS | 72 | 66 | 6 | 92% | +0.264 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | 72 | 38 | 34 | 82% | +0.050 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 72 | 2 | 70 | 100% | +0.164 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | 72 | 7 | 65 | 57% | +0.007 |

全体 288 个切片×目标中 175 个没有任何可评对照抽取（四分位档在本片目标外可测基因中候选不足，逐次记录在 `subset_draws.tsv`）。**每个目标单独读；不给跨目标合并的差，也不把该差当显著性、生物参照或 HF-01。**单次对照内两半不重复取同一基因（`control_shared_genes` 全为 0）。

每个互斥子集至少 3 个可测基因（合计 ≥6），否则该片该目标为 NOT_ESTIMABLE（见 `not_testable.tsv`）；对照档位候选不足只让**该次抽取**不可评，逐次计数保存在 `subset_draws.tsv`，不得改档位或补抽。

## 5. 单位与分母（完整单位等权；任一片不可评则该单位该指标不进均值）

- n_measurable：完整单位 260/260，单位均值 149.033
- effective_gene_count：完整单位 260/260，单位均值 69.052
- top_gene_share：完整单位 260/260，单位均值 0.090
- proxy_technical_r2：完整单位 260/260，单位均值 0.380
- pair_raw.median：完整单位 260/260，单位均值 0.044
- pair_raw.q10：完整单位 260/260，单位均值 -0.021
- pair_raw.q90：完整单位 260/260，单位均值 0.163
- pair_residual.median：完整单位 260/260，单位均值 0.024
- pair_residual.q10：完整单位 260/260，单位均值 -0.026
- pair_residual.q90：完整单位 260/260，单位均值 0.121
- technical_r2.median：完整单位 260/260，单位均值 0.031
- deletion_corr_vs_measurable_mean：完整单位 260/260，单位均值 0.993
- deletion_jaccard_q75：完整单位 260/260，单位均值 0.944
- subset_corr_median：完整单位 260/260，单位均值 0.775
- subset_jaccard_q75_median：完整单位 260/260，单位均值 0.534
- control_corr_median：完整单位 100/260，单位均值 0.485

### 6 个可辨来源（仅 ELIGIBLE；乳腺/肝脏为器官外推）

| 来源 | 目标 | 单位数 | 完整 | 两半相关中位 | 匹配对照中位 |
|---|---|---:|---:|---:|---:|
| https://www.biorxiv.org/content/10.1101/2022.11.09.5（器官外推） | B_AXIS | 2 | 2 | 0.460 | 0.201 |
| https://www.biorxiv.org/content/10.1101/2022.11.09.5（器官外推） | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 2 | 2 | 0.857 | NA |
| https://www.biorxiv.org/content/10.1101/2022.11.09.5（器官外推） | HALLMARK_INTERFERON_GAMMA_RESPONSE | 2 | 2 | 0.541 | 0.601 |
| https://www.biorxiv.org/content/10.1101/2022.11.09.5（器官外推） | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 2 | 2 | 0.816 | NA |
| https://www.biorxiv.org/content/10.1101/2023.01.23.5 | B_AXIS | 6 | 6 | 0.609 | 0.269 |
| https://www.biorxiv.org/content/10.1101/2023.01.23.5 | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 6 | 6 | 0.911 | NA |
| https://www.biorxiv.org/content/10.1101/2023.01.23.5 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 6 | 6 | 0.732 | 0.665 |
| https://www.biorxiv.org/content/10.1101/2023.01.23.5 | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 6 | 6 | 0.820 | NA |
| https://www.biorxiv.org/content/10.1101/2024.06.04.5 | B_AXIS | 2 | 2 | 0.723 | 0.316 |
| https://www.biorxiv.org/content/10.1101/2024.06.04.5 | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 2 | 2 | 0.927 | NA |
| https://www.biorxiv.org/content/10.1101/2024.06.04.5 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 2 | 2 | 0.828 | 0.665 |
| https://www.biorxiv.org/content/10.1101/2024.06.04.5 | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 2 | 2 | 0.915 | NA |
| https://www.cell.com/cell/fulltext/S0092-8674(21)013 | B_AXIS | 2 | 2 | 0.565 | 0.433 |
| https://www.cell.com/cell/fulltext/S0092-8674(21)013 | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 2 | 2 | 0.952 | NA |
| https://www.cell.com/cell/fulltext/S0092-8674(21)013 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 2 | 2 | 0.772 | 0.825 |
| https://www.cell.com/cell/fulltext/S0092-8674(21)013 | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 2 | 2 | 0.880 | NA |
| https://www.cell.com/cell/fulltext/S0092-8674(23)012 | B_AXIS | 34 | 34 | 0.544 | 0.295 |
| https://www.cell.com/cell/fulltext/S0092-8674(23)012 | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 34 | 34 | 0.904 | 0.682 |
| https://www.cell.com/cell/fulltext/S0092-8674(23)012 | HALLMARK_INTERFERON_GAMMA_RESPONSE | 34 | 34 | 0.779 | 0.740 |
| https://www.cell.com/cell/fulltext/S0092-8674(23)012 | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 34 | 34 | 0.892 | 0.886 |
| https://www.nature.com/articles/s41467-023-43458-x（器官外推） | B_AXIS | 1 | 1 | 0.856 | 0.371 |
| https://www.nature.com/articles/s41467-023-43458-x（器官外推） | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 1 | 1 | 0.947 | NA |
| https://www.nature.com/articles/s41467-023-43458-x（器官外推） | HALLMARK_INTERFERON_GAMMA_RESPONSE | 1 | 1 | 0.859 | NA |
| https://www.nature.com/articles/s41467-023-43458-x（器官外推） | HALLMARK_OXIDATIVE_PHOSPHORYLATION | 1 | 1 | 0.916 | NA |

标题型 `source_key`（`title:...`）属来源身份不足的 19 片，**只留在内部描述**，不得读作第 7 个来源；肺切片全部在其中，故肺来源不可评。

完整分母与分层见 `unit_metrics.tsv`（`__n`、`__state`）、`unit_macro.tsv`、`unit_strata.tsv`；逐片见 `section_metrics.tsv`。

## 6. 不可评价记录

- 共 0 条：{}
- 结构化不可评价（不是生物学阴性）：无 panel 基因；可测基因 <3（协同/留一）；<6（两半各需 ≥3）；常量图；技术残差退化；对照档位候选不足（逐次记录）；条码/坐标不一致（硬失败）。

## 7. 边界与第二段

- 描述性：无 PASS 阈值、无显著性、无独立病理或生物真值；片内同向 ≠ 生物协同 ≠ HF-01。
- 两半互斥、另配匹配对照，避免同一批基因既构图又自证。
- **第一段不选定第二段的基因成员、权重或候选数**；第二段未授权，本报告不提供拆分结论。

