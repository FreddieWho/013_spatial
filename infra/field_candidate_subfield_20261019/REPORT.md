# 两分量候选场探针（D-279）：结果报告

## 口径与结论上限

- 四目标全报；候选为逐折固定 K=2 的两分量非负分解，A/B 按基因名与种子 0 拆分。
- A、B 图共用“排除本目标整个面板”的外部文库分母；训练单位冻结可测 A、等量抽样与缩放；测试 B 仅事后评价。
- 普通参照（A 均值、K=1、单基因、全 A Ridge）与候选同折同 A 同训练抽样；全 A Ridge 是更强参照。
- **即使稳定、非冗余、受 B 支持，也只称候选场**：不是独立病理真值，不是独立生物场，不是 HF-01；
  旧 2026-10-06 FC-03 分数与氧化磷酸化优势不计入本轮。

## 运行状态

- 目标×划分单元：44（期望 44）
- 状态分布：{'OK': 44}

## 各目标：候选与参照（完整单位等权，B 基因逐单位汇总）

### B_AXIS

| 臂 | 可评单位数（范围） | 基因×划分行 | Pearson 中位 | Jaccard 中位 | Pearson>0 比例 |
|---|---|---|---|---|---|
| A_mean | 1–34 | 172 | 0.131 | 0.105 | 0.98 |
| K1 | 1–34 | 172 | 0.105 | 0.101 | 0.96 |
| single_gene | 1–34 | 171 | 0.051 | 0.088 | 0.92 |
| K2 | 1–34 | 172 | 0.106 | 0.109 | 0.97 |
| full_A_ridge | 1–34 | 172 | 0.140 | 0.110 | 0.97 |

### HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION

| 臂 | 可评单位数（范围） | 基因×划分行 | Pearson 中位 | Jaccard 中位 | Pearson>0 比例 |
|---|---|---|---|---|---|
| A_mean | 1–34 | 721 | 0.283 | 0.230 | 0.91 |
| K1 | 1–34 | 721 | 0.283 | 0.228 | 0.90 |
| single_gene | 1–34 | 721 | 0.214 | 0.205 | 0.88 |
| K2 | 1–34 | 721 | 0.294 | 0.235 | 0.95 |
| full_A_ridge | 1–34 | 721 | 0.286 | 0.235 | 0.98 |

### HALLMARK_INTERFERON_GAMMA_RESPONSE

| 臂 | 可评单位数（范围） | 基因×划分行 | Pearson 中位 | Jaccard 中位 | Pearson>0 比例 |
|---|---|---|---|---|---|
| A_mean | 1–34 | 705 | 0.127 | 0.165 | 0.97 |
| K1 | 1–34 | 705 | 0.128 | 0.166 | 0.98 |
| single_gene | 1–34 | 705 | 0.070 | 0.133 | 0.86 |
| K2 | 1–34 | 705 | 0.136 | 0.162 | 0.96 |
| full_A_ridge | 1–34 | 705 | 0.146 | 0.165 | 0.98 |

### HALLMARK_OXIDATIVE_PHOSPHORYLATION

| 臂 | 可评单位数（范围） | 基因×划分行 | Pearson 中位 | Jaccard 中位 | Pearson>0 比例 |
|---|---|---|---|---|---|
| A_mean | 1–34 | 712 | 0.261 | 0.193 | 0.99 |
| K1 | 1–34 | 712 | 0.264 | 0.194 | 0.99 |
| single_gene | 1–34 | 651 | 0.055 | 0.132 | 0.68 |
| K2 | 1–34 | 712 | 0.256 | 0.193 | 0.99 |
| full_A_ridge | 1–34 | 712 | 0.230 | 0.190 | 1.00 |

## K=2 相对参照的配对差（同目标同划分同 B 基因）

| 目标 | 对比 | 配对数 | 中位差 | K2 更高比例 |
|---|---|---|---|---|
| B_AXIS | K2 − A_mean | 172 | -0.014 | 0.30 |
| B_AXIS | K2 − full_A_ridge | 172 | -0.018 | 0.25 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | K2 − A_mean | 721 | 0.011 | 0.70 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | K2 − full_A_ridge | 721 | 0.003 | 0.53 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | K2 − A_mean | 705 | 0.002 | 0.55 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | K2 − full_A_ridge | 705 | -0.008 | 0.40 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | K2 − A_mean | 712 | -0.001 | 0.45 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | K2 − full_A_ridge | 712 | 0.022 | 0.87 |

## 内部与来源留出（完整单位等权，Pearson 中位）

| 目标 | 臂 | 内部 5 折 | 来源 6 屏 |
|---|---|---|---|
| B_AXIS | A_mean | 0.134 | 0.129 |
| B_AXIS | K1 | 0.106 | 0.104 |
| B_AXIS | K2 | 0.106 | 0.106 |
| B_AXIS | full_A_ridge | 0.153 | 0.129 |
| B_AXIS | single_gene | 0.050 | 0.053 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | A_mean | 0.271 | 0.306 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | K1 | 0.273 | 0.310 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | K2 | 0.280 | 0.326 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full_A_ridge | 0.276 | 0.303 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | single_gene | 0.198 | 0.248 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | A_mean | 0.129 | 0.125 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | K1 | 0.128 | 0.129 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | K2 | 0.138 | 0.133 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | full_A_ridge | 0.153 | 0.139 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | single_gene | 0.075 | 0.065 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | A_mean | 0.267 | 0.256 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | K1 | 0.268 | 0.257 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | K2 | 0.265 | 0.248 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | full_A_ridge | 0.242 | 0.215 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | single_gene | 0.066 | 0.029 |

## 分量稳定性与技术相关

| 目标 | 折内重抽样配对余弦（中位，按分量均） | 两分量载荷余弦（中位） |
|---|---|---|
| B_AXIS | 1.000 | 0.273 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 1.000 | 0.510 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | 1.000 | 0.644 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | 1.000 | 0.650 |

| 目标 | 分量 | 深度相关（中位） | 检出率相关（中位） |
|---|---|---|---|
| B_AXIS | 0 | 0.152 | 0.160 |
| B_AXIS | 1 | 0.227 | 0.205 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 0 | 0.188 | 0.200 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 1 | 0.608 | 0.587 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | 0 | 0.472 | 0.476 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | 1 | 0.364 | 0.349 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | 0 | 0.562 | 0.546 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | 1 | 0.703 | 0.712 |

## 失败与不可评

- 逐片行 133750，逐单位行 119840，完整单位汇总行 11488。
- 不可评（单位级）：10580。
- 稀疏/并列导致的逐片不可评见 `summary_per_section.tsv` 的 `status` 列。

## 边界

- 片内同向与 B 支持不等于生物协同、独立病理真值或 HF-01；片数不是患者数。
- 本探针不检验场间空间恢复；隐藏场恢复需在候选冻结后另做，仍保留同点与整片背景对照。
- 人工 TLS 标签、旧人工评价、HookNet 全轮封存。

