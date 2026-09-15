# R-16 设计方案 v2.1（2026-09-15，D-118 尺子工厂版）

依据 D-112/D-113/D-114/D-115/D-118。本文档冻结参数；改动需追加 decision。

**v2.1 变更（D-118）：** Tier-1 新增数据驱动分支 v2 工厂——邻域增强（λ∈{0,0.2,0.5,0.8}）→逐片逐λ PCA top-10 候选轴→LISA tile 相干（999 次值置换＋BH-FDR）→新旧对照（|r|≥0.7 KNOWN）→跨病人 loading |cos| 平均链接匹配→复发按有 tile 结构病人数。结论为可信阴性：NEW 复现轴为 0。

**v2.0 变更（D-115）：** 普查单层结构 → **两层结构**：Tier-1 连续组成轴尺子厂（地基）＋Tier-2 残余 cluster 发现。第一性原理修正：硬聚类把"程序"与"混合邻居"焊死（L1）；离散合法但旋钮必须披露可扫（L2 用户修正）；跨患者证据可用连续共享度（L3）；患者独有 pattern 升一等公民但注明可能污染（L4）；Tier-1 轴已实证（L5）。

## 1. 范围与架构

```
Tier-1（地基，已实证）：marker 模块轴 → 轴分场 → 相干检验 → 等高线几何
Tier-1 v2（尺子工厂，已执行，结论阴性）：邻域增强×λ → PCA 候选轴 → LISA tile → 新旧对照 → loading 匹配
Tier-2（发现）：对不上轴的 cluster 普查（平均链接＋纯度＋轴对齐）
Lane A（下一轮）：gene~f(轴分) 或 gene~f(到等高线有符号跳数)
Lane B（下一轮）：回归掉全部轴/距离协变量后的残余空间结构筛查
```

首轮交付：Tier-1 全部＋Tier-2 全部＋registry v1（已完成 2026-09-14）。Lane A/B 待 registry 过目后启动。

## 2. 输入

- **发现集**：Vanderbilt training 47 片/30 患者（单 lineage 单平台）。
- **复现引擎（首轮不触碰）**：ST-CRC 14 片/7 患者、USZ 8 片/8 患者（带 TLS 注释）；registry 行升级时在复现引擎检验，跨平台按 D-113 降半级。
- **gene universe**：整合 HVG top-10000（D-114）。cache：`infra/r16/census_cache_hvg10k/`（本地不入 git，manifest 入 git）。
- 坐标图跳数，标称 ×100µm/跳标记近似。

## 3. Tier-1：组成轴尺子厂（已实证，`r16/axes.py`＋`scripts/r16_axis_rulers.py`）

| 参数 | 取值 |
|---|---|
| 轴 | **工作集六轴（D-117）**：T/B/Mye/Epi/Stromal/Plasma。ILC 退役（稀少、定义不稳、实证最弱）。「+x」仅浆细胞已测；内皮/增殖/缺氧未测，加入须另批。 |
| 轴分 | 模块基因 z 分 log1p 均值（片内 z-score） |
| 空间相干 | Moran's I＋值置换 null 200 次；axis 级判据：≥2 患者片中位 p≤0.01 |
| 等高线几何 | 分位 0.6/0.7/0.8 三档，相邻 Jaccard 入 registry notes |
| 阳性对照 | CXCL13 对 B 轴 q0.7 等高线的有符号跳数谱；判据：多数患者峰在 ≤+1 跳 |

实测（`infra/r16/axis_rulers_20260914.json`）：相干 27–30/30 患者；Jaccard ≈0.71；CXCL13 峰 ≤+1 跳：25/30（峰分布 −4…+5，主峰深部）。

## 3b. Tier-1 v2：卷积×tile 尺子工厂（D-118，已执行，`r16/axis_factory.py`＋`scripts/r16_axis_factory.py`）

| 参数 | 取值 |
|---|---|
| 邻域 | 空间 kNN k=15 行归一化均值（自身排除），对 log1p 求均值再片内 z-score；k=30 只做半径复核 |
| 混合 | Z=√(1−λ)·X_z＋√λ·M_z，λ∈{0,0.2,0.5,0.8} 全档披露 |
| 候选轴 | 每片每λ PCA top-10（seed 固定）；loading 即基因权重 |
| tile | 逐点 LISA＋BH-FDR（α=0.05）；热点=连通分支≥20 点；**抽样 999 次**（200 次的置换地板在 BH 下已证伪） |
| 新旧对照 | 候选场与六轴场 max\|Pearson r\|≥0.7 记 KNOWN，否则 NEW |
| 匹配 | loading \|cos\| 平均链接，切 0.25（=cos 0.75）；复发=有 tile 结构成员的病人数；纯度<0.5 强制 DESCRIPTIVE |
| 升级 | 复发≥2 病人 → EXPLORATORY_REPRODUCED；registry tier `1b`，kind `discovered_axis` |

实测（`infra/r16/axis_factory_20260915.json`，47 片×4λ×10PC=1880 候选）：pilot 证实 λ 为已知→新奇旋钮（tile_ok 22→41；known 比率 0.18→0.10）；半径复核 k15-vs-k30 中位 0.975 通过；全量 857 组仅 2 组复现且均为 KNOWN（R16B_001 持家/Epi 程序 29 病人；R16B_002 基质程序 3 病人），**NEW 复现轴=0**；跨病人 loading 相似度中位 0.306，0.75→0.5 全程无阈值膝点；1299 个 NEW 且有 tile 结构的候选均匀散在 30 病人（新奇不少，共享为零）。registry v1.1 共 1085 行。

## 4. Tier-2：残余 cluster 普查（`scripts/r16_composition_census.py`）

| 参数 | 取值 |
|---|---|
| 聚类 | 每片 log1p→z→PCA30→kNN15→Leiden 0.25/0.5/1.0（seed 20260914），min 20 spots |
| 匹配 | **平均链接层次聚类**（cosine 距离，切 0.25=cos 0.75）——杀单链接链式雪崩 |
| 纯度 | 组内最小两两余弦；**<0.5 → LOW_PURITY_MIXED_GROUP，强制 DESCRIPTIVE**（毛发团如实登记） |
| 轴对齐 | 组均值签名在 7 轴上的得分；best_axis ≥0.5 记 `axis_explained`，否则 `UNEXPLAINED`（D-112 残余维度的内建测量） |
| 跨档稳定 | 相邻档 Jaccard ≥0.5；resolution_support 入 registry |
| 空间相干 | mask Moran's I＋值置换 null 200 次 |
| 普遍性 | **按患者计数**；患者独有组一等公民，等级注明可能污染（L4） |

实测（`infra/r16/composition_pattern_census_20260914.json`）：995 cluster → 221 组；2 个毛发团（纯度 0.19/0.28，已强制 DESCRIPTIVE）；13 个干净多患者组（EXPLORATORY_REPRODUCED）。

## 5. Registry（`infra/r16/field_registry.tsv`，schema r16.field_registry.v1）

列：`pattern_id, tier(1|2), kind, axis_name, lineage, n_sections, n_patients, resolution_support, top_markers, geometry, effect_size, null_type, null_draws, null_p, fdr_q, cross_patient_shape_corr, residual_dimension, uncertainty_ci95, evidence_grade, notes`。
等级：`DESCRIPTIVE_SINGLE_PATIENT | EXPLORATORY_REPRODUCED | CLAIM_CANDIDATE`；CLAIM_CANDIDATE 需另行用户批准；复现引擎检验是升级前提。
v1：228 行（7 Tier-1＋221 Tier-2；20 行 EXPLORATORY_REPRODUCED）。

## 6. 零假设规格（沿用 v1.1 §6，I-021 细化版）

自相关检验（Tier-1 相干、Tier-2 mask、Lane B）→ 值置换合法；关联检验（Lane A 信号~轴/距离）→ **空间保持 null**（60° 旋转＋凸包内平移＋重映射，≥95% spot 落 0.5 跳内接受）——Lane A 开工前须落地并验收（I-021 仍对 Lane A 阻塞）。抽样 200 次。

## 7. Lane A / Lane B（下一轮，本版不执行）

- Lane A 双接法：gene~f(轴分)（组成坐标，无需几何）与 gene~f(到等高线跳数)（几何坐标）；函数族三类（等调单调/样条 df=4/壳层单峰）；升级判据 ≥2 患者同族同向＋形状相关 ≥0.5＋效应地板 0.5 SD。
- Lane B：回归掉全部轴＋距离协变量后的残余，Moran＋值置换；跨患者一致性按患者。
- 工程：每片预计算旋转/平移索引映射全基因复用。

## 8. 计算预算与产出

已完成：materialize（47 片/10k 基因）、Tier-1（CPU 分钟级）、Tier-2（CPU ~10 分钟）、registry v1、Tier-1 v2 工厂全量（CPU ~5 分钟）＋半径复核、registry v1.1（1085 行）。产出物全部 `EXPLORATORY_*_NOT_CLAIM`。

## 9. 已知限制（随结果一起交付）

单 lineage；坐标跳数近似；毛发球吞弱特异性 pattern（L-006 marker-overlap 候选解）；ILC 轴弱（12 基因，相干 27/30）；轴集合覆盖受 marker 库限制（D-105 失效条件沿用）；mask 视为固定（不确定性传播留敏感性检查）；**PC 匹配只做符号对齐，旋转混合歧义未解——同一程序在不同片拆分到不同 PC 上无法匹配，复发检验对此类程序灵敏度不足（见 LEADS L-008）**。
