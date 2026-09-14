# R-16 设计方案 v1（2026-09-14，待用户批准）

依据 D-112/D-113。本文档冻结参数；批准后除"校准条款"（§9）外不再改动，改动需追加 decision。

## 1. 范围

本方案覆盖 R-16 全管线（普查→距离尺→Lane A→Lane B→registry），但**首轮只执行第 ① 步（组成 pattern 普查）＋registry 骨架**。②–⑤ 在 ① 的产出经用户过目后再启动。

## 2. 输入

- 训练 manifest 47 行（`infra/r04/role_manifests/training_manifest.json`），30 患者，**全部 HTAN_VANDERBILT_CRC 单 lineage 单平台**——本方案不产任何跨队列/跨平台结论。
- counts：panel cache（4000 基因，CSR npz＋barcode json）。
- 坐标：源 h5ad `obsm/spatial`，按 barcode 精确对齐（已验证：obs `_index` 与 cache barcode 一致，spot 数一致）。
- 坐标单位不可靠（中位最近邻间距 ≈44.8，非 µm）——**距离一律用图跳数**，标称换算 ×100µm/跳，标记为近似值。

## 3. 第①步：组成 pattern 普查（首轮执行）

| 参数 | 取值 | 理由 |
|---|---|---|
| 预处理 | log1p → 片内逐基因标准化 → PCA 30 维（sklearn） | 与 R-04 惯例一致 |
| 聚类 | Leiden（leidenalg），表达 kNN 图 k=15 | 空间聚集作为**评分**而非约束——聚不聚交给 null 判 |
| 分辨率档位 | 0.25 / 0.5 / 1.0，三档，seed=20260914 | 跨档稳定性入分 |
| 最小 cluster | ≥20 spots（665–5000 spot 片的 0.4%–3%） | 更小即 Visium 尺度噪音 |
| cluster 签名 | 全 4000 基因 z-score 均值谱（匹配用）＋ top-10 t 统计 marker（registry 备注用） | 匹配与解释分离 |
| 跨片匹配 | 签名余弦 ≥ 0.75，连通便为 pattern 组（单 lineage 内） | 阈值按 §9 校准 |
| 跨档稳定性 | 相邻档位 Jaccard ≥ 0.5 视为同一 cluster 延续；pattern 组记 3 档中获支持档数 | 防分辨率采樱桃 |
| 空间一致性 | mask 二值场 Moran's I，**空间保持 null**（§6），200 次 | 禁用标签置换（I-021） |
| 普遍性 | 拥有匹配 cluster 的**患者数**（非片数） | a29 教训内建 |

产出：pattern 组清单（每组：患者数、片数、档位支持、top marker、Moran's I 及 null p），全部进 registry，`evidence_grade` 最高只到 `EXPLORATORY_REPRODUCED`。

## 4. 第②步：距离尺（下一轮）

- 图：坐标 cKDTree k=6（Visium 六边形拓扑），BFS 跳数距离。
- 几何：**有符号到边界距离**（mask 内为负，主用）；到质心无符号距离（辅，标记）。跳数封顶 ±8。
- 分辨率地板写入登记行：直径 <3 spot 的 pattern 不支持壳层类检验；<2 跳的梯度判混叠。

## 5. 第③④步：Lane A / Lane B（下一轮）

- 函数族（预选，三类）：(a) 等调单调（升/降）；(b) 自然样条 df=4；(c) 壳层（单峰位置自由，两侧单调）。
- Lane A 单元 = pattern × 基因 × 族：跳数分箱谱拟合加权 R²；显著性 = 空间保持 null 200 次经验 p；片内 BH-FDR q<0.05 仅为修剪；**升级判据 = ≥2 患者同族同向且形状相关 ≥0.5**，效应地板 = 分箱谱极差 ≥0.5 片内 SD。
- Lane B：逐片 ridge 回归掉本组全部 pattern 指示变量 → 残余场 Moran's I＋空间 null → 跨患者一致性。**这是 D-112 的残余评分维度，不是门禁。**
- 关键工程：每片预计算 200 个旋转/平移索引映射，全部基因复用（单基因单次 null ≈ 50µs，整屏 CPU 分钟–小时级）。

## 6. 空间保持 null 规格（I-021 换代件）

1. **主用**：场相对组织 mask 做 k×60° 旋转（绕质心）＋凸包内随机平移，最近 spot 重映射；≥95% spot 落在 0.5 跳内才接受，否则重抽。保持自相关、只破坏对齐。
2. **抽检**（v2，非首轮）：经验变异函数匹配 GP 模拟，验证主用 null 的 p 值校准。
3. 禁用：naive 值置换、标签置换（连续信号场景）。
4. 抽样次数 200（项目惯例，p 分辨率 0.005）。

## 7. Registry schema（`r16.field_registry.v1`）

每行：`pattern_id, kind(composition_mask|distance_function|residual_structure), lane(A|B), lineage, n_sections, n_patients, resolution_support(k/3), top_markers[10], geometry, effect_size, null_type, null_draws, null_p, fdr_q, cross_patient_shape_corr, residual_dimension, uncertainty_ci95, evidence_grade(DESCRIPTIVE_SINGLE_PATIENT|EXPLORATORY_REPRODUCED|CLAIM_CANDIDATE), notes`。CLAIM_CANDIDATE 需另行用户批准，本管线不自动升级。

## 8. 计算预算与产出

- ①：47 片 × 3 档 Leiden + 200 null × 复用索引 → CPU 十分钟级；无 GPU、无新数据、无外部请求。
- 产出：`infra/r16/composition_pattern_census_20260914.json`（+TSV）＋ `infra/r16/field_registry.tsv` 骨架＋测试。
- 验收：所有 47 片处理成功或明确记失败原因；registry 每行可追溯到片/患者/参数。

## 9. 校准条款（允许不回批的调整）

首轮跑前用 5 片 pilot 校准：余弦匹配阈值（看分布取膝点，报告原值与新值）、最小 cluster（若档位产出全 <20 则降到 10）。**只许收紧不许放宽**；调整值写回本文档 v1.1 备注。

## 10. 已知限制（随结果一起交付）

单 lineage 无跨队列结论；panel 盲区（I-020）；坐标跳数近似；mask 视为固定（不确定性传播按 D-113 留作敏感性检查）；发现即登记、登记≠claim。
