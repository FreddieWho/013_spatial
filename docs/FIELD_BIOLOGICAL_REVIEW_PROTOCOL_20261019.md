# 实测场资格与冻结预测恢复：生物学复核协议（2026-10-19）

**性质：** 描述性复核协议。不含训练、不含新模型选择、不含 PASS/FAIL 生物学阈值、不产生 p 值。
**依据：** `infra/paseo_field_review_20261019/implementation_task.md`（主代理任务书）与用户 2026-10-19 补充边界。
**执行合同：** `infra/field_biological_review_20261019/contract.json`（机器可读、先于运行冻结）。
**实现：** `scripts/field_biological_review_20261019.py`。
**输出根目录：** `infra/field_biological_review_20261019/`。

本文件只描述**已冻结输入**上的复核口径。它不修改 `docs/plan.md`，不新增 HF 编号，不改变既有冻结脚本、模型、清单或结果。

---

## 1. 要回答什么

在本地已经冻结的四目标预测（B 细胞轴、氧化磷酸化、IFN-γ 响应、EMT）之上，形成一份**描述性清单**：

1. **分子检测覆盖与多基因支持**：每个切片的每个目标场，其分子读出由多少基因实际支撑、检测稀疏到什么程度（连续报告，不设通过线）。
2. **富集位置、范围、空间尺度**：实测场的高值在组织内落在哪里、是否成片、在多大尺度上衰减；单位一律为**该切片自身的经验点间距**（pitch），不换算微米、不主张分辨率。
3. **与可见场的相对分布关系**：实测场与同一视图内各可见模块场的相关、共定位与位置/尺度差异；以及冻结预测是否复现了这组关系（而非只复现整体水平）。

同时记录**哪些格子不可评价**以及原因。不可评价只来自"该统计量在数学上无定义"，不来自任何生物学百分比门槛。

**本阶段不做的判断：** 不判定"某场恢复成功/失败"的二分结论，不给关系排序筛选，不声称发现独立生物场，不声称 TLS 结构成立或不成立。

---

## 2. 输入边界

### 2.1 允许读取

| 类别 | 具体来源 |
|---|---|
| 冻结预测 | `infra/field_source_holdout_20261018/{internal_0..4, source_0..5}/*.npz`（1000 份包，5 个模型 × 2 个视图） |
| 冻结划分与身份 | `infra/field_source_holdout_20261018/{splits.json, input_manifest.json, contract.json, definitions.json, vocabulary.json}` |
| 冻结预测误差表（仅交叉核对） | `infra/field_source_holdout_20261018/section_metrics.tsv` |
| 冻结几何与实测读出 | `infra/field_discovery_20261018/cache/{sid}.npz`、`{sid}_near.npz`、`{sid}_outer.npz` |
| 冻结可见基因矩阵 | `infra/field_source_holdout_20261018/genes/`（small5）、`infra/field_followup_20261018/genes/`（full15） |
| 表达与基因身份 | `data/other_sources/hest1k/st/{sid}.h5ad` 的 `X`、`var/_index`、`obs/_index`、`obs/in_tissue` |
| 模块定义 | `definitions.json`（按 target × view 的模块→基因表）、`vocabulary.json` |
| 代码 | `scripts/field_discovery_20261018.py`（作为库复用）、本脚本 |

### 2.2 禁止读取（硬边界）

- 人工 TLS 标签、`annotations/`、`incoming/`、任何历史人工标签评价产物；含 `tls_annotation_pack`、`tls_label_field_20261001`、`tls_four_cancer_labels_20261008`、`zenodo_st_crc_cms/Pathology_SpotAnnotations`、HookNet 输出与权重。上述子串写在合同的 `denied_path_substrings`，脚本安装第二道审计钩子：任何 `open()` 路径命中即抛错；同时启动时断言**白名单内没有任何路径命中拒读表**。
- 外部网络下载、外部联系、任何新数据获取。
- 非本协议白名单内的仓库文件；实现用 `InputGuard`（仓库内 Python `open` 审计钩子 + HDF5 显式字段白名单）拦截。**该白名单是 Python 级访问审计，不冒称操作系统沙箱。**
- metadata 的 `disease_state`/`oncotree_code` 不用于自动合并 Cancer 与 Tumor；文件名不补患者或诊断身份（本阶段不使用这些字段）。

### 2.3 资源与数值分辨率

- 仅 CPU（`CUDA_VISIBLE_DEVICES` 非空即拒绝启动），线程数受合同限制。
- 运行前/中检查可用空间 ≥ 200,000,000,000 bytes，另需合同声明的 `estimated_max_new_bytes` 余量。
- **方差下限 `sd_floor = 1e-6`（计算规则，非生物阈值）**：冻结场以 float32 存储，低于该下限的起伏无法与表示噪声区分。任一方 sd ≤ floor 时，相关、自相关、连通片、质心与尺度指标为 `NOT_ESTIMABLE`（`constant_field` 或 `degenerate_variance`），**不填 0、不填 1**；RMSE/MAE/bias 仍保留为诊断。已知事实：来源留出 1000 份 `train_mean` 包中 364 份 sd=0、其余约 3e-8，均低于下限；抽查切片的实测读出 sd 约 1e-2。逐场 `sd_band` 列报告 `zero` / `below_floor` / `measurable`。

---

### 2.4 三本账（不得共用“成功”一词）

| 账本 | `evidence_role` | 内容 | 不得当作 |
|---|---|---|---|
| 测量资格 | `measurement_eligibility` | 目标基因是否在面板、检出覆盖、基因支撑、实测场分布/几何 | 恢复成功 |
| 输入回声关系 | `input_echo_relation_not_spatial_increment` | 场与**已进入该模型输入**的可见模块场的相关/共定位；模块基因与目标基因已断言不相交 | 独立生物关系或空间增量 |
| 探索性再分析 | `exploratory_reanalysis_of_developed_predictions` | 只重读已存预测并重算指标，不是重训、不是独立标志基因验证、不是独立患者确认 | 新证据 |

三本账写在不同表且带 `evidence_role` 列；预测行另有 `cross_budget_verdict='none'`，本阶段不产生任何“哪个模型赢”的字段。

### 2.5 信息预算标注（不做胜负）

预测行带 `model_class`、`uses_neighbourhood_input`、`is_posthoc_smoothing`：

| 模型 | 视图 | `model_class` | 邻域输入 | 后处理平滑 |
|---|---|---|---|---|
| `gene_ridge` | small5/full15 | `raw_gene_ridge` | 否 | 否 |
| `gene_ridge_smooth_half` | small5/full15 | `raw_gene_ridge_posthoc_neighbourhood_smoothing` | 是（后处理） | 是 |
| `fixed_point_hgb` | small5/full15 | `module_point_hgb` | 否 | 否 |
| `fixed_local_hgb` | small5/full15 | `module_local_hgb` | 是 | 否 |
| `train_mean` | 任意 | `constant_prior` | 否 | 否 |

small5 与 full15 改变基因数与可见归一化分母；因此**跨预算只描述，不判胜负，也不宣布空间增量**。

## 3. 记号与几何约定

- 切片 `s`：组织内点 `n_s` 个（`obs/in_tissue == 1`），坐标为源像素 `xy`（`cache/{sid}.npz: xy_source_pixels`）。
- **经验点间距 pitch** `h_s` = 坐标最近邻距离的中位数（源像素）。所有长度除以 `h_s` 后以 "pitch" 为单位报告。**不报告微米、不主张超分辨率**。
- **冻结邻域图**：`near`（≤ 1.1 h_s 的行归一化权重图）与 `outer`（1.1–3.1 h_s），直接读取 `cache/{sid}_near.npz` / `_outer.npz`，与冻结模型使用的邻域完全一致。
- **邻域均值** `nm(v) = W_near v`（孤立点取自身值，与冻结代码一致）。
- **细节分量** `d(v) = v − nm(v)`：扣除邻域平均后的局部组分，用来区分"整体水平"和"切片内空间位置"。
- **位置权重** `w(v) = v − min(v)`：非负、对整体平移不变；`w ≡ 0`（常量场）时位置类指标不可评价。
- **模块场**：模块 `m` 在视图 `V` 下的场 = 该模块在 `definitions[target][V][m]` 中的基因、取自冻结矩阵 `genes/{sid}_t{ti}.npy` 的列均值（该矩阵已是 log1p(10⁴×counts/library) 归一化读数）。该定义与冻结模型的 `point` 设计逐列一致，并在检查中核对：**模块场使用与冻结模型相同的可见信息预算，目标基因全体已从可见词汇中剔除，且脚本在启动时断言可见词汇与目标基因集不相交**，因此不存在目标泄漏。
- **实测场** `y`：`cache/{sid}.npz: y[:, ti]`，即冻结目标读出（目标基因集内可用基因的 log1p 归一化计数均值）。复核从 `X` 重算并与冻结值逐点比对，最大绝对差写入 `input_checks.json`。**本复核只读取已存预测并重算描述量，不是独立训练重跑，也不产生新的模型比较。**
- **目标基因来源**：目标基因表来自冻结面板文件（`panel_controls.json`，哈希与冻结 discovery 收据一致）；`definitions.json`/`vocabulary.json` 只描述**可见**模块与可见基因并集，不作为目标定义。面板外缺失（`sup_genes_absent_from_panel`）与面板内零计数（`sup_genes_in_panel_zero_counts`）分开报告，不合成同一状态。
- **可见基因矩阵的来源核对**：`genes/{sid}_t{ti}.npy` 的列必须与 `vocabulary.json` 一一对应、行必须与条码顺序对应。检查阶段在真实切片上用 `field_followup.gene_input` 从 `X` 独立重算并与冻结矩阵逐元素比较（容差 1e-5），同时用一个行置换对照确认行序不对时能被检出。QC 列（`*_qc.npy`）在本阶段**不使用**。

**支撑与边界：** 评价对象只有 `obs/in_tissue == 1` 的组织条码；不把阵列上非组织点补进分母；边缘点与孤立点**保留**、不丢弃也不补零，邻域均值对孤立点返回自身值（与冻结 `neighbor_mean` 完全一致，并在检查中用稠密独立实现比对）；坐标重复即报错（不静默去重）；条码重复即报错。逐片报告 `geo_near_degree_min/median/max` 与 `geo_near_isolated_spots`，使边界/空洞影响可见。

**距离尺度的锁定：** 坐标列固定为 `obs/pxl_col_in_fullres` / `obs/pxl_row_in_fullres`，经冻结 cache 存为 `xy_source_pixels`；pitch = 组织内点的**最近邻距离中位数**（重复坐标已报错排除），重算值必须与冻结 manifest 的 `empirical_pitch_source_pixels` 相等（容差 1e-9）；邻域半径固定 1.1 / 3.1 pitch，**直接复用冻结 `{sid}_near.npz`/`{sid}_outer.npz`**，因此与模型邻域同一支撑，不存在 `NOT_COMPARABLE` 情形。不换算微米。

---

## 4. 指标

三组指标对**实测场与每个冻结预测场使用同一实现**，保证可比。所有计算用 float64。

### 4.1 分子检测覆盖与多基因支持（`sup_*`，仅实测场）

从原始计数（`X` 中目标基因列）计算，不用读出分数：

| 指标 | 生物学含义 | 量程 |
|---|---|---|
| `sup_genes_defined` | 目标定义基因数 | ≥1 |
| `sup_genes_in_panel` | 本切片基因面板中存在的定义基因数 | 0…defined |
| `sup_genes_detected` | 至少 1 个组织内点计数非零的基因数 | 0…in_panel |
| `sup_det_frac_min/q25/median/q75/max` | 单基因检出点占比的分位数：该基因在多少比例的点上被检出 | 0–1 |
| `sup_spots_with_any_gene_frac` | 至少检出 1 个目标基因的点占比 | 0–1 |
| `sup_genes_per_spot_median/q90` | 每点检出的目标基因数分位数 | 0…detected |
| `sup_effective_gene_count` | 表达占比的参与比 `(Σs_g)²/Σs_g²`：读出实际由"多少个基因等效支撑" | 1…detected |
| `sup_top_gene_share` / `sup_top3_gene_share` | 单一/前三位基因占总读出表达占比 | 0–1 |
| `sup_target_count_share_median` | 目标基因计数占该点总库中位比例 | 0–1 |
| `sup_mean_if_detected` / `sup_sd_if_detected` | 仅在"检出任意目标基因"点上的均值/标准差（读出稀疏度的敏感性描述） | 实数 |
| `sup_readout_max_abs_diff` | 复核读出与冻结 `y` 的最大绝对差（完整性检查，非指标） | ≥0 |

**全部分子支持量连续报告。** 本阶段不使用"检出率 ≥ X%"之类的通过线；低覆盖不触发 NOT_TESTABLE，只作为数值与诊断标记一起呈现（见 §6）。

### 4.2 分布与几何（`dist_*`、`geo_*`，实测与预测同口径）

| 指标 | 生物学含义 | 量程 / 单位 | 不可评价条件 |
|---|---|---|---|
| `dist_mean/sd/min/max/q05/q25/median/q75/q95` | 场的取值分布（读出单位） | 实数 | — |
| `dist_frac_zero` | 读出恰为 0 的点占比；**表示该点未检出可支持转录本，不表示生物学阴性** | 0–1 | — |
| `dist_frac_negative` | 预测为负的点占比（读出自构造非负，负值为模型伪影） | 0–1 | — |
| `dist_skew` | 偏度 | 实数 | sd=0 |
| `geo_centroid_shift_pitch` | 位置权重质心到几何中心（等权质心）的距离 / pitch：高值偏在哪一侧 | ≥0，pitch | 常量场、h=0 |
| `geo_centroid_box_x/y` | 位置权重质心在组织外接矩形内的归一位置 | 0–1 | 同上 |
| `geo_centroid_hull_radial` | 位置权重质心到组织凸包边界的最短距离 / (0.5·√凸包面积)：0=贴边，≈1=居中 | ≥0 | 凸包退化（共线/点过少） |
| `geo_autocorr_near` | 细节分量 `dv = v − nm(v)` 在冻结近邻图上的空间自相关指数 `Σ_ij W_ij·dv_i·dv_j / Σ_i dv_i²`：相邻点是否同向偏离邻域均值 | ≈[−1,1] | 常量场或 `dv ≡ 0` |
| `geo_autocorr_outer` | 同上，1.1–3.1 pitch 环带 | 同上 | 常量场 |
| `geo_corr_band1..6` | 距离分带指数 `Σ_pairs dv_i·dv_j / (\|pairs\|·mean(dv²))`，`dv = v − nm(v)`；带边见合同（单位 pitch），配 `geo_band_k_mean_dist_pitch` 报告该带实际平均距离 | ≈[−1,1] | 常量场 / 该带无点对 |
| `geo_decay_half_pitch` | 相关衰减到首带一半时的距离（线性内插，单位 pitch）：空间尺度描述 | ≥0，pitch | 首带 ≤0，或 10 pitch 内未衰减 |
| `geo_patch75_n` | 阈值 q75 下高值点在冻结近邻图上的连通片数。阈值定义为**严格大于**该场 q75 分位（`v > quantile(v,0.75)`）：零膨胀或大量并列值时集合自动退为"高于该值"的点，而不是被并列值填满 | ≥0 | 常量场；严格集合为空（`empty_threshold_set`） |
| `geo_patch75_largest_share` | 最大连通片占高值点比例：集中还是分散 | 0–1 | 常量场 |
| `geo_patch75_largest_radius_pitch` / `mean_radius_pitch` | 连通片半径（片内点到片质心的最大距离 / pitch） | ≥0，pitch | 常量场 |
| `geo_patch75_area_frac` | 严格高值点占比（连续场≈0.25；零膨胀场即非零点占比，逐场读出实际值） | 0–1 | — |
| `geo_patch75_at_threshold_frac` | 恰好等于阈值分位的点占比：并列质量大小，用于判断该阈值是否有区分度 | 0–1 | — |
| `geo_patch90_*` | 同上，阈值 q90 | 同上 | 同上 |
| `geo_patch75_*_shuffled_median/p10/p90` | 同一切片内把该场取值随机重排后（保持取值边际）同名统计量的中位与分位 | 同上 | 常量场 |

**分位阈值的性质：** q75/q80/q90 是**操作阈值**，用于把连续场切成可比较的点集，**不是生物真值、不是阳性定义**。阈值取**严格大于**分位值，避免零膨胀字段（大量未检出点读出恰为 0）被并列值填满；每场同时报告 `area_frac` 与 `at_threshold_frac`，使实际集合大小与并列程度都可追溯（不同场之间点集大小因此**可能不同**，比较的是空间组织而不是阳性数量）。协议同时报告 q75/q80/q90 以显示阈值敏感性；打乱参考按合同仅在 q75 计算（成本参数，非科学门槛）。

### 4.3 与可见场的相对分布关系（`rel_*`，实测与预测同口径）

对每个（切片，目标，视图，模块）组合：

| 指标 | 含义 | 量程 |
|---|---|---|
| `rel_pearson` / `rel_spearman` | 场与该模块场逐点相关（线性 / 秩） | [−1,1] |
| `rel_jaccard_q75/q80/q90` | 场与模块场在操作阈值下的高值点交集 Jaccard | 0–1 |
| `rel_jaccard_q75_shuffled_median/p10/p90` | 保持场的取值边际、打乱其空间位置后与模块场的 Jaccard 分布（**描述性扰动参考，非零模型、无 p 值**） | 0–1 |
| `rel_centroid_distance_pitch` | 场与模块场位置权重质心距离 / pitch | ≥0，pitch |
| `rel_module_n_genes` / `rel_module_mean` / `rel_module_sd` / `rel_module_autocorr_near` | 模块场自身属性（基因数、均值、标准差、近邻自相关） | — |

**注意：** 预测场是可见模块（及可见基因）的函数，因此"预测场与可见模块的相关"部分由模型输入结构决定，不能读作独立生物学证据。本表只描述"预测是否复现了实测场与可见场的关系图样"，该限制写入结果报告。

**关系轮廓（profile）：** 每个（切片，目标，视图，场）对应一个长度 = 可见模块数的相关向量。定义：

| 指标 | 含义 | 不可评价条件 |
|---|---|---|
| `prof_relation_pearson` | 实测关系向量与预测关系向量在模块间的相关：预测是否按同样的模块排序 | 任一向量常量；模块数 < 3 |
| `prof_relation_l2` | 两向量欧氏距离 | 模块不可用 |
| `prof_abs_delta_mean` | 逐模块相关差的绝对值均值 | 模块不可用 |
| `prof_top_module_agreement` | `argmax|corr|` 的模块是否一致（并列即不可评价） | 并列/常量 |

### 4.4 恢复指标（`rec_*`，实测 vs 冻结预测）

| 指标 | 含义 | 量程 | 不可评价条件 |
|---|---|---|---|
| `rec_rmse` | 均方根误差（读出单位） | ≥0 | — |
| `rec_centered_rmse` | 去均值后 RMSE：只比空间图样，不受整体水平影响 | ≥0 | — |
| `rec_mae` | 平均绝对误差 | ≥0 | — |
| `rec_bias` | 预测均值 − 实测均值 | 实数 | — |
| `rec_sd_ratio` | 预测标准差 / 实测标准差；<1 表示压平（过度平滑的常见标志） | ≥0 | 实测 sd=0 |
| `rec_calibration_slope` | 以预测为自变量、实测为因变量的 OLS 斜率；1 为理想 | 实数 | 预测常量化 |
| `rec_pearson` / `rec_spearman` | 逐点相关 | [−1,1] | 任一方常量 |
| `rec_local_pattern_pearson` | 细节分量 `d(y)` 与 `d(p)` 的相关：是否恢复了扣除邻域均值后的局部位置信息 | [−1,1] | 任一方细节分量常量 |
| `rec_jaccard_q75/q80/q90` | 实测与预测高值点集在操作阈值下的 Jaccard | 0–1 | 任一方常量 |
| `rec_jaccard_q75_shuffled_median/p10/p90` | 保持预测取值边际、打乱其空间位置后的 Jaccard 分布（描述性扰动参考） | 0–1 | 同上 |
| `rec_jaccard_q75_excess` | 观测 Jaccard − 打乱中位；**描述性差值，不是 p 值、不用于筛选关系** | [−1,1] | 同上 |
| `rec_centroid_displacement_pitch` | 实测与预测位置权重质心距离 / pitch：位置是否复原 | ≥0，pitch | 任一方常量 |
| `rec_autocorr_near_delta` / `rec_autocorr_outer_delta` | 自相关差：预测是否更平滑/更碎 | ≈[−2,2] | 任一方常量 |
| `rec_decay_half_pitch_delta` | 尺度衰减差（pitch），配两侧绝对值报告 | 实数，pitch | 任一侧不可评价 |
| `rec_patch75_n_delta` / `64_largest_share_delta` / `_largest_radius_delta`、q90 同 | 连通片结构差 | 实数 | 任一方常量 |

**关于自相关的使用限制：** 纯平滑会人为抬高 `geo_autocorr_*` 与连通片聚集度。因此任何"更平滑 ⇒ 更好"的读法都是无效的；`geo_autocorr_*` 在结果中必须与 `rec_local_pattern_pearson`、`rec_jaccard_*`、`rec_centroid_displacement_pitch` 一起读，**不得单独用更高的自相关判优**。

### 4.5 单位内一致性（`within_unit_*`）

同一冻结分组单位含多个切片时（72 片 / 65 单位），比较同单位切片的**实测关系轮廓**一致性：两切片关系向量的 Pearson 相关（`within_unit_relation_pearson`）与关键描述量的差（`within_unit_autocorr_near_absdiff`、`within_unit_centroid_shift_absdiff_pitch`）。这是"目标不稳定"的描述手段之一，不构成独立重复，也不升级为患者级确认。

---

## 5. 分层与聚合

1. **切片内**：所有指标在 `n_s` 个组织内点上计算，不在单位内混池。条码顺序逐点断言；预测包 `barcodes`/`y` 与冻结 cache 必须完全相等，不一致即终止（不静默重排）。
2. **切片 → 单位**：单位内切片等权平均（沿用冻结协议；切片数逐单位记录）。**单位身份只由既定 `unit_id` 与分析键（`split, arm, target, view, model, model_class, cross_organ_extrapolation`）确定；切片级元数据（`source_key`、`organ`、`platform`、`patient_known`、`confirmed_unit`）作为元数据保留，同一单位内不一致时标 `mixed`（并给 `*_n_distinct`），不得因标题或 source_key 把一个既定单位重复赋权。** **某指标若在该单位任一切片不可评价，则该单位的该指标为 `NOT_ESTIMABLE_UNIT`（不取部分切片的均值、不让剩余切片获得整单位权重），并记入 `partial_unit_metrics`；`complete_for_all_metrics` 标出单位是否全指标完整。**
3. **单位 → 宏平均**：`(arm, target, view, model)` 内单位等权；报告 `n_units`、`n_estimable`、`n_units_partial`、均值、中位、p10/p25/p75/p90。**分层（strata）在某一层全不可评价时仍保留行，`mean` 为空并给 `reason='no_estimable_units_in_stratum'`，不静默消失。**
4. **敏感性口径**（同表并列）：
   - `section_weighted`：切片等权（切片多的单位权重更大）；
   - `source_balanced`（仅来源臂）：先在各来源内做单位等权，再对来源等权。**来源归组使用划分级 `held_source`（`source_split_key`），不用切片级 `source_key`（后者可为 `mixed`）。**
5. **缺失值**：任何指标不可评价时该单位该指标记 `NOT_ESTIMABLE_UNIT`，并保留原因计数；**绝不静默删除行或改分母**。宏平均只对可评价单位求值，同时报告可评价单位数。
6. **臂不合并**：内部臂与来源臂分别出记录；同一片同时出现在两臂时是两行，不合并为一条恢复结论。`summarize` 断言 `(arm, section, target, view, model)` 唯一且两臂均存在。**臂级宏平均的来源外推标记为 `true`/`false`/`mixed`（来源臂同时含器官外推与非外推来源时为 `mixed`）；划分级行保留单一布尔值（同一划分内一致）。**
7. **身份与器官外推**：`patient_known`/`confirmed_unit` 取自修正身份表（`patient_confirmed_within_study`），仅用于分层计数：确认 14 片 / 8 单位，其余为未知单位，不升级为独立患者；来源臂再带 `cross_organ_extrapolation`（留出器官无同器官训练时为 True，乳腺与肝脏属于此类）。来源身份不足的 19 片保留在实测清单与内部臂，不进入来源分母。
8. **分母与顺序**：所有交集/覆盖统计基于冻结条码顺序，读取时逐点断言。
9. **命名审计**：`summarize` 检查输出列名中不出现 `pass`/`success`/`significant`/`p_value`/`truth` 等字样，避免把描述性差值写成判定。
10. **跨阶段一致**：`run` 阶段重算的实测描述量必须与 `prepare` 冻结的实测表逐字段相等（`cross_stage_consistency.json`），否则中止；检查阶段另外验证实测描述量对打乱开关不变（`*_measured_descriptor_shuffle_invariance`）。

---

## 6. 数值规则的性质区分（必须写明）

| 类别 | 规则 | 性质 | 是否影响结果筛选 |
|---|---|---|---|
| 计算资格 | 场常量（sd=0）、`0 < sd ≤ sd_floor`（浮点噪声级）、凸包退化、该距离带无点对、模块数 <3、`n_s` 过少致连通片无定义、严格阈值集合为空（`empty_threshold_set`） | **数学可计算性/数值分辨率**，不涉及生物学 | 只决定该格子是否出现数值，不决定"成功/失败" |
| 报告状态列 | `state`：`DESCRIBED` / `DESCRIBED_WITH_NOT_ESTIMABLE_METRICS` / `NOT_TESTABLE_NO_TARGET_GENES` / `NOT_TESTABLE_CONSTANT_READOUT`（实测）与 `NOT_ESTIMABLE_CONSTANT_PREDICTION`（预测） | 仍只是计算资格的可读标签，不是科学判定 | 不改变任何指标 |
| 诊断标记 | `diag_sparse_detection`（任意目标基因检出点占比 < 0.02）、`diag_few_effective_genes`（有效基因数 < 2）、`diag_dominant_gene`（首位基因占比 > 0.5） | **报告便利标记，非门槛**；取值在见数据前写死在合同；同时报告相邻切点（0.01/0.02/0.05、1.5/2/3、0.4/0.5/0.75）的计数敏感性 | **不参与任何筛选、不改变任何指标、不产生状态** |
| 无 | 任何"成功率""覆盖率""相关性"百分比 PASS 线 | 本阶段**不引入**（用户 2026-10-19 边界） | — |

**打乱扰动的性质：** §4.2/§4.3/§4.4 的 `*_shuffled_*` 是**保持取值边际、破坏空间排列**的描述性扰动。它**不保持空间支撑与密度**，未经校准，因此：不报告 p 值、不做显著性判断、不用于筛选关系、不命名为零模型。只作为"观察值与随机重排相比处在什么位置"的参考量级。

### 6.1 数学可算 ≠ 有足够生物学参照（2026-10-19 用户补充边界）

上表的"计算资格/报告状态"**只描述统计量能否计算**。它既不是生物学可用性，也不是评价资格的授予。必须同时明确：

1. **数学上可计算不等于有足够生物学参照。** 一个格子出现数值，只说明该统计量在该数据上存在，不说明该目标场在生物学上已被充分测量、可被独立比较或已有独立参照。
2. **连续证据仍会限制解释，且逐项报告、不合成单一分数**：任意目标基因检出点占比（`sup_spots_with_any_gene_frac`）、单基因/少数基因主导程度（`sup_top_gene_share`、`sup_top3_gene_share`、`sup_effective_gene_count`）、检出基因数与面板外缺失（`sup_genes_detected`、`sup_genes_absent_from_panel`）、目标基因读数占总库比例（`sup_target_count_share_median`）、零点比例（`dist_frac_zero`）、方差档（`sd_band`）、边界/孤立点与邻域度数（`geo_near_degree_*`、`geo_near_isolated_spots`）、经验点间距（`pitch_source_pixels`；HD 的物理尺度未确认，不得据此声称分辨率）、可见模块数量（`prof_n_modules`）。这些量在同一张表内并列，供解释时逐一考量。
3. **不使用统一的百分比规则删除材料**：不因检出覆盖率、有效基因数或任何单一比例删除切片、目标或行；`diag_*` 标记与 `diag_*_sensitivity` 只作报告便利与切点敏感性展示。
4. **不把数学上可算的格子升级为"已具生物学评价资格"**；同样，不把 `NOT_TESTABLE`/`NOT_ESTIMABLE` 读作"该场不存在"或"恢复失败"。
5. **本阶段产出仅是描述性清单**：生物学可用性（哪些场、哪些关系值得进一步作为生物对象处理）仍需主代理结合上述连续证据与项目背景判断，不由本清单代替。下游任何自定判据都必须单独说明理由、适用范围与敏感性，并保留本清单的原始逐项证据。

---

## 7. 输出

`infra/field_biological_review_20261019/` 下：

| 文件 | 内容 |
|---|---|
| `contract.json` | 冻结口径（本协议的机器可读版） |
| `input_checks.json` | 输入哈希、pitch 重算、读出重算比对、条码对齐、基因集隔离、模块定义一致性、资源余量 |
| `prepare_receipt.json` / `splits/*/split_receipt.json` | `InputGuard` 访问记录（含 `human_tls_files_opened`，必须为空） |
| `measured_fields.tsv` | 逐（切片 × 目标）：分子支持 + 分布 + 几何 + 诊断标记 + `state` + 不可评价清单 |
| `measured_summary.tsv` | 实测场资格的分层中位/极值（按目标 × 器官/平台/来源） |
| `measured_module_relations.tsv` | 逐（切片 × 目标 × 视图 × 模块）：实测场与可见模块场的关系 |
| `prediction_metrics.tsv` | 逐（划分 × 切片 × 目标 × 视图 × 模型）：恢复指标 + 几何差 + 不可评价清单 |
| `prediction_module_relations.tsv` | 逐（划分 × 切片 × 目标 × 视图 × 模型 × 模块）：预测场与可见模块场的关系 |
| `unit_metrics.tsv` | 单位级聚合（含可评价计数与原因；元数据冲突标 `mixed` 并给 `*_n_distinct`） |
| `macro_summary.tsv` | 单位宏平均与分位（含 section_weighted / source_balanced 敏感性口径；后者按 `source_split_key` 归组） |
| `macro_strata.tsv` | 按来源划分（`source_split`）、切片级来源元数据（`source_metadata`）、器官、平台、身份、外推标记的分层宏平均；全不可评层保留行并给 reason |
| `within_unit_consistency.tsv` | 同单位多切片的实测与预测关系轮廓一致性。读数使用各表真实列名：实测表 `rel_pearson`，预测表 `measured_rel_pearson` / `predicted_rel_pearson`（`profile_column` 标明）；列缺失时失败关闭（`KeyError`），不静默产生空值。仅描述，不升级为独立重复 |
| `cross_stage_consistency.json` | `prepare` 实测表与 `run` 阶段重算的实测描述量逐字段一致性（7 个字段，最大绝对差必须 ≤1e-12，否则中止） |
| `diagnostic_sensitivity.tsv` | 诊断标记相邻切点的计数敏感性（证明标记不门槛化） |
| `checks/` | 阶段一合成检查与小样本真实检查的原始记录（含旧版本失败日志的 `checks/archive/`） |
| `REPORT.md` | 人类可读描述性小结（由 `summarize` 生成） |
| `freeze_manifest.json` | 全部产物哈希（`summarize` 后生成） |

每个 TSV 行携带身份列（`section_id, unit_id, organ, source_key, platform, patient_known, n_spots`）与 `not_estimable` 原因串，便于追溯。

---

## 8. 本阶段不主张的内容

- 不主张空间增量：本阶段不重做同点/邻域基线对照，也不构造新的同信息预算比较；三组指标只描述既有冻结预测与实测场的关系。
- 不主张任何"恢复成功/失败"的二分判定；不给出通过线。
- 不主张独立 marker 验证：既有模型用完整目标基因集训练，事后拆基因不构成独立留出标志基因验证；本阶段只报告**实测分子参照一致性**，真正的独立 marker 评价需另行隔离并重训（未运行）。
- 不主张真实患者独立确认：已知患者仅 8 个研究内单位，其余未知；切片数不冒称患者数。
- 不主张物理尺度或超分辨率；不主张视野外恢复；不主张因果、通讯或机制。
- 不主张打乱参考构成统计显著性。
- 旧模型之间不新增调参、不训练、不改动既有冻结产物。

**与 PLAN 的关系（探索性）：** 本清单为 HF-01（条件空间恢复）与 HF-04（迁移与可信范围）提供**描述性对照材料**：哪些实测场在分子支持上可描述、预测在场图样上与实测有多接近、以及这种接近是否随目标与来源变化。它不构成 HF-01 的空间增量证据（没有同信息基线对照），也不构成 HF-04 的校准证据。HF-02/HF-03/HF-05 不在本轮范围。**且按 §6.1，本清单不判定生物学可用性，也不把"可计算"当作"可评价"。**

---

## 9. 复核容差（固定于看结果之前）

| 比对 | 容差 | 来源 |
|---|---|---|
| 重算读出 vs 冻结 `y` | 1e-6 | 冻结 `y` 为 float32 存储 |
| 重算 RMSE / 去均值 RMSE / Pearson vs 冻结 `section_metrics.tsv` | 1e-5 | 冻结表由 float32 算术写出 |
| 可见基因矩阵 vs 从 `X` 独立重算 | 1e-5 | 同上 |
| pitch vs 冻结 manifest | 1e-9 | 纯 float64 统计 |

容差不为了通过而放宽；检查中若出现超出容差的差异，作为实现问题处理并保留失败日志。已知案例：冻结表对常数 `train_mean` 基线给出 ±1e-6 的数值噪声相关，本复核按 `sd_floor` 判 `NOT_ESTIMABLE`，并在检查中记录该差异。

## 10. 未运行项（NOT_RUN）

- 全量 `prepare`/`run`/`summarize`（未获授权；阶段一仅跑过抽样 `check` 与限定切片的端到端烟测 `infra/field_biological_review_20261019/smoke/`）。
- 全量训练、重训、独立 marker 隔离与重训。
- 人工 TLS 标签读取、结构验证与人工评价。
- 复杂零模型、密度/支撑保持校准、正式显著性。
- 外部数据获取、GPU 计算。
- 复杂不确定性校准。
