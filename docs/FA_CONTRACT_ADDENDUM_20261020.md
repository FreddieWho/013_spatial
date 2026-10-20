# FA 观测合同补充条款（2026-10-20，基于 Grok preflight 的 GAP 关闭）

**性质：** 本文件是 `docs/field_analyzer_killtest_20261020/05_COMMON_PROTOCOL.md` 与 `07_ALGORITHM_ROUTES.md` 的项目侧补充条款。施工包目录是只读来源快照，不修改其正文；凡本文件与包内文字冲突，以实现本文件为准，并在 `infra/field_analyzer_20261020/` 的报告中注明冲突点。

**依据：** `infra/field_analyzer_20261020/review/preflight_fa00.json`（审核代理 `a5bcb012`，Provider=pi，Model=`xai/grok-4.7`，Thinking=high，Revised=文档与少量源码入口的只读预审）。结论：`allow_fa00_implementation = true`，无 BLOCKER，7 项 GAP 可按下列固定断言在实现中关闭。审核范围未覆盖实现、未运行毒化反例、未核验本地原始文件是否存在。

以下编号沿用预审的 items id。

## 1. 掩膜覆盖的入口清单（item 1，GAP）

05 §3 已点名可测性、方向、全片 PCA、预平滑、查询资格、隐藏 0、孤立点自我回退；**必须再显式覆盖**：

1. 文库大小／归一化分母（含"排除全部核验 B 的并集"的外部参考分母）；
2. 逐点深度与检出率（及其任何分位或阈值）；
3. 整片均值、方差、分位等全片统计（旧 `feature_views` 会把整片均值/标准差广播进每一行，先算后掩膜即泄漏）；
4. **按表达加权的邻接或图权重**（必须只由可见值构造，权重归一化后记录支持量与缺失标记）；
5. 查询点资格／`ok_rows`（不得由 H 内表达决定）。

判据：先掩膜再聚合。任何"先在完整片上算一个数、再把这个数发给模型"的路径都算泄漏，即使那个数看起来无害。

## 2. 核验组的强制方式（item 2，GAP）

布尔开关不足以强制"排除全部核验 B 的并集"。接口必须要求**显式的基因角色表**（`gene_role`：`A_FEATURE / B_WITNESS / EXCLUDED / BACKGROUND`），并：

- 任一 `B_WITNESS` 基因出现在任何输入层的基因轴上即报错停止；
- 禁止整包读取旧 `cache_section` 产物（A 图、B 图与全片统计同处一份缓存）；如需复用，只能按显式列白名单重读并逐列核对角色；
- 归一化参考基因集与 A、B 目标集合的交集必须为空，并在报告中打印该交集的实际大小（必须为 0）。

## 3. 轴顺序必须可失败（item 3，GAP）

`ObservationBundle` 必须带显式轴声明（例如 `axis_order: "obs_by_gene"`、`n_obs`、`n_genes`），读入时校验：

- `shape` 与 barcode 数、gene 数双向一致；
- 若 `shape` 两种解释都自洽，**报错**并要求显式声明，禁止任何形状启发式转置；
- 转置、barcode 重排、重复 symbol、缺基因各自有定向测试。

## 4. 单位角色表（item 4，GAP）

FA-00 必须交付 `infra/field_analyzer_20261020/unit_role_assertions.tsv`，列至少：`unit_id, section_id, source_id, patient_status, n_points, role(training/development/evaluation), mask_group, rotation_group, is_replicate_of`；并附断言：

- 同一 `unit_id` 的多片同角色；同一 `base_id` 的 mask/旋转/生成器变体同角色；
- 重复 mask 不计入"额外患者"；统计单位是单位，不是 spot；
- 未知身份不得补成独立患者，只能按既定单位分析并标注。

## 5. 指纹对象的显式清单（item 6，GAP）

毒化反例必须声明"哪些对象参与逐位比较"，至少：

`values_visible`（数组+形状+dtype）、`visible_mask`、`tissue_support`、`barcodes`、`xy`、`gene_role` 解析结果、`input_provenance`、特征变换参数、邻接/权重矩阵、查询点列表与资格掩膜、模型参数、预测数组。

同时必须证明两件事：
- (i) 只改 H 内表达／只置换 B／只移动历史标签文件时，上述对象全部逐位不变；
- (ii) 评价结果**确实改变**（否则测试是空的，例如把 H 值改掉却没有任何评价量随动，说明评价没读到 H）。

指纹用数组内容哈希或 `max|Δ|=0` 报告，不得只写"未报错"。

## 6. 资产状态必须是三态（item 7，GAP）

`reuse_resolution.tsv` 的 `availability` 在包内四态之上，对"存在但不可足量使用"单列状态，禁止混入 `NOT_PRESENT`：

`VERIFIED_LOCAL / PRESENT_BUT_INSUFFICIENT / NOT_PRESENT / REFERENCED_LOCAL_NOT_VERIFIED / TRACKED_SOURCE_ONLY`

`PRESENT_BUT_INSUFFICIENT` 必须写出实测值（实际单位数、实际片数、缺失片 ID、barcode 对不齐的行数占比）。不得用旧评分把不足的队列补满，也不得把缺片写成阴性。

## 7. B 路线精度式的错误修正（item 8，GAP；FA-03B 前置，非 FA-00 前置）

07 §2.3 的能量

`E(f)=1/2||Mf−y_O||²_{R⁻¹} + τ/2||f−μ0||² + 1/2 Σ_l ||W_l^{1/2}(L_l f − b_l)||²`

对应梯度 Hessian 为 `MᵀR⁻¹M + τI + ΣLᵀWL`。包内写出的 `Q = τI + ΣLᵀWL` **漏掉观测项 `MᵀR⁻¹M`**，与已写对的右端项 `MᵀR⁻¹y_O + τμ0 + ΣLᵀWb` 不自洽：用该 Q 求解，观测噪声不约束解。

修正口径（实现时二选一，均需记录）：

- **A（带观测项）：** `Q = MᵀR⁻¹M + τI + ΣLᵀWL`，右端项不变；`Q` 仍正定（`MᵀR⁻¹M ⪰ 0`）。
- **B（块条件化）：** 真正对可见块做条件化（Schur 补），不能直接把观测项删掉后当等价。

实现要求不变：float64、稀疏/矩阵自由、共轭梯度 + 小规模直接解对照、记录残差与迭代；不做全片稠密求逆。

## 8. FA-00 范围微调（item 9 结论 + scope 建议）

保留原三份交付，另加：

- `unit_role_assertions.tsv`（第 4 条）；
- 小几何样本中**包含至少一个被掩膜的场值**，用于检验"可见区算统计、隐藏区不算"的路径；
- 毒化反例覆盖第 1 条的五个入口，不只覆盖表达矩阵。

不加：能量式求解、旧 p 值修复、全库普查、72 片清点。

## 9. 不适用的延伸

预审第 5、10 项为 PASS，无需补充。预审未评估：实现正确性、真实文件存在性、`NOT_ASSESSABLE_FROM_DOCS` 项（如"图乘在实现中是否真的按可见权重归一化"）——这些留给实现后的独立复核。
