# FA-13 等价性核对：是否已有“目标自身邻域 + 另一场邻域”的公平嵌套比较

核对日期：2026-10-23。核对对象：FA-12、FA-05、旧 HF 的 M0/M1/M2、旧 FC 的 X/WX/VX。
判定口径（合同 `docs/FA13_CONDITIONAL_RELATION_20261023.md`）：只有同时满足下列五项才算**完全等价**，否则一律 `NOT_EQUIVALENT`，不得只因名称相似而引用为已验证结果：

1. 同一个目标关系（主关系 EMT ← IFN-γ；对照 IFN-γ ← OXPHOS）；
2. 相同观测预算（同一批 14 片 / 8 个研究内单位、同一固定查询点 Q、Q 的表达/深度/QC 不进入输入）；
3. 同一独立基因半边（特征半边与核验半边分离，且核验半边不进入任何输入）；
4. 单位级留出（留一单位，单位不跨训练/评价）；
5. **模型包含目标自身的可见邻域项**，再与“自身邻域 + 另一场邻域”做嵌套比较。

## 逐项核对

| 来源 | 目标/关系 | 观测预算与查询点 | 基因拆分 | 单位留出 | 含目标自身邻域？ | 判定 |
|---|---|---|---|---|---|---|
| FA-12 `field_analyzer/spatial_lag.py`（`compare_pairs`，提交 `06e5d59`） | EMT←IFN、IFN←OXPHOS 等四对 | 14 片 / 8 单位；10% 掩膜（同 seed 0）；但 `_prepare` 在**掩膜前**对全片 residualise/zscore，`section_scores` 由 `build_inputs(...support,support...)` 在全片评分（即 Q 参与标准化/深度系数） | 有 A/B 半边（`roles.build_field_splits`），特征只在 A 半边 | 有（留一单位） | **没有**：`same_neighbour`、`cross_neighbour`、`same_spot` 是三个**各自单独**的预测器，各自单独评价；不存在 M0（自身邻域）M1（+跨场）嵌套模型 | **NOT_EQUIVALENT**（单独比较 ≠ 条件增量；且预处理含 Q） |
| FA-05 `field_analyzer/block_recovery.py::predict_multifield` | 目标场恢复（任一目标） | 阻塞中心 2 pitch 核心；可见预算相同 | 其它场特征 | 单元级 | **没有**：只用“其它场”的平滑值，目标自身邻域不在模型中 | **NOT_EQUIVALENT**（缺少 M0 项，无法回答条件增量） |
| 旧 X/WX/VX（`docs/FIELD_RAW_GENE_SPATIAL_CONTRACT_20261019.md`，脚本 `scripts/field_raw_gene_spatial_20261019.py`，提交 `0604dfc`/`0ca330f`） | 单目标逐基因可恢复性（EMT/B_AXIS 等），非跨场关系 | 冻结可见矩阵、512 点/单位抽样 | 无 A/B 半边（同一目标自身的基因与 QC） | 有（单位五折与来源留出） | 有 **同场** 邻域（WX 近邻、VX 外环），但**没有另一场的邻域块** | **NOT_EQUIVALENT**（无跨场项，不是 M0 vs M0+M1） |
| 旧 HF M0/M1/M2（`docs/IMPLEMENTATION_REPAIR_20260921.md`、`docs/IMPLEMENTATION_AUDIT_20260921.md`） | 旧 GO 读数场预测（ST-CRC / USZ 队列） | 不同队列（7 患者/14 片、8 患者/8 片），不是本 pilot 的 8 单位/14 片；旧 M1/M2 被审计记录为**隐藏表达泄漏**（标准化在全片拟合） | 旧口径 | 旧留一 | M0/M1 是旧模型代际，不是“自身邻域 vs 自身+跨场” | **NOT_EQUIVALENT** 且旧数值不可用于本节点 |
| 旧 FC01/FC02/FC03（`scripts/fc01_field_recovery.py`、`fc02_input_ladder.py`、`fc03_compare_fixed.py`） | 输入阶梯/字典探针；源码检索无 WX/VX 联合块 | 与 FA-13 不同预算（4213 面板外基因 vs 23 固定分数等） | 无 A/B 半边 | 旧划分 | 无目标自身邻域的嵌套项 | **NOT_EQUIVALENT** |

## 结论

**本地不存在与 FA-13 完全等价的联合条件比较**，因此本轮必须按合同实现并运行；不引用任何旧结果作为本检验的答案。

需要特别记录的三点：

- FA-12 的 0/8 只说明“跨场邻居**单独**没有超过同场邻居单独”，不能推出“已知自身邻域后跨场增量为零”。这正是 FA-13 补的缺口。
- FA-05 的多场参照缺少目标自身邻域项，不能承担条件增量的判断。
- 旧 HF 的 M0/M1 既不同队列也有已记录的隐藏表达泄漏（审计已定性），本轮不复用其数值，也不把它当作“已做过”的证据。
