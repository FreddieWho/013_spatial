# FA-00 交付报告：最小观测适配与资产清点

**日期：** 2026-10-20
**节点：** FA-00（基础设施；不产生 AF 证据）
**实现：** 主代理直接实现（D-283），非分派
**口径：** `docs/FA_CONTRACT_ADDENDUM_20261020.md` 优先于施工包正文

## 1. 交付物

| 文件 | 内容 |
|---|---|
| `field_analyzer/observation.py` | `ObservationBundle`、`GeneRoles`、掩膜统计、指纹、可见/hidden 物理分离 |
| `field_analyzer/adapters.py` | 显式轴声明的 h5ad 读取、重复 symbol/零深度/派生 obs 字段处理、header 快速探测 |
| `field_analyzer/geometry.py` | 经验 pitch、kNN 邻接（含最大边距）、距离对、局部方向、空洞/掩蔽区分 |
| `field_analyzer/roles.py` | 按基因名的确定性 A/B 拆分、witness 并集、角色注册表 |
| `field_analyzer/pipeline.py` | **唯一信息入口** `build_inputs`（五类入口全部在掩膜后重算）、`build_truth`（仅评价）、`ridge_predict`（确定性参照） |
| `field_analyzer/poisoning.py` | 毒化反例 harness（P1/P2/P3 + 反向断言 + 静态源码扫描） |
| `field_analyzer/cli.py` | `assets` / `units` / `poison` / `geometry-smoke` 四个子命令 |
| `infra/field_analyzer_20261020/reuse_resolution.tsv` | 12 项资产，五态可用性 |
| `infra/field_analyzer_20261020/section_asset_check.tsv` | 14 张片的头部核验 |
| `infra/field_analyzer_20261020/unit_role_assertions.tsv` | 8 单位 / 14 片的分折与角色断言 |
| `infra/field_analyzer_20261020/poisoning_report_{ZEN36,NCBI642}.json` | 毒化结果（含逐数组 max\|Δ\|） |
| `infra/field_analyzer_20261020/geometry_smoke/*` | 带洞/窄颈合成几何的数值结果 |
| `tests/test_fa_{observation,adapters,geometry,roles,poisoning}.py` | 32 个定向测试 |

## 2. 资产清点（实测，不是引用）

- 身份表实际 72 行，其中 `patient_confirmed_within_study == True` 者为 **8 个研究内单位 / 14 张片**，与施工包文字一致；其余 **58 行为 unknown**，不升级为独立患者。
- 14 张片全部存在，`X` 均为 `csr_matrix`，行数 = barcode 数 = `obsm['spatial']` 行数，**14/14 轴一致**；基因轴统一 36,601。片大小 328–2,385 spots，单文件 86–110 MB。
- 复用裁定：`reuse_resolution.tsv` 12 项全部 `VERIFIED_LOCAL`。其中 `scripts/tls_multiscale_gobp.py` 标 `SEALED`（标签端点在发现分支禁用），三个 legacy 脚本标 `ADAPT`/`REBUILD`，两个历史结果为 `HISTORICAL_ONLY`。
- 本轮**未**全量逐字节 hash；只做存在性、字节数、头部形状与 barcode 对齐核验。

## 3. 观测合同实现要点（对应补充条款）

1. **显式轴**：`AxisSpec` 是必填语义；形状与 obs/var 不符时若恰好等于转置只报错并提示"声明 `gene_by_obs`"，**绝不静默转置**（旧 `read_matrix` 会转）。定向测试覆盖转置文件、声明不符、方阵歧义。
2. **角色注册表**：`GeneRoles` 强制 witness 基因不能出现在任何输入基因轴；A/B 按基因名 + 固定种子拆分（`split_field_genes`），与表达读取无关，且对输入顺序不变。
3. **外部文库分母**：`library_scope = all_genes_excluding_all_analyzed_field_genes`，即四个场的 A∪B 全部排除，仅用其余基因求每 spot 文库；实测 `n_library_genes = 36601 − 333 = 36268`（四个场 49+200+200+200 定义基因中实际存在 48+200+196+200）。这比"只排除 B"更严格，也避免 A 基因进入自己的分母。
4. **掩膜先于统计**：文库、深度、检出率、整片均值/方差、场评分、邻接权重全部只由可见行计算；隐藏行在 `values_visible` 中为 NaN，`ObservationBundle.__post_init__` 拒绝隐藏行出现有限值。
5. **禁止 self fallback**：孤立点邻居值为 NaN，绝不复制自身值（旧 `neighbor_mean` 会复制）。邻接有最大边距 `3.0 × pitch`，断开组织不会被拉成邻居。
6. **尺度不硬编码**：pitch 由数据的中位最近邻距离估计，单位随 `coordinate_unit` 记录；旧 130.01 μm 未进入新代码。
7. **空洞 vs 掩蔽**：`region_classes` 分 VISIBLE / MASKED_POSSIBLE_TISSUE / NON_TISSUE；`assert_edges_avoid_non_tissue` 拒绝对穿真实空洞；不使用凸包补组织。

## 4. 毒化反例（本节点核心证据）

同一几何掩膜（以中位组织坐标为中心、半径 6×pitch 的圆盘，纯几何、不用表达）下，把原始矩阵改造成三种毒化版本后**从头重算全部模型输入**：

| 情形 | ZEN36 | NCBI642 |
|---|---|---|
| 隐藏 spot 数 | 122 | 128 |
| P1 隐藏区表达置为 +1000×max | 22 个数组全等，max\|Δ\|=0 | 同 |
| P2 隐藏行互相置换 | 全等，0 | 同 |
| P3 witness(B) 基因跨 spot 置换 | 全等，0 | 同 |
| 参照模型预测 | max\|Δ\|=0 | 0 |
| 评价量是否随动 | 是（B_AXIS 隐藏区真值 max\|Δ\|=0.513；EMT 1.01） | 是 |
| 静态扫描：本包无标签/TLS 读取器 | clean | clean |

被比较的 22 个数组覆盖五类入口：`library`、`depth`、`detection`、`slide_mean/slide_std`、`neighbor::*`（表达加权邻接）、`neighbor_struct::*`、`query_eligible`、四个 `field::*`、`degree`/`isolated`/`pitch`、`score_genes`。反向断言（评价必须变）单独检查，避免"评价根本没读隐藏区"的空测试。

反向对照（证明断言有检出力）：`tests/test_fa_poisoning.py::test_leaky_slide_statistic_would_be_caught` 用"掩膜前算整片均值"的写法，毒化后数组出现显著差异而被判失败；同一数据按可见行计算则完全一致。

## 5. 几何冒烟（合成，已知答案）

64×64 规则格，中心圆洞 + 一条穿过洞的窄颈被标为支持组织；另有一个与场重叠的掩膜圆盘。

- 4096 点，支持 3802，pitch = 1.0，中位度 6，**孤立点 294**（全部为不支持的空洞点），孤立点邻居值 NaN。
- 距离对 12 个距离箱，全部端点都在支持集内；**穿越非组织的边 = 0**。
- 局部方向：全可见时 92.8% 点可测；加掩膜后 88.5%，且**掩膜区内角度全为 NaN**。
- 掩膜统计：掩膜数组的可见均值与直接对可见值求均值完全相等，隐藏位置无有限值。

## 6. 测试与失败

- FA 定向测试：**32 passed**（observation 11 / adapters 8 / geometry 6 / roles 5 / poisoning 4，含负对照）。
- 全库测试：438 passed，4 failed，5 subtests passed（746 s）。4 个失败全部在 legacy 套件：`tests/test_r16_graphst_bridge.py`(2)、`tests/test_tls_multiscale.py`(1)、`tests/test_tls_patient_validation.py`(1)，错误为 numpy 2.x 的 `solve` 签名不兼容（`gufunc signature (m,m),(m,n)`），与 `field_analyzer` 无关；`git status scripts/` 为空，legacy 代码未被我改动。**未修复**，记为环境既有问题。

## 7. 资源

| 项 | 实测 |
|---|---|
| 墙钟 | assets 3 s；units 1 s；geometry-smoke 0.7 s；单张片毒化 27–30 s |
| 峰值内存 | 毒化 1.29–1.40 GB（主要来自 36601 列稀疏子集与 1983×333 稠密化） |
| 线程 | 本节点无大矩阵分解；BLAS 线程未固定（无 BLAS 重负载）。FA-01 起需固定并记录 |
| 磁盘余量 | 运行后 697.9 GB（≥200 GB 要求满足） |
| GPU | 未使用 |
| 下载 | 无 |

## 8. 未做与限制

- 未实现 FA-01 的尺度曲线（半变异函数）；本节点只到"可测性接口 + 几何正确性"。
- 未做任何真实拟合、掩膜恢复或路线比较；因此**不产生 AF-01～AF-04 的任何证据**。
- 58 个 unknown 身份单位未纳入；pilot 只用 8 个已确认单位。
- 未运行 MISTy / LSGI / NicheCompass / TrimNN 等公共方法（I-054 仍 OPEN）。
- 毒化只在前两张片（ZEN36、NCBI642）与合成数据上跑通，未覆盖全部 14 片；FA-04 前应至少补到每来源一张。
- 本报告不声称合同无缺陷；独立复核仍待执行（D-283：节点完成复核，Grok `xhigh`）。

## 9. 对 AF 假设的影响

**无直接影响。** 本节点是基础设施：它把"模型能看到什么、看不到什么"变成可执行且可被毒化检验的接口。按 roadmap 规则，FA-00 标记为基础设施，不计入 AF 支持。

## 10. 复现

```bash
cd /home/huyudi/013_spatial
/opt/anaconda3/bin/python -m field_analyzer.cli assets
/opt/anaconda3/bin/python -m field_analyzer.cli units
/opt/anaconda3/bin/python -m field_analyzer.cli poison --section ZEN36
/opt/anaconda3/bin/python -m field_analyzer.cli poison --section NCBI642
/opt/anaconda3/bin/python -m field_analyzer.cli geometry-smoke
/opt/anaconda3/bin/python -m pytest tests/test_fa_*.py -q
```
