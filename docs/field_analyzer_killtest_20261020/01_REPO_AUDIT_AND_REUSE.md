# 仓库审阅与有条件复用

## 1. 核查范围与当前事实

核查提交为 `8b174e2`，2026-10-20 04:59:42 UTC发布。读取了AGENTS、STATUS、plan、roadmap及本包来源表列出的源码片段和报告；这是**静态审阅，不是独立重训、全库代码审计或服务器文件验收**。[R01–R04]

现行无人工标签轮次为72片、65个既定单位、213,304点；研究内身份确认的只有8个单位，其余57个不能变成“独立患者”。来源留出53片/47单位；部分来源缺失，乳腺/肝脏包含器官外推。四个代理读数为B_AXIS、IFNG response、EMT、OXPHOS。固定邻域对照没有形成稳定的额外空间证据；K=2候选也没有稳定超越同信息普通参照。旧结果仍是当前有效的描述性证据。[R03–R05]

**工作范围纠正：** 旧完整通道隐藏下的同点参照是合法强基线；区域内全表达隐藏下，它只能作为获得额外本地信息的参照。不能靠改名评价任务，把旧失败重判为成功。

## 2. 复用优先级

| 资产 / 精确入口 | 决定 | 可复用的部分 | 必须改写或重新检查的部分 |
|---|---|---|---|
| `scripts/field_discovery_20261018.py`：`strings/read_matrix/assert_grouped/normalized_subset` | ADAPT | HDF5稀疏读取、字符串转换、单位跨折检查、逐点归一化基本运算 | `read_matrix`含形状启发式转置；新适配器必须明确轴语义和barcode映射，不静默猜。重复基因、零深度由新合同处理 |
| 同文件：`InputGuard` | ADAPT | 允许路径和HDF5字段记录思路 | 仅Python访问审计，不是OS沙箱；当前 `record()` 会逐字节hash大文件，不要照搬全池hash；缓存命中也不能跳过掩膜合同 |
| 同文件：`neighbor_graphs/neighbor_mean/feature_views` | REBUILD_FEATURE_PATH | cKDTree/CSR图乘的低层实现可参考 | 全片均值/标准差、全片邻域、孤立点self fallback不适合隐藏区域；禁止隐藏值经图乘或fallback回流 |
| `scripts/field_candidate_subfield_20261019.py`：`split_genes/unit_sections` | ADAPT | 基因名确定性划分、完整单位汇总；目标外分母避免核验基因直接进A | `cache_section`同时保存A/B及全片非零统计、`ok_rows`；不直接给模型完整缓存。跨多个程序须排除**所有核验B的并集**；旧目标专属分母不等于通用场坐标系 |
| `scripts/field_biological_review_20261019.py`：`weighted_centroid/patch_stats` | ADAPT | 点集重合、组件、质心等描述函数；正确单位聚合规则 | 质心不是峰/源；凸包不是组织支持；`largest_radius_pitch`实现取各组件半径最大值，并不明确等于“点数最大组件半径”，新名称与公式需对应 |
| 同文件：`describe_field/band_autocorr/decay_half_pitch` | REBUILD_SCALE | 距离对枚举和有限值处理作为辅助 | `dv=v-neighbor_mean(v)`的相关尺度不作原场尺度；缺口/不相连组织不可被ConvexHull补为实测组织 |
| `scripts/tls_multiscale_gobp.py`：`components` | ADAPT | 在给定合法图上取连通组件 | `tissue_graph`硬编码130.01 μm；仅适用原设计，不可搬到HD。`build_endpoints/prepare`读取TLS/TUM标签，当前分支禁止调用 |
| `scripts/tls_multiscale_forward.py`：`path_coordinates/trace_path` | ADAPT | 已给定折线的弧长/侧距、前驱路径还原 | 弯折、自交、近接平行线、零长度段、非连通图需反例；“最近顶点相邻两段”不是任意曲线的全局精确投影 |
| 同文件：`corridor_weights/endpoint_adjust/operators` | HISTORICAL_ONLY | 旧对照失败和宽度/长度诊断可参考 | 五段四段同号不是连续连接；对称侧带可能抹掉偏侧晕；含标签端点/全路径和原mask语义。不得原样当新桥分析器 |
| `infra/r16/t2p2_20260928/`、`scripts/r16_t2p2_variogram.py` | NEGATIVE_CONTROL_SOURCE | 旧指数+噪声变异函数拟合与模型失配案例 | 正式关联p未校准；不能复用为新显著性。长length但sill极小的结果不解释为大生物场 |
| `docs/TLS_PATIENT_VALIDATION_20261003.md`、患者crosswalk | REUSE_METADATA_IF_PRESENT | 旧GSE18片/17位患者及c_2合组信息，来源区分 | 不与72片HEST单位机械拼接；物理别名需可追踪。监督AUC不能进入新模型选样/调参 |
| 旧72片原始表达、逐点未平滑计数/评分 | CONDITIONAL_LOCAL | 原始观测、几何、冻结身份/来源表 | 当前运行时是否存在未验证；全片目标QC、跨点处理、目标选择过的GO子集不能直接继承 |
| 旧26片标签端点/桥库、74片HookNet图 | SEALED_FOR_DISCOVERY | 历史说明与最终验证计划 | 不是本轮无标签训练、选ROI、定向、选mask、桥真值来源 |

`reuse_register.tsv`提供机器可读版本。ADAPT并非旧脚本重新跑一遍；优先把纯函数移入新适配层，保留来源路径和Git版本。不要为了省几十行代码导入一个会自动读取标签、安装全局hook或改写旧目录的模块。[R06–R10]

## 3. GitHub和服务器不是同一层资产

本轮GitHub API访问 `infra/field_source_holdout_20261018` 返回404，但源码明确引用它。故该目录记为 `REFERENCED_LOCAL_NOT_VERIFIED`，不是已丢失，也不是已在包内。原始数据默认在 `data/other_sources/hest1k/st/{section_id}.h5ad`，只是源码路径证据。

FA-00在服务器只做一次定向清点：相对路径、存在性、字节数、已有manifest版本、数据来源、掩膜前处理范围。优先使用已有freeze指纹，不全盘逐字节hash。确有身份冲突时再对少数相关资产做必要完整核验；头尾摘要只能用于快速比对，不能声称证明两份物理样本独立。

缓存缺失时按顺序：已有原始计数→最小新缓存；仅有可证明逐点且无监督选择的旧评分→限评分层实验；均无则报告对应真实数据分支不可执行，合成路线照常推进。不得自动外部下载，也不重新构建全部历史产物。

## 4. 已明确的新任务风险

- 旧模型已使用同点X，不应再以“邻域模型没有同点信息”解释旧阴性。
- 完整A/B缓存包含核验信息，而且该候选缓存主要保存ok_rows和数值，不直接携带完整barcode/xy契约。新分离器必须重建并核验原始行→barcode→坐标对应，再按观测范围和证据角色拆分，而非只在最终loss上打mask。
- HEST原始 `patient_known` 非空字符串不是确认患者；以修正表 `identity_review.tsv` 的证据字段为准。
- 原 `neighbor_graphs`以经验pitch做邻接，物理距离不能跨平台直接通用；沿组织边界的连线也须检查是否跨真实空洞。
- 旧 `patch_stats`在空高值集返回None，不等于生物阴性；q75会在无结构噪声里强行产生约四分之一“热点”，不能作为桥存在分类器。
- 旧桥旋转参考失败限制显著性，不阻塞重新定义的描述性连接、合成真值恢复和预测实验。[R05,R09–R12]

## 5. 旧结果如何进入新报告

只用于提出机制假设、选择强基线、定义已知失败边界。旧数值原样保留并注明不同归一化/目标/分组/监督角色。新表不与旧RMSE直接相减，除非观测合同、读数、病例、权重和模型信息全部一致。凡已看过的切片，无论新造多少mask，都不是全新外部盲验。
