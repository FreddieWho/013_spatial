# 技术路线：空间结构反演

本文档的角色(需写入文档开头): 本文档只回答"为了回答 plan.md 中的问题,需要依次做什么"。它是可执行的,会频繁变动。本文档不重述科学问题的动机,不论证其重要性——需要时引用 plan.md 中的条目编号。本文档也不记录"为什么当初这样选",那属于 decisions.md。

## R-00｜节点图与关键路径

### 主关键路径

`R-01 → R-02 → R-04 → R-05 → R-06 → R-07 → R-08 → R-09 → R-10`

### 邻近平面分支

`R-03 → R-11`

### 主动观测分支

`R-06/R-11 → R-12 → R-13`

### 条件探索分支

`R-07 或 R-09 暴露明确失败模式 → R-14`

- **最早可否证节点：R-05。**
- **不可逆投入：R-13。**
- **允许并行探索：R-08、R-09、R-14；但不得越过各自依赖节点。**
- **当前数据范围约束：仅使用已经收录的数据；CROST、Synapse/HTAN 原始数据与受控数据维持 HOLD。确切核心队列名单为 `[待定]`，需要完成 R-01 后依据独立患者/组织块、重复收录、结构 GT 和专用验证价值确定。**

## R-01｜现有数据去重、注册与角色冻结 `[基础设施｜COMPLETE_WITH_EXCLUSIONS]`

1. **要做什么：** 对当前已收录肿瘤空间数据建立唯一物理样本注册表，恢复 `study_id / patient_id / block_id / section_id / z_position`，识别跨 HEST、GEO、10x、HTAN、STOmicsDB 等来源的重复样本，并把数据标记为 discovery、training、internal validation、external validation、high-resolution validation、serial-section validation 或 archive。
2. **服务的假设：** `[基础设施]`，不直接检验科学假设。
3. **完成判据：** 每个进入分析的数据单元可追溯到患者和真实组织块，或满足下述 fail-closed 条件的 patient-linked physical specimen；同一物理样本不承担相互冲突的角色；当前下载集合完成去重；核心工作集不超过已讨论的 6–10 个逻辑数据单元。
4. **对假设信心的影响：** 不直接改变任何假设的可信度；若无法恢复患者或组织块身份，则相关数据不能承担确认性证据。

**投入上限：** 不恢复当前 HOLD 的数据，不进行全量 RAW 重处理；人工投入上限为 `[待定]`，需要项目可用人力和可接受延迟才能确定。

**R-01 注册与冻结规则：**

- 一个逻辑数据单元定义为 `canonical source study lineage × acquisition/platform protocol × 可审计物理身份层级`。同一 study 的不同下载入口、聚合镜像、同一 block 的不同 section 或同一 section 的重处理版本不得拆成独立单元。
- 任何进入 claim-bearing 核心集的单元都必须由数据自带或经批准获取的官方 metadata 支持 patient 与真实 block 身份，或满足窄定义的 `patient-linked physical specimen` 例外：官方 patient 映射、实体组织描述和稳定 specimen locator 至少形成 E2 互证；`block_id` 必须保持空，并单独登记 specimen identity、equivalence basis 与失效条件。GSM、BioSample、sample、slide、capture area、section、目录或文件名都不得写入或静默等同于 block。
- patient、block、serial、重处理或 possible-match 任一关系可能相同时，在切分和角色分配上采用保守 leakage group；不同 checksum 不能证明物理独立。
- 角色冻结只使用身份完整性、本地资产存在性、来源谱系、平台、癌种和预注册 capability 等 outcome-blind 信息。TLS 阳性率、结构数量、场强或任何模型性能不得用于挑选或更换 validation 单元。
- 真实 patient/block 与 specimen-equivalent 都不能仅凭 physical row 自报进入 gate：patient/block 或 patient/specimen/basis evidence 必须引用校验和通过的 metadata asset，normalized value 与 raw key/value 必须闭合；canonical source lineage 必须在 duplicate group 中登记，角色 leakage group 必须与该 lineage 一致。
- 节点允许 `COMPLETE_WITH_EXCLUSIONS`：核心 6–10 个逻辑单元满足全部判据，其余库存明确退出确认性用途。若穷尽本地自带 metadata 后仍不足 6 个 patient/block 可追溯单元，状态为 `BLOCKED_IDENTITY`；若去除镜像和保守 leakage group 后无法形成独立 external validation lineage，状态为 `BLOCKED_INDEPENDENCE`。

**2026-07-31 完成证据：** 用户批准后仅获取官方 GEO sample/series metadata、PubMed/OA provenance metadata 和一个 patient—GSM crosswalk supplement；未保留表达矩阵、图像或结果表。GSE211956、GSE226997、GSE274103、GSE274557 分别通过 E2 patient-linked physical specimen 审计，各只增加 1 个 canonical logical unit，并与既有 Atlas/本地镜像合并为同一 source lineage。加上 HTAN Vanderbilt CRC 与 10x Breast Block A，956 条注册记录中有 121 条合格 physical records，折叠为 6 个逻辑单元。

**冻结结果：** GSE274557 为 discovery，HTAN Vanderbilt CRC 为 training，GSE274103 与 GSE226997 为 internal validation，GSE211956 为 external validation，10x Breast Block A 为 serial-section validation。角色只使用身份、来源、平台、本地资产和预注册 capability；未使用 TLS、结构结果或模型表现。机器 gate 为 `COMPLETE_WITH_EXCLUSIONS`，R-01 对 plan.md 各假设仍无直接影响；R-02 的结构 GT 审计可开始，但不得把 R-01 的 specimen equivalence 当成结构 GT。

**2026-08-03 增补注册（D-037 三路径）：** 用户批准三条公开来源路径后，新增三个 validation 逻辑单元并完成注册与角色冻结：GSE175540（Meylan ccRCC，24 个 GSM physical units，一人一样本，external validation）、ST_CRC_CMS（Valdeolivas CRC，7 患者 × 2 连续切片共 14 个样本，internal validation）、TLS_VISIUM_USZ（Zenodo 14620362，肾癌 3 + 肺癌 5 共 8 个样本，external validation）。注册记录 1002 条中合格 physical rows 由 121 增至 167，冻结逻辑单元由 6 增至 9；机器 gate 保持 `COMPLETE_WITH_EXCLUSIONS`。角色分配只使用身份、谱系、平台与预注册 capability；GSE226997、GSE274103、GSE211956 的既有冻结角色不变（保持 GT-less，不承担确认性主张），GSE274557 维持 discovery。本节点仍为基础设施，对 plan.md 各假设无直接影响。

## R-02｜结构本体、GT 重叠和切分规则冻结 `[基础设施｜PARTIAL_GT_READY]`

1. **要做什么：** 为已讨论的 TLS、血管、坏死和肿瘤—基质边界建立结构实例注册，作为已知结构的验证与定位锚点；记录几何类型、证据等级、GT 定义模态、需要排除的直接输入通道、边界不确定性、同一物理结构在相邻切片中的对应关系和用途场景；不把这些已知结构当作潜在空间场的完整候选集合；按 `block_id` 冻结外层切分。
2. **服务的假设：** `[基础设施]`，不直接检验科学假设。
3. **完成判据：** 每个进入命名结构确认性定位的实例具有唯一 `instance_id`，GT 生成信息和允许输入已登记；任何同一结构的跨切片化身不会跨折；没有明确结构锚点的潜在空间场可以进入发现与重建路径，不因缺少 GT 被排除；若要提出第二个命名结构 claim-bearing 主张，再按结构级患者数、实例数、分辨率和独立 GT 确定。
4. **对假设信心的影响：** 不直接提高假设可信度；若 GT 与输入无法形成可审计边界，则该结构只保留发现性用途。

**投入上限：** 在第二结构确定前，仅处理上述四类已讨论结构，不扩展新的结构本体；人工标注上限为 `[待定]`，需要病理方可用工时确定。

**R-02 GT 与切分规则：**

- 确认性 GT 必须关联 R-01 physical unit，并提供内容寻址、可重放的 mask、polygon、centroid 或等价空间 locator；presence/count、结构 ID、类别、面积或样本级阳性不能替代空间几何。
- GT 生成模态、方法、坐标系、分辨率、边界不确定性和输入依赖必须闭合。表达或同一拟用输入派生的标签只能用于 discovery；来源未知时 fail closed。GT 本身、到 GT 的距离、GT 决定的遮蔽形状及生成依赖的传递闭包全部禁止进入预测输入。
- outer group 是同患者、同真实 block、重复/重处理、serial-section、显式跨切片同结构关系和 R-01 leakage edge 的保守连通分量。真实 block 优先；block 未知的 E2 specimen-equivalent 只按患者级 envelope 冻结，`block_id` 保持空，不声明 block-level holdout。
- 同一 block 不自动证明切片相邻或同一结构跨切片；只有 metadata 明确给出对应关系时才建立 cross-section structure link。未标注区域保持 unknown，不能静默当作 negative。

**2026-07-31 审计结果：** 6 个冻结逻辑单元中共识别 57 条 Table S4 真实 TLS ID 候选，但均缺少可重放空间几何，且现有上游流程包含 marker scoring、phenotype inference 与 TLS segmentation，无法建立独立 GT—输入边界。本地另有 19 个文件名含 `TLS_annotation` 的非核心候选资产，但 provenance 与 physical-unit 映射未闭合。血管、坏死和肿瘤—基质边界的可审计实例均为 0，第二结构不能冻结。四个 specimen-equivalent GEO lineage 的 72 条 physical rows 只形成 30 个患者级保守 envelope；加上 HTAN 与 10x 后，所有 121 条合格 physical rows 共冻结为 61 个 patient-wide outer groups，并单独保留 32 个真实 block groups。

**2026-07-31 解阻证据（范围更新获批后）：** 用户批准扩大外部数据范围（D-032）后，获取并注册了唯一公开合格的 GT 来源：Heiser et al. Cell 2023（PMC10756562）commit 钉死仓库的 42 个逐 spot 病理标注 CSV 与 sample key。审计发现身份粒度必须是 piece 而非 capture area（5 个 capture 为 TMA 多芯、分属不同患者；Visium barcode 跨 capture 复用，仅 capture 内唯一，D-034），47 个合格 piece 全部经 h5ad spot 索引重放坐标（D-033）。结果：42 个可审计 GT source（1 个 piece 因本地缺 h5ad 挂 REPLAY_BLOCKED，1 个 capture 因 R-01 排除挂 UNLINKED）；195 个确认性实例——TLS 4 个（3 患者）、肿瘤—基质边界 191 个（19 患者、30 piece）；另有 28 个癌前/正常上下文实例仅作 discovery。第二 claim-bearing 结构冻结为 TUMOR_STROMA_BOUNDARY（D-035），血管与坏死的可审计实例在所有来源中仍为 0。

**2026-08-03 解阻证据（validation 谱系，D-036/D-037/D-038）：** 三条公开沉积路径全部通过审计并注册为带可执行 verifier 的 GT 来源：GSE175540 作者沉积逐 spot TLS 标注 CSV（18 个可审计 source、35 个 TLS 组件/18 患者；3 个 T_agg-only 文件登记 NOT_AUDITABLE，1 个样本无沉积标注文件，2 个零 TLS 样本为已审计阴性无行）；TLS_VISIUM_USZ 沉积 h5ad 内 `obs/ground_truth` 人工标签（8 个 source、108 个 TLS 组件/8 患者；D-038 授权仅读 `_index`/`x_array`/`y_array`/`ground_truth` 列）；ST_CRC_CMS 病理医生逐 spot 类别 CSV（12 个 source、551 个肿瘤—基质边界组件/6 患者；`tumor&stroma` 家族 4 标签入范围，`IC aggregate*` 类标签 fail-closed 排除出 TLS 范围）。实例定义与 Heiser 一致（六边形网格连通组件），几何全部经同源坐标资产重放并索引钉死。确认性实例总数 889：TLS 147（HTAN 4、GSE175540 35、USZ 108），肿瘤—基质边界 742（HTAN 191、STCRC 551）；可审计 GT source 共 80 个。KIRC 谱系另贡献 30 条 Table S4 真实 TLS ID 候选，scoped inventory 由 57 增至 87。GSM↔Atlas R_P 逐样本 crosswalk 未公开，GSE175540 与 Atlas KIRC 行的重复链接保持 study 级；重算 35 组件与 Atlas 30 条的差异为定义粒度差异，正式计数以可重放重算为准。

**门禁与假设影响：** 机器状态为 `PARTIAL_GT_READY`：training 谱系（HTAN Vanderbilt CRC）与三个公开沉积 validation 谱系（GSE175540、TLS_VISIUM_USZ 为 external，ST_CRC_CMS 为 internal）均已具备可审计实例级 GT，TLS 与 TUMOR_STROMA_BOUNDARY 对这四个谱系冻结为 claim-bearing。该 GT 是已知结构的验证锚点，不是发现未命名潜在空间场的前置条件；血管与坏死也不是当前主路径的必需 GT。R-04 至 R-07 的训练侧工作（分子-only 输入、冻结外层切分内）可以启动；R-04 完成判据中的队列间方向一致性与 R-05 的两来源一致性判据恢复可执行（限上述四个有 GT 的谱系）。R-02 是基础设施节点，本结果不直接支持或否定 plan.md 的任何假设，但使 H-01/H-03/H-05 的确认性检验首次可在 training 与 validation 双侧执行，H-04 因第二结构冻结恢复可检验路径。剩余缺口：GSE226997、GSE274103（internal）与 GSE211956（external）仍无可审计 GT，保持 GT-less 不承担确认性主张；血管与坏死在所有已查公开来源中可审计实例仍为 0（GSE274103 的作者索取邮件降为可选，其独占剩余价值仅为未冻结结构的标注；跟踪于 I-001）。

## R-03｜相邻切片与配准可测性普查 `[基础设施]`

1. **要做什么：** 在当前库存中统计连续组织块、切片顺序、切片间距、切片厚度、可用 WSI 或正交 GT、可达到的配准误差，以及当前平面与目标结构的相交关系。
2. **服务的假设：** `[基础设施]`，为 H-06 和 H-07 判断是否可测提供前提。
3. **完成判据：** 形成组织块级清单；能够区分同平面遮蔽、邻近平面预测、抽中间平面和栈外外推；若关键字段缺失，明确哪些任务不可测。
4. **对假设信心的影响：** 不直接改变 H-06；若没有足够独立组织块或配准误差接近目标误差，H-06 保持未判定而不是被视为失败。

**投入上限：** 只普查当前已下载数据，不为增加普通横断面样本恢复外部下载；是否申请独特 serial-section 数据为 `[待定]`，需要本节点结果后另行决策。

## R-04｜潜在空间场发现、距离场复现与组成—状态拆分

1. **要做什么：** 不预设所有空间场都有明确结构锚点，先从分子空间数据发现候选场；对有独立 GT 的 TLS 与肿瘤—基质边界，再估计其周围的多通道距离响应，分别表示细胞组成、固定细胞类型内状态、可解释 program 和数据驱动生态信号。已知结构只作为候选场的验证与解释子集。
2. **服务的假设：** H-01、H-03、H-05。
3. **完成判据：** 至少一个候选空间场在患者或组织块层面显示可重复的空间残余信息；对有结构锚点的场，报告其与独立 GT 的对应关系；对未命名场，不要求先有 GT，但要用独立患者、队列或留出空间检验其可重复性，并给出组成与状态分量及不确定性。
4. **结果如何改变信心：** 若潜在场在独立队列方向一致，则提高对 H-01、H-03 和 H-05 的初始信心；若只有已知结构附近的复现而没有新的潜在场，结论收缩为已知结构场描述；若只存在组成变化或只在单队列出现，则削弱 H-03/H-05，并阻止进入强反演结论。

**R-04 内部科学顺序：**

1. 先从分子空间数据发现候选空间场，不用已知结构标签定义候选集合。
2. 再拆分细胞组成、细胞类型内状态和其他生态信号，判断候选场是否仍有残余信息。
3. 在独立患者或队列中检查候选场是否重复，避免把单一队列现象当成可迁移空间场。
4. 对 TLS 和肿瘤—基质边界等已有结构锚点，检验它们是否能解释或验证候选场的一部分；没有对应结构的场保留为潜在空间场。
5. 只有经过上述检验仍稳定的场，才进入 R-05 的捷径对抗和 R-06/R-07 的重建、定位评估。

**2026-08-07 工程状态：进行中。** 已建立 molecule-only 输入合同、双模型连续场适配器、按 factor 跨 section 的候选聚合、组成 cross-fitting、候选冻结、GT 隔离、checkpoint/hash 运行时和合成 smoke；一个 HTAN 单 section pilot 已执行但仍为 `FIT_COMPLETE_NOT_VALIDATED`。尚未完成正式多 section input manifest、五次 restart、独立 lineage 复现或 R-04 科学 gate，因此本节点仍未完成，不能据此提高或降低 H-01/H-03/H-05 的信心。

**2026-08-10 输入与冻结推断接口完成，科学节点仍进行中。** 已生成 71 个 eligible section 的显式 locator 和 molecule-only manifest，完成 71-section 预检并冻结 HTAN-only 4,000-gene universe；ST_CRC、USZ 和乳腺只承担预先冻结的验证角色。6 患者双模型 50-step CPU pilot 已完成并落盘，仍标记为 `FIT_COMPLETE_NOT_VALIDATED`；未执行 restart、patient/block bootstrap 或独立验证。训练后的 loading、gene-level nuisance/profile、输入 hash 和环境信息现可写入 `r04.frozen_model.v1`，新 section 只在固定模型下优化 section-specific GP 潜变量；冻结推断 CLI 已通过同一 pilot manifest 的 wiring smoke，但这不是独立验证。该结果解除 I-011 输入阻塞并提高 H-01/H-03/H-05 的可测性，不改变科学信心。
**2026-08-10 候选发现实现与门禁执行，科学节点仍未完成。** 本轮把 71-row manifest 拆成 role-pure partitions，training 的 4,000-gene panel 已物化为带坐标、full-library offset 和 hash 的稀疏 cache；fit/infer 现在 fail closed 检查角色，validation 缺失 gene 只允许显式 common-panel 条件推断，不做 zero-fill。双模型加入 gene-minibatch 的无偏 likelihood 缩放、未缩放的 global GP-KL、deterministic checkpoint step、factor effective-rank/energy/cosine 诊断；跨模型匹配改为 sign-invariant 的一对一 Hungarian conjunctive gate，仍不把场离散成 cluster。40 项 R-04 相关回归测试通过；三个 seed 的 20-step planted-factor 诊断仍显示 loading cosine 0.9964–0.9978 和有效秩约 1.00，三类不规则支持的 10-step calibration 也未恢复 planted fields，所有结果均保持 `CALIBRATION_COMPLETE_NOT_VALIDATED`；修正 offset 后的 47-section cache 单 mNSF step 约 69.5 秒、峰值 RSS 约 13.4 GiB。因此正式五次 restart 尚未启动，I-015 记为 CPU 时间硬阻塞，factor gate 也未通过；该结果不提高或降低 H-01/H-03/H-05 的科学信心。

**2026-08-13 目标函数与可识别性校准，科学门禁仍未通过。** 发现并修正 gene-minibatch 目标中把 spot-gene 联合均值直接乘以 `G/B` 的缩放错误；修正后 dense 与所有 3-gene batch 组合的 objective oracle 误差约为 `1.9e-6`，并新增连续场子空间、principal-angle、平台期和 `K_model/K_eff` 诊断。新增无场、单场、双独立场、共线场与 crescent/branch/disconnected 连续控制；正确指定双场控制在 300 steps 仍未平台，1000-step disconnected 代表性控制也未平台，三个几何的第二 canonical correlation 均未稳定超过 coordinate-permutation null。本轮定向回归测试共 54 项通过。一次 2-fold、单 seed、20-step 的 K=1/2/4 held-out 比较仅产生 K=4 探索性候选，因未达到平台且没有跨 seed/null 证据，最终 K 保持 `K_NOT_IDENTIFIABLE`。该结果说明校准仍不充分，不能判定真实场只有一个维度，也不能授权正式 restart；修正目标后的真实 47-section cache 单 mNSF step 为 93.63 秒、峰值 RSS 14,713,848 KiB，资源门禁仍保持关闭。

**2026-08-13 可识别性修复后科学校准通过，R-04 仍未完成。** 进一步发现并修正 inducing-value 映射、fit/infer objective 不一致、K 依赖初始化、合成 gene effect 共线和 K 评价泄漏。修正后，disconnected、crescent、branch 三种不规则支持的双独立场在三个 seed 均达到平台并恢复两个 spot×gene effect directions；单场和共线双源只支持一个方向。三 seed 的 section 留出 + gene cross-fit 在 planted `K_eff=2` 上比较 K=1/2/3/4，平均 held-out NB 分数约为 -29.38/-28.65/-29.10/-29.25，one-SE 唯一选择 K=2；K=4 不是有效最优候选，未扩到更大 K。无场控制比较 K=0/1/2，平均分数约为 -21.98/-22.14/-22.46，唯一选择 K=0。新版 18-control gate 为 `PASS_SCIENTIFIC_CALIBRATION`，说明方法选择规则在可控条件下可区分 0/1/2 个连续效应维度；不构成患者级发现，也不设定真实数据最终 K，因此 H-01/H-03/H-05 信心不变。修正参数化后的 47-section、4,000-gene、K=2 单步 CPU 用时 95.55 秒、峰值 RSS 16,277,456 KiB；五次 600-step mNSF restart 线性下限约 79.6 小时且未含 signed/validation，资源门禁仍为硬阻塞，正式 restart 未授权。
**2026-08-13 二次审计与证据替换，R-04 仍未完成。** 审计发现上一版 gene split seed 间接包含 K，导致不同 K 使用不同 evaluation genes；因此撤回 `k_calibration_aggregate_v1_20260813.json`、`k0_gate_aggregate_v1_20260813.json` 和 `scientific_gate_20260813_v5.json` 的判定用途。修正后，data seed、optimization seed 与 K-independent split seed 分离，聚合器强制同一 data seed 下的输入 hash、truth、执行参数、evaluation folds 和 K=0 baseline 完全一致，并以 data seed 而非同源 gene fold/section fold 作为独立标准误单位。重新完成 24 个 K=1/2/3/4 双场单元后，平均 held-out NB 分数为 -29.273/-28.688/-28.912/-29.177，三个 data seed 均以 K=2 最佳，K=4 支持比例为 0 且不在搜索边界，因此不扩到 K=5/6/8；12 个无场 K=0/1/2 单元的平均分数为 -21.975/-22.191/-22.524，三个 data seed 均选择 K=0。`scientific_gate_20260813_v7.json` 在 18 个控制、24 个正 K 单元和 12 个无场单元上同时通过嵌套来源哈希与底层语义重算，状态为 `PASS_SCIENTIFIC_CALIBRATION`。这只恢复 R-04 的方法可测性，不改变 H-01/H-03/H-05 信心；正式患者级 restart、组成拆分和独立验证仍未运行，资源硬阻塞不变。

**2026-08-13 独立审查修复与最终重跑，R-04 仍未完成。** 审查进一步发现：跨数值环境重放 truth 时浮点逐字相等会误阻断；mNSF 的 loss 在 Adam 更新前计算却保存更新后参数；非有限 canonical correlation 可绕过比较。现改为离散 truth 严格相等、浮点 truth 在 `1e-12` 容差内且必须有限；最佳参数在更新前写入 shadow checkpoint；总门禁对 objective、控制和 K shards 的所有必要数值做有限性与执行合同检查。修复后重新运行 18 个控制、24 个双场 K 单元和 12 个无场单元：K=1/2/3/4 平均 held-out NB 分数为 -29.273/-28.689/-28.912/-29.177，三个 data seed 均以 K=2 最佳，K=4 支持比例为 0；K=0/1/2 无场分数为 -21.975/-22.191/-22.523，三个 data seed 均选择 K=0。`scientific_gate_20260813_v10.json` 从底层来源重算后为 `PASS_SCIENTIFIC_CALIBRATION`。当前实现的 20-step 真实 training 资源测量显示 optimizer 为 92.80 秒（4.64 秒/步），计入 179.51 秒固定开销后，600-step 单 restart 规划值约 0.82 小时、五次 mNSF 约 4.12 小时；CPU 时间不再是硬阻塞。该结果只提高 R-04 的可执行性和可测性，仍不改变 H-01/H-03/H-05 信心；下一节点是实现并运行真实数据 K=0 向上、患者/组织块留出的非边界 K 粗筛，然后执行五次 restart。
**2026-08-13 真实 K 搜索入口开始执行，R-04 仍未完成。** 真实 K=0 已改为与正 K 共享 `nonspatial_rank=1`、offset、NB 目标、患者分组和 gene cross-fit，只移除坐标效应；新增入口固定 30 患者/47 section 的 grouped folds、K=0 向上候选和 checkpoint。缩小输入的一折 K=0/1、2-step smoke 已贯通 manifest/cache、留出推断和 scoring，但两个模型均未平台，结果保持 wiring-only，不能写成真实 K 选择。下一步仍是完整 5-fold K=0,1,2,3,4 粗筛；若最佳或 one-SE 选择触及上界，再扩 K=5,6。该工程结果不改变 H-01/H-03/H-05 信心。
**2026-08-13 完整真实粗筛已启动，R-04 仍未完成。** 完整 training 数据的 5-fold、患者/组织块留出、K=0–4 任务已在固定 CPU 环境启动，持久入口记录为 `RUNNING`；在聚合文件完成前，真实 K 仍为 `NOT_RUN`。空间置换 null、五次 restart、组成拆分和独立 lineage 复现尚未启动。该运行只服务 R-04/H-01/H-03 的可测性，不改变 H-05 的单 lineage 限制。
**2026-08-14 收敛门禁结果，R-04 仍未完成。** 真实 5-fold、K=0–4 粗筛已完成 25 个单元，30 个患者和 47 个 section 均有留出分数；但 25/25 个 fit platform 均未通过 600-step 平台期门禁，因此机器状态为 `K_NOT_IDENTIFIABLE_NOT_CONVERGED`，`selected_k=null`。K=3 只有原始预测分数最高这一筛选信号，不能写成真实最优 K，也不改变 H-01/H-03/H-05 信心。已保存 loss trace 后确认尾部仍有明显变化；独立持久 CPU 诊断正在以同一数据划分测试 K=0/2/3、fold 0、1200 steps，空间 null 和 restart 尚未启动。该结果对 plan.md 的影响是：H-01/H-03 的真实候选发现仍未判定，R-04 内部顺序保持不变。
**2026-08-16 低学习率诊断启动，R-04 仍未完成。** 先前 1200-step 哨兵中 K=0/2/3 的最佳步均接近末尾且 fit 未平台，因此新建独立持久 CPU 任务，使用完整输入、K=0/1/3、fold 0 和 fold 4、1200 steps、学习率 0.01。fold 0 代表原 K=0 轨迹的中位区域，fold 4 代表最差尾部；该任务只诊断共同优化协议，明确标记为 wiring-only，不改变患者折、gene cross-fit、NB 目标或收敛门槛。空间 null、restart、组成拆分和独立 lineage 复现仍未启动。

**2026-08-20 checkpoint 收敛审计启动，R-04 仍未完成。** 在再次改变 optimizer 前，先对既有 K=0/K=3 冻结 checkpoint 做完整 gene panel、固定 GP 抽样的 objective 重算，并记录同一 checkpoint 的重复误差及 baseline—nonspatial nuisance 的可能规范漂移。若固定 objective 已平台而在线 minibatch 门禁失败，先校准测量方式；若 objective 仍下降，再只测试一个经标准化更新量筛出的 nuisance 参数块。该审计不改变目标函数、不读取 GT、不选择 K；空间 null、restart、组成拆分和全量 25-cell 继续暂停。

**2026-08-20 checkpoint 审计完成：** K=0 的 dense objective 在 step 6600→7200 仅变化约 0.002%，且同一 checkpoint 的固定抽样重复误差为 0；原始 minibatch 平台门禁与固定 objective 不一致。K=3 在 step 3600→4800 仍下降约 3.1%，尚未平台。结果支持先把 deterministic objective 纳入训练监测，再延长 K=3 的同协议 continuation；不放宽原门禁、不解释 K、不启动全量搜索。

**2026-08-20 K=3 continuation 已启动：** 从固定 joint、fold 0、step 4800 checkpoint 继承 Adam state，保持 batch=512、lr=0.01、K=0 公平基线定义不变，续训至 step 7200；每 300 步记录一次固定随机数的 dense objective。首次提交因误用系统 Python 缺少 SciPy，未进入模型；已按此前成功审计的 Anaconda 环境重新启动，当前仍为 wiring-only，结果待定。

**2026-08-21 continuation 完成：** K=3 的 dense objective 从 5352.24 降至 5300.02；最后 600 步仅下降约 0.057%，但原始 fit 门禁仍为 false，inference 未运行。发现并修复顶层诊断 JSON 将 K 写死为 0 的元数据 bug，新增回归测试并通过 11 项 targeted tests。当前对 step 7200 checkpoint 做同一输入的重复审计；在该审计完成前不修改科学门禁。

**2026-08-21 确定性平台审计完成：** 用相同 dense objective、固定抽样和相同输入重算 step 6600 与 7200，目标从 5303.035 降至 5301.126（仅 0.036%）；step 7200 重复评估差为 0。确定性目标已满足“后段平台”的诊断证据，但原始 minibatch fit 门禁仍为 false。R-04 在“是否把确定性完整 objective 纳入正式 fit 收敛判据”处暂停；该选择会改变科学判定规则，不能由工程实现自动替代。

**2026-08-21 双门方案获用户批准：** 将 dense objective 平台作为优化稳定性的主要记录，将 minibatch 门禁保留为警告；两者同时写入结果，不互相覆盖。已补充计划定义、诊断函数和回归测试，并启动 K=0/K=3 fold-0 的 1200-step held-out inference。该阶段仍是诊断，不选择 K、不运行全量 25-cell。

**2026-08-21 首轮 inference 完成：** K=0 的 1200-step inference 通过原有平台门禁；K=3 的 inference 仍下降约 0.36%，其 fold-0 留出均值分数暂低于 K=0，因此不能区分“未收敛”与“真实劣势”。已启动 K=3 的 2400-step inference，并从 K=0 step 7200 补足至 7800 以形成同样的 dense-objective 平台窗口；随后按同一 2400-step 协议复跑 K=0。

**2026-08-23 fold-4 continuation 完成：** K=0 与 K=3 均从 fold-4 的 step 1200 checkpoint 继承 optimizer state 续训至 step 7200；两者 dense objective 均通过双向相对变化平台判定，原始 minibatch 门禁仍作为警告。两个 checkpoint 均为 wiring-only，尚未产生 fold-4 留出分数；该结果继续服务 R-04 的优化稳定性诊断，不改变 H-01/H-03/H-05 信心。

**2026-08-23 fold-4 inference 启动：** 已按与 fold-0 完全相同的 2400-step、gene cross-fit、GP 和 K=0 公平基线协议启动 K=0/K=3 inference；两个 CPU 持久任务均处于运行状态。结果落盘后先做 fold-0/fold-4 配对比较，再决定是否扩展到完整 5-fold；不提前选择 K 或启动空间 null、restart、组成拆分。

**2026-08-23 fold-4 inference 完成：** K=0 与 K=3 的 inference 均达到平台。fold-4 六名患者的均值分数为 K=0 `-2812.28`、K=3 `-2752.92`，K=3 在 6 名患者中的 5 名更高；但 fold-0 的对应差值为 `-28.04`，fold-4 为 `+59.36`，方向跨 fold 不一致。因此当前结果只证明两种模型都能在代表性 fold 上稳定完成诊断，不能选择 K=3 或提高 H-01/H-03/H-05 的科学信心。下一步需在相同协议下扩大独立 fold 证据，再决定是否进入完整 5-fold。

**2026-08-23 fold-1/2/3 当前协议重跑启动：** 由于库存中的旧 step-600 checkpoint 不属于当前双门训练协议，已对 fold 1、2、3 分别从 step 0 启动 K=0/K=3，使用 lr=0.01、gene batch=512、7200 fit steps 和 2400 inference steps。三个持久 CPU 任务均已进入运行；本轮只补足跨 fold 异质性证据，不执行 K 选择、空间 null、restart 或组成拆分。

**2026-08-24 fold-1/2/3 的 K=0 阶段完成：** 三个 fold 的 K=0 均完成 7200-step fit 和 2400-step inference；dense objective 与 inference 均达到平台，原始 minibatch fit 门禁仍为警告。三个任务随后进入 K=3：当前 checkpoint 为 fold1 step 1200、fold2 step 1800、fold3 step 1800。尚未汇总 K0/K3 的五-fold 科学比较。

**2026-08-24 K=3 continuation 进展：** 三个 fold 仍在同一持久 CPU 任务中运行，最新已保存 checkpoint 为 fold1 step 5400、fold2 step 5400、fold3 step 5400；尚未开始各自的 2400-step inference。该阶段仍只服务跨 fold 稳定性诊断。

**2026-08-25 五-fold K0/K3 诊断完成：** fold 1/2/3 的 K=3 fit 与 2400-step inference 已完成；加上 fold 0/4，10 个 K0/K3 诊断 cell 均完成，inference 平台全部通过，dense objective 按当前双向变化规则全部通过，原始 minibatch fit 门禁全部保留为警告。K3 在 30 名患者中 24 名高于 K0，患者加权平均差值为 `+22.80`；但 fold 0 的差值为 `-28.04`，其余四个 fold 为正，存在明确 fold 异质性。因此 R-04 仍未完成，`selected_k=null`，空间 null、restart 和正式 K 搜索均未运行。
**2026-08-25 极端 fold 重启诊断启动：** 为区分初始化敏感性与患者/fold 异质性，已在 fold 0（唯一负向 fold）和 fold 4（最强正向 fold）各启动一条独立 CPU 任务；每条任务依次从 step 0 重跑 K=0 与 K=3，并使用 `restart-index=1`、同一数据划分、gene batch=512、lr=0.01、7200-step fit 和 2400-step inference。该节点只服务 R-04/H-01/H-03 的稳定性判断，不选择 K；空间置换 null、正式 K 扩展和组成拆分继续暂停。
**2026-08-27 极端 fold 重启诊断完成：** restart 1 的四个 cell 均完成，dense objective 与 inference 平台均通过，gene-minibatch fit 门禁仍作为警告。fold 0 的 K3-K0 从 `-28.04` 变为 `-33.36`，fold 4 从 `+59.36` 变为 `+6.72`；两个 fold 的方向均保持，但效应幅度对初始化敏感。因此该结果支持 fold/患者异质性解释方向，不足以建立正式 seed 稳定性，R-04 仍不选择 K。
**2026-08-27 五 seed 稳定性面板启动：** 为补齐既定五次 restart，已启动 restart indices 2、3、4，覆盖 fold 0 和 fold 4 的 K=0/K=3；每条任务保持同一 split、gene cross-fit、lr=0.01、gene batch=512、7200-step fit 和 2400-step inference。该面板只判断极端 fold 内的初始化稳定性；空间置换 null、正式候选 K 比较和组成拆分继续暂停。
**2026-08-27 重启面板协议切换：** 原 restart 2–4 任务在完成 fold0/K=0 的部分 fit 后停止，部分 checkpoint 保留但不纳入科学汇总。为减少只读监测和写盘开销，在新目录从 step 0 重启相同的 restart indices；只调整 dense 监测间隔、诊断记录间隔和 checkpoint 间隔，不调整模型、目标、随机种子、fit/inference 步数或数据划分。线程池不做未经验证的限制；本轮仍只服务五 seed 稳定性判断。
**2026-08-31 五 seed 面板汇总完成：** canonical 与 restart 1-4 的 fold 0/4、K=0/3 共 20 个 cell 已完成审计汇总。fold 0 的 K3-K0 在五个 seed 均为负，均值范围 `-39.70` 至 `-23.56`；fold 4 均为正，范围 `+2.06` 至 `+59.36`，说明 fold 方向可重复但正向幅度对整体优化随机性敏感。optimization seed 同时控制初始化、gene-minibatch 顺序和变分采样，因此该面板不是纯初始化实验。跨 seed canonical correlation 显示第一方向较稳定、第二方向中等、第三方向很不稳定；无 collapse warning 不能替代有效秩证据。20/20 inference 达到平台；canonical 外置 audit 按当前规则重放后，canonical 与 restart 1 共 8/20 dense-objective 平台，restart 2-4 的 12/20 未平台。因此 R-04 保持 `K_UNDECIDED_OPTIMIZATION_NOT_RESOLVED`，不改变 H-01/H-03/H-05 的科学信心。
**2026-08-31 配对收敛补跑准备：** continuation 与 inference 工具已显式接入 `restart_index`，确保恢复 checkpoint 后仍使用原 optimization seed；12 项定向回归测试通过。下一步仅对 restart 2-4 的 12 个 cell 做 step 7200→8400 的 fit-only continuation，每 400 步评估 dense objective；只有配对 K0/K3 均通过后才重跑 inference 和评分。空间 null、K=1/2/4+、组成拆分和场命名继续暂停。

**2026-08-31 step-8400 配对审计完成并停止：** restart 2-4 的 12 个 K0/K3、fold 0/4 cell 均完成 7200→8400 fit-only continuation；冻结 checkpoint 的来源、分片、重复评估和 K0/K3 患者配对均通过 provenance 检查。用冻结的更新后 checkpoint 终点替换在线更新前终点后，10/12 cell 通过完整门禁。`restart 3 / K=3 / fold 0` 的 7600→8000 dense-objective 相对变化为 `0.1006526%`，严格高于 `0.1%`；`restart 2 / K=3 / fold 4` 的相同 checkpoint 两次评估有一个 draw 相差 `0.000488`，未达到逐值完全重复。依 D-081，状态为 `STOPPED_OFF_PLATFORM`，held-out inference、评分、K 扩展、空间 null、组成拆分和场命名均保持 `NOT_RUN`。该结果只确认优化/数值平台仍未完全解决，不改变 H-01/H-03/H-05 的科学信心。

**2026-08-31 严格审计与推理门禁复核：** 独立代码审查发现首版聚合未持久化 `ridge`、未逐项验证 repeat 数值结构，且未来 resume 路径理论上可在已审计终点后继续 fit；这些漏洞均已 fail-closed 修复。没有手工补写旧结果，而是在新目录按同一 CPU、checkpoint、seed 和四 draw 协议重跑 12 个 repeat，显式记录 `ridge=1e-5`。严格版 12/12 repeat 逐值一致，11/12 cell 通过；唯一正式失败仍是 `restart 3 / K=3 / fold 0` 的 `0.1006526%` 区间。首版 `restart 2 / K=3 / fold 4` 的单 draw 微小不一致保留为间歇性数值风险，不被删除或解释成已彻底解决。全局仍为 `STOPPED_OFF_PLATFORM`，不运行 inference；该工程与优化结果仍不改变 H-01/H-03/H-05 信心。

**2026-08-31 实际接受并启动冻结推断：** 用户依据唯一失败的边缘幅度、完整端点净变化和 12/12 精确 repeat，接受 step-8400 面板作为进入下一诊断阶段的实际充分结果；严格审计仍保留 `STOPPED_OFF_PLATFORM`、11/12，不改写为全通过。新增的独立授权绑定 strict audit、source、repeat、checkpoint bundle、输入和固定推断协议，只允许 12 个冻结端点运行 held-out inference 与 NB scoring，禁止继续 fit，保持 `selected_k=null`。三个 CPU worker 已分别启动 restart 2、3、4，每个 worker 顺序运行 fold 0/4 的 K0/K3；一次错误 locator manifest 启动在模型加载前失败，改为 47-section training role manifest 并验证 fold 0/4 输入 hash 后已正常重启。该节点服务 R-04/H-01/H-03/H-05 的优化误差与留出稳定性判断；结果完成前不改变这些假设的科学信心，不运行空间 null、K 扩展、组成拆分或场命名。

**2026-09-01 冻结推断面板完成并回收：** restart 2–4 的 12 个固定 step-8400 checkpoint cell 全部完成 2400-step、两次 gene cross-fit 的留出推断与评分，三份 panel 均为 `ALL_CELLS_COMPLETE`，失败 0；27 个输出文件已回收并逐项 SHA-256 一致。结果继续标记 `FROZEN_HELDOUT_DIAGNOSTIC_NOT_K_SELECTION`、`selected_k=null`，严格审计保持 `STOPPED_OFF_PLATFORM`；CPU/GPU 等价性对照按用户为避免租用 GPU 空转而中止，记为 `NOT_RUN`，不宣称跨设备数值等价。该结果完成 R-04 的冻结推断计算子阶段，但不完成 R-04 科学节点，不改变 H-01/H-03/H-05 信心。下一步按五 seed 聚合患者级分数、fold 方向和 top-1/2/3 子空间稳定性，再决定 K=2 或 K>4；空间 null、组成拆分和场命名仍不启动。

## R-05｜核心遮蔽与组成捷径对抗

1. **要做什么：** 对 R-04 的目标结构执行结构核心遮蔽、直接定义性 marker 限制、composition-only ceiling、同切片基础区室匹配、cross-fitted residualization，并在可用高分辨率数据中检验固定细胞类型后的状态变化。
2. **服务的假设：** H-01、H-03。
3. **完成判据：** 患者或组织块级比较明确报告完整信息相对 composition-only、marker-only 和空间平滑基线的增量；至少两个独立数据来源给出一致结论；匹配变量的因果角色已登记。
4. **结果如何改变信心：** 若完整信息在上述限制下仍有稳定增量，则第一次实质性支持 H-01/H-03；若增量消失，则否证 H-01 的当前强形式，并把项目收缩为结构—生态场描述或直接检测。

**为什么不能更早否证：** 在 R-01 与 R-02 完成前无法排除重复样本、同一物理结构跨折和 GT—输入循环；在 R-04 完成前也没有独立于测试集定义的候选场。更早得到的阳性或阴性都无法区分真实结构信息、组成捷径和数据泄漏。

## R-06｜概率三角定位与可识别性校准

1. **要做什么：** 把多个空间点的结构距离信息转化为结构位置概率密度，建立多通道、有符号和可非单调的响应表示；同时评估真实操作点的经验误差地板、核误设、源数量、空间自相关和后验覆盖。
2. **服务的假设：** H-02。
3. **完成判据：** 对每个进入定位任务的结构输出概率密度、可信区域、后验熵和不可识别区域；在模拟与真实留出数据中报告定位误差、覆盖率和相对最强单点或最近邻基线的增益。
4. **结果如何改变信心：** 若多点后验稳定收缩且校准优于基线，则支持 H-02；若后验持续弥散、错误多峰或只靠强先验收缩，则否证该结构在当前操作区间的定位 claim。

## R-07｜同平面隐藏结构定位 benchmark

1. **要做什么：** 在当前平面完全遮蔽目标结构本体和缓冲区，比较概率反演、composition-only、marker/program-only、morphology-only、空间平滑和其他适用基线。
2. **服务的假设：** H-01、H-02。
3. **完成判据：** 以患者或组织块为单位报告存在概率、位置误差、概率图覆盖和校准；模型增益必须大于 GT、配准和经验误差地板中的适用项。
4. **结果如何改变信心：** 若反演在严格遮蔽下仍稳定优于最强基线，则提高 H-01/H-02 信心并解锁切片内 hidden architecture 叙事；若只能在核心可见时有效，则该叙事被否证。

## R-08｜有结构锚点空间场的独立证据链 `[探索可并行]`

1. **要做什么：** 只有在需要提出第二个命名结构或更强通用性主张时，才从血管、坏死、肿瘤—基质边界或其他候选中选择有结构锚点的第二类空间场，并独立重复 R-04 至 R-07 的证据链；潜在空间场的发现不依赖本节点。
2. **服务的假设：** H-04、H-05。
3. **完成判据：** 第二结构具有独立患者或组织块测试、独立 GT 审计、结构级性能和适用词汇；若没有候选达到要求，明确记录失败而不替换成事后挑选的“漂亮结构”。
4. **结果如何改变信心：** 若第二种机制或几何不同的结构成立，则支持 H-04 并提高方法级通用性；若所有命名结构候选失败，则收缩结构通用性主张，但不自动否定未命名潜在空间场的发现结果。

## R-09｜重叠结构的多源可分解性 `[探索可并行]`

1. **要做什么：** 在人工控制和真实重叠区域中比较单一共同场、独立单源模型和联合多源后验，量化结构对之间的残余场相似性和实例归属。
2. **服务的假设：** H-04。
3. **完成判据：** 对每个结构对报告是否可分、只能联合报告或完全不可识别；联合模型必须在组织块留出中改善实例归属或概率校准。
4. **结果如何改变信心：** 若至少一组重叠结构可稳定分解，则增强 H-04；若所有结构对均只能归为共同场，则否证多源分解 claim，但可保留联合生态场输出。

## R-10｜跨队列、癌种与平台迁移

1. **要做什么：** 对 R-04 至 R-09 中成立的结构场执行 leave-one-cohort-out、leave-one-cancer-out 和跨平台评估，分别比较零样本迁移、少量目标域适配和目标域重训；确切队列组合为 `[待定]`，需要 R-01 的去重与角色冻结。
2. **服务的假设：** H-05。
3. **完成判据：** 报告共享成分、上下文特异偏移、迁移预算、校准变化和失败域；外部测试数据不参与 program、核或阈值选择。
4. **结果如何改变信心：** 若至少部分结构场在独立队列或平台保留方向和定位增益，则支持 H-05；若全部性能只存在于原队列，则否证跨域共享主张并把结论限制在队列内。

## R-11｜邻近平面与三维竞争 benchmark

1. **要做什么：** 根据 R-03 可用性，在同平面 core-masked、邻近平面 out-of-plane、抽中间平面或栈外外推中选择可测任务，并在相同输入预算下比较连续性、H&E、组成、UniST/NTF 式插值或其他适用三维基线。
2. **服务的假设：** H-06。
3. **完成判据：** 当前输入平面不包含目标结构核心；组织块级留出；模型增益超过最强基线和配准/GT 误差地板；概率区间达到预先冻结的经验覆盖要求。具体可执行子任务与阈值为 `[待定]`，需要 R-03 的 serial-section 普查结果。
4. **结果如何改变信心：** 若分子场反演提供稳定增量，则支持 H-06 并允许相应的 nearby 3D 或 out-of-plane 词汇；若没有增量，则否证三维扩展，但不影响已经成立的同平面反演。

## R-12｜主动观测的回顾性删观测评估

1. **要做什么：** 在已有连续切片或多模态数据中隐藏部分观测，根据 R-06/R-11 的不确定性选择下一切片、ROI 或验证 marker，并与固定间距、随机和可获得的专家选择对照。
2. **服务的假设：** H-07。
3. **完成判据：** 以新增确认结构数、后验熵下降或单位成本信息增益为指标；选择规则在独立组织块评估，并报告校准依赖。
4. **结果如何改变信心：** 若回顾性选择稳定优于对照，则增加 H-07 的初步信心；若不优于对照，则削弱 H-07，并停止前瞻投入。

## R-13｜主动观测的前瞻小队列 `[不可逆投入]`

1. **要做什么：** 在新的组织块上按模型建议选择切片、ROI 或验证 marker，并与预先规定的对照策略比较。样本数、伦理路径、实验平台和预算均为 `[待定]`，需要可用组织来源、成本和审批信息。
2. **服务的假设：** H-07。
3. **完成判据：** 前瞻结果能够计算单位实验投入的信息增益或结构确认率，并且选择和终点在获取结果前冻结。
4. **结果如何改变信心：** 若前瞻选择提高单位成本信息量，则支持 H-07 和 active sampling 叙事；若无改善，则否证其转化价值，不以回顾性结果替代。

**投入性质：** 样本消耗、染色和高分辨率测量一旦发生难以撤回；在 R-12 未通过前不得进入本节点。

## R-14｜表示与模型复杂度条件探索 `[探索可并行]`

1. **要做什么：** 仅当 R-07、R-09 或 R-10 暴露明确失败模式时，比较线性 program/NMF/topic/ICA 基线与图模型、masked spatial model、隐式场或其他 `[待定]` 非线性表示；候选由失败模式决定，不按流行度决定。
2. **服务的假设：** H-02、H-04、H-05。
3. **完成判据：** 新表示在患者或组织块外层留出中解决预先指定的失败模式，并报告相对统计基线的净增益、校准和额外复杂度。
4. **结果如何改变信心：** 若复杂表示恢复了可识别性、可分解性或迁移性，则提高对应假设信心；若没有稳定增益，则保留简单基线并停止模型扩张。

## R-15｜假设覆盖检查

- H-01：R-04、R-05、R-07。
- H-02：R-06、R-07、R-14。
- H-03：R-04、R-05。
- H-04：R-08、R-09、R-14。
- H-05：R-04、R-08、R-10、R-14。
- H-06：R-11。
- H-07：R-12、R-13。

没有未被非基础设施节点引用的 plan.md 假设。
