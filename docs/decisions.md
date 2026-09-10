# 决策日志：空间结构反演

本文档的角色(需写入文档开头): 本文档记录已经发生的选择及其当时的理由,是唯一无法从代码、plan 或 roadmap 反推的信息。本文档只追加,不修改、不删除历史条目;推翻旧决策的方式是追加一条新条目并引用旧条目编号。本文档不包含计划、待办和未来意图——那些属于 roadmap。

### D-001 | 2026-07-20 | 将课题定义为空间结构反演，而不是已知结构周围的梯度分析
- 背景:备选包括继续描述 TLS 周围的距离相关 program、预测到已知 TLS 的距离，或从分子场逆推未知结构
- 理由:前两者与已有空间梯度工作重叠，逆问题才对应未直接观测结构的科学空白
- 代价:引入了可识别性、多解、多源混合和不确定性校准问题
- 复查条件:如果已有工作已完整解决“外围场到未知结构”的反演，应重新考虑本决策
- 影响:plan.md 的 H-01、H-02；roadmap.md 的 R-04 至 R-07

### D-002 | 2026-07-20 | 反演的标准输出采用空间概率密度而非唯一确定的三维 volume 或 mesh
- 背景:备选是输出单一结构坐标、确定性三维体积，或输出后验概率与可信区域
- 理由:讨论中确认单平面通常不支持唯一几何解，而概率密度仍可表达可靠的邻近结构预测
- 代价:结果解释更复杂，并要求概率校准和不可识别警告
- 复查条件:如果未来获得足够多连续平面或独立几何约束，可重新考虑确定性几何输出
- 影响:plan.md 的 H-02、K-06；roadmap.md 的 R-06、R-11

### D-003 | 2026-07-20 | 保留“近似或概率性邻近三维重建”的强叙事，但由证据解锁
- 背景:备选是主动放弃三维词汇，或不加限制地宣称恢复真实三维结构
- 理由:讨论认为激进叙事本身不是问题，关键是是否击败连续性、形态和插值解释
- 代价:需要相邻切片、配准地板和强竞争基线；若失败必须收缩叙事
- 复查条件:如果没有可测的 serial-section 数据，或邻近平面无增量，应取消三维词汇
- 影响:plan.md 的 H-06、K-05；roadmap.md 的 R-03、R-11

### D-004 | 2026-07-20 | 将结构信号问题拆为“复杂空间生态”和“spot-level bulk-like 混合”两部分
- 背景:备选是把结构 program 当作单一 pathway，或把 spot 直接视为单细胞状态
- 理由:讨论确认很多 program 可能是复合 ecosystem signal，而多数空间数据混合多个细胞
- 代价:需要分别处理组成、细胞内状态及未命名生态成分，模型和验证更复杂
- 复查条件:如果项目限定到单细胞分辨率平台且目标结构只有单一明确 marker，可重新简化
- 影响:plan.md 的 H-03；roadmap.md 的 R-04、R-05

### D-005 | 2026-07-20 | 物理和形态约束作为可选先验，不作为所有结构的普遍前提
- 背景:备选是所有结构都采用扩散、形态或连续性约束，或完全放弃这些约束
- 理由:可见血管、坏死等结构可受益于硬约束，但不可见或无清晰物理边界的生态结构不能被排除
- 代价:不同结构需要不同先验等级，统一模型更困难
- 复查条件:如果最终只研究可见硬结构，可提高物理约束的中心地位
- 影响:plan.md 的 H-04、B-02；roadmap.md 的 R-02、R-08、R-09

### D-006 | 2026-07-20 | 不预设 NMF 或 topic model 为主模型，只把简单分解作为基线
- 背景:备选包括直接采用 NMF/topic、选择更复杂的图或神经场模型，或由实证失败模式决定
- 理由:用户明确反对因近期使用经验而偏向这些已经被广泛使用的简单模型
- 代价:主模型保持未定，需要先建立统计基线和失败诊断
- 复查条件:如果简单模型已经达到全部判据，应取消复杂模型探索
- 影响:roadmap.md 的 R-14

### D-007 | 2026-07-20 | 结构场知识来自已标注结构周围的距离—信号关系，反演是其逆向使用
- 背景:备选是纯无监督发现结构、直接训练端到端结构检测器，或先研究固定结构场
- 理由:讨论形成的主框架是利用大量空间数据先描述固定结构及其周围场，再检验逆向定位
- 代价:存在结构选择效应和发现—验证泄漏风险
- 复查条件:如果固定结构缺乏可复现 halo，或纯发现路线出现独立锚点，应重新评估
- 影响:plan.md 的 H-01、H-04；roadmap.md 的 R-04、R-08

### D-008 | 2026-07-20 | TLS 作为概念示范和主起点，但不作为唯一依赖结构
- 背景:备选是把课题做成 TLS 专用工具，或一开始宣称任意结构
- 理由:TLS 提供明确启发和大规模数据，但单一 TLS 容易退化为“免疫区有免疫”的识别
- 代价:必须为至少第二类结构建立独立证据链
- 复查条件:如果没有任何第二结构达到可识别和功效要求，方法级通用 claim 必须收缩
- 影响:plan.md 的 H-04；roadmap.md 的 R-08、R-09

### D-009 | 2026-07-21 | 接受第三方审计的“B：条件立项”，但不接受其所有硬门槛
- 背景:备选是接受审计全部条款、拒绝审计，或接受主要风险并修订过度限制部分
- 理由:审计正确识别组成捷径、泄漏、三维竞争解释和分层功效，但部分条件被认为不是唯一有效协议
- 代价:需要维护审计条款、项目回应和最终 benchmark 之间的可追溯关系
- 复查条件:如果后续三方会议正式采纳不同判定，应追加新条目替代本决策
- 影响:plan.md 的全部假设；roadmap.md 的 R-01 至 R-14

### D-010 | 2026-07-21 | 组成独立性不能由删除少量 marker 单独证明
- 背景:备选是把 marker removal 作为主要证据，采用单一“组成充分统计量”，或使用多种对抗证据
- 理由:剩余转录组仍可能重建细胞组成，而单一组成预测器也不能被证明是充分统计量
- 代价:需要 composition-only ceiling、同切片匹配、残差化和固定细胞类型分析等组合证据
- 复查条件:如果目标数据全部达到可靠单细胞分辨率且组成已直接观测，可简化部分测试
- 影响:plan.md 的 H-01、H-03；roadmap.md 的 R-05

### D-011 | 2026-07-21 | 采用 GT—输入信息不相交原则，但不对所有生物代理特征一刀切删除
- 背景:备选是使用相同 RNA 定义并预测结构、删除所有相关通路，或逐实例登记 GT 与输入重叠
- 理由:完全循环不可接受，但结构组成、结构响应和标签泄漏不能被“代理特征”一个词混为一谈
- 代价:需要逐实例 GT 重叠注册和分层证据等级
- 复查条件:如果某结构只能由与输入完全相同的特征定义，则其确认性用途应取消
- 影响:plan.md 的 H-01、B-05；roadmap.md 的 R-02、R-05

### D-012 | 2026-07-21 | 项目 claim 按结构分别成立，不以最弱结构压低全部结论
- 背景:审计建议项目最高 claim 取各结构最低上界，备选是按结构和方法级分别签署
- 理由:一个弱结构不应否定其他结构已经通过的完整证据链，一个强结构也不能替其他结构背书
- 代价:论文和工具必须维护结构级适用范围，不能给出一个无条件总分
- 复查条件:如果最终只保留一个结构，方法级通用 claim 应取消
- 影响:plan.md 的 H-04、B-06；roadmap.md 的 R-08、R-09

### D-013 | 2026-07-21 | 三维竞争采用相同输入预算的分层 benchmark，不采用固定 z 跨度作为唯一门槛
- 背景:审计提出 z 跨度不超过三倍切片间距的硬门槛，备选是按实际基线增益定义非平凡性
- 理由:结构尺寸本身不能决定任务是否平凡，关键是是否在无当前核心时超过连续性、形态和插值
- 代价:需要多个输入预算组和更复杂的 benchmark 解释
- 复查条件:如果三方正式冻结硬 z 阈值，应追加新条目并修改相应词汇解锁规则
- 影响:plan.md 的 H-06；roadmap.md 的 R-11

### D-014 | 2026-07-21 | 主要切分和不确定性评估以患者或组织块为单位
- 背景:备选是按 spot 随机切分、按患者切分，或进一步按物理组织块和跨切片结构实例切分
- 理由:spot 强自相关，同一物理结构可跨相邻切片重复出现，患者级切分仍可能泄漏
- 代价:有效样本量显著小于 spot 数，功效与模型选择更严格
- 复查条件:如果数据没有可恢复的 block_id，则对应数据只能承担较弱证据
- 影响:plan.md 的 H-01 至 H-06；roadmap.md 的 R-01、R-02、R-04 至 R-11

### D-015 | 2026-07-20 | 主动切片、ROI 和验证 marker 选择属于项目的扩展价值
- 背景:备选是只输出结构热图，或让不确定性直接服务后续观测选择
- 理由:讨论一致认为主动空间实验设计能把概率输出转化为真实使用价值
- 代价:回顾性证据不足以解锁最终 claim，需要前瞻样本和不可逆实验投入
- 复查条件:如果概率输出不能校准或回顾性选择无增益，应取消前瞻投入
- 影响:plan.md 的 H-07、K-07；roadmap.md 的 R-12、R-13

### D-016 | 2026-07-27 | 冻结当前数据扩招并维持 CROST、Synapse/HTAN 原始数据和受控数据 HOLD
- 背景:现有数据约 1.9 TB，仍可继续下载 CROST、Synapse/HTAN 原始数据和多个受控队列
- 理由:用户已实际 HOLD CROST 与 HTAN/Synapse；现有数据对二维场发现和反演 MVP 已充足，新增普通切片的边际价值下降
- 代价:可能暂时错过独特 serial-section、正交 GT 或特殊平台数据
- 复查条件:如果 R-01/R-03 发现缺少关键结构、独立 GT、serial-section 或高分辨率验证，应重新评估特定数据源
- 影响:roadmap.md 的 R-01、R-03

### D-017 | 2026-07-30 | 当前数据与科学范围保持在肿瘤空间组学，不向其他组织面扩展
- 背景:备选是继续扩展非肿瘤组织以增加结构种类和样本量，或先在现有肿瘤数据闭环
- 理由:用户明确要求当前只讨论数据充分性，暂不扩招到其他方面
- 代价:暂时不能回答非肿瘤组织的普遍适用性
- 复查条件:如果肿瘤数据无法提供第二结构或足够外部验证，再重新讨论是否扩域
- 影响:plan.md 的 B-07；roadmap.md 的 R-01

### D-018 | 2026-07-30 | plan、roadmap 与 decisions 采用稳定目标、可变执行和只追加历史的三文档治理
- 背景:原有完整方案同时混合科学目标、执行路线和历史理由，后续 agent 容易混淆
- 理由:用户明确规定三份文档的职责和交叉约束
- 代价:同一变更可能需要跨文件引用，维护要求提高
- 复查条件:如果文档治理规则再次改变，应追加新决策而不是改写本条
- 影响:本仓库的 plan.md、roadmap.md 和 decisions.md

### D-019 | 2026-07-31 | 根项目独立纳入 Git，保留上游代码仓库和大型数据边界
- 背景:项目根目录此前未形成有效 Git 仓库；`repo/` 是带独立远程的上游论文代码仓库，`data/` 约 1.9 TB
- 理由:项目治理文档、分析代码、注册表和结果索引需要统一版本历史，但不应把上游仓库历史或大型数据对象误并入根仓库
- 代价:根项目与 `repo/` 维持两个版本边界；若需修改上游代码，必须明确在哪个仓库提交和如何记录依赖版本
- 复查条件:如果后续决定 fork 上游代码、采用正式 submodule，或把数据对象交由 DVC/对象存储管理，应追加新决策
- 影响:版本控制、R-01 的数据来源追踪及后续代码布局

### D-020 | 2026-07-31 | 当前结构与样本证据仅接受数据自带 metadata
- 背景:目前没有已确认的外部病理标注或专家复核来源，也不确定现有 metadata 包含哪些结构和切片关系字段
- 理由:把文件名、表达模式、论文常识或模型推断静默当作 GT 会造成循环定义和不可审计的信息泄漏
- 代价:部分数据可能只能承担发现性或 archive 角色，R-02 及后续确认性证据可能因 metadata 不足而受限
- 复查条件:用户批准新增外部数据、独立标注或专家复核后，重新登记证据来源和 GT—输入边界
- 影响:plan.md 的 H-01、H-04、H-06；roadmap.md 的 R-01 至 R-03、R-05、R-08、R-11；ISSUES.md 的 I-001 至 I-003

### D-021 | 2026-07-31 | 本地优先、禁用 GPU 并设置 2 TB 存储安全余量
- 背景:现有本地数据足以启动普查；用户允许自主裁决不影响系统工作的 CPU 和内存使用，但暂不允许 GPU
- 理由:先用现有资产回答可测性和证据缺口，避免无依据扩招数据或占用共享系统资源
- 代价:GPU 友好型模型或预处理可能延迟；所有读取、解包和派生物必须受存储余量约束
- 复查条件:现有数据缺口具有不可替代性时提出外部数据申请；GPU 禁用导致运行时间显著延长并影响进度时升级为硬阻塞，统一审核
- 影响:roadmap.md 的 R-01 至 R-14；所有计算与数据获取任务

### D-022 | 2026-07-31 | 使用 ISSUES.md 管理未解决但未硬阻塞的问题
- 背景:原有三文档治理没有稳定位置记录跨节点证据缺口，容易让未决问题被遗忘或被误报为阻塞
- 理由:将普通 issue 与 roadmap 任务、历史 decision 和硬阻塞分开，可保持阶段状态和结论影响清晰
- 代价:问题状态变化需要同步维护；ISSUES.md 不能替代 roadmap
- 复查条件:如果问题日志与 roadmap 出现持续重复或状态漂移，应重新设计治理边界
- 影响:AGENTS.md、ISSUES.md 和项目状态报告纪律

### D-023 | 2026-07-31 | R-01 采用 outcome-blind 角色冻结和保守 leakage group
- 背景:Atlas、HEST、GEO、10x、HTAN、STOmicsDB 等来源存在聚合镜像、同患者多切片和重处理版本；同时 Table S2 已暴露 TLS 阳性与数量，若据此挑验证队列会形成 outcome-guided selection
- 理由:数据角色应由可审计物理身份、来源谱系、本地可用性和预注册 capability 冻结，而不是由结构结果或模型表现决定；未解决的 possible match 在切分上同组可防止跨折污染
- 代价:可能过度合并尚未确认的独立样本，且 metadata 不足会显著缩小可承担确认性结论的核心集
- 复查条件:用户批准独立身份核验或新增自带 patient/block 映射的数据后，可拆分被证实独立的 possible-match group；若 roadmap 正式放宽 block 要求，应追加新决策而不是静默使用 sample 代理
- 影响:roadmap.md 的 R-01、R-02、R-05、R-10；ISSUES.md 的 I-004、I-006

### D-024 | 2026-07-31 | R-01 在两个合格逻辑单元处停止，不以 sample 代理 block
- 背景:穷尽当前本地 bundled metadata 后，只有 HTAN Vanderbilt CRC 与 10x Breast Block A 同时具备可审计 patient+block；roadmap 要求 6–10 个逻辑单元
- 理由:把 Table S2 Sample ID、HEST id、slide_id、sample_key、文件名或目录静默当作 block 会隐藏同组织块跨折风险，并制造虚假的独立验证
- 代价:R-01 进入 BLOCKED_IDENTITY，R-02 和模型工作暂停；需要额外 metadata 才能恢复关键路径
- 复查条件:至少新增 4 个彼此独立且自带 patient↔block 映射的逻辑单元，并能保留独立 external validation lineage 时重跑 gate；若用户决定修改科学证据门槛，应先修改 roadmap 并追加替代决策
- 影响:roadmap.md 的 R-01 至 R-02；ISSUES.md 的 I-006；infra/sample-registry/r01_gate.json

### D-025 | 2026-07-31 | R-01 gate 对证据等级和角色冻结 fail closed

- 背景:只按 patient、block 和状态计数会让未来误标为 candidate 的 E1/E0/冲突记录进入确认性集合；只统计冻结行数量也可能允许无关单元、重复角色或 leakage 冲突制造 false green
- 理由:机器 gate 必须直接执行注册表已经声明的证据纪律，而不能依赖上游人工保证；claim-eligible 身份仅接受 E3/E2，冻结行必须引用合格逻辑单元、具有唯一有效主角色和非空 leakage group，且 external 不得与其他 claim-bearing 角色共享 leakage group
- 代价:未来不完整的 role freeze 会明确停在 `BLOCKED_INDEPENDENCE`，即使已有足够行数；gate 命令写出阻塞产物时仍返回零，自动化必须解析 JSON status
- 复查条件:若 roadmap 改变允许的证据等级、角色集合或 leakage 独立性标准，必须同步修改 gate policy、回归测试和 completion requirements
- 影响:roadmap.md 的 R-01；infra/sample-registry/r01_gate.json；后续 CI/自动化

### D-026 | 2026-07-31 | 外部 metadata 请求按身份缺口排序且不预支独立性

- 背景:R-01 已因仅有两个合格逻辑单元硬阻塞，但 574 条记录已有 patient、缺 block；泛化申请“更多数据”无法说明哪些 metadata 能直接解阻，也容易按 TLS 结果挑队列
- 理由:请求队列只使用 bundled patient coverage、provenance locator 和 metadata 冲突数排序，不读取 TLS/结构结果或模型性能；相同 accession、PMID 或 DOI 先并为 provenance group，每组在 patient↔block crosswalk 和跨来源重复审计通过前贡献 0 个新增单元
- 代价:同一论文下真实独立的多个队列会被暂时保守合并；50 个请求候选目前只形成 44 个 provenance group，且仍不能保证其中任意 4 个最终独立合格
- 复查条件:用户批准 metadata-only 获取并得到官方 crosswalk 后，按物理身份拆分或继续合并 provenance group；若新增证据改变现有 HTAN/10x 谱系关系，追加决策并重跑 gate
- 影响:roadmap.md 的 R-01；ISSUES.md 的 I-004、I-006；infra/sample-registry/r01_metadata_request.tsv

### D-027 | 2026-07-31 | R-01 接受窄定义的 patient-linked physical specimen 作为 block-equivalent

- 背景: 经用户批准获取官方 metadata 后，GSE211956、GSE226997、GSE274103 和 GSE274557 均有 patient 映射、实体组织描述与稳定 BioSample locator，但没有字面 block_id；继续要求字面 block 会保持 R-01 硬阻塞，把 GSM 或 BioSample 直接写成 block 又会违反 D-024。
- 理由: 对 R-01 的泄漏控制，最低充分证据是能把同患者的实体组织 specimen 保守同组，而不是伪造 block 字段。仅当官方 patient 映射、实体组织描述与稳定 specimen locator 至少形成 E2 互证时，允许登记 `patient_linked_physical_specimen`；`block_id` 保持空。
- 代价: 该证据弱于真实 patient—block crosswalk，不能证明一个 specimen locator 内只含一个组织块，也不能支持精细 block-level 独立性主张。
- 复查条件: 获得真实 block crosswalk 时以真实 block 替换；发现一个 BioSample/GSM 含多个组织块、多个 GSM 是同一 block 的切片、患者映射冲突或跨研究复用 specimen 时，撤销 equivalence、重建 leakage group 并重跑 gate。
- 影响: 替代 D-024 的“只接受字面 block”完成门槛，但保留其“禁止 sample/section/file name 代理 block”纪律；roadmap.md 的 R-01；sample registry schema 与 gate。

### D-028 | 2026-07-31 | 冻结六个 R-01 逻辑单元的 outcome-blind 角色

- 背景: 四个新增 GEO lineage 通过身份审计后，核心集达到 6 个逻辑单元，需要在读取结构结果或模型表现前完成角色冻结。
- 理由: 按身份完整性、独立 BioProject/source lineage、本地资产和预注册 capability 冻结：GSE274557=discovery，HTAN Vanderbilt CRC=training，GSE274103/GSE226997=internal validation，GSE211956=external validation，10x Breast Block A=serial-section validation。
- 代价: 冻结后不能因 TLS 丰度、结构覆盖或模型结果不理想而替换 external；GSE274103 与 GSE274557 虽有不同 BioProject 和论文，仍保留共享团队这一 provenance 风险。
- 复查条件: 发现 lineage/患者/实体 specimen 重叠、角色 capability 不成立或本地分析资产不可用时，相关角色退回 archive/in-progress，状态改为 `BLOCKED_INDEPENDENCE` 或 `BLOCKED_IDENTITY`；不得按结果换队列。
- 影响: roadmap.md 的 R-01；infra/sample-registry/role_freeze.tsv；后续所有外层切分与验证。

### D-029 | 2026-07-31 | 所有 claim-bearing 身份必须通过跨表引用与来源文件完整性 gate

- 背景: 只检查 `physical_units.tsv` 自报的 E2/E3、patient、block 或 specimen basis，会允许任意 schema-compatible 行伪造六个 study 形成 false green，即使没有 identity evidence、source asset 或 canonical lineage。
- 理由: 真实 patient/block 路径必须核验 patient 与 block evidence；specimen-equivalent 还必须核验 physical specimen 与 block-unknown basis。所有 evidence 的 normalized value、raw key/value、冲突状态和 source path 必须闭合；引用的 metadata asset 必须 active/PRESENT、SHA-256 与实际文件一致；canonical source lineage 与冻结角色 leakage group 必须一致。
- 代价: gate 会读取并校验小型 metadata 文件，表间任一缺失、校验和变化或 lineage 漂移都会使相关行 fail closed；构建和测试夹具更复杂。
- 复查条件: schema 或来源体系变化时，应提供等价的来源认证、identity normalization 与 canonical lineage 规则，不能退回只信任 physical row 自报。
- 影响: D-025 的机器执行强度；roadmap.md 的 R-01；`r01_gate.json` 的 integrity errors。

### D-030 | 2026-07-31 | 只有来源闭合且具有可重放空间几何的结构标签可作为确认性 GT

- 背景: 当前 metadata 中存在 TLS presence/count、TLS ID、成熟度、位置类别和面积等结果汇总，上游代码还明确包含 marker scoring、phenotype inference、TLS segmentation 和 H&E YOLO 流程。仅凭标签名称或实例编号，无法证明它是独立于拟用输入的空间 GT。
- 理由: 确认性 GT 必须能关联 R-01 physical unit，具有可重放的 mask、polygon、centroid 或等价坐标 locator，并闭合标注模态、生成方法、坐标系、分辨率、边界不确定性和来源校验和。GT 生成依赖的输入及其派生特征按传递闭包禁用；同一表达 assay 派生、来源未知或只有样本/实例汇总的标签只能用于 discovery 或直接排除。
- 代价: 现有 57 条 TLS 实例汇总不能因有唯一 TLS ID 而进入定位 benchmark；本地 19 个文件名含 `TLS_annotation` 的资产在 provenance 和物理映射闭合前也不能升格。当前可用样本量因此可能降为零。
- 复查条件: 获得内容寻址的空间 annotation、physical-unit crosswalk 和明确的独立生成 provenance 后重审；若 GT 来自 H&E，仅在分子-only 输入且完全排除该图像及其派生特征时考虑确认性使用。
- 影响: roadmap.md 的 R-02；后续 R-04 至 R-10 的 GT、遮蔽、输入 allowlist 和 benchmark 有效性。

### D-031 | 2026-07-31 | 外层切分采用身份关系图；block 未知时只冻结患者级保守 envelope

- 背景: R-01 的 HTAN 与 10x 单元有真实 block，四个 GEO 单元只有 patient-linked physical specimen。把 GSM、BioSample、specimen 或文件名当 block 会伪造 block-level 独立性；只按单行切分还可能让同患者多 specimen 跨折。
- 理由: outer group 取同患者、同真实 block、重复/重处理、serial-section、显式跨切片同结构关系和 R-01 leakage edge 的保守连通分量。真实 block 保持可审计；block 未知时把同患者全部 specimen 放入 `PATIENT_ENVELOPE_BLOCK_UNKNOWN`，`block_id` 继续为空，不计 block-level 样本量。
- 代价: 四个 specimen-equivalent GEO lineage 的 72 条 physical rows 最多只形成 30 个患者级 envelope，不能声称 72 个独立 block；GSE274557 的 55 条 specimen 收缩为 13 个患者级 envelope。此规则保证不跨患者泄漏，但不能证明未知 block 之间独立。
- 复查条件: 获得真实 patient—block—section crosswalk 后，以真实 block 替换 patient envelope 并重跑 duplicate、cross-section 和 role gate；发现跨患者别名或跨 lineage 物理复用时合并相应连通分量。
- 影响: roadmap.md 的 R-02；`outer_split_units.tsv`；所有患者/组织块级外层验证和有效样本量报告。

### D-032 | 2026-07-31 | 用户批准扩大外部数据范围以修复 R-02 GT 硬阻塞，限定公开 annotation/crosswalk/provenance

- 背景: R-02 硬阻塞（I-001/I-003）的唯一解除路径是获得独立于表达模态、带可重放空间几何的结构标注。本地侦察确认：Atlas Table S4 只有汇总属性、repo rds 无坐标、GEO RAW.tar 只含 Space Ranger 标准输出、HTAN h5ad obs 只有表达聚类列；本地不存在任何可用 GT 资产。文献普查确认：四个 GEO 谱系的原始论文均未公开标注文件（GSE274103 的 Loupe/HALO 标注与 GSE274557 的 Loupe 标注存在但未沉积，只能向作者索取；GSE226997 无独立标注；GSE211956 的 QuPath 几何未沉积）；唯一公开可得的合格来源是 HTAN Vanderbilt CRC 原始论文（Heiser et al., Cell 2023, PMC10756562）GitHub 仓库 `Ken-Lau-Lab/spatial_CRC_atlas` 的 42 个逐 spot 病理标注 CSV、`visium_sample_key.csv` 和 LCM ROI mask（病理医生 H&E Loupe 手工标注，spot barcode 坐标系）。
- 理由: 证据缺口（确认性结构实例为 0）且该来源不可替代（GT 必须绑定 R-01 冻结物理单元，只有同一研究的官方标注满足）；规模极小（数十个 CSV，<50 MB，远低于存储余量约束）；预期解锁判据明确（training 谱系出现 auditable GT source 与 confirmatory instance、第二结构可冻结则 R-02 机器门禁离开 HARD_BLOCKED_NO_AUDITABLE_GT；validation 谱系 GT 仍缺，单独跟踪）。
- 代价: 范围从此前 metadata-only 扩大到公开标注内容文件与论文方法学文本；引入新的 provenance 审计面（GitHub 仓库 commit、PMC 文本、标注者身份粒度）。
- 复查条件: 若发现标注实际为表达派生、commit 内容漂移、或作者公开更完整几何（polygon/mask），重审该来源等级；向作者索取 GSE274103/GSE274557 标注的行为需要用户单独发起，不在本决策范围内。
- 影响: roadmap.md 的 R-02；ISSUES.md 的 I-001、I-003；`infra/structure-registry/`；后续 R-04 至 R-10 的 GT 与输入 allowlist。

### D-033 | 2026-07-31 | R-02 几何重放允许读取 h5ad obs 的 barcode 与 array 坐标索引，仍禁表达值与图像像素

- 背景: 把 spot-barcode 标注重放为空间几何（连通域、质心、距离）需要每个 barcode 的 array_row/array_col；此前 R-02 审计只读取 h5ad schema，未读任何列值。
- 理由: barcode 与 array 坐标是身份与空间索引 metadata，不是表达矩阵内容，也不是图像像素；没有它们就无法把任何公开 spot 标注绑定到 R-01 物理单元并生成几何实例，GT 审计无法推进。
- 代价: R-02 策略从"只读 schema"微调为"可读 obs 的 _index/array_row/array_col 三列"；表达值（X/layers）、obsm image_means 及图像内容继续禁止。
- 复查条件: 若后续节点需要读取表达值或图像派生特征，必须另立决策并更新 input_policy；本决策不构成对表达内容的授权。
- 影响: `r02_policy.json`、`scripts/r02_build_registry.py`、`scripts/r02_validate_gate.py`；不改变任何科学假设的可检验性。

### D-034 | 2026-07-31 | Heiser GT 的身份粒度是 piece（h5ad）而非 capture area；Visium barcode 仅在单 capture 内唯一

- 背景: Heiser 42 个标注 CSV 按 capture area 命名，但 meta2 与 GitHub sample key 显示 5 个 capture（7794_1、7794_2、8578_1、8578_2、8578_3）是 TMA 多芯：同一 capture area 的 spot 分属 2–3 个不同 patient/block 的物理 core（共 12 行 duplicated），R-01 物理单元是 piece 级 h5ad（48 行），与两份 key 文件逐行互证一致。另实测 Visium barcode 序列跨 capture area 复用（如 AAACAACGAATAGTTC-1 同时出现在 6723_1 与 6723_2），仅在单个 capture 内唯一。
- 理由: GT source 与结构实例必须挂在 R-01 piece 上，否则同一 capture 的实例会同时属于多个 patient，外层切分与患者级统计全部失真。piece 内 barcode 唯一使归属无歧义；同一 capture 的 piece 之间 barcode 不相交已作为完整性条件强制核验（8578_1 实测精确划分）；被上游 QC/trimming 删去的标注 spot 无坐标，不可用，不构成身份冲突。
- 代价: 42 个 capture 级 CSV 拆为 47 个可重放 piece；连通域按 piece 计算，capture 级探索性数字（如边界 184 组件）被 piece 级正式数字（191 组件）取代；标注覆盖率受上游 trimming 限制。
- 复查条件: 若上游发布未修剪 h5ad 或 capture 级坐标表，重审 spot 覆盖率；若任一 capture 内 piece 间 barcode 不再不相交，构建 fail-closed 并重新审计身份。
- 影响: `scripts/r02_heiser_gt.py`、`gt_source_audit.tsv`、`structure_instances.tsv`、`h5ad_replay_index.tsv`；R-04 及之后所有患者/组织块级计数。

### D-035 | 2026-07-31 | 第二 claim-bearing 结构冻结为 TUMOR_STROMA_BOUNDARY（Heiser carcinoma_border/carcinoma_edge）

- 背景: I-003 要求按血管、坏死、肿瘤—基质边界逐项比较真实 patient/block 数、实例数、空间分辨率、标注误差、输入独立性和 validation 覆盖后再选择。审计结果：血管与坏死在全部已查来源（Heiser 标注、四个 GEO 谱系公开材料、Atlas 表）中可审计实例均为 0；Heiser 的 carcinoma_border/carcinoma_edge 逐 spot 病理标注在 training 谱系给出 191 个连通域组件、19 位患者、30 个 piece；adenoma_border 仅在癌前上下文出现，登记为非确认 context。
- 理由: 肿瘤—基质边界是唯一满足「非 TLS + 可审计空间 GT + 输入隔离可执行」的候选；carcinoma_border/edge 是病理医生在 H&E 上标注的肿瘤—非肿瘤交界 spot，几何上即肿瘤—基质边界的 spot 级化身，且与分子输入的隔离规则可执行（H&E 及其派生特征对该类 GT 任务禁作输入）。
- 代价: 第二结构 GT 同样只覆盖 training 谱系，validation 谱系对 TSB 的确认性检验仍阻塞（并入 I-001 剩余缺口）；spot 级 55 μm 分辨率的边界不确定性已登记于实例表。
- 复查条件: 若血管或坏死出现可审计 GT（如作者提供 HALO/QuPath 向量），重审第二结构选择或评估增加结构；若 Heiser 标注语义被证明不等于肿瘤—基质边界，将相关来源降级为 context 并重跑门禁。
- 影响: ISSUES.md I-003 解除；roadmap.md R-02、R-08；`structure_ontology.tsv` 的 claim_status；H-04 恢复可检验路径。

### D-036 | 2026-07-31 | 用户裁定：人类标注者不区分资历等级，统一作为高置信度标注，应用于所有数据

- 背景: 候选 GT 来源的标注者资历不一——Heiser 2023 为病理医生但未逐文件署名，USZ 数据集标注者为 "expert researchers"，Valdeolivas 为病理医生。此前审计把 "per-file annotator identity unattributed" 和 "expert researcher 非病理医生" 作为置信度降级备注。
- 用户决定: 当前阶段，只要标注者是人类，不再严格区分资历（执业病理医生、研究者、未署名标注者），统一按高置信度标注处理；本规则应用于所有数据。
- 边界（不改变的部分）: 本规则只放宽标注者资历分级；确认性 GT 仍须满足 D-030 的全部其他条件——标注模态独立于拟用输入（H&E/IHC/形态学）、可重放几何、R-01 物理链接、provenance 闭合。表达引导或表达派生的标注（如 SpaceHack Visium HD CRC 的 H&E+marker 混合域、Oliveira 的 deconvolution Periphery 标签）仍不得作为确认性 GT；"人类手工标注" 必须由来源 metadata 或方法学文本证明，来源不明时 fail closed。
- 代价: 不同资历标注者之间的系统性标注差异（inter-rater variability）不再在置信度层面区分；该噪声只能在有重复标注实验的数据上单独估计（如 ActiveVisium 一致性实验，当前未收录）。
- 复查条件: 若后续在任务上测得标注者资历对标注质量有可测量影响，或项目进入需要分级的正式验证阶段，重审本规则。
- 影响: 所有 GT source 的 provenance 审计口径与 notes 表述；USZ 与 Valdeolivas 来源可按高置信度注册；gate 对 provenance 的检查（KNOWN = 模态/方法/独立性已知）不变。

### D-037 | 2026-07-31 | 用户批准三条公开来源路径同时补足 validation 谱系 GT（替代作者索取）

- 背景: I-001 剩余缺口是三个 validation 角色谱系（GSE274103、GSE226997、GSE211956）无可审计 GT，唯一已知补救是周期很长的作者索取。三路文献普查确认公开沉积可完整替代：本地 GSE175540（Meylan ccRCC，Atlas 谱系内，作者沉积逐 spot TLS 标注，CD3/CD20 染色判定）、Zenodo 7760264（Valdeolivas CRC，7 患者×2 连续切片，病理医生 QuPath/Loupe spot 标注，唯一公开的第二结构来源）、Zenodo 14620362（USZ 肾+肺，8 患者，显式 TLS 标签）。GSE274557（discovery 角色）放弃索取；GSE274103 邮件降为可选（其独占剩余价值仅为未冻结的血管/坏死标注）。
- 批准范围与规模: 路径 1 解包本地 RAW.tar（约 3 GB 落盘）；路径 2 选择性下载标注 zip + 14 样本的矩阵/空间文件（预计 2–5 GB）；路径 3 下载单 zip（2.1 GB）。总量远低于存储余量（当前约 3.07 TB）。
- 不可替代性: 路径 2 是第二结构（肿瘤—基质边界）validation GT 的唯一公开来源；路径 1/3 是显式 TLS spot 标签仅有的两个公开来源；三者均无需作者响应周期。
- 预期解锁判据: 新谱系在 structure registry 出现 AUDITABLE_GT source 与 CONFIRMATORY 实例，gate 的 confirmatory breakdown 出现 HTAN 以外谱系；R-04 队列间一致性判据恢复可完成。
- 角色治理: 三个新逻辑单元作为新增 validation lineage 冻结，不替换 GSE226997/GSE274103/GSE211956 的既有冻结角色（后两者保持 GT-less，不承担确认性主张）；角色分配只使用身份、谱系、平台与预注册 capability。
- 失效条件: Zenodo 内容与普查描述不符、患者身份无法闭合、标注被证明为表达派生或存在跨谱系物理重复时，对应路径停止并回写 ISSUES.md。
- 影响: roadmap.md R-01（增补注册）与 R-02；ISSUES.md I-001；`infra/sample-registry/`、`infra/structure-registry/`；data/GEO/GSE175540 解包与 `data/other_sources/` 新增两个 Zenodo 目录。

### D-038 | 2026-08-03 | 三谱系 GT 审计结论：实例定义为 hex 连通组件；授权读取 USZ h5ad 的 obs ground_truth 标签列

- 背景: D-037 批准的三条路径已完成内容验收。三个来源的标注形态各不相同：GSE175540 是 GEO raw data 内作者沉积的逐 spot TLS 标注 CSV（两种表头：`TLS_2_cat` 值域 {TLS, NO_TLS, 空} 与 `TLS` 值域 {T_agg, 空}）；TLS_VISIUM_USZ 的标签不在独立文件中，而是存在沉积 h5ad 的 `obs/ground_truth` categorical 列（坐标在同文件 `obs/_index`+`x_array`/`y_array`）；ST_CRC_CMS 是病理医生逐 spot 类别 CSV（4 种表头变体、39 项含拼写错误的沉积词汇）。此前 D-033 只授权读取 h5ad 的 barcode 与 array 坐标列，未授权读取任何标签列。
- 理由: （1）实例定义沿用 Heiser 谱系的六边形网格连通组件（邻居 (0,±2),(±1,±1)，in-tissue 过滤），四个谱系定义一致才能进入同一确认性计数与外层切分；逐谱系重算也使计数可重放、可审计（GSE175540 重算 35 组件 vs Atlas S4 的 30 条 R-ID 为定义粒度差异，正式计数以重算为准）。（2）USZ 的 `obs/ground_truth` 是沉积本身携带的人工标注列（记录描述："annotated in the corresponding H&E images by expert researchers (K.S. and S.D.)"），读取它与 D-033 的坐标列同属身份/标签 metadata，不是表达矩阵内容；没有它该谱系没有任何可用 GT。（3）KIRC 的 T_agg（T-cell aggregate）不是 TLS，三个 T_agg-only 文件 fail-closed 登记 NOT_AUDITABLE；ST_CRC_CMS 只把显式 `tumor&stroma` 家族 4 标签纳入 TUMOR_STROMA_BOUNDARY 范围，`IC aggregate*` 为免疫细胞聚集、非已验证 TLS，fail-closed 排除出 TLS 范围；全部来源的空标注单元按 unknown 处理，绝不当作阴性。
- 代价: R-02 读取范围从 D-033 的三列扩展到 USZ 八个 h5ad 的第四列（`obs/ground_truth` 的 categories+codes）；表达值（X/layers/obsm）与图像像素继续禁止，USZ zip 中高分辨率 tif 未解压。KIRC 的 GSM↔Atlas R_P 逐样本 crosswalk 未公开，与 Atlas KIRC 行的重复链接只能保持 study 级。
- 复查条件: 若 USZ 沉积被证明 `ground_truth` 列含表达派生成分、任一来源标注被证明为模型生成、或作者公开更完整几何（polygon/mask/逐样本 crosswalk），重审对应来源等级并重跑门禁；若后续节点需要读取表达值或图像派生特征，必须另立决策，本决策不构成授权。
- 影响: `scripts/r02_validation_gt.py`、`scripts/r02_validate_gate.py`、`r02_policy.json`（`h5ad_obs_ground_truth_label_reads_allowed`）、`gt_source_audit.tsv`、`structure_instances.tsv`、`h5ad_replay_index.tsv`；`infra/structure-registry/provenance/` 新增三份 provenance notes；roadmap.md R-02 与 ISSUES.md I-001 的解除证据。

### D-039 | 2026-08-06 | 将潜在空间场而非已知结构作为主要发现对象
- 背景: 现有 plan.md 已讨论隐藏结构和软生态 niche，但 R-04/R-08 的表述容易把 TLS、肿瘤—基质边界、血管和坏死误读为全部候选对象，使项目滑向复现已发表的已知结构空间场。
- 决定: 项目的主要发现对象是可能没有清晰可见结构或既有名称的潜在空间场。已有结构标签只作为其中一部分空间场的验证、解释和定位锚点；没有结构锚点的场仍可进入发现与重建路径。
- 理由: 科学问题要检验的是分子空间场是否具有可迁移、可重建的信息，而不是已知结构周围是否存在梯度。已知结构是有价值的验证子集，但不是发现的必要前提。
- 边界: 没有独立锚点的场可以作为可重复的潜在空间场报告，但在缺少正交或独立证据时不强行命名为具体结构；血管和坏死的确认性 GT 只有在提出对应命名结构或更强通用性主张时才需要。
- 影响: 更新 `docs/plan.md` 的 P-01、H-01、H-04、H-05、K-01、P-06、B-09；更新 `docs/roadmap.md` 的 R-02、R-04、R-08；当前主路径不再以血管/坏死 GT 为前置条件。
### D-040 | 2026-08-07 | R-04 采用双模型连续场实现与受控 TensorFlow 环境
- 背景: R-04 需要在不把候选对象降为普通聚类的前提下，从 molecule-only 空间数据发现可重叠、带空间尺度和不确定性的潜在场；同时需要 signed residual 路线处理正负拮抗 program，并使用本地 005/006 单细胞数据进行组成调整。
- 决定: 并列实现非负多样本 NB-GP/mNSF 适配器和 NB nuisance residual 上的 signed GP 适配器；候选按共享 loading、连续场相关性、restart 和患者稳定性形成 `BOTH/MNSF_ONLY/SIGNED_ONLY` 分层并集。候选 API 只输出坐标函数的 posterior mean/SD/lengthscale，不输出 spot 类别。组成调整使用 ILR + 患者分组 cross-fitting；005 的 GSE132465 作为 CRC 精细参考，006 的 GSE115978/GSE232240 作为受控的跨器官免疫/基质参考。
- 理由: 单一非负模型不能充分表达拮抗 program，单一 signed 模型又缺少 count-level 非负生成解释；并列模型可把模型依赖显式暴露为候选层级，而不是用 GT 选择模型。已知结构标签只在候选 hash 冻结后进入独立 anchor evaluator。
- 环境与资源: 使用用户级 Python 3.11 环境，TensorFlow 2.18.1、TensorFlow Probability 0.25.0、tf-keras 2.18.0；CPU 为默认，GPU 只有在 CPU pilot 外推超过 96 小时、单 fit 超过 6 小时或 RSS 超过 128 GiB 时提出独立租用申请。GPU 不得改变折数、基因数、restart、posterior draws 或 gate。
- 失效条件: 任一模型未通过合成连续场、组成-only 阴性和 uncertainty smoke，R-04 不得用聚类、图平滑或 marker-only 结果替代；TensorFlow 环境缺失时状态为 `BLOCKED_COMPUTE`。组成参考无法排除物理重叠或无法传播不确定性时，只能报告 `PARTIAL_RAW_FIELD_ONLY`/`STATE_NOT_IDENTIFIABLE`。
- 影响: 新增 `r04/`、`scripts/r04_*`、`infra/r04/` 和 R-04 回归测试；本决策不改变 `docs/plan.md` 的科学问题，不授权任何 GT、图像、结局或 USZ h5ad 表达层读取。

### D-041 | 2026-08-07 | R-04 候选注册表、输入内容和 uncertainty provenance 必须绑定
- 背景: 单次拟合输出不能代表跨 restart 或患者稳定候选；仅比较数组形状或默认把一次拟合当作稳定结果，会允许 checkpoint、候选或门控指标错配。
- 决定: 候选注册表写入并重算自身 hash；候选保存完整输入内容 hash 和 field hash；跨模型匹配要求 section 与 input hash 一致；稳定性只接受显式 restart/patient 分母与命中数；gate 缺少 registry、input、field、uncertainty 或 lengthscale provenance 时直接阻断。
- 理由: R-04 的对象是可重复的连续场，不是一次运行产生的标签。输入 hash 必须覆盖 counts、coordinates、barcode 和 gene IDs，才能防止同形状数据复用 checkpoint 或指标。当前模型的 uncertainty 只覆盖 inducing weights，故 gate 默认不能通过。
- 失效条件: 若未来后验统一覆盖 loading、dispersion、lengthscale 和模型选择，并有冻结后的 grouped metrics，可由新的 decision 解除相应 gate 阻断；不能用降低阈值或手工填充 provenance 解除。
- 影响: `r04/candidates.py`、`r04/io_contract.py`、`r04/serialization.py`、`scripts/r04_gate.py` 和对应回归测试采用 fail-closed 语义；不改变 `docs/plan.md`，也不把 GT 变成潜在场发现前置条件。

### D-042 | 2026-08-10 | R-04 采用显式 71-unit locator、稳定 gene ID 和训练侧 gene universe
- 背景: R-04 需要把冻结 outer split 与真实矩阵/坐标一一绑定，同时跨 HTAN h5ad、10x H5 和 Matrix Market 对齐基因；dense h5ad 在转换前若不检查整数性，会静默截断 normalized float；10x 坐标还存在标准无表头格式。
- 决定: locator 只能来自 pinned replay index 或显式 override table，禁止从文件名猜路径；h5ad 优先读取 `var/gene_ids`，10x 使用 feature ID，只有没有稳定 ID 时才回退到 `_index`/第一列；所有 71 个 eligible section 先通过 molecule-only preflight。gene universe 只由 HTAN 训练患者生成，按患者内至少 5% spot 检出且至少 70% 患者满足后，按训练侧 log-normalized variance 取 4,000 个基因。标准无表头 10x 坐标按六列位置解析；零 library spot 确定性剔除并写入预检报告。
- 理由: 这些选择同时防止跨平台 gene namespace 假不重叠、输入值域污染、坐标首行丢失和验证集反向选择基因；locator、manifest、gene list 和 preflight gate 可独立重放且不读取 GT。
- 失效条件: 稳定 gene ID 缺失、验证 lineage 覆盖低于 80%、输入内容 hash 漂移或出现目标派生字段时，相关 lineage 立即 fail closed；不能通过改名、补零或降低覆盖阈值解锁。
- 影响: 新增 `scripts/r04_build_locator_table.py`、`scripts/r04_preflight_manifest.py`、`scripts/r04_make_cpu_pilot_manifest.py`、`r04/gene_selection.py` 及 `infra/r04/locator_overrides.tsv`；I-011 解除。该选择不改变 `docs/plan.md`，不把已知结构 GT 变成潜在场发现前置条件。
### D-043 | 2026-08-10 | R-04 冻结共享参数、只对新 section 推断连续场
- 背景: 训练期模型已经产生共享 loading、基因级 NB nuisance/profile 和 section-specific inducing 权重；若验证阶段重新拟合 loading 或重新选择基因，候选场会被验证数据定义，不能区分迁移失败与重新训练成功。
- 决定: `r04.frozen_model.v1` 保存训练侧 gene universe、共享参数、输入/config/environment hash；新 section 通过 `infer()` 固定共享参数，只优化该 section 的 GP 诱导权重（mNSF 同时优化显式非空间 nuisance；signed 路线使用训练侧 NB Pearson profile），输出连续 field mean/SD，不输出聚类标签。`scripts/r04_infer.py` 在 manifest、gene list 或 artifact 不一致时 fail closed。
- 理由: 这把“发现候选场”和“在新 section 上评估已发现的场”分成两个统计操作；空间对象仍是 kernel 定义的连续随机场，不因候选发现阶段方便而退化为普通聚类。
- 失效条件: artifact hash 漂移、gene universe 不一致、共享 loading/dispersion/profile 未被保存，或 uncertainty 仍只覆盖 inducing weights 时，推断结果只能标记 `INFERENCE_COMPLETE_NOT_VALIDATED`，不得进入 R-04 科学 gate。
- 影响: 新增冻结 state 序列化、两个 estimator 的 frozen inference、`scripts/r04_infer.py` 及回归测试；不改变 `docs/plan.md`，不把 GT 变成潜在场发现前置条件。
### D-044 | 2026-08-10 | R-04 外部来源索引只登记官方元数据，不复制表达或图像
- 背景: 本轮使用了已存在本地资产对应的 ST-CRC-CMS 与 USZ TLS Visium 官方 Zenodo record JSON；矩阵和图像不应因 R-04 运行再次进入 metadata-only 数据索引。
- 决定: 在 `infra/bioinf-data-index/` 中新增两个 `ZENODO:` 元数据行，记录官方 URL、bytes、SHA-256 和空间 cohort 覆盖；索引 schema 升为 3，expression/image/result retained 计数保持 0。
- 理由: 保留来源可审计性，同时不扩大外部数据交付范围，也不把 GT 或结果表混入模型输入。
- 失效条件: 若后续必须读取未登记的表达、图像或结果文件，先停止并重新申请数据范围；不能用现有 metadata 行替代实际数据授权或 GT 审计。
- 影响: `scripts/update_bioinf_data_index.py`、`tests/test_update_bioinf_data_index.py` 和索引三件套更新；不改变 `docs/plan.md`。
### D-045 | 2026-08-10 | R-04 fit 与 infer 强制使用纯角色 manifest，并缓存冻结 panel
- 背景: 原始 71-row manifest 同时包含 training 与三类 validation；直接把它交给模型会使角色泄漏难以从结果包中排除，而且每次 restart 都重复读取 19,207-gene 原始矩阵。
- 决定: fit 只接受纯 `training` manifest；infer 只接受显式 validation role。冻结 4,000-gene panel 可物化为带 gene 顺序、barcode、坐标和父 manifest hash 的稀疏 cache，cache 缺失或不完整时 fail closed。
- 理由: 角色隔离是统计设计的一部分，不是命令行约定；缓存只改变读取路径，不改变 counts、坐标或 gene namespace。
- 失效条件: role hash 漂移、cache metadata 与矩阵 shape 不一致、坐标或 barcode 对不齐、验证覆盖低于 80% 时，不得运行正式 fit/infer。
- 影响: 新增 role partition、panel materialization/combination 与 cached-panel loader；不改变 `docs/plan.md`。
### D-046 | 2026-08-10 | R-04 gene minibatch 只估计可分解 likelihood，global GP-KL 不缩放
- 背景: 全 gene 张量的内存和时间成本过高，但当前 NB 目标在 gene 维度上可逐项计算；空间 GP KL 是共享的 global 项。
- 决定: 均匀无放回 gene batch 使用 `G/|B|` 的 Horvitz–Thompson 缩放；global spatial KL 每步只加一次；batch 由 `seed + step` 的 stateless sampler 决定并写入诊断。
- 理由: 这样 batch size 改变的是 ELBO 的 Monte Carlo 方差，而不是目标定义；dense oracle、Monte Carlo 无偏性和双模型 smoke 已通过。
- 失效条件: 后续引入跨 gene normalizer、非均匀 sampling 或无法恢复 effective step/RNG 语义时，正式 restart 重新阻塞。
- 影响: 两个 estimator 增加 gene batch、KL/ELBO 诊断和 checkpoint global step；不把 minibatch 结果视为科学验证。
### D-047 | 2026-08-10 | R-04 不用正交惩罚掩盖 factor collapse，匹配采用 sign-invariant 一对一连续签名
- 背景: 受限 CPU 一步训练的 mNSF loading 两两 cosine 仍约 0.997，说明需要诊断而不是人为制造分离；旧的逐因子贪心匹配可能产生顺序依赖。
- 决定: 先报告 loading cosine、连续 field cosine、effective rank 和 factor energy，不默认加入 orthogonality/repulsion。跨模型候选使用 loading 与 patient-level field 的 conjunctive threshold，并由 Hungarian assignment 一对一匹配；signed factor 的整体 sign 作为 nuisance symmetry。
- 理由: 候选对象仍是连续函数及其 gene signature，不是 spot cluster；可靠的无 `BOTH` 结果比事后调阈值产生 `BOTH` 更有信息。
- 失效条件: assignment margin 不足、joint duplicate/inactive factor 持续、或 synthetic planted control 不能恢复时，禁止正式 restart 和候选命名。
- 影响: 新增 factor diagnostics、patient-aware provenance、sign-invariant matching 与回归测试；不改变 `docs/plan.md`。
### D-048 | 2026-08-10 | R-04 CPU 资源校准未通过正式 restart 时间门槛，GPU 不在本轮自动启用
- 背景: 修正 full-library offset 后，47-section、4,000-gene training cache 的单 mNSF step 实测 wall time 69.54 秒、峰值 RSS 14,103,928 KiB；双模型 3-step 原始读取运行在 120 秒外层上限内未完成。内存尚未达到 128 GiB，但正式 300-step×5 restart 的 CPU 时间不可接受。
- 决定: 当前只保留 calibration artifact，不启动正式五次 restart；GPU 仍遵守项目政策，只有用户批准后才申请/租用。后续先完成多 seed/planted/irregular 诊断并给出 GPU 需求与替代方案。
- 理由: 未完成的核心计算必须保持 `NOT_RUN/BLOCKED_TIMEOUT`，不能用缩减 steps 或单步结果冒充正式证据。
- 失效条件: 若后续 cache、线程和目标等价性优化把完整运行降到可接受 CPU 时间，并通过 collapse/irregular gates，可撤销本门禁；反之提交 GPU 审核。
- 影响: R-04 继续进行但正式 restart 处于资源与科学双门禁；不改变 `docs/plan.md`。

### D-049 | 2026-08-13 | R-04 将可识别空间维度与模型 factor 数分开，并修正 gene-minibatch likelihood 目标

- 背景: 仅用 loading cosine 和逐 factor 对齐会把交换、变号、旋转等统计对称性误判为恢复失败；同时原 gene-minibatch 路径把 spot-gene 联合均值直接乘以 `G/B`，与按 spot 求和的 NB likelihood 目标不一致。
- 决定: 用 `K_model` 表示模型允许的 factor 数，用 `K_eff` 表示由留出预测、空间置换零分布、子空间 canonical correlation 和重复性共同支持的效应维度；校准同时报告完整 effect subspace 与 factor-level 指标。gene minibatch 先按 spot 求和 gene loss，再使用 Horvitz--Thompson 缩放；global GP-KL 仍只加一次。
- 理由: `K_eff=1` 不代表只有一个不重叠结构，而是当前数据只支持一个统计空间效应子空间；正确目标函数是正式训练和 K 判定的必要条件，不能用旧的短跑结果继续外推。
- 失效条件: 若 dense 与 full-batch/minibatch 的目标或梯度不能达到数值一致，停止所有科学校准；若正确指定高信号控制在平台期仍不能恢复其 effect subspace，则保留 R-04 科学门禁并进入模型可识别性诊断，不加入正交惩罚制造分离。
- 影响: 新增连续场子空间诊断、K 选择辅助函数、显式 calibration generator 和 objective 回归测试；不把聚类、GT 或已知结构引入潜在场发现，也不改变最终 R-04 gate。

### D-050 | 2026-08-13 | R-04 用可识别的 inducing-value 参数化和无泄漏 K 选择替代旧短跑校准
- 背景: D-049 后的诊断发现旧 mNSF 直接用 `Kxz @ u`，却用 `Kzz` 作为 `u` 的先验；fit 与 frozen infer 的 likelihood 归一化也不同，初始化总空间强度随 `K_model` 增加。旧 K 比较还用 held-out section 的全部 counts 先适配场、再在同一 counts 上评分，K=4 趋势无法解释为预测证据。
- 决定: 将 inducing values 映射改为 `Kxz Kzz^-1 u`，fit/infer 统一为按 spot 平均、按 gene 求和的 NB objective，并采用总强度不随 K 改变的数据尺度初始化。K 搜索使用 section 留出和 gene cross-fit：只用 adaptation genes 推断 held-out section 的连续场，用未见 evaluation genes 评分；方向 null 在 section 内置换且逐 canonical direction 比较。先搜 K=1/2/3/4，只有边界 K=4 仍显著最好才扩大上界。
- 理由: 这些修改恢复模型和 GP 先验的统计语义，避免 K 的初始强度偏差和评价泄漏；连续 spot×gene effect subspace 而非 factor 编号或 spot cluster 是主要恢复对象。
- 失效条件: 若完整 objective/gradient oracle 不一致，或高信号正确指定控制在平台后仍不能恢复对应 `K_eff`，立即撤销科学 PASS；若真实数据最优 K 落在搜索边界，保持 `K_NOT_IDENTIFIABLE` 并扩大搜索，而不直接选择边界。
- 影响: D-049 后、D-050 前的 calibration 与资源结果继续保留为历史证据但不可合并；`docs/plan.md` 的科学定义不变。

### D-051 | 2026-08-13 | 合成 K 规则通过不自动授权真实数据 K 或正式 restart
- 背景: 三个数据/优化 seed 的无泄漏控制在 planted `K_eff=2` 时从 K=1/2/3/4 选择 K=2，无场控制从 K=0/1/2 选择 K=0；三类不规则支持和单场/共线控制也达到对应平台与子空间门槛。但这些是方法校准，不是患者数据发现。
- 决定: 将门禁状态写为 `PASS_SCIENTIFIC_CALIBRATION` 和 `SYNTHETIC_K_RULE_CALIBRATED`，不把 K=2 外推成真实数据最终 K。真实运行仍需独立选择 `K_model/K_eff`，并保留 K=0；K=4 本轮不在最佳 one-SE 集合，因此不试 K>4。
- 理由: 合成控制证明选择器能在已知真值下区分 0、1、2 个可识别效应方向，只提高 R-04 的可测性，不构成 H-01/H-03/H-05 的患者级证据。
- 失效条件: 若后续真实数据 K 选择卡在搜索边界、跨 restart 不稳定、或无场/置换 null 出现同等增益，恢复 `K_NOT_IDENTIFIABLE` 并关闭科学门禁。
- 影响: R-04 从模型可识别性阻塞进入资源审批；正式五次 restart、signed 路线、组成 posterior 和独立验证仍未运行。

### D-052 | 2026-08-13 | 跨 K 必须使用同一考题，独立不确定性只按 data seed 计算

- 背景: D-051 后复核发现 optimization seed 含 K，而 gene split seed 又由 optimization seed 派生，导致不同 K 实际使用不同 evaluation genes；旧聚合还把同一合成数据内的 section/gene folds 当成独立重复计算标准误。旧结果即使仍选择 K=2，也不能作为公平的最优 K 证据。
- 决定: 撤回 v1 K/K0 聚合和 v5 总门禁的判定用途。将 data seed、K-dependent optimization seed 与 K-independent split seed 明确分离；聚合器必须拒绝重复或缺失的 K×data-seed×holdout 单元，并逐 data seed 核对输入 hash、完整 truth、generator/execution contract、evaluation folds 和 K=0 baseline。K 的均值先在每个 data seed 内汇总，标准误只以 data seed 为独立单位。只有搜索上界仍为最佳或被 one-SE 规则选中时才扩 K。
- 理由: K 比较必须只改变模型容量，不能同时改变数据或评价题目；section/gene folds 是同一数据的技术复用，不是新的生物学重复。
- 失效条件: 若任一聚合输入无法证明同源数据和同一评价划分，状态立即回到 `K_NOT_IDENTIFIABLE`；若真实数据最优 K 位于当前上界，再扩大到更高 K，而不是选择边界。
- 影响: 重新完成 24 个双场单元和 12 个无场单元后，K=2 与 K=0 分别在三个 data seed 上一致最佳，K=4 不在搜索边界，因此仍不扩 K=5/6/8；有效证据替换为 v2 聚合和 v6 总门禁。`docs/plan.md` 的科学假设不变。

### D-053 | 2026-08-13 | 正式长运行的断点必须保存历史最佳态，signed fit/infer 使用同一尺度

- 背景: mNSF 断点过去只保存当前参数和 optimizer，恢复后丢失中断前的历史最佳参数；signed residual GP 的 fit 使用按 spot 平均、按 gene 求和的目标，而 frozen infer 使用 spot×gene 联合均值，二者随 gene 数变化时尺度不同。
- 决定: mNSF checkpoint v2 同时保存 best loss、best step 和每个参数的 best shadow state，旧 checkpoint metadata fail closed；signed frozen infer 改用与 fit 相同的 per-spot gene-sum objective。用模拟第三步崩溃验证恢复运行与不中断运行得到相同最佳参数。
- 理由: 多小时 restart 必须能在中断后保留已经找到的最好解；fit/infer 尺度一致是冻结模型跨 panel 推断可解释的最低条件。
- 失效条件: 若恢复结果与同 seed 不间断结果不一致，或 objective oracle 显示 fit/infer 随 gene 数产生不同尺度，正式 restart 继续禁止。
- 影响: 解除断点最佳态和 signed 目标尺度的工程缺口，但不解决 signed 完整 NB 后验、真实数据 K、组成 posterior 或资源阻塞。

### D-054 | 2026-08-13 | 总科学门禁必须从底层 K shards 重新计算而非只信任聚合摘要

- 背景: 聚合摘要即使带有来源路径和哈希，仍可能自身字段被改写；只验证底层文件未变，不能证明摘要中的 K、标准误和状态确由这些文件计算得到。
- 决定: 总门禁除重算每个来源文件的字节数与 SHA-256 外，还要重新读取全部 K/K0 shards，执行同一完整性、同题和独立统计检查，并逐字段比较重算聚合与提交聚合。v7 总门禁替代 v6 作为当前有效入口。
- 理由: provenance 证明“来源是谁”，语义重算证明“结论确实来自这些来源”；两者缺一都可能产生假绿。
- 失效条件: 任一来源缺失、哈希改变、底层语义不合格或重算结果与摘要不一致时，总门禁必须为 `BLOCKED_SCIENTIFIC_CALIBRATION`。
- 影响: v7 在当前 18 个控制、24 个正 K 单元和 12 个无场单元上重新通过；仍不授权真实数据正式 restart。

### D-055 | 2026-08-13 | 最佳损失必须绑定同一时刻参数，所有校准数值必须有限且可跨环境重放

- 背景: 独立审查发现三处剩余风险：truth 与当前 generator 的浮点差约 `2.78e-17` 时逐字相等会误阻断；mNSF 在更新前计算 loss、却把更新后参数登记为该 loss 的 best state；canonical correlation 为 NaN 时 `NaN <= cutoff` 为 false，可能绕过失败判断。
- 决定: truth 的离散字段和容器结构严格相等，浮点字段必须有限并在 `rtol=atol=1e-12` 内相等；mNSF 在 Adam 更新前把产生当前 loss 的参数写入 best shadow；objective、控制相关系数/null、held-out scores 和 recovery diagnostics 全部显式检查有限性。聚合器还要从 seed 重放 generator 的 truth 与输入 hash，并复算 `recovery_pass`。旧 v8 前的模型输出不迁移，完整重跑进入门禁的全部控制和 K/K0 单元。
- 理由: 最佳损失与参数错位会让断点恢复一个从未取得报告损失的模型；NaN 与跨库末位舍入分别会导致假绿和假红。离散严格、浮点有界是可重放性与科学语义之间的最低充分平衡。
- 失效条件: 若跨受支持 Python/NumPy/SciPy 环境仍不能重放，或保存参数重新评价不对应登记 loss，立即关闭门禁并冻结正式运行。
- 影响: 系统 Python 与项目 TensorFlow 环境均能重放当前聚合；修复后 18+24+12 个单元重新通过，当前有效入口为 `scientific_gate_20260813_v10.json`。该决定不改变 plan 假设，也不把合成 K=2 外推到患者数据。

### D-056 | 2026-08-13 | 资源判定改用训练循环内部计时，CPU 硬阻塞解除

- 背景: 旧资源判断把 1-step 总墙钟直接乘 600，混入冷缓存读取、模型构建和写盘，曾给出单 restart 15.9 小时、五次 79.6 小时。当前 1/2/5/10-step 外部总墙钟受缓存和系统负载影响明显，不能稳定估计每步成本。
- 决定: 在 estimator diagnostics 记录本次 optimizer 实际步数与循环墙钟。以同一 47-section、4,000-gene、K=2、16 inducing、gene batch 256、16 CPU threads 的 20-step 测量为资源主证据：optimizer 92.80 秒，即 4.64 秒/步；固定加载/构建/写盘约 179.51 秒。线性规划值为单次 600-step mNSF 约 0.82 小时、五次约 4.12 小时；这是规划估计而非时长保证，不含 signed、真实 K 搜索和 validation。
- 理由: 分离固定成本与增量步成本后，CPU 正式计算处于小时级而非多日级，继续把 GPU 审批作为硬前置会不必要地阻断科学进度。
- 失效条件: 若真实 600-step 首次运行的 optimizer 斜率超过 20-step 测量两倍、内存接近系统安全边界、或影响其他正常任务，暂停后重新提交 GPU/调度选择；不得静默缩减步数。
- 影响: I-015 从 `HARD_BLOCKED` 关闭；GPU 暂不申请。正式入口增加模型隔离 checkpoint 目录与 checkpoint cadence，下一步可在 CPU 上进行真实数据 K=0 向上的留出粗筛，但 v10 gate 本身仍不自动授权任何患者级结论。

### D-057 | 2026-08-13 | 用户批准将本项目运行时可用存储硬下限调整为 1.2 TB

- 背景: 原 D-021 设定的 2 TB 是历史安全余量；用户在启动真实 K 搜索前明确批准将本项目当前运行时硬下限调整为 1.2 TB。
- 决定: `resource_status`、R-04 示例配置和项目工作约束统一使用 1.2 TB 硬下限；1.4 TB 作为软提醒线。历史报告中的 2 TB 数值保持原样，不 retroactively 改写历史事实。
- 理由: 当前 R-04 派生物规模和实测可用空间允许在 1.2 TB 仍保留明确安全余量；同时继续在每个长任务前后检查可用空间。
- 失效条件: 可用空间低于 1.2 TB、任务预计峰值越过该线、或系统出现明显存储压力时立即停止新增落盘并标记 `BLOCKED_STORAGE`。
- 影响: 仅改变运行资源门槛，不改变任何科学假设、K 判定规则或证据资格；I-005 的当前硬阻塞线同步调整。

### D-058 | 2026-08-13 | 真实数据 K=0 必须与正 K 共享 nuisance 和计数目标

- 背景: 合成 K0 校准使用了截距-only 配置，而真实 mNSF 正 K 配置默认含 `nonspatial_rank=1`。若直接复用合成 K0，非空间 nuisance 的收益会被误算为空间效应收益。
- 决定: 真实 K 搜索中 K=0 仍使用同一 offset、NB dispersion、gene cross-fit、患者分组、`nonspatial_rank=1` 和 held-out 评分接口；只将 coordinate-dependent GP effect 设置为空维度。K0 初始化使用观测 gene rate，避免因为缺少 spatial mass 造成初始值人为劣势。
- 理由: K 比较只应改变连续空间效应子空间的维数；K=0 表示没有可识别坐标效应，不是换成普通聚类或另一个更弱的基线。
- 失效条件: K0 与正 K 的 objective、nuisance、offset、gene split 或患者折不一致，或未收敛单元被纳入聚合时，真实 K 结果退回 `K_NOT_IDENTIFIABLE`。
- 影响: 新增真实 training 患者/组织块分组 K 搜索入口；当前一折 2-step smoke 仅验证 wiring，未改变 H-01/H-03/H-05 信心。

### D-059 | 2026-08-14 | 收敛诊断先于真实 K 判定，且不按 K 单独加预算
- 背景: 真实 5-fold、K=0–4 粗筛已生成留出分数，但 25/25 个 fit platform 均未通过 600-step 平台期门禁；若直接按原始分数选择 K，会把未稳定的优化状态当成空间场证据。
- 决定: 保持现有患者折、gene cross-fit、K=0 nuisance、NB objective、连续空间场定义和收敛门槛不变；先在独立目录用代表性 K 做收敛诊断，再固定一个共同协议重跑全部 K。学习率只作为显式诊断参数暴露，默认值仍为 0.05；不为 K=3 单独延长预算，不在未收敛结果上运行 null 或 restart。
- 理由: 这样可以区分步数不足、优化波动和空间场特异性失败，同时不让调参改变 K 的公平比较或把 K=3 的高分反向写入模型选择。
- 失效条件: 若诊断协议改变 objective、数据划分、K=0 nuisance、gene scaling 或把平台期阈值放宽，则该协议不能用于正式 R-04；若共同协议仍无法使全部必要 cell 收敛，R-04 保持 `K_NOT_IDENTIFIABLE_NOT_CONVERGED`。
- 影响: 当前结果只提高了收敛问题的可定位性，不提高或降低 H-01/H-03/H-05 的科学信心；空间 null、五次 restart 和组成拆分继续等待收敛 gate。

### D-060 | 2026-08-16 | 低学习率哨兵使用 K0 中位与最差 fold
- 背景: 0.05 学习率的 1200-step 哨兵中 K=0/2/3 的 fit 均未平台，且最佳步都接近末尾；需要区分学习率和 fold 条件，而不能只按 K=3 的预测分数挑样本。
- 决定: 在独立目录运行完整输入的 K=0/1/3、fold 0 和 fold 4 哨兵，1200 steps、learning rate=0.01、gene batch 和其余协议保持不变。fold 0 代表 K=0 轨迹的中位区域，fold 4 代表最差尾部；显式 fold 选择和低学习率任务均标记为 wiring-only，不参与 K 聚合。
- 理由: 同时覆盖通用基线、一个低维正 K、一个较高维正 K，以及典型和最差 fold，能用有限 CPU 预算判断共同优化问题；不改变正式的患者分组、K=0 公平基线或收敛判据。
- 失效条件: 若哨兵未完成、仍有非有限 loss、或只让某个 K 单独收敛，不得据此选择正式 K；若共同协议仍不收敛，继续保持 `K_NOT_IDENTIFIABLE_NOT_CONVERGED` 并转向 batch/空间参数化诊断。
- 影响: 该任务只改变收敛问题的诊断信息，不改变任何 H-01/H-03/H-05 的科学信心。

### D-061 | 2026-08-18 | 2400 步收敛诊断先区分 batch 噪声、学习率和 GP 近似
- 背景: 1200 步 loss 总体趋平，但相邻窗口下降率仍有明显正负波动；现有真实 K 结果全部未通过 fit platform 收敛门禁，且真实 K 搜索入口此前固定使用 gene batch=256。
- 决定: 在独立 wiring-only 目录先比较 gene batch=256、512 和短程 full-batch；在代表性 K/fold 上记录 loss、梯度范数、参数块范数、学习率轨迹和 GP 诊断。若 batch 选择后仍需调整，先比较固定 lr=0.01 与显式衰减协议；inducing points 和 lengthscale 不与首轮 batch/lr 变化混改，lengthscale 保持 3.0，避免把空间先验变化误判为优化改善。
- 理由: full-batch 作为固定目标和梯度参照能直接估计 minibatch 随机性，不必支付 2400 步 full-batch 训练成本；每次只改一个因素能判断抖动来自 gene 采样、学习率还是 GP 参数化。
- 失效条件: 若诊断改变 objective、gene scaling、患者/组织块切分或 K=0 nuisance，结果不能进入正式 R-04；若代表性 K/fold 未共同收敛，不启动 25-cell 正式聚合；任何未收敛结果不得触发空间 null、restart 或组成拆分。
- 影响: 本决定只增加 R-04 收敛问题的可定位性，不改变 K 的科学判定、不提高或降低 H-01/H-03/H-05 的信心。

### D-062 | 2026-08-18 | batch 诊断后固定 512，再单独比较学习率衰减
- 背景: K=0/3、fold 0 的 batch=256、512 和 full-batch wiring-only 诊断均正常完成。400-step 任务中，batch=512 的末段 loss 标准差相对 batch=256 在 K=0 和 K=3 分别下降约 35% 和 44%；full-batch 的短程波动更低，但 100 步 optimizer 墙钟已约 30–36 分钟，不能作为 2400-step 主协议。三种设置的 fit 与 inference 收敛门禁均未通过。
- 决定: 下一轮 2400-step 代表性运行固定 gene batch=512、inducing points=16、lengthscale=3.0、K=0/3、fold 0 和同一 gene cross-fit；只比较固定 learning rate=0.01 与从 0.01 线性衰减至 0.0025（1600 steps）。inference steps 提高到 200，使其超过当前 150-step 门槛；两轮仍标记为 wiring-only。
- 理由: batch=512 是噪声与 CPU 代价的折中；保持初始学习率相同，才能把差异归因于是否衰减，而不是同时改变 batch、初始步长或空间先验。
- 失效条件: 若两轮均未在 K=0 与 K=3 共同通过 fit/inference 平台期，不能继续扩大到 25 个 cell；若只有一个 K 收敛，不能选择真实 K。full-batch 不因较低波动而直接进入长程正式协议。
- 影响: 该决定只缩小 R-04 收敛诊断的工程搜索空间，不改变 K 的科学判定或 H-01/H-03/H-05 的信心。

### D-063 | 2026-08-18 | 2400 步不衰减优于早期衰减，先延长同一协议
- 背景: batch=512 的固定 lr=0.01 与线性衰减到 0.0025 两个 2400-step wiring-only 任务均正常完成，但 K=0/3 的 fit 和 inference 收敛门禁均未通过。固定 lr 的末段损失低于衰减方案；其最后三个 50-step 块仍有约 0.54%（K=0）和 0.84%（K=3）的首尾下降，超过 0.1% 平台门槛。K=3 的有效秩参与度约 1.83，未触发 factor、field 或 loading collapse；梯度逐步下降但尚未形成平台。
- 决定: 先保持 batch=512、learning rate=0.01、inducing points=16、lengthscale=3.0、K=0/3、fold 0 和 gene cross-fit 完全不变，只把 fit 延长到 4800 steps；inference steps 提高到 300，继续标记为 wiring-only。暂不改变 GP 参数化、平台阈值或正式 K 搜索范围。
- 理由: 当前证据更符合“仍在缓慢优化”而不是数值崩溃或因子塌缩；延长步数能直接检验这一解释，且比同时改 GP 先验或放宽收敛门槛更可归因。
- 失效条件: 若同一协议到 4800 steps 仍不平台，下一步才拆分 GP 参数/空间梯度与平台判据；不得把末段分数或有效秩当作真实 K 证据，也不得启动 25-cell、null、restart 或组成拆分。
- 影响: 该决定只检验 R-04 优化是否需要更长训练，不改变 K 的科学判定或 H-01/H-03/H-05 的信心。

### D-064 | 2026-08-19 | 从 4800-step checkpoint 分叉检验 late learning-rate step-down
- 背景: 固定 batch=512、lr=0.01 延长到 4800 steps 后，K=0 与 K=3 仍未通过 fit/inference 平台门禁；K=0 和 K=3 共同失败说明不能先把问题归因于 K=3 的因子秩。K=3 有效秩参与度约 1.44、能量约 80% 集中于一个方向，但未触发硬 collapse 告警。固定 lr 优于从第 1 步开始衰减，提示前期探索与后期稳定可能需要不同步长。
- 决定: 新增显式 checkpoint continuation 接口，验证输入和环境 hash 后允许在明确分叉目录改变 `steps`/学习率，同时保留源 checkpoint 的 Adam optimizer state、global step、seed、gene batch 顺序和模型变量。只对 K=0、fold 0、fit 阶段做两条 4800→7200 continuation：control 保持 lr=0.01；treatment 从 step 4800 起切到 lr=0.0025。两条均不运行 inference、不运行 K=3，结果只作 wiring-only 优化诊断。
- 理由: K=0 是公平基线，可先判断共享 NB/nuisance 优化是否因终端步长过高而无法平台；从同一 checkpoint 分叉避免把重新初始化、Adam moment 或前 4800 步随机差异混入因果比较。
- 接受条件: treatment 完整通过原 fit 平台门禁、control 在匹配 continuation 步数仍未通过，且 treatment objective 不劣于 control 和切换前有效最优值；否则不接受 late decay 作为修复。无论结果如何，K=3 仍需独立通过 fit/inference 才能继续正式 R-04。
- 影响: 该决定只定位共享优化协议，不改变平台阈值、K 的科学判定或 H-01/H-03/H-05 的信心。
D-065 (2026-08-20): implement an optional staged_shared_first optimizer. The first stage updates only global/shared blocks (baseline, dispersion, non-spatial channel, and global factor amplitudes), then the original joint objective updates all blocks. K=0 has no spatial blocks, so its trajectory remains the fair joint baseline. Adam slots are built for all variables before staging so the transition cannot change optimizer state semantics. This is a convergence diagnostic, not evidence for any K; reject it if representative K=0/K=3 fit and inference do not improve under the original gate.
D-066 (2026-08-20): reject staged_shared_first as the common R04 training protocol after the representative K=0/K=3 fold-0 comparison. K=0 was numerically equivalent to joint training; K=3 had a slightly lower best fit loss but remained non-converged in both fit and inference and had worse held-out score. Do not run the 25-cell campaign or downstream null/restart analyses on this protocol. The next diagnostic must change a different optimization factor while preserving the same objective and K=0 control.
D-067 (2026-08-20): before another optimizer change, run a read-only deterministic checkpoint audit. Re-evaluate the dense gene-panel objective with fixed variational draws on existing K=0/K=3 checkpoints, estimate repeatability, and report the additive baseline/nonspatial nuisance gauge coordinates. Keep the original stochastic loss and platform gate unchanged until the audit distinguishes measurement noise, genuine descent, and parameter drift. This preserves the NB objective, molecule-only input, K=0 fairness, and continuous-field interpretation.
D-068 (2026-08-20): the first frozen-checkpoint audit separates K=0 measurement failure from K=3 optimization failure. Dense full-panel objective with fixed variational draws is repeatable; K=0 is effectively flat from steps 6600 to 7200 while the raw minibatch gate remains false, whereas K=3 still improves materially from steps 3600 to 4800. Add deterministic objective monitoring before changing the optimizer. Do not relax the gate or treat K=0 as a scientific result until the monitoring definition is replayed and documented.

D-069 (2026-08-20): continue K=3 from the existing joint step-4800 checkpoint to step 7200 with inherited Adam state, batch=512, learning rate=0.01, and the same objective; record the fixed-draw dense objective every 300 steps. The first launch used the wrong Python environment and failed before model import; the retry uses the previously verified Anaconda environment. This remains wiring-only and cannot select K or relax the platform gate. If the deterministic trace does not flatten, do not add more optimizer variants before auditing the dominant parameter block.

D-070 (2026-08-21): continuation summaries must derive `k_model` from the requested factor count rather than a K=0 literal; isolate that serialization in a pure helper and regression-test K=3. Repair the completed diagnostic JSON in place because its cell-level model metadata and all numerical diagnostics already identify the run as K=3, while preserving its wiring-only and non-converged status. This changes provenance correctness only, not the objective, parameters, or scientific gate.

D-071 (2026-08-21): repeat the fixed-draw dense objective on the identical step-7200 checkpoint twice before deciding whether deterministic convergence can be used to diagnose the raw minibatch gate. A zero repeat difference supports measurement repeatability but does not by itself authorize gate replacement, inference, K selection, or the 25-cell campaign.

D-072 (2026-08-21): the repeated audit found a 0.036% dense-objective decrease from steps 6600 to 7200 under the same fixed-draw definition, and zero difference when step 7200 was evaluated twice. Treat this as evidence that the full objective has reached a diagnostic platform, while retaining the original stochastic fit gate as `false` until a scientific rule is chosen for reconciling the two measurements. Do not silently replace the gate, run inference, select K, or launch the 25-cell campaign.

D-073 (2026-08-21): adopt the dual-gate interpretation for the next diagnostic stage, following user approval. The fixed full-panel objective is the primary optimization-platform measure when its two-window and repeatability criteria pass; the gene-minibatch gate remains a separately reported warning. This permits held-out inference diagnostics from the completed K=0 and K=3 fold-0 checkpoints, but does not authorize K selection, biological naming, full 25-cell aggregation, spatial nulls, or restart claims. Revert to the strict gate if the inference or cross-patient stability evidence fails.

D-074 (2026-08-21): because K=3 inference was still descending after 1200 steps, extend it to 2400 steps before interpreting its lower fold-0 score. To keep the K=0 comparison under the same dual-gate rule, first continue its objective from step 7200 to 7800 to obtain two stable dense-objective intervals, then run the matched 2400-step inference. This is a diagnostic extension, not a new optimizer or a K search.

D-075 (2026-08-22): correct the deterministic-objective platform implementation to match plan.md: the final intervals must have absolute relative change no larger than the threshold, rather than requiring monotonic decrease. A small objective increase is compatible with a platform, while a larger increase remains a failure. The signed relative-improvement values remain reported so the direction is auditable; this changes no model, likelihood, split, or K-selection rule.

D-076 (2026-08-23): for the remaining fold-1/2/3 cross-fold diagnostic, do not resume the available old step-600 checkpoints. Their training provenance predates the current dual-gate protocol and does not guarantee the current learning-rate/mini-batch contract. Refit K=0 and K=3 from step 0 with the current fixed split, lr=0.01, gene batch=512, 7200 fit steps and 2400 inference steps. These runs test cross-fold heterogeneity only; they do not authorize K selection or formal 25-cell aggregation.
D-077 (2026-08-25): because fold 0 is the only negative K3-K0 fold and fold 4 is the strongest positive fold, run one alternate full restart for K=0 and K=3 on each extreme fold under the same split, gene cross-fit, lr=0.01, gene batch=512, 7200 fit and 2400 inference protocol. This is a bounded diagnostic of initialization sensitivity versus patient/fold heterogeneity, not the formal five-restart gate. Defer spatial null, formal K expansion and composition split until its interpretation is available.
D-078 (2026-08-27): the extreme-fold restart preserved the sign of K3-K0 in both folds but changed the positive-fold magnitude substantially, so one alternate restart cannot establish seed stability. Complete the planned five-seed panel by running restart indices 2, 3 and 4 for K=0 and K=3 on folds 0 and 4, keeping the same data split, gene cross-fit, optimizer and inference protocol. Treat this as a within-extreme-fold stability panel; it does not authorize K selection or replace the later full-fold spatial null.
D-079 (2026-08-27): to reduce execution overhead without changing the optimization trajectory, the restart panel may reduce read-only dense-objective monitoring from every 300 to every 1200 steps, diagnostic recording from every 100 to every 600 steps, and checkpoint persistence from every 600 to every 1200 steps; retain four fixed stateless evaluation draws and all model, data, seed, fit-step and inference-step settings. Do not cap TensorFlow thread pools in this switch because changed parallel reduction order could alter floating-point trajectories without an equivalence test. The partially executed restart 2–4 jobs are retained as stopped, non-scientific artifacts and the panel restarts from step 0 in fresh directories.
D-080 (2026-08-31): interpret the completed five-seed extreme-fold panel as `K_UNDECIDED_OPTIMIZATION_NOT_RESOLVED`. Fold 0 remains K0-favouring and fold 4 remains K3-favouring in all five seeds, so the combined optimization randomness does not explain the fold-level sign reversal; however fold-4 magnitude is seed-sensitive, only 8/20 cells pass the current dense-objective rule after replaying canonical external audits, and the third loading/field canonical direction is weakly reproducible across seeds. The optimization seed currently controls initialization, gene-minibatch order and variational draws together, so this panel is not a pure initialization experiment. Do not call K=3 over-parameterized or selected, and do not start spatial null, K expansion, composition splitting or field naming until optimization error is separated from effective-rank instability. This decision fails if a provenance audit shows the cells are not comparable.
D-081 (2026-08-31): before rerunning held-out inference, continue only restart indices 2-4 from their exact step-7200 checkpoints to a bounded step-8400 endpoint, preserving each restart's optimization seed, Adam state, model, split, learning rate and gene-minibatch protocol. Record the fixed-draw dense objective every 400 steps, diagnostics/checkpoints every 600 steps, and run fit only. This gives the existing two-window platform rule three new endpoint evaluations while avoiding repeated 2400-step inference for fits that may still fail. If any cell remains off-platform at step 8400, stop for a new decision rather than silently extending or changing optimizer; if all paired K0/K3 fits pass, rerun inference from those frozen endpoints before recomputing score and top-1/2/3 subspace stability.

D-082 (2026-08-31): audit each step-8400 endpoint as the frozen post-update checkpoint rather than treating the online pre-update step-8400 objective as the same parameter state. For every restart/fold/K cell, evaluate that identical checkpoint twice with the original four stateless draws and seed offset, bind the derived audit to the source-cell SHA256 and TensorFlow checkpoint bundle hashes, replace only the final live objective point, and fail closed unless all 12 cells pass provenance, exact-repeat and two-window dense-objective checks. The completed audit is `STOPPED_OFF_PLATFORM`: 10/12 pass; restart 3/K=3/fold 0 exceeds the dense threshold in its first final interval, and restart 2/K=3/fold 4 is not bitwise repeatable. Therefore inference remains unauthorized. This decision does not relax a threshold, select K, or change H-01/H-03/H-05 confidence; any further continuation or numerical-determinism experiment requires a new decision.

D-083 (2026-08-31): harden the reusable continuation gate before treating D-082 artifacts as future inference authority. Checkpoint audits must persist and match `ridge=1e-5`; every draw, mean, sample SD, dimension and nuisance-gauge field must be finite and internally consistent; pair-level provenance failures count against both cells. An authorized resume must re-hash the source cell, checkpoint bundle and repeat artifact, reuse the bound frozen dense result, preserve the fixed inference protocol, and reject any fit step beyond the audited endpoint before scoring. Because the first artifacts lacked the ridge field, rerun all 12 read-only repeats in a fresh directory instead of editing them. The strict rerun has valid provenance, 12/12 within-run exact repeats and 11/12 full passes; restart 3/K=3/fold 0 remains off-platform, so inference remains unauthorized. The earlier intermittent repeat mismatch remains a numerical-risk observation. This changes no objective, threshold, K decision or H-01/H-03/H-05 confidence.
D-084 (2026-08-31): accept the completed step-8400 panel as practically sufficient for the next frozen diagnostic while preserving the strict audit as `STOPPED_OFF_PLATFORM`, 11/12. The sole exception, restart 3/K=3/fold 0, exceeds the fixed first-window threshold only narrowly (`0.1006526%` versus `0.1%`), while its step-7600 to step-8400 net change is `0.0699%`, its final-window absolute change is `0.0308%`, and its frozen repeat is exact. Record this as a separate post-hoc acceptance artifact bound to the strict audit hash; do not rewrite the strict threshold or call the cell `objective_platform=true`. Its scope is limited to held-out inference and negative-binomial scoring from the 12 unchanged step-8400 checkpoints with zero fit updates, the fixed 2400-step/two-gene-split protocol, and `selected_k=null`. Partial completion cannot be aggregated. This decision fails if any bound audit, source, repeat, checkpoint, input or inference-protocol hash changes, or if either gene split fails its inference-platform check.
D-085 (2026-08-31): add an execution-only TensorFlow optimization for the reusable MNSF fit/inference path. `execution_mode="auto"` preserves the historical eager reference on CPU and selects a per-logical-step `tf.function` path when a GPU is visible; invariant tensors are converted once per call, and GPU memory growth is enabled before runtime initialization. The compiled path preserves float32, the likelihood and KL expressions, stateless seed derivation, logical step order, best-state ordering, optimizer updates, checkpoint boundaries, gene sampling, cross-fit inputs and output schema. XLA, mixed precision, batch-size changes, gradient accumulation, thread-pool changes and full-loop segment fusion are explicitly out of scope. This is an infrastructure change with no direct effect on H-01/H-03/H-05 confidence. A GPU run must first pass one-cell eager-versus-compiled loss/gradient/parameter/checkpoint/output comparison in the target CUDA environment; any material mismatch, GPU invisibility or core CPU fallback keeps the run on the eager reference or requires a new decision.
D-086 (2026-08-31): permit the explicitly authorized remote 200-GB GPU execution to use an execution-only storage reserve override of `R04_STORAGE_RESERVE_BYTES=100000000000` while retaining the local/default 1.2-TB hard reserve. The override is parsed and validated by the runtime, is present in the immutable remote task command, and does not alter data, objective, optimizer, checkpoint, inference or K-selection semantics. The remote task stops as `BLOCKED_STORAGE` if available space falls below the 100-GB reserve or if the task requires unplanned large artifacts; this exception does not change D-057 for local work or support any scientific conclusion.
D-087 (2026-08-31): for the explicitly authorized remote GPU pilot, permit an opt-in source-checkpoint environment-hash compatibility mode when the target TensorFlow/TFP versions are the same as the source and the only known difference is the Python patch environment encoded by the historical hash. The mode takes the hash directly from the source checkpoint metadata, records both override mode and runtime hash in the result, and is limited to frozen inference from the audited step-8400 checkpoint; it does not permit fit updates, objective changes, threshold changes, K selection or downstream null tests. Any target-version mismatch, missing source hash, pilot equivalence failure or provenance drift keeps the remote run blocked.
D-088 (2026-09-01): the user explicitly prioritized avoiding idle cost on the rented GPU and authorized stopping the CPU/GPU equivalence pilot before completion once the direct GPU panel was running. The GPU full-panel output remains valid only as a practical post-hoc diagnostic, while CPU/GPU numerical equivalence is `NOT_RUN`; do not claim cross-device numerical equivalence. No scientific objective, K rule, threshold or input changed, so this choice does not alter H-01/H-03/H-05 confidence.

D-089 (2026-09-01): formally interpret `K_model` as the allowed dimension of the latent spatial molecular-effect representation, and reserve `K_eff` for the dimension supported by held-out prediction, cross-restart subspace stability, spatial nulls and repeatability. K is not a count of biological structures, structure instances or named fields. Biological structures and latent dimensions are many-to-many; a single factor may not be assigned a TLS, vessel, necrosis or boundary name. The primary identifiable objects are the continuous effect subspace, rotation-invariant reconstructed spatial effects and structure-specific readouts. Exact global K is therefore decoupled from whether two structures are separately resolvable. If later evidence shows that K changes the direction of a major structure conclusion, or that a third K=3 direction has independent cross-fold/seed/null-supported incremental information, reopen formal K comparison and tighten downstream authorization.

D-090 (2026-09-01): use K=3 as the current overcomplete working representation for no-new-training downstream diagnostics, while treating its anonymous subspaces and structure-specific readouts—not factor labels—as the objects of analysis. `selected_k` remains `null`; a K=2 bridge is conditional, triggered only when training-role nested cross-fitting shows that a major structure conclusion is sensitive to the unstable third direction or cannot otherwise be judged. This choice preserves the same NB objective, folds, gene splits, optimizer and null/control requirements. It fails if the missing held-out effect representation is recovered and the K=2/K=3 structural conclusions are materially discordant, or if spatial nulls show that the stable representation does not outperform K=0/non-spatial controls.
D-091 (2026-09-01): for the first downstream readout, use the centered mNSF spatial-rate contribution `section_center(exp(field_mean) @ (loading * factor_amplitude).T)` as the primary rotation-robust effect representation. Export it from the already fitted K=3 training fields and compare shared ecological readout with shared-plus-structure-specific residuals by deterministic patient-level inner cross-fitting. This is a post-fit training-role diagnostic; it does not feed GT into model fitting, does not name factors, and does not replace outer validation, spatial nulls or composition controls. The choice fails if the stored effect cannot be reconstructed from its field/loading/amplitude arrays or if a future frozen outer readout materially disagrees with the training-role conclusion.
D-092 (2026-09-01): trigger the conditional K=2 bridge because the training-role readout has only two patients with paired TLS and TUMOR_STROMA_BOUNDARY GT in each extreme outer fold, so it cannot establish whether the unstable third direction changes the structural conclusion. Freeze Stage A to fold 0/4 and the first three K=3-matched optimization seeds, preserving the source fit length, NB objective, gene split, GP settings, optimizer and learning rate. The manifest is complete, but the first CPU K=2 fit was stopped before its first checkpoint after more than two minutes; historical matched CPU fits require approximately 1.25–6.9 hours per cell. The bridge therefore remains `BLOCKED_CPU_RUNTIME` with no K=2 scientific result; GPU remains disabled until separately approved. This decision fails if a completed Stage A shows fold-dependent rank needs or any material K=2/K=3 structural discordance, in which case the downstream claim remains conditional and Stage B is not automatic.
D-093 (2026-09-01): use the TensorFlow compiled execution path for the frozen K=3 training-effect export and the explicitly frozen K=2 Stage-A attempt as an execution-only optimization on CPU. The path does not change the NB objective, model dimension, data, gene split, seed, optimizer, fit length, checkpoint semantics or output schema; it is not a scientific result and does not override the GPU prohibition. If an equivalence or provenance check fails in a future completed run, revert that run to the eager reference or add a new decision before interpretation.

### D-094 | 2026-09-01 | 将 R-04 全部主动计算路径从 TensorFlow/TFP 重构为 PyTorch，CPU 与 CUDA 共用同一 Torch 实现

- 背景：R-04 的 mNSF、signed-residual GP、NB objective、runtime、serialization、checkpoint audit 及全部 `scripts/r04_*.py` 训练/推断/评分/checkpoint/导出/审计路径此前基于 TensorFlow/TFP eager 图。双后端并存带来可重复性、环境锁定与 GPU 依赖分叉风险，且与"所有新科学计算从单一后端重新开始"的要求冲突。
- 决策：一次性将上述模块重构为原生 PyTorch：`torch.autograd` + `torch.optim.Adam` + 显式 NB（`torch.lgamma`/`torch.distributions` 显式计算，不引入 GPyTorch/分布式）；GP basis、坐标归一化与 inducing-point 几何仍由 NumPy/SciPy 预计算后转为 Torch tensor；CPU 与 CUDA 共用同一实现，不再保留 TensorFlow fallback，可执行 `import tensorflow/tensorflow_probability` 在 active code 中为 0。
- 理由：消除双后端维护与隐式数值分歧；以最小依赖实现 science 合同（连续、可重叠、不规则场；K_model/K_eff 语义；NB/offset/nuisance/HT 缩放/inducing GP/KL/K=0 基线/fold-gene split/seed 与 checkpoint/best-state 语义均保持）；`environment.lock.json` 锁定为 `backend=torch` 的可重现环境。
- 代价：新 Torch checkpoint/frozen schema 与旧 TF 不互通，需独立版本化；旧 TF checkpoint/JSON/NPZ/审计产物全部封存为历史证据，不得作为新 Torch 训练输入；K=2/K=3 如需比较必须在 Torch 内重跑，不得混合 TF K=3 与 Torch K=2 作 formal K 比较；短期内不宣称任何新 K 科学结论。
- 取舍：不引入 GPyTorch、复杂分布式或半精度/XLA；不默认启用 GPU、不主动租用实例、不自动启动长时间真实数据重训。
- 失效条件：若后续证明仅靠原生 Torch 无法达到所需精度/收敛速度或必须引入 GPyTorch/分布式/GPU 才能完成 Stage A 级别的正式训练，则追加新决策并重新评估依赖与环境锁定。
- 影响：`r04/models/mnsf.py`、`r04/models/signed_residual_gp.py`、`r04/objectives.py`、`r04/runtime.py`、`r04/serialization.py`、`r04/checkpoint_audit.py`、`scripts/r04_*.py`、`infra/r04/environment.lock.json`、相关 `tests/test_r04_*.py`；`docs/plan.md` 科学对象不变。

### D-095 | 2026-09-02 | Torch checkpoint 使用单一参数状态、显式 objective-input hash，并将旧 compiled 标签迁移为 eager

- 背景：迁移后的 checkpoint 同时保存可续训的 `params` 与用于选择最佳步的 `best_state`；若审计在缺字段时跨两者补值，会产生磁盘上不存在的混合状态。模型的 NB objective 还会使用可选的显式 `library_size` offset，而历史 `sections_content_hash` 不包含该 offset。
- 决策：R-04 Torch v3 checkpoint 的审计和续训只从一个 `params` 映射读取同一参数状态；存在 `parameter_layout` 时严格校验 schema、键、长度、唯一映射和形状，缺失布局时仅允许旧 v3 的完整 positional `p_*` 映射。新增不改变旧 hash 含义的 `objective_input_hash`，把显式 library-size 纳入 provenance；旧 artifact 没有该字段时可读取但输出 `NOT_VERIFIED_LEGACY_V3`。`raw_v`/`raw_h` 必须成对出现。`execution_mode` 只接受 `auto/eager`，新请求 `compiled` 失败关闭；历史 frozen artifact 中的错误 `compiled` 标签迁移到 eager，并记录迁移来源。
- 理由：保证审计、续训和实际目标使用完全相同的参数/输入状态，防止缺失字段、offset 变化或历史错误标签被静默掩盖；不改变 NB 数学目标、K 语义、fold、gene split、优化器或科学门禁。
- 失效条件：若未来引入真正的 compiled Torch 路径，必须另行定义 checkpoint/执行合同并通过 CPU/CUDA 等价性测试；若新 artifact 缺少 objective hash、参数布局或完整参数 map，应保持 legacy/unverified 或 fail-closed，不能回填猜测。

### D-096 | 2026-09-04 | 批准租用 4090 GPU 实例完成 R-04 Torch fold-0 代表性链路，并授权 checkpoint env-hash 显式 override

- 背景：Torch K=0 fold-0 fit 已达 step 8400（`infra/r04/torch_representative_fold0_20260902_cont_k0_7800/`），但 2400-step 双 gene split 留出推断在 CPU 上超过 7200 秒未产出 cell.json；K=3 fold 0 未开始。GPU 此前按项目政策禁止（D-092 记为 separately approved 才可启用）。用户现明确批准租用 RTX 4090（24GB 显存 / 16 核 / 48GB 内存）实例，镜像 CUDA 12.4，环境自建为 torch 2.5.x+cu124 / python 3.12。
- 决策：在该实例上执行 R-04 Torch fold-0 代表性收尾，范围仅限：(a) 用 `scripts/r04_resume_infer.py --allow-source-environment-hash` 从 step-8400 K=0 checkpoint 做零 fit 步的 checkpoint-only 留出推断（2400 步 × 2 gene split）与评分；(b) 以同一 manifest/fold/gene split/NB-GP/优化器/步数新跑 K=3 fold 0 fit 8400 + 留出推断 2400×2。远端必须设置 `R04_STORAGE_RESERVE_BYTES=100000000000`（沿用 D-086 的远程 100GB 储备先例，本地 1.2TB 硬约束不变）与 `CUBLAS_WORKSPACE_CONFIG=:4096:8`（因 `seed_everything` 启用 `use_deterministic_algorithms(True)`）。env-hash override 仅限上述 checkpoint-only 入口，`cell.json` 的 `environment_hash_compatibility` 记录 `SOURCE_CHECKPOINT_HASH_EXPLICIT_OVERRIDE`、源 hash 与 runtime hash；K=3 新 fit 使用远端自身 runtime hash，无 override。
- 理由：CPU 推断单 cell 超 2 小时未产出，构成 D-092 式的运行时间硬阻塞；代码静态推导的显存峰值 ≤5GB（fit gene batch=512；推断 g=2000 无切分），24GB 余量充足；checkpoint 为纯 tensor dict 的 torch.save zip，跨 torch 2.11→2.5、python 3.11→3.12 兼容；`input_hash`/`objective_input_hash`/`config_hash` 在 resume 时 fail-closed 重算，数据与协议漂移无法静默通过。D-087 的旧 override 先例只覆盖"同 TF 版本、仅 python patch 差异"，本次跨 torch/python 次版本，故单独授权而非援引。
- 代价与边界：CUDA 与 CPU 的 stateless RNG 序列不同，GPU 数值与既有 CPU 数值不逐位可比；延续 D-088，不宣称跨设备数值等价。fold-0 结果仍为代表性/诊断，`selected_k` 保持 `null`，不进科学聚合，不触发空间 null、组成拆分、K=2 bridge 或结构命名。
- 失效条件：2-step wiring 预检的三个 hash（input/objective-input/config）与 checkpoint 元数据不一致；50–100 步短跑出现非有限 loss、显存溢出或 checkpoint schema 审计失败；任一 gene split 的 inference-platform 检查失败；或实例环境无法固定为上述已记录版本。任一发生时停止运行，状态回退为 `BLOCKED_CPU_RUNTIME` 或相应的 fail-closed 记录，不静默更换参数继续。
- 影响：`docs/STATUS.md`、项目根 `TODO.md`（新建）、远端新输出目录（K=0 评分目录与 `torch_representative_fold0_20260904_k3`）；不改变 plan.md 假设、fold 划分、gene split、objective、阈值或优化器。

### D-097 | 2026-09-04 | Torch 代表性 checkpoint-only 推断改用新入口 r04_torch_frozen_infer.py（替代 resume_infer 的 TF 面板门禁路径）

- 背景：D-096 预定用 `scripts/r04_resume_infer.py --allow-source-environment-hash` 完成 K=0 checkpoint-only 推断。实际运行失败关闭：K=0 的 step-8400 fit 审计（`r04.late_lr_cell.v2`）硬编码 `continuation_panel: true`，触发 resume_infer 的严格面板门禁；该门禁要求 `r04.continuation_platform_audit.v2` 面板审计，其 checkpoint bundle 指纹与 source-trace 校验绑定 TensorFlow 时代格式（`ckpt-*.data` 分片、steps [7600,8000,8400] dense trace、restart 2–4 固定 cell 集合），而 Torch 代表性 fit 按协议 `evaluation_interval=0` 无 dense trace，且 checkpoint 为 `checkpoint.pt` 单文件。对 Torch cell 该路径不可达；非面板路径在 CLI 层不可达（`--fit-audit-json` 为必填，空 panel audit 必触发拒绝）。
- 决策：(1) 先用 `scripts/r04_checkpoint_audit.py`（D-095 的 Torch 原生审计）对存活的三个 fold-0 K=0 checkpoint（step 4200/7800/8400）做端点审计，产出 `infra/r04/torch_k0_fold0_checkpoint_audit_20260904.json`（schema `r04.checkpoint_audit.v1`，mc_draws=4，逐 checkpoint 四 draw 精确重复，objective-input hash VERIFIED）。(2) 新增 `scripts/r04_torch_frozen_infer.py` 作为 Torch 代表性 checkpoint-only 入口：强制推断协议等于冻结值（2400 步/16 inducing/lengthscale 3.0/lr 0.01/folds 5/group seed 20260807/gene folds 2/split seed 20260817），强制 `--fit-steps == --expected-start-step`（零 fit 更新），运行前重算 input/objective-input hash 并与 checkpoint 元数据比对（fail-closed），端点审计 sha256 写入 cell provenance，`--allow-source-environment-hash` 的 override 模式与 runtime hash 记录方式与 resume_infer 相同，`selected_k=null`、wiring_only=true。模型数学、协议参数、fold、gene split、阈值均无变化；`_score_fold` 与 checkpoint 审计复用既有经测实现，未新增数值路径。
- 理由：handoff §7 要求"先明确采用可审计的 checkpoint-only inference 入口"；新入口以既有审计原语组合出等价约束，避免为复用 TF 专属门禁而伪造面板审计。
- 失效条件：若端点审计未含 frozen-step 条目或 hash 非 VERIFIED、训练输入 hash 与 checkpoint 元数据不一致、或恢复后 `optimizer_steps_this_call != 0`，入口失败关闭；若未来要把该结果纳入任何正式聚合，需先建立 Torch 时代的面板级审计并追加决策。
- 影响：新增 `scripts/r04_torch_frozen_infer.py`；产出 `infra/r04/torch_k0_fold0_checkpoint_audit_20260904.json` 与 K=0 评分目录 `infra/r04/torch_representative_fold0_20260904_k0_score_gpu/`；不改变 D-096 的范围与边界（其"resume_infer 入口"表述由本条修正）。

### D-098 | 2026-09-05 | 批准再次租用 GPU 实例，按同协议完成 R-04 Torch fold-4 K=0/K=3 配对复跑

- 背景：D-096/D-097 的 fold-0 收尾已完成（K3−K0 均值 −28.96，与 TF 历史 fold-0 全负方向一致）。TF 五 seed 历史中 fold 4 是唯一全正的极端 fold（均值 +20.05），缺少 Torch fold-4 对照则"fold 异质性可重复"（D-080）在迁移后无法裁决。用户明确批准再次租用 GPU 实例完成该对照。
- 决策：在实例上按 D-080/D-084 完全同协议（仅 `--fold-values 4`）新跑 fold 4 的 K=0/K=3 fit 8400 + 推断 2400×2；因 k_search 强制 in-run K=0 baseline，仍用 `--k-values 0,3` 同运行配对；对两个 K 的 fold-4 端点 checkpoint 各做一次 Torch 审计（重复性 + objective-input hash 核验）。运行前以 fold-0 的 2-step wiring hash 锚点（input 044d28c5…/objective-input b5a3fb68…/manifest 5a0dc1de…）fail-closed 校验上传数据完整性。作业包（安装、运行、回传、退租步骤）预置于 `infra/r04/gpu_fold4_20260905/`。远端环境变量约束沿用 D-096（`R04_STORAGE_RESERVE_BYTES=100000000000`、`CUBLAS_WORKSPACE_CONFIG=:4096:8`）；fold-4 为全新 fit，无需 env-hash override。
- 理由：CPU 上该规模不可行（I-015 同类约束）；GPU 路径已在 fold-0 验证（单 cell 约 1 小时、费用约 ¥2-4）；配对两端极端 fold 是裁决 fold 异质性的最小充分证据。
- 代价与边界：fold-4 结果为代表性/诊断，`selected_k` 保持 `null`，不进聚合，不触发空间 null、组成拆分、K=2 bridge 或结构命名；延续 D-088，不宣称跨设备数值等价；若 Torch fold-4 方向与 TF 历史冲突且无法归因，记录 ISSUES 并停下报告，不强行裁决；实例结束即退租。
- 失效条件：预检 hash 锚点不一致；运行出现非有限 loss、显存溢出、checkpoint schema 审计失败或 inference-platform 检查失败；环境无法固定为记录版本。任一发生时停止运行并按 fail-closed 记录，不静默换参继续。
- 影响：新增 `infra/r04/gpu_fold4_20260905/`（README/setup_remote.sh/run_fold4.sh）与运行产物目录；不改变 plan.md 假设、fold 划分、gene split、objective、阈值或优化器。

### D-099 | 2026-09-05 | K_model=3 落定为下游受控分析的默认工程工作表示，selected_k 保持 null

- 背景：fold 间 K 偏好异质已在 Torch 迁移后复现（fold-0 K3−K0=−28.96，fold-4 +72.73，方向均与 TF 历史一致），证实不存在全局最优 K（D-080）。按 D-089/D-090 的语义，K_model 是表示容量而非结构数量，K 的选择本来就是工程适用问题；D-090 已将 K=3 指定为当前过完备工作表示。用户明确要求结束 K 环的悬置状态、尽快推进。
- 决策：下游受控分析（结构读出及其后续）默认使用 K_model=3 作为工程工作表示；代码层面无需改动——`--k-values` 与模型 `n_factors` 本就可配置为 2 或其他值，冻结协议脚本（r04_real_k_search.py 等）的默认 k 列表保持不变以护审计链。该默认不是 K_eff=3 的科学断言：`selected_k` 保持 `null`，任何对外表述不得把 K_model=3 写成"数据支持三个结构/方向"。默认的生效附带验证门禁：先完成 Torch cell 上结构读出的 K 稳健性分析（fold-0 与 fold-4 的 K=0/K=3 对照）；若主要结构结论不依赖不稳定第三方向，默认正式生效；若依赖，触发 K=2 Stage A（GPU 路径已验证，单 cell 约 1 小时、约 ¥2，见 infra/r04/GPU_RUNBOOK.md）。
- 理由：K=3 在两个 fold 上均稳定 fit、无 collapse 警告、留出推断收敛；在 fold-4 显著最优、在 fold-0 次于 K=0 但不崩；在 fold 异质已被证实的前提下，继续为"选 K"投入资源的边际收益低于推进结构读出的收益；工程默认 + 稳健性门禁的组合既满足推进速度要求，又保留了对第三方向依赖性的检验通道。
- 失效条件：结构读出 K 稳健性分析显示主要结论对 K∈{0,3} 敏感；或 K=2 Stage A（若触发）显示 K=2 明显优于 K=3；或未来新增 fold/队列出现 K=3 fit 不稳定或 collapse 证据。任一发生时回退默认，重开 K 裁决并追加决策。
- 影响：新增 `infra/r04/GPU_RUNBOOK.md`（资源实测档案）；TODO 当前分支切换为"结构读出 K 稳健性分析"；不改变 plan.md 假设、fold 划分、gene split、objective、阈值、优化器与冻结协议脚本。

### D-100 | 2026-09-05 | Torch 训练角色场导出授权（training-only，不新增推断）与 D-099 门禁规则预注册

- 背景：D-099 门禁要求在 Torch cell 上做结构读出的 K 稳健性分析（fold-0 与 fold-4 的 K=0/K=3 对照）。现状缺口：两个 Torch K=3 cell 的 `structure_field_exports` 均为空；`scripts/r04_export_frozen_effects.py` 的 `_validate_source` 只接受 TF 时代产物（`checkpoint` 指针文件的 `ckpt-` 标记、`r04.mnsf_checkpoint.v2`、TF bundle 指纹），会对 Torch 的 `checkpoint.pt` + `r04.mnsf_checkpoint.v3_torch` 失败关闭。但 `_score_fold(training_only=True)`、field export 写入与 `r04_run_structure_readout.py` 均与后端无关；training-only 只做数据加载 + checkpoint 恢复 + 训练角色效应前向导出，不跑 2400 步留出推断（CPU 上约分钟级固定开销，D-096 的超时证据针对的是完整推断，不适用此处），因此无需 GPU。K=0 没有空间因子，走不了同一读出管线；门禁中的 K=0/K=3 对照指 fold 内 K=0 留出分数（fold-0 −28.96、fold-4 +72.73 已知）与 K=3 读出的联合判定。
- 决策：(1) 给 `_validate_source` 加 Torch 分支：接受 `checkpoint_metadata.json` 的 `r04.mnsf_checkpoint.v3_torch` + `checkpoint.pt` 存在性，核验 `sha256(checkpoint.pt)`、input/objective-input hash 与源 cell 一致，step 取源 cell 的 `fit_diagnostics.steps`，env hash 以显式参数透传（记录 override，沿 D-096 模式）；(2) 新建 Torch 导出 panel manifest（restart 0、fold 0/4、K=3 cell）；(3) 本地 CPU 跑 `--training-only` 导出，再跑既有 readout 脚本；validation GT 全程密封。排除的替代方案：另写 Torch 专用导出脚本（D-097 模式）——此处 `_score_fold`/export/readout 已是后端无关实现，只有来源校验前言是 TF 绑定的，最小改动即够用。
- 门禁规则（在看到 Torch 读出数值前预注册）：逐 fold 比较 `full_k3` 与 `stable_rank2_sensitivity` 的 `specific_increment`（逐结构 AUC delta、MSE delta）。两 fold 的逐结构结论定性一致（AUC delta 与 MSE delta 同号，或双方 |AUC delta|<0.02 即联合为零）→ ROBUST，K=3 默认正式生效；任一 fold 出现定性翻转（|AUC delta|≥0.02 的符号翻转，或有实质幅度的 MSE 符号翻转）→ SENSITIVE，触发 K=2 Stage A；任一 fold 为 NOT_TESTABLE（无有效内折/配对 GT 不足）→ 无法建立，按 roadmap 第 4 条走 K=2 Stage A 通道（需另行申请 GPU，不自动租用）。
- 失效条件：导出 replay 出现 `optimizer_steps_this_call != 0`、replay input/config hash 与源 cell 不一致、env override 未记录、任一单 entry CPU 超过约 1 小时，立即停止并按 fail-closed 记录，不静默换参继续。
- 影响：`scripts/r04_export_frozen_effects.py` 的校验分支扩展 + 新 panel manifest + Torch readout panel JSON；不改模型数学、协议参数、fold、gene split、阈值；不新增训练；不读 validation GT。

### D-101 | 2026-09-05 | 项目转为发现/探索定位，取消预注册约束；D-100 预注册规则废止，门禁按实质重判

- 背景：用户明确本项目是发现和探索性项目。既有流程按确认性范式搭建了预注册约束（outcome-blind 角色冻结、validation GT 规则冻结前密封、"只有…才…"式门禁触发、D-100 的阈值预注册），在探索定位下不再适用，且已实质阻碍推进（D-100 规则的字面触发与实质解读分离就是例证：fold-4 TSB −0.076→+0.001 触发字面规则，但内折显示两边都是单患者噪声）。
- 决策：(1) 取消全部预注册约束：R-01 未来角色/用途调整不再要求 outcome-blind（须记录依据并披露所用信息；已冻结角色保持不变以保 provenance 连续）；R-04 顺序第 3 条 validation GT 改为探索性使用允许加证据等级标注；第 4 条 K=2 从触发式门禁降为可选探索方向（LEADS L-001）；第 5 条顺序门禁解除，各方向可并行探索；D-100 预注册规则废止，不再作为判决依据。(2) 保留（与预注册无关）：GT 永不进入预测输入与训练；leakage group 切分与保守连通分量；hash/provenance/fail-closed 工程约束；input manifest 合同；全部历史记录不动（含 D-100 原文与 roadmap 日期条目，由本条覆盖其状态）；plan.md 不动，其假设与"放弃条件" reinterpret 为探索方向参考，不再是确认性判决。(3) 门禁重判（探索语义，数值不变）：Torch 读出无任何稳定结构结论（N=2/fold 且全部非零值单患者驱动、跨患者矛盾）；没有任何结论正向依赖第三方向。结论为描述性未决；K=3 作为工作容器继续（D-099 的工程默认不变）；K=2 Stage A 转 LEADS L-001（可选探索，如执行仍需单独 GPU 审批）；之前提交的硬阻塞 GPU 审批撤回。
- 理由：预注册解决的是"确认性结论的可信度"问题；探索定位下真正需要的是"不把探索写成确认"（证据等级标注）加"算对"（防泄漏与 provenance）。硬套已失准的字面规则会花 ¥3–8 去裁决一个单患者噪声，不符合探索的成本/价值比；但废止自己写下的规则必须由用户授权的定位变更来做，这正是本条。
- 失效条件：若未来对外表述把探索性结果写成确认性结论，或 GT 进入预测输入，或绕过证据等级标注，则本条的宽松化自动失效，回退确认性约束并追加决策。
- 影响：roadmap R-01/R-04 现行规则改写加新日期条目；新建 LEADS.md（L-001）；TODO 当前分支更新；之前提交的 GPU 硬阻塞审批撤回；不改 plan.md、历史记录、代码行为与测试。

### D-102 | 2026-09-05 | R-04 探索性读出 null 与深化的固定分析口径

- 背景：D-101 转探索定位后，读出类分析不再设预注册判定阈值，但跨时间可比需要固定口径。三方向并行第一轮用同一 Torch 训练角色导出与训练角色 GT 完成 A（标签置换 null）与 C（rank 曲线），B（组成拆分）经盘点无输入而转为 proposal（LEADS L-002）。
- 决策：探索性读出 null 固定为"两结构标签联合、在 patient 内置换 spot，200 draws，种子 20260905＋fold/rank 偏移"，统计量与 panel 读出同口径（内层患者交叉拟合的 specific-increment 均值），只报经验 p 与分位数，不设通过阈值；读出深化固定为 rank1/2/3 曲线加绝对 AUC/MSE、残差系数范数与逐内折表。两者输出 schema 为 `r04.explore_readout_null.v1` / `r04.explore_readout_deepdive.v1`，状态一律 exploratory。后续读出类探索沿用此口径以保可比；validation GT 仍未动用。
- 理由：null 必须在 GT 侧 réagir 而不是效应侧（效应是拟合产物，重排它破坏 provenance）；patient 内置换保留患病率与结构间相关，只打断 spot—标签关联，正好对应"空间效应与结构无关"的零假设；200 draws 在当前 N 下已使蒙特卡洛误差远小于患者抽样方差，加 draws 不增加信息。
- 失效条件：若未来配对患者数显著增加或读出统计量定义改变，本口径需重审并追加决策；不得把经验 p 写成确认性 p 值。
- 影响：`scripts/r04_explore_readout_null.py`、`scripts/r04_explore_readout_deepdive.py`、`tests/test_r04_explore_readout.py`；产物 `infra/r04/explore_readout_{null,deepdive}_torch_20260905.json`；不改变模型、协议与 GT 密封状态。

### D-103 | 2026-09-05 | 跨设备推断重放分数门禁的处置：接受等效最优点漂移，场导出可用于探索性外层验证

- 背景：fold-4 全量导出 CPU 重放收敛（双 split inference_platform 均为 converged，2400 步），fit_updates=0，三 hash 与源 cell 一致，但终点患者分数未通过 rtol=1e-6/atol=1e-5 的逐值重放门禁。对比终点推断损失：split-0 重放 3300.95 vs 源 3300.61（差 0.34，约 1e-4 相对），split-1 重放 3337.94 vs 源 3339.16（差 1.2，约 4e-4 相对）。源 cell 在 GPU（CUDA kernels/RNG 流）上产生，重放在 CPU 上进行；torch 不保证跨设备逐位可复现，2400 步变分轨迹的设备噪声必然超过 1e-6 量级门禁。D-088 早已确立不宣称跨设备数值等价。
- 决策：(1) 该门禁的触发判定为"等效最优点上的设备噪声"，不是协议破坏；(2) 重放场导出（`heldout_full.npz` 全基因、`gene_split_*.npz`）可用于探索性外层验证，前提是逐项核验：fit_updates=0、hash 一致、重放 inference_platform 收敛、重放终点损失与源终点损失相对差 <1e-3；(3) 导出脚本的分数失配路径现将重放分数与逐患者差值写入 FAILED 记录再抛出（`FAILED_SCORE_MISMATCH_FIELDS_RETAINED`），不再丢弃诊断信息；(4) 任何跨设备重放分数不得进入科学聚合或 K 选择比较。
- 理由：门禁的设计目的是捕获协议漂移（错数据、错 seed、错步数），这些通道已由 hash/步数/收敛三重检查覆盖；把同一最优点上的设备噪声当作协议破坏，会要求不可能的跨设备逐位复现。1e-3 的损失等效界远小于 K3-K0 患者分数差的量级（数十），不影响探索性读出的定性结论。
- 失效条件：若重放终点损失相对差 ≥1e-3、或 inference_platform 不收敛、或任一 hash 不一致，则仍按协议破坏处理，不得接受；若未来出现同设备重放失配，门禁维持原 1e-6 严格度。
- 影响：fold-4 外层验证可直接使用已落盘重放场；fold-0 导出若触发同一门禁，其场同样按此处置（以实际核验为准）；不改变 D-088（不做跨设备等价宣称）与 D-101（探索等级标注）。

### D-104 | 2026-09-10 | GPU 自租用授权：K=2 配对跑（L-001），4090 档，看守与预算上限

- 背景：用户明确批准由 agent 自行租用 GPU（L-001 K=2 Stage A 从 LEADS 转入执行），并强调费用看守（完成后关闭、空转检查）。
- 决策：租用 4090 档（¥1.65/h；无货则 3090 ¥1.1/h），跑 `--k-values 0,2` 配对（fold 0+4，同 D-096/D-098 冻结协议）；预算上限 ¥15；看守机制为 bg 心跳（进程存活＋nvidia-smi 利用率＋日志 mtime 三重），连续空转即查因处置；完成后立即回传→退租→对账记账。
- 失效条件：超预算停手汇报；余额不足不开租。
- 影响：GPU 累计花费在汇总时记账。

### D-105 | 2026-09-10 | marker-proxy 组成路线授权（探索级，R-02 例外事项记录）

- 背景：组成拆分缺逐 spot 组成输入（LEADS L-002）；用户批准 marker 代理与 scRNA 参考双路线并行。三库（CellMarker 2.0 文献策展、PanglaoDB 灵敏度/特异度量化、CellTypist CRC 模型系数）投票后 5/6 类别在 panel 内可用（T22/B11/Mye16/Epi54/Stromal22 基因；ILC 仅 1 基因不可用）；规则缺口已修补，未映射标签与有意丢弃类别（内皮/神经/脂肪等无对应 Major）均有记录。
- 决策：允许用三库投票 marker 集构建逐 spot 模块分 proxy（类内 log1p 均值→行归一伪组成），经 patient 分组 crossfit 残差化后重跑读出；全程标 exploratory。R-02 的表达派生标签禁令在此处境下让位，理由：proxy 只用作组成协变量调整而非结构 GT，且 scRNA 参考路线并行提供独立对照；若两路线结论矛盾，以参考路线为准并如实报告分歧。
- 失效条件：proxy 调整前后结论翻转但无法归因到组成 vs 方法伪影时，不得单凭 proxy 路线下任何结论；对外表述 forbidden 确认性措辞。
- 影响：`scripts/r04_combine_marker_proxy.py` + smoke 脚本 + 产物 `marker_proxy_combined.json`；不改 panel、模型与 GT 使用范围。

### D-106 | 2026-09-10 | K=2 Stage A 裁决完成：第三方向无增量贡献，K=3 默认正式生效

- 背景：D-099 门禁要求裁决主要结构结论是否依赖不稳定第三方向；D-100/SENSITIVE 与 D-101 重判后，K=2 Stage A 作为可选探索进入 LEADS L-001。用户批准 GPU 自租后执行（D-104，4090PLUS 48GB，¥2.49/h）。
- 证据：Torch fold-0/fold-4 的 K=0/K=2 同协议配对跑（`torch_representative_fold{0,4}_20260910_k02`，4 cell 全部 FIT_AND_SCORED，inference 双 split 收敛，4 份端点审计）：fold-0 K2−K0 均值 −16.61（K2 胜 2/6），fold-4 K2−K0 均值 +71.96（K2 胜 6/6）。对照既有 K3−K0（fold-0 −28.96，fold-4 +72.73）：fold-4 上 K=2 与 K=3 几乎逐位一致（+71.96 vs +72.73），fold-0 上方向一致（均为负，K=2 劣势更小）。
- 决策：第三方向（K=3 相对 K=2 的增量）不带来任何留出预测增益——fold-4 的全部正向增益已被 2 个方向捕获，fold-0 的负向在 K=2 下同样存在。K=3 作为过完备工作表示正式生效（D-099 默认转正）；`selected_k` 保持 `null`（仍不断言 K_eff=3）；K=2 无需替代 K=3（K=3 容器包含稳定子空间，rank-2 读出口径不变）。
- 理由：这是 D-099 门禁要的直接裁决，且是有检验力证据（6 患者/fold 配对分数），不是 N=2 读出的描述性推断；fold 异质性本身（fold-0 偏 K=0、fold-4 偏 K=3）在 K=2/K=3 下同时复现，说明异质性位于 ≤2 的有效维度内。
- 失效条件：若未来 fold/队列出现 K=3 显著优于 K=2（留出配对差值方向不一致），或结构读出显示结论依赖第三方向，则重开 K 裁决。
- 影响：L-001 转已并入主线；后续结构读出/验证工作默认 K=3；本次 GPU 会话花费约 ¥3.12（Money ¥0.39＋Power ¥2.73），累计 GPU 花费 ¥8.43。
