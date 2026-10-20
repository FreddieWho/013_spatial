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

### D-107 | 2026-09-11 | K closure 精确语义：D-092 six-cell 范围被 D-104/D-106 取代，K 搜索关闭

- 背景：D-092 曾冻结 six-cell（fold 0/4 × restart 0/1/2）K=2 Stage A 并记为 `BLOCKED_CPU_RUNTIME`；D-101 转探索定位后，D-104 授权更小的 GPU panel（fold 0/4 × restart 0，`--k-values 0,2` 两次配对跑），D-106 以该 panel 关闭 K 问题。但 `scripts/r04_finalize_k_semantics.py` 与 `infra/r04/k2_bridge_stage_a_20260901/stage_a.json` 仍停留在旧状态，机器入口与文档"已关闭"矛盾。
- 决策：(1) 对 D-106 的精确解释：第三方向没有稳定正增量（fold-4 K2≈K3：+71.96 vs +72.73；fold-0 K2−K3 均值约 +12.35，即 K2 劣势更小，但 K2/K3 均不优于 K0）；采用 `working_k_model=3 + primary_readout_rank≤2`，但不声称 `global K_eff=2`，`selected_k` 保持 `null`。(2) D-104/D-106 的 exploratory GPU closure 正式 supersede D-092 的 six-cell 执行范围：不再补跑缺失 seed；压缩的依据是 K=2/K=3 配对分数已对第三方向做出有检验力裁决（6 患者/fold），额外 seed 只能改变幅度估计，不能改变"第三方向无稳定正增量"的定性结论；该压缩只支持"停止 K 搜索/选定工作表示"，不升级为确认性全局 K 结论。(3) 旧 bridge driver 与 stage artifact 保留为历史（legacy），`r04_finalize_k_semantics.py` 升级为 canonical K closure 入口，直接绑定 K=2/K=3 evidence artifact 与 hash。
- 失效条件：若未来 fold/队列出现 K=3 显著优于 K=2，或读出依赖第三方向，则重开 K 裁决（沿 D-106 失效条件）。
- 影响：`scripts/r04_finalize_k_semantics.py` + 测试；`infra/r04/k_closure_canonical_20260911.json`；旧入口标 legacy。

### D-108 | 2026-09-11 | Composition adjustment 改为 nested patient-level cross-fitting，旧 smoke 产物保留为探索历史

- 背景：`scripts/r04_explore_marker_smoke.py` 先在全体 paired 患者上 pooled 拟合 `crossfit_composition_adjustment`，再做内层 patient split 读出。内层训练患者的 adjusted 特征经由含内层测试患者数据的 residualizer 计算，构成 transductive/nested-CV 泄漏（2 患者 LOPO 下：测试患者的 adjusted 值虽只用训练患者拟合，但训练患者的 adjusted 值用了测试患者拟合，读出映射仍被污染）。
- 决策：final audit 把 residualizer 嵌套进每个内层 split：只用内层训练患者拟合 `effect ~ ilr(proxy)`（单训练患者时用普通 ridge、无 crossfit），同一映射变换训练/测试；测试患者不参与任何 residualizer 参数估计。新增泄漏回归测试（向测试患者注入 composition-only 信号，旧路径必须检出、新路径必须干净）。旧 `explore_marker_smoke_fold0_20260910.json` 保留不动；final audit 写新产物。另记录仓库事实冲突：现行 roadmap/phase-report 写"full-TSB +0.127→−0.068"，但产物中 a29-test-fold 的 full-TSB unadjusted 实际为 −0.011（+0.127 属 rank2-TSB unadjusted）；新审计以产物值为准，历史文字不改写（见 final report）。
- 失效条件：若 nested 调整后训练/测试映射无法构造（如单患者无 target variation），记 NOT_TESTABLE 而不放宽嵌套；不得回退到 pooled 拟合。
- 影响：`r04/composition.py` 新增 train-only 拟合接口；新脚本 `scripts/r04_composition_final_audit.py`；新测试。

### D-109 | 2026-09-11 | 最终审计规则预承诺（outcome-neutral，结果出来前冻结）

- 背景：R-04 最终门禁必须允许阴性关闭，不能事后调规则。本条在看到 final-audit 数值前冻结两条审计的判定规则。
- 决策：(1) Composition final audit（nested patient-level crossfit，D-108）：rank1/rank2（主）/full-K3（sensitivity），unadjusted vs nested-adjusted，训练角色导出（fold-0/fold-4，不新增训练）；uncertainty 为内层 LOPO 经验 spread＋test-spot bootstrap CI（200 draws，seed 20260911）。residual-signal 记录 adjusted 效应＋CI，不设通过阈值，由最终门禁综合判定。(2) Unnamed-field audit（旋转不变、无标签）：两 fold K=3 loading 的 SVD 内禀 rank-r 子空间 canonical correlation，基因行置换 null（200 draws，seed 20260911）；rank-1/rank-2 候选规则：min(obs cc) ≤ null q975 → DOES_NOT_SURVIVE；> max(null q975, 0.5) → SURVIVES；之间 → NOT_IDENTIFIABLE；rank-3 只描述。(3) fold-4 TSB rank3 0.544：若 rank1/rank2、fold-0 或 nested-adjusted 任一不复现，默认记 ISOLATED_EXPLORATORY_SIGNAL，不追第三方向。
- 失效条件：规则一旦冻结，看到数值后不得修改阈值；若输入 artifact 被证伪，审计作废重跑而非调参通过。
- 影响：`scripts/r04_composition_final_audit.py`、`scripts/r04_unnamed_field_audit.py`＋测试；产物 `composition_final_audit_20260911.json`、`unnamed_field_audit_20260911.json`。

### D-110 | 2026-09-11 | R-04 最终门禁结果：阴性完成，无可重复残余场，R-05 不触发

- 背景：D-109 预承诺规则下跑完 nested composition final audit 与 unnamed-field audit，`scripts/r04_finalize_phase.py` 从 5 份登记产物重算（hash/患者独立性 fail-closed）。
- 证据：nested 调整后唯一非平凡残余是 a29 单患者的 TSB rank2 信号（adjusted +0.157，CI [0.095,0.212]）；独立患者 0bd3（+0.005）与 286667（−0.002）均为零，且 a29 同时在两 fold 训练集中、不能自我复制；TLS rank2 最大 reviewer 内效应 +0.034/+0.006，未过项目 null 等价带；rank1 读写全零（单维度下特异残余空间坍缩所致的结构性零，非证据）；外层验证双 fold null；unnamed rank-1 loading 子空间共享（cc 0.582）但无组成残余/外层支撑，不等于残余场。
- 决策：最终状态 `R04_COMPLETE_NO_REPRODUCIBLE_RESIDUAL_FIELD`（阴性完成，不是 fishing 未果）。合成规则（D-109 把阈值制定委托给最终门禁，此处冻结）：hit 需 adjusted AUC delta ≥ 0.02（项目自 D-100 时代沿用的 null 等价带：|.|<0.02 联合读零；bootstrap CI 只覆盖 spot 抽样噪声，信号必须明显越过该带）且 q025 > 0；SURVIVES 需 ≥2 个不相交患者同号命中。4 候选全部 DOES_NOT_SURVIVE；uncertainty PASS；coordinate null 与 independent lineage 按 §6 记 NOT_TRIGGERED（无 surviving candidate，不是 open TODO）；R-05 记 `NOT_TRIGGERED_NO_SURVIVING_R04_FIELD`，R-06/R-07 不触发。
- 失效条件：若未来独立患者/队列出现跨患者复制的组成残余信号，或外层验证转阳，则重开 R-04 证据评估（新决策，不改写本门禁）。
- 影响：`infra/r04/r04_final_gate_20260911.json`；roadmap R-04 完成判据改为 outcome-neutral 三状态并记最终状态；I-019 关闭，I-012/I-013/I-014/I-017 按门禁记 resolved/limitation/deferred；旧 smoke 结论（a29 信号被吃掉）被 nested 审计推翻，以新产物为准（旧产物保留，历史文字不改写）。

### D-111 | 2026-09-11 | 用户批准 serial-section 数据申请，启动外部搜寻（只侦察不下载）

- 背景：R-03 普查确认已下载数据中零对已知切片间距，R-11 邻近平面三类任务不可执行（I-002 对 R-11 已升级硬阻塞）。用户批准数据申请，并指示用 research agent 搜寻合适数据。
- 决策：批准启动外部搜寻，本轮只做候选清单与逐项核实，不下载任何数据；下载前另凭清单报批。搜寻标准（缺口→不可替代性→规模→解锁判据）：(1) 同一组织块≥2 切片且片间距有书面记录（论文方法/数据集说明/Zenodo 记录均可，文件名暗示不算）；(2) 每切片有 counts＋坐标＋H&E（配准用）；(3) 优先 Visium 人肿瘤（与现有 HTAN/ST-CRC/Block-A 可比），平台不同则须论证可比性；(4) 规模：先求 1–2 个块共≥3 片（含一组已知间距邻片）做可行性验证，不求大队列；(5) 解锁判据：拿到间距后 R-11 邻近平面预测任务由 NOT_MEASURABLE 转为可执行（做不做得成另判）。
- 失效条件：找不到满足 (1) 的数据则如实报告，I-002 硬阻塞维持，不放宽"记录间距"为"推测相邻"。
- 影响：新增侦察任务（TODO #32）；不改变数据边界（下载仍需另批）。

### D-112 | 2026-09-14 | 场发现策略修正：残余降级为评分维度，守门只问方法不问内容，field registry 立项

- 背景：用户连续质疑 R-04 门禁语义：(1) TSB/TLS 为何零存活；(2) 组成本身就是场的一部分，要求"组成之外的残余"作为参赛资格，等于把本质即组成的真实结构挡在门外。TLS 阳性对照（`infra/r04/tls_positive_control_20260912.json`，13,726 标注点／170 阳性／4 片／3 患者，exploratory）证实：经典 TLS 基因（CXCL13、MS4A1/CD20、CD3D/E/G、CD19、CCL19、CCL21、LTB 等）**全部不在** frozen 4000-gene panel 内；可用信号在样本内可测（PTPRC AUC 0.78、组合 logistic 0.93），但最强项恰是组成型信号（白细胞富集＋非上皮）；浆细胞模块在 0bd3 为正（0.81）、在 a29 为负（0.34–0.38，符号翻转）。用户立场：理想路线是不设内容门禁、发现所有场（组成解释的和不可见的），再统一打分解释；agent 防线：取消复现检验，"发现"一词失去真值（a29 单患者信号曾是近误发反例），且 fixed-cell-type state 与"发现即新组成"在缺参考系时不可判定。双方达成路线共识。
- 决策：
  1. "组成之外的残余"从**参赛资格**降级为**证据评分维度之一**：所有候选场（组成可解释的和不可见的）都可以被发现和登记；评分表必须包含残余维度，但残余不单独否决任何场。
  2. 统一证据标准（守门只问方法，不问内容）：升级为 claim 必须跨患者/队列复现、赢过相应 null、不确定性量化；标准统一，但 null 的具体实现允许因场制宜（离散结构与连续场的零假设可不同——统一门槛高度，不统一刻度）。
  3. 新增 **field registry**（scored list）为持久资产：每行一个候选场，列为跨患者复现、过 null、残余维度、效应量、不确定性、证据等级；账本永不关闭，新数据追加行；停止规则内建（本轮候选打完分即收工）。
  4. 工具方向保留并加两条修正：预测器（新切片→场强度/到场距离）必须双头输出（组成解释量＋残余量），且输出自带证据等级（弱证据的场不得裸奔出去预测）。
  5. R-04 阴性关闭（D-110）**不改写不回溯**：它回答的是"残余定义下是否有可重复场"（答：没有）；本决策改变的是**后续**场发现与登记方式。
- 理由：H-03 纠缠问题的工程化解法——防平庸（残余检验）与防错杀（内容放开）不再焊死；TLS 阳性对照证明残余资格会过滤真实生物学；同时 a29 的教训证明无复现检验的全发现会放任噪音过关。此拆分同时满足"全发现"诉求与复现防线。
- 失效条件：统一评分制下出现系统性水货（大量无复现场被当 claim 引用）则收紧内容侧并追加决策；registry 长期回填不出任何可复现场，按 plan 否证逻辑整体收缩。
- 影响：roadmap（R-04 判据补注＋registry 资产要求）；TODO；I-020（panel 覆盖缺口登记）；后续先建 registry schema＋回填；R-11 smoke 与数据申请维持 hold 待批。

### D-113 | 2026-09-14 | 发现引擎定型：组成 pattern 为地基的距离函数框架＋双车道分步 DAG；连续信号 null 机器换代

- 背景：用户质疑"R-04 空转后只发现 4 个信号"。核实确认：仓库**从未存在**空间 pattern 挖掘代码（无聚类/hotspot/Moran's I）；R-01/R-02 结构清单继承自文献注释（atlas Table S4 的 87 个 TLS ID、STOmics 19 份注释），R-04 是锚定读出评估——"4 候选"是评估顺路碰到的对象，不是发现扫描的产出。用户提出发现原理：看空间聚集 gene pattern 在样本间的普遍性；并进一步主张：组成型 pattern 是**地基**，目标是找到与组成 pattern 的**距离构成函数关系**的 pattern（不论细胞组成还是信号组成）；同时指出 Moran's I/置换等开关变量检验对连续信号不适用，要求判断分步还是整合。
- 决策：
  1. **距离函数框架为主路线**：每个登记的组成 pattern 自动生成一把距离尺（TSB 尺的推广），在全基因×全尺上筛连续信号的距离函数关系。函数族预选：等调单调＋限自由度样条；壳层（中间距离峰）以单峰约束单列一类——只筛单调会漏壳层，全放开非参数则多重检验失控。几何选择（到边界/到质心/有符号距离）与分辨率地板（Visium 点距约 100µm，窄于 2–3 点距的梯度物理不可见）写入每个 pattern 的登记行。
  2. **连续信号 null 换代**：现有 patient 内标签置换 null（`r04_explore_readout_null.py`）是开关变量检验，对其原任务（AUC vs 标签）正确，对连续信号不适用；naive 值置换破坏空间自相关使 p 值反保守。连续信号一律使用**空间保持 null**（场相对 mask 旋转/平移、按经验变异函数模拟 GP null、块置换）；检验对象是信号与距离协变量的关联（带空间相关误差的回归/GP 对比），不是全局自相关。Moran's I 本是连续变量检验，但回答的是"是否空间结构"，答错了"是否锚定于某 pattern"这个问题。
  3. **残余维度列保留**：组成本身沿距离平滑变化，"信号随距离变"可能只是"组成随距离变"的重述——D-112 的残余评分维度在距离框架下继续保留，不做参赛资格，只做解释约束。
  4. **分步 DAG，不做一步整合**：没有 mask→没有距离→没有筛选，依赖链强制分步；registry 要求逐 pattern 逐行记分，整合模型把各层不确定性糊在一起（H-03 第三次出现），其唯一真实优点（不确定性传播）降级为敏感性检查保留。
  5. **双车道**：Lane A（锚定）组成 pattern 普查→距离变换→全基因×全尺筛距离函数，跨患者一致性**按病人计数**（a29 教训内建）；Lane B（无锚）回归掉全部距离协变量后筛残余空间结构——H-01 在 D-112 后的操作化，防止框架系统性漏掉不锚定任何已知 pattern 的场。
  6. **执行序列**：①组成 pattern 普查（跨分辨率稳定性入分）→②距离变换＋几何登记→③Lane A→④Lane B→⑤全部进 field registry。①②纯 CPU 复用 panel cache；③④复用改造后的 null 机器。
- 理由：用户"全发现"诉求与复现防线的第二次拆分（继 D-112 后）：发现引擎（普查＋距离筛选）与证据引擎（null＋复现＋记分）解耦；普查只产组成型 pattern 是合法资产而非失败。
- 失效条件：Lane A 的距离函数全部不可跨患者复现→按 plan 否证逻辑收缩；Lane B 在 Lane A 充分后仍只产不可复现残余→同样收缩；普查 mask 跨分辨率崩解→不静默换方法，追加决策。
- 影响：roadmap 新增 R-16 节点；R-15 假设覆盖更新；I-021（连续信号 null 缺口）登记；TODO（先交一页设计方案）；plan.md 不动（H-01 不变，此为其操作化）。

### D-114 | 2026-09-14 | R-16 gene universe 改用整合 HVG top-10000；样本范围定为 30 发现＋16 复现

- 背景：用户追问 4000 基因来历，机理查明（I-020）：D-042 两阶段过滤（prevalence ≥5% spot 且 ≥70% 患者 → log1p 方差前 4000）系统性杀灭"小、稀有、患者间不均"结构的 marker——14 个经典 TLS 基因全灭（12 死于 prevalence、2 死于方差排名），方差榜被 IG/持家/上皮霸占。用户提出用经典 Seurat 整合 HVG 策略重新选基因。pilot（`infra/r16/hvg_recurrence_pilot_20260914.json`）：每片 seurat_v3 top-5000（loess frac=0.3）→ 患者 recurrence 合并（多片患者一票）→ **14/14 经典 TLS 基因进 top-8000**（rec 13–25/30）；天然截点 **recurrence ≥10/30 ≈ top-10000**（10,383 基因）。frozen 4000 在新排名中位数 14,722，两套哲学近乎正交。
- 决策：
  1. **R-16 census universe = 整合 HVG top-10000**（每片 seurat_v3 top-5000 → 患者 recurrence 合并 → 总排名前 10000；≈"≥1/3 患者中复发的高变异基因"）。从源 h5ad 重新 materialize census 专用 cache（含坐标），不改动 R-04 任何冻结资产（D-042 仅作 R-04 历史范围决策保留，不回溯）。
  2. **样本范围（用户批准选项 C）**：发现引擎只用 Vanderbilt training（47 片 / 30 患者，单 lineage 单平台）；ST-CRC（14 片 / 7 患者）与 USZ（8 片 / 8 患者，带 TLS 注释）作为**独立复现引擎**——registry 行欲升级证据等级须在这两个队列中复现，跨平台复现按 D-113 自动降半级记录；Block-A 维持挂起；HTAN WUSTL 55G 在 R-01 注册表之外，不并入。
  3. HVG 排名产物（`infra/r16/hvg_final_rank_20260914.npy`）与 pilot 脚本入库，全确定性；seurat_v3 为 Python 复刻非逐字节 Seurat，回收率（14/14）即其验证记录。
- 失效条件：若发现复刻实现与参考实现（scanpy seurat_v3）系统性偏离，重推排名并审计受影响产物；若跨队列复现因平台效应系统性不可行，按 D-113 降级规则如实记录而非强推。
- 影响：设计方案 v1.1；I-021 细化（自相关检验用值置换合法，关联检验才需空间保持 null）；TODO；随后 materialize＋普查第①步。

### D-115 | 2026-09-14 | R-16 升级为两层结构：Tier-1 连续轴尺子厂＋Tier-2 残余发现（用户批准）

- 背景：普查 v0 暴露逻辑层问题（用户第一性原理追问）：(L1) 硬聚类把"程序"与"混合邻居"焊死，全签名相似度无法区分"程序不同"与"邻居不同"，且跨平台放大；(L2) 用户保留——地基可离散，但旋钮必须披露可扫（Leiden resolution 是藏起来的旋钮：pattern 数 198→488）；(L3) 连续共享度试一下；(L4) 患者独有 pattern 是 tradeoff（好货 vs 污染）；(L5) 用户质疑连续轴能否更好服务下一步连续场发现。实证（`infra/r16/axis_rulers_20260914.json`）：7 根 marker 模块轴（T99/B45/Mye57/ILC12/Epi125/Stromal55/Plasma6 基因）在 **27–30/30 患者**空间相干（Moran 值置换 null，p≤0.01）；等高线分位 0.6/0.7/0.8 三档相邻 Jaccard ≈0.71（旋钮诚实）；端到端阳性对照 **CXCL13 在 25/30 患者的 B 轴等高线内侧/边缘达峰**（教科书位置）。
- 决策：
  1. **Tier-1（地基）**：组成轴尺子厂——轴分=模块基因 z 分均值，跨患者可比性由构造保证（同基因同公式），无匹配步骤；几何 ruler 由等高线切出（离散合法，但分位旋钮披露并扫三档）；轴全部登记 registry，EXPLORATORY_REPRODUCED 起步。
  2. **Tier-2（发现）**：对不上任何轴的 cluster 进普查（平均链接层次聚类杀链式雪崩＋纯度列，纯度 <0.5 强制 LOW_PURITY_MIXED_GROUP 且压回 DESCRIPTIVE）；每组记 best_axis 对齐分（≥0.5 记 axis_explained）；患者独有组升一等公民，等级注明可能污染。
  3. Lane A 双接法：gene~f(轴分)（组成坐标）与 gene~f(到等高线有符号跳数距离)（几何坐标）都允许；L3 连续共享度留作 Tier-2 跨患者测量（R-04 rank-1 方法的复用）。
  4. 已知限制入 LEADS：毛发球会吞弱特异性真 pattern（B 细胞 6 cluster 被吞），候选解 marker-overlap 匹配（L-006）。
- 理由：轴把"跨患者 cluster 匹配"（病态）换成"每患者测同一根轴"（良态）；离散 mask 保留在 Tier-2（无轴可打的区域）合法地位；D-113 连续信号原则落实为"地基连续＋离散几何旋钮诚实"。
- 失效条件：轴集合覆盖不足（新结构全部 UNEXPLAINED 且 Tier-2 全孤证）→扩轴需追加 decision；等高线分位敏感性证伪下游结论→结论随旋钮作废；marker 库本身偏差（D-105 失效条件沿用）。
- 影响：设计 v2.0；`r16/axes.py`＋`scripts/r16_axis_rulers.py`＋census 升级（平均链接/纯度/轴对齐）＋`scripts/r16_build_registry.py`；`infra/r16/field_registry.tsv` v1（228 行：7 Tier-1＋221 Tier-2；20 行 EXPLORATORY_REPRODUCED）；I-021 对 Tier-1 不再阻塞（自相关用值置换），Lane A 关联检验仍需空间保持 null。

### D-116 | 2026-09-15 | Tier-1 轴在 ST-CRC 与 USZ 上独立复现，不混发现集

- 背景：用户看到 registry 全是 Vanderbilt CRC，担心单来源/单平台/单疾病偏见，询问是否应把多来源混入发现。判断：泛用性要，发现阶段混合不要（批次不可分、烧掉验证集、冲掉癌种特异信号）。正确梯子是同一套 marker 模块在每个谱系上独立打分。
- 决策：ST-CRC（同病 7 患者）与 USZ（肾+肺 8 患者、TLS 真值）只做独立轴复现；不用 Vanderbilt HVG-10k 过滤（用各切片实际存在的 voted 基因）；不做跨谱系联合聚类/匹配；不写入发现 registry 行，只写 replication 产物。发现集仍仅 Vanderbilt（D-114 C 维持）。
- 失效条件：若独立复现因基因符号对不齐或坐标读入错误导致系统性零相干，先修 IO 再解释生物学；若用户明确批准「混合发现」才允许改 D-114 C。
- 影响：`scripts/r16_axis_replication.py`；`infra/r16/axis_replication_{ST_CRC_CMS,TLS_VISIUM_USZ}_20260915.json` 与 comparison JSON。

### D-117 | 2026-09-15 | 退役 ILC 轴；工作集定为 T/B/Mye/Epi/Stromal/Plasma 六轴

- 背景：用户指出 ILC 在组织切片（尤其单细胞）里通常稀少、定义不稳。实证上 ILC 也是最弱轴（Visium Moran I 0.19–0.26；D-105 旧 panel 仅 1 基因、viable=False）。「7 轴+x」中唯一额外测试过的 x 是浆细胞模块；内皮/增殖/缺氧等未测。
- 决策：ILC 移出工作集（`r16.axes.RETIRED_AXES`），后续 Lane A/普查轴对齐不再使用。历史产物保留 ILC 数字作档案，registry 行标记 `RETIRED_D117`。工作集六轴：T、B、Mye、Epi、Stromal、Plasma。新 x 轴须单独批准并做与 CXCL13–B 同级的阳性对照，不得默默加入。
- 失效条件：若未来有稳定、跨队列可复现的 ILC 空间程序（且不与 T/NK 模块共线），可追加决策恢复。
- 影响：轴脚本与 census 对齐跳过 ILC；不重跑已完成的三队列实验。

### D-118 | 2026-09-15 | Tier-1 v2 尺子工厂：卷积×tile 轴发现，结论为可信阴性

- 背景：用户要求把 tile 化与卷积化前移到 Tier-1 尺子生产环节（此前六轴全是已知 marker，D-115 L5 遗留"尺子只能量已知"问题）。设计 v2.1：邻域增强表示（λ∈{0,0.2,0.5,0.8}，自表达排除）→ 逐片逐λ PCA top-10 候选轴 → LISA tile 相干（999 次值置换＋BH-FDR）→ 新旧对照（与六轴场 |r|≥0.7 记 KNOWN）→ 跨病人 loading |cos| 平均链接匹配（切 0.25）→ 复发数按有 tile 结构的病人计。
- 决策：新轴升级判据 = 复发 ≥2 病人（须有 tile 结构）＋纯度 ≥0.5；LISA 抽样取 999（200 次的置换地板 0.005 在 BH 下对 256 点测试集已证伪：45 点阳性团全灭，999 次全回收）；半径复核 k=30 只在 pilot＋已升级组所在片执行；旋转混合（非符号）歧义记为已知限制，不静默处理。
- 实测（`infra/r16/axis_factory_20260915.json`，47 片×4λ×10PC=1880 候选）：pilot 证实 λ 即"已知→新奇"旋钮（tile_ok 22→41、known 比率 0.18→0.10）；半径复核中位 0.975 通过；全量 857 组，仅 2 组复现且均为 KNOWN（R16B_001 持家/Epi 程序 29 病人；R16B_002 基质程序 3 病人），**NEW 复现轴为 0**；跨病人 loading 相似度中位仅 0.306，且 0.75→0.5 全程无阈值膝点（降阈只会造毛发团）；1299 个 NEW 且有 tile 结构的候选均匀散在 30 个病人身上——新奇不少，共享为零。
- 理由：这是可信阴性而非方法失效——(1) 无阈值膝点，阴性不依赖 0.75；(2) 与 Tier-2 独立方法收敛到同一结论（可复现结构稀少且薄）；(3) pilot 预注册了旋钮行为，实测吻合。含义：病人特异性主导，未命名共享程序若存在，其 loading 形式在本队列不可见。
- 失效条件：若未来引入旋转感知匹配（子空间/Procrustes）后出现 NEW 复现轴，本结论按方法升级修订（见 LEADS L-008）；若复现引擎（ST-CRC/USZ）上出现本队列没有的共享轴，同样修订。
- 影响：`r16/axis_factory.py`＋`scripts/r16_axis_factory.py`（＋半径复核脚本）＋6 项测试；registry 升 v1.1（1085 行：7 Tier-1＋857 Tier-1b＋221 Tier-2）；设计 v2.1；I-021 对 Tier-1 v2 不再阻塞（自相关用值置换），Lane A 关联检验仍阻塞。

### D-119 | 2026-09-15 | HTAN 内部联合嵌入发现轨：专用方法 + 独立复现

- 背景：用户问为何不按单细胞方法把所有片子整合起来分析。判断：发现集内部联合可以，验证集绝不能混入（D-116 防火墙）；且整合的目标函数就是"抹掉样本间差异"，拿对齐完的数据检验一致性是循环论证——本轨只能当望远镜（发现候选），不能当法官（复现仍须独立检验）。
- 环境勘误：此前"harmonypy、scVI、R 全没有"的调查结论**三条全错**，特此更正——R 一直在（/usr/bin/Rscript，R 4.5.1）；`sc` env 有 harmonypy；`sc_inte` env 有 scvi；R 891 个包中 Seurat 4.4.0/sctransform/harmony/Signac/SpatialExperiment 俱全。教训：环境断言必须以实测命令为准，不凭记忆。
- 决策：联合轨限定 HTAN Vanderbilt 47 片/30 患者；四臂分工——A numpy 手写 ComBat 基线（scanpy import 链在本环境损坏，HVG 缓存+确定性种子已验证行为；跑全量中）、S Seurat v4 SCT anchors（reference 取最大 5 片）、H Harmony（merged log-normalize + RunHarmony）、P PRECAST（待装； intrinsic CAR，本轨唯一空转专用模型）。候选定义与 Tier-2 同口径（Leiden 三档、min 20 点、Moran 200 次）；内部裁判为 split-half 病人对半独立重跑（ComBat 在半内重 fit，匹配用 raw-space 签名）；升级线不变（≥2 病人）；registry 记 tier 2、kind joint_pattern_group。
- ComBat 适用性质疑（用户提出，当面记录）：ComBat 假设各 batch 装的是同一生物学，切片间真实组成差异会被当批次效应吃掉——这是假设层问题，调参修不好。本轨内将其后果转为假阴性而非假claim（split-half 对不上就如实报不复现）；另提议补 Arm 0 无矫正对照以量化 ComBat 足迹，待用户批准。
- 失效条件：任一臂出现系统性单病人纯 cluster（矫正失败）或全混浆糊（过度矫正），该臂结果降级为方法学阴性对照，不进 registry；PRECAST 若装不上或跑不动，P 臂取消，不阻塞 S/H 臂。
- 影响：`r16/joint.py`＋`scripts/r16_joint_embed.py`＋`tests/test_r16_joint.py`；`scripts/r16_export_seurat_bridge.py`＋`scripts/r16_seurat_integration.R`＋`scripts/r16_score_r_labels.py`（待 R 产物）；`scripts/r16_build_registry.py` 接 tier2j/tier2r。

### D-121 | 2026-09-15 | GraphST GPU 臂准备：逐片训练＋全 epoch＋V100 券优先

- 背景：用户批准 GPU 路径并指定 GraphST，要求全 epoch 不打折。事前核实（非记忆）：GraphST 源码无 DGL/torch-geometric 依赖（纯 torch）、聚类默认走 R-mclust（rpy2＋R 4.0.3）、`preprocess` 原生 seurat_v3 top-3000 HVG、`construct_interaction` 对 Visium 建稠密 O(n²) 距离矩阵。
- 决策：逐片训练，不做 merged vertical——113k 点 pooled 的稠密距离矩阵约 100GB，24/48GB 卡放不下；跨片复现由 registry 现有机器（匹配/纯度/Moran/split-half）检验，不强迫对齐（D-119 原则沿用）。原生 seurat_v3 top-3000（输入池为 HVG-10k，披露近似）；random_seed 改 20260914（上游默认 41，披露）；聚类用 Leiden 三档＋min 20（不上 mclust，免 R 桥、与 Tier-2 同口径）；epochs=600 全量；labels 逐片增量落盘（防抢占作废）；torch 2.1.2+cu118（sm_70/V100 有文档支持，不追新版）；V100 6 小时券先烧，不够切 4090（¥1.65/时，先前 ¥15 封顶继续有效）。
- 失效条件：任一片 embedding 出 NaN/全零则该片记失败不入 registry；若逐片 embedding 跨片匹配率系统性低于 Tier-2 同口径基线，记方法学阴性对照；V100 上 torch/DGL 行为异常则停租转 4090，不烧券硬扛。
- 影响：`scripts/r16_export_graphst_bridge.py`＋`scripts/r16_graphst_run.py`＋`scripts/r16_graphst_setup.sh`＋`tests/test_r16_graphst_bridge.py`（3 项本地全绿，torch 路径 stub 覆盖）；scorer 接 arm G（纯加法，缺文件跳过）；bridge 数据 git 忽略可重建。

### D-122 | 2026-09-16 | GraphST 全量不限时跑完即停（用户指令变更计费纪律）
- 背景：dry-run 全绿（GPU 真训、labels 落盘 665 行、列齐）；全量预估 1.5–2.5h。用户明确：不用管 6 小时券，跑完即停。
- 决策：GraphST 全量（47 片×600 epoch）不限 6 小时券窗口，以完成为准；结束后立即退租（用户另有明确"跑完帮我及时退租关闭"指令）。券外部分的计费如实记录，不另设硬封顶——D-121 的 ¥15 封顶仅适用于 4090 fallback，不适用于本次 V100 全量。
- 影响：`graphst_out_full` 远端增量落盘；完成后回收 labels→scorer→registry，落盘推送，然后执行退租。

### D-123 | 2026-09-16 | GraphST 臂标签语义修复：逐片 id 禁止跨片池化

- 背景：GraphST 全量 47 片跑完（113345 点，1 片用 HVG fallback）。评分前发现 scorer 按原设计把同值标签跨片直接池化——对 Seurat/Harmony（联合嵌入，全局 id）正确，对 GraphST（逐片训练，片内 id）是错的，会凭空捏造跨病人候选（16 个全局取值→约 16 个假复现组）。
- 决策：arm G 走 Tier-2 同构匹配（逐片 cluster 为单元，片内 z 签名，平均链接余弦≥0.75，纯度<0.5 强制 DESCRIPTIVE）；匹配逻辑抽成 `r16.census.match_units_to_groups` 并加回归测试（同值异质必分、同质必合、raw id 无关）。另修 numpy 定宽字符串截断组 id 的静默碰撞（`<U2` 下 "g123"→"g1"），赋值前一律转 object。
- 失效条件：若匹配后干净组过少导致 arm G 无信息量，如实报阴性，不降阈捞数。
- 影响：`r16/census.py` 新增 helper、`scripts/r16_score_r_labels.py` 改写 arm G 路径、`tests/test_r16_g_matching.py` 5 项全绿；registry 行待打分完成后合并。

### D-124 | 2026-09-17 | Seurat 臂 S 判死刑、H 落盘：Harmony 83 候选 78 复现入库

- 背景：R 集成脚本连续四次后台失败，逐一修好——(1) future 全局变量 50GB 上限（plan sequential 无用，在 future() 创建时收集，改 maxSize=Inf）；(2) Arm S 注释掉后工具函数 strip_prefix 一起被注释（提到 block 外）；(3) v4 FindClusters 无 cluster.name 参数（v5 API，静默忽略，改手动拷 Idents）；(4) suffix 地雷：as.character(1.0)=="1" 导致 HARM_1 ≠ HARM_10（改位置后缀）。随后 IntegrateData 在 540 万 anchors 上 kNN 内部长度不一致报错——Seurat v4 在 47 片规模的工程上限，非调参可解。
- 决策：Arm S 按 D-119 失效条件降级为方法学阴性对照，不抢救（用户确认）；Arm H（merged log-normalize＋RunHarmony）全量跑通：47 片、108,775 点（Seurat QC min.cells=3/min.features=200 吃掉 4,570 点，覆盖度 0.9597）。scorer 改 QC-aware：QC 掉点挂哨兵位保 pooled 对齐、不进候选不进 Jaccard、coverage 进 artifact（附带修了 H/G 单独跑时 xarm 只有臂内的如实状态）。
- 实测（`infra/r16/joint_embed_H_20260917.json`）：83 候选，78 EXPLORATORY_REPRODUCED（≥2 病人），5 单病人描述。markers 生物学面孔与 G 臂高度一致：持家（PABPC1/NPM1/HSPD1）、杯状（MUC2/FCGBP/TFF3/SPINK4）、浆（IGKC/IGHA1/IGLC3/IGHJ6）、基质/Epi 程序；Moran 中位普遍 0.5–0.75。H 复现率（94%）远高于 G（9%）——预期内：Harmony 联合嵌入全局 id 直接可比，G 是逐片 id 经匹配，口径不同不可比强弱。
- 失效条件：H 臂结论同样只到 exploratory，复现仍须独立检验（ST-CRC/USZ）；若外部引擎上 H 候选系统性对不上，按 D-119 修订。
- 影响：registry 1740 行（7 Tier-1＋857 Tier-1b＋221 Tier-2＋655 Tier-2j：A 99＋B 117＋G 356＋H 83）；`scripts/r16_seurat_integration.R`（S 封存＋H 位置后缀）＋`scripts/r16_score_r_labels.py`（QC 哨兵）已提交（b80fd1e）。

### D-125 | 2026-09-17 | PRECAST 全量跳过：能跑但不跑（用户批准）

- 背景：PRECAST smoke 已验证可跑（2 片/2000 基因/分钟级）；全量是 47 片 11 万点联合 CAR＋K 网格＋MBIC，与 smoke 不是一个量级。A/B/G/H 四臂已收敛到同一结论（可复现的全是组成程序；新 bench NEW 复现为 0，可信阴性 D-118）。
- 决策：P 臂记 `NOT_RUN_LOW_MARGINAL_VALUE`，smoke 产物与脚本留档当方法资产，不删除。与 S 臂死法不同分开记账：S 是技术上跑不动（方法学阴性对照，D-124），P 是判断性跳过。
- 重开触发器（预先写死）：只有当外部复现引擎（ST-CRC/USZ）出现四个臂全漏掉的空间相干程序时，才请出 PRECAST 追——现在跑是锦上添花，那时才是雪中送炭。
- 影响：机器让给 Lane A（空间保持 null＋距离函数筛查）。

### D-126 | 2026-09-17 | Lane A 全线阴性：6 轴×10000 基因无跨病人一致距离函数

- 背景：Lane A 筛 gene~f（到轴等高线有符号跳数），6 工作轴 × 全基因 × 47 切片 × 20 draws 空间保持 null（180°旋转＋小平移，D-125 机器）。升级线冻结：null p≤0.05＋峰一致≥0.5＋峰高≥0.5z。
- 实测（`infra/r16/laneA_{B,Stromal,Epi,T,Mye,Plasma}_20260917.json`，60000 测试）：B 0、Stromal 0、T 0、Mye 0、Plasma 0（p≤0.05 者 0–16 个，chance 期望 500/轴，null 偏保守）；Epi 19 个但全擦线（p 0.03–0.05、agree 0.51–0.57、峰 19/19 全在 -4 区最深处、面孔全是持家/骨架/Epi 自家基因）——同义反复，非距离函数。top hits 跨轴一致模式：峰全钉 -4、曲线全是"区内高区外低"单调台阶，无一例边缘峰/区外峰/单峰涟漪。
- 阳性对照澄清：CXCL13~B 中位曲线形状正确（区内 0.54→区外 -0.2）但正式口径下 p=0.20、agree=0.23——"区内高往外降"成立，"峰在哪"病人间不一致。旧"25/30 峰≤+1"结论是宽门（-4…+1 全算）下的描述，不收回但降级为门内描述，不作跨病人一致峰解读。
- 决策：Lane A 在 6 轴距离尺上判全线阴性（exploratory）。这是实质阴性非技术失败：null 保守、阳性对照形状成立、top 面孔生物学自洽（自家基因台阶）。H-01 未命名场分支按否证逻辑收缩：组成区周围无跨病人一致的连续距离函数。
- 失效条件：若未来引入旋转感知/多尺度距离尺后出现复现涟漪，本结论按方法升级修订；ST-CRC/USZ 上若出现本队列没有的轴周涟漪，同样修订。
- 影响：registry 待补 19 个 Epi 擦线行（标 EXPLORATORY_REPRODUCED＋擦线披露）；Lane B（残余筛查）成为未命名场唯一剩余入口。

### D-127 | 2026-09-18 | Lane B 全量完成：10000 基因残余筛查，191 复现入库（exploratory）

- 背景：Lane B 回归掉 6 轴＋距离壳（OLS，设计列 25–37）后筛残余空间结构（Moran＋值置换 null 50 draws，I-021 细则自相关检验合法）。判据：病人中位 p≤0.05 且相干病人≥2＋残差 Moran 中位≥0.1（修复 smoke 期 p≤0.01 地板不可达 bug）。
- 实测（`infra/r16/laneB_residual_20260917.json`，10000 测）：191 EXPLORATORY_REPRODUCED（残差 I 中位 0.138，最高 SERF1A 0.467；相干病人中位 24/30）。面孔：持家样（核糖体/骨架/代谢）＋自家地盘基因（EPCAM/CEACAM5/CLDN4 上皮、MUC2/FCGBP 杯状、FN1/DCN 基质、FOS/FOSB 应激）——轴回归吃掉了大头但同类基因仍有残余相干，符合"组成内状态/亚型"解释，不作新场解读。
- 决策：191 行以 kind=laneB_residual_structure 入库（D-126 的 19 个 laneA 行同批在库），registry 1950 行（7＋857＋221＋655＋19＋191）。Lane B 关闭；未命名场两条入口（Lane A 距离函数、Lane B 残余）均已执行完毕。
- 失效条件：残余相干≠新场；任何 laneB 行升级为 claim 前须先排除"同类基因亚型/状态"解释并经复现引擎检验。
- 影响：`scripts/r16_laneB_screen.py`＋产物＋registry 行已提交；`scripts/r16_build_registry.py` 计数器补 laneA/laneB 分项。

### D-128 | 2026-09-18 | T1 PC 符号修复完成：分组 857→772，复现仍为 2 组且均为六轴联合可解释

- 背景：恢复包 A02（signed cosine 分组＋raw 平均）。修复限 PC 语义：回归测试（`tests/test_r16_recovery_pc_matching.py`，3 绿，含旧代码必失败项）＋`r16/recovery/pc_matching.py`（|cos| 分组＋medoid 定向平均＋conflict 上报）＋`scripts/r16_recovery_T1_rematch.py` 从 65M loadings 直接重算。
- 实测（`infra/r16/recovery_20260918/audit/pc_matching_reanalysis.json`）：857→772 组；8 旧组分裂、91 新组合并（均为两旧组合一，符号翻转对）；复现组仍为 2（R16B_001 持家/Epi 29 病人 184 成员；R16B_002 基质 3 病人 9 成员），conflict 均为 0.0。联合解释力（A03 正确新奇性）：两组 joint_axis_R²=0.02/0.08——不是"单轴 max|r|<0.7 的 NEW"，而是旧口径误标：它们本就 KNOWN（Epi/Stromal），联合口径下 `new_reproduced=0` 维持。
- 子空间诊断（L-008 有边界检查）：R16B_002 三成员片内 top 匹配集中单 PC（cos 0.78–0.95），无跨 PC 拆分证据——逐向量匹配对此组够用，不升级旋转方法。
- 几何诊断（A10）：q-Jaccard 降级为 `nested_quantile_overlap_descriptive`；q0.7 等高线保留定位功能，仅"稳定性"解读作废。
- 决策：D-118 的 NEW=0 在符号修复＋联合口径下维持（772 组，0 NEW 复现）。I-023 可 RESOLVED（遗留旋转限制转 LEADS L-008）。产物路径见上；历史 JSON 不改写。
- 失效条件：若未来对未解释候选的子空间诊断显示跨 PC 拆分，按 L-008 升级方法后重判。

### D-129 | 2026-09-18 | T2 机械修复验证通过：患者先汇总＋同算子闭环；remap 近似性质首次实测

- 背景：恢复包 A05/A06/A09。`r16/recovery/laneA_stats.py`（唯一 patient-first＋同算子）＋回归测试 3 绿＋诊断集（6 患者×先一片×10 基因×199 联合 draws，`scripts/r16_recovery_T2_laneA_small.py`）。
- 实测（`infra/r16/recovery_20260918/laneA/`）：SYNTH_NEG 双 p>0.1（机械无假阳性）；SYNTH_BUMP 固定对比检出（tcon+1.53/p=0.005），峰高统计量对宽 bump 不敏感（p=0.215，符合预期）；CXCL13~B tcon+1.08（区内高成立，不答峰位）。
- remap 实测（1194 draws）：接受率 1.00（诊断片几何规则）；frac 中位 0.71–0.84；唯一源比例 0.72–0.84（多对一量化约 20–28% 复用）；Moran 保持比 0.63–8.96，小 Moran 片爆炸——remap 不保变异函数尺度，正式降为诊断对照，主证据等变异函数匹配 surrogate。
- 决策：A05 机械部分关闭；全部占位 p 记 `NOT_CALIBRATED`；I-022 保持 OPEN（待主 surrogate 落地）。历史 D-126 数字不动。
- 失效条件：主 surrogate 校准失败则关联 p 永久记 NOT_CALIBRATED，转 T3/T4 独立预测为主证据（02 T2 完成条件已预告）。

### D-130 | 2026-09-18 | T3 六程序冻结：191＋联合臂→556 union→共表达→6 程序，全员带四标签

- 背景：02 T3＋04-B 分层（身份/状态/技术/空间信息分开评价；D-127"排除亚型/状态"改为分层）。
- 实测（`infra/r16/recovery_20260918/programs/`）：union 556（191＋365；191 中 140 在 A/B、48 在 Tier-2 十三组——同批东西，多方法不计独立证据）；共表达 1 巨模块 M213（159，|r|中位 0.59，回归六轴后残余 0.18——搭便车为主但残余非零）＋Ig M93＋基质 M132＋平滑肌/ MHC2/应激小模块；368 单例无协同。
- 六程序：P-epi（16，QC残余0.34，30/30病人）/P-stromal（9，0.22，非线性quad砍半如实记）/P-plasma（11，与Plasma轴loo 0.815 不可分，降级只验分子再现）/P-smmhc（3，小而干净）/P-mhc2（2，唯一非丰度候选，残余最弱0.14）/P-stress（2，FOS系，技术梯度未排除双标签）。MT模块按技术伪影排除；HIST1H簇归增殖亚态不单列。
- 决策：6 程序定义冻结结构已写（`program_definitions.json`，frozen时间戳空，T4前填）；readout基因T4划分并机器检查分离；患者效应全表180行无截断（A18关闭）。I-024待T4后降级。
- 失效条件：T4 无增量则程序收缩为描述性；P-plasma 若分子亦不复现则关闭。

### D-131 | 2026-09-18 | T4 独立复现＋增量：M1 全员成立，M2 仅 mhc2/stress 成立

- 背景：02 T4＋03 §6（分子/空间/增量三层分开；input/readout分离；留一患者；局部校准口径）。定义冻结 2026-09-18T16:33:33Z（hash 675df82fde39b306）；readout 三队列全覆盖已验；P-plasma 只验分子（NOT_SEPARABLE）。
- 实测（`infra/r16/recovery_20260918/{external_replication.tsv,incremental_prediction.tsv}`）：
  - 分子再现：ST-CRC 全员强（corr 0.49–0.89）；USZ 分化——stromal 0.81/mhc2 0.72/epi 0.58 强，plasma 0.66 中，stress 0.33 弱（含 0.08–0.23 失败病人）。
  - 空间再现：input/readout Moran 双双>0.2（除 USZ stress-input 0.23 及个别肺片），但"人人Moran阳性≠结构相同"（03 §6.1 已警告），只作必要条件。
  - 增量 dM1（M1-M0，程序活性超组成）：10/10 队列-程序中位为正，其中 8 格 bootstrap CI 不含零；mhc2 两队列 CI 跨零（含 5399a83c1999 -80% 反例：M0 已很好 0.08，input 外推翻车）。
  - 增量 dM2（M2-M1，邻域空间项）：仅 mhc2（9.4/8.2%，CI 不含零）与 stress/ST-CRC（15.8%）成立；epi/stromal/smmhc 中位≈0（epi/ST-CRC -0.1%）——活性有信息、空间排布无额外信息，诚实记为场活性复现而非空间组织发现。
- 决策：T4 完成。P-epi/stromal/smmhc：组成外活性信息成立（H-01/H-03 支持），空间组织不成立；P-mhc2：三层全成立但最弱、CIs 宽，保留为候选；P-stress：ST-CRC 强、USZ 分化（技术梯度标签未摘）；P-plasma：分子再现成立（ST-CRC 0.60/USZ 0.66），增量未测。I-024 降级为方法备注（残差 p 仅筛查用，主证据为独立预测）。
- 失效条件：T5 遮蔽若显示程序分可被最近邻插值替代，M2 结论再收缩；USZ 肺片系统性失败若经组织学确认为缺靶组织，则不记为生物学反例。

### D-132 | 2026-09-18 | T5 遮蔽关门：场预测 m1 碾压 NN，m2 归零；mhc2 超 B 轴成 TLS 最强线索

- 背景：02 T5＋03 §7（窗与 GT 无关；隐藏派生不可见；泄漏探针；ST-CRC IC aggregate 不进 B，D-033）。
- 实测 A（场预测，396 行）：m1 对 NN 改善 45–86%（全程序×两队列）——程序分携带场信息，与 T4 dM1 一致；m2 对 m1 改善≈0（-1.8%…+5.2%）——遮蔽尺度下邻域项无增量。m2 全 NaN bug（隐藏点特征被填 NaN）已修并重跑，泄漏探针 22/22 通过。
- 实测 B（USZ-TLS 隐藏核心，56 行）：阳性窗仅 0–3/片（窗 8×8 对 30–344 TLS 点太粗，KC1/KC2 零阳性窗）；6 可评片中 mhc2 5 片 ≥B 轴（LC3 0.67>0.54、LC2/LC1 1.0 持平、LC5 0.57>0.43、LC4 0.56>0.50），stress 全弱于 B 轴。mhc2 与 T4-M2 独立双命中，但窗粒度只给 0.1–0.2 精度，记"最强线索"不记结论。
- 决策：T5 完成。主张上限：场预测成立（A）；结构定位仅 mhc2 线索（B）。H-02 未触及（无位置后验）；H-01/H-03 由 T4 承担，T5 只补"场信息可遮蔽预测"的独立证据。

### D-133 | 2026-09-19 | T6 恢复包关门：plan 澄清为说明不修文；旧否证外推维持 METHOD_LIMITED

- 背景：恢复包 04 要求在 P-01/H-03/H-05 处澄清身份/活动/空间信息分层。核对 plan 原文：H-01 放弃条件（两队列 CI 含零＋不优于基线）与 T4 证据结构一致；H-03 放弃条件（组成-only 等价＋高分辨无固定类型后状态变化）与分层评价一致；H-05 条件性偏移与 USZ 分化结果一致。plan 写的是"什么算数"，恢复包做的是"这次算出来什么"——无冲突，故澄清只做说明，不修改 plan 正文（AGENTS.md：plan 改动需科学判断改变）。
- 实际执行改变：D-128（T1）/D-129（T2-small）/D-130（T3）/D-131（T4）/D-132（T5）；I-023 关闭，I-024 降 WONT_FIX；I-022/I-025 保持 OPEN（T2P2/深度确证未做）；旧产物路径见 old_vs_new_conclusions.tsv；历史决策原文未改写。
- 失效条件：若 T2P2 主 surrogate 落地或深度确证完成，按新证据追加 D-13X 并更新 I-022/I-025，不追改本条。

### D-134 | 2026-09-19 | L-009 mhc2 TLS验证失败：细窗963配对下B轴显著更好，退回场活性描述

- 背景：T5B"mhc2超B轴5/6片"基于1–3阳性窗/片（粒度0.1–0.2）。L-009改TLS-foci自适应细窗（3/5/8 hops titration）+3 seeds，同窗同seed配对head-to-head（`scripts/r16_recovery_L009_fine.py`，963行）。
- 实测（`infra/r16/recovery_20260918/{masked_structure_fine.tsv,L009_verdict.md}`）：mhc2 hit 0.606 vs B轴 0.672，差-0.065；不一致对34 vs 97，McNemar精确双侧p=3.3e-08。三档窗全输；8片7片输或持平；LC5小focus窗覆盖不足属设计局限。
- 决策：T5B的"最强线索"修订为不成立；mhc2退回场活性程序（T4的M1/M2证据不受影响）；H-02仍未测试。L-009关闭。
- 失效条件：若未来有独立TLS队列+更密标注再验出mhc2稳定超B轴，可重开；本轮不再投入。

### D-135 | 2026-09-20 | GO线开工：GOBP全量7583通路当先验程序，T3线归档不删

- 背景：用户要求用全部 GOBP 重跑一遍，不自己发现 pattern（T3 共表达 556→6 线保留归档）。方向切换：数据驱动发现→先验通路普查；H-03（复合生态）下 GOBP 是现成词典，P-07 B-09 允许无锚点软生态场。
- 决策：GMT=c5.go.v2025.1.Hs.symbols.gmt（sha256 见 build_provenance.json）；coverage=三队列交集（cache10k∩STCRC∩USZ=9012基因），<3基因→NOT_TESTABLE（用户拍板，596条）；half-split seed 20260920 shuffle+ceil-input（input≥2，零重叠机器检查）；covered集相同的117个重复ID只跑一份（id_map保留，原始全留，家族归并只做report层）。
- 范围：testable 6987 GO ID → unique 6870 份；验证机器逐字复用 T4（发现集QC三层Moran＋外部三层＋留一患者M0/M1/M2）；smoke20 已通（119s，rep/inc各300行）。
- 失效条件：dM1/dM2 口径与 T4 完全一致；多重检验按有效独立通路计，不按6987计。

### D-136 | 2026-09-20 | GO线全量完成：6870通路×三层验证，820双队列dM1>20%，家族归并只做report层

- 背景：用户要求全部 GOBP 重跑，不自己发现 pattern（D-135 开工）。GMT v2025.1，coverage=三队列交集9012基因，<3基因596条NOT_TESTABLE，half-split seed 20260920，117重复ID跑一份。
- 实测（`infra/r16/recovery_20260918/programs_go/go_summary_full.tsv`，6987行含重复展开）：
  发现集QC残余中位0.040（p90=0.102，722通路>0.1）；ST-CRC dM1中位1.2%（12.9%>20%）、USZ dM1中位15.1%（44.0%>20%）；双队列dM1>20%共820通路。
  家族归并（GO ID前两段词根，2727家族）：top为有丝分裂/细胞周期/剪接/代谢家族（两队列均值40–76%），bottom为神经/昼夜/转运家族（≈0%）。
  T3旧6程序对照：上皮分化/抗原呈递/肌肉家族在GO线同样强（USZ 55–91%），应激家族分化（ST≈0/USZ 34%），与T4结论同向。
- 决策：GO线与T3线并列，不互相覆盖；多重检验按有效独立通路计（report层），原始6870份全留；NOT_TESTABLE 596条诚实标记。
- 失效条件：家族归并词根粗糙，大通路家族内部异质性高时以单通路行为准，不以家族均值杀单通路。

### D-137 | 2026-09-20 | Spline 距离曲线：mhc2 翻盘失败，但发现 peri-TLS 晕环形态

- 背景：M2 升级候选1（抄灵感论文 Science 2026）：程序分 ns(df=3)~到TLS距离，
  USZ 8片TLS GT锚点，TUM spots，距离rescale 0–1，形态=Spearman单调+15%变化地板。
- 实测（`infra/r16/recovery_20260918/spline_{summary,curves}_full.tsv`，64 fits）：
  形态分布 flat 32 / nonmonotonic 22 / low-to-high 6 / high-to-low 4——干净单调极少。
  mhc2：仅KC1 high-to-low（r²=0.07），LC3/LC4 nonmonotonic（r²=0.20/0.22），其余flat；
  距离梯度翻盘失败，L-009结论维持（mhc2退回场活性）。
- 新形态（非预设发现）：LC3/LC4的mhc2曲线呈"近TLS峰→谷→远端回升"晕环形
  （LC3: 1.70→0.66→1.12；LC4: 1.16→0.59→0.55），B轴同片曲线平坦（0.43→0.26→0.30）。
  mhc2振幅（中位0.34）数倍于B轴（0.07）——mhc2在TLS周围有B轴没有的空间结构，
  只是不是单调梯度，而是晕环。这是M2线性均值和L-009 hit/miss都看不见的形态。
- 决策：spline纳入M2形态库（与线性M2并列，不替代）；晕环形态记LEADS L-011
  （需多片复现+机制解释，当前仅2片描述性）；mhc2"距离梯度翻盘"关闭。
- 失效条件：若晕环在ST-CRC（无TLS GT，需B轴代理锚点）或更多USZ片复现失败，
  退回单片描述，不升级。

### D-138 | 2026-09-20 | Spline×820通路：形态由切片主导，通路身份被抹平；L-011晕环降级

- 背景：用户要求把 spline 铺到 GO 线 820 个双队列 dM1>20% 通路（USZ 8 片 TLS 锚点）。口径同 D-137。
- 实测（`spline_summary_go820_test5.tsv` 6560 fits）：形态几乎完全由切片决定（LC3 100% nonmonotonic，KC3/LC2 99.9% flat，KC1 71% high-to-low）；通路多数票 nonmonotonic 725 / flat 95 / 单调 0；没有任何通路在 ≥3 片上稳定单调。LC3 随机通路对曲线 Spearman ρ≈1.0。
- 原因：spline 拟合未残差 log1p 均值，深度/组成沿 TLS 距离的共变被所有通路共享。L-011 peri-TLS 晕环修订为切片几何，非 mhc2 特异。
- 决策：820 通路 spline 作为方法阴性（切片主导）归档；不把形态当通路身份。继续 spline 必须先 Q+C 残差。mhc2 梯度翻盘（D-137）维持。
- 失效条件：QC 残差后通路形状去共线且 ≥3 片同向，可重开通路特异形态主张。

### D-139 | 2026-09-20 | Spline三路线：残差拆共线，秩分数不能替代；残差后无迁移梯度

- 背景：用户要求秩分数（Wilcoxon AUC，AUCell同类最快）、仅Q+C残差、秩+残差三条路线分开做，铺到820双队列通路。
- 实测（每路6560 fits）：raw LC3 ρ≈1.00；rank LC3 ρ=0.90且因AUC振幅小被15%门槛打成全flat；resid LC3 ρ=0.10、r²中位0.012；rank+resid ρ≈0、r²中位0.008。残差后≥3片h2l 47/l2h 42，但r²~0.02，抗原呈递家族2/2/3混合。C相对B无增量。
- 决策：形态主张关闭。残差是拆切片共线的必要操作；秩分数是稳健赋值、不能替代残差。47条弱同向不升级。M2线性均值和spline距离曲线都不再作为通路身份证据。C6轴名必须用 `Mye` 不是 `Myeloid`（残差0-fit bug）。
- 失效条件：换独立锚点（非TLS）或残差后r²跨片中位>0.10且形态同向，可重开。

### D-140 | 2026-09-21 | 审计修复与重算口径冻结

- 授权：用户要求修复审计问题并报告最终结论。范围为 F01–F11；不启动新走廊/几何字典，不重训 R-04，不恢复 GPU 或外部数据。
- 输入与产物：历史结果保留；新输出写 `infra/repair_20260921/`。六程序沿用已有 input/readout；GO 沿用 6870 代表的划分，将别名定义统一到代表并钉死 hash。本地实际共同基因表用于完整性检查；没有覆盖或没有变化均显式 NOT_TESTABLE。
- T4/GO：主评分 mean(log1p counts)，组成轴同尺度，不逐片 z-score；每程序将 input+readout 从组成轴、QC 深度/检出数的基因集合中排除。任何轴无剩余基因，不填零冒充组成控制，增量记 NOT_SEPARABLE。训练/误差先切片后患者等权，留一患者；M0=QC+六轴、M1=+input、M2=+可见 input 的邻域项。六程序另做排除集合外深度归一化敏感性及训练患者内选择的单 input marker 基线；GO 主流程全量重算，计数仅描述，不声称独立发现或 H-03 成立。深度归一化敏感性限六程序，不暗示已覆盖全 GO。
- T5：隐藏窗口内全部分子值不可见；QC/组成/程序及读出训练只使用窗外。隐藏位置特征来自可见邻居，完整预测入口必须通过隐藏值干预。所有模型共享可见分子预算，与最近邻和邻域平滑读出基线比较；这是局部校准的场预测，不是零样本外推。
- 结构窗：unknown 保持 NaN；阳性窗要求已知覆盖≥80%且已知区域 TLS 占比≥10%；阴性窗要求全部点有明确非TLS标签。否则 UNKNOWN/PARTIAL，不参与二分类。细窗正例仍按已知 TLS 焦点居中，只能称靶向邻域对照；不得称无 GT 的盲定位。程序/B轴分数都只从窗外重建。按焦点/窗口先汇总到患者，患者级区间；不把 seed/窗宽增加当生物学 N。
- 几何：10x array 坐标明确转为等距六边形 pitch 单位；焦点使用注册时的六边邻接，不跨缺失点用 kNN 接桥。TLS 距离采用到已标注 TLS spot 集合的欧氏距离；核心负距离到集合外最近 spot，不叫连续组织边界距离。
- spline：真正施加边界二阶导为零，knots 只来自拟合数据；形态用范围内曲线变化/观测分数 SD，门槛0.2 SD，导数≥95%同号才记单调，去掉拟合曲线箱值显著性。原820通路及修复后新阈值命中的并集，全8片、raw/rank/resid/rank+resid重算；六程序也重算。残差与T4使用同一排除规则。普通F检验不再输出正式p；空间关联p继续 NOT_CALIBRATED。分块留出距离增量与注入场检出用于描述能力，不能把无检出等同于否证。
- T1/R04：T1用实际spot分数对组成分数的回归重判两复现组，loading投影只保留为几何诊断；R04分别记录共享场/特异信息，缺少共同残余检验时为未判定，不自动触发后续节点或K搜索。
- 失效/停止条件：输入 hash/标签对齐失败即停相应分支；不能满足已知标签、组成覆盖或数值可识别性则明确不可测，不放宽门槛追求结果。bootstrap区间为患者层、条件于固定CV预测的探索性区间；没有空间null和新独立队列，不升级确认性空间场。

## D-141｜修复期间补充的基因覆盖与缓存缺失约束（2026-09-21）

- 来源：对冻结定义与源 feature/var metadata 逐片核对。USZ 的 epi 缺 CEACAM6/KRT19、smmhc 缺 MYH11；stress 的 DEPRECATED 标识在外部源不存在；plasma 也缺多项。不能把静默删基因后的评分当完整冻结程序。ST-CRC 的 epi/smmhc 源文件覆盖完整，不等于这批数据不可用，但它们不满足本轮跨队列共同面板的完整程序比较。
- 选择：主重算保持冻结 input/readout，不猜测 DEPRECATED 标识、不改成部分基因代理。六程序的匹配面板增量只对完整覆盖者计算，其他项显式 NOT_TESTABLE；源数据覆盖与当前任务资格分表记录。不新增数据，也不把少基因评分冒充原程序。
- 另发现：发现集 cache 的 10k 列是并集，部分源文件没有的基因被补零。9012 分析面板是两外部队列实际共同基因与该缓存面板的交集，不代表每个发现切片的源 var 都完整。现在另外读取 47 份源 var：发现集中每程序逐片检查完整性，缺基因记 NOT_TESTABLE；组成均值排除源文件未测基因，不能把缺失当实测零。外部完整读出不会仅因某张发现片缺基因而被取消资格。
- 失效条件：若源 symbol 与 stable ID 存在明确的可核验 crosswalk，另行冻结映射后重算；当前 metadata 未提供可直接跨队列采用的 DEPRECATED 映射。没有完整 readout、足够患者或可分组成时保留不可测，不以零填补或阈值调整恢复数量。
- 影响：这是覆盖与资格修正，不改变 plan.md 的科学假设。新的完整 GO 尝试仍覆盖全部 6870 代表，但成功计算数、缺失数与常数目标数分别报告。

## D-142｜修复结果采用边界与收尾（2026-09-21）

- 依据：`docs/IMPLEMENTATION_REPAIR_20260921.md` 与 `infra/repair_20260921/validation_receipt.json`。33 项测试、45 个显式加权最小二乘误差复算、44 次真实完整预测干预通过；全量 GO 与曲线重算完成，缺失/常数目标单列。
- 选择：采用修正结果，撤回旧遮蔽优势、窗口群体显著性及普遍残余场阴性。保留基质的分子预测价值，不能升级为隐藏结构定位。两个复现 PC 组的实际 spot 组成解释 R² 较高，仅支持组成相容诊断。
- 标注与基因不足分别触发 I-027/I-028 的局部硬阻塞；不把 NOT_TESTABLE 写成阴性，不为恢复原样本量放宽已冻结条件。R04 的共同残余场仍未测试，科学状态为可识别性阻塞；操作停止、K 搜索关闭保持。
- 后续边界：本次不启动 T2P2、走廊/字典、旧模型重训或新数据获取。先依据完整基因与独立标注确定能回答的问题，再安排后续计算。未修改科学假设本身，plan.md 保持不变。
- 失效条件：若数据自带 metadata 支持新的完整映射、独立标签覆盖改变，或合格留出试验出现超过最强基线的稳定增量，应重新评估；当前分子共变、单片距离关联和无噪声线性注入不足以跨越这些条件。

## D-143｜2026-09-22｜修复后资格普查与有限方向性描述
- 授权：用户要求继续工作。先完成 STATUS 指定的本地数据资格核查；不恢复旧模型训练、K 搜索、GPU、外部数据或确认性定位主张。
- 范围：检查 69 片修复输入的完整程序覆盖；注册表用于核对独立标注来源。逐 spot 重放限当前 22 片外部比较数据，不能把修复缓存中的空标签解释为原始数据没有标注。USZ 按源 h5ad 标签验证；ST-CRC 仅清点原始病理类别，不从类别名称猜造 TLS 标签。
- 在看分子方向结果前冻结：仅 P-stromal/P-mhc2，USZ 八片；TLS 六边连通域质心指向最近的明确 TUM spot。距离 4–24 pitch；半径为该距离的 0.25–1，分三等宽径向带；朝向和左右 90 度对照各取半角 30 度。仅明确非 TLS 的已标注 spot，每带三扇区各至少 3 点、每焦点至少 2 带。每带先对照左右均值，再跨带/焦点/切片/患者等权。资格只看坐标与标签，不看表达；不根据结果调宽扇区或补 unknown。
- 如存在合格焦点，计算 readout 对 QC+六轴的片内残差（排除完整 input/readout；raw 与深度归一化两路线），以片内残差 SD 为单位报告定向差异。片内回归是描述性组成调整，不是留出预测；GT 只定义分析几何。完整空间扫描、抗原呈递 GO 家族扩展与正式空间 p 均不在本次计算中。
- 排除方案：不以旋转 90 度或随机打乱直接冒充保留空间自相关的零模型；不把扇区差异当作隐藏定位，也不把无差异当作没有空间场。该描述不替代 T2P2、患者外推和邻域平滑竞争。
- 失效条件：源标注与缓存不一致、hash 失配、空组成轴或几何资格不足，则停止对应单元并记录不可测。少于 3 位合格患者不作群体汇总区间；即使满足也只报告描述性效应，不输出关联 p。新结果仅用于评估 H-01/H-05 后续可测性，不升级假设支持。

## D-144｜2026-09-22｜有限方向描述结果采用与停止扩展
- 依据：D-143 冻结规则及 infra/directional_eligibility_20260922/receipt.json；108焦点中11个/3患者合格，两个完整程序四条评分组合均完成。
- 选择：只采用患者层描述；不以原始抗原呈递3/3同向作为复现，归一化后该方向不稳定。基质2/3同向不升级稳定走廊。当前不扩到GO家族，不调阈值增加N；下一步优先空间误报率/检出能力校准。完整空间场字典未执行。
- 补充解释：源标注缺失条形码按未知处理，与修复缓存一致；ST-CRC已有病理GT不能因修复缓存未加载而丢弃。当前片内组成调整不是留出验证，径向带对照不保证消除带内距离/区室/形态差异。
- 失效条件：预先设定的校准显示方法有可用误报率与检出能力，且合格独立患者/队列支持稳定增量后，方可重新评估后续扩展。无差异不能替代对所有空间场的否证；GPU、外部数据、K搜索和旧模型重训边界不变。

## D-145｜2026-09-22｜全 GO BP 阶段规则与结果采用
- 授权与记录时点：用户要求完成GO BP测试，并要求阶段完成后统一落报告。因此执行前固定 infra/gobp_stage_20260922/run_contract.json 与 calibration_contract.json，本条在阶段结束统一登记；期间未改动基因划分或筛选阈值。
- 范围：6870冻结代表、22片/15患者、raw与背景深度归一化；分队列LOPO、12个固定遮蔽窗、D-143方向几何、全量USZ组成调整分块距离预测。input/readout继续从组成与QC排除；不可测目标保留状态。原始流程重新全量执行以核对加速计算。
- 选择：采用稀疏逐点评分和分批回归，完整预测仅接触可见行；拒绝先拟合全数据再减隐藏矩的数值泄漏路径。132次真实隐藏值干预和直接计算对照通过。原103050条分子误差与7208条距离增量复算一致。
- 推断：方向的3患者双侧精确检验最小p=0.25。患者级符号翻转仅在患者独立与零假设符号对称下解释；三个完整端点家族分别Holm校正。已知高斯空间模型上的120条件校准不替代真实表达协方差/变异函数匹配，T2P2及真实空间关联p仍未完成。
- 结果采用：分子双队列>20%筛选629→45，两种评分共同45。遮蔽M2胜平滑共同2497，但同时胜组成/最近邻/平滑、两队列两评分均为正者仅92，其最佳最弱条件增量0.3671%，与45交集为0。计数为描述性，不是独立发现；接近数值零的正号不作生物增益。空间方向与距离保留局部描述，没有升级隐藏定位、完整H-03或零样本H-05。
- 未运行与失效条件：无GPU、新数据、旧模型重训、完整字典或全GO单marker验证。若独立合格患者/标签或预先限定端点的校准支持实质增量，可重审定位路线；不得通过扩大基因集数量替代患者数，不把未检出当作普遍空间场阴性。plan.md保持不变。

## D-146｜2026-09-24｜GO BP候选重叠与GSE175540验证层级边界
- 背景：全GO阶段筛出的45项分子预测候选中，input/readout并集有42/990对的 Jaccard≥0.50；本地R-02审计显示GSE175540的18条TLS来源都映射到冻结外部角色和患者，但没有块编号。
- 选择：保留45项作为可追溯的预测候选分别列示，同时报告成对重叠；不把它们当作45条独立机制，也不在事后按任意重叠阈值合并术语。GSE175540只按R-02冻结的患者级外部角色使用，不作块级独立重复主张。
- 理由：尚无预注册的GO术语折叠规则，Jaccard清点不能决定生物机制独立性；空白 `block_id` 不支持块级推断，但不撤销冻结患者级验证角色。
- 复查条件：若事先冻结了GO折叠与增量判据，或R-02获得可核验组织块映射，重新评估解释层级。当前H-01/H-03/H-05状态不变；单基因预测基线与GSE175540样本/基因资格核对仍未运行。

## D-147｜2026-09-24｜45个GO BP候选的单基因基线采用嵌套患者留一
- 背景：45个候选是用两队列、两种评分的既有结果筛出，input 集合为50–512个基因；旧双程序对照只在外层训练患者内按训练误差选marker，没有内层留出选择。
- 选择：对每个候选，只在冻结 input 基因中选单基因；沿用 M0（QC+六组成轴）、M1（完整 input 程序）及同一患者级外层 LOPO 和 raw/深度归一化评分。单基因选择在每个外层训练集内再做患者级 LOPO，以内层留出 MSE 均值最低者为准；随后在整个外层训练集拟合，并评估未见外层患者。并列时按基因符号排序决胜。
- 理由：不让外层患者参与marker选择，并避免以训练内误差直接挑单基因而低估择优偏差。45个候选本身仍是用全队列结果后选，本轮对照只作条件性探索，不作独立确认性证据或多重检验结论。
- 复查条件：若将候选选择也嵌入外层训练或取得独立冻结队列，重做确认性比较。当前不改变H-01/H-03/H-05结论；完整空间推断仍需另行验证。

## D-148｜2026-09-27｜GSE175540外部评分：冻结映射规则、端点与两臂并行
- 背景：roadmap要求先冻结重复feature符号处理与患者级TLS端点，再定16×31还是聚合后45。用户指示两条都跑，扩展路线的两聚合规则都试。
- 选择：Arm1严格（16严格患者×31无歧义程序）；Arm2同16患者×45程序，两路线：AGG-SUM（重复符号跨feature行先加总计数再打分）与AGG-DROP（剔除歧义符号，<2有效基因记NOT_TESTABLE）。lib_excl两臂一致，只让分子不同。评分raw主+深度归一化敏（GO惯例，D-140）；端点主SMD+敏AUC；跨患者中位+符号计数+双侧符号p描述性，每臂每评分SMD家族内Holm。患者为单位，无块级、无跨癌种主张。冻结先于任何表达分数，见 infra/gobp_external_score_20260927/run_contract.json。
- 理由：SUM把切分注释当同一基因合拢，DROP不作同一性假设；两者分叉即报告对歧义符号的敏感性，不事后二选一。
- 复查条件：若来源metadata给出stable-ID級crosswalk，另冻映射重算；当前D-141约束不变。

## D-150｜2026-09-28｜T2P2变异函数匹配不可用，Lane A关联p永久未校准
- 依据：D-149 门槛与 `infra/r16/t2p2_20260928/summary.json`。六张诊断片均为六边形 pitch。一张远处点只有2个，统计量不可算，校准用其余5张（实现500次、错规格各400次、功效200次）。
- 结果：正确高斯场误报6.0%，在[0.02,0.10]内；图平滑20.5%、零膨胀平滑18.3%，出窗。近区1 SD注入检出90.5%。CXCL13 5/6、COL1A1 6/6 片通过拟合门槛。判定 `PERMANENTLY_NOT_CALIBRATED`，不是 `BLOCKED_IMPLEMENTATION`。
- 选择：Lane A 距离关联 p 永久不可用于发现或定位。不放宽门槛，不改用旋转或值置换补救。已知模型校准不替代。旧 D-126 数字保留，外推仍为 `METHOD_LIMITED`。
- 对假设：H-01 无新支持或否证。没有可用 p 不等于没有空间场。H-02 仍未测。`docs/plan.md` 不变。I-022 改为 `SEALED_NOT_CALIBRATED`，项目不关闭。
- 失效条件：只有新的、事先冻结的非平稳或非高斯零模型在同类错规格场上把误报打回 2%–10%，才能重开关联 p。不能用这次失败的 p 值筛选基因后再称复现。

## D-149｜2026-09-28｜T2P2变异函数匹配可用性门槛先于结果冻结
- 授权：用户要求空间null的两条路线都试：跑T2P2主surrogate，并在不可用时永久标NOT_CALIBRATED。
- 选择：只试一个主方法（平稳指数+块金，经验变异函数，Cholesky surrogate）。几何用T2已选的6个效应盲诊断片。坐标优先源obs的array_row/col转六边形pitch；对不齐则显式回退像素坐标，不假装已转换。统计量是B轮廓内（hop≤0）减远处（hop≥4）。
- 可用当且仅当：正确模型实现FPR落在[0.02,0.10]；图平滑与零膨胀平滑两个错规格家族的FPR也落在该窗；1 SD近区注入检出率≥0.40；CXCL13或COL1A1至少一个基因在≥5/6片上拟合相对RMSE≤0.35且surrogate变异函数相对RMSE≤0.35、与原场绝对相关≤0.20。否则Lane A关联p记PERMANENTLY_NOT_CALIBRATED。实现FPR出窗则记BLOCKED_IMPLEMENTATION，不把代码错误写成几何否证。
- 排除：不在看到结果后放宽门槛；已知模型校准不替代本门；本轮不输出生物学发现p，不扩GO。门槛文件 infra/r16/t2p2_20260928/run_contract.json。
- 失效条件：若后续换锚点或非平稳协方差模型，须新决策，不能改本条数字。

## D-151｜2026-09-28｜GSE175540患者级评分采用：主终点不显著，深度归一化敏感性不升级
- 依据：D-148 合同与 `infra/gobp_external_score_20260927/receipt.json`。3872 行患者端点，收据哈希复核通过。
- 主终点：严格臂 raw，16 患者 × 31 项。22/31 中位差为正，但中位差仅 0.045，中位 AUC 0.510。无一项未校正符号 p<0.1，Holm 后为 0。适应性免疫 11/16 为正，含一位负向患者。
- 敏感性：深度归一化下适应性免疫 16/16 为正，中位差 1.37，三臂 Holm 后仍 <0.05。加总与丢弃相对严格臂无实质差异。该命中不改变主终点。
- 选择：不把敏感性命中升为主结论，不称跨癌种或块级验证，不把已标注 TLS 上的免疫分数差写成隐藏定位或 H-03 确认。H-01/H-03/H-05 均不升级。`docs/plan.md` 不变。
- 失效条件：若事先改主评分为深度归一化并在另一独立患者队列同向，可重审该一条程序的标签关联；当前不能用这次 Holm 命中回改主终点。

## D-152｜2026-09-28｜组成残差对照先于结果冻结
- 授权：用户要求规划并执行“去掉细胞组成之后还剩什么”。
- 选择：只在 D-148 的 16 位严格患者、31 项无歧义程序上做。组成轴为 B/T/Mye/Epi/Stromal/Plasma，不用已退役的 ILC。每个程序先从各轴删掉自己的 input+readout；剩不到 2 个基因的轴删掉，不补零；剩不到 3 条轴则 NOT_SEPARABLE。Q 为截距、排除本程序基因后的 log10 文库量、以及同样排除后的检出基因数。raw 与 depth 两路线的轴分数和程序分数用同一尺度、同一 lib_excl。回归只在明确 TLS 与 NO_TLS 点上拟合，空白不进回归也不当阴性。
- 主问题只问适应性免疫。一个模式算 REMAINS，必须残差中位差>0.2、保留比例≥0.5、且至少 14/16 位患者残差为正。保留比例<0.5 或正号患者<12 算 COLLAPSED，其余 ATTENUATED。只有 raw 和 depth 都是 REMAINS，才可以说组成之外的标签优势还在。31 项的 Holm 只是上下文，不能改这条规则。
- 排除：不重复加总/丢弃臂；不把这次结果升成定位或跨癌种。门槛文件 infra/gobp_composition_residual_20260928/run_contract.json。
- 失效条件：若事后放宽 0.2、0.5 或 14/16，本条作废，必须另立决策重跑。

## D-153｜2026-09-28｜组成残差对照：适应性免疫未通过双路线规则
- 依据：D-152 与 `infra/gobp_composition_residual_20260928/receipt.json`。16位患者均完成，992行。
- 结果：raw 残差中位差0.48，保留72%，16/16为正，REMAINS。深度校正残差中位差0.57，保留39%，16/16为正，但保留比例低于0.50，COLLAPSED。总体 `DOES_NOT_SURVIVE`。组成轴对 raw 分数的患者中位 R² 为0.93，深度路线0.66。
- 选择：不放宽50%保留线，不把“残差仍为正”改写成通过。不把上下文里两条都 REMAINS 的细胞活化升为主目标。H-01 无新定位证据，H-03 不升级。`docs/plan.md` 不变。
- 失效条件：只有事先改规则并在另一批独立患者上重跑，才能重开“组成之外仍有标签优势”。不能用本次深度残差的正号单独主张。

## D-154｜2026-09-30｜TLS场算子在看曲线之前冻结
- 依据：目标 muniy54y-ukz62g。合同 `infra/tls_field_20260930/run_contract.json`。此时没有新的曲线结果。
- 选择：八名 USZ 病人按编号排序，前五发现、后三留出。瘤内 TLS 只由标签几何定义。曲线只用肿瘤点，按单个 TLS 把径向距离缩到 0–1，并记下 0、0.75、1 的微米。分数是 AUCell。近高远低用本地名称对应的 I 型干扰素、II 型干扰素、MHC-II 和炎症；近低远高用 G2/M 和 EMT。MYC 与 KRAS 本地没有词表，不测、不手补、不算通过。B 轴用同一算子，不得呈近低远高。场不稳则反推不启动。
- 排除：不拼接缩放距离，不扫 200/300/500，不用已封存的空间 p，不改 plan.md。
- 失效条件：若在出曲线后改基因集、留出名单、形状门槛或 IT 定义，本条作废，必须另立决策并重跑。

## D-155｜2026-09-30｜TLS场门未过，反推不启动
- 依据：D-154 与 `infra/tls_field_20260930/discovery_curves.tsv`。21 条曲线，未改门槛，未拼接缩放距离，未计算空间 p。
- 结果：五名发现病人里只有 KC3、LC1 有可评瘤内 TLS，少于合同要求的三名。三份可评曲线全部 flat/nonmonotonic，没有预先指定的近高远低或近低远高。scaled 0 为 100 μm，scaled 1 为 3643–8945 μm，已写入曲线表。
- 选择：场不稳定。留出反推不启动，不制作成功对照表，不加半径、不改门槛。目标按“反推不启动”收束。
- 对假设：H-01 无新支持。这个算子没有找到稳定场，所以反推没有被检验，也不能据此说文献中的 TLS 场不存在。H-02 仍未测。`docs/plan.md` 未改。
- 失效条件：只有另立决策并事先改算子，才能重开反推。不能用这次曲线挑选新半径后再称同一合同通过。

## D-156｜2026-09-30｜TLS atlas Fig-4 复现：正常分支方向可复现但被 B 轴解释，瘤内分支不可评
- 依据：`scripts/tls_pancancer_repro.py` 与 `infra/tls_repro_20260930/`（21 万条 TLS-分箱行、4400 条曲线、220 条汇总）。
- 口径：论文没有用滑动 bin；它用 Semla 径向距离（TLS 中心）+ 0–1 缩放 + Monocle 样条 df=3。本复现额外输出滑窗 bin（用户要求），分类用论文端点规则。第一种草案把 AUCell 秩方向写反，已修复并重跑；以修复后结果为准。
- 结果：正常组织分支（400 μm 内无肿瘤的 TLS，用 NOR 点）干扰素/MHC-II/抗原呈递 7/8 片为正、分类 5 片 high-to-low，方向与论文一致；但 B 轴对照同样 6/8、5 片 high-to-low，因此不能当组成外证据。瘤内分支（Fig. 4 主结果）USZ 只有 3 片 15 个 TLS，方向互相矛盾（KC3 正、LC4 负、LC1 混），不可评，不构成否证。每片最小最大缩放使 distance=1 达 3.6–11.5 毫米，曲线斜率主要来自大尺度梯度而非 50–200 μm 短程场。
- 选择：正常分支方向复现记为方法验证；不把免疫梯度写成新发现；瘤内主结果记为 USZ 不可评；不改 plan.md。
- 失效条件：若获得更多含瘤内 TLS 的队列（新数据需审批），或事先另立口径做瘤旁敏感性，可重审瘤内分支。不能把本次正常分支结果当作图 4 的复现。

## D-157｜2026-09-30｜论文式计算 TLS 检测：可跑，但扩出的是未验证热点
- 依据：`scripts/tls_computational_detect.py`、`scripts/tls_computational_gradient.py`、`infra/tls_repro_20260930/detect_contract.json` v2 与 detect/comp_anchor 三表。报告见 docs/TLS_REPRO_20260930.md 追加节。
- 选择：实现论文后段（签名分数 + 六边空间聚类 + 分割），跳过 iStar（USZ 无 H&E）。检测分数必须用 AUCell 秩分数；第一版均值 z 分数与文库深度秩相关高达 0.88，已按冻结合同 v2 弃用。锚点自己的点从读出中排除，保持输入/读出不相交。
- 结果：检测器对病理 TLS vs 正常点 AUC 0.61–0.93（KC1/KC3 0.93），LC4/LC5 0.29/0.13。top-10% 下瘤内组件 42 个，41 个不含病理 TLS 点、尺寸 2–4 点、距最近已知 TLS 200–2400 μm。用这些锚点复算：免疫侧 6–7/7 片为正（循环，检测用它），G2M/MKI67 5/7 片为负（探索级，非循环），MYC/EMT/COL1A1 混杂。
- 决定：不把计算锚点当已验证 TLS，不据此宣布复现论文主结果；记录为探索性扩样本与敏感性。LC4/LC5 病理 TLS 为低深度点（文库中位 260/278 vs 正常 655/6953），两片 TLS 导向的分析需单独标注。
- 失效条件：若能取得论文原始 TLS 签名与 iStar（或 USZ 的 H&E），或引入独立组织学验证（FDC/GC 网络标记），可重审计算锚点的 TLS 资格；在此之前不得把锚点写进确认性结论。

## D-158｜2026-09-30｜跨队列/平台/组织 TLS-signature 池扩到 55 section、2859 组件（只计数）
- 依据：`infra/tls_pool_expand_20260930/run_contract_v3.json`（FROZEN_BEFORE_SCORING）、`scripts/tls_pool_expand.py`、`pool_diagnostics.tsv`（55 行）、`pool_components.tsv`（2859 行）、`pool_summary.json`（55/55 ok）。报告见 docs/TLS_POOL_EXPAND_20260930.md。
- 选择：同一签名（B 49 + 浆 7 + 趋化因子 4，共 57 符号）、同一 AUCell 算子、同一 top-10% 阈值、同一六边组件规则跑全部本地队列；平台只换几何（Visium 系 100 μm；HD/Xenium 先聚 55 μm 伪格）。重复符号跨行加总。不拟合梯度、不反推、不算 p。
- 结果：USZ 8（430）＋GSE175540 18（1008）＋ST-CRC 14（364）＋HTAN 8（452）＋Cervilla 对 2（76）＋10x 公开 CRC 2（250）＋HD 1（125）＋Xenium 2（153）。标签验证：USZ 0.63–0.94（LC4/LC5 异常沿用 D-157 标注）、GSE175540 0.56–0.96、HTAN 0.66–0.999。置换 null：Visium 系 obs/null 约 1.3–3.5（HD 12.5）；Cervilla CytAssist v2 0.95（约等于 null）；Xenium 0.6–0.7（大团块方向，panel 受限）。实现层修了两处不改算子的 bug（脚本入口缺失、tar 二进制成员被当文本解码、老 numpy 无 trapezoid 改手工梯形），重跑后 55/55 全过。
- 决定：池子只用于样本量与诊断；组件仍是 B-签名热点而非已验证 TLS；数量不跨平台直接比；ST-CRC 的 IC aggregate 类维持 FAIL-CLOSED；不改 plan.md。
- 失效条件：任何把池内组件直接当 TLS 证据用的分析，必须另立合同并先过独立组织学或签名外验证；不得用本池数量论证生物学结论。

## D-160｜2026-09-30｜场筛查改用 Hallmark（用户推翻 D-159 小 panel）＋P0/P1 先行门控
- 背景：用户要求必须下 MSigDB Hallmark，不退回小 panel；乐观才切 GOBP（归并减规模）；去掉新增外部数据锁；同意 P0–P3 门控意见。GMT 文件已在本地（`data/geneset/h.all.v2026.1.Hs.symbols.gmt`，50 集，32–200 基因/集，中位 199，sha256 `eecaf6…f596`），本次无下载动作；D-111 对未来拉取仍然有效。
- 选择：读数＝Hallmark 50＋B 轴（对照兼晕问题本身）；mask＝组件＋一圈隔离；距离主口径绝对微米公共 bins，rescale 只做论文仿真敏感性；P0 角度分箱＋P1 走廊对照先行，P2/P3 以阳性为开门条件；Xenium 冻结；GOBP 归并只在乐观后另立冻结。合同 `infra/tls_field_20260930/run_contract_v4.json`。
- 失效条件：出结果后改 panel、改门槛、把 Xenium 拼入计数、或用免疫读数单独宣布发现，本条作废重冻。

## D-161｜2026-10-01｜Hallmark 场筛查 P0/P1：B 晕有、走廊无；P2/P3 不开
- 依据：v4 合同、`scripts/tls_field_hallmark.py`、`field_curves/angular/summary_v4{,r}.tsv`。
- 结果：径向 B 轴 52/53 片为正（中位 drop 0.012）；免疫集 72–92% 为正但量级千分之几且循环；G2M/MYC 平坦；EMT 反方向（基质共定位解释）。角度走廊 51 集全部归零（中位 ±0.002）。
- 限制：锚点 top-10% 选择效应无梯度 null，只能描述；汇总指标因掩膜改过一次（首个可用 bin），门未动。
- 决定：P0 径向对 B 描述性通过，角度与 P1 走廊为阴性 → P2/P3 不启动；H-01/H-03/H-05 不升级；plan.md 不变。
- 失效条件：任何把 B 晕或免疫方向写成淋巴结构定位证据的用法，必须另立合同并先过独立验证与梯度 null。

## D-162｜2026-10-01｜场筛查审计＋选择效应校准：B 晕是选择效应＋局部 B 场追踪，无独立场（H-01 不升级）
- 审计修正：（a）Cervilla loader 坐标错位——组件位置全错，池 2859→2900、Cervilla 76→117、obs/null 0.95→1.58/1.90（场筛查在此修复后运行，未受影响）；`components()` 加长度断言；D-158 正文的 2859 以本条为准。（b）角度走廊口径错：只有 USZ 的标签值是 `TUM`（HTAN=`carcinoma`、STCRC=`tumor`），v4 把 30 片当有标签，实际 22 片扇区 0 为任意方向；USZ-only 重算仍阴性。（c）HD 巨组件（4612/3304/1789 点）使质心距离失效。（d）循环性表述更正：B 轴 49/49 与检测签名重合（真循环），IFN 系仅 3–4/200（非循环），31 集零重合。（e）池表 STCRC 364→365。验证：53 片 Visium 上场锚点与池组件完全一致（2747/2747）。
- 主臂（条件化残差，v5.1）：原始近-远对比为正（B 轴 53/53、+0.0092；IFN-γ/α 44/53）→ 扣掉“本点 B 水平”后归零（B −0.0000、IFN-γ +0.0007）→ 再扣“周围 300 μm B 场”转负（过度扣减，不算反晕）。事前规则要求两条条件化口径都 ≥60% 为正，六个读数全部不通过。
- 次臂（匹配自相关代理场，200 抽/片）：178 个可评估检验中 2 个 p ≤ 0.05（1.1%，零假设期望 5%），p 中位 0.49，B 轴实测中位 0.105 < 代理场 0.170；单 σ 拟合对长程成分失配（实测 ρ(700–1000) 中位 0.24，拟合预测 ≈0），方向偏保守。
- 幂验证：注入 0.25 倍标准差人造晕即可在 9/9 片检出（统计量 +0.007，实测仅 −0.0000/+0.0007），零结果不是管线无电。
- 结论：B 晕＝峰值选择＋读数随局部 B 场水平变化；无独立于 B 场之外的分子晕/场；走廊阴性（USZ-only）。
- 失效条件：若将来引入**组成对照**（C6 协变量）或独立组织学锚点后重新检出方向一致的晕，本条需重开；任何引用“无晕”结论时必须同时说明组成对照尚未执行。报告见 docs/TLS_FIELD_NULL_20261001.md。

## D-163｜2026-10-01｜GOBP 压缩规则（用户选定 keep-specific）＋compressed panel 晕扫掠（screen-only）
- 背景：用户要求用 Hallmark 之外的 GOBP 扫一遍分子晕；先做压缩。本地已有 v2025.1 的 c5.go GMT（`/home/huyudi/006/data/pathway/`）与 GO 本体（`012_conference/…/go-basic.obo`，releases/2026-06-15），零下载。
- 数据：GOBP 7583 集；名称匹配本体 7331（96.7%），252 未匹配（旧名/弃用）；仅尺寸 ≥20 基因保留 3808，仅层级去后代保留 867。用户选定 **keep-specific（每条谱系留最具体者，floor 20 基因，Jaccard 0.7）→ 1691 集**，面板文件 `infra/gobp_compress_20261001/survivors_specific_floor20_jac70.tsv`，sha256 `512842dd…f7ac`。失效条件：若换用“保留泛化”或加 depth 窗，需另写一条。
- 扫描口径：与 D-162 条件化臂一致（边缘距离、每锚点近-远差、raw/own/local 三口径）；实现改为“每片只排一次秩 + 一次稀疏乘法打 all sets”，开跑前用 5 个随机集与 `aucell` 逐点比对 ≤1e-6 通过；53 片 × 1691 集 919 秒。合同 `infra/gobp_halo_20261001/run_contract_v6.json`。
- 结果：raw 中位为正的集 52.3%、own 后 38.5%；中位 raw +0.00005 → own −0.00021。片内标签置换零假设：own ≥0.7 实测 13 对零假设中位 4（上限 10），own ≥0.65 为 59 对 32，raw ≥0.7 为 101 对 16；集大小与 own_pos 无关（Spearman −0.025）。
- 结论：**筛选结果，NOT_CALIBRATED，不作声明**。top 名单里补体经典途径、肽抗原-MHC 组装、抗菌体液反应与 v5 方向一致；心肌形态发生、剪接位点识别等项生物学上无法解释。下一步（合同已预登记）：只对 own 前 50 集跑匹配自相关代理场 null。
- 报告：docs/GOBP_HALO_SCREEN_20261001.md。

## D-164｜2026-10-01｜GOBP 前 30 集空间校准：全部归零，筛选超额未转成显著性
- 对象：v6 按 own_median 前 30 的 GOBP 集 + 3 对照（B 轴＝循环阳性参考、MYC_V1＝干净对照、随机 32 基因＝噪声底线），冻结于 `infra/gobp_halo_20261001/target_sets_v7.tsv`（sha256 `dd114df5…`）。
- 零假设设计与 v5 不同（合同 v7 跑前说明理由）：v5 的“混合＋重挑锚点”是为 raw 统计量造的；v7 校准的是条件化残差，它在零假设下已与锚点几何无关，故**固定观测锚点、只随机场**（同 σ 平滑高斯场 + 分位数映射），200 抽/片。用 v5 设计会把显著性抬高。
- 结果：33 个集 claim 全为 0。补体经典途径（v6 own_pos 0.80）校准后切片级显著比 0.12，肽抗原-MHC 组装（0.71）为 0.03。对照全在预期位置：B 轴 0.045、MYC 0.057、随机集 0.0（观测中位 −0.001）。整体 644 个切片级检验 8.4% p≤0.05（名义 5%），超额集中在只有 2 片的 Cervilla（17%）与 Parent（20%），大队列 5–9%；目标集 8.9% vs 对照 4.3%。
- 结论：v6 筛出的“扣掉本点 B 场后仍存的小额正差”未能转成空间显著，与 D-162 在 Hallmark 上的判断一致；B 热点周围无超出平滑场与局部 B 场的独立分子晕。
- 限制：可评估切片少（3–35），“无证据”而非“有效证明无”；VEGF 负向调控集 σ 不可拟合（n=0）；组成对照仍未执行。
- 失效条件：若获得独立组织学锚点或补上组成对照后出现跨切片一致、能过固定锚点零假设的信号，本条重开。报告：docs/GOBP_HALO_NULL_20261001.md。

## D-165｜2026-10-01｜独立锚点（病理 TLS 标签）场检验：无场，循环性反驳被排除
- 动机：v4–v7 的锚点全部是签名热点（B/浆/趋化因子 top-10%），与读数同源；改用病理标签定义几何。
- 设计变更（跑前探测驱动，v8.1）：逐锚点径向廓形不可行——26 片中仅 5–6 片能达到 ≥8 个可用锚点（远带 400/500/700/600+ 四种口径分别是 5/5/6/6），因为标注点少且碎、Voronoi 领地达不17远处；合并组件只会更差。在未算任何读数统计量前改为主统计量＝**片级** TLS 距离廓形（近 200–300 vs 远 700–800，各≥30 点），更贴近 Cho 2026 RadialDistance。v8.2 增补：σ 拟合不出时用同边际白噪声 fallback（只补回被丢弃的检验，不改已跑的，已逐行标注，方向偏激进）。
- 结果：USZ 8 + GSE175540 18 片，20 片进入分析；968 个切片级检验 6.5% p≤0.05（名义 5%）；**无读数通过事前判定**（最高 0.20，需 0.50）。对照：B 轴 0.05、MYC 0.00、随机集 0.00。剔除弱标签片后 7.0%；单队列 6.1–6.8%。
- 意义：用与读数独立的病理标签重做后结论不变，“结论由签名锚点循环性造成”这一最强反驳被排除；B 轴在独立锚点下降到名义水平，印证当年的“B 晕”是签名挑锚点的产物。
- 限制：“无证据”而非“证明无”；半数检验用白噪声 fallback；标签是病理区域注释非分子 GT。
- 失效条件：若换成分子级 TLS GT（如空间 TCR/浆细胞原位验证）或加入组成对照后出现能过固定掩膜零假设的跨片一致信号，本条重开。报告：docs/LABEL_ANCHOR_FIELD_20261001.md。

## D-166｜2026-10-01｜组成对照补齐：残余超额不被组成解释，也无集通过；场线封口
- 目的：v4 门里唯一未执行的“beats composition controls”。
- v9.0 无效：R-04 冻结 NNLS 参考是 **CRC 单细胞**（GSE236581），套到本轮的肾/肺（USZ）与 ccRCC（GSE175540）上退化——usz-LC3 均值 T=0.00/B=0.63/Epi=0.00，gse…c_2 的 Epi=0.06；参考自身 profile 重建完美（各自 1.0），故是**领域不匹配**而非代码 bug。v9.0 结果不承担结论，仅存档为 `*_v90_nnls.tsv`。
- v9.1 有效：改用跨组织 marker-proxy 模块（CellMarker2.0 + PanglaoDB + CellTypist CRC，加 r16.PLASMA_GENES）的 T/Mye/Epi/Stromal/Plasma 分数（AUCell，转十分位）作协变量；排除 B（锚点相邻轴，已由签名分位代表）与 ILC（非 viable）。跑前有效性检查通过（T/B/Plasma 在 TLS 升、Epi 在肿瘤升）。
- 结果：own 臂 968 个切片级检验 7.4% p≤0.05、comp 臂 **6.8%**（USZ 7.1%→6.1%，GSE 7.6%→7.3%）；**无读数通过事前判定**（最高 comp 臂 0.25，需 0.50）；对照 MYC 0.00、随机集 0.00。
- 合并结论：签名锚点（v5/v7）、独立病理标签锚点（v8.2/D-165）、组成对照（v9.1）三个方向一致——无独立分子晕；“锚点循环性”与“组成”两个最强反驳均已被排除。残余 ~1.8pp 超额散在单切片，不聚集于任何通路。
- 限制：无分子级 TLS GT；Visium 100 μm 分辨率；20/26 片可评估；半数检验用白噪声 fallback；横断面无串行切片（I-002）。
- 失效条件：分子级 TLS GT 或更高分辨率数据出现能过同一固定掩膜零假设的跨片一致信号时重开。报告：docs/COMPOSITION_CONTROL_FIELD_20261001.md。

## D-159｜2026-09-30｜B 热点场筛查方案冻结（含三条反驳）：内圈一环、冻结小 panel、P0/P1 先行、Xenium 冻结
- 背景：用户确认内圈隔离一环、读数取免疫＋B＋肿瘤内在、P0–P3 全开、Xenium 先冻结，并要求主动反驳。
- 反驳一（不跑全 GOBP）：本地 GO 定义 6987 个、中位仅 6 基因，AUCell 在小集合上噪声大；2859 锚点 × 6987 集 ≈ 2000 万条曲线，不可行；45 候选中已有 42/990 对 Jaccard≥0.5，全量跑等于把高度相关的检验数几千倍，且通不过多重检验。Hallmark 本地无文件，下载即新增外部数据，需另行批准。冻结小 panel：免疫（IFN_I/II、MHC_II、炎症、抗原呈递）＋B 轴＋干净肿瘤内在（G2M/EMT/MYC/MKI67/COL1A1）＋少量单基因（CXCL13/CXCL10/CCL19/CCL21/TNFSF14 等）。
- 反驳二（不四条全开）：Visium 点距 100 μm、点径 55 μm，200 μm 内只有 1–2 圈，4 扇区×多径向带在近端天然缺格；P2 梯度向量与 P3 定向滤波在该分辨率下先天不足。P0 角度分箱＋P1 走廊对照先行；P2/P3 以 P0/P1 阳性为开门条件，不并行铺开。
- 反驳三（Xenium 冻结，同意用户）：Xenium 是亚细胞单细胞 regime，55 μm 伪格只是投影；其 panel 仅覆盖签名 20/57，阈值与连通含义不同。冻结 Xenium，只做 Visium 系（含 HD 伪格敏感性）。
- 选择：mask＝组件＋一圈隔离；读出按上条冻结 panel；距离用绝对微米公共 bins 为主口径（各平台按分辨率取 bin 宽、分段报告缺格），rescale 0–1 只做论文仿真敏感性；场判定线＝掩膜外有信号＋随距离有结构（单调或有方向）＋跨锚点中位数一致＋赢过组成对照（B 轴形状/幅度＋肿瘤内在干净读数）；热点仍只是漏斗入口，不升格。
- 失效条件：出结果后改 panel、改门槛、把 Xenium 拼入 Visium 计数、或用免疫读数单独宣布发现，本条作废重冻。

## D-167｜2026-10-02｜结论审计：先撤回过强否证，修复数值算子并隔离历史输出
- 依据：当前 v6–v9 的 `residual_operator` 对秩亏设计使用未截秩 QR；两个 sigma 拟合器的 x 已除以 4、反解又除以 4；v8/v9 把无法拟合相关长度改作白噪声。解析回归测试先复现前两项错误。
- 选择：使用 SVD 的有效列空间投影，修正 sigma 反解；模拟读数也经过同一残差投影；相关长度不可拟合时不生成 p、不作为阴性。修复后空间模拟仍为探索性（NOT_CALIBRATED），正式 claim 留空，另报原阈值的探索性通过状态。保留原合同与结果，用独立 audit 输出登记变更，不事后改写原冻结定义。
- 排除：不把修好公式等同于真实空间 null 已校准；不靠白噪声 fallback 补齐分母；不把 marker 分数代理等同真实细胞比例；不把无显著性等同等价性或普遍不存在。原 v7 top30 来自受影响的 v6，重算只作固定历史面板敏感性，不能声称修复后的全 GO 确认。
- 执行：CPU、本地数据、独立目录，保留 1.2 TB 空间；完整重算最新独立标签与组成代理两臂，追踪排除原因。其他受影响旧 null 的结果明确作废，不用局部重算冒充全线恢复。
- 对假设：撤回 D-162–D-166 中“无独立场／排除组成／场线封口”的普遍结论；H-01/H-03/H-05 尚无本轮合格确认或否证。H-02 仍未测，R-05/R-06/R-07 保持 NOT_TRIGGERED，不启动新科学路线。plan.md 不变。
- 失效条件：只有经同操作点误报率/功效验证、无读数自扣减的组成对照和合格患者级证据，才可重新评估正式结论；当前修复不自动解锁。
- 开跑前补充：本地 GSE175540 positions 与 USZ x_array/y_array 均为 Visium 双倍列坐标，`hex_xy` 却再次加奇偶偏移且未除以 2，六个一阶邻居距离不等。修复为 (col/2, row*sqrt(3)/2)*pitch；影响使用该 loader 的旧组件与距离结果，全部降为历史无效几何输出。最新标签队列按源标签重算，不重用旧组件；旧池表仅用于冻结样本名单，不能再作几何证据。R-04 与未调用该转换的历史路径不据此一概作废。

## D-168｜2026-10-02｜用户确认保留外围组成信号，组成独立性降为辅助问题
- 来源：本项目用户明确提出“TLS 附近免疫细胞多本身就是有价值的结论，不应该严格扣除”，并认可主分析保留外围组成信号、组成拆分作为辅助解释的调整。
- 科学理由：核心不可见时，外围细胞丰度及状态都可能携带结构信息。组成解释不能自动否决这种信息；“外围有可用信息”和“同类细胞存在组成之外的状态变化”是不同问题。
- 选择：修改 plan P-01/H-01/K-01/P-06，H-03 保留为次级机制假设。主分析以独立病理标签定义 TLS、遮蔽核心与隔离带，保留外围免疫/B轴/组成相关读数，不预先回归掉检测签名或组成代理；比较组成-only/简单marker模型用于说明复杂信息是否有额外价值，不作为所有外围信息的否决门。
- 保留门槛：技术混杂、空间误报率/功效、患者级重复性、标签泄漏和预测留出仍需处理；只靠表达不声称细胞数增加或 TLS 因果招募。已知锚点的距离曲线不等于定位成功，未来定位输入不得由答案生成可见 mask 孔洞。
- 本次范围：准备下一节点与分析设计，尚未冻结具体统计合同、未运行新分析；旧 own/comp 结果不改名为新主分析，不恢复旧强阴性。现有23片资格只属于旧距离口径，新的可评样本须重新由metadata和冻结mask核定。
- 门禁变化：I-032 不再阻塞外围信息主分析，仍限制“组成之外”的解释；I-031 的空间校准阻塞保留。R-05/R-06/R-07、GPU及既有停止路线不自动重启。
- 失效条件：若今后把组成独立性、因果机制或定位成功重新纳入主张，必须补齐对应证据并另记决策；不能以本次目标调整追认历史未校准结果。

## D-169｜2026-10-02｜R-17执行：原始计数、独立标签mask与空间校准先行
- 来源：用户授权按D-168开始下一步。新评分尚未计算。
- 数据选择：26片既有名单；USZ改读随数据提供的Space Ranger filtered原始计数，h5ad只读取条码/病理标签；GSE175540读原始filtered矩阵与官方positions/annotation。不使用表达平滑，不用隐藏核心做跨spot预处理，重复基因符号按同符号原始count求和后逐spot排序。全部Hallmark50与既有B-axis组成51读数，冻结清单，不使用旧筛选top30。
- 几何：已修正的100μm双倍列六角坐标，主mask为到任一TLS点距离≤130μm，敏感性≤230μm。固定近带[200,400)、远带[700,1000)，每带至少30点、核心至少20点；两mask分别报告可评资格与可见点。全距100μm分箱只作描述；不按结果调整半径。AUCell在可见spot的全基因背景逐行计算，不回归B/组成。主要效应为近远均值差及可见spot SD标准化差；0.2 SD为本轮预先规定的有用效应参考，不是已知生物学阈值。
- 校准候选：固定五分量（nugget及相关长度100/300/800/1600μm的高斯核）非负变异函数拟合，空间相关性保留于组织实际可见坐标；拟合失败不回退白噪声。用实际mask上多个已知单/多尺度、高斯/非高斯与趋势场检验误报及注入功效，先冻结具体模拟设置。模型内校准通过不自动外推真实数据。失败则仍完成全数据描述，但禁止正式空间p及阴性证明，不事后替换null救结果。
- 身份：现有USZ登记patient源实际仅为SAMPLE编号，不能作明确患者映射；GSE atlas有24个不同患者的队列级旁证，但原GSM与R_P编号交叉表需核实。当前严格patient_id保持未知，不将样本汇总当患者级推断。详列源字段与历史登记冲突，旧表不被静默当作权威。
- 选择排除：不直接平移病理mask当可交换阴性（组织支持/unknown存在问题）；不打乱spot；不把样本bootstrap当患者CI；不新增数据或启GPU。R-17结果只回答本地已知标注核心外关联，H-02不启动。
- 失效条件：看新信号后改panel、mask、效应方向、统计量、校准门槛或患者认定，需要另立决策，不能称本合同通过。

## D-170｜2026-10-02｜R-17全量执行完成：B轴外围候选保留，校准失败不升级
- 依据：`infra/tls_peripheral_20261002/contract.json`（sha256 0ef70069…95ccf6）、qualification/effects/calibration与validation receipt；报告 `docs/TLS_PERIPHERAL_RESULTS_20261002.md`。
- 结果：26片全部核查，24片×2mask×51读数=2448效应行、24480分箱行；51读数均满足冻结覆盖阈值。B轴主/扩圈均USZ8/8、GSE13/16近高远低，中位效应主0.436/0.340、扩圈0.445/0.359 SD。保留为描述性外围关联，不将组成相关性提前消除，也不声称细胞数或因果机制。
- 校准：完整24片×2mask×6类场=288个null设置，每设置400次；另1152个功效设置，每设置200次。56/288未过误报上界门，510/576的0.5SD设置未过功效下界门。候选检验FAILED，不输出正式p/正式多重检验声明，不以阴性证明无场。
- 身份：源metadata不支持当前严格逐例患者映射，I-033硬阻塞患者级确认；汇总仅按切片。USZ旧patient_id推断已在新分析与registry overlay隔离，历史冻结表不重写。
- 选择：R-17执行结束、状态DESCRIPTIVE_ONLY/CONFIRMATION_BLOCKED；H-01仅描述性候选，H-05非独立迁移确认；H-03机制与H-02定位未运行。R-05/R-06/R-07不启动。plan保持D-168版本。
- 失效条件：后续只有另立预登记的适用校准与身份依据才能升级；不在本合同下调null、半径或panel抢救结果。

## D-171｜2026-10-02｜全量slim GO BP与多尺度桥：用户新目标授权
- 来源：用户要求独立完成全量slim GO BP扫掠、检验TLS↔分子场/晕双向关系，并重点研究TLS–肿瘤区域、TLS–肿瘤边缘、TLS–TLS桥；随后明确“并不是说桥一定是毫米级，只是要求在多个尺度测试”。
- 科学变化理由：既有短程径向统计不能回答方向性连接及跨尺度形态；新增H-08，桥长度与宽度同时测试，不以“必须毫米级”或“必须远大于TLS”作为资格先验。相对TLS尺寸只作为结果维度。
- 全量定义：既有keep-specific/floor20/Jaccard0.7面板1691项（项目压缩GO BP，非声称官方GO-slim本体）；保留全部项的覆盖/不可测记录，禁止换成top30或子集交付。沿用原始counts逐spot评分、保留组成信息；本地数据，无GPU或外部下载。
- 数据资格：先完整处理26片；与前轮不同，小TLS不是统一按20点排除，按目标实例与尺度的几何支持分别判断。缺失已知TLS的片几何确认仍受限；未知标签不作阴性。USZ有TUM，可支持肿瘤相关目标；GSE只有TLS/NO_TLS，不虚构肿瘤区域/边界。
- 前向分析：晕、边、桥分开登记。桥包括直线及组织图上的几何路径，分别检查沿程连续性、平行侧带、端点遮蔽与端点径向叠加解释；表达驱动路径若使用，必须把路径搜索也纳入校准，不能将画出的连线当桥。
- 反向分析：按切片留出、内部验证；每个候选查询点使用同一规则遮蔽，不暴露答案生成的mask孔洞、TLS距离或真实端点。全1691项的筛选/组合仅在训练折内；简单B/免疫、深度与空间基线一起比较。患者独立性不明时不宣称患者级或外部确认；双向指结构→分子关联与分子→结构预测，不指双向因果。
- 范围：本目标授权开展新桥与反向检验，不因前轮R-17未过而只做短程B轴；但前轮失败的p不能复用，正式结论仍需适用校准。新统计合同按metadata/几何先冻结，未计算新GO效应。
- 失效条件：任何事后改尺度、panel、端点定义、训练/测试分组或用数据驱动路径而不重做选择校准，都不得称原合同确认。

## D-172｜2026-10-02｜多尺度前向算子细化（首次GO效应查看前）
- 几何清点：USZ三类桥共5025对端点（TLS–TLS605、TLS–肿瘤区域4159、TLS–肿瘤边缘261），GSE TLS–TLS23对；这是标注连通域/固定地理tile的候选对数，不是已经存在的桥，也不是独立生物重复。所有候选和宽度/长度保留资格记录。
- 走廊：中心半宽w，两侧各(w,3w]；去除所有TLS核心和隔离带及目标端点disk后，将剩余中心走廊的几何投影范围等分5段，每段中心和两侧各需≥5点。端点遮蔽导致短尺度不可测时记录，不合并为阴性。直线与组织图最短路径都计算；不使用GO表达选择路径。
- 桥效应：每段中心均值减两侧均值的平均；5段均值为整体效应，至少4段同号仅作连续性描述。两侧差异单独记录；宽走廊缺少侧带是不可测，不放宽阈值。
- 端点晕基线：在中心走廊外的可见spot上，拟合两端距离的固定指数基（100/250/500/1000/2000/4000μm）与截距，预测走廊/侧带的端点叠加贡献；报告原始及剩余对比。使用训练外mask的线性算子，不读取核心/端点disk表达。基线无足够支持时明确不可测，不能用空缺当“超出晕”。
- 晕与边：每张片、两mask的完整预定距离区间统计，至少30点的bin才可评；相邻径向bin对比以及TUM与明确非TUM界面分带另列，不冒充桥。身份、真实空间p与两方向确认限制保持。
- 失效条件：若事后换路径代价、分段或侧带宽度以获得连续信号，需要新合同并重做全选择校准。

## D-173｜2026-10-02｜反向验证合同：固定网格、统一mask、逐切片留出
- 查看前向GO效应前冻结。查询网格由双倍列%4=0、行%2=0确定（约200×173μm），覆盖整张源组织；不把TLS中心挑出来另加查询点。核心遮蔽半径230/500/1000μm逐档完整测试。
- 真值：查询中心为源TLS标签且所属标注连通域≥2点、整个域再加130μm隔离均在统一圆mask内时，为可评正例；中心明确NOR/TUM/INFL/LN/NO_TLS为中心标签负例。负例指中心不是TLS，不宣称整窗口没有TLS。UNKNOWN不改负例；未完全藏住目标域、缺失已知TLS的片不进入该操作点评价。
- 输入：只读取查询圆外；距离壳为R→R+250→R+750→R+1750→R+3750→组织可用极限。每壳8扇区均值；5个径向均值与一、二阶角向幅度各5个，描述方向性/线状信息而不输入真实TLS位置或桥端点。径向少于20点、角向任一扇区少于5点时该特征中性补0并报告几何缺失；前两壳均≥20点、总可见≥100才可评价。按每查询可见表达均值/SD标准化，不能用被遮蔽表达算分母。
- 模型：每个GO与51固定对照分别使用15维特征的固定ridge（协方差标准化后λ=1）预测中心标签，按切片/类别等权训练。再在每个训练折内按训练得分选前10个GO组成联合模型；这是训练内选择，外折标签不用于选集/调参。B轴、Hallmark固定对照与几何/深度基线同输入预算。全1691项都尝试，缺覆盖/常数保留不可测。
- 验证：全部逐切片留出；另做USZ↔GSE来源留出。AUC/AP仅按明确中心标签计算并报告各折正负数；缺少一类时AUC不可测。置信/独立单位限制为I-033，不把重叠查询窗当患者。固定panel在旧数据上有探索历史，当前为内部验证，不称全新独立盲验。
- 失效条件：加入答案形状mask、TLS距离/端点、在留出数据选GO/尺度、未知当负、在结果后换查询网格，即不能沿用本合同结果。

## D-174｜2026-10-02｜反向模型实现细化：方向性增量与训练内选择
- 首次训练/查看反向效果前补充：对全1691项及固定对照分别报告径向5特征、角向10特征、联合15特征三臂，便于区分“晕信息”和方向性信息。角向信息不是桥的充分证据，仍需前向完整路径与对照支持。
- 按切片/类别等权计算训练一、二阶矩，ridge用训练协方差标准化、λ=1。联合模型前10 GO按训练切片AUC的宏平均排序，仅在全部训练来源具备基因覆盖的项中选择；留出片不可测选中程序时联合模型记录不可测，不用留出覆盖信息重选。
- 简单几何/深度基线同样只用查询圆外信息。缺单一类的留出片保留预测，但AUC/AP不冒充有效二分类评价；来源留出两方向均执行，有资格限制时明确记录。
- 失效条件：引入留出标签决定模型、top10或方向臂，以及把角向增量直接命名为桥，均不合法。

## D-175｜2026-10-03｜多尺度桥的真实场旋转参考与校准
- 尚未按GO效果挑选桥或通路。新增完整原位参考：每条已可评走廊绕端点中点旋转4.5°至175.5°，共39方向，长度/宽度/几何路径形状不变，沿程5段使用真实走廊在端点mask后剩余的同一投影范围。
- 主参考只匹配几何支持：每段中心/两侧仍各≥5点，各自点数为真实走廊的0.5–2倍。实际TLS核心与目标端点mask保持相同；不按表达选择参考。TUM与UNKNOWN比例差≤0.2单列为次级组织区室敏感性，不通过这一条件否决允许组成贡献的主信号。
- 所有1691项与固定对照在全部合格参考上执行同一中心−侧带、沿程符号一致性、端点晕基线调整，不抽GO子集。长宽比单列：长度不足两倍半宽的结果称宽片/宽带，不仅因对比高就命名桥；仍保留完整结果。
- 旋转只能提供参考，不自动证明在不规则组织上的方向交换性。正式校准须在实际支持上检验完整方向/尺度/路径选择，且真实与参考使用相同候选集合；参考不足保持不可测。旧R-17失败p不复用。
- 计算实现：端点基线采用与原算子代数等价的“可见域低秩矩+中心局部修正”，避免为每个参考存储整片稠密回归算子；须先通过独立最小二乘等价测试。全部原始和调整后参考分数留档。
- 失效条件：按表达筛参考、事后只保留最强桥尺度、忽略参考不足或将旋转参考直接当独立患者证据，均不能据此确认桥。

## D-176｜2026-10-03｜原位旋转参考的全量汇总规则
- 尚未按GO效果挑选通路或桥。每个已可评几何保留实际合格旋转数，全部1691项及51对照均汇总调整后效应相对旋转参考的均值差、标准化差和经验分位；这些分位不是已校准p值。
- 用于主要方向性比较的参考支持要求：39个候选方向中至少19个可用，且0–45、45–90、90–135、135–180度四个区间各至少2个。19个参考只提供约5%尾部的粗分辨率；此门槛不能代替空间交换性或多尺度选择校准。不满足者仍保留所有得分，明确为参考不足，不作阴性。
- 长度小于两倍走廊半宽的几何保留为宽片/宽带；桥候选汇总另列长宽合适者。正、负两种沿程方向均保留；至少4段同号仍仅为描述性连续性。主要参考保留组成信息，TUM/UNKNOWN匹配子集另列为敏感性分析，不作为主结果必须通过的条件。
- 先在源TLS实例内汇总，再在切片内汇总，不把大量端点或旋转当独立患者。全部长度、宽度、mask和两种几何路径均交付；不从汇总中选择一个尺度冒充预先指定检验。
- 失效条件：以程序结果决定参考支持门槛、丢掉不利尺度或将经验分位写成正式显著性，需要撤回对应结论；本条不解除I-031/I-033。

## D-177｜2026-10-03｜真实旋转方法的必要校准门：共同候选与完整几何搜索
- 新GO效果尚未用于挑选桥。真实旋转不是天然可交换检验；先验证其必要条件，不能因有39个参考就直接报告显著性。
- 几何集合：对每片，取39个方向全部满足原支持匹配、且长度至少为走廊全宽的几何。只在这一共同集合上，让真实方向及39个旋转方向执行完全相同的路径、mask、宽度、长度和三类端点搜索。其他几何保留D-176完整描述，并明确此校准不覆盖，不能称阴性。
- 统计量：每个几何的40个调整后对比作对称标准化；仅保留原始五段中至少四段同号、调整后效应及标准化差同向的读数。每个方向取所有共同几何的绝对标准化差最大值；真实方向在40个值中的保守尾部比例作为诊断秩，含真实方向、相等值计入尾部。真实数据逐1691项及51对照完整输出，另记录1691项联合最大值，不把事后最强GO或尺度当预定单检验。
- 必要校准：每片在实际组织点与mask上生成六类明确已知且无TLS定向连接的场，各400次。四类为长度150/500/1500/4000μm的高斯谱基场（512个固定随机频率），第五类混合上述四长度（方差权重0.4/0.3/0.2/0.1），第六类为混合场的偏斜变换、位置依赖方差与固定平面梯度。用固定seed复现，不拟合真实GO表达，不把这些场声称为已验证的真实GO零分布。
- 每次模拟都执行上述完整共同几何搜索；诊断秩≤0.05计一次误报。单侧95%二项上界≤0.08为该设置通过。任何必要设置失败，旋转交换性不能通过，不再输出正式空间p；全部设置通过也仅表示必要门通过，不自动证明1742维联合表达场适用性。
- 功效：在共同集合按端点类型、路径、mask、长度档、半宽分组，固定取几何编号中位数，不按表达选择。在对应中心走廊注入0.5/1.0可见域SD，使用前200个模拟背景；要求注入目标本身相对所有方向/所有共同几何足够极端，避免将其他位置误报计成功效。0.5SD的单侧95%二项下界≥0.8才允许该设置支持有用效应的阴性解释；不达标保留不确定。
- 这是新的旋转检验适用性审计，非旧R-17 p重用。校准失败、共同参考不足或患者身份不明，均不得被总体平均掩盖。晕和边的适用对照仍须单独完成。
- 失效条件：根据新GO结果或模拟失败事后更换共同集合、零场、门槛或功效定义，不得沿用本条通过状态。

## D-178｜2026-10-03｜完整反向与桥校准的解释边界
- 84个反向留出任务及全1691 GO、51固定对照、三种信息臂已完成。500μm条件存在内部预测增益；1000μm跨来源时B轴优于训练内选择的GO组合，复杂组合不作普遍优势主张。保留全部mask和两方向来源留出，不在结果后只报告最好尺度或重新挑组合。
- 一些几何/深度基线低于随机水平。全量汇总增加对0.5 AUC及固定比较项逐片最高值的差，防止将“胜过失败基线”称为有用增益；后者仅为保守描述，不是测试集选择的部署模型。不改变已完成的模型或预测。
- 原位参考全量完成，必要旋转校准35/48误报设置失败；767/888个0.5SD功效设置失败。禁止将原位分位/诊断秩改写成正式p，禁止普遍否证桥，不继续调这个合同的门槛抢救。实际空间图仍可作为带明确限定的探索性证据。
- 反向统一查询mask允许读取圆外其他可见TLS；前向排除全部已知TLS核心。两方向证据不完全对应同一背景，不强归因于目标TLS独有晕，也不把中心标签分类改称位置后验或桥检出。
- 完成结果与全量表见两份20261003结果报告；五张Python图及源数据索引已导出SVG/PDF/600dpi TIFF并视觉核查。R-18仍需晕/边背景对照、全GO证据整合及全目标验收，plan假设不因局部结果擅改。
- 失效条件：若后续获得可追溯患者映射、独立验证或新的适用空间检验，可按新合同升级相应证据；本条的内部结果不能追认成独立确认。

## D-179｜2026-10-03｜全量晕与边的原位几何背景及适用性审计
- 保持已完成的前向距离分箱与所有1691 GO、51固定对照，不按反向预测或单个GO效果选择距离范围。晕对全部已知TLS的最近距离计算；边只用既有明确TUM/非TUM界面，GSE无肿瘤标注保持不可测。
- 原位背景将整套TLS点或边界点平移，保留原表达、真实组织区室标签与真实TLS核心遮蔽。平移基向量为六个100μm六角格方向，以及六个长度约173.2μm的格方向；倍数1/2/3/5/8/13/20，共84个候选，覆盖近至远的位移。它们是合成几何参考，不新造TLS或肿瘤真值。
- 每个平移至少90%的锚点须落在现有组织点51μm以内。每个原本可测分箱，参考仍需≥30点且点数为真实箱的0.5–2倍；逐箱保留所有合格与不合格原因。晕参考另遮蔽平移后TLS及同宽隔离带；真实TLS核心始终遮蔽。边参考不改变TUM/非TUM标签，只改变到几何边界的距离。
- 效应包括全部相邻径向箱对比、边界两侧各自的相邻距离箱对比，以及相同距离箱的TUM−明确非TUM对比。标准差只用对应位置和mask的可见域。对比所需的任一箱不足支持即不可测，不越过空箱连接非相邻区间。
- 为共同选择诊断，使用对全部已可测对比都满足支持条件的同一平移集合；至少19个平移才有约5%尾部分辨率。每个位置都执行相同的距离、mask及正负方向最大值搜索，并逐GO与1691项联合最大值保留诊断秩。参考不足仍完整输出逐箱/逐对比背景；秩不当正式p。
- 必要校准沿用D-177六类明确已知场，各400次，但针对本晕/边算子重新执行；不是重用桥或R-17的p。误报单侧95%上界≤0.08为设置通过，任何失败保留。通过必要门不自动验证真实联合GO场或患者独立性。
- 功效模板为100/250/500/1000/2000/4000μm指数晕；边分别用界面峰和带TUM/非TUM符号的界面变化。每个模板按几何算子最大真实对比定标，使注入对比为0.5或1.0个可见域SD，背景用前200次模拟。判定要求模板最大对比本身胜过所有参考位置的完整搜索，避免把其他位置误报计为检出。0.5SD检出率单侧95%下界≥0.8才支持对应条件的有用效应阴性解释。
- 平移不保证不规则组织上可交换；必要门失败或参考不足时，结论保持描述性。主分析不回归掉细胞组成，区室差异可有价值，但不能仅凭区室差异声称TLS特异机制。
- 失效条件：按效果挑平移、放宽失败设置的支持/误报/功效门槛、变更分箱后沿用原合同，均不合法；进一步方法属于新研究合同，不能改写本次完整扫掠结果。

## D-180｜2026-10-03｜全GO证据整合与展示规则
- 保留每个GO在全部来源、晕/边距离箱、桥长度/半宽/路径/mask/端点类型及反向操作点的完整数值与不可测状态。数值汇总只作描述，不生成一个混合各尺度的“确认分数”。桥在源TLS实例内先汇总，再汇总至切片，保持既定重复单位。
- 晕/边分别保留原始箱均值、相邻/界面对比、原位参考均值与分位、有效参考数及对应校准状态。桥补入原位参考、负向连续性和双侧不对称；参考不足不是零效应。GO压缩面板1691项与51固定对照有共同索引，不能以榜单替代全量。
- 探索性反向榜单按每个mask和信息臂分别给出：以两个来源留出方向中较低的切片平均AUC排序，只给达到该操作点完整可评切片覆盖的GO排名，其余仍列出未排名原因。0.5为随机排序参考。此榜单来自已查看的内部验证，不能叫新的独立确认；不跨mask挑一个最高值代表全程序。
- 代表空间图使用固定B轴，以及500/1000μm条件联合信息臂的各一个榜首GO，均明确是展示而非新验证。连线按几何资格选取每张片、每类端点中mask130/半宽250/组织图路径的中位长度候选，不按表达选路径或最亮连线；无合格几何时不补线。全部1691项另有完整前向曲线图集与数值索引。
- 失效条件：把后验挑出的榜首/代表图当预先指定发现、隐藏未通过尺度或把描述性参考分位当p，均应撤回对应主张。整个目标仍以完整验收而非好看图例结束。

## D-181｜2026-10-03｜全量本地交付完成，正式确认不改判为成功
- 完成原冻结范围的全量评分、前向晕/边/三类桥、84个反向任务、适用性与功效审计、全程序整合及完整图表；执行验收见 `docs/TLS_MULTISCALE_ACCEPTANCE_20261003.md`。14项针对性测试通过，机器验收核对输入/缓存/预测及1742页图集。
- 原目标要求正式双向关系与桥确认；现有内部关联/预测、必要校准失败与功效不足不满足该目标，不将“文件齐全”替代科学确认，也不修改plan中的成功判据。长度与宽度跨尺度、组成信息保留的用户要求不变。
- 患者对应缺口在连续推进的本地计算、反向/桥校准、整合交付阶段始终未解除；最新26行patient_id全空、94个原始输入SHA不变。可执行的本地交付已完成，身份/独立真值无法由重复计算补出，正式确认记阻塞。空间方法缺口另列I-031，不把方法困难单独当作停止全部工作的理由。
- 已准备两既有来源官方患者对应元数据的≤10MB补证方案，尚待用户按AGENTS.md授权；未下载外部数据。选择先查已有来源的身份依据，不直接获取新表达队列；理由是它是患者级重复判据缺少的最小证据。若官方metadata仍没有交叉表，保持未知，不由样本命名或表达推断。
- 失效/重开条件：获得可追溯患者对应、用户明确确认独立性，或授权与新证据足以启动独立验证时重新评估；身份补齐不自动解除真实空间校准、功效和结构标注限制。进一步方法研究另立合同，不事后放宽本轮失败门槛。

## D-182｜2026-10-03｜执行已批准的有限身份元数据补查
- 用户本轮明确允许补查，授权限现有USZ Zenodo14620362和GSE175540及直接关联的身份metadata，总量≤10MB；不取表达矩阵或影像、不联系他人。先查官方记录与样本字段，保留URL、时间、大小、SHA与逐样本结论，并更新数据索引。
- 仅明确患者字段/交叉表可解除对应样本的I-033；样本序号、肿瘤数与队列总患者数不足以推断一一对应。若仍缺失，记录检索范围与失败原因，不扩展为未授权数据获取。

## D-183｜2026-10-03｜用户提供GSE175540补充表的显式ID匹配与重复组风险
- 表头含cohort、ID、pT及Spatial transcriptomics。以GEO官方Sample_title中的样本别名拆出Frozen/FFPE制备方式与原始ID，再与Excel的ID和制备字段联合核对；保留GSM、原始标题、Excel行号和来源SHA，不从本地文件顺序推身份。
- 全24个GSM对应23个表内ID；本轮18片对应17个ID。c_2同时列Frozen / FFPE，对应GSM5924030与GSM5924050，登记为同ID保守分组候选。正式patient_id暂不赋值，待原表图注明确ID单位；不把这两片默认当独立重复，旧逐片留出不能冒充按此ID分组的留出。
- 分期和检测开展标记可作为表内metadata；x/X与空白不是检测阳性/阴性，不能当TLS成熟度、丰度或治疗响应。来源为用户提供文件，原出版链接/图注待确认，非自行验证的官方附件。
- 失效条件：原表ID不是跨制备同一单位，或来源/ID规则出现冲突时撤回对应分组；确认ID为患者后才升级患者身份，并重新进行适用的分组留出。当前只登记与标记受影响结论，不静默覆盖已冻结结果。

## D-184｜2026-10-03｜采纳用户对GSE补充表患者ID的明确确认
- 用户在明确询问ID是否代表患者后答复“id确认”，据此将表内ID登记为患者分组；证据等级是用户确认，不冒充独立核验的出版社图注。全24片为23位患者，本轮18片为17位，c_2两制备片同组。
- 新建版本化患者交叉表与26源当前身份overlay，保留旧冻结身份与预测结果；GSE身份门解除，旧逐片留出未自动变为患者级验证。USZ8片仍未知。失效条件为用户修正或来源出现矛盾；后续训练/评估须使用同患者一起留出的分组并重新验收。

## D-185｜2026-10-03｜GSE患者分组验证，固定模型并复用等价折
- 用户要求继续推进。GSE以用户确认的17位患者分组，USZ继续按8个肿瘤样本，不认定患者独立或连续切片。全1691 GO、51对照、230/500/1000μm及三信息臂不变；不按新结果重新选尺度。
- 固定此前特征、岭模型、训练内top10选集及训练时切片/类别权重，以隔离验证分组的影响。16位单片患者的训练/测试集合与旧逐片折完全相同，复用其48个折；c_2两片同时留出，在三个尺度重新拟合，禁止两片分别出现在训练和测试。训练权重仍是原切片权重，非患者均衡训练，明确披露，不同时改变估计目标。
- 来源留出6个模型训练/测试未变化，复用完整预测；GSE按每片先算指标、同患者片均值、再患者等权汇总，USZ仍样本等权。所有GO保留，单类别/覆盖不足按原规则不可测；每项保留有效患者数与片数。患者内部分片不可测时使用可评片并标记，另列全片可评支持，不把缺测补零。
- 原模型与冻结成果不覆盖，复用折核对预测SHA及精确训练/测试集合。新的分组结果仅内部验证，不是独立盲验或空间校准/定位确认；不重排榜首充作新发现。
- 失效条件：身份表修订、特征/训练方法/选集算法变化或复用折集合不一致时必须重新计算相应折；本次复用结论仅在上述输入不变时成立。该步骤检验H-01内部信号是否受同患者重复影响，并收紧H-05的证据单位；不单独确认H-08桥。

## D-186｜2026-10-03｜前向患者汇总与c_2成对诊断
- 用户要求继续。保持全1742项、所有距离/路径/半宽/mask及原可测规则，原始评分/几何/空间参考不变。晕、边、桥先保留原切片值，再同一患者同条件的可测片取均值，最后在患者间报告中位数和正负方向计数；USZ每个肿瘤样本单独作为描述单位。不得跨不同尺度或端点类型补齐缺测，不把没有候选几何的患者计为阴性。
- 同步输出各单元的贡献切片数、完整片支持标记与来源行号；原切片中位数和患者中位数并列。参考z/分位及沿程比例只作描述性汇总，不合成正式p或宣称空间校准通过。原桥先TLS实例后切片汇总保持不变，患者层只是新增一层。
- c_2诊断限定已冻结B轴及两项既有展示GO，检查原始barcode/标签/坐标、原始counts重算得分、查询资格和可见特征是否一致；分解1000μm固定B轴模型的径向/角向贡献，并保留230/500/1000μm支持。比较其他可见TLS和组织采样范围只作解释线索，不据此新选读数、删除负例、翻转方向或改mask。
- 不把同患者两种制备视为物理相邻或同位置切片，不从一位患者的差异推出制备方法效应。失效条件：身份修订、层级/尺度混合或缓存与源不一致；若发现计算错误，另行记录并重算受影响结果，不能用描述掩盖。服务H-01/H-05的重复单位与适用条件，H-08仍须独立空间证据。

## D-187｜2026-10-03｜冻结查询设计的TLS尺寸与组织支持适用范围
- 用户要求继续。仅核查既有230/500/1000μm查询设计，不改mask、网格、标签或每环/可见点数门槛，不重选程序。所有26片均登记，无肿瘤/患者身份的信息仍未知。
- 每个metadata TLS连通域记录点数、点坐标直径、以该域spot为中心的最小完整遮蔽半径（加130μm隔离带）、固定查询网格内可行中心数量及实际可评正查询；区分singleton、无网格中心、核心无法完整隐藏、可见组织不足、源标签缺口。理论直径上界2×(mask−130)仅为必要条件，不等同于完整可测资格。
- 患者/样本层保留总查询与正负支持、不同TLS连通域贡献数；GSE按已确认患者分组，USZ按肿瘤样本。跨片同患者的TLS连通域没有对应依据，不把相加的切片实例称为独立生物结构，不用查询点数冒充有效患者数。
- 各尺度完整列出，不用表现好的尺度覆盖失败尺度。该节点说明H-01/H-02可识别前提及H-05适用范围，不自行验证H-08桥。失效条件：源metadata或几何设计改变后应重算适用性，不能将当前门槛直接外推到另一设计。

## D-188｜2026-10-03｜共同查询上的500/1000μm配对内部比较
- 固定D-185的分组模型预测，不重训、不筛最好半径、不翻转失败方向。比较已保存的全1742项联合特征模型、训练内GO组合、几何/深度基线；径向/角向单臂预测未保存，本次不冒充覆盖这些单臂。
- 每张片以源spot_index交集配对两mask的合格查询，要求标签一致；共同正查询天然对应同片同TLS连通域，共同负查询也固定。以metadata资格匹配，不用预测值筛样本。GSE用患者整体留出和来源留出，USZ只报来源留出的肿瘤样本结果。
- 逐片在同一查询集合上计算两档AUC/AP，患者内对两档均可测的相同切片取均值，再患者等权；USZ样本等权。缺任一类别/程序训练资格的比较保留不可测，不能用不同支持人群的均值相减。原完整人群结果并列保留，共同子集不能替代全量。
- 同时报全量到共同子集的患者/切片/TLS实例/查询变化，查询不是独立重复。该比较减少测试对象组成变化，但两档输入范围、归一化可见域、训练标签与训练选集仍不同，因此不能称纯距离效应或因果衰减。用已保存float32预测同时重算两档，明确与历史float64指标的小数精度区别。
- 服务H-01/H-05的操作范围比较；失效条件为身份或模型SHA不一致、标签冲突或未配对即平均。仍不确认H-02位置概率或H-08桥。

**D-188存储精度补注：** 实际归档中单程序prediction为float32，geometry/combined为float64；本轮均按各自已保存精度配对重算，没有量化或更改基线/组合预测。

## D-189｜2026-10-03｜桥对照失效条件诊断，不更改检验
- 使用D-177已保存的48个零场模拟结果、1776个功效设置及全部26源几何资格，分来源/场类型报告误报，分端点类型/路径/长度/半宽报告功效。已有零场统计只保存全几何搜索最大值，不能事后把误报率归到某一桥宽或某个端点类型。
- 对每份40方向×400次最大统计量，逐方向轮换为目标计算相同诊断秩，比较实际方向与其余方向的误报频率；复核原方向k/400与原表一致。方向平均只是诊断，不是新的真实p。
- 另计算所有共同几何、40方向的线性端点晕调整算子在单位独立噪声下的方差（权重平方和），量化原方向和参考方向的算子差异；用显式稠密权重核对低秩公式。该几何诊断不替代真实相关场校准，亦不用于重标准化既有结果。
- 功效通过必须同时标明匹配零场设置和整源六类场必要门是否通过；选择过后看见的通过子集不作为新确认集。原失败门槛、模拟次数和全GO诊断不改。失效条件：重定义统计量/尺度后沿用旧校准，或把独立噪声方差当真实场方差，均应撤回相应解释。服务H-08可检验性，保留不可判定而非生物阴性。

**D-189确定性注入补充诊断：** 对原148个几何选定的功效模板计算不加噪声时的端点调整增益和同一完整搜索下的目标诊断秩，检查宽模板是否也激活参考方向。复用原模板定义，不新增效应选择，不把无噪声极限当真实数据功效或正式p。


## D-190｜2026-10-06｜用户授权转向空间分子场条件恢复
- **背景：** 用户要求解包 `docs/013_spatial_field_completion_20261006.zip` 并依据包内内容更新项目上层文件。来源为本次明确指令及 [科学方案](field_completion_20261006/01_SCIENCE.md)、[路线](field_completion_20261006/02_ROADMAP.md)；实际 HEAD 为 `5a2d97aa8d1ef3315b32773dabf89f6d717b3ce7`，已跟踪工作区原先干净。
- **理由／选择：** 现行目标改为学习场间条件空间关系，从可见场任意可用子集恢复未观测场的身份、分布及不确定性。使用 HF-01～HF-05 与 FC-00～FC-05；TLS 是有独立结构参照的首个正案例，easy / harder 为核心，extreme / 视野外为探索。替换旧入口而非附加竞争主线。
- **代价：** 旧中心分类、外围描述及桥计算不等于新整场恢复结果；FC 实验节点本轮均未运行。原 plan、roadmap 及其他入口保存在 [实际 HEAD 的快照](archive/pre_field_completion_5a2d97aa/README.md)，历史 H/R/D 含义保留。
- **复查条件：** 若用户修订科学目标，或实际任务只剩同点共表达／直接复制，应重新说明科学范围，不能沿用“空间条件恢复”的贡献主张。
- **影响：** 重写 plan、roadmap、STATUS 与 TODO；沿用可审计资产和已有失败。此次目标变化已获授权；本轮仅文档迁移与接手准备，无新科学结论。

## D-191｜2026-10-06｜固定字典先检验恢复，分子组成发现随后接入
- **背景：** 新目标允许固定代理字典与学习字典互换；旧独立场发现与 K 搜索未证明全局唯一分解。
- **理由／选择：** FC-00→FC-01→FC-02 为主路径，先用已有 GO/Hallmark/B-axis 与独立 TLS 参照构建完整通道缺失任务。FC-01 有基本复核的恢复信号后即可开启 FC-03，不等待完整 FC-02；TLS 弱标签扩样并行，FC-04 随结果整理轻量关系表。排除先求“正确 K”、完成全库发现或大规模扩样才恢复的方案。
- **代价：** 固定字典只是代理读数，不能声称独立生物场或库外身份发现；发现阶段需要跨独立单位的身份对应与共同任务比较。旧 R-04 和 K 搜索维持停止，FC-03 不自动恢复旧运行。
- **复查条件：** 若固定字典失败有可复核的字典局限证据，可做有限学习字典探针；若新组合不稳定或不改善共同任务，保留固定字典／部分候选，不无限调 K。
- **影响：** FC-01/02 不被场发现阻塞；具体可逆实现由执行者选择并记录影响，不增加逐项模型审批。GPU 禁用、外部获取审批及 1.2 TB 存储余量继续适用。

## D-192｜2026-10-06｜分层证据并收窄旧空间校准阻塞范围
- **背景：** 旧旋转／场零模型已有误报与功效失败；GSE 18片对应17位患者并完成分组内部验证（D-184/185），USZ 8个肿瘤样本的患者对应仍未知。新问题以合法留出上的恢复比较为主要证据。
- **理由／选择：** I-031 继续硬阻塞具体桥及相应空间显著性／普遍不存在主张，不再前置于全部患者级恢复探索。预测增益、经验不确定性校准、结构显著性与因果解释分别评价。自动标签可带来源、置信度和版本用于弱监督／探索；独立病理或正交参照承担结构测试真值，不通过 RNA 签名或模型一致性自证。
- **代价：** 退出主路径不等于校准通过；旧中心 AUC 不升级整场定位，已反复用于开发的来源不重新称盲验。USZ 身份未知仍限制其患者级主张；短窄桥、稀疏边界与无视野外真值按具体任务保留限制。
- **复查条件：** 若主张涉及正式桥 p、因果作用、库外发现、三维／视野外精确定位，须提供对应独立证据；若测试输入读到目标原始读数、全图 latent 或循环弱标签，相应恢复结论不可用，先纠正信息边界。
- **影响：** ISSUES 按来源与主张更新范围；新目标风险先作为设计关注记录，不虚构已发生故障或新增统一审批门。旧失败、不可测和未知记录继续保留。

## D-193｜2026-10-06｜FC-00 固定首批面板、目标与可见输入
- **背景：** FC-00 核对本地 26 片表达、坐标、标签和既有评分缓存。GSE 18 片对应已确认 17 位患者，USZ 仍是 8 个患者未知的肿瘤样本。评分缓存为逐点 AUCell，没有跨点平滑。
- **理由／选择：** 首批固定字典用已冻结的 50 个 Hallmark 加 B-axis，不把全部 1691 个 GO 同时作为输入。首个结构目标是独立 TLS 标注；首个非 TLS 分子目标是氧化磷酸化 Hallmark，依据是它与 B-axis 无共有基因、不是淋巴细胞签名，且 26 片覆盖合格。B-axis 只作为 TLS 相关的次要分子目标。选择依据是字典身份和基因重叠，不是测试集表现。
- **代价：** 51 项面板可能信息不足，不能代表全部 GO。氧化磷酸化仍有 100/200 个基因出现在其他 Hallmark 中，通道隐藏不等于基因隔离，必须另做去重叠对照。b_18 没有明确阴性，USZ 不能作患者级确认。
- **复查条件：** 若空间增益只来自重叠基因，相应关系不能支持 HF-01。若 51 项面板失败且训练折显示是字典覆盖不足，只允许在训练信息内做一轮有限修正，不能扩大成全 GO 扫描或按测试图重选目标。身份表或缓存与原始计数不一致时，受影响结论撤回。
- **影响：** FC-01 按此输入边界开跑。摘要见 [FC-00 资格](FC00_ELIGIBILITY_20261006.md)。此选择不改变 plan，也不把旧结果改称为新恢复。

## D-194｜2026-10-06｜冻结 FC-01 留出、基线与继续规则
- **背景：** D-193 已固定面板和目标。FC-01 需要在看到留出结果之前写死模型、划分和“什么算空间增益”。
- **理由／选择：** GSE 按 17 位患者留出，USZ 按肿瘤样本单独留出，两来源不混训。基线为训练均值、同点线性模型、同点梯度提升，以及仅用于 TLS 的 B-axis 逻辑回归。空间模型使用同一梯度提升，只多两圈可见通道的邻域均值。氧化磷酸化和 B-axis 必须同时报告去掉共有基因程序的对照。打开 FC-02 的操作规则是患者方向一致且达到预先写下的实用下限；这不是显著性门槛。没有增益时只允许合同中的一轮修正。
- **代价：** 圈层、惩罚和迭代是可逆的工程选择，可能不是最优模型。实用下限可能放过很小的增益，也可能挡住不稳定的小数差异；两种情况都要报告原始数值，不能只报是否开门。
- **复查条件：** 若在看到测试结果后改目标、圈层、划分或增益规则，相应比较撤回。计时探针的分数不得当作评价结果。若正式运行前只因墙钟超过 12 小时而把迭代降到 30，必须先追加说明，不能事后用探针分数挑选。
- **影响：** 合同见 [FC-01 合同](FC01_CONTRACT_20261006.md)。尚未产生恢复结果，不改变 HF-01 或 HF-04 的证据。

## D-195｜2026-10-06｜计时探针未触发迭代下调
- **背景：** D-194 允许在正式比较前用一个患者折估计墙钟。若外推超过 12 小时，只能把梯度提升迭代从 50 降到 30，并在开跑前记录。
- **理由／选择：** 探针完成，单折各模型拟合均在 1 秒内。17 个患者折加 8 个肿瘤样本折的外推远低于 12 小时，因此正式比较保持迭代 50。探针输出单独放在 `probe_not_evaluation/`，不进入评价，也不用来改目标、圈层或增益规则。
- **代价：** 短探针不能证明全部折都同样快；若正式运行中途触及存储下限，仍须停止。
- **复查条件：** 若事后用探针分数选择模型或阈值，相应正式比较撤回。
- **影响：** 正式留出比较按原合同开跑。此条不提供 HF-01 或 HF-04 的证据。

## D-196｜2026-10-06｜FC-01 未达门槛，只做预先规定的容量修正
- **背景：** 正式留出比较完成，泄漏核对通过。GSE 上 TLS 的患者平均 AUC 增益为负，氧化磷酸化去重叠对照的 RMSE 相对降低约 1.9%，低于预先写下的 2% 下限。两者都不打开 FC-02。
- **理由／选择：** 同点模型在留出患者上已经明显好于训练均值：TLS 的 AUC 约 0.85 对 0.50，氧化磷酸化去重叠对照的训练 RMSE 比先验低约 40%。这不是“训练接近先验”的字典不足，因此不加入 GO，不改目标或圈层。唯一允许的修正是把梯度提升迭代从 50 改为 120，其余合同不变。不因 1.9% 接近 2% 而降低门槛。
- **代价：** 多迭代可能只减少训练误差，不产生空间增量。若修正后仍不达标，FC-02 与 FC-03 关闭。
- **复查条件：** 若用测试图重选目标、改圈层或改下限，相应结论撤回。
- **影响：** 修正结果单独存放，不覆盖原 50 次迭代的正式比较。此条仍不是 HF-01 的成立证据。

## D-197｜2026-10-06｜容量修正后仅氧化磷酸化去重叠臂达到操作门槛
- **背景：** 50 次迭代的正式比较未达 D-194 的继续规则。D-196 的唯一修正把迭代改为 120，其余合同不变。
- **理由／选择：** 120 次迭代后，GSE 17 位患者上氧化磷酸化去重叠臂的 RMSE 相对同点梯度提升平均降低 4.1%，12/17 位同向，达到预先的 2% 下限和半数方向规则。TLS 的平均 AUC 差仍约 +0.001，不达标。含共有基因的输入上空间模型仍不优于同点模型。按冻结合同，这一结果打开 FC-02，但不把 HF-01 写成已经成立，也不把 B-axis 或 TLS 写成有空间增量。
- **代价：** 均值受少数患者影响，b_18 空间模型明显更差。操作门槛不是显著性或校准。USZ 仍不是患者级证据。
- **复查条件：** 若事后放宽 2% 下限、删掉 b_18，或用含重叠基因的臂声称空间增益，相应继续决定撤回。
- **影响：** 进入有限的输入组合比较。FC-03 可以做有限探针，但不能重开 K 搜索。见 [结果](FC01_RESULTS_20261006.md)。

## D-198｜2026-10-06｜FC-02 只比较三条预先规定的输入阶梯
- **背景：** D-197 打开输入组合比较。不能在看过子集表现后再挑选“最小有用输入”。
- **理由／选择：** 沿用 120 次迭代、原留出和原邻域。氧化磷酸化只在无共有基因的程序里，按训练折绝对相关取 1 个、5 个和全部。TLS 的单输入固定为 B-axis，五输入是 B-axis 加训练折相关最高的四个其他程序，再与全部 51 项比较。组合增益沿用 2% RMSE 或 0.02 AUC，且至少一半患者同向。不枚举幂集，不扫描 GO。
- **代价：** 训练相关选出的程序可能不是生物上互补的组合。阈值仍是操作规则，不是机制证明。
- **复查条件：** 若用测试目标挑选阶梯，或把未达标的更大集合写成高阶机制，该比较撤回。
- **影响：** 合同见 [FC-02](FC02_CONTRACT_20261006.md)。比较完成前不声称组合增益。

## D-199｜2026-10-06｜FC-03 只用固定 K=8 的有限字典探针
- **背景：** D-197 给出一个狭窄的初步信号。路线允许此后做分子组成探针，但不能把 K 搜索或旧 R-04 重开。
- **理由／选择：** 在看拟合结果前固定 K=8、非负矩阵分解和种子 0。氧化磷酸化探针去掉目标基因。与固定面板使用同一留出和同一读出。不稳定因子只称候选。不扫描其他 K，不使用全部转录组做无界发现。
- **代价：** 8 个因子可能太少或不可对齐。探针失败只限制学习字典，不推翻固定面板上已经记录的同点关系。
- **复查条件：** 若比较多个 K 后挑选最好的，或用因子自重建宣称优于固定字典，该探针撤回。
- **影响：** 合同见 [FC-03](FC03_CONTRACT_20261006.md)。探针完成前不声称发现了新的分子场。

## D-200｜2026-10-06｜FC-03 探针的归一化与不稳定阈值
- **背景：** D-199 固定了 K=8，但没有写表达尺度和“不稳定”的操作界。探针尚未运行。
- **理由／选择：** 使用全部基因表达特征的文库大小归一到一万后再取对数，只保留 51 项面板里的基因。矩阵分解用随机初始化、种子 0、最多 100 次迭代。训练患者对半拟合后，载荷配对的中位绝对相关低于 0.5 视为不稳定，只能称候选。
- **代价：** 这一尺度和阈值是可逆的工程选择，可能让因子更平滑或更容易被判成不稳定。
- **复查条件：** 若在看过探针结果后改归一化、K 或 0.5 阈值，相应比较撤回。
- **影响：** 这些细节补入 [FC-03 合同](FC03_CONTRACT_20261006.md)。仍不授权扫描其他 K。

## D-203｜2026-10-06｜启动图像 TLS 弱标签探针
- **背景：** 用户批准用非 RNA 方式补充 TLS 位置标注。本地已有 YOLO 两步权重（`repo/data/TLS_classification_HE_YOLO_step*.pt`，泛癌 TLS 管线）与 GSE 每片配对组织图，无需下载。ultralytics 本地 pip 安装，CPU 运行。
- **理由／选择：** 先在 1 片上做计时探针，再核对检测框与测序点配准、与独立标注的偏差。标签只作弱监督／扩样，不升级为确认性真值；未标注区不补阴性。YOLO 训练集与考试样本的重叠核查在扩样前完成。
- **代价：** CPU 推理可能慢；若超时则报阻塞，不租 GPU。USZ 本地无原图，图像路线覆盖不到。
- **复查条件：** 若用弱标签作关键测试真值、或凭模型一致升级标签，该结论撤回。
- **影响：** 探针输出单独存放，不污染现有独立标注。

## D-204｜2026-10-06｜用户批准下载外部TLS定位权重并外包算法调研
- **背景：** 本地YOLO只是切块分类器，不能定位；HookNet-TLS等定位链缺权重。用户明确批准下载外部权重，并要求researcher子代理调研热门、可靠、基准优秀的TLS预测定位算法，找一个公认方法先试；A路切块分类放到最后。
- **理由／选择：** 只下载researcher推荐的定位模型权重，不批量下载外部数据集。下载后记录来源、版本、许可证与SHA-256，再做训练重叠核查。GPU仍禁用。
- **代价：** 外部权重可能与考试样本重叠；重叠部分必须隔离，不能当独立验证。
- **复查条件：** 若下载整套外部数据、租GPU或把弱标签升级为真值，该批准不覆盖，需另行确认。
- **影响：** researcher子代理已启动，只做调研不下载；A路暂缓。

## D-205｜2026-10-06｜researcher子代理启动失败，未产生结果
- **背景：** 已授权的子代理调研在启动时失败，未开始实质工作。
- **理由／选择：** 失败原因是该 researcher 角色请求的 `fetch_content`、`get_search_content`、`source_check` 在当前子代理工具白名单中不可用；属同一协议的通道故障，不是研究结论。未改用其他隐蔽通道继续。
- **代价：** 算法调研暂无结果；A路仍暂缓，B路权重也未下载。
- **复查条件：** 需用户明确选择后才可继续：由我直接做调研，或修正工具配置后重派 researcher。
- **影响：** 等待用户指示。

## D-206｜2026-10-06｜TLS定位算法首选HookNet-TLS
- **背景：** 直接调研（researcher通道故障后由执行者直接做）：HookNet-TLS（Communications Medicine 2023）在约1019张TCGA切片（肺鳞癌、膀胱癌、肾透明细胞癌）上穷举标注训练，与Faster R-CNN和阅片者间差异做过基准比较，并用IHC引导标注验证过；权重在Zenodo公开（10614942），代码在GitHub，Grand Challenge有可部署版本；输出TLS加生发中心分割掩膜，可配准到测序点。备选EfficientNet分割（Dice 0.91，内部集）权重可得性不明，TLScope太新。
- **理由／选择：** 公认程度、基准完备性、权重可下载、输出可配准四项都满足，且与本仓库已有管线引用一致。先试HookNet-TLS。
- **代价：** 训练集含大量肾透明细胞癌，与GSE同癌种；论文还用过苏黎世USZ作独立测试集，与我方USZ队列可能重叠。因此USZ不能作为HookNet标签的独立验证，只能用GSE或另行隔离；重叠核查在扩样前完成。需Docker与openslide，CPU推理先计时。
- **复查条件：** 若权重无法下载、CPU超时或重叠核查发现考试样本在训练集中，改用A路或另选算法。
- **影响：** 按D-204下载权重，记录版本与哈希。弱标签仍不升级为真值。

## D-207｜2026-10-06｜训练重叠核查结论与A路并行启动
- **背景：** HookNet-TLS训练于TCGA（LUSC/BLCA/KIRC），另用USZ档案15例（7 KIRC、8 LUSC）做基准。我方GSE为Meylan独立队列，与TCGA无共有切片证据，但同癌种，记同癌种迁移。
- **理由／选择：** 我方USZ（Zenodo 14620362，2025年发表，3肾5肺Visium配对）与论文USZ集（15例，IHC连续切片）研究不同、数量不同，直接重叠可能性低，但同医院同癌种，无法证伪独立性。因此USZ隔离，不承担HookNet标签的独立验证；只用GSE做偏差检查。用户批准滞后任务并行，A路切块分类探针与HookNet环境搭建并行启动。
- **代价：** A路切块误报需偏差检查把关；两路并行占CPU，已限8线程。
- **复查条件：** 若发现考试切片确在任一训练集中，相应标签作废并隔离。
- **影响：** 弱标签仍只作训练扩展。

## D-208｜2026-10-06｜A路切块探针失败：YOLO是分型器不是检测器
- **背景：** c_57整片1083点逐点切块分类，CPU约251秒。
- **理由／选择：** 全部判Mature，置信度TLS全1.0、非TLS中位0.999，无鉴别力。该YOLO是候选patch成熟度分型器，密铺超出其设计。此形态判失败，不扩样。
- **代价：** 探针成本约4分钟CPU。
- **复查条件：** 换权重或换用法后重探需另立决策。
- **影响：** A路此形态关闭。

## D-209｜2026-10-06｜B路HookNet探针失败：本地H&E分辨率不够
- **背景：** 权重已验哈希，独立conda环境（TF 2.13）跑通，两种分辨率假设下输出通道1全片饱和、TLS与阴性无区分、无空间结构。
- **理由／选择：** 本地GSE组织图是约5µm/px的hires PNG，上采样无法恢复模型所需的40×级别细胞核纹理；Space Ranger产物不含真全分辨率图，GEO也下不到。因此不是模型坏，是输入端先天不足。在Visium队列上判B路失败，不跑全片。
- **代价：** 已花时间建环境与探针；环境保留但不继续投入。
- **复查条件：** 若未来拿到真高分辨率配对图像，另立决策重探。
- **影响：** 图像弱标签扩样两条路都关闭；弱标签分支回暂缓。

## D-210｜2026-10-06｜外部TA-TLS图谱adata_tls_full.h5ad验看结论
- **背景：** 用户提供figshare私链并授权下载单个文件。经浏览器下载`adata_tls_full.h5ad`（606MB，CC BY 4.0），只读验看，未改写。
- **理由／选择：** 该文件23641点×38732基因，20个癌种（含KIRC 3982点），带成熟度（LA/PFL/SFL）、细胞富集类、TLS内距离字段；但点已预选为TLS富集区、无空间坐标（obsm仅UMAP）、病人ID为1～89数字编号无对应表。因此不能直接扩充整场恢复队列（缺整片上下文、缺阴性、缺坐标），也不作独立验证。可用作成熟度参照与组成参照。真正可能有用的是该数据集两个原始队列文件夹，需另行评估与归属核查。
- **代价：** 已落盘606MB，存储仍高于1.2TB下限。KIRC点的来源重叠尚未核查，在用作任何训练前必须先查。
- **复查条件：** 若将其点直接当整片恢复输入或当独立验证，该结论撤回。
- **影响：** 本文件记为参照级外部数据，不进入恢复队列。

## D-212｜2026-10-06｜按引用ST流程扩充：先监督5片，再探弱标签
- **背景：** 引用方法（ST签名打分→平滑→k-means→连通域→成熟度分层）可在本地复用。6个闲置GSM身份均为用户已确认病人，其中5片有作者TLS标注（034无）。
- **理由／选择：** 先做监督扩充（5片进核心，不做弱标签冒充）；再在有独立标注片上探ST检测法的灵敏度/特异度，通过才用于无标注队列；成熟度E/P/S作为新标签维度。引用iStar像素插值不做（需训练神经网络，超出现场资源），用点级分数加邻域平滑代替并如实记录偏差。
- **代价：** 034无标注暂不进；ST检测标签永不进关键测试。
- **复查条件：** 若把ST检测标签当独立真值，或身份表修订，相应扩充撤回。
- **影响：** 核心26→31片（病人17→22），见资格补遗。

## D-214｜2026-10-06｜闲置6片重定级：2片直扩，3片隔离，1片无标注
- **背景：** 6个闲置GSM身份虽确认，但标注模式不同：c_21/b_13用TLS_2_cat（与核心一致）；a_1/a_17/b_7只用TLS列标T_agg阳性、无阴性、其余空白；c_10无标注文件。
- **理由／选择：** 只有c_21/b_13进监督核心（26→28片，17→19病人）。T_agg语义不明且无阴性，不当TLS用，按未知处理，其表达坐标可留作无监督用途；c_10同样无标签。T_agg与TLS_2_cat不互认、不合并。
- **代价：** 监督扩充比预期小；成熟度只在TLS_2_cat实例上做。
- **复查条件：** 若查清T_agg与TLS对应关系并经独立核验，可另立决策升级；空白永不补阴性。
- **影响：** D-212的数量目标按本条修正。

## D-213｜2026-10-06｜ST检测探针失败，只做监督扩充加成熟度
- **背景：** 在c_57/c_2上有独立标注的两片上探引用ST流程（点级签名均值→邻域平滑→k-means取最高簇→连通域），灵敏度0.23/0.00，特异度0.97/0.92；成熟度分组在仅2个域上退化。
- **理由／选择：** 该形态漏检过高，用作训练会把真TLS教成阴性，判失败，不用于任何队列。扩充改为监督5片（身份与作者标注现成）加全量独立TLS实例的成熟度E/P/S标注（边界来自独立标注，只在域内算表达分型，不循环）。
- **代价：** 无标注队列仍无TLS标签；成熟度是表达分型，只作新标签维度。
- **复查条件：** 若把成熟度分型反过来论证TLS存在，该部分撤回。
- **影响：** 写扩展评分与成熟度脚本，结果单独存放。

## D-211｜2026-10-06｜新数据文件夹无法自助下载，转作者索取
- **背景：** 6例新数据文件夹成员旧ID全部404，文件夹打包URL服务器端502；12个大文件的新ID有效但只是处理好的子集。跨命令引用页面元素曾导致一次误下重复文件，已删除。
- **理由／选择：** 停止反复撞下载墙。回map验证与重叠核查都依赖原始文件夹（含条码、坐标、图像），在拿到之前两道前置都不算通过，不启动全部病人下载。
- **代价：** 扩容暂停，等待作者回复。
- **复查条件：** 拿到有效链接或样本元数据后恢复验证。
- **影响：** 起草给作者的索取邮件待用户确认发出。

## D-201｜2026-10-06｜FC-02 保留非单调的输入组合，不写成机制
- **背景：** 预先规定的 1 个、5 个和全部输入比较已经完成。选择发生在训练折内。
- **理由／选择：** 氧化磷酸化的全部去重叠输入优于 1 个和 5 个，但 5 个不如 1 个，因此不把前五个相关程序当成有用组合。TLS 上五个优于只有 B-axis，全部 51 项反而不如这两个较小集合，因此保留较小输入。两者都不命名为高阶机制。USZ 只作样本描述。
- **代价：** 操作下限会把小而同向的组合增益记为有用，也会忽略未达半数方向的平均差异。
- **复查条件：** 若事后删掉 b_18、改阶梯或把全套 TLS 输入说成优于 B-axis，相应结论撤回。
- **影响：** 结果见 [FC-02](FC02_RESULTS_20261006.md)。HF-02 只获得有限、非单调的支持；HF-04 仍无校准或患者级外部证据。

## D-202｜2026-10-06｜FC-03 只部分保留学习因子，不称发现
- **背景：** 固定 K=8、种子 0 的探针已完成，34 个折×目标的载荷稳定性都高于预先的 0.5 界。
- **理由／选择：** 氧化磷酸化上学习因子优于固定去重叠面板，达到预先下限；TLS 上不如固定面板。只保留氧化磷酸化方向的部分候选，不扫描其他 K。学习因子接触了固定面板之外的基因，增益不能直接解释成更好的分子身份。
- **代价：** 信息预算不完全相等，收敛上限是已知限制。
- **复查条件：** 若用其他 K 挑选最好结果，或把因子自重建写成胜利，该探针撤回。
- **影响：** 结果见 [FC-03](FC03_RESULTS_20261006.md)。HF-03 只获得一个目标上的部分候选证据。

## D-215｜2026-10-06｜主线不按癌种划分，起步队列按“独立参照×可用量”选
- **背景：** 用户指出肾癌不是主线，主线不区分癌种，只是训练模型；应选可用样本最多的癌种先行，后期扩多癌种。核查确认 plan.md 从未定肾癌为主线（通用框架，TLS 只是首个正案例）；肾癌起步是执行惯性（Meylan 有最好的逐点 TLS 标注），USZ 本身就是肾 3＋肺 5，从不限定单癌种。
- **理由／选择：** 训练起步队列的选择规则定为“有独立参照前提下的最大可用量”，当前恰好是肾癌，不是因为肾癌特殊。按纯样本量排，乳腺是第一候选（HEST 140、10x 与 STOmicsDB 均以乳腺居多），但乳腺缺独立 TLS 参照，只能先做无监督／弱标签训练，不能做法验证；参照缺口补上前不进验证。跨癌种主张仍按 HF-04 逐项凭证据，不 blanket 推广。
- **代价：** 现有恢复结论仍只覆盖肾为主的内部队列；乳腺扩量前需先过身份去重与平台差异审计。
- **复查条件：** 若把某癌种结果直接推广到全癌种，或把无参照队列当验证，该结论撤回。
- **影响：** STATUS 中“肾癌切片”措辞改为中性表述；后续扩量优先乳腺，但以参照到位为前提。

## D-216｜2026-10-06｜四癌扩充可行性：ST全覆盖，视觉只看真全片
- **背景：** 用户要求按引用方法在肝/乳腺/肠/肺扩样本，并判定视觉标注可行。核查：HEST四癌Visium的h5ad与wsi全在本地（肠78/肺18/乳腺11/肝5）；10x四癌43条目目录在（部分spatial包未解）；HTAN49、ST_CRC14、GSE23标注俱在。真全片如TENX13达24240像素边长，GSE小图仅2000像素边长。
- **理由／选择：** ST扩充四癌全可行（沿用冻结面板打分），但输出只是候选弱标签，先在28片核心上做灵敏度/特异度审计，不合格不用。视觉标注只看真全片：人工可行；HookNet自动标注需在HEST真全片上重探一次（GSE小图已判失败，不代表全片）。弱标签永不进关键测试。
- **代价：** 引用iStar像素插值与病理复核环节本地没有，ST版是降级复刻，只能当候选。解包与WSI落盘前重估存储（现余1.6TB，下限1.2TB）。
- **复查条件：** 审计不达标而强行使用者，相应标签作废；USZ身份未知仍限患者级主张。
- **影响：** 优先乳腺人工标注试点＋四癌ST打分；HookNet在真全片上重探一次。

## D-217｜2026-10-06｜四癌打分完成，弱标签审计因方法已死而终止
- **背景：** HEST四癌99片冻结面板打分完成（肠78/肺18/乳腺11/肝5，中位有效项1642）。其中TENX28/TENX23是千基因级靶向小panel，覆盖约0，已标不可用。
- **理由／选择：** D-213已判ST检测探针失败（灵敏度0.23/0.00），28片审计的前提不存在，不再为产生TLS标签而审计；97片有效分数转作训练背景扩展（无监督输入），不作标签。CSC/转置/小panel三类格式问题已修并记录，manifest合并覆盖事故已从npz重建。
- **代价：** 四癌仍零TLS标签；标签缺口只能靠人工/半自动视觉标注补。
- **复查条件：** 若把这97片分数当TLS标签用，相应结论撤回。
- **影响：** 待乳腺人工标注试点与HookNet真全片重探。

## D-218｜2026-10-06｜D-213死刑越权，改全量28片审计后再判
- **背景：** 用户指出D-213仅凭2片肾切片判ST检测法死刑，换大队列应重测。核查属实：探针只跑了c_57/c_2，方法失败模式（只抓最亮核心）是否跨癌种成立并未检验；且新癌种虽无标签无法直接验，但28片核心（含USZ肺）足以做更强审计。
- **理由／选择：** D-213降格为“2片探针失败”，不作为方法死刑。在全部20片GSE＋8片USZ上重跑审计：GSE报灵敏度/特异度，USZ报TLS覆盖率＋TUM+NOR误报率（UNASSIGNED/空白永不参评）。通过标准事前冻结：GSE中位灵敏度≥0.5且中位特异度≥0.9，否则方法仍判失败。
- **代价：** 即使通过，新癌种标签仍是未校准弱标签，永不进关键测试。
- **复查条件：** 审计不达标而使用者，标签作废。
- **影响：** D-217中“审计终止”改为“审计扩大到28片”；D-213原文保留，仅限定适用范围。

## D-219｜2026-10-06｜28片审计确认ST检测法不达标，关闭该线
- **背景：** D-218要求20片GSE＋8片USZ重审，冻结通过线为中位灵敏度≥0.5且中位特异度≥0.9。
- **理由／选择：** 实测GSE中位灵敏度0.22（20片中4片为0，最高1.0，异质性极大）、中位特异度0.95；USZ覆盖率0~0.52，两片在肿瘤/正常区误报超0.1。灵敏度远不达标，方法判失败——这次是28片的结论，不再是2片探针的问题。失败模式与D-213一致（只抓最亮核心），证实是方法内禀局限，非肾癌运气差。
- **代价：** ST弱标签线彻底关闭；四癌97片分数只作训练背景。
- **复查条件：** 换方法（如完整iStar＋病理复核）后可另立决策重探，本 maladapted 流程不再启用。
- **影响：** 标签缺口只剩视觉标注一条路（乳腺人工试点＋HookNet真全片重探）。

## D-220｜2026-10-07｜iStar能跑，但只在真全片上值得跑
- **背景：** iStar代码公开（daviddaiweizhang/istar），逐样本训练、无需预训练权重，输入（H&E＋counts＋坐标＋像素尺寸＋组织掩膜）本地齐全；官方称GPU高度推荐。实测GSE hires约3.3µm/px（离0.5µm/px差约6.6倍），HEST真全片约0.27µm/px（达标且有富余）。
- **理由／选择：** 能跑≠跑了有用。GSE小图缺细胞核纹理，学出的组织学→表达映射近似平滑，超分增益存疑；只批准一片HEST＋一片GSE对照试跑，以留出点重建误差和TLS分数对比验收，不达标不铺开。推理只存签名基因，不存全转录组。
- **代价：** CPU训练单片分钟级，0.5µm/px全片推理单片数小时（后台跑）；clone外部代码需记录版本commit。
- **复查条件：** 对照试跑无增益则关闭iStar线；弱标签永不进关键测试。
- **影响：** 见TODO对照探针任务。

## D-221｜2026-10-07｜四癌iStar/HookNet可跑口径冻结
- **背景：** 用户要四癌可跑样本数。核查HEST全目录（非抽样清单）：本地h5ad+wsi双全肠78/肺18/乳腺11/肝5；像素尺寸以HEST_v1_3_0.csv为准。
- **理由／选择：** iStar口径＝有H&E＋counts＋坐标＋像素尺寸＋组织掩膜可推导，HEST本地112片全达标，GSE24片可跑但降级（约3.3µm/px）。HookNet口径＝像素尺寸≤0.6µm/px：肠75/肺18/乳腺11/肝2（仅2例HCC；3例健康肝约0.9不达标）。10x各集spatial包未解，解包验收后另计；USZ无tif，两者皆不可，先取图。
- **代价：** 口径是工程判定，不是生物学结论；Visium HD先按Visium计，差异另行审计。
- **复查条件：** 像素尺寸元数据被证伪，或用超标样本出标签，相应计数与结论撤回。
- **影响：** 后续探针与标注任务按此口径取数。

## D-222｜2026-10-07｜批准直连下载HIPT特征权重
- **背景：** iStar对照探针已备好HEST与GSE输入，特征提取需要HIPT的vit256/vit4k权重；用户明确批准下载并要求尽量不用代理。
- **理由／选择：** 只下载该探针所需的两个公开特征权重，不下载演示数据与整库权重；使用直连，显式绕过代理环境变量；校验文件大小与SHA-256后落盘并记录来源、版本与用途边界。
- **代价：** 外部权重仅作图像特征提取器，不进入标签真值层；若下载中断或校验失败，先停后报，不继续训练。
- **复查条件：** 若权重被用于确认性真值、因果主张或超出对照探针范围，相关结论撤回。
- **影响：** 权重就绪后启动一片HEST与一片GSE的400轮CPU对照训练与推理。

## D-223｜2026-10-07｜HIPT直连下载失败，备选需用户定
- **背景：** D-222批准直连下载HIPT两个权重。实测`upenn.box.com`可达但二次跳转到`upenn.app.box.com`后连接超时（curl 28，131秒无响应）；GitHub Releases无.pth；HuggingFace API不可达。GitHub主站可达（慢）。
- **理由／选择：** 直连通道判定不可用，暂停下载。备选：①用户批准后仅对这两个文件走代理（校验SHA-256）；②用户线下提供权重；③跳过HIPT改用替代特征（偏离原文，可比性下降，不推荐默认走）。
- **代价：** 在用户拍板前，iStar对照探针停在输入就绪阶段；HookNet臂不受影响（weights.h5已在本地）。
- **复查条件：** 任一备选获批后恢复；若用替代特征，必须另立决策说明偏离。
- **影响：** todo #21保持in_progress（阻塞），不标记完成。

## D-224｜2026-10-07｜用户特批走代理下载HIPT权重
- **背景：** D-223确认直连不可用后，用户批准仅对vit256/vit4k两个文件走代理。
- **理由／选择：** 代理仅用于这两个公开特征权重，落地后验SHA-256并记录哈希；不扩大到其他下载。
- **代价：** 与“尽量不用代理”的偏好冲突，但直连已证伪，这是代价最小的合规路径。
- **复查条件：** 哈希不符则删除重下，不进入训练。
- **影响：** 权重就绪后启动HEST＋GSE对照探针。

## D-225｜2026-10-07｜Box链路已死，改走Git LFS官方仓库
- **背景：** 代理下Box返回404（文件已不在），直连跳转主机超时；HIPT官方仓库用Git LFS托管同名检查点，GitHub可达。
- **理由／选择：** 改从MahmoodLab/HIPT仓库`git lfs pull --include "*.pth"`获取同名权重，落地后验大小与可用性（无官方哈希可对，以能加载为准）。仍仅限这两个特征权重。
- **代价：** 多一次clone；若LFS对象也被删则转入线下索取。
- **复查条件：** 文件名/大小与iStar期望不符时停用。
- **影响：** 探针输入就绪后，双片对照训练因主环境skimage二进制冲突失败（numpy 2.x vs 1.25编译的skimage），需换干净环境或隔离venv重跑；GPU租用估算已交付用户（见D-226，若获批）。

## D-226｜2026-10-07｜iStar的GPU租用估算（单片约5-15分钟，全量约半天）
- **背景：** 用户提议租V100。从iStar实际代码路径估算：HEST真全片缩到0.5µm/px后约13k像素边长、16个4096瓦片、约67万16px子令牌；训练5状态×400轮的小网络；推理只跑组织掩膜内约23万令牌、批次50、5状态取中位。
- **理由／选择：** 保守估计单片约5～15分钟（含特征提取、训练、推理和读写），GSE小图同量级；112片HEST全量约11～28小时，单卡V100一天内可完成；四癌探针4片约23～60分钟。CPU路径当前因环境二进制冲突未跑通，修复后估慢一个数量级。
- **代价：** 估算基于V100与官方代码路径，实际吞吐受磁盘IO和批次大小影响；租用决策与费用由用户定。
- **复查条件：** 若实际单片超过估算的三倍，重估后再铺开。
- **影响：** 若租用获批，优先跑对照探针2片验证增益，再决定是否全量。
- **实测更正（2026-10-07，2080 Ti）：** HEST_TENX13 墙钟约 32 分钟，不是 5–15 分钟。主因不是卡，见 D-229。

## D-229｜2026-10-07｜iStar 单片 32 分钟，慢在代码结构不在 2080 Ti
- **背景：** D-226 按 V100 估单片 5–15 分钟。实测 RTX 2080 Ti 上 HEST_TENX13 约 32 分钟。用户问为何比 V100 估算更慢。
- **理由／选择：** 不把这 32 分钟解释成卡慢。2080 Ti 与 V100 在 FP32 上同级（约 13.4 vs 14 TFLOPS），显存带宽 V100 更高（616 vs 900 GB/s）；本任务是 FP32，用不上 2080 Ti 的游戏光栅优势。时间拆开：CPU 缩放 4.6 分钟；特征 16 次位移各约 30 秒共 8 分钟（`get_embeddings_shift` 对整张图做 4×4 位移并每次重载 ViT，D-226 只算了 1 次 16 瓦片）；小网络训练 5×186 秒共 15.5 分钟（30 batch/epoch、62 it/s，是启动开销不是算力）；写出 66 个场约 4 分钟。更快的同级卡省不掉这 16 次重载和 2000 次小步。
- **代价：** 112 片若都按这张片的大小，这张卡约 60 小时，不是 D-226 的 11–28 小时。正式“超过上限三倍”（45 分钟）尚未触发，但半天排期已失效。
- **复查条件：** 若改成权重只加载一次、或减少位移，再重测单片；未改代码前不用卡的标称速度外推。
- **影响：** 探针验收仍按原规则。全量是否在这张卡上铺开，按约 2 天墙钟重新决定，不按半天。

## D-230｜2026-10-07｜探针验收后关实例，不在这张 2080 Ti 上铺全量
- **背景：** 用户改口：探针验收后即关实例，并问 iStar/HookNet 是否吃显存、RTX 2060 6G 能否跑。
- **理由／选择：** 覆盖 D-228 的“通过就全量”。本机只跑完双片探针并取回验收产物，然后释放。显存实测：官方 HIPT 单张 4096 瓦片前向峰值 1839MB（保留 2232MB），开不开梯度一样，因为权重不要求梯度。两个 ViT 教师权重仅 168MB+94MB。训练是小网络，不吃显存。HookNet 权重 49MB、批大小 1、输入 1244，按结构应在 2–3GB 内，尚未在这张卡上实测。
- **代价：** 2060 算力约为 2080 Ti 的一半，CPU 阶段省不掉，全量墙钟会比 32 分钟/片更长；6G 卡要打开 TF 显存按需增长，避免启动时一次性占满。
- **复查条件：** HookNet 若在 6G 上 OOM，只降批大小或切块，不改标签规则。
- **影响：** 验收完成后释放当前实例。全量改到 6G 卡在显存上可行。

## D-231｜2026-10-07｜iStar 双片探针未过冻结线，不铺全量
- **背景：** 2080 Ti 上 HEST_TENX13 与 GSE c_57 官方管线跑完（各 66 基因场）。验收按 README：90/10 种子 0 留出重建，以及 c_57 超分辨 TLS 分数对独立标注的 AUC，对照 B-axis 0.8439 / 同点 HGB 0.8497 / 空间 HGB 0.8464。
- **结果：** HEST 留出平均 Pearson 0.233（中位 0.203）；TLS 相关基因接近零（MS4A1 0.019，CXCL13 0.0005）。c_57 留出平均 Pearson 0.084。c_57 TLS AUC 0.682（1084 点，69 个 TLS），低于三条基线约 0.16。
- **理由／选择：** 低于基线即无增益（README 规则 3，D-220）。关闭 iStar 路线，不铺 112 片。弱标签与超分场不进真值。HookNet 未再占用这张卡。
- **代价：** 只否定了「图像超分能超过现有逐点 RNA 分数」这一用法，不改 HF-01 的文字。
- **复查条件：** 若换不放大的真全片、且留出相关与独立标注同时超过冻结基线，才重开。本次 c_57 是放大后的 hires，HEST 无独立 TLS 标注。
- **影响：** 对 HF-01：图像插值没有提供可写入方案的空间恢复证据。实例在结果入库后释放。

## D-227｜2026-10-07｜GPU对照探针代码就绪，待算力登录方式
- **背景：** 用户要求先备代码再给算力。已交付 `infra/gpu_probe_20261007/`：环境脚本（numpy 1.25 钉死避开skimage冲突）、iStar官方管线脚本（跳过select_genes保护66基因名单）、HookNet探针脚本、两个验收脚本（留出重建90/10种子0、c_57 TLS AUC对入库基线、RNA高低区通道分离）。
- **理由／选择：** 全套自包含：iStar代码（commit 3cb0e53）＋HIPT权重（Git LFS官方，672M＋378M）＋HookNet权重＋探针输入目录已列清单。墙钟估算已交付（单片5～15分钟，对照合计约1小时）。验收规则事前冻结在README，不因结果改。
- **代价：** eval脚本未在真实GPU上验跑；登录后先跑环境脚本与单片即可发现适配问题。
- **复查条件：** 远端路径或CUDA驱动不同时，只改setup_env.sh的torch轮子，不改科学规则。
- **影响：** 待用户提供算力登录方式后执行。

## D-228｜2026-10-07｜测试通过直跑全量，完工即关机，绝不空置
- **背景：** 用户明确：对照探针通过后直接跑全量；跑完用ai-galaxy MCP关实例；任何情况下不允许空置服务器。
- **理由／选择：** 探针→全量无需二次确认（以冻结验收为准）；全量结束（或任一分支永久停止）后立即走plan_release→release释放；kept_disk同样及时清理。MCP当前可见1 running／0 kept_disk。
- **代价：** 释放后环境与中间产物丢失；验收产物必须先取回本地入库。
- **复查条件：** 若释放前验收产物未入库，不得释放。
- **影响：** 关机需要实例名（MCP无列表接口），待用户从控制台确认当前SSH机是否为MCP名下那台running。





## D-232｜2026-10-08｜四癌TLS区域标签与分层使用操作方案
- **背景：** 用户要求将本会话讨论的标签形式与下游用途落盘为具体操作方案。允许有误差的标签用于训练；需要区分结构参照、分子代理、置信度与成熟度。
- **理由／选择：** 采用实例区域＋spot阳性/阴性/未知两层标签；自动方法提供候选，独立病理复核提供结构参照。先在四癌本地材料上做有限标注试点，依据实际阳性预测值、漏检和位置误差判断候选用途，不以“必须超过RNA预测基线”作为标签可用性的统一条件。细节见 [操作方案](TLS_FOUR_CANCER_LABELING_PLAN_20261008.md)。
- **代价：** 独立标注人员和可用染色尚未落实；自动高分不等于正确概率。HEST112片包含健康与Tumor未细分材料，不能称112片癌组织或患者；10x43条目须与HEST去重。现有表达派生成熟度不升级为独立成熟度真值。
- **复查条件：** 患者/来源重叠、标注参考读取同一RNA定义、坐标失配、或跨癌种误差无法外推时，收窄相应标签用途并重估；未检查区域不得补阴性。旧方法失败原样保留，不能改阈值后追认旧验收通过。
- **影响：** 为FC-00增加四癌标签准备工作，服务后续HF-01/02/04；本次仅完成方案（基础设施），没有新增标签或实验。具体候选生成不自动重开iStar全量/GPU；资源限制与外部获取授权保持，plan科学假设不变。

## D-233｜2026-10-08｜四癌资格清单不把病人号或14片试点锁成样本
- **背景：** 方案要求先做身份、来源重叠和图像资格，再选最多14片。本地112片Visium/Visium HD文件都在，但病人字段大多为空。
- **理由／选择：** 病人号只在同一研究内使用。同一研究里“Patient 1”对应多个不同样本编号时，身份改回未知，不合并成一人。捕获区后缀A1/B1/C1/D1不当成另一个样本。估计像素尺寸在0.05到0.6微米且短边至少8000像素、坐标99%落在图内，才记为全片候选；嵌入像素尺寸为0或1不算。10x按页面或标题标可能重复，不加入112。目视配准未做，试点不锁定，不导出标注图像包，不生成标签。
- **代价：** 能对上的癌组织病人只有一个肠癌研究的6位，加两片肝癌的两位。肺和乳腺在这批材料里没有可数病人。外部测试集不能从现有metadata形成。
- **复查条件：** 若把未知身份凑进14片并称独立病人，或把P1–P5的“Patient 1”合并，相应名单作废。目视配准或新的病人对应表可以重开试点选择。
- **影响：** 清单在 `infra/tls_four_cancer_labels_20261008/eligibility/`。这是基础设施，不支持或否定HF-01。

## D-234｜2026-10-08｜四癌数据重检：8片配准错位、HD尺度不可用55换算、独立单位被高估
- **背景：** 用户要求重新核对数据本身的问题，不只看身份清单。
- **理由／选择：** 逐片重检全片头、h5ad结构、坐标落位、面板覆盖、物种与尺度来源。像素尺度改以元数据为准，并用 scalefactors 交叉验证：99 片普通 Visium 的 55 微米换算与元数据一致；13 片 HD 隐含常数为 128 微米，55 微米换算无效，不用它推尺度。溢出判定以一个点直径 55 微米为容差。
- **结果：** 112 片全为人，无重复条码、无缺失坐标、行数一致、全片尺寸一致。8 片坐标超出 91–1524 微米（TENX68/72/128/155/168/169/170/171）判配准错位不可标注；4 片溢出 21–38 微米记边界约定。NCBI828 标 Healthy 但注释为肝炎/急性肝衰竭。TENX28/23 面板可用 0/3 项。6 片像素大于 0.6 微米。6 组捕获区（A1/B1/C1/D1）与 ZEN 的 Rep1/Rep2 使独立单位远少于 112。`spots_under_tissue` 有 9 片与实际不符。目录内另有 66 片未下载（33 非人物种、33 为 Treated/Diseased），10x 清单实际 73 条而非 43 条。
- **代价：** 可进入配准的候选少于 112；试点仍不锁定。肝的背景组织不可作阴性。
- **复查条件：** 若用 55 微米换算 HD 尺度，或把捕获区/重复当独立病人，或把溢出片当可标注，相关结论撤回。43 片全片图需走 PIL 而非 tifffile。
- **影响：** 问题清单见 `eligibility/DATA_ISSUES.md`。仍属基础设施，不改变 plan。

## D-235｜2026-10-08｜配准核验完成，更正坐标越界判据，并锁定试点提案
- **背景：** 用户要求继续推进。配准检查需回答点坐标是否落在组织上，并导出可给人看的叠加预览。
- **理由／选择：** 判据写在脚本开头：落在组织轮廓内的点占比低于 0.90 判配准错位；非本身变换且高出 0.05 以上判需要变换。判定用“落点内占比”而不是“是否有任一点越界”。试点入选规则同样先写后跑：配准与面板门槛通过、一个标本只取一张、优先有病人号、病人号不足才用身份未知补位且不计入独立样本。
- **结果：** 96 张配准通过、18 张错位、无一需翻转或旋转。目视抽查四张与自动判定一致。合并门槛后可标注候选 71 张（肠53、乳腺10、肺6、肝2），对应 38 个标本、其中仅 8 个有可数病人。试点提案 10 张、独立病人 6 位，留出 2 个标本。
- **更正：** D-234 把 12 片坐标越界列为阻断，过重。除 TENX170（4.10%）与 TENX171（2.83%）外，其余越界点占比均低于 0.7%，属单个边缘点。质心偏移量在通过与错位两组间不分开，不得作判据。
- **代价：** 外部锁定测试集仍不存在；留出的 2 个标本不构成验证集。独立参照因无病理标注者而 NOT_RUN。
- **复查条件：** 若把身份未知的 4 张计入独立病人，或把 2 个留出标本称作外部验证，或在无标注图册与标注者的情况下开始画 TLS 区域，相关结论作废。
- **影响：** 包在 `annotator_pack/`（10 张切片、50 块 1 mm 视野、纯 H&E 无点圈）。仍为基础设施，不支持也不否定 HF-01。

## D-236｜2026-10-08｜判读者为项目用户本人，专家级人工判读，单阅并如实登记
- **背景：** 用户要求由自己完成病例标注，并要求提供数据与操作方法。
- **理由／选择：** 先写死判读规则与输出格式（D-232 要求规则先于判读），再提供工具。证据等级用 D-232 原定的 `independent_reference`，因为判读只依据 H&E 形态、不读取任何 RNA 读数或模型输出。强制记录 `reviewer_credential=expert_human_review`、`double_read=false`、`annotation_protocol=blind_first_pass`。资质按 D-036 与 D-237 处理，不区分资历等级。
- **盲法依据：** 这 10 张从未跑过任何自动 TLS 检测（HookNet 未运行、iStar 已按 D-231 关闭、ST 检测法已按 D-219 判死），因此不存在可被污染的自动候选答案。如实登记的一点：判读者同时是本项目分析的设计者，不是外部第三方读者。
- **代价：** 单阅，无分歧记录与裁定；4/10 切片身份未知；独立病人 6 位；外部锁定测试集不存在。
- **复查条件：** 若出现第二判读者，须重做双阅并记录分歧，现有 v1 标注不得追认为双阅。若把这批标签用于确认性结构主张，须同时给出单阅这一限制。
- **影响：** 工具 `scripts/tls_annotate_gui.py`（自建，PyQt5＋PIL 金字塔，无外部依赖）、校验 `scripts/tls_annotation_check.py`、说明 `docs/ANNOTATION_GUIDE_20261008.md`。标注结果当前为空，independent_reference 仍 NOT_RUN。

## D-237｜2026-10-08｜用户裁定：人工审核内容一律按专家级建议登记
- **背景：** 我在 D-236 中把判读者写成“非执业病理医生，证据等级下调”。这是未经询问的臆断，且与 D-036（人类标注者不区分资历等级、统一按高置信度处理）相悖。用户明确纠正。
- **用户裁定：** 任何涉及人工审核的内容，都应当按专家级（病理医生/博士）建议登记，除非用户特别说明其非专业。
- **理由／选择：** 取代先前的臆断，并把 D-036 从“不降级”推进为“默认按专家级登记”。本规则只管人工审核内容的登记口径，不放宽其它条件：标注模态仍须独立于拟用输入（D-036 边界）、几何可重放、provenance 闭合；由 RNA 签名生成的标签仍不得作为确认性真值；模型一致性不得升级为真值。
- **代价：** 登记口径不再反映资历差异，inter-rater 差异只能靠重复标注实验单独估计。单阅／双阅、判读者是否外部仍必须逐项如实记录，不受本规则影响。
- **复查条件：** 若某次任务的标注质量确实受资历影响，或项目进入需要分级的正式验证阶段，重审。
- **影响：** 已据此改写 `docs/ANNOTATION_GUIDE_20261008.md`、`TODO.md`、`STATUS.md`、`annotation/README.md`、`scripts/tls_annotate_gui.py`、`scripts/tls_annotation_check.py` 与 D-236。

## D-238｜2026-10-08｜交付离线标注包（无图形界面环境）
- **背景：** 用户所在环境无图形界面，要求把待判读数据与标注工具打包，Windows/Ubuntu 开箱即用。
- **理由／选择：** 本机无可用病理标注软件（无 QuPath/napari，装 QuPath 需外部下载且版本漂移），因此改为自建纯静态浏览器查看器：零依赖、无构建步骤、文件结构由本项目完全掌控。影像按切片预渲染成 JPEG 金字塔块，最高一层 1.1 µm/px。启动脚本只依赖 Python 3 标准库；无 Python 时可直接打开 `viewer/index.html`。
- **分辨率取舍：** 原生 0.274–0.548 µm/px，全原生打包约 1.4 GB；1.1 µm/px 打包为 194 MB，体积是前者的 1/7。经实测 1.1 µm/px 下隐窝、上皮与固有层淋巴细胞形态清晰可辨，足以判聚集是否成形与密度；不足之处是生发中心内精细结构（套状区、滤泡树突网）。若需该尺度，按 `--target-um 0.55` 重建。
- **验证：** node 桩环境对真实 app.js 做端到端测试（真实 HTTP 取块），33 项全通过，含取数链路、瓦片寻址、多边形落点与坐标换算回原始像素、覆盖率算法、导出与导入回环。另以同一套寻址数学在软件里合成渲染，确认瓦片拼接连贯、无错位。瓦片寻址随机 4000 次抽样零缺失；全 10 片各层瓦片数与理论值一致；zip 完整性测试通过。
- **过程中修掉的真实缺陷：** `dh = dw` 使边缘块被拉成正方形（右侧点圈明显变形）；`$("slideSel").onchange()` 在严格模式下 `this` 为 undefined 会中断初始化；选层方向取反导致按原生而非目标分辨率打包（体积 4 倍）。
- **代价：** 交付件分辨率低于原生；单阅限制不变。
- **复查条件：** 若判读认为 1.1 µm/px 不足以判断生发中心内结构，按 `--target-um 0.55` 重建后重做相应切片，已有标注按坐标保留不受影响。
- **影响：** 包 `dist/tls_annotation_pack.zip`（179 MB，sha256 `43f4718…`，2026-10-09 重建后为 `95662ed472a404092d56f9f620f64b0e09e5332eefdea85236cbe2acdf19f748`）；源 `infra/tls_annotation_pack_src/`；构建 `scripts/tls_build_portable_pack.py`；测试 `scripts/tls_viewer_test.js`。不改变任何标签状态。

## D-239｜2026-10-09｜修正标注器黑屏：相对路径误解析，并修测试盲区
- **背景：** 用户在真实浏览器中报告 UI 正常但画布黑屏。
- **根因：** `fetch("data/index.json")` 是相对路径，文档位于 `/viewer/`，浏览器解析为 `/viewer/data/index.json` → 404 → 取数失败 → `S.meta` 为空 → 画布永远停在底色。附带缺陷：瓦片解码完成后无重绘回调（画面需用户操作才显影）；画布尺寸仅在 `load` 时取一次；文档链接同样指向错误目录。
- **修法：** 根路径改为从 `document.currentScript.src` 推导，并保留 `../`、`./` 候选回退；所有数据与瓦片路径经该根拼接；瓦片 `onload`/`onerror` 触发重绘；画布接 ResizeObserver；失败时在画布上直接画诊断文字并在顶栏报错，不再静默。
- **测试盲区（更重要）：** node 桩原先把相对路径拼到站点根，因此**恰好掩盖了这个 bug**，所以 40 项断言全绿而真实浏览器黑屏。改为按浏览器规则用 `new URL(u, documentUrl)` 解析后，旧代码立即报 `index.json 404` 失败、新代码通过。测试从“证明当前代码可用”升级为“能抓住这一类回归”。
- **验证边界：** 无头 `firefox --screenshot` 在本机是按首次绘制抓取而非等 `load`，无法作为画布渲染证据；Marionette 在本版本不回命令。已验证的是：URL 解析（浏览器等价）、瓦片寻址（4000 随机视图零缺失）、拼接几何（同一套数学的软件合成图像连贯无拉伸）、解码触发重绘、只绘制真实瓦片。**未在真实浏览器中亲眼看到画布渲染**。
- **复查条件：** 若在有图形界面的机器上仍见黑屏，直接看画布左上角的诊断文字与浏览器控制台，据此定位，不再猜测。
- **影响：** 服务已在本机 `127.0.0.1:8848` 运行；包已重建并通过 zip 完整性测试。

## D-240｜2026-10-09｜标注包扩至全部71张候选；修选层容差与覆盖率缩放，规则改发HTML
- **背景：** 用户要求把剩余切片全部纳入在线标注；并报告判读规则打开是乱码。
- **乱码根因：** `http.server` 对 `.md` 不发字符集，浏览器按本地编码解 UTF-8。改为构建时用 `scripts/md_to_html.py` 生成带 `<meta charset="utf-8">` 的网页，侧栏改链。转换后逐项核对：标题 13/13、表格 3 个无裸管道符漏入段落、自查勾选框 19/19、抽查中文全部保留、无乱码特征串。
- **选层容差 bug：** 原判据要求分辨率 ≤ 目标；41 张 MISC 肠癌片原生 0.2757 µm/px，逐层 0.551、**1.103**、2.205，而目标写 1.10，于是 1.103 > 1.10 被卡在 0.551，线性差一倍、像素差四倍。这 41 张吃掉约 3 GB。改为允许 8% 容差。
- **两档分辨率：** 试点 10 张加全部 14 张有病人号切片给细档 1.10 µm/px（316 MB），其余 53 张给粗档 2.20 µm/px（370 MB），合计 759 MB。粗档中有病人号的切片为 0 张，能支撑病人级结论的材料未被降分辨率。全细档为 3987 MB，不可行。
- **覆盖率缩放 bug：** 覆盖率算法用组织轮廓的**原始全分辨率**坐标去测**基准层**网格上的点，未做缩放，因此返回空值或错误数字。直接影响到标注者判断“第一遍覆盖是否走完”。已按 `base_w/full_w` 缩放，并加回归测试：未缩放版本会使测试报两项失败。
- **测试补强：** 新增断言“组织轮廓缩放后落在基准层网格内”与“整片覆盖多边形应报约 100%”。
- **代价：** 粗档 2.20 µm/px 下淋巴细胞核约 3 像素，够看成形与密度，不够看细胞形态；该层切片不用于细胞层面判读。
- **留出集消失：** 原预留的两块肠癌标本（ZEN38/ZEN42、ZEN36/ZEN40）已含在 71 张内。全部候选被标注后，批内不再有任何按病人留出的切片。后续评价必须用**预先固定的病人级划分**，划分口径须在看到标注结果之前定。
- **复查条件：** 若粗档切片需细胞层面判读，按 `--bulk-um 1.10 --no-two-tier` 重建（约 4 GB）。若覆盖率读数仍与预期不符，先看顶栏 µm/px 与侧栏覆盖百分比是否同时可信。
- **影响：** `dist/tls_annotation_pack.zip`（646 MB，sha256 `50204b22b36823843e99819e57eae97e1994e20ba5d53d4b9ef1b43a62ca9f69`），71 张：肠 53、乳腺 10、肺 6、肝 2。不改变任何标签状态。

## D-241｜2026-10-09｜修“始终在加载”：零尺寸画布使 fit 产出 scale=0 并永久卡死
- **背景：** 用户报告有 GUI 但读不出图片，页面始终在加载。
- **根因：** `fit()` 无零尺寸保护。画布尚未布局时 `clientWidth` 为 0，算出 `scale = 0`。后续 `pickLevel()` 的 `log2(1/0)` 得 `Infinity` 被兜底为 0（最细层），而可见范围除以 0 得 `NaN`，循环条件恒假——**不请求任何瓦片**，画布永远空。叠加 `fitted` 标志在改写时丢失，首次 fit 失败后不再自动重算，于是不会自愈。
- **次因：** 单帧无瓦片数量上限，任何退化视图都可能一次性请求整张金字塔，使标签页永不完结。
- **修法：** `fit()` 在画布小于 8 像素时置 `needsFit` 并返回，不写 scale；`draw()` 在无画布尺寸时直接返回，不绘制；`size()` 发现 `needsFit` 且画布已有尺寸时重算；新增 `MAX_TILES_PER_PAINT=400` 与 `chooseLevel()`，可见瓦片超预算即自动降层；退化视图在顶栏报出画布尺寸与缩放值。
- **测试：** 新增 `[3d]` 零尺寸场景与 `[3e]` 极端缩放下单帧瓦片数上限两组断言，共 50 项全通过。人为去掉零尺寸保护后测试报 `nothing painted while the canvas has no box` 失败，确认可拦回归。
- **交付核对：** 71 张全部 `meta.json`、`tissue.geojson` 齐全，各层瓦片数与理论值一致（0 失败）；经 HTTP 实取每张每层的四角与中心共 **1168/1168 全部 200**；14 张有病人号切片全部在细档。
- **影响：** 包已重建（646 MB，sha256 `50204b22b36823843e99819e57eae97e1994e20ba5d53d4b9ef1b43a62ca9f69`）。仍未在真实浏览器中亲眼确认渲染（D-239 已说明本机无头抓取不可用）；若仍异常，请取顶栏文字，它现在会直接报出画布尺寸与缩放值。

## D-242｜2026-10-10｜等待专家判读期间：建成映射管线，并发现自动门槛误排了人工可判的片
- **背景：** 判读已交外部专家。等待期做不依赖标签的工作。
- **映射管线：** 新建 `scripts/tls_map_labels.py`。冻结的保守三态规则：bin 的全部采样点落在实例内为阳性；全部落在已检查区内且不碰任何实例为阴性；其余（部分交叠、冲突、未检查、配准存疑）为未知。**未知绝不并入阴性。** 实例即使没有 bin 被完整覆盖也保留，故小结构不丢。输入支持浏览器导出的 `annotations.json` 或每片 GeoJSON。用合成标签跑 8 项规则断言全通过（含“碰实例的点绝不判阴性”“未检查区判未知”“三态划分穷尽所有点”）。一次试跑：NCBI642 得阳性 1、阴性 1432、未知 554（未知占 27.9%，符合保守规则的预期）。
- **门槛错排（新发现）：** 被排除的 20 张里，**7 张是按自动管线门槛排掉的，人工判读并不需要**。TENX28／TENX23 只因 1742 项分子面板可用 0–3 项被排除，而其配准为 0.997／0.999；另 5 张只因像素大于 0.6 µm/px（0.69–0.91）被排除，而 0.6 这条门槛是给 HookNet 定的，其配准为 0.987–0.9995。人工在 0.9 µm/px 下判读聚集成形与密度无碍。
- **配准错排：** 18 张配准阻断中 17 张落在 0.80–0.90，仅 MEND36 为 0.545 且已目视确认为真错位。其余更像局部偏移或组织轮廓不全，而非全局坐标系错误，因此存在通过逐点校正或剔除越界点挽回的可能。此项尚未执行。
- **代价／待决：** 加回 7 张需同步告知正在工作的专家，属其决定；身份问题不因标注推进而缓解。
- **复查条件：** 若把未知当作阴性、或按阳性直接调下游阈值，相应结果作废。若配准修复采用逐点校正，阈值必须在修复前定，且不得据修复后结果回改。
- **影响：** 映射管线为基础设施，不支持也不否定任何假设。标签一到即可运行。

## D-243｜2026-10-10｜人工门槛与自动门槛分离；发现一对字节级重复；HookNet CPU 可行
- **背景：** 用户要求把 7 张人工可判切片加入标注网页，并完成身份核查、配准修复、HookNet 计时三项。
- **门槛分离：** 原候选池按自动管线需求筛选（0.6 µm/px 是 HookNet 的纹理要求，分子面板是打分要求），二者对病理判读均不成立。新增 `gate_status_human`（配准通过 + 像素 ≤1.0 + 瘤组织语境）。多出 4 张癌切片（TENX23/28/29/49）；3 张肝（NCBI828/829/830）标 Healthy 且 NCBI828 说明写肝炎/ALF，留在背景类、不得当阴性，查看器直接显示该矛盾。
- **字节级重复：** TENX28 与 TENX29 的 TIFF 逐字节相同（md5 `b63fbb04b50c9fe85d7424769bac2cad`，359,956,737 字节，轮廓 IoU=1.0000），HEST 用两个 id 登记同一张图。全池 75 张候选用「尺寸+头尾各 1 MB 哈希」指纹得 74 个不同指纹，唯一冲突即此对，方法用已知重复对验证。**故标注集为 74 张（原 71 + 3）**。
- **身份核查（否定）：** 本地 h5ad 无 donor/patient 字段，无 cxg 目录；CZ CELLxGENE 397 个公开集合无 COLON MAP；网页检索未得到带标本表的补充材料。COLON MAP 为 36 张切片、7 个编号（6723/7003/7319/7794/8270/8578/8899），横跨两篇 Cell（S0092-8674(21)01381-7 与 (23)01222-9）；其余 25 张多为 10x 演示集，病人身份通常不公布。解开需补充表或联系作者，属外部获取，待批准。
- **配准修复（否定）：** 对 17 张内点 0.80–0.90 的切片做平移搜索（±1000 px，40 再 8 px），判据先写后跑。全部 17 张在 (0,0) 即全局最优，增益 0；两张看似有增益的算出 12,370 px 与 20,610 px 平移，落出图像被 R4 否决。故**不是整体坐标偏移**。机制进一步未能确定：MISC 系列（有白背景峰）上组织外点大多确在背景（轮廓大致正确）；MEND/FFPE 系列与 TENX170 无背景峰（近白像素 0–3%，中位 192–206），固定阈值会全判成组织、逐图 Otsu 会把阈值压到 6–14，**全局阈值分不开组织与背景**。**处置：17 张继续排除**，因为无可信组织掩膜时无法保证区域映射到正确的点。
- **HookNet CPU 与 GPU 估算（同一条数字改了两次，以本节为准）：** 首在真全片计时（hooknet 环境 TF 2.13.1，CUDA_VISIBLE_DEVICES=-1）。ZEN39 用第 0 层 0.4595 µm/px。**更正一：** HookNet 细圈 1244 px × 0.5 µm = **622 µm 视野/patch**，我先前声明的 32 µm 栅格密了约 380 倍。**更正二：** 我把 0.587 s/patch 写作“单核”但未验证；强制单线程重测为 **9.231 s/patch**，即实际并行加速仅 **15.7 倍**而非 64 倍，先前的墙钟估算都乐观了约 5 倍。全池 74 张染色组织约 3136 mm²，按满机 0.587 s/patch：622 µm（无重叠）8,106 patch ≈ **1.3 h 推理**；311 µm（50% 重叠）≈ 5.3 h；64 µm ≈ 124.8 h；32 µm ≈ 499.4 h。另加 TIFF 解码 20–40 分钟（可与推理重叠）。**按原生视野跑全池：约 1.5–2 小时，主要是推理；GPU 只把它压到约 25–45 分钟。** V100 只在刻意铺密栅格时决定性（32 µm 从约 21 天降到 1–2 天，按 10–30 倍估计，无 GPU 实测）。只测吞吐，未测检出能力。
- **代价：** 标注集扩到 74 张后包为 799 MB；身份问题不因标注推进而缓解。
- **复查条件：** 若有人把 Healthy 切片当阴性、或把 TENX29 当第二份材料、或声称 17 张配准失败已定位到轮廓，相应结论作废。HookNet 铺格必须按 622 µm 原生视野计算，用更密的栅格前需说明理由。
- **影响：** 详见 `infra/tls_four_cancer_labels_20261008/eligibility/PARALLEL_FINDINGS.md`；HookNet 计时原始产物存 `infra/tls_four_cancer_labels_20261008/hooknet_probe/`（JSON + 日志）。仍属基础设施。

## D-244｜2026-10-10｜HookNet 输入范围错了一个数量级；先前探针结论作废；GPU 包就绪
- **背景：** 用户提供 1 小时 V100 窗口，要 50% 重叠。开窗前必须证实脚本可用。
- **发现的错误（关键）：** HookNet 模型内部有 `Rescaling(scale=1/255)`，它要**原始 0–255** 输入。我的运行器又除了一次 255，网络看到近全黑图，输出退化为常数：纯白、随机噪声、真实组织给出几乎相同结果，TLS 通道**恒为 0.0000**。改为原始范围后 TLS 通道最大 0.968、空间标准差 0.15、类别和恰为 1.0000。
- **连带作废：** `scripts/tls_hooknet_crop_probe.py` 同样除了 255（`/ 255.0`），所以 D-209“HookNet 在 GSE 小图上饱和、无鉴别力”的**实测基础不成立**。D-209 的分辨率论据（需 0.5 µm/px 纹理）与输入范围无关，仍然成立；但“实测饱和”这一条不能再引用为证据。HookNet 线是否可行，需用正确输入重测。
- **采样约定（核过源码，不是推断）：** 视野 622 µm（1244×0.5）；官方步长 512 µm（overlap 220）；本包按用户要求 311 µm（50% 重叠，比官方密 2.7 倍）；输出 1030×1030×3；放置取输出**中心 1024** 于 `info.x−512`（`writing.py` 的 `TILE_SIZE=1024` 与 `apply.py` 坐标式）；**TLS＝通道 2**（`postprocess.py` 的 `inv_label_map={2:'tls',3:'gc'}`）。输出中心对输入中心的依据：`apply.py` 对标注做中心裁剪，否则模型学不出来。
- **包与脚本：** `scripts/tls_hooknet_stage.py`（每张源层重采样到恰好 0.5 µm/px + 位置表，约 3–5 GB，而非 649 GB）、`scripts/tls_hooknet_gpu_run.py`（自盯时间、可续跑、逐片原子落盘、信号处理）、说明 `infra/hooknet_gpu_20261010/README.md`。
- **本地已验证：** 推理路径、续跑（第二次跳过已完成片）、预算耗尽不新开片、mosaic 形状与覆盖率、类别和=1。另修两个自造 bug：`np.savez_compressed` 会在非 `.npz` 结尾补后缀；同进程重复构造 HookNet 会抛 `_build` 类属性异常。
- **时间：** 本地 CPU 实测 0.563–0.691 秒/位置 × 32,425 位置 ≈ 5.4 h；V100 估 11–32 分钟（**估算，无 GPU 实测**）。
- **代价：** 本轮只产出概率图，不判定对错；产物不得给正在盲读的专家看。
- **复查条件：** 若有人引用 D-209 的“实测饱和”，须同时引用本条。若 GPU 上 TLS 通道仍恒为 0，先查输入范围而非模型。
- **影响：** 待用户提供登录后开跑。

## D-245｜2026-10-10｜HookNet 在 V100 上跑通；两个自造错误；实测 0.378 s/位置
- **背景：** 用户提供 V100（SXM2-32GB，驱动 550.127.05；内存 15 GB、8 核）。开机、上传 6.4 GB、装 TF 2.13.1，但 TF 起初看不到 GPU。
- **自造错误一（决定性）：** 运行器里有 `os.environ.setdefault("CUDA_VISIBLE_DEVICES", "-1")`——为本地测试写的，机上未设该变量时默认生效，**每次都强制走 CPU**。之前所有的"GPU 时隐时现"都是这一行造成的，GPU 自始至终是好的。改为仅在 `--cpu` 时禁用。另外系统 CUDA 是 12.4 而 TF 2.13 要 cu11，机上到 PyPI 不通，改为本机下载 cu11 的 10 个 wheel（1.8 GB）上传离线安装；并补了 `libcuda.so` 无版本号符号链接。
- **自造错误二：** 把早前 CPU 回退时的 `std::bad_alloc`（batch 16 的激活内存在 15 GB 上爆掉）误判为图像过大，曾删掉上下文层反而让取块变慢。改为整片解码一次、上下文层只算一次。
- **实测（非估算）：** 首个片 0.707 s/位置（含 XLA 预热），热机 0.378 s/位置。全量 32,090 位置约 **3.4 小时**，不是先估的 11–32 分钟。先前所有 HookNet 的工时估算作废，以本节为准。
- **运行安排：** 运行器自盯时间（--budget-min 55）、可续跑、逐片原子落盘、收信号先写完当前片、每片失败不炸掉整个运行、缺的切片跳过并稍后重扫。运行顺序改为**有病人号切片优先**（约 14 张、约 2300 位置、约 15 分钟），保证窗口耗尽前最有价值的部分先完成。
- **代价：** 预算内大概率只跑完一部分；未完成部分靠续跑补齐，或后续 CPU/别的卡补。
- **复查条件：** 若有人引用先估的“V100 几分钟跑完全集”或“单核 0.587 s/patch”，须同时引用本条。若运行中 TLS 通道恒为 0，先查输入范围（应为原始 0–255）。
- **影响：** 运行进行中；监视器已挂，完成后拉取产物。关机在产物取回后，按用户指示走 ai-galaxy（可能因非 MCP 创建实例被拒，届时用账号释放接口按 SSH 端口 24432 匹配）。

## D-246｜2026-10-10｜HookNet V100 全量完成 74/74，产物入库，实例已释放
- **结果：** 32,090 个位置全部完成，0 跳过，`stopped_by_signal: false`。74 张 mosaic（prob uint8×3 通道 + coverage）全部满足：3 通道、覆盖非空、行和≈1。本地抽查 6/74 通过。产物在 `infra/hooknet_gpu_20261008/`（277 MB），另有每片 sidecar JSON 记载间距、形状、放置约定与耗时。
- **释放：** 用户明确指示完成后经 ai-galaxy 关机。`plan_release` 对非 MCP 创建实例拒绝（上次已验证），故走同后端账号接口 `release_instance`：`17750662501_lyg2144_a6e886d22d9248579c5094df9f536c53_u1_c0` 返回 Status 0，退款 1.30。复查：实例列表中该实例为 0、全 44 个实例无一运行中、保留磁盘 0。**未从机器内部关机**（那不停费）。
- **通知：** Bark 普通消息两条，均 code=200：完成后一条（74/74、32,090 位置、本地验证通过），释放确认后一条（已停、0 运行、0 保留盘）。
- **代价：** 计费按 1.254/小时自动续费多跑了约两小时；接受，因用户明确要求跑完全部、不设时间预算。
- **复查条件：** 概率图是 `automatic_candidate`，不得直接当标签；与专家标签的比较（实例级查全/查准、点级 PR-AUC、定位误差）待标签到手后做，见 `infra/hooknet_gpu_20261010/README.md` 第 9 节。

## D-247｜概率图叠加层提前上线，之后所有导出记 assisted 复核（v2）
- **背景：** 用户决定不等盲读锁定，现在就把 HookNet TLS 概率叠加进标注页。
- **理由／选择：** 技术上热力块与 H&E 共用同一坐标系（mosaic 格 4µm → base 像素仿射换算，逐片记录），叠加开/关、阈值、不透明度可调，阈值只改显示不改数据。但代价是证据等级：之后所有导出（实例与覆盖区）记 `reviewed_candidate` + `assisted_review_single_reader` + `annotation_version v2`，不再是盲读。叠加是否开过按切片记录（`overlay_used`），随导出一起走。
- **实现与验证：** 热力块只取 TLS 通道（索引 2），mosaic 原生分辨率起层、不做上采样虚增（全 74 片约 109 MB，而非放大到 base 的 469 MB）。viewer 侧灰度块经阈值+红色 ramp+透明度的 LUT 实时渲染，阈值/透明度变化只清 LUT 缓存。node 桩端到端测试更新断言并通过（含叠加块 HTTP 可达、热点数值 247/252 在块内完整保留）。
- **代价：** 之前若已有盲读导出（v1），与之后的 v2 不得混用；引用时必须注明版本。包体积 646→763 MB。
- **复查条件：** 若有人把 v2 标注当独立参照引用，相应结论作废。叠加阈值的选择不进入标签，只影响显示。

## D-247 附记｜叠加层零绘制 bug：未定义变量被外层 try 吃掉
- **现象：** 用户在任何片子都看不到叠加，但服务端一切正常（新版 viewer 在线、74 张全带热力元数据、块 HTTP 200）。
- **根因：** `heatTile` 里 `out._heatUrl = u` 引用的 `u` 在该作用域不存在 → ReferenceError → 被 `draw()` 外层的 try/catch 静默吞掉。URL 照发（所以"块可达"测试是绿的），但一张 canvas 都没画出来。
- **教训与补救：** 沿用 D-239 的教训——只断言"请求发出"不够，必须断言"画上去了"。测试新增 [7c]：开叠加、等解码、重画，断言处理过的 heat canvas 确实进入 `drawImage`；另修了测试自身的竞态（`exportSlide` 第二个文件经 350ms 定时器落地，改成立刻按名捕获 envelope）。回归验证：把该行改回 bug 版，测试报 `processed heat canvases painted to screen n=0` 失败。
- **另：** 全池热力信号最强的是 TENX23（>128 格占 17.7%，max 253），给用户做验证用。TENX23/TENX24 数值高度相近（不同 panel，疑似同块组织），TENX13/14 为相邻切片，均属正常相似，不作重复处理。

## D-248｜2026-10-14｜跨模态三角验证（HookNet 图 × RNA TLS 签名）+ 通道身份实证 + 划分冻结 v1
- **做了什么：** 在 65 片（74 图 ∩ 99 评分）上，把 HookNet TLS 概率在 spot 足迹（55µm 盘≈6.9 格@4µm）内平均，与 14 基因 TLS 签名（log1p(CPM) 逐基因 z 均值）做 Spearman 相关 + 100 次环面平移零假设。坐标链：pxl_*_in_fullres × pixel_size_um_estimated ÷ 4 = mosaic 格（注意 geometry npz 的 xy 是 hex 阵列坐标，不能用于图像映射）。脚本 `scripts/tls_hooknet_rna_triage_20261014.py`，产物 `infra/tls_hooknet_rna_triage_20261014/`。
- **通道身份实证：** 官方 io.py 的 heat1=tls/heat2=gc 命名与保存通道的对应关系不明；用生物学约束裁决——GC 必须嵌在 TLS 内。高信号片上 ch1⊂ch2 邻域比例 0.79-0.81（反向 0.10-0.26），且全池 rho 排序 tls>gc>other（other 中位 -0.050 为负）。**ch2=TLS、ch1=GC、ch0=肿瘤/其他**，与此前所有叠加工作一致。
- **结果（弱标签级，非验证）：** 全池中位 rho_tls=0.029，13/63（20.6%）p<0.05。信号集中在**已知患者结肠**（ZEN36 0.349 / ZEN38 0.282 / ZEN44 0.256 / ZEN42 0.244，均 p=0.01，6/12 ZEN p≤0.01）；乳腺/肝/肺中位≈0，无显著片。TENX23 弥散巨块 rho=0.005（p=0.49）——弥散激活不追踪 RNA，支持"弥散≠TLS 生物学"的读法。
- **解读边界：** rho 0.25-0.35 是 spot 级稀疏衰减后的中等相关；这是两模态互证，不是真值验证。乳腺等器官为零可能是 HookNet（乳腺训练集）迁移差、RNA 签名稀疏、或真实无 TLS，三者暂不可分。
- **划分冻结 v1（草案待用户确认）：** `docs/SPLIT_FREEZE_20261014.md` + `split_folds_v1.tsv`（sha256 03348878…）。66 单元（8 已知患者 + 58 单例；TENX23/24、TENX13/14 强制同单元），5 折分组 CV，seed=20261014。标签被查看前生效。
- **评估管线：** `scripts/tls_eval_hooknet_20261014.py` 自测全绿（合成标签 recall/precision/PR-AUC=1，打乱图 PR-AUC≈流行率，unknown 不计入 negative）；真实路径等 annotations。
- **分诊产出：** `triage_spots.json` 每片 top-10 图像高/RNA 低冲突点（fullres 坐标）；标注优先级建议：ZEN36/38/42/44（双模态一致+已知患者）→ 冲突区。

## D-249｜2026-10-14｜划分冻结 v1 生效；叠加层人眼确认；enter-wsl 暂缓
- **用户确认 SPLIT_FREEZE v1 生效**（`docs/SPLIT_FREEZE_20261014.md`，fold 表 sha256 03348878…）：66 单元 5 折分组 CV，自此刻起任何标签相关评估按此执行。修订只走文件中声明的失效条件。
- **用户人眼确认 viewer 热力叠加可见**：D-247 附记的零绘制 bug 修复经真实浏览器验证关闭。该修复路径（服务侧元数据→LUT 渲染→回归测试 [7c]）定型。
- **enter-wsl GPU（100.109.41.95）暂搁置**，疑网络问题，用户指示以后再说。当前工作均为 CPU 可行，不构成阻塞。

## D-250｜2026-10-14｜三角验证 NaN 修正（更正 D-248 数字）+ 实例库分层 + COLON MAP 解阻草稿
- **自造错误（已修）：** 环面零假设在 rho_tls=NaN（常数通道=全零 TLS 图）时假性显著——NaN 参与比较使所有 null 都不 ≥ 观测，p 恒为 0.0099。5 张全零图（MISC45/46/47/49、TENX152）被误判显著。**D-248 的"13/63 (20.6%) 显著"作废**，修正为 **9/58（15.5%）**，且 9 张全是结肠：ZEN36/38/44/42/45/40（已知患者）+ MISC54/63/34。中位数不变（tls 0.029 / gc 0.048 / other -0.050）。脚本已加 CONSTANT_CHANNEL 单列，per_slide.tsv/summary.json 已重建。教训同 D-239/D-247 附记：断言必须落到最终数值的有效性，而非流程跑通。
- **实例库与分层（`infra/tls_instance_library_20261014/`，6902 实例）：** 分层规则冻结（diffuse=最大块≥0.79mm²；focal=有≥0.03mm² 块；其余 sparse）。结果：diffuse 4 张**全是乳腺**（TENX13/14/23/24，4.3–12.8mm²，rho≈0），与三角验证互证"弥散≠TLS"；focal 19；sparse 51。**sparse≠差**：ZEN42/44 仅 4–6 块但显著追踪 RNA（少而准）。评估时 diffuse 层必须单列。GC 嵌套率在 ZEN 上 0.14–0.47（乳腺训练模型在结肠的 GC 识别偏弱，记录备查）。
- **COLON MAP 解阻草稿**（`docs/COLON_MAP_UNBLOCK_20261014.md`）：路径 A=PMC 补表自查（PMC8941949/PMC10756562 附表检索 7 token），路径 B=作者邮件草稿；发送/下载均需用户批准。解阻后按 SPLIT_FREEZE 失效条件重新生成 fold。
- **对 PLAN 的影响：** HF-04 方向获得第一个跨模态正信号（图像模型 × RNA 签名在已知患者结肠一致，9/58），但仅结肠、仅中等强度；乳腺等三器官为零（迁移差/签名稀疏/真无 TLS 不可分）。不构成 HF-04 成立，是评估框架就绪 + 信号存在性的证据。

## D-251｜2026-10-18｜用户将存储硬下限调整为200GB
- **背景：** 接手核对实测可用空间约1.054 TB，低于原1.2 TB；用户明确要求调整到200GB。
- **理由／选择：** 现行硬下限为200 GB（200,000,000,000 bytes，十进制），显著落盘前继续估算峰值并核对实际余量。原1.2 TB及历史2 TB仅保留历史含义。
- **代价：** 可用安全缓冲减少；不得据此默认复制全库或删除数据。
- **复查条件：** 预计或实际余量低于200 GB时，暂停新增计算落盘并报告存储硬阻塞。
- **影响：** 当前存储阻塞解除；GPU禁用、外部获取及对外发送授权规则不变。

## D-252｜2026-10-18｜算法先发现，人工TLS标签仅用于最后验证
- **背景：** 用户明确要求本阶段不使用人工真TLS标签，直接用算法开始发现，最后再用人工标签验证。这覆盖D-232等旧方案中人工标注先行和本阶段监督训练的安排。
- **理由／选择：** 本阶段允许原始表达、图像、坐标、非目标metadata和有来源的固定算法代理；禁止读取人工TLS区域/类别/成熟度用于训练、选样、通道判定、阈值调节、选模型或提前评价。人工标签可独立制作并封存，发现不等待回传；冻结候选、模型、预测、阈值和评价合同后才进入最终结构验证。
- **代价：** 既有FC监督比较属于已看标签的历史开发证据，不能重命名为新盲法发现。预训练HookNet含外部监督先验，不是从零无监督的生物场发现；与RNA一致也不等于TLS真值。已看模型叠加的人工判读必须标assisted，不能承担该模型的独立确认。
- **复查条件：** 发现当前候选依赖人工标签、评价信息被用于再次选模、预训练来源与测试重叠或分组不独立时，隔离受影响结论。必要时需新的独立参照；不能靠重命名文件解除暴露。
- **影响：** plan补充科学证据边界；roadmap/TODO改为无人工标签的发现→结果冻结→最终验证。独立核对先检查现有资产与入口，不自动重跑旧TLS监督训练或启动GPU；HF假设尚未因路线变化获得新证据。

## D-253｜2026-10-18｜独立核对后限定机器候选与旧评价用途
- **背景：** 按用户要求完成无人工TLS标签读取的产物/划分/代码核对及合成反例；证据见 `infra/tls_discovery_audit_20261018/` 与 [审计报告](TLS_DISCOVERY_AUDIT_20261018.md)。
- **理由／选择：** MISC39的sidecar为complete=false、16/578位置，旧run_summary的74完成不能作为逐片完成证据；当前HookNet完整候选输入只取其余73片，74片图像池保留，MISC39记PARTIAL而不删除。6902是跨阈值/通道记录，不是唯一TLS数。分组v1重放一致且SHA不变，但未知组不构成已证独立患者。
- **代价／结论收窄：** 原9/58名义p计数保留为旧方法描述；平移覆盖算子失配使其不能确认空间显著性，弥散/低RNA一致性不证明假阳性。TLS通道语义仍需训练标签映射来源，嵌套不能独立证明通道生物身份。人工标签映射和旧评价的实现缺陷升级为最终验证硬阻塞，不再宣称标签一到即能出可信报告。
- **复查条件：** 只有新版本完成记录和覆盖核查通过后才纳回MISC39；旧评价修复后须以同分、覆盖区、冲突、尺度及暴露元数据的合成反例验证，再冻结最终规则。真实标签在发现冻结前不读取、不用于修正算法。
- **影响：** 保留原代码/预测/表/划分作为历史证据；本轮不重跑推理、不覆盖结果。主线改为当前合同约束下的无人工标签探索，最终结构验证另行解封；HF未因此被确认或否证。

## D-254｜2026-10-18｜先运行四目标、两种可见集合的无人工标签恢复
- **背景：** 用户授权“开始吧”。本地74片冻结池同时具备表达与坐标；按基因覆盖检查，TENX23/TENX28的靶向面板不能覆盖四个预定目标。人工标签继续封存。存储仍按用户本任务明确指定的200GB（D-251）执行，通用AGENTS文本中的1.2TB不覆盖该专项授权。
- **理由／选择：** 计算前冻结 `infra/field_discovery_20261018/contract.json`：B-axis、氧化磷酸化、IFN-γ、EMT四个固定分子代理；5与15个可见模块；隐藏完整目标并从全部输入及其归一化/QC中删除目标基因。复用原五折和单位，不重分样本；资格只用基因覆盖、可读性和坐标，不用目标强度或人工TLS。所有留出点保存预测。
- **排除的方案：** 首轮不用人工监督TLS、不把HookNet作真值、不扩全GO扫描、不下载数据、不调K。先用同点、含无空间排列的整片背景、简单平滑与同容量邻域模型，避免只与弱均值模型比较。固定两种容量、训练子采样和一个重复种子，不按结果追加搜索。
- **代价／边界：** 原始HD文件已经聚合，metadata的native间距不能当聚合后的坐标间距；用每片实测最近邻间距归一化，只解释相对邻域，不声称统一物理尺度或超分辨率。未知患者身份、来源重复和已用于开发的数据使结果仅属内部探索；分子代理可恢复不等于TLS结构正确，也不确认HF-03。
- **判据／失效条件：** 合同预定相对每单位最强普通基线的保守门槛、失败保留、去均值误差和重复子采样稳定性，不生成空间p值。若输入审计发现目标复制、单位跨折、人工标签读取、余量不足200GB或目标定义/坐标不可读，隔离受影响运行。跨来源/患者独立确认、学习字典与最终人工验证均另列未完成，不冒充本轮已经检验。
- **影响：** 本轮服务HF-01/HF-02的有限恢复比较；HF-04只描述现有单位/来源的适用范围。结果尚未产生，不能提前判定任何HF成立。

## D-255｜2026-10-18｜身份分层采用已有资格表，保留首轮模型与预测
- **背景：** 首轮脚本的 `patient_known` 实际只检查metadata的patient字符串是否非空，得到10个单位；TENX152/TENX154虽均写Patient 1，既有资格表已因样本标签冲突将其列unknown，不能升级为两个已确认患者。
- **理由／选择：** `infra/field_discovery_review_20261018/` 按原资格表的known_within_study且patient_key与冻结unit相符才确认为研究内已知患者，明确8个已知单位、57个未知单位；输出修正的逐片身份、评价表及身份分层。原v1字段仅保留“原始患者字符串存在”含义，不能用于患者独立性主张。
- **代价／边界：** 首轮分组、样本资格、训练权重及模型输入均不使用该布尔字段，整体分子指标与候选门槛不受影响；原分组继续冻结，未知关系不擅自并成同一患者。历史v1产物保留哈希，当前身份解释与修正表必须一并引用。此次核对是保存产物的独立指标重算，不冒充独立重训或新验证队列。
- **失效条件：** 后续若有证据证明未知单元跨折同患者，对相应独立性结论隔离并重新分组评价；现在不把未知当独立。未来复用必须采用修正身份口径，不能从首轮非空字符串字段生成患者N。
- **影响：** HF-04仍无患者/来源独立确认；这项元数据修正不改变HF-01/HF-02的整体结果。人工TLS标签继续封存，未参与纠错。

## D-256｜2026-10-18｜限定背景/残差诊断与等基因预算学习字典探针
- **背景：** 用户要求继续。首轮整片背景模型弱于同点模型，邻域补回部分损失仍不胜强基线；固定字典的失败不能区分表示限制和模型失配。
- **选择与理由：** 新运行写入 `infra/field_followup_20261018/`，拟合前冻结合同。沿用72片/65单位及原五折、四目标和完整通道隐藏。先核对保存模型的训练/留出误差，比较去除整片背景的局部模型；用训练单位内三折交叉预测形成残差，分别用同点和邻域修正，系数固定0.5。外层测试不参与残差训练。
- **学习表示：** 各目标仅使用原full15可见基因并集、剔除目标基因，再取72片var身份的交集，所有对照共享相同基因、归一化分母和QC。训练折等单位取样拟合标准化PCA，固定8个带正负权重的分子组合；比较固定模块、PCA同点/局部模型与原始基因Ridge。此探针检验数据驱动表示，不把有符号轴叫非负TLS实例或独立生物身份；暂不并列搜索NMF/K/多种字典。
- **判定与代价：** 保存全留出预测、基因载荷、模型及失败；按合同评价空间/表示增量和跨折载荷匹配。跨折训练重叠，载荷稳定不等于独立重复。共同基因交集仅用检测面板身份，新的固定模块对照必须重算，不能与信息预算较大的旧分数直接判胜负。普通基线继续包含同点/平滑与原始基因模型。
- **失效条件与边界：** 原始目标进入输入/归一化、字典读到测试表达、内层单位跨折、共同模块少于10基因或存储余量不足200GB则停止受影响计算。身份采用D-255修正表。全部结果是已见来源上的后续开发证据；人工TLS、独立迁移和校准不确定性仍未测，不因新模型过门槛升级确认性结论。
- **影响：** 诊断服务HF-01，等预算学习表示探针服务HF-03；不修改顶层科学目标。

## D-257｜2026-10-18｜重复基因符号按实际输入范围检查
- **背景：** 新准备阶段在TENX128检测到TBCE、HSPA14、TMSB15B重复符号并停止；三者均不属于本轮模块或目标，尚未训练模型。
- **选择／理由：** 保留逐片重复符号清单，对实际可见基因和目标符号要求唯一；无关列的重复不排除整片。可见集合不能以字典索引静默取重复项的最后一列。共同可见基因、目标、样本和模型规则不因该修复改变。
- **失效条件：** 重复符号落入可见集合或目标时停止对应准备，需有来源的身份处理后才能计算，不擅自求和或去重。原失败日志保留。

## D-258｜2026-10-18｜保留局部改进，将原始基因回归/平滑纳入后续主基线
- **证据：** D-256完整五折已运行并冻结。去整片背景的局部模型相对原最强模块级普通模型，B-axis/IFN-γ平均RMSE降低3.18%/5.20%；不能再把首轮结果概括成所有空间输入均无收益。相同基因预算下，PCA局部模型仍比原始基因Ridge或其固定平滑高19.93%～60.49%误差。全部预定候选均未通过保守单位门槛。
- **选择／理由：** 保存局部正向结果和失败，后续主基线包括原始可见基因Ridge、固定简单平滑与同点/模块对照。PCA改善固定模块不能单独证明新生物身份或空间关系；不因均值改善而放松原单位门槛，也不按当前结果追加K搜索。
- **边界／失效条件：** 8维PCA是限定基因集合内的线性有符号表示，非所有学习字典的代表；跨折载荷稳定不是独立重复。原始基因对照适用于可见基因可获取的条件；若实际只能获取汇总场值，必须在相应信息约束下重新定义合法对照，不能跨信息入口比较。
- **后续：** 在现有本地、可追溯来源与更少可见信息的条件下检验边界，再决定复杂空间分支的价值。独立来源/患者验证未运行，人工TLS不解封。此结论影响HF-01/HF-03的当前实现判断，不否定全部条件恢复问题。

## D-259｜2026-10-18｜来源分组留出与5/15输入的强基线比较
- **授权：** 用户再次明确沿用200GB存储余量；CPU运行，人工TLS继续封存，不新增外部数据。
- **拟合前固定：** 沿用四个目标、72片/65单位、共同检测基因身份与原五折，比较small5/full15两种输入。原始基因Ridge、固定0.5邻域平滑、模块同点HGB、模块局部HGB及训练均值共享可见基因集合；每种输入单独按可见基因归一化，目标基因不进入分母或QC。不追加PCA/K或超参数搜索。原始表达重新生成small5，full15复用已冻结缓存。
- **来源资格：** 53片/47单位具备明确研究链接，按6个研究分别留出，训练也仅使用这些已确认研究来源；19片来源不足标为NOT_TESTABLE_SOURCE_IDENTITY，仍参加全72片内部比较。肺癌6片均属来源不足。肝脏和乳腺各只有一个明确来源，留出时无同器官训练，单列为器官与来源联合迁移。研究分组无已知单位跨组；研究链接不保证公共材料没有跨研究复用，既有开发暴露不能追认为新盲验。
- **判据：** 保存全量留出预测和模型，按单位等权RMSE、中心化RMSE比较每个基线，同时列出每来源结果，避免34片大来源掩盖异质性。局部候选须对每个普通基线平均误差改善至少3%、中心化误差改善为正、至少60%单位改善；这只形成开发候选，不构成患者独立确认。不以测试最优选择可部署模型。
- **失效条件：** 目标泄漏、已知单位跨训练测试、输入哈希变化、参与输入的基因符号不唯一或资源余量不足时停止对应计算。来源资格不足限制来源验证，不删除其内部结果。结果服务HF-01/HF-02/HF-04，不改变科学问题。

## D-260｜2026-10-18｜来源留出后保留普通恢复基线，限定空间模型失败范围
- **证据：** D-259的11个划分、1000份预测和5000条误差已完成；16个目标/输入/划分条件均未同时通过所有普通基线门槛。5→15输入使Ridge内部平均误差降低10.83%～25.66%，来源留出降低3.41%～23.43%；并非高阶交互证明。相同样本按六来源等权，full15 Ridge误差相对内部增加3.07%～36.50%，但乳腺/肝脏同时跨器官，不能当纯来源效应。
- **选择／理由：** 当前模块局部HGB不进入空间关系候选库；保存全部普通模型恢复结果及来源边界，不追加其超参数/字典数搜索。后续若检验空间增量，优先在同一原始基因表示下构建有限空间对照，区分模块压缩损失与空间信息缺失；不能据当前比较停止所有空间路线。
- **范围／失效条件：** 本轮是已有开发样本的新分组，患者与来源资格限制继续有效（I-044）；不升级为独立四癌验证或TLS发现。后续同表示空间模型在预定留出与强基线下出现可靠增量时，可更新当前实现判断。人工TLS继续封存。
- **复现说明：** 新旧内部Ridge抽样/标准化一致，预测最大差4.60e-5、单片RMSE最大差4.12e-6；原严格数值一致性断言失败日志保留。后续审计分开记录数据/分组正确性与数值复现，未改模型结果或科学判据，不声称逐位一致；底层差异原因未独立定位。
- **数值精度更正：** 上述来源留出Ridge增量下界按完整产物四舍五入应为3.40%，其余范围与判断不变；精确逐目标值见结果表。

## D-261｜2026-10-19｜用户指定Paseo/Pi实现与独立审核分工
- **授权与选择：** 用户明确要求代码工作交由Paseo管理的Pi代理 `opencode-go/deepseek-v4.1-flash`、推理 `max`；独立代码审核交由Pi代理 `xai/grok-4.7`、推理 `xhigh`。主代理负责科学要求、调度、约束记录与结果分析，不自行实现计算代码或替代审核。
- **执行边界：** 每次核对实际provider/model/thinking及工作目录，记录agent ID。实现和审核分离，审核绑定代码/合同哈希；影响结果的问题修复并复核前不用于科学结论。此要求不能理解成绝对无缺陷保证。
- **先行任务：** 先准备四目标实测场资格和现有预测生物学核验的实现/合同，经独立审核后再运行全量描述性分析；不预设空间形状。旧3%规则保留历史，不能代替当前生物学价值判断。旧模型完整目标基因训练后，事后拆分基因不作为独立验证。
- **资源与失效条件：** 人工TLS封存、本地数据、CPU、200GB下限不变。指定provider/model/档位不可用则报告阻塞，不静默换模型或由主代理代写。当前已确认宿主Paseo运行，但Pi provider配置disabled；启用前未创建执行/审核代理，不声称已派工。
- **影响：** 分工配置是基础设施，对HF假设无新增实证。任务书保存在 `infra/paseo_field_review_20261019/`；实际运行状态以同目录登记为准。
- **启用与派工实录：** 用户随后明确允许启用Pi。仅将宿主Paseo的 `agents.providers.pi.enabled` 从false改为true，热加载成功且无须重启。实时模型列表和运行检查均确认两个精确模型/档位可用；执行代理 `41f406bd-bd3a-44e6-8dd8-2fa5c0ea93a6` 为Pi/`opencode-go/deepseek-v4.1-flash`/max，审核代理 `a5bcb012-6618-488e-8f88-f700773fd4c4` 为Pi/`xai/grok-4.7`/xhigh，工作目录均为本项目。尚未完成代码审核时不得宣称审核通过。

## D-262｜2026-10-19｜Paseo调度改用已验证MCP，等待不做高频轮询
- **核实：** 用户修复后当前进程可读取PASEO_PASSWORD；不输出密钥，直接调用MCP initialize/get_agent_status/send_agent_prompt成功。当前会话原生工具清单仍未注入Paseo，不能把协议可用说成自动完成回调已接通。
- **运行事实：** 指定DeepSeek上轮报错为上游模型不可用，已有代码继续保留。用相同模型/max档位通过MCP恢复；审核代理的合同预审已完成，代码审查未完成。
- **等待选择：** 本机MCP send_agent_prompt(background=false)约30秒后返回running，不是任务失败。未具备agent-scoped回调时使用一次无超时的paseo wait等待结束，不反复发状态查询或短超时轮询；用户主动询问或需要干预时可做一次检查。
- **边界：** 仍按D-261分离实现与审核，当前修改后的代码不得沿用旧哈希检查结果。此节点是基础设施，对HF无新增实证。

## D-263｜2026-10-19｜独立审核拒绝后修复汇总语义并冻结交接版本
- **依据：** Grok首轮代码审核NOT_ACCEPTED，R-01至R-04能导致汇总中断、单位权重变化、外环指标偏差或器官外推误读。自测通过不作为全量放行依据，见I-045。
- **选择：** 交回原DeepSeek代理修复，要求真实表模式的汇总贯通检查。一个多切片单位若某指标部分不可算，保留部分单位及原因，不冒充完整单位进入该指标的完整单位等权均值；全空分层保留记录和分母。混合来源的外推状态不得取任意首行代表整臂。
- **交接：** 代码、协议、合同、检查与handoff哈希必须一致；交付后停止修改直至审核反馈。早期快照或旧版验收不能替新版本背书。首轮审核记录另存，不覆盖。
- **失效条件／影响：** 新版本未关闭影响结果的问题前，全量描述性运行保持硬阻塞。本次工程修复本身不产生HF-01/02/04科学证据；旧模型产物继续冻结。

## D-264｜2026-10-19｜独立复核关闭4项缺陷，按绑定版本进入全量描述性运行
- **证据：** Grok round2 `ACCEPTED_FOR_DESCRIPTIVE_RUN_ONLY`，未关闭结果相关问题为空。脚本52d11904…、协议6062a19d…、合同a71eb3d8…与交接/磁盘一致；独立小反例及5片实际汇总贯通产物关闭R-01至R-04，1550/0自测仅作辅助，不替代审核。原不通过记录保留。
- **执行：** 仅按上述哈希运行prepare/run/summarize，全72片及两臂1000份冻结预测，CPU、本地数据、人工TLS封存和200GB余量不变。实际运行若需改代码，旧许可失效并须重新复核，不静默续用。
- **含义：** 解除I-045对该版本的执行硬阻塞。新输出仍是已开发材料的描述性再分析，数学可算不等于生物学参照资格；尚无新的HF实证结论，不把审批本身作为科学节点完成。

## D-265｜2026-10-19｜使用 Git 管理审核版本，取消重复代码哈希记账

- **依据：** 用户明确要求沿用项目既有 Git 版本管理，不另行维护代码哈希。
- **选择：** 后续执行交接与独立审核以 Git 提交及必要的工作区 diff 标明范围；取消手工代码哈希和重复版本清单要求，替代 D-263、D-264 中相关交接要求。未提交改动如实记录，不冒充已经提交或审核。
- **边界：** 保留历史审核记录与既有产物，不为此次流程调整修改正在运行的计算代码。结果相关变更仍交由指定执行代理修改、指定审核代理复核。
- **影响：** 仅调整版本管理流程，不改变科学合同和结果。

## D-266｜2026-10-19｜关闭单位重复计权后再解释全量描述性输出

- **依据：** 全量交接暴露 `singleton::TENX13` 的两张切片被不同标题/source_key拆成两个单位组，违反完整单位等权合同；已由指定DeepSeek代理在 Git `5fee80f` 修复，内部恢复65单位，来源47单位保持不变。
- **选择：** 复用逐片预测，仅重做受影响汇总；完整单位由既定unit_id确定，来源等元数据冲突如实保留mixed，不推断独立患者。交由指定Grok代理复核Git提交及实际产物，复核前不使用受影响汇总作科学结论。保留旧版对照，不为此重训或重跑完整预测。
- **解释边界：** 先形成四目标的分子支持、位置、范围、局部细节和关系恢复表，不引入生物学PASS百分比。当前细节分量自相关不充当原场自相关；几乎全不可评的半衰尺度不充当生物尺度失败证据。
- **失效条件：** 若复核发现其他结果相关问题，仍由原执行代理修复并复核，未关闭产物不进入结论。工程修复本身不支持HF；科学影响须依据复核后的结果单列。

## D-267｜2026-10-19｜复核通过后按具体空间特征解释恢复能力

- **证据：** Grok对Git `5fee80f` 的复核允许描述性结果分析，单位重复计权已关闭；交接旧分母和有效基因中位数勘误，不改计算产物。72片/65单位、来源53片/47单位的四目标能力与边界见 `FIELD_BIOLOGICAL_REVIEW_RESULTS_20261019.md`。
- **判断：** 普通模型已恢复相当部分整体图样，B轴与IFNG的局部细节和高值区域仍较弱；平滑改善部分重合度却可能压低起伏、扩大范围或增加质心偏移。5→15输入在图样/重合/细节指标上均有有限收益，支持HF-02的有限输入增量，不证明高阶组合。未建立HF-01同输入空间增量；HF-04获得描述性边界，未获得校准/独立确认。HF-03/05未测。
- **下一步选择：** 优先围绕B轴/IFNG的已知缺口设计有限的同原始基因输入空间对照，OXPHOS作为普通模型较强的参照、EMT作为范围与位置取舍参照。暂不扩大模型搜索或外部数据，不将分位重合度定义成生物真阳性门槛。具体新计算须先落盘比较合同。
- **失效条件：** 若提出原场尺度、独立身份或TLS真值主张，I-047按升级条件硬阻塞该主张；若新对照不在同一输入预算或使用开发集调出的规则，不能据此宣称空间增量或盲验。

## D-268｜2026-10-19｜固定同原始基因输入的邻域回归对照

- **依据：** 用户继续授权执行D-267；B轴/IFNG局部细节与高值区域仍有缺口，先排除输入表示差异而非扩大模型搜索。
- **选择：** 四目标、full15、原内部5折及来源6划分；固定alpha=10。唯一候选为同点/近邻/外环原始可见基因特征的Ridge，比较同点、固定平滑、整片背景、固定置换图对照。三个标准化块按1/sqrt(3)缩放，小检查验证复制同点输入时不改变有效惩罚。具体合同见 `FIELD_RAW_GENE_SPATIAL_CONTRACT_20261019.md`。
- **理由：** 固定线性模型便于直接检验邻域信息，避免将模块压缩、额外基因、背景均值或正则强度变化误称空间增量；不启动复杂网络或超参数扫描。块缩放不保证不同特征秩完全匹配，置换仅为机制破坏，不作显著性。
- **判定：** 局部细节、高值区域、位置、范围绝对偏差与幅度分别解释，不设置RMSE百分比成功门槛；主要服务HF-01，HF-04限开发范围。失败只限制本次固定线性对照。
- **失效条件：** 泄漏、分组/归一化/条码错误、求解不收敛或资源无法覆盖完整范围时停在对应阶段并报告，修复后指定审核代理复核；不静默删目标/划分。

## D-269｜2026-10-19｜训练前明确输入复用与评价口径

- **依据：** Grok合同预审B-01至B-06指出可导致实现不同义的边界；尚未运行新全量评价。
- **选择：** D-268合同补充指定复用原full15矩阵/QC/划分/抽样，真实/背景/置换均从全片可见X构造，正式求解流程验证复制块等价，固定0.5近邻输出平滑；绑定Git5fee80f评价定义及逐片绝对范围误差。内部和来源评价均不用于调参或选择种子。
- **数值约束：** 根据已记录旧Ridge复算误差量级，预测最大绝对差容差1e-4、每片RMSE差1e-5，超出停止查因，不以新评价结果放宽。具体科学合同补充节和预审记录可追溯。
- **失效条件：** 无法证明同输入/同分组、数值复现或有效惩罚等价时保持运行阻塞，交原执行代理修复和Grok复核；这些工程约束不构成生物学成功标准。

## D-270｜2026-10-19｜统一器官/平台辅助输入的正则惩罚

- **依据：** Grok对Git7ba553f审核NOT_ACCEPTED。正式三块设计把仅出现一次的器官/平台指示列E一起乘1/sqrt(3)，相对同点模型改变其惩罚；复制整份点设计的检查未覆盖真实实现。
- **选择：** X限定基因/QC，三块X按既定缩放，E仅拼接一次且不参与块缩放，标准化/惩罚与point一致。原DeepSeek修复正式设计及对应复制块检查，Grok复核；全量继续阻塞。
- **记录勘误：** 复制块数值界限确实由1e-6放宽到至少1e-5，不以“不是放宽”掩盖；审核认可其数值分辨率依据而非事先固定性。修复后界限固定1e-5，保留求解器地板和错缩放负对照，超出不自动放宽。原冻结point界限不变。交接预计5000行勘误为2500。
- **失效条件：** 同输入或有效惩罚仍不一致、数值检查失败则不能作空间增量比较，修复和复核前不运行全量。本次无新增HF证据。

## D-271｜2026-10-19｜停止定时读取等待会话，使用MCP回复

- **依据：** 用户明确要求不要轮询，使用MCP等待subagent回复。
- **选择：** 停止周期性get_agent_status、write_stdin读取CLI等待会话，以及短超时MCP循环；只使用一次MCP同步回复或真正的完成回调。已派发修复任务继续，不停止代理。
- **实测限制：** 当前工具目录无原生Paseo工具；直接MCP的send_agent_prompt(background=false)在30秒后返回running超时，notifyOnFinish按工具schema仅支持agent-scoped调用。当前Codex会话未接入该回调，因此不能承诺自动唤醒继续审核/全量运行；不能用轮询掩盖这一限制。
- **影响：** 不改变实验合同。R-01修复由指定代理继续，全量仍须独立复核放行；不能把等待或运行中视为完成。

## D-272｜2026-10-19｜复核关闭辅助输入问题，启动完整空间对照

- **证据：** Grok第二轮 `ACCEPTED_FOR_FROZEN_DESCRIPTIVE_RUN`，允许Git `0e6e5f5` 全量运行，R-01关闭，结果相关未关闭问题为空。原冻结point复现界限及固定等价界限保持；审核不等于科学成功。
- **执行：** 按既定full_run_task由原Pi/DeepSeek max执行四目标×11划分、full15的fit/evaluate/summarize；预期500预测包、2500逐片模型行。原数据、参数、分组、标签封存及CPU/200GB约束不变。
- **边界：** 只确认运行许可，HF尚无新结果。运行若触发数值界限或发现代码问题，停止相关流程、保存失败并交指定代理修复复核；不静默跳过。等待遵守D-271，不轮询。

## D-273｜2026-10-19｜基线复现超界后停止评价并独立核查数值原因

- **证据：** Git1df03a3全量交接：44/44拟合、500预测包完成，11/44冻结基线复现超出原数值界限，审计返回1；evaluate/summarize未运行。
- **选择：** 保留全部训练产物和失败，原容差不改；由Grok核查复现审计、设计/抽样/权重/列与数值环境一致性，给出必要修复及重算范围，再交DeepSeek执行。
- **边界：** 执行者提出环境敏感性解释，但受控证据是否充分仍待审核；同环境不保证模型差值无数值偏差。本轮不能产生HF-01的新性能结论。
- **解除条件：** 独立审核关闭结果相关问题并明确允许后续评价，必要修复与重算完成；不以训练完成或数组有限替代。

## D-274｜2026-10-19｜先用单个受控探针验证线程原因

- **依据：** Grok数值审核确认11/44基线复现失败，认为BLAS原因仅为推测，N-01仍阻断。
- **选择：** 仅source_0/EMT，保存脚本，固定训练设计/行序/权重/列mask/dtype/参数，对比实际1与8线程及冻结/存盘预测；交Grok判断，不先改正式启动器或重算全部44组合。
- **失效条件：** 若线程不能复现差异，不宣称原因已定，转向最小设计/保存路径定位；若确认原因，仍需明确修复及复核，再确定重算。原数值界限、标签封存、CPU/200GB不变。


## D-275｜2026-10-19｜解除旧预测逐位复现对描述性评价的阻塞，改查结论敏感性

用户指出线程一致性排查偏离科学问题，并明确要求继续推进。旧冻结预测差异仍按原1e-4逐点/1e-5逐片RMSE界限记录失败，不事后改为通过；该工程界限没有被证明对应生物学或方法比较有效性，故不再作为读取现有500包、计算描述性空间指标的先决条件。本轮使用同批运行的point/smooth/candidate/background/permuted作配对比较，不将旧轮次成绩混入主比较，不更换线程/求解器择优、不全量重拟合。

DeepSeek负责执行已有evaluate/summarize及必要的最小代码修改；Grok独立审核结果链与数值误差对结论的实际影响。可先生成未审核描述性结果，但审核未关闭的结果不能支持科学结论。停止扩展线程扫描；如主要优势极小、排序或空间指标对数值扰动敏感，只将对应结论标为未确定，再提出最小必要验证，不自动重开全量工程排查。不能把一个point案例的误差上界推广到所有模型，也不能因共同运行环境假定误差抵消。

失效条件：出现泄漏、条码/分组/归一化错误、非有限预测或其他足以影响比较的具体证据，应阻止对应科学解释并修复。I-049保留为结果解释约束；其原全局评价/汇总硬阻塞由本决定取代。科学目标与HF-01/HF-04判据不变。
## D-276｜2026-10-19｜先核查四个代理场定义，再决定候选场拆分

用户要求沿已讨论路线推进。现有真实邻域线性对照没有提供稳定的空间排列增量，且 B_AXIS 的有效支持基因数中位仅约 2.87、首位贡献中位约 0.535；这说明代理定义可能是限制因素，但尚未证明。因此先按 [代理场检查合同](FIELD_PROXY_AUDIT_CONTRACT_20261019.md) 检查四集合的片内基因空间协同、主导性、独立基因子集一致性与跨单位／来源重复性，再决定是否在训练单位内提出少量候选子场。原四代理与单基因图保留为参照，留出数据不用于挑选场定义。当前不使用人工 TLS 标签、不以旧恢复分数选新场。

排除立即扩大空间模型或全基因无约束搜索：前者无法分辨代理混杂与关系模型不足，后者在已开发的小来源集合中易于把技术梯度或测试表现选成“场”。失效条件：若检查显示预设代理已稳定且无明显混合，应把重心转回场间空间关系表示；若可测基因不足，停止对应集合拆分并报告不可评。该节点首先服务 HF-03；代理可靠性本身不证明 HF-01。

## D-277｜2026-10-19｜关闭代理场检查预审指出的评价口径缺口

Grok 对 D-276 合同的预审 `infra/paseo_proxy_audit_20261019/review/preflight.json` 提出 P-01～P-05，认为技术梯度、原读数可比性、全数据描述与后续留出隔离、单位汇总和可测基因下限会改变解释，第一段暂不放行。故在任何第一段指标计算前补充 [合同执行口径](FIELD_PROXY_AUDIT_CONTRACT_20261019.md#预审后的执行口径grok-p-01p-05d-277)：同一全基因文库归一化；点级非目标深度/检出率敏感性对照；每子集至少3个可测基因；10次固定种子拆分及按检出/表达四分位匹配的目标外对照；完整既定单位等权；第一段不锁第二段基因和候选数。技术控制可移除真实生物信号，不能用它单独判生物/伪影。

这些固定值是可评性与工程预算规则，不是“场存在”的成功阈值。若固定对照候选不足，保留不可评，不放宽匹配；若同一归一化/条码顺序无法复现原读数，停止图间比较并交执行代理修复。经 Grok 对补充口径复核放行后才完整运行第一段；第二段仍待第一段结果与另行范围判断。此预审不产生 HF-03 新证据。

## D-278｜2026-10-19｜四目标统一进行训练内两分量候选场探针

第一段代理场检查已在 Git `d184488` 完成，R-01 对照抽样问题由 `a46b50a` 修复并经 Grok 在 `51703bf` 关闭，允许描述性结果。B_AXIS 原始计数加权有效基因数中位约 2.87、首位份额 0.535，但本轮相同读出尺度的留一图相关中位 0.985；两者口径不同，既不能说 B_AXIS 只是一个基因，也不能据此确认它是单一干净的场。IFNG 两半基因图有部分同向；EMT/OXPHOS 匹配对照分别只在 2/72、7/72 片可评，不能用其对照差挑选目标。

为避免看完全数据描述后挑有利目标，下一有限探针四目标均运行，固定 K=2，构场基因 A／核验基因 B 按基因名与固定种子在表达读取前拆分；候选定义与全部可学习参数只在既定训练单位拟合，留出 B 图只作评价。同 A 基因预算比较简单均值、K=1、单基因和全 A 基因普通 Ridge；逐基因、逐单位、逐来源保留失败与不可评。合同见 [两分量候选场探针](FIELD_CANDIDATE_SUBFIELD_PROBE_20261019.md)。这是 HF-03 的候选身份探索，不是 HF-01 场间空间恢复；使用过本批材料的目标范围选择仍属于开发，不能称独立外部确认。

排除只做 B_AXIS 或根据第一段匹配对照差删掉 EMT/OXPHOS，因为会产生看后选目标且后两者对照覆盖不足。固定 K=2 仅测试最小拆分，不是最优 K 搜索。失效条件：Grok 预审发现 A/B 泄漏、训练/测试分组不独立、技术梯度或强参照不可比，应先修订合同/实现再运行；若候选不稳定或不优于普通参照，如实限制本固定表示，不自动扩 K 或换方法找好看分数。

## D-279｜2026-10-19｜候选场预审后消除核验基因的分母入口

Grok 在 `infra/paseo_candidate_subfield_20261019/review/preflight.json` 提出 S-01～S-06，认为第二段拟合前存在会改变解释的缺口，尤其原全基因文库分母含核验 B 基因，B 可代数进入构场 A 图。第二段尚未运行。故在任何候选拟合前修订 [两分量合同](FIELD_CANDIDATE_SUBFIELD_PROBE_20261019.md#预审后固定口径s-01s-06d-279)：A/B按排序基因名与固定NumPy种子精确拆分；A和B统一用排除本目标全部基因的外部文库分母，明确与旧代理图尺度不同；训练单位冻结可测 A、抽样和缩放；K=2非负模型及载荷一一对齐、不可评条件固定；同A信息入口的普通Ridge强参照与B稀疏/q75规则固定；技术协变量排除整个目标集合。

这是排除核验泄漏和比较不公平的必要选择，不是根据留出效果调参。共同外部文库仍会带来技术相关，需报告而不能称独立真值。若外部文库不足、训练可测基因不足、模型不收敛或Grok复核仍不放行，保留 `NOT_TESTABLE`/阻塞，不更改K、基因拆分或分母以求有利结果。预审及可行性检查不产生HF-03的新证据。

## D-280｜2026-10-19｜独立代码审核 Grok 推理档位改为 high

用户明确将审计用Pi `xai/grok-4.7` 的推理档位由 `xhigh` 改为 `high`。主代理已通过 Paseo MCP 更新现有审核代理 `a5bcb012-6618-488e-8f88-f700773fd4c4`，随后确认 `effectiveThinkingOptionId=high`、model=`xai/grok-4.7`、cwd=`/home/huyudi/013_spatial`。执行代理仍为 `opencode-go/deepseek-v4.1-flash` `max`。后续审核均按 `high` 记录实际档位；已完成的 `xhigh` 审核保留历史事实，不追溯重写。若运行时该模型/档位不可用，按项目约束报告阻塞，不静默替代。

## D-281｜2026-10-19｜冻结固定 K=2 探针，暂不按留出分数扩大 K 搜索

Git `054a16b` 完成四目标×11划分候选场探针；Grok `high` 在 `db75404` 独立复核，允许描述性结果，无必修代码问题，但不允许独立场科学主张。固定K=2对留出核验B图的配对Pearson中位，相对A均值为B轴−0.014、EMT+0.011、IFNG+0.002、OXPHOS−0.001；只有EMT有有限改善，且相对全A Ridge仅+0.003。OXPHOS相对全A Ridge+0.022但不如A均值，且技术因素相关最高。折内重抽样载荷稳定不能代替跨单位身份确认。完整解释见 [结果](FIELD_CANDIDATE_SUBFIELD_RESULTS_20261019.md)。

选择冻结本轮K=2代码与产物，不因看到EMT/OXPHOS小差距就继续扫描K、重分A/B或换归一化挑更好成绩。后续若要探索新的场表示，先提出能区分技术梯度、组成混合和场间关系的证据与同信息强参照，再设新合同。失效条件：若独立审核发现实际产物有影响结果的错误，或出现真正新的独立数据/分子证据，可按新版本重评；不能把本次结果追认成独立生物真值。该节点未给HF-03新增正支持，未测试HF-01。
