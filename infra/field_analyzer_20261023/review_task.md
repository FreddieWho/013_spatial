# FA-13 独立完成复核任务

- 派发日期：2026-10-23。
- 实现版本：`e7affab`；计算前合同：`3a07f06`。
- 实现代理：Pi `opencode-go/deepseek-v4.1-flash` / `max`，`41f406bd-bd3a-44e6-8dd8-2fa5c0ea93a6`。
- 独立审核代理：Pi `xai/grok-4.7` / `xhigh`，`82d66e60-3f71-4310-8840-14777ad3151f`。
- 工作目录：`/home/huyudi/013_spatial`。
- 经 Paseo MCP 派发；完成等待采用一次事件订阅，不循环查询状态。

请完成本轮 FA-13 独立完成复核。工作目录 /home/huyudi/013_spatial。实际配置必须 Pi / xai/grok-4.7 / xhigh；不符则报告阻塞，不换档。本轮实现代理是独立 Pi DeepSeek max，主代理负责科学合同与结果解释。先读 STATUS.md、docs/plan.md、docs/FA13_CONDITIONAL_RELATION_20261023.md 和 D-284；合同在计算前由 3a07f06 冻结。

你不是独占仓库，不得回退或修改实现、现行主控文档或他人文件。你只拥有 infra/field_analyzer_20261023/review/ 下审查 JSON/说明。请用当前 Git 提交及必要 diff 指明实现版本，无需手工代码哈希。不得读封存人工 TLS 标签，不能下载数据/权重，CPU，存储至少留 200 GB。

审查新增 field_analyzer/conditional_relations.py、scripts/field_analyzer_conditional_20261023.py、tests/test_fa_conditional_relations.py 及 infra/field_analyzer_20261023/ 的等价性核对、结果和资产清单。只做影响结果的必要检查，不默认跑全库测试。

重点：
1. 比较是否确为 M0 自身邻域 vs M1 自身+跨场邻域，固定 EMT←IFN 与 IFN←OXPHOS，目标邻居只用特征半边，Q/训练行/共同评分点一致；正反两个半边都交付，不筛选较好结果。FA-05/FA-12/旧 X-WX-VX 的等价性核对是否准确。特别注明：FA-12 same_neighbour 用目标 witness 在 O 的值，新合同为了分子内部核验改用 feature 半边，这个输入预算差异必须在报告中明示；不能把新联合结果直接称为超过旧 FA-12 强参照或超过全基因可见预算下的最优模型。
2. 遮挡是否从原始观测开始？逐基因缩放、场评分、深度残差、标准化有没有在 Q 上拟合？隐藏 Q 的计数/深度/QC 能否改变输入或预测？不要只读测试名字，必要时重放一个真实切片的毒化，注明哪些检查是重跑、哪些仅审读存档。离线 witness 目标可以用于评价，但不能回流输入。
3. 所有场的共享基因与 witness 并集排除、背景文库、交换两半是否隔离；邻域不含自身，不从隐藏点借表达，不将缺失当观测零。坐标/条码/行顺序和 metadata 独立单位是否对应，留一单位是否把同单位切片一起留出。
4. 模型拟合是否只读训练单位，留出单位无监督校准是否仅用许可的 O；测试误差绝对和、均值、单位等权/逐片等权有没有混淆。按同一查询集核算主效果及两半/单位一致性。
5. 合成例子是否真正展示“cross 单独较弱而 joint 有增量”，无增量反例是否合理；不要将其当真实生物学确认。
6. 本地 Xenium/Visium HD 数据资产检查是否仅基于 metadata/轴；ROI、连续切片、患者没有互相冒充，未知保持未知。不能沿用旧 Base/Addon 同患者/连续切片叙述；不因有 COL1A1 等单个基因就宣称结构真值。
7. 无空间匹配零分布或独立结构证据时，即使条件增量为正也不能直接升级 AF-01，更不能声明桥/非局部连接/因果或患者级泛化。本轮随机 10% 整点留出也不是连续区域缺失实验，不能直接升级 AF-03。阴性同样仅限具体模型和观测条件。

本次交付有两项请特别核对：报告声称“≥6/8 是预先固定”，但主代理冻结合同并未写这个通过阈值，不能把执行者新增的描述性规则冒充合同预注册。报告另称 Xenium 基因轴为 10,756–16,915 个，但历史已解包结肠轴为 541；请用少量实际 gene_panel/features 交叉确认提取的是基因名而非 JSON 字符、非基因实体或混杂字段。说明归档可读范围和未完成资产检查，不要一概宣称资产核验完整。

输出 review/fa13_review.json：actual_runtime、git_commit、scope.read/ran/not_checked、findings（严重性、文件行/原因/影响/修复建议）、effect_checks、allow_node_conclusion、limitations。发现结果性错误则 allow=false，交回实现者修复，你不要代写代码。主代理会修复后再请你复核。最终简短说明主要发现和未检查范围。
