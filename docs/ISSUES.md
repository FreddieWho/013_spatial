# 问题日志：空间结构反演

本文档记录尚未解决、但当前未造成硬阻塞的问题和证据缺口。它不是任务清单：执行顺序属于 `roadmap.md`，历史技术理由属于 `decisions.md`。每个问题必须写明影响、当前证据、下一次判定点和升级为硬阻塞的条件；状态变化时更新原条目。

允许状态：`OPEN`、`RESOLVED`、`HARD_BLOCKED`、`WONT_FIX`。

## I-001｜validation 谱系缺少可审计的结构 GT（training 与三个公开沉积 validation 谱系已解除）

- **状态：** `OPEN`
- **发现日期：** 2026-07-31
- **影响节点：** R-04 完成判据（队列间一致性）、R-08、R-10，以及所有 internal/external validation 主张的确认性证据等级。
- **当前事实：** training 谱系（HTAN Vanderbilt CRC）已于 2026-07-31 解除（Heiser 42 个逐 spot 病理标注 CSV，42 个可审计 source、195 个确认性实例）。2026-08-03 经用户批准的三条公开沉积路径（D-037）再解除三个 validation 谱系：GSE175540（Meylan ccRCC，作者沉积逐 spot TLS 标注，18 个可审计 source、35 个 TLS 组件/18 患者，external）、TLS_VISIUM_USZ（Zenodo 14620362，h5ad `obs/ground_truth` 人工标签，8 个 source、108 个 TLS 组件/8 患者，external）、ST_CRC_CMS（Zenodo 7760264，病理医生逐 spot 类别，12 个 source、551 个肿瘤—基质边界组件/6 患者，internal）。确认性实例合计 889（TLS 147、TSB 742），可审计 GT source 合计 80，机器门禁保持 `PARTIAL_GT_READY`。剩余缺口按角色分层：internal validation 的 GSE274103 有 Loupe/HALO 标注（覆盖 LA/TLS、血管、坏死、侵袭前沿四类结构）存在但未沉积，只能向作者索取（lead contact Linghua Wang）——其独占剩余价值已降为未冻结结构（血管/坏死/侵袭前沿）的标注，索取邮件由必需降为可选；external validation 的 GSE211956 有 QuPath 几何未沉积，同样只能索取（可选）；internal validation 的 GSE226997 原论文无独立标注流程，没有任何可索取对象；discovery 角色的 GSE274557 已明确放弃索取，无主张损失。此前记录的其他事实不变：87 条 Table S4 TLS ID 仍只是汇总属性，19 个本地 STOmics 文件仍在 R-01 冻结核心之外。
- **污染与泄漏：** Table S2 的 TLS presence/count、Table S4 的 TLS ID/maturation/location/size 和 `repo/data/ST_*_maturation_location.csv` 的 Cluster/Location 都是目标或目标派生字段，进入输入会造成直接 target leakage；表达派生标签反过来验证表达模型会造成循环验证。GSE211956 的 sample title 编码 treatment response，属于 outcome leakage；patient/site/treatment/TNM/sample/file identifiers 均是 shortcut 风险。GSE274103/GSE274557 的官方 metadata 说明预处理已删除 acellular stromal-interface spots，因此肿瘤—基质边界还存在选择偏差。上述字段全部禁止作为普通模型输入。
- **当前影响：** R-04 至 R-07 的工作可在分子-only 输入与冻结外层切分内启动；R-04 的队列间方向一致性判据与 R-05 的两来源一致性判据恢复可执行，限四个有 GT 的谱系（HTAN Vanderbilt CRC、GSE175540、TLS_VISIUM_USZ、ST_CRC_CMS）。GSE226997、GSE274103、GSE211956 保持 GT-less，不得承担任何确认性主张；血管与坏死无公开 GT，第二结构以外的结构主张仍冻结。R-01 的身份 gate 仍有效。
- **下一判定点：** R-04 启动时确认队列间一致性的谱系组合（四个有 GT 的谱系）；若未来需要血管或坏死的确认性检验，再评估是否向 Wang 实验室索取 GSE274103 的 HALO/Loupe 血管/坏死标注（唯一一封仍有独占价值的邮件）。GSE211956 的 QuPath 索取仅在外部验证层需要第三独立谱系时再考虑。
- **解除条件：** 至少一个结构在至少一个独立 validation lineage 中具有内容寻址、可重放空间几何、physical-unit 映射、已知标注模态与生成方法、可执行输入排除规则——2026-08-03 起已对三个 validation 谱系满足；本条目对 GSE226997/GSE274103/GSE211956 三个冻结 GT-less 单元保持 OPEN。

## I-002｜serial-section 与配准相关 metadata 覆盖未知

- **状态：** `OPEN`
- **发现日期：** 2026-07-31
- **影响节点：** R-03、R-11、R-12，以及 H-06/H-07 的可测性。
- **当前事实：** 数据库存量大，但尚未确认自带 metadata 是否足以恢复同一组织块的切片顺序、切片间距、厚度、相邻关系和配准依据。HEST 中 29/492 条记录包含统一为 3.0 的 `z_step_size`，主要来自 Xenium；该字段描述采集或分割的 z-step，不能污染性地解释为组织切片间距或顺序。
- **当前影响：** 不阻塞主关键路径 R-01；H-06/H-07 保持未判定，不得因字段缺失提前判为失败或默认可测。
- **下一判定点：** R-03 形成组织块级 serial-section 清单。
- **硬阻塞条件：** 若现有数据均无法建立可审计的相邻切片关系，R-11 的邻近平面检验不可执行；届时提出针对性外部数据申请或取消三维相关叙事。

## I-003｜第二个 claim-bearing 结构尚未确定

- **状态：** `RESOLVED`
- **发现日期：** 2026-07-31
- **影响节点：** R-02、R-08、R-09，以及 H-04/H-05 的方法级通用性。
- **当前事实：** 按既定判据逐项比较后（D-035），第二 claim-bearing 结构于 2026-07-31 冻结为 TUMOR_STROMA_BOUNDARY：Heiser 标注的 carcinoma_border/carcinoma_edge 在 training 谱系给出 191 个可审计连通域组件、19 位患者、30 个 piece，输入隔离规则可执行；血管与坏死在所有已查来源中可审计实例仍为 0，不能入选；adenoma_border 仅癌前上下文，登记为非确认 context。此前记录的 GSE274103/GSE274557 acellular-interface 选择偏差事实不变，validation 侧的边界主张仍需原始空间审计。
- **解除方式：** 范围更新获批后获得 Heiser 逐 spot 病理标注，经 piece 级 crosswalk 与坐标重放建立可审计 GT；`structure_ontology.tsv` 中 TLS 与 TUMOR_STROMA_BOUNDARY 均为 FROZEN_CLAIM_BEARING，机器检查 `second_structure_frozen=true`。
- **剩余限制：** 第二结构的 GT 只覆盖 training 谱系；validation 谱系的 TSB 确认性检验并入 I-001 跟踪，未解除前 R-08/R-10 不得声称跨队列完成。
- **再次硬阻塞条件：** 若 Heiser 标注语义被证明不等于肿瘤—基质边界、commit 内容漂移或 crosswalk 失效，第二结构回退为未冻结并按 D-035 复查条件重审。

## I-004｜聚合来源之间的物理样本重复范围未知

- **状态：** `OPEN`
- **发现日期：** 2026-07-31
- **影响节点：** R-01、R-02 及所有患者或组织块级外层验证。
- **当前事实：** HEST、论文 atlas、GEO、10x、HTAN/OSF、STOmicsDB 和 SpatialDB 可能收录同一患者、组织块、切片或其重处理版本。相同校验和只能确认资产级重复；不同校验和不能证明物理样本独立。
- **已确认或冲突证据：** `repo/data_meta/ST_CRC_cohort_meta2.csv` 的 48 行中有 12 行自带 `duplicated=Yes`；本地 HTAN 49 个 h5ad 与该表 48 个资产名之间存在两个本地独有和一个 metadata 独有文件。HEST 有 6 条 `patient` 与 `subseries` 的显式 P 编号冲突。HEST 与本地 10x manifest 之间有 90 条记录共享 35 个数据源页面：47 条仅确认同一来源数据集谱系，15 条保留为可能物理相关并保守同组，28 条经显式映射审计后没有逐条关系证据且不建立 leakage edge。四个新增 GEO accession 均与既有 Atlas 行及本地镜像合并为同一 canonical source lineage，未重复计数；GSE274103 与 GSE274557 共享研究团队，但因 BioProject、论文与队列均不同，当前不合并患者，保留 provenance 风险。
- **当前影响：** 不阻塞 metadata 来源盘点；在物理身份澄清前，不得因来源数据库、目录或文件格式不同而把样本分配到相互独立的 claim-bearing 角色。
- **下一判定点：** 针对性外部 metadata 获批后，用 patient↔block crosswalk 复核 `r01_metadata_request.tsv` 的 provenance group；同一 PMID、accession 或 DOI 目前只代表来源重叠，不自动证明或否定物理独立。
- **硬阻塞条件：** 若候选核心队列无法形成彼此物理独立且身份可审计的训练与验证组，则 R-01 升级为 `BLOCKED_INDEPENDENCE`，后续确认性 benchmark 不得启动。

## I-005｜存储安全余量有限

- **状态：** `OPEN`
- **发现日期：** 2026-07-31
- **影响节点：** 所有需要解包、复制或生成派生数据的节点。
- **当前事实：** 2026-07-31 在 R-02 结束前复查，项目所在文件系统有 3,098,826,780,672 bytes（约 3.10 TB）可用，满足至少保留 2 TB 的要求，安全余量约 1.10 TB。2026-08-03 三路径数据落盘（GSE175540 解包 1.1 GB、两个 Zenodo 数据集 5.6 GB）后复查：2,982,261,653,504 bytes（约 2.98 TB）可用，仍满足要求，安全余量约 0.98 TB。
- **当前影响：** 不阻塞小型 metadata 扫描和注册表输出；禁止无估算的大包解压、全量复制或大规模派生物落盘。
- **下一判定点：** 每个预计产生显著落盘的新任务开始前检查可用空间，并记录预计峰值。
- **2026-08-13 D-057 更新：** 用户批准将当前运行时硬下限从历史的 2 TB 调整为 1.2 TB，1.4 TB 作为软提醒线；历史复查数字不改写。
- **硬阻塞条件：** 任务预计峰值会使可用空间低于 1.2 TB，或运行中实际余量接近该阈值时，升级为 `BLOCKED_STORAGE` 并停止新增落盘。

## I-006｜patient/block 可追溯核心单元不足

- **状态：** `RESOLVED`
- **发现日期：** 2026-07-31
- **影响节点：** R-01 的完成判定和角色冻结，并影响 R-02、R-05、R-10 的独立性。
- **当前事实：** 用户已批准 metadata-only 获取。956 条 physical-unit 记录中，121 条来自 6 个合格逻辑单元：HTAN Vanderbilt CRC 与 10x Breast Block A 有真实 patient+block；GSE211956、GSE226997、GSE274103、GSE274557 以 E2 patient-linked physical specimen 条件计入，`block_id` 全部保持空。四个 GEO 研究各只计 1 个单元，PDX、scRNA、Atlas/GEO 镜像和同患者多 specimen 未增加逻辑计数。
- **当前影响：** 6 个 outcome-blind 角色已冻结，独立 external lineage 为 GSE211956，机器 gate 为 `COMPLETE_WITH_EXCLUSIONS`。R-01 身份阻塞解除，但这不证明任何结构 GT，也不提高 plan.md 假设可信度。
- **解除方式：** 使用官方 GEO sample/series metadata、稳定 BioSample locator 和一个官方 patient—GSM supplement crosswalk建立 E2 互证；没有下载或保留表达矩阵、图像、结果表，也没有联系作者或访问受控/付费数据。
- **下一判定点：** R-02 审计结构 GT；若真实 block crosswalk 后续出现，替换 specimen equivalence 并重跑 leakage audit。
- **再次硬阻塞条件：** 任一 equivalence 的 patient、实体组织或 locator 证据失效，发现跨 study specimen 复用，或 external lineage 与其他 claim-bearing 角色共享 leakage group 时，R-01 回退为 `BLOCKED_IDENTITY`/`BLOCKED_INDEPENDENCE`。

## I-007｜三个 metadata extractor 的可选输出路径未统一加固

- **状态：** `OPEN`
- **发现日期：** 2026-07-31
- **影响节点：** R-01 的非默认复跑和未来复用这些 extractor 的工程安全性。
- **当前事实：** Atlas、HEST、HTAN extractor 的默认输出均位于可重建的 `infra/sample-registry/staging/`，但自定义 CLI 输出路径尚未统一采用项目根边界和原子替换；中断时可能留下部分多文件输出。
- **当前影响：** 不影响当前注册表、机器 gate 或科学结论；默认路径已通过两次确定性构建。未加固前不得把自定义输出指向不可重建资产或项目外路径。
- **下一判定点：** R-01 解阻后、首次修改 extractor 或使用非默认输出路径前，统一复用已有的根目录校验与原子写入实现。
- **硬阻塞条件：** 后续任务必须覆盖不可重建文件、写到项目边界之外，或多文件原子性成为交付前提时，先修复再运行。

## I-008｜Atlas 与官方 GEO metadata 存在来源污染和错误归因

- **状态：** `OPEN`
- **发现日期：** 2026-07-31
- **影响节点：** R-01、R-02 及任何依赖 Atlas 患者数、论文或 accession 的队列选择。
- **当前事实：** GSE242311 的官方 GEO metadata 对应 5 位患者、6 份乳腺肿瘤 specimen 和 16 个 GSM，并引用 PMID 39456890；bundled Atlas 将 16 行近似当作 16 位患者，且关联 PMID 36674951。GSE246011 的官方 GEO 当前包含 7 个 GSM（4 STAD、2 BLCA、1 LUAD），而 Atlas 记录为 20 条 STAD。两组均未进入 R-01 核心集。
- **当前影响：** 这是明确的 metadata provenance 污染风险：若直接信任 Atlas，会夸大独立患者数、错配论文/癌种并制造虚假验证独立性。当前 gate 未受污染，因为这两组保持排除。
- **下一判定点：** 任何后续节点尝试使用 GSE242311 或 GSE246011 前，必须以官方 accession metadata 重建行级 crosswalk并解释 Atlas 差异。
- **硬阻塞条件：** 若计划中的核心/验证 lineage 只能依赖无法解释的 Atlas—官方来源冲突，则相关 lineage 先排除，核心集不足时 R-01 回退为 `BLOCKED_IDENTITY`。

## I-009｜外部补充包曾短暂越过 metadata-only 范围

- **状态：** `RESOLVED`
- **发现日期：** 2026-07-31
- **影响节点：** R-01 的数据边界、外部数据索引与污染审计。
- **当前事实：** Europe PMC 的 PMC11508537 `supplementaryFiles` 端点曾返回一个含 7 JPG、7 GIF 和嵌套 ZIP 的混合包；未解读图像，发现后立即删除。随后仅解包无图像 supplement，发现两份 XLSX 都是表达结果表而非身份 crosswalk，也立即删除。另有 4 份核验来源时下载的 PMC 全文 XML 超出交付边界，已删除。保留的 PMC10991508 官方 crosswalk 工作簿含未使用的测序统计列；提取器只读取 `Sample matrix` 中 patient 与 Visium GSM 两列，未将统计列写入证据或角色选择。混合包和两份结果表的 SHA-256 已记录；四份全文 XML 删除前没有计算 SHA-256，manifest 诚实记录路径、删除状态和 `not_computed_before_deletion`，不补造校验和。
- **当前影响：** 保留索引中表达文件、图像文件和结果文件均为 0；这些内容未进入注册表、角色选择或科学判断。污染已清除，但事件作为范围防护回归证据保留。
- **下一判定点：** 下次调用返回混合 supplement 的端点前先检查清单和 MIME/成员列表，仅提取 allowlisted metadata/crosswalk。
- **再次硬阻塞条件：** 若无法在不读取表达或图像内容的情况下提取所需身份信息，则停止并重新申请范围，不得静默扩大。

## I-010｜会话上下文含与本项目无关的 WSDM 标签

- **状态：** `RESOLVED`
- **发现日期：** 2026-07-31
- **影响节点：** 全项目的任务边界和上下文污染审计。
- **当前事实：** 本次会话注入的 agent context 中出现 `$CMEM wsdm2027` 标签，同时明确写着没有历史 session；该标签与本项目无关。项目工作区排除大型 data 和独立上游 repo 后检索不到 WSDM 字符串。
- **当前影响：** 这是会话级上下文污染，不是样本或训练/验证泄漏。该标签未用于来源选择、角色冻结、代码、文档内容或科学判断，本项目也不与 WSDM 建立任何关系。
- **下一判定点：** 后续会话初始化时继续忽略不属于 `/home/huyudi/013_spatial` 的竞赛/项目标签，并以本仓库 `STATUS.md`、`docs/plan.md` 和 `docs/roadmap.md` 为控制面。
- **再次硬阻塞条件：** 若外部上下文开始改变数据选择、指标、角色或科学问题，立即停止并清除污染后重审相关产物。

## I-011｜R-04 真实数据 input manifest 与冻结 gene universe
- **状态：** `RESOLVED`
- **发现日期：** 2026-08-07
- **解决日期：** 2026-08-10
- **影响节点：** R-04 的真实 discovery、HTAN training 和独立 validation。
- **解决事实：** 已从冻结 `outer_splits.tsv`、HTAN replay index 和显式 override table 生成 71 个 eligible molecule-only locator：HTAN training 47、ST_CRC_CMS internal validation 14、TLS_VISIUM_USZ external validation 8、TENX breast serial validation 2。ST_CRC 只抽取 filtered H5，乳腺只抽取坐标；没有读取图像、原始矩阵或 GT。
- **预检证据：** `infra/r04/preflight_gate.json` 为 `READY`；71 个 section 均通过整数 counts、有限坐标、barcode/坐标交集、重复 ID 和患者分组检查。HTAN 30 患者训练数据冻结 4,000 个 gene universe；ST_CRC/乳腺覆盖率为 99.375%，USZ 为 90.825%。一个 HTAN section 的 1 个零 library spot 被确定性剔除并记录。
- **当前影响：** 真实 CPU pilot 和后续冻结 discovery 已解除输入阻塞；该结果只证明输入可审计，不提高或降低 H-01/H-03/H-05 的科学信心。当前双模型 pilot 仍为 `FIT_COMPLETE_NOT_VALIDATED`。
- **下一判定点：** 对当前冻结接口做五次 restart、患者/组织块 bootstrap、lengthscale sensitivity，并在未参与拟合的 validation rows 上只推断场值；同时补齐共享 loading、dispersion/profile 和 model-selection uncertainty。
- **再次硬阻塞条件：** 若任一后续 lineage 的 gene coverage 低于 80%、输入 hash 与 manifest 不一致、或出现非整数/目标派生输入，立即回退为 `BLOCKED_INPUT_CONTRACT`，不得用聚类或图平滑替代。

## I-012｜组成参考的 posterior uncertainty 尚未接入正式估计器
- **状态：** `OPEN`
- **发现日期：** 2026-08-07
- **影响节点：** R-04 的 composition/state 拆分与 H-03。
- **当前事实：** 005/006 参考 adapter 已实现 donor/type cap、gene intersection 和 h5ad allowlist；组成调整函数支持 posterior draws，但 cell2location/RCTD 的正式运行和多重插补尚未接入。当前 deterministic reference 只能产生 `POINT_REFERENCE_NO_POSTERIOR`。
- **当前影响：** raw latent field discovery 可以继续；组成之外的强结论和 cell-intrinsic state 结论必须保持 `PARTIAL_RAW_FIELD_ONLY` 或 `STATE_NOT_IDENTIFIABLE`。
- **下一判定点：** 在冻结候选后运行 cell2location，保存 provenance、cell-type vocabulary、posterior draws 和 gene-fold cross-fit manifest，再运行 RCTD sensitivity。
- **硬阻塞条件：** 若参考无法排除物理重叠或不能提供可审计 cell-type uncertainty，停止组成/状态 gate，但不回滚无标签 raw field discovery。

## I-013｜双模型当前的后验不确定性主要覆盖诱导权重
- **状态：** `OPEN`
- **发现日期：** 2026-08-07
- **影响节点：** R-04 候选 uncertainty、患者 bootstrap 和最终 gate。
- **当前事实：** TensorFlow 实现已对 section-specific inducing weights 做 variational draws，并输出 field mean/SD；共享 loading、dispersion、长度尺度和模型选择不确定性尚未由同一后验完整覆盖。合成 smoke 不等于完整科学不确定性。
- **当前影响：** 当前输出可用于工程 smoke 和方法调试，不能单独签署 `PASS_R04_STABLE_FIELD`。
- **下一判定点：** 加入 loading/dispersion 的变分或外层 bootstrap、5 restart、patient/block bootstrap 和 lengthscale sensitivity，并把 uncertainty coverage 写入 gate。
- **硬阻塞条件：** 若最终 uncertainty 只能依赖单次点估计，R-04 只能保留连续 field 的探索性输出，不能进入强确认性结论。

## I-014｜signed residual adapter 仍是 pilot likelihood，不是完整 NB-GP 后验
- **状态：** `OPEN`
- **发现日期：** 2026-08-07
- **影响节点：** R-04 双模型一致性、signed-only 候选解释和最终 gate。
- **当前事实：** signed 路线仍先用固定过度离散参数构造 NB Pearson nuisance residual，再对残差做 MSE 拟合；diagonal variational inducing-weight GP-KL 已纳入训练，fit 与 frozen infer 已统一为按 spot 平均、按 gene 求和的目标尺度，但还没有把完整 NB observation likelihood、dispersion uncertainty 和 loading uncertainty 统一纳入后验。
- **当前影响：** signed-only 或 `BOTH` 候选可用于方法调试和探索性发现，不能把该适配器的输出写成已完成的 signed NB-GP 科学证据；gate 继续保持关闭。
- **下一判定点：** 在真实 manifest 冻结后比较完整 signed NB likelihood、当前 residual pilot 和非空间残差基线，并传播 dispersion、loading、lengthscale 与 inducing weights 的不确定性。
- **硬阻塞条件：** 若 signed 路线无法在 grouped held-out likelihood 上超过 nuisance/非空间基线，删除 signed 模型贡献，不得用普通聚类或调参后的 marker program 替代。
## I-015｜R-04 正式 restart 的 CPU 时间可行性
- **状态：** `RESOLVED`
- **发现日期：** 2026-08-10
- **影响节点：** R-04 五次 restart、后续 bootstrap 和 validation。
- **当前事实：** 47 个 training section 的 4,000-gene panel 已分段物化并通过 47/47 loader 验证。修正 full-library offset 后，缓存单 mNSF、单 step 实测 wall time 为 69.54 秒，峰值 RSS 为 14,103,928 KiB；双模型 3-step 原始读取校准在 120 秒外层上限内未落盘，状态为 `NOT_RUN/BLOCKED_TIMEOUT`。
- **当前影响：** 正式 300-step×5 restart 暂不启动；这不是科学阴性结果。GPU 当前仍被项目政策禁止，需用户明确批准后才可提交租用申请。
- **下一判定点：** 先完成等价性与多 seed 诊断，给出 CPU cache/线程优化后的完整时间估计；若仍超出可接受窗口，提交 GPU 需求、预计时长和成本供用户审核。
- **解除条件：** 完整训练在可接受资源窗口内可重复、可恢复，并通过 factor collapse、ELBO 和 checkpoint preflight。

**2026-08-13 更新：** D-049 修正了 gene-minibatch likelihood 目标，因此 2026-08-10 的 69.54 秒单步历史测量不能直接作为当前正式训练时间。修正目标后，对同一 47-section、4,000-gene cache 的单 mNSF step 为 93.63 秒、峰值 RSS 14,713,848 KiB；该 run 仍为 `FIT_COMPLETE_NOT_VALIDATED`，没有解除本条硬阻塞。当前只保存资源重测，不据此授权正式 restart。

**2026-08-13 D-050 更新：** 科学校准已通过，但 inducing 参数化与初始化再次改变，旧 93.63 秒结果不能作为最终版本资源依据。对同一 cache、K=2、16 inducing points、lengthscale=3.0 的新版单步实测为 95.55 秒、峰值 RSS 16,277,456 KiB。按线性下限，600-step 单 restart 约 15.9 小时、五次 mNSF restart 约 79.6 小时，且不含 signed 模型和验证。CPU 时间硬阻塞仍成立；下一解除条件为用户明确批准 GPU，或明确接受多日 CPU 运行并确认调度方式。

**2026-08-13 D-056 解决：** 旧估计错误地把包含冷缓存读取、构建和写盘的 1-step 总墙钟逐步相乘。当前代码新增 optimizer 内部计时；同一 47-section、4,000-gene、K=2 配置的 20-step run 总墙钟 272.31 秒、峰值 RSS 16,264,992 KiB，其中 optimizer 92.80 秒（4.64 秒/步），固定成本约 179.51 秒。按固定成本加实测斜率，600-step 单 mNSF restart 规划值约 0.82 小时、五次约 4.12 小时；不含 signed、真实 K 搜索和 validation，且不是运行时保证。CPU 时间已从多日硬阻塞降为可调度的小时级工作，本 issue 关闭，暂不申请 GPU。若真实首个 600-step 的每步斜率超过当前两倍或影响系统正常工作，按 D-056 重新打开并提交资源选择。
## I-016｜factor collapse 与不规则场恢复尚未通过科学校准
- **状态：** `RESOLVED`
- **发现日期：** 2026-08-10
- **影响节点：** R-04 候选场发现、跨 restart 稳定性和 H-01/H-03/H-05 信心。
- **当前事实：** 47-section 单步 mNSF calibration 的 loading cosine 为 0.9974、effective rank participation 约 1.00；这只能触发告警，不能代替平台期、多 seed 和 planted-factor 判定。三个独立 seed 的 20-step planted-factor 诊断仍有 loading cosine 0.9964–0.9978、effective rank participation 1.0022–1.0036，matched field absolute correlation 约为 [0.248, 0.396]、[0.487, 0.214]、[0.567, 0.103]；这说明重复方向不是单一 seed 的偶然现象，但尚未达到平台期判定。crescent、branch、disconnected 三类支持的 10-step synthetic calibration 也均未恢复 planted fields：mNSF 平均 matched absolute correlation 分别约 0.107、0.147、0.096，signed 分别约 0.054、0.080、0.053；所有 artifact 均明确为 `CALIBRATION_COMPLETE_NOT_VALIDATED`。
- **当前影响：** 不加入强制正交惩罚，不命名重复方向，不把低相关结果写成“模型不支持不规则场”；正式 restart 保持科学门禁。
- **下一判定点：** 完成 tiny dense oracle、至少 3 个诊断 seed、planted-factor recovery、训练平台期和 crescent/branch/disconnected sensitivity；若有效秩稳定低于设定 factors，比较降低 K 后的 held-out 预测损失。
- **解除条件：** planted control 可恢复且平台期无未解释 duplicate/inactive factor；否则保留较小 K 或收缩候选 claim。

**2026-08-13 更新：** objective oracle 已通过，dense 与所有 3-gene batch 组合误差约 `1.9e-6`；但修正后的正确指定双场控制在 crescent、branch、disconnected 上 300 steps 均未达到平台，1000-step disconnected 仍未平台。三种几何的第二 canonical correlation 分别约为 0.019、0.100、0.197，而各自 coordinate-permutation null 的 95% 上界约为 0.197、0.193、0.243；K=0、K=1 和共线控制已生成但尚未到平台。一次 2-fold、单 seed、20-step 的 K=1/2/4 held-out 比较给出 K=4 的探索性 one-SE 候选，但最终状态仍为 `K_NOT_IDENTIFIABLE`。当前状态仍为 `OPEN`，不能将其解释为真实数据只有一个 factor，也不能加入正交惩罚掩盖问题。

**2026-08-13 D-050/D-051 解决：** 根因包括 inducing-value 参数化、核尺度、K 依赖初始化、合成 gene effect 共线及 held-out 评价泄漏；均已修正。三个数据/优化 seed、三种不规则几何的双独立场达到平台并恢复两个 effect directions；单场和共线双源只支持一个方向。无泄漏 K=1/2/3/4 聚合选择 K=2，无场 K=0/1/2 聚合选择 K=0；18-control 科学门禁为 `PASS_SCIENTIFIC_CALIBRATION`。本 issue 关闭。该结果只证明方法校准通过，不代表真实数据 K=2，也不提高 H-01/H-03/H-05 的科学信心。

**2026-08-13 D-052/D-054/D-055 最终审计：** 上述首版 K 聚合因 split seed 间接包含 K 而被撤回；随后又修复浮点重放、best loss/参数错位和 NaN 绕过。所有进入门禁的 18 个控制、24 个双场 K 单元和 12 个无场单元均用最终数值语义重跑。K=2 与 K=0 仍分别在三个独立 data seed 上一致最佳，K=4 支持比例为 0；v10 总门禁核验 hash、重放 generator、重新聚合 shards 并检查有限性后通过，I-016 继续保持 `RESOLVED`。若同题、有限性或重放约束失效则立即重新打开。本结果仍只恢复方法可测性，不改变患者级假设信心。
## I-017｜缺失 gene 的条件推断接口已实现，但真实 USZ 可测性未完成
- **状态：** `OPEN`
- **发现日期：** 2026-08-10
- **影响节点：** R-04 external validation 与最终跨队列结论。
- **当前事实：** infer 现在可在所有 validation section 的 common observed panel 上冻结模型参数并推断，禁止 zero-fill，并记录 panel hash、coverage、retained loading energy、rank loss 和 offset source；但尚未在真实 USZ 上完成 full-versus-mask technical control，也未证明每个 factor 的 retained energy 与最差 patient 阈值。
- **当前影响：** USZ 结果只能保持 `NOT_TESTABLE`/`INFERENCE_COMPLETE_NOT_VALIDATED`，不能被写成独立验证阴性或阳性。
- **下一判定点：** 用完整 gene 技术样本施加 USZ 同一 mask，比较 full/masked field correlation、factor loading energy 和 rank；通过后再运行真实 external validation。
- **解除条件：** common observed panel 不造成 factor-specific 信息丢失，且 full/masked 判定一致；否则仅保留覆盖充分的 validation lineage。

## I-018｜D-049 之前的 minibatch 校准结果与当前目标不再可直接比较
- **状态：** `RESOLVED`
- **发现日期：** 2026-08-13
- **影响节点：** R-04 历史 calibration、资源外推和 factor gate。
- **当前事实：** 旧实现把 spot-gene 联合均值直接乘以 `G/B`；D-049 已改为先按 spot 求和 gene loss 再缩放，并由 objective oracle 验证。旧的 20/80/300-step 结果仍保留为历史诊断，但不能与修正目标下的结果合并或用于正式资源外推。
- **当前影响：** 需要以修正目标重新生成进入 gate 的 calibration 和资源基准；旧结果不支持真实数据结论，也不自动构成回归失败。
- **下一判定点：** 完成修正目标下的 dense/full-batch/minibatch 梯度与 checkpoint 一致性复核，并确认所有正式输入的 config hash 已更新。
- **解除条件：** 新目标下的 oracle、平台期控制和资源基准均有独立 provenance，历史 artifact 被明确标记为不可合并。

**2026-08-13 D-050 更新：** objective v3 同时通过 objective 和平均 minibatch gradient oracle；fit/infer 已统一 dense 缩放，D-050 后的合成校准与 47-section 资源重测均有独立 artifact。D-049 之前及 D-049 到 D-050 之间的 artifacts 继续保留但不合并。本 issue 关闭。

## I-019｜R-04 全局表示秩、fold 异质性与下游 K 稳健性尚未闭合
- **状态：** `OPEN`
- **发现日期：** 2026-08-13
- **影响节点：** R-04 正式候选发现及 H-01/H-03/H-05。
- **当前事实：** 合成校准只证明方法在受控数据上能识别效应维度，不能给真实数据指定 K。真实数据已有五 seed 的 K=0/K=3 极端 fold 诊断和 restart 2–4 的冻结留出推断；fold 方向可重复，但幅度存在异质性，前两类子空间方向较稳定，第三方向较弱且不稳定。严格训练审计仍为 `STOPPED_OFF_PLATFORM`，冻结 practical acceptance 仅覆盖 held-out inference/scoring；当前冻结 cell 没有 GT 对齐的 held-out spot×gene 空间效应矩阵。
- **当前影响：** `selected_k=null` 保留，且不单独构成永久硬阻塞；现有证据可以支持 K 作为模型表示容量的诊断，但还不能计算受完整门禁支持的 `K_eff`，也不能判断主要结构读出是否依赖第三方向。R-04 科学节点仍未完成。
- **下一判定点：** 已完成最新面板的 provenance-checked 聚合；下一步在 training-role GT 的嵌套 cross-fitting 中使用旋转稳健的空间效应表示，比较共享成分与结构特异残余以及完整 K=3 与稳定子空间版本。只有下游结论对第三方向敏感或无法判定时，才触发 K=2 bridge；随后再进入空间 null、组成拆分和独立 lineage 复现。
- **硬阻塞条件：** 只有当合理 K、seed 或 fold 的变化导致主要结构结论方向不稳定，或空间置换 null 后稳定表示不优于 K=0/非空间基线时，才将 R-04 的科学路径判为硬阻塞；`selected_k=null` 本身不是硬阻塞。
- **历史记录边界：** 下方带日期的条目保留当时的诊断状态和停止决定；不将历史的“尚未运行”描述改写成当前事实。
- **2026-08-13 执行更新：** 已实现真实 training manifest/cache 的患者分组 K 搜索入口和同 nuisance 的 K=0；一折 K=0/1、2-step wiring smoke 产出 47-section/30-patient 路径结果，但 fit/inference 均未平台，已明确标为 `WIRING_ONLY_NOT_SCIENTIFIC`。正式 5-fold K=0 向上搜索、空间置换 null、五次 restart 和组成拆分仍未运行。
- **2026-08-13 当前状态：** 完整 5-fold、K=0–4 粗筛已启动为持久 CPU 任务（PID 2623526），但在 `r04_real_k_search.json` 完成前仍保持 `NOT_RUN`；运行中断或任一 K 未平台时不得从部分 cell 推断真实 K。
- **2026-08-14 当前状态：** 完整 5-fold、K=0–4 粗筛已完成 25 个 cell，但 25/25 个 fit platform 均未通过 600-step 平台期门禁；机器结果为 `K_NOT_IDENTIFIABLE_NOT_CONVERGED`、`selected_k=null`。K=3 的较高原始预测分数仅作诊断信号，不能进入科学候选。当前正式 R-04 gate 因训练收敛处于硬阻塞；一个独立的持久 CPU 哨兵任务正在以 K=0/2/3、fold 0、1200 steps 诊断步数不足与优化波动，默认学习率仍为 0.05，尚未改变正式协议。空间置换 null、五次 restart、组成拆分和独立 lineage 复现继续保持未运行。
- **2026-08-16 诊断更新：** 0.05 学习率的 1200-step 哨兵已正常结束，但 K=0/2/3 的 3 个 fit 仍全部未平台；最佳步为 1182、1184、1200，支持“训练仍未走完/优化参数需要诊断”的判断，但该单 fold 任务仍为 wiring-only。已启动下一轮独立低学习率哨兵：完整输入、K=0/1/3、fold 0/4、1200 steps、learning rate 0.01；结果完成前不改变 `K_NOT_IDENTIFIABLE_NOT_CONVERGED`，不进入空间 null 或 restart。
- **2026-08-18 诊断更新：** 修复了 fold 子集最终 JSON 的 `wiring_only` 标记并通过回归测试；真实训练入口新增 gene batch、学习率衰减和梯度/参数轨迹记录。当前已持久启动 K=0/3、fold 0 的 batch=256、512 和短程 full-batch 诊断，均为 wiring-only，尚未改变 `K_NOT_IDENTIFIABLE_NOT_CONVERGED`。
- **2026-08-18 诊断完成更新：** 三个 batch 任务均以 systemd `success` 正常结束，共生成 6 个 `FIT_AND_SCORED` cell；但 K=0/3 的 fit 与 inference 收敛门禁全部为 false，因此没有科学 K 结果。batch=512 的末段 loss 波动较 batch=256 低约 35%–44%，full-batch 只作为短程噪声参照。下一轮固定 batch=512，单独比较固定学习率与线性衰减的 2400-step wiring-only 运行。
- **2026-08-18 2400-step 更新：** 固定 lr=0.01 与早期线性衰减两个任务均以 `success` 完成，生成 4 个 cell；K=0/3 的 fit 和 inference 仍全部未平台。固定 lr 的末段损失优于衰减方案，K=3 有效秩参与度约 1.83 且未出现 collapse 告警，但末段仍未达到 0.1% 平台门槛。下一步保持同一协议延长到 4800 steps，仍不进入正式 K 聚合。
- **2026-08-19 continuation 更新：** 4800-step 延长任务也正常完成，但 K=0/3 fit 与 inference 仍未平台；K=0 和 K=3 的共同失败优先指向共享优化时程。已实现并通过回归测试的显式 checkpoint continuation 将从同一 K=0、fold 0、4800-step checkpoint 分叉 control/lr-step-down 两条 fit-only 诊断，保留 Adam state；不运行 K=3、inference 或正式聚合。
- **2026-08-20 staged optimizer:** implemented and regression-tested, but no real-data evidence yet. The first representative K=0/K=3 run remains required; until it passes the existing fit and inference gates, R04 stays `K_NOT_IDENTIFIABLE_NOT_CONVERGED`.
- **2026-08-20 staged run started:** K=0/3, fold 0, 2400 steps, batch 512, learning rate 0.01, and 600 shared-first steps is `RUNNING`; it remains wiring-only and cannot change the R04 status before fit and inference gates are evaluated.
- **2026-08-20 staged run completed:** both cells finished normally, but fit and inference platform gates were false for K=0 and K=3. K=0 matched joint training; K=3 had lower best fit loss but worse held-out score and no convergence recovery. Staged training is rejected as the common protocol; R04 remains `K_NOT_IDENTIFIABLE_NOT_CONVERGED`.
- **2026-08-20 checkpoint audit started:** deterministic full-panel objective audit is `RUNNING` for existing K=0/K=3 fold-0 checkpoints. It is read-only diagnostic evidence; until repeatability and checkpoint-to-checkpoint drift are assessed, no optimizer change or gate relaxation is authorized.
- **2026-08-20 checkpoint audit completed:** K=0 dense objective was effectively flat from step 6600 to 7200 with zero repeated-evaluation spread, despite the raw minibatch gate being false. K=3 still declined materially from step 3600 to 4800. The old gate is therefore unreliable for K=0, while K=3 remains genuinely unresolved; deterministic monitoring must be added before further optimizer comparisons.
- **2026-08-20 K=3 continuation:** the first service attempt failed before model import because it used a system Python without SciPy; no scientific computation or checkpoint was produced. A second service was relaunched with the verified Anaconda environment and is currently `RUNNING`; its output remains diagnostic-only until completion and replay checks.
- **2026-08-21 K=3 continuation completed:** the inherited-state fit reached step 7200 and produced a final checkpoint. The dense objective was nearly flat in the last 600 steps, but the original fit gate remained false and inference was not run. A top-level summary serialization bug incorrectly wrote `k_model=0`; it was fixed, the existing artifact was corrected, and 11 targeted regression tests passed. The identical-checkpoint repeat audit is now `RUNNING`; no scientific gate has been changed.
- **2026-08-21 deterministic platform audit completed:** under the same dense objective and fixed draws, step 6600→7200 decreased only 0.036%, and repeating step 7200 gave exactly the same objective. The unresolved issue is now a scientific gate-definition conflict, not a runtime failure: the deterministic full objective is stable while the stochastic minibatch gate is false. Escalate to a hard decision before inference, K selection, or the 25-cell campaign.
- **2026-08-21 dual-gate decision resolved:** user approved recording the deterministic full objective as the primary optimization-platform measure while retaining the minibatch gate as a warning. K=0 and K=3 fold-0 held-out inference (1200 steps, inherited 7200-step checkpoints) is now running; results remain diagnostic-only and cannot select K.
- **2026-08-21 first dual-gate inference completed:** K=0 inference platformed at 1200 steps; K=3 remained descending by about 0.36% and had a temporarily worse fold-0 held-out mean than K=0. This is not yet evidence against K=3 because inference length was insufficient. K=3 is being extended to 2400 steps; K=0 receives a 7200→7800 objective continuation before its matched 2400-step inference.
- **2026-08-31 five-seed restart audit:** all 20 cells are complete and all held-out inference runs platformed. Fold 0 is K0-favouring in 5/5 seeds and fold 4 is K3-favouring in 5/5, but fold-4 effect size varies from `+2.06` to `+59.36`; one patient in each fold changes sign across seeds. Cross-seed subspace diagnostics show a recurrent first direction, a weaker second direction and an unstable third direction. Replaying canonical external audits under the current two-sided rule gives 8/20 dense-objective passes; restart 2-4 account for the 12 failures. Because optimization seed jointly changes initialization, minibatch order and variational draws, this is overall optimization-randomness stability rather than pure initialization stability. The hard R-04 state is `K_UNDECIDED_OPTIMIZATION_NOT_RESOLVED`; no K, spatial null, composition split or field naming is authorized.
- **2026-08-31 next decision point:** run a bounded, seed-preserving, fit-only continuation for restart 2-4 from step 7200 to 8400 with dense evaluation every 400 steps. If all paired K0/K3 cells platform, rerun inference and reassess top-1/2/3 subspaces; if any remain off-platform, stop and decide between further optimization diagnosis and a matched dense/full-batch tail. A later K=2 bridge is permitted only after optimization error is resolved and the third direction remains unstable.
- **2026-08-31 step-8400 hard block:** all 12 bounded fit-only continuations and frozen-checkpoint repeat audits completed with valid source/checkpoint provenance, but only 10/12 passed the combined gate. Restart 3/K=3/fold 0 remains above the `0.1%` dense-objective threshold in one final interval; restart 2/K=3/fold 4 has a one-draw bitwise repeat mismatch of `0.000488`. D-081 therefore triggered its stop condition. Held-out inference, score recomputation, K expansion, spatial null, composition split and field naming are `NOT_RUN`; the next action requires a new decision and cannot be silently replaced by a relaxed tolerance or an extension of only the failing K=3 cell.
- **2026-08-31 strict-audit update:** after fail-closed schema, ridge, checkpoint-bundle and frozen-inference hardening, all 12 repeats were rerun rather than patched. The strict artifact set has no provenance or structural errors, records `ridge=1e-5`, and is exactly repeatable in 12/12 within-run pairs; 11/12 cells pass the full gate. Restart 3/K=3/fold 0 still fails because 7600→8000 changes by `0.1006526%`, above the fixed `0.1%` threshold. The first audit's intermittent `0.000488` repeat difference remains a numerical-risk observation even though it did not recur. The hard block and all downstream `NOT_RUN` states remain unchanged.
- **2026-08-31 practical-acceptance update:** the user accepted the single marginal restart 3/K=3/fold 0 exception for the limited purpose of frozen held-out inference. The strict artifact remains `STOPPED_OFF_PLATFORM`, 11/12; a separate hash-bound artifact authorizes no fit updates and no K selection. Restart 2-4 frozen inference/scoring is now `RUNNING` under three CPU workers, so the optimization hard block is lifted to `OPEN`; it returns to `HARD_BLOCKED` if any cell or gene split is incomplete/off-platform, provenance changes, or the completed five-seed score/subspace assessment remains non-identifiable. Spatial null, K expansion, composition splitting and field naming remain `NOT_RUN`.
- **2026-09-01 frozen inference completion:** the authorized restart 2–4 GPU panel completed 12/12 cells with no failed cell; remote results were recovered as 27 files and all SHA-256 checks passed. The outputs remain `FROZEN_HELDOUT_DIAGNOSTIC_NOT_K_SELECTION`, `selected_k=null`, and strict audit remains `STOPPED_OFF_PLATFORM`. The CPU/GPU equivalence pilot was intentionally stopped at user request to avoid rented-GPU idle cost and is `NOT_RUN`; the practical diagnostic can be used, but no cross-device numerical-equivalence claim is allowed. The next unresolved issue is aggregation of patient-level K3-K0 scores and cross-seed subspace stability; K expansion, spatial null, composition splitting and field naming remain `NOT_RUN`.
- **2026-09-01 semantic and evidence-integration update:** the five-seed panel was aggregated with the recovered frozen cells, explicit historical endpoint exceptions and separate endpoint/environment strata. The new semantics artifact records `K_model` as representation capacity, leaves `K_eff` uncomputed because the frozen cells lack GT-aligned held-out spatial-effect matrices, and records structure readout/downstream K robustness as `not_tested`; this is an evidence gap, not a negative result. No K=2 bridge is triggered in this round because no downstream structure readout has yet shown sensitivity to the third direction. The next issue-closing evidence is the training-role nested cross-fitted readout; spatial null, composition split and independent lineage replication remain pending.
- **2026-09-01 training-role readout update:** all 10 K=3 training-effect exports (five seeds × folds 0/4) were recovered from frozen checkpoints with the primary centered spatial-rate contribution, zero fit updates, complete array/hash checks and validation GT sealed. The shared-versus-specific readout ran by patient-level inner cross-fitting, but the paired TLS/boundary support is only two patients per extreme outer fold. Therefore `downstream_k_robust=not_tested`; the mixed descriptive rank-2/rank-3 increments are not evidence that the third direction is biologically meaningful or meaningless. Single-factor naming remains forbidden.
- **2026-09-01 conditional bridge update:** because the training-role evidence cannot judge third-direction sensitivity, the six-cell K=2 Stage-A manifest is frozen. Its first K=2 CPU fit was stopped before the first checkpoint; no K=2 score, readout or K conclusion was produced. The stage artifact is `BLOCKED_CPU_RUNTIME`, with historical matched CPU fit estimates of approximately 1.25–6.9 hours per cell. The current linked handoff is `k_semantics_and_downstream_robustness_20260901_final.json`. The next decision is a bounded GPU approval or explicit acceptance of the corresponding CPU runtime; spatial null and composition split remain downstream and are not silently substituted.
- **2026-09-01 Torch backend migration:** R-04 active compute has been refactored from TensorFlow/TFP to PyTorch (D-094); `environment.lock.json` `backend=torch` is now the single enforced backend, CPU/CUDA share one `torch.autograd`/`torch.optim` implementation, and no `import tensorflow`/`tensor_probability` remains in active code. Legacy TensorFlow checkpoints/artifacts are sealed as historical evidence and are not consumed by Torch training. Current `selected_k=null` is preserved; no new K=2/K=3 scientific result has been produced under Torch, and formal K comparison requires matched Torch reruns.
- **2026-09-02 Torch migration audit update:** the new checkpoint contract is now fail-closed against mixing `params` with `best_state`, malformed parameter layouts, missing parameter tensors, one-sided nuisance tensors and explicit library-size changes. The isolated Torch environment (including the existing h5ad reader dependency) passes the full 178-test R-04 suite and the objective oracle; this is engineering/provenance evidence, not a new real-data result. The issue remains `OPEN`: matched Torch K=0/K=3 real-data reruns are still required before the old fold-heterogeneity signal can be reassessed, and the K=2 bridge remains conditional.
