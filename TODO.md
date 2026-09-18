# TODO（项目共读入口）

当前唯一活跃工作面：R-16 局部恢复与推进（恢复包 20260918，分支 `recovery-20260918`，基准 `7399652`=D-127）。旧 Lane A/B 已完成（D-126/D-127），其否证外推标记为 `METHOD_LIMITED_PENDING_REANALYSIS` 待 T1/T2 复核。
科学判断见 `docs/plan.md`，节点状态见 `docs/roadmap.md`，技术选择见 `docs/decisions.md`。

## 下一步（恢复包 T0–T6；单 agent 按依赖串行）

- [x] T0 现状与输入契约：preflight 通过（HEAD=基准，registry 1950 行，cache 47 片，磁盘余量 1.62TB）；toy 审计 12 测试全过；`input_contract.json` 落盘（源 X=int64 原始 counts，r16 路径无深度校正=分析选择非数据缺失）
- [x] T1 PC 匹配符号修复＋联合解释力＋嵌套 Jaccard 降级（D-128：857→772组，NEW=0维持；I-023关闭）
- [x] T2-small Lane A 机械修复验证（D-129：同算子闭环+NEG干净+BUMP检出；remap降诊断对照；占位p记NOT_CALIBRATED）
- [ ] T2P2 变异函数匹配主 surrogate（硬骨头：非平稳+零膨胀+不规则几何无精确解；失败只封关联p，不阻塞T3/T4）
- [ ] T3 191 基因→≤10 多基因程序（执行中：candidate_union＋共表达模块＋program_definitions＋triage_log）
- [ ] T4 冻结程序 ST-CRC/USZ 独立复现＋增量预测（M0/M1/M2）
- [ ] T5 遮蔽任务接口（无适用目标则记 NOT_TESTABLE）
- [ ] T6 阶段总结：old_vs_new 对照＋程序卡＋final_decision＋canonical 文档更新
- [ ] hold 中（需逐项批准）：R-11 smoke、Xenium CRC、dbGaP、两封数据邮件、CRDC/NCI 凭据；GPU/新外部数据/付费资源本轮一律不动

## 当前执行分支（三线并行执行中；2026-09-10 用户批准开工）

- [x] A marker代理路线：三库集齐（CellMarker2.0 CRC子集2490行＋PanglaoDB GI489行＋CellTypist-CRC模型36类top50；HPA/TISCH因墙/无批量接口放弃，已记录）；投票5/6类可用（ILC判死）；smoke完成——a29信号被组成调整吃掉、0bd3不动（D-105，roadmap已记）
- [x] B scRNA参考路线：GSE236581参考构建完成（94k细胞×3854基因，marker抽查全对）；NNLS两版（全基因/DE限定）免疫份额均与经典marker零/反相关，按预定判据判死（LEADS L-003已放弃）；组成调整只走marker-proxy
- [x] C GPU租用成功 09-10：4090PLUS 48GB（¥2.49/h，6h到期），实例lyg1076；环境SETUP_DONE（torch2.5.1+cu124，确定性自检过）
- [x] C K=2完成 09-10：4090PLUS约1h跑完（K2−K0：fold-0 −16.61，fold-4 +71.96；K=3默认生效，D-106）；实例已退租、无残留磁盘；本次¥3.12，累计¥8.43
- [x] 汇总判读与记账完成 09-10（含本轮阶段总结推送）

## 当前执行分支（fold-0 外层验证；B 在 LEADS L-002 hold）

- [x] 输入确认：K=3 fold-0 checkpoint 存在；留出 11 section / 6 患者；GT 覆盖 TSB 8 section（5 患者），TLS 仅 1 section（不做统计，只报覆盖）
- [x] 外层验证脚本 + 4 测试全绿（`scripts/r04_explore_outer_validation.py`：训练冻结共享分→留出场重建→冻结映射→TSB AUC + section内置换 null；重建自洽门禁 fail-closed）
- [x] 患者扩容方案确认（无 GPU）：fold-4 全量导出并行启动，TSB 留出 8+7=15 section / 9 患者；TLS 两 fold 各仅 1 section，只报覆盖不做统计
- [x] fold-0 heldout 补跑收尾：heldout 版与拼装版一致（TSB 差 <0.005，同判 null），结论确认，无分歧
- [x] fold-4 导出触发分数门禁的处置：核验为等效最优点设备噪声（终点损失相对差≤4e-4），追加 D-103 接受重放场用于探索性外层验证（完成 09-05）
- [x] 外层验证脚本支持 direct 模式（sidecar 驱动加载 + hash/收敛门禁，不依赖 panel 状态；完成 09-05）
- [x] fold-4 外层验证完成 09-05（TSB rank3 弱阳性 0.544，rank1/2 null；TLS 单 section 只报数；`explore_outer_validation_fold4_20260905.json`）
- [x] fold-0 超时处置：双 split 已落盘验完（hash/收敛/fit_updates=0 全过），缺的只是 heldout_full；外层验证改走 split 拼装模式（脚本已支持+D-103 同等门禁）
- [x] fold-0 外层验证完成 09-05（TSB 三 rank 全 null 0.48–0.50；TLS 单 section 只报数；`explore_outer_validation_fold0_20260905.json`）
- [x] 双 fold 合并判读：外层零发现（fold-4 rank3 的 0.544 在 fold-0 无复现）；H 不变，R-04 未完成（完成 09-05）
- [x] fold-0 heldout 补跑收敛落盘 + heldout 版重跑完成 09-06：与拼装版一致（TSB 差 <0.005，同判 null），结论确认，无分歧

## 当前执行分支（TLS 共享生态稳健性深挖；B 在 LEADS L-002 hold）

- [x] 稳健性深挖 CPU 运行完成 09-05：fold-0 双向一致（rank1 基因余弦 1.0）+ rank3 一致性崩塌（第三方向患者特异噪声的直接证据）+ fold-4 双向反号 + 跨 fold 正交
- [x] 判读与记录：roadmap 日期条目 + 记忆；结果标 exploratory，基因只列 ID 不做生物学命名（完成 09-05）

## 当前执行分支（三方向并行探索，D-101 探索定位；子项可并行，无输入依赖）

- [x] A 空间 null：置换 null 完成 09-05（200 draws；唯一超带信号 fold-0 rank2 TSB p_up=0/200，单患者驱动，记探索性候选）
- [x] B 组成拆分：盘点确认无输入，不编造数字；转 scoped proposal（LEADS L-002：外部参考+审批，或 marker 代理+单独决策）
- [x] C 读出深化：rank 曲线 + 绝对 AUC + 范数 + 逐内折表完成 09-05（fold-0 首方向主导；fold-0 TLS shared AUC 0.86 vs fold-4 <0.5，患者异质又一证据）
- [x] 汇总：D-102 口径决策 + roadmap 日期条目 + 记忆；结果全部标 exploratory（完成 09-05）

## 暂缓分支

- K loop：已关闭（D-107/D-110，`k_search_closed_for_r04=true`）。
- marker-proxy 扩样本 / TSB 0.544 追查 / 未命名场检验：R-04 阴性关闭后不再作为开放任务；仅当新决策重开证据评估时恢复。
- 独立 lineage 复现：`NOT_TRIGGERED_NO_SURVIVING_CANDIDATE`（无存活候选，不再是 open TODO）。
- 强 deconvolution（cell2location/RCTD）：`NOT_TRIGGERED`（无存活候选）。

## 停止分支

- 旧 TensorFlow 计算路径：已封存（D-094），不与 Torch 数值混合。
- CPU 上的 2400-step 推断：超时不可行，遗留 `RUNNING` 标记目录不作为结果（handoff §3）。

## 分支记录

- fold-0 GPU 收尾分支（09-04 启动，D-096）：全部完成——环境/上传/双预检/K=0 端点审计（4200/7800/8400 平台，hash VERIFIED）/K=0 零 fit 推断评分（新入口 r04_torch_frozen_infer.py，D-097）/K=3 配对跑（K3−K0 均值 −28.96，fold 0 偏 K=0，与 TF 历史方向一致）/产物回传/实例退租（花费 ¥3.88）。09-05 起整组从待办移除，审计线索见本条与 `docs/roadmap.md` 2026-09-04 条目。

- 2026-09-05 午：三方向并行第一轮完成——A 置换 null（`explore_readout_null_torch_20260905.json`）+ C 深化（`explore_readout_deepdive_torch_20260905.json`）+ B 转 proposal（LEADS L-002）；新增 D-102 固定探索口径；roadmap 已记数值与假设影响。
- 2026-09-05 午（二）：TLS 共享稳健性深挖完成（`explore_shared_robustness_torch_20260905.json`；脚本+4测试）；结论见 roadmap；B（L-002）继续 hold。
- 2026-09-05 午（三）：D-103 跨设备重放处置＋外层验证脚本 direct 模式；fold-4 外层验证完成（见 roadmap）；NaN-AUC bug 修补＋回归测试。
- 2026-09-05 午（四）：fold-0 外层验证（split 拼装）+ 双 fold 合并判读完成；外层验证整体零发现。
- 2026-09-06 晨：heldout 版复核与拼装版一致，外层验证收尾；结论维持零发现。
- 2026-09-10：用户批准三线并行（marker代理＋scRNA参考＋GPU K-2自租）；B从LEADS L-002转入执行（L-002记已并入主线）；C（L-001）GPU由我租用，看守＋预算上限¥15。
- 2026-09-10 晚：K=2 Stage A 完成（D-106），K=3 默认生效，L-001 关闭；GPU 已退租对账。
## 变更记录

- 2026-09-04：建立本文件。起因：接手 R-04 Torch 交接（docs/HANDOFF_R04_TORCH_2026-09-04.md），用户批准租用 4090 实例（D-096），当前执行分支从"CPU 等待"切换到"GPU 收尾"。
- 2026-09-04 晚：环境/上传/双预检/K=0 端点审计/K=0 推断评分全部完成；K=0 原定 resume_infer 入口因 TF 专属面板门禁不可达，改用新 Torch 入口（D-097）；新增"退租实例"项（用户指示）。
- 2026-09-04 深夜：K=3 配对跑完成（K3-K0 均值 −28.96）；产物回传本地；roadmap/STATUS 已更新。
- 2026-09-05 凌晨：4090 实例已退租（CLI release-plan+release，MCP 拒放非 MCP 创建实例），花费 ¥3.88。fold-0 分支整组完成并移除（见分支记录）。新增当前分支 fold 4 配对复跑（D-098）：作业包预置完成，等实例地址；暂缓分支中的 fold-4 条目随之移除。
- 2026-09-05 早（三）：D-101 项目转探索定位，取消预注册约束（R-01/R-04 现行规则改写、D-100 规则废止）；门禁按实质重判为描述性未决，K=3 容器继续，K=2 转 LEADS L-001，GPU 审批撤回；新建 LEADS.md。
- 2026-09-05 早（二）：Torch 读出门禁判 SENSITIVE（roadmap 有完整数值与内折记录），K=2 Stage A 通道待 GPU 审批；此期间空间 null、组成拆分继续暂停。
- 2026-09-05 早：D-099 门禁开工——定位 Torch 场导出缺口（cell 内 exports 为空、TF 校验不可达），追加 D-100（含门禁规则预注册），CPU training-only 导出运行中。
- 2026-09-05 凌晨（三）：K_model=3 落定默认工程工作表示（D-099，用户指示尽快推进）；资源峰值与开销固化到 infra/r04/GPU_RUNBOOK.md；当前分支切换为结构读出 K 稳健性分析。
- 2026-09-05 凌晨（二）：fold 4 配对跑完成——K3−K0 均值 +72.73、K=3 胜 5/6，方向与 TF 历史一致，无冲突无需记 ISSUES；端点审计 VERIFIED、近精确重复；产物已回传；roadmap/STATUS 已更新。fold-4 分支整组完成并移除（见分支记录）。新增待确认分支：结构读出 K 稳健性分析。
- 分支记录补充：fold-4 GPU 复跑分支（09-05 启动并完成，D-098）——作业包预置（修正两处脚本问题：rsync 目录层级、预检字段位置）→ 实例 2 执行 → 预检 PASS → 配对跑 + 审计完成 → 产物回传 → 实例 2 已退租（本次花费 ¥1.43，含退回 ¥1.89；累计 GPU 花费 ¥5.31）。
- 2026-09-10：三线并行收官——marker-proxy smoke（a29 信号被组成调整吃掉，D-105）、NNLS 路线判死（LEADS L-003 已放弃）、K=2 配对跑完成 K=3 默认生效（D-106，L-001 关闭）；移除两个已完成的旧待办（全量场导出运行项、双 fold 合并判读项，原因：对应工作 09-05/09-06 均已完成）；D-105 重号修正为 D-106；GPU 已退租对账（本次¥3.12，累计¥8.43）。
- 2026-09-11：立项 R-04 完成缺口清单（对照完成判据逐项盘点，最大缺口为外层零发现）；首项为路线三选一，需用户拍板。
- 2026-09-11：独立节点报告 docs/PHASE_REPORT_R04_EXPLORATORY_2026-09-10.md 落盘并推送；STATUS 末段同步。
- 2026-09-11：执行 R04 finalization prompt，R-04 阴性关闭（D-107…D-110，最终门禁落盘）；K 控制面修复；nested 终审推翻旧 smoke 结论；文档原位更新；final phase report 落盘；313 测试全过并推送（db625ed）。
- 2026-09-11：R-03 普查完成（22 多切片块零已知间距；I-002 对 R-11 升级硬阻塞）；是否申请新数据待用户决策。
- 2026-09-11：R-03 工作冻结为交接文档 docs/HANDOFF_R03_2026-09-11.md（普查＋侦察＋Cervilla 落地＋A/B/C/D 待决策＋恢复步骤）；恢复时只读该文档即可。
- 2026-09-12：数据补充阶段关闭——W0（100文件16G）＋P0/P1（312文件，W0去重后新增224，落盘52G）＋HT480B1两对CANDIDATE（32文件3.55GB）全部md5/文件头校验落盘；W1分五档（P0 48f/6.5GB、P1 264f/35.6GB、P2_multi 216f/29GB、P2_single 432f/58GB、P3 96f/13GB），P2_multi验完（仅HT480B1两对判CANDIDATE，HT339B1判BORDERLINE不动，168文件判不同块REJECT）；P2_single＋P3跳过（~71GB零价值）；CRDC/Gen3/NCI、dbGaP、OHSU保持外部阻塞；流量配置为API走代理＋S3/GCS bulk直连（download_env.sh）。后续只按触发条件补数，不再囤积。（本条目曾重复登记，09-14 清理为单条）
- 2026-09-12：用户定方向为关闭 R-04（收缩，不重开证据评估）；R-11 相关新数据（Xenium CRC／dbGaP／作者邮件）一律暂缓，需另行明确批准。
- 2026-09-14：场发现策略修正（D-112）——残余降级为评分维度、field registry 立项、TLS 阳性对照结果与 panel 盲区（I-020）落盘；文档已原位更新（roadmap/decisions/ISSUES/LEADS/STATUS）。
- 2026-09-14：发现引擎定型（D-113）——组成pattern普查＋距离函数双车道（R-16 新节点）；连续信号null换代登记（I-021）；TODO下一步改为先交一页设计方案。registry schema并入设计方案，不再单列。
- 2026-09-14：R-16 设计方案 v1 落盘（docs/R16_DISCOVERY_DESIGN_20260914.md）；侦察确认 bash 环境 sklearn/leidenalg/igraph 齐备、坐标在源 h5ad obsm/spatial 可按 barcode 对齐、坐标单位不可靠故距离用图跳数。待用户批准后执行第①步普查。
- 2026-09-14：R-16 两步走全部落地——D-114 基因宇宙（整合 HVG top-10000，14/14 TLS 基因救回）、D-115 两层结构（Tier-1 轴尺子厂＋Tier-2 残余普查）、registry v1 228 行交付；用户第一性原理追问驱动的逻辑层重构（L1-L5）全部处置。
- 2026-09-15：D-116 独立复现完成——ST-CRC 7/7、USZ 8/8 轴相干；CXCL13 峰≤+1 为 6/7 与 7/8；USZ TLS~B 中位 AUC 0.71（肾高肺杂）；不混发现集。
- 2026-09-15：D-117 退役ILC轴；工作集改为T/B/Mye/Epi/Stromal/Plasma；未测的x（内皮/增殖/缺氧）不得默默加入。
- 2026-09-15：补Tier-2人读摘要——221组漏斗摊开：13个干净多患者组，其中7个六轴解释不了。
- 2026-09-15：D-118 Tier-1 v2 尺子工厂交付——邻域增强×λ→PCA候选轴→LISA tile（999 次置换）→新旧对照→loading 匹配；1880 候选→857 组，仅 2 组复现且均为 KNOWN，NEW 复现轴=0（可信阴性）；registry v1.1 共 1085 行；LEADS L-008 记录旋转歧义限制。
- 2026-09-15：D-119 联合嵌入轨全量完成——Arm A/B 各99/117候选，split-half双半匹配；67个≥2病人候选与Tier-2 marker交集<0.3（杯状/Paneth样/Wnt/MHC-II等）；registry共1301行；R专用臂进行中。
- 2026-09-15：D-121 GraphST GPU 臂准备完成——逐片训练（merged稠密矩阵约100GB不可行）、全600 epoch、Leiden代mclust、V100券优先；代码＋数据＋测试本地全绿，待用户给登录方式开跑。
- 2026-09-16：用户指令全量不限时跑完即停（D-122）；GraphST 47片×600epoch已开跑，跑完回收＋退租。
- 2026-09-16：GraphST GPU臂全量完成——47片×600epoch，356组32复现入库；registry共1657行；退租关闭。
- 2026-09-16：V100 实例 lyg2143 已退租（0 在跑，0 残留磁盘，未用时长退款 ¥1.10）；SSH 密码文件已销毁。
- 2026-09-18：TODO 顶部清理——旧 R-04/R-16 完成项归入历史正文；当前工作面切换为恢复包 T0–T6（分支 recovery-20260918）。起因：外部审计包指出 Lane A/B 与 PC 工厂的方法实现问题（A02/A05/A06/A09/A10/A13），旧否证外推待复核。
