# R-04 最终收尾 Agent Prompt｜2026-09-11

## 0. 角色与总目标

你正在维护仓库 `FreddieWho/013_spatial`。你的任务不是重新设计 R-04，也不是继续寻找一个哲学上唯一的 K，而是在**严格继承现有 plan / roadmap / decisions / issues / STATUS 与 2026-09-10 exploratory phase report 的前提下，完成 R-04 的最终科学收尾、控制面修复和机器可审计闭环**。

本轮必须把本 Prompt 中经仓库事实验证后成立的规则和状态，**原位更新回现有项目文档**；不得另起一套 plan/roadmap 造成断裂。允许新增一份最终 phase report 作为历史快照，但它不能替代现有文档体系。

开始前先读并遵守：
- `AGENTS.md`
- `STATUS.md`
- `TODO.md`
- `docs/plan.md`
- `docs/roadmap.md`
- `docs/decisions.md`
- `docs/ISSUES.md`
- `docs/PHASE_REPORT_R04_EXPLORATORY_2026-09-10.md`
- `infra/r04/README.md`
- R-04 当前脚本、测试、运行产物与 provenance/hash 记录

若本 Prompt 中任何事实与当前仓库 HEAD 冲突，以可审计仓库事实为准，但必须在最终报告中指出冲突，不能静默改写历史。

---

## 1. 本轮不可改变的科学语义

### 1.1 K 的语义已经关闭，不再开启 K 搜索循环

继续沿用当前 `docs/plan.md` 中已经确立的定义：

- `K_model`：模型允许的潜在空间分子效应表示容量；
- `K_eff`：由留出预测、重复性、空间 null 与子空间稳定性共同支持的有效维度；
- K **不是** TLS、血管、坏死、肿瘤—基质边界等生物结构的数量；
- 生物结构与潜在维度允许多对多映射；
- 单个 factor 的编号、符号和旋转不是可唯一识别的生物对象，主要解释对象应为稳定子空间、旋转不变效应和结构特异 readout。

结合截至 2026-09-10 的实际 K=2/K=3 证据，本轮默认采用以下项目级工作状态，除非仓库中出现更新的反证：

```text
working_k_model = 3
primary_readout_rank = 2
global_k_eff = NOT_IDENTIFIABLE
selected_k = null
k_search_closed_for_r04 = true
third_direction = NO_STABLE_POSITIVE_INCREMENT
```

解释必须精确：
- fold 0 中 K=2 相对 K=0 仍为负，但比 K=3 更不差；
- fold 4 中 K=2 与 K=3 的 held-out 增益近似；
- 因此第三方向**没有显示稳定的正增量价值**；
- 这支持把 K=3 当作 overcomplete working container，并把主要下游 readout 限制在稳定 rank≤2；
- 这**不等于**证明真实组织的全局 `K_eff=2`，更不能把 `selected_k` 改成 2 或 3。

### 1.2 禁止本轮继续做的 K 工作

默认禁止：
- 新一轮 K=0/1/2/3/4/更大 K 搜索；
- 为了补齐旧 six-cell Stage A 协议而机械补跑 K=2 额外 seed；
- 继续延长 K=3 优化步数；
- 因为 fold 4 某个 rank3 指标看起来较高而专门追第三方向；
- 用新的 K 计算推迟 R-04 收尾。

只有当本轮最终 readout 对 rank≤2 与 full K=3 的结论产生**明确、可重复、方向相反的科学结论**时，才允许升级为新的明确决策；不得自行启动计算。

---

## 2. 先完成一次“现状—实现—文档”审计

在修改代码前，制作一个短的内部审计表，至少核对：

1. `docs/PHASE_REPORT_R04_EXPLORATORY_2026-09-10.md` 的每条关键结论能否追溯到实际 artifact；
2. K=2 GPU 运行是否来自当前 canonical `r04_real_k_search` 路线，输入 hash、fold、gene split、fit/inference step、checkpoint audit 是否齐全；
3. 当前五 seed K=3 稳定性证据与 K=2 paired evidence 的统计语义是否一致；
4. marker-proxy composition smoke、NNLS dead-end、outer-validation、label permutation null 是否各自有对应代码/产物；
5. 哪些脚本、README、ISSUES 仍停留在 2026-09-01 或 CPU-blocked 旧状态；
6. 当前测试覆盖哪些关键语义，哪些只是历史协议测试；
7. 是否存在任何“报告已关闭但机器入口仍显示 blocked/not-run”的控制面冲突。

审计结果必须进入最终 phase report 的“provenance / continuity audit”部分，但不要把审计表变成新的 roadmap 节点。

---

## 3. 必须修复的控制面断裂

### 3.1 K 状态统一

重点检查并修复：
- `scripts/r04_finalize_k_semantics.py`
- `scripts/r04_k2_bridge_stage_a.py`
- `infra/r04/README.md`
- `docs/ISSUES.md` 中与 K 未闭合相关的旧条目
- 对应 tests

已知风险：旧 finalizer/Stage-A 逻辑可能仍硬编码 `BLOCKED_CPU_RUNTIME`、CPU-only、`COMPLETE_STAGE_A_NOT_CLOSED` 等旧状态，而真正的 K=2 证据已经通过 2026-09-10 GPU `real_k_search` 路线完成。

要求：
- 历史 artifacts 不删除、不伪造、不回写；
- 旧 bridge driver 如果仍需要历史可重放性，可以保留，但要明确标成 legacy/historical，不得继续作为当前 K 状态权威入口；
- 当前 K 状态必须有**一个** canonical、可重放、机器可读的聚合入口；
- 该入口必须直接绑定实际 K=2/K=3 evidence artifact、hash/provenance，并输出第 1.1 节的工作语义；
- 更新相关测试，禁止未来再次出现文档“已关闭”而机器结果“CPU blocked”的矛盾。

建议优先复用并升级 `scripts/r04_finalize_k_semantics.py`，而不是再造多个 K 总结脚本。

### 3.2 append-only 决策澄清

不要修改历史 decision。追加新 decision，至少澄清两点：

1. 对 D-106 的精确解释：
   - 第三方向没有稳定正增量；
   - fold 4 K2≈K3；
   - fold 0 K2 明显优于 K3，但 K2/K3 都不优于 K0；
   - 因而采用 `working_k_model=3 + primary_readout_rank≤2`，但不声称 `global K_eff=2`。

2. 如果 D-092 曾冻结 six-cell K2 Stage A，而实际 D-101/D-104/D-106 在 exploratory scope 下以更小 GPU panel 关闭 K 问题，则明确记录：
   - 哪个后续 decision 正式 supersede 了哪个旧执行范围；
   - 为什么不再补跑旧 Stage A 的缺失 seed；
   - 该压缩只支持“停止 K 搜索/选择工作表示”，不升级为确认性全局 K 结论。

---

## 4. 修正 R-04 的最终完成定义：允许阳性，也允许可信阴性

当前 R-04 的完成判据如果仍要求“至少一个候选空间场必须成功”，会让一个科学上有效的阴性结果永远无法完成 R-04，并诱导继续 fishing。

因此，在**不改变核心科学问题和 H-01/H-03/H-05 含义**的前提下，更新 `docs/roadmap.md` 的 R-04 completion gate 为 outcome-neutral：

### A. `R04_COMPLETE_WITH_REPRODUCIBLE_RESIDUAL_FIELD`
至少一个候选满足：
- 患者/组织块层面重复；
- 在 rank≤2 主表示下成立；
- 组成控制后仍有增量残余；
- 有可审计的不确定性；
- 必要时通过空间 coordinate-permutation/null；
- 对命名结构，结构 GT 只用于解释/验证，不反向定义候选场。

此时冻结 surviving candidate(s)，允许进入 R-05；后续 R-06/R-07 按原 roadmap 条件推进。

### B. `R04_COMPLETE_NO_REPRODUCIBLE_RESIDUAL_FIELD`
在预先限定、现有数据的审计范围内，没有任何候选同时通过重复性、组成残余和不确定性要求。

此状态是**R-04 已完成的阴性结果**，不是“继续找更多 K/更多结构/更多参数”的理由。
- 对 H-01/H-03/H-05 的信心按 `plan.md` 的否证逻辑下调或记为当前数据不支持；
- R-05/R-06/R-07 的强反演主线记为 `NOT_TRIGGERED_NO_SURVIVING_R04_FIELD` 或等价状态；
- 不得用外部数据、第三方向或新 marker 做事后救援，除非另有用户批准的新项目决策。

### C. `R04_BLOCKED_IDENTIFIABILITY`
仅在**连阴性结论都无法可靠判断**时使用，例如关键输入缺失、composition control 无法构造且会决定主结论、或 provenance 破坏。

必须明确：
- 缺什么；
- 为什么现有数据无法给出阳性或阴性判断；
- 最小解锁数据/计算是什么；
- 不得泛化成“再找更多数据可能更好”。

`docs/plan.md` 当前关于 K 和多结构语义已经基本正确，**默认不修改**。只有科学问题/假设本身改变时才修改，并先在 `decisions.md` 留下理由。

---

## 5. 用现有数据完成 R-04 最后一次有界科学审计

本轮默认：**不下载新外部数据、不申请 GPU、不训练新的 K 模型**。

### 5.1 Composition closure：先修 nested cross-fitting

检查当前 marker-proxy composition adjustment 的实现。已知需要重点排查：composition residualization 是否先在全体患者上 cross-fit，然后再做结构 readout 的内层 patient split；如果是，这会让一个最终 test patient 参与其他 train patients 的预处理映射，形成轻度 transductive/nested-CV 泄漏。

必须把 composition adjustment 嵌套到每个结构 readout 的 outer/inner patient split 中：

- 只用当前 train patients 拟合 `latent/readout feature ~ composition proxies` 的关系；
- 用这一个 train-fitted mapping 同时变换 train/test；
- test patient 不得参与任何 residualizer 参数估计；
- 添加专门单测，能在故意注入 test-only composition signal 时捕获泄漏；
- 旧 exploratory artifacts 保留，不覆盖；重新生成 final-audit artifacts。

然后在所有**当前已经可用且不需要新模型训练**的 training/held-out patient/section 上做 bounded composition audit：

- rank1；
- rank2（主分析）；
- full K=3 只作 sensitivity；
- unadjusted vs nested composition-adjusted；
- marker proxies 必须继续称为 proxy，不得写成真实 cell fraction；
- ILC 等已知 unusable 类别不得强行补齐；
- NNLS 失败路线保持 abandoned，不继续调 NNLS。

统计单位优先 patient/block，输出：
- patient-level effect；
- bootstrap/empirical CI 或同等级不确定性；
- 调整前后方向、效应量和不确定性；
- 是否存在 composition adjustment 后仍稳定的 residual signal。

如果没有候选 surviving，**不要**自动升级到 cell2location/RCTD 或重新构造 scRNA deconvolution；先进入阴性 closure 判断。

### 5.2 Unnamed field reproducibility closure：回到 R-04 原始主对象

R-04 的主对象不是 TLS/TSB，也不是 K，而是“无需结构标签定义的潜在空间分子场”。必须用现有 fitted artifacts 做一次最终、rotation-invariant 的 unnamed-field audit。

要求：
- 不用 TLS/TSB label 参与候选发现、对齐或阈值选择；
- 不按 factor 1/2/3 名称硬匹配；
- 基于稳定 rank≤2 子空间、gene-effect/loading 子空间、principal angle/canonical correlation 等旋转不变量做重复性判断；
- 不把不同患者的绝对空间坐标强行对齐；空间形状只使用患者内可比较的 summary，不制造跨患者伪几何对应；
- 跨 fold、restart、患者/组织块报告稳定性；
- 每个候选最终只允许：
  - `SURVIVES`
  - `DOES_NOT_SURVIVE`
  - `NOT_IDENTIFIABLE`
- 不命名未锚定场，不用外部常识补标签。

优先复用现有 null/稳定性框架。若当前已有证据已足以判定某 candidate 不稳定，不为它新增计算。

### 5.3 Anchored readout synthesis：TLS/TSB 只做解释与验证

整合：
- fold 0 / fold 4 existing outer validation；
- composition-adjusted final audit；
- rank1/rank2/full-rank sensitivity；
- label permutation null；
- patient-level uncertainty。

特别规则：
- fold 4 中曾出现的 TSB full/rank3 AUC≈0.544，如果 rank1/rank2、fold0 或 composition-adjusted 结果不复现，默认记为 `ISOLATED_EXPLORATORY_SIGNAL`；
- 不为了“解释 0.544”而追第三方向；
- TLS 若每 fold 只有一个 section，只能保留 descriptive，不升级为 patient-level reproduction。

---

## 6. Conditional work：只有 surviving candidate 才触发

只有第 5 节出现至少一个 `SURVIVES`，才评估以下任务是否为最终主张不可替代：

### 6.1 更强 composition control
仅当 marker-proxy residualization 无法排除 composition shortcut 且该问题决定 candidate 是否 surviving 时，才提出更强 deconvolution。
- 不继续失败的简单 NNLS；
- 优先当前已有本地 scRNA reference 与经过验证的工具；
- 若需要新大规模计算/GPU，先形成最小申请，不自动执行。

### 6.2 Fixed-cell-type state decomposition
仅当现有数据可以在不循环定义、不泄漏 GT 的条件下完成时执行。
若无法可靠完成，允许输出：
`STATE_COMPONENT_NOT_IDENTIFIABLE_CURRENT_DATA`
并收缩 H-03 的措辞；不要为了勾选 roadmap 而伪造 cell-state 证据。

### 6.3 Coordinate-permutation spatial null
注意：结构 label permutation null 不等于“空间场存在性”的 coordinate-permutation null。
仅对 surviving unnamed/residual field 跑真正的坐标置换/null，以验证空间坐标贡献不是普通低维表达结构或平滑伪影。

### 6.4 Independent lineage reproduction
只有 surviving candidate 且主张需要跨 lineage 支撑时触发。
若没有 surviving candidate，记：
`NOT_TRIGGERED_NO_SURVIVING_CANDIDATE`
而不是继续作为 R-04 open TODO。

任何新外部数据获取、GPU 租用或不可逆投入仍按仓库现有用户批准规则执行。

---

## 7. 本轮明确不要求完成的旧分支

除非第 6 节被 surviving candidate 正式触发，否则不要为了“完成 R-04”去完成：

- signed NB-GP adapter 的完整后验统一；
- 失败 NNLS 的继续优化；
- USZ/common-panel 更深技术补强；
- 新的 K2 seed panel；
- K>3 搜索；
- rank3-only 结构信号救援；
- 新 external cohort 下载；
- R-05/R-06/R-07 的实质工作。

这些应根据 final R-04 gate 被标为：resolved / deferred / not-triggered / limitation，而不是永远悬挂为“还没做完”。

---

## 8. 创建真正的 R-04 最终机器门禁

当前 synthetic `r04_scientific_gate.py` 只证明方法校准，不应被误用成真实 R-04 closure gate。

新增或复用一个明确的最终聚合入口，例如：

`scripts/r04_finalize_phase.py`

它必须：
- 只读取已登记的 final-audit artifacts；
- 校验 evidence paths、hash/provenance、patient/fold independence；
- 不重新训练模型；
- fail closed；
- 生成类似：
  `infra/r04/r04_final_gate_20260911.json`

至少输出：

```text
phase = R04
status = one of:
  R04_COMPLETE_WITH_REPRODUCIBLE_RESIDUAL_FIELD
  R04_COMPLETE_NO_REPRODUCIBLE_RESIDUAL_FIELD
  R04_BLOCKED_IDENTIFIABILITY

working_k_model = 3
primary_readout_rank = 2
global_k_eff = NOT_IDENTIFIABLE
selected_k = null
k_search_closed_for_r04 = true
surviving_candidates = [...]
composition_gate = ...
unnamed_field_reproducibility_gate = ...
anchored_readout_summary = ...
uncertainty_gate = ...
coordinate_null = PASS / FAIL / NOT_TRIGGERED
independent_lineage = PASS / FAIL / NOT_TRIGGERED
next_nodes_authorized = [...]
blocked_or_not_triggered_nodes = [...]
evidence = [{path, sha256, role}, ...]
```

最终 gate 不能通过手工改 JSON 得出，必须由代码从来源 artifact 重算。

---

## 9. 文档必须原位更新，保证方案连续性

### `docs/plan.md`
- 默认不修改；
- 保留当前 K_model/K_eff、多对多结构语义和 factor 不可唯一命名规则；
- 只有核心科学假设真的改变时才改，并在 decisions 先说明。

### `docs/roadmap.md`
必须修改：
- 把 R-04 完成判据改成第 4 节 outcome-neutral 三状态；
- 明确 R-05/R-06/R-07 是否被 surviving candidate 触发；
- 更新 R-04 最终状态；
- 不新增一个平行“R-04b”或新的 K 节点。

### `docs/decisions.md`
只追加：
- K closure 精确语义；
- D-092 six-cell protocol 与后续 exploratory GPU closure 的 supersession 关系；
- nested composition crossfit 修复；
- final R-04 gate 与最终 closure 结果。

### `docs/ISSUES.md`
逐项清理至少 I-012–I-019 中与 R-04 相关的状态：
- K 未关闭的旧 I-019 必须与当前 K closure 对齐；
- composition/state/uncertainty/signed/common-panel 等 issue 要按 final gate 写成 resolved / limitation / deferred / not-triggered；
- 不允许已经不再阻塞 R-04 的旧问题继续伪装成 hard blocker。

### `STATUS.md`
按 `AGENTS.md` 四段规则整份重写：
1. 我们在回答什么问题；
2. 目前知道了什么；
3. 今天做了什么；
4. 接下来做什么。

必须用普通语言明确：K 搜索已结束；R-04 是阳性完成、阴性完成还是 identifiability blocked；下一节点为什么被允许/不被允许。

### `TODO.md`
- 删除/关闭已经完成的 K loop；
- independent lineage、强 deconvolution 等按触发条件处理；
- 最终只保留真正尚需执行的节点。

### `infra/r04/README.md`
必须更新旧 CPU/GPU/K2 状态，确保与当前 canonical evidence 和 final gate 一致。

### Phase report
保留：
`docs/PHASE_REPORT_R04_EXPLORATORY_2026-09-10.md`
作为不可改写的历史快照。

R-04 最终关闭后新增：
`docs/PHASE_REPORT_R04_FINAL_2026-09-11.md`

它是报告，不是新 roadmap。

---

## 10. 测试与审计要求

至少补充/更新测试覆盖：

1. K finalizer 不再把当前状态误报成 CPU blocked；
2. historical K2 bridge 与 current canonical closure 的角色不会混淆；
3. nested composition residualization 不读取 test-patient 信息；
4. factor permutation/sign/rotation 不改变 unnamed-field reproducibility 结论；
5. final R-04 gate 对缺失 artifact/hash mismatch/fold overlap fail closed；
6. negative outcome 可以合法生成 `R04_COMPLETE_NO_REPRODUCIBLE_RESIDUAL_FIELD`；
7. no surviving candidate 时 external lineage / coordinate null / stronger deconvolution 正确标为 NOT_TRIGGERED；
8. surviving candidate 时必要 conditional gate 不能被静默跳过。

运行当前完整测试集。不要把 phase report 中历史记录的测试数量作为必须固定的数字；报告实际 HEAD 下的测试结果和失败数。

---

## 11. 停止规则

本轮的优先级是**完成 R-04，而不是把所有可能分析做完**。

只要能够从现有数据可靠进入第 4 节 A 或 B，就立即关闭 R-04。

禁止因为：
- AUC 接近 0.5 但略高；
- 某一 fold/rank3 有孤立改善；
- “也许换一个 deconvolution 更好”；
- “多跑几个 seed 会更放心”；
- “外部数据可能救回来”

而无限延长 R-04。

若只能进入 C，必须给用户一个**最小、明确、可否决**的解锁申请，而不是自行扩张项目。

---

## 12. 最终交付

完成后必须提交：

1. 代码与测试修改；
2. canonical K closure machine-readable artifact；
3. final R-04 machine-readable gate；
4. `docs/PHASE_REPORT_R04_FINAL_2026-09-11.md`；
5. 原位更新后的 roadmap / decisions / issues / STATUS / TODO / infra README；
6. 一个非常短的最终回复，格式如下：

```text
R-04 final status: <三状态之一>

K:
- working K_model = 3
- primary readout rank = 2
- global K_eff = NOT_IDENTIFIABLE
- selected_k = null
- K search = CLOSED

Scientific result:
- surviving residual fields: <n + names/IDs or none>
- composition-residual evidence: <一句>
- unnamed-field reproducibility: <一句>
- anchored TLS/TSB result: <一句>

Drift/repairs:
- <只列真正修复的 2–5 项>

Next:
- <被授权的下一 roadmap 节点，或为什么不触发>

Tests:
- <passed/failed>
```

不要用“基本完成”“接近完成”“还可以进一步优化”这类模糊结论。必须给出明确的 R-04 final state。
