# 阶段报告：R-04 探索阶段收官（2026-09-10）

本文件是 R-04 探索阶段（D-101 转探索定位之后）三线工作的交接快照：C 线 K=2 裁决配对跑、A 线 marker-proxy 组成路线、B 线 scRNA 参考＋NNLS 路线，以及此前完成的外层验证合并判读。面向接手研究的读者：只靠本文件与所引用的产物，即可理解每条线做了什么、得到什么数、下什么结论、证据等级是什么。本文是时间点快照，不是持续维护文档；项目的五份维护文档见 `AGENTS.md`（`STATUS.md`、`docs/plan.md`、`docs/roadmap.md`、`docs/decisions.md`、`docs/ISSUES.md`）。R-04 节点本身**仍未完成**，缺口见 §7。

证据等级总声明：本阶段全部新结果均为探索性/描述性（exploratory），小样本（6 患者/fold 配对分数、N=2 读出对照），不改变 `docs/plan.md` 中任何假设的信心；`selected_k=null` 保持。

## 1. 范围与状态

- **科学背景**：检验"无清晰边界区域的分子信号能否形成可重复的连续空间场"（H-01…H-07 见 `docs/plan.md`）。R-04 内部顺序为候选场发现→组成拆分→独立重复→结构锚点验证→交 R-05（见 `docs/roadmap.md` R-04 节）。
- **本阶段三线**（2026-09-10 用户批准并行开工）：
  - A marker-proxy 组成路线：三库投票构建逐 spot 伪组成，经 crossfit 调整后重跑读出（smoke，fold-0）。
  - B scRNA 参考＋NNLS 路线：GSE236581 参考构建＋两版 NNLS 反卷积，按预定判据判死。
  - C GPU K=2 Stage A（LEADS L-001 转入执行，D-104/D-106）：fold-0/fold-4 的 K=0/K=2 同协议配对，裁决第三方向。
- **此前已完成并在本阶段合并判读**：Torch 结构读出 K 稳健性门禁、外层验证双 fold（fold-4 TSB rank3 弱阳性 0.544 孤立、fold-0 全 null，合并判零发现）、heldout-full 复核一致（差 <0.005）。

## 2. 关键数值结果

### C 线：K=2 裁决（`torch_representative_fold{0,4}_20260910_k02`，4/4 FIT_AND_SCORED，inference 双 split 收敛，4 份端点审计）

| fold | K2−K0 均值 | K2 胜 | K3−K0 均值（既有） | 方向 |
|---|---|---|---|---|
| fold-0 | −16.61 | 2/6 | −28.96 | 一致为负 |
| fold-4 | +71.96 | 6/6 | +72.73 | 几乎逐位一致 |

结论（D-106）：第三方向零增量贡献；fold 异质性位于 ≤2 有效维度内；K=3 作为过完备工作表示正式生效（K=3 解空间包含 K=2，切换无收益；读出仍用 rank-2 口径，第三方向不参与判决）。

### A 线：marker-proxy smoke（`explore_marker_smoke_fold0_20260910.json`，fold-0，N=2 患者对照）

- 三库投票可用基因：T22 / B11 / Mye16 / Epi54 / Stromal22；ILC 仅 1 基因，判不可用。
- 组成调整后关键变化全在 a29 患者：full-TSB +0.127→−0.068、rank2-TSB +0.127→−0.004（门禁讨论的单患者信号被吃掉）、TLS 转明显负值；0bd3 患者小幅正向或基本不动。
- 解读（描述性）：a29 特异信号与组成共线，"相容于组成驱动"，不定案（反向因果警告：回归分不清组成驱动与真信号恰共线，即 H-03 纠缠问题）。

### B 线：NNLS 判死（LEADS L-003，已放弃）

- 两版（全基因 3854 / DE 限定 512）免疫份额均与经典 marker 零/反相关：LYZ–Mye −0.05→−0.19，PTPRC–T −0.02→−0.04；Epi/Stromal 始终 +0.55–0.65 正常。
- 根因：T 与 ILC 均值 profile 余弦 0.991（数学上不可分）＋scRNA/Visium 平台 gap。按预定判据两版皆不过，不再投入；脚本与参考保留为负结果资产。

### 外层验证合并（exploratory，9 患者 TSB）

- fold-4：TSB rank3 0.544（N=4，孤立弱阳性），rank1/2 null；fold-0：TSB 三 rank 0.486/0.500/0.480 全 null；TLS 留出两 fold 各仅 1 section，只报数。
- 合并判读：零发现；fold-4 rank3 不升级。heldout-full 版与拼装版差 <0.005，无分歧。

## 3. 规范产物（下游契约表面）

### 脚本（`scripts/`，均有对应测试，298 项全绿）

| 脚本 | 用途 |
|---|---|
| `r04_build_composition_reference.py` | GSE236581 参考构建（含 provenance） |
| `r04_deconvolve_nnls.py` | 两版 NNLS 反卷积＋预定判据 |
| `r04_combine_marker_proxy.py` | 三库投票＋伪组成构建 |
| `r04_explore_marker_smoke.py` | crossfit 组成调整＋读出重跑 |
| `r04_explore_outer_validation.py` | 外层验证（含 direct/split-direct 双模式＋重建自洽门禁） |
| `r04_explore_readout_null.py` | 空间置换 null（200 draws） |
| `r04_explore_readout_deepdive.py` | rank 曲线＋绝对 AUC＋范数＋逐内折 |
| `r04_explore_shared_robustness.py` | TLS 共享生态稳健性深挖 |
| `r04_export_heldout_only.py` | heldout-only 导出 |

### 数据产物（`infra/r04/`，大文件与 checkpoint `.pt` 按 `.gitignore` 排除在外）

- `marker_proxy_combined.json` — 投票 marker 集（含未映射标签与有意丢弃类别记录）。
- `composition_ref/provenance.json` — GSE236581 参考构建 provenance（`.npz` 矩阵 443MB 本地保留，不进仓库）。
- `explore_*.json` — null / deepdive / shared_robustness / outer_validation（fold-0、fold-4、heldout-full）/ marker_smoke 结果。
- `torch_k{0,2}_fold{0,4}_checkpoint_audit_20260910.json` — 4 份 K=2 端点审计。
- `gpu_k2_20260910/run_k2.sh` — K=2 作业脚本（预检 hash 锚点 input 044d28c5 / objective-input b5a3fb68 / manifest 5a0dc1de）。
- `training_structure_readout_torch_20260905.json`、`torch_effect_export_panel_20260905.json` — 读出与导出面板。

## 4. 决策链（`docs/decisions.md`）

- D-099：K=3 为默认工程工作表示（`selected_k=null`）；D-101：项目转探索定位，预注册约束取消（D-100 规则废止），GT 永不进入输入与训练永久保留。
- D-102：探索口径（结果标 exploratory，不写确认性结论）；D-103：跨设备重放场接受条件（终点损失相对差 ≤4e-4）。
- D-104：GPU 自租用授权（4090 档，看守＋¥15 预算上限）；D-105：marker-proxy 路线授权（R-02 例外事项记录：proxy 仅作组成协变量，不作结构 GT）。
- D-106：K=2 裁决完成，K=3 默认转正；失效条件为未来某 fold 出现 K=3 显著优于 K=2 或读出依赖第三方向。

## 5. GPU 记账

| 会话 | 实例 | 花费 | 状态 |
|---|---|---|---|
| 09-04 fold-0 | 4090 24GB | ¥3.88 | 已退租 |
| 09-05 fold-4 | 4090 24GB | ¥1.43 | 已退租 |
| 09-10 K=2 | 4090PLUS 48GB（¥2.49/h，约1h） | ¥3.12（Money 0.39＋Power 2.73） | 已退租，无残留磁盘 |
| **累计** | | **¥8.43** | running=0，kept disk=0（双重核验） |

经验：4090PLUS 的 torch cu124 原生可用（本次首次验证）；其 quote TTL 仅 60 秒，下单须 plan→rent 连打。作业手册见 `infra/r04/GPU_RUNBOOK.md`。

## 6. 对假设的影响

H-01/H-03/H-04/H-05 信心均不变。K 环关闭是表示容量问题，不是生物学发现；外层零发现意味着训练侧信号尚未拿到独立重复证据；a29 的组成共线提示 H-03 纠缠真实存在，但 N=2 不定案。

## 7. R-04 完成缺口（见 TODO.md"下一步"分支）

1. 路线三选一（需用户拍板）：外部数据补充 vs 结论收缩 vs 现有数据深挖。
2. marker-proxy 扩样本验证；TSB fold-4 rank3 的 0.544 追查（两者可并行）。
3. 组成与状态分量＋不确定性量化（完成判据要求，目前缺）。
4. 未命名场的可重复性检验（目前只验过有锚点的场）。
5. 独立 lineage 复现（暂缓，需 GPU 审批）。

## 8. 阅读顺序

1. `STATUS.md` — 两分钟通俗现状。
2. `docs/plan.md` — 科学问题与判定标准。
3. `docs/roadmap.md` R-04 节 — 节点判据与日期条目（含 09-10 三条）。
4. 本文件 — 三线接口与产物细节。
5. `docs/decisions.md`（D-099…D-106）、`LEADS.md`（L-001 已并入主线、L-003 已放弃）。
