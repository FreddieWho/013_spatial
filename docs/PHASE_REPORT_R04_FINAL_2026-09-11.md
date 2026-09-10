# 阶段报告：R-04 最终关闭（2026-09-11，阴性完成）

本文件是 R-04 最终收尾的历史快照，记录关闭 R-04 的证据、规则、修复与判定。`docs/PHASE_REPORT_R04_EXPLORATORY_2026-09-10.md` 保留为不可改写的上一阶段快照。最终机器门禁：`infra/r04/r04_final_gate_20260911.json`（由 `scripts/r04_finalize_phase.py` 从登记产物重算，fail-closed），最终状态：**`R04_COMPLETE_NO_REPRODUCIBLE_RESIDUAL_FIELD`**。

## 1. Provenance / continuity audit（Prompt §2 审计结论）

1. Exploratory phase report 的结论全部可追溯到产物，但发现一处转录冲突：roadmap 09-10 条目与 phase report §2 写"full-TSB +0.127→−0.068"，而产物 `explore_marker_smoke_fold0_20260910.json` 中 a29-test-fold 的 full-TSB unadjusted 实际为 −0.011（+0.127 属 rank2-TSB unadjusted）。新审计以产物值为准；历史文字不改写（D-108 已记录）。
2. K=2 GPU 运行走 canonical `r04_real_k_search` 路线：4/4 cell `FIT_AND_SCORED`，fit 8400 步，inference 双 split 收敛，input/objective hash 锚点（fold-0 预检锚；fold-4 为跨 cell 一致＋审计 VERIFIED 绑定），4 份 checkpoint 审计齐全（D-106 数值逐位复算吻合）。
3. 五 seed K=3（TF，restart 0–4）与 K=2 配对（Torch，restart 0）统计语义不一致处已如实记录：K=2 证据是单 restart；跨 seed 第三方向不稳定性引自 TF 面板。两者方向一致（第三方向无稳定正增量），闭合结论只做到"停止 K 搜索"，不升级为确认性全局 K 结论。
4. smoke / NNLS / outer-validation / null 各有对应代码与产物（见 exploratory report §3）。
5. 过期状态已修复：`r04_finalize_k_semantics.py`（旧输出仍称 `BLOCKED_CPU_RUNTIME`）、`k2_bridge_stage_a`（记 legacy）、`infra/r04/README.md`（CPU-blocked/GPU-禁止旧段落）、ISSUES I-012–I-019。
6. 测试：旧 finalize 测试断言 BLOCKED 状态（保留为历史语义回归）；新增 closure/nested-leak/invariance/gate 测试 15 项。
7. 控制面冲突即第 5 条，已修复；修复后机器入口与文档一致。

## 2. 控制面修复

- K 统一（D-107）：`scripts/r04_finalize_k_semantics.py --mode closure` 为唯一 canonical 入口，输出 `infra/r04/k_closure_canonical_20260911.json`（16 个证据绑定，D-106 均值 fail-closed 重算）；旧 handoff builder 保留可重放；`r04_k2_bridge_stage_a.py` 头部标 legacy，行为不动。
- 修复中发现的仓库事实：fold-4 cell 的 input/objective hash 与 fold-0 不同属正常现象（不同 split 的输入本就不同）；闭合逻辑对 fold-4 采用跨 cell 一致＋审计 VERIFIED 绑定，而非 fold-0 前缀锚。

## 3. Composition closure（D-108/D-109）

旧 smoke 的 pooled 拟合构成 transductive 泄漏（内层训练特征经由含测试患者的 residualizer 计算）。Final audit（`scripts/r04_composition_final_audit.py`，产物 `composition_final_audit_20260911.json`）把拟合器嵌套进每个内层 split（只用训练患者拟合，同一映射变换训练/测试），rank1/2/full，200-draw bootstrap CI。结果：fold-0 rank2 a29-TSB nested-adjusted +0.157（CI [0.095,0.212]），旧 smoke 的"被吃掉"（−0.004）是泄漏伪影；但独立患者 0bd3（+0.005）/286667（−0.002）均为零，且 a29 横跨两 fold 不能自我复制；TLS 最大效应 +0.034 未过项目 null 等价带（0.02）；rank1 全零是单维度残余空间坍缩的结构性零。

## 4. Unnamed-field closure（D-109）

`scripts/r04_unnamed_field_audit.py`（产物 `unnamed_field_audit_20260911.json`）：两 fold K=3 loading 的 SVD 内禀子空间 canonical correlation，基因行置换 null。rank-1 子空间 cc 0.582（null q975 0.089）→ SURVIVES；rank-2（0.903/0.474，最小 0.474 < 0.5 floor）→ NOT_IDENTIFIABLE；rank-3 只描述。共享的 rank-1 loading 方向真实存在，但无组成残余/外层支撑，不计为残余场。

## 5. Anchored synthesis 与最终判定（D-110）

4 个锚定残余候选（TLS/TSB × rank1/2）全部 DOES_NOT_SURVIVE；fold-4 TSB rank3 0.544 记 `ISOLATED_EXPLORATORY_SIGNAL`（rank1/2、fold-0、nested-adjusted 均不复现，不追第三方向）；uncertainty PASS；coordinate null 与 independent lineage 按规则记 NOT_TRIGGERED；R-05 记 `NOT_TRIGGERED_NO_SURVIVING_R04_FIELD`，R-06/R-07 不触发。对 H-01/H-03/H-05：当前数据不支持（非证伪）。

## 6. 条件工作处置（§6/§7 清单）

更强 deconvolution、fixed-cell-type state、coordinate-permutation null、independent lineage：全部 NOT_TRIGGERED（无 surviving candidate）。signed NB-GP 统一、NNLS 继续优化、USZ/common-panel 补强、新 K2 seed、K>3、新 external cohort、R-05/06/07 实质工作：记 deferred / not-triggered / limitation（见 ISSUES 处置与 roadmap 条目），不再悬挂为 open TODO。

## 7. 交付与审计线索

- 代码：`r04/composition.py`（fit/apply 接口）、`r04_finalize_k_semantics.py`（closure 模式）、`r04_k2_bridge_stage_a.py`（legacy 标记）、`r04_composition_final_audit.py`、`r04_unnamed_field_audit.py`、`r04_finalize_phase.py`（新）。
- 产物：`k_closure_canonical_20260911.json`、`composition_final_audit_20260911.json`、`unnamed_field_audit_20260911.json`、`r04_final_gate_20260911.json`（新）。
- 决策：D-107（K 精确语义＋D-092 取代）、D-108（nested 修复＋转录冲突）、D-109（审计规则预承诺）、D-110（最终判定＋合成规则）。
- 文档：roadmap（完成判据 outcome-neutral 化＋09-11 关闭条目＋R-05 触发态）、ISSUES（I-019 关闭；I-012/I-014 WONT_FIX；I-013 有限制关闭；I-017 保持开放备跨平台）、STATUS（重写）、TODO（K loop 关闭、条件工作记 NOT_TRIGGERED）、infra README。
- 测试：全量 313 项通过＋5 subtests（2026-09-11，约 10.5 分钟；上一轮基线 298 项）。新增覆盖：K closure 重算/防篡改/fold-4 交叉绑定、nested 泄漏捕获（新旧路径对照）、旋转/符号/置换不变性、门禁 fail-closed/阴性关闭/双患者复制/效应地板/未测试不可判。
