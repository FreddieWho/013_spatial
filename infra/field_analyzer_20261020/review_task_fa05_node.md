# FA-05 逐节点复核（xai/grok-4.7，xhigh）

只审 FA-05：**缺口可恢复范围**。注意：你上一轮已发现 G1（汇总口径）与 G2（措辞），本轮确认修复并检查**其余**实现，不要重复要求已修的两点。

## 范围

- 代码：`field_analyzer/block_recovery.py`、`field_analyzer/cli.py` 的 `cmd_block_recovery` / `_unit_level` / `_skill_units`
- 产物：`infra/field_analyzer_20261020/block_recovery/{recovery_metrics.tsv,summary.json}`
- 报告：`infra/field_analyzer_20261020/REPORT_FA05.md`
- 历史：`review/fa01_06_review.json`（G1）、`review/g1_fix_review.json`（G1 关闭 + G2）

## 必答

1. **G1/G2 是否落实**：主口径是否单位等权（先 unit 内中位，再 unit 间中位）？报告/STATUS/roadmap 是否还有"8 pitch 可恢复"或"多场从未超过均值"的残留？
2. **遮挡是否真的遮住全部表达**：`evaluate_section` 里 hidden 与 visible 的构造；core Q 是否与半径无关；区域内任何派生量（深度、检出、QC）是否可能进入预测？
3. **多场岭回归是否公平**：特征只用"其他场"的平滑值？拟合只用可见行？有没有偷看目标场自身？`alpha=1` 是否固定且与单场参照同预算？
4. **参照是否可比**：背景均值用的是可见区均值（不是含隐藏区的全片均值）？IDW/核的 k 与带宽是否对四个场一致？
5. **中心选取**：`geometric_centers` 用坐标哈希、且要求最小间隔；是否真的与表达无关？3 个中心会不会系统避开组织边缘从而高估性能？
6. **数字一致性**：从 `recovery_metrics.tsv` 复算单位等权的 kernel 三档（应约 +0.129/−0.080/−0.204）与 multifield（+0.137/−0.004/−0.082），并核对 `n_units=8`、行数 2500。
7. **结论强度**："可恢复范围约 4 pitch"是否被产物支持？有没有把"在这批材料/这个观测量级下"的限定漏掉？

## 交付

写 `infra/field_analyzer_20261020/review/fa05_node_review.json`，`scope.read`/`scope.ran` 必须非空。字段同上模板（`node` 填 `FA-05`，findings 前缀用 `K`）。

不要修改仓库里除该 JSON 以外的任何文件。
