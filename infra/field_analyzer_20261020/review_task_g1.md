# G1 修复复核（xai/grok-4.7，xhigh）

你上一轮审 `9465644` 给出 `allow_phase_conclusion: true` 与 1 项 MAJOR：

> G1（FA-05）8 pitch 仍优于均值只在把中心和场当成独立行时成立；按单位等权 kernel 在 r=8 为 −0.08。

## 已做的修改（提交 `6b26f28`）

- `field_analyzer/cli.py`：新增 `_unit_level` / `_skill_units`——**先在每个单位内取中位 RMSE，再在 8 个单位间取中位**，作为主口径；原行级口径改名 `skill_vs_background_rows_secondary`。
- `summary.json`：新增 `n_units`、`skill_vs_background_units_primary`，notes 首条声明"主口径按单位等权，行不是独立观测"。
- 重跑 `block-recovery`，主口径结果：kernel +0.129 / −0.080 / −0.204；idw +0.156 / −0.116 / −0.270；kernel_bw2 +0.159 / −0.038 / −0.256；multifield_ridge +0.137 / −0.004 / −0.082。
- `REPORT_FA05.md`：结论改为"可恢复范围约 4 pitch"，并写明行级口径把 8 pitch 虚高到 +0.08 的原因；该修正归因于 G1。
- `REPORT_FA06.md`、`docs/roadmap.md`、`TODO.md`、`STATUS.md` 同步为 4 pitch 口径。
- `recovery_metrics.tsv` 实际 2500 行，报告原写 336 行，已改。

## 请回答

1. 新的主口径是否真的做到"单位等权"（检查 `_unit_level` / `_skill_units` 的实现与调用）？
2. 主口径数值是否与 `recovery_metrics.tsv` 可复算一致？行级口径是否仍被清楚标为次要？
3. 报告与状态文档里是否还残留"8 pitch 可恢复"或等价的过强表述？
4. `REPORT_FA06` 的"AF-03 负结果"在 4 pitch 口径下是否仍然成立、是否被写得更强或更弱了？

## 交付

写 `infra/field_analyzer_20261020/review/g1_fix_review.json`：

```json
{
  "reviewer": {"provider": "xai", "model": "grok-4.7", "thinking": "xhigh", "cwd": "/home/huyudi/013_spatial"},
  "commit_reviewed": "<git rev-parse --short HEAD>",
  "g1_closed": true,
  "items": [{"question": 1, "verdict": "...", "evidence": "..."}],
  "new_findings": [],
  "allow_phase_conclusion": true,
  "limitations": "..."
}
```

不要修改仓库里除该 JSON 以外的文件；只报实际读到/跑到的内容。
