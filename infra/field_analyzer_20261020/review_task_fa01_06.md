# FA-01～FA-06 节点完成复核（xai/grok-4.7，xhigh）

FA-00 已由你两轮复核（`fa00_completion.json`、`fa00_fix_review.json`）。本轮复核**其余六个节点**的已交付工作，重点只有一件事：**有没有会让结果解释变错的错误**。

## 复核范围（提交 9465644 及之前的 FA-01～06 交付）

| 节点 | 代码 | 产物 |
|---|---|---|
| FA-01 尺度/方向 | `field_analyzer/scales.py` | `infra/field_analyzer_20261020/scales/*.tsv`、`REPORT_FA01.md` |
| FA-02 合成 | `field_analyzer/synthetic.py` | `synth/scenario_manifest.tsv`、`pair_matching_diagnostics.tsv`、`REPORT_FA02.md` |
| FA-03A 路线 A | `field_analyzer/route_a.py` | `route_a/*.tsv`、`REPORT_FA03A.md` |
| FA-03B 路线 B | `field_analyzer/route_b.py` | `route_b/*.tsv`、`REPORT_FA03B.md` |
| FA-03C 路线 C | `field_analyzer/route_c.py` | `route_c/*.tsv`、`REPORT_FA03C.md` |
| FA-04 真实核验 | `field_analyzer/real_validation.py` | `real_validation/*.tsv`、`REPORT_FA04.md` |
| FA-05 缺口恢复 | `field_analyzer/block_recovery.py` | `block_recovery/*.tsv`、`REPORT_FA05.md` |
| FA-06 汇总 | 无新计算 | `REPORT_FA06.md` |

## 必答问题

1. **A/B 隔离是否真的对称**：`real_validation.half_specs` 的两套 spec，其外部文库是否都排除了 A∪B？有没有半边通过分母或 witness 进入另一半？
2. **FA-02 的配对是否仍可能被"匹配手段"污染**：改版后的 S2（纯位移 + 噪声标定）是否真的只改位移、不改同点相关以外的结构？S1 的径向匹配、S5 的可见逐位相同是否仍成立？
3. **模型选择的口径**：FA-03A 的"最优 CV 5% 内取最简"是否会掩盖真实差异（例如双峰种子 8）？报告是否已经把这种情形写成"不可识别"而不是"模型更好/更差"？
4. **FA-03B 的能力扫描与 S2 结论是否自相矛盾**：扫描 15/18 恢复 vs S2 全部判错，报告是否把差别归因于观测（并有潜在真值对照），而不是同时宣称"有能力"和"没能力"？
5. **FA-03C 的指标是否有鉴别力**：C0/C1 的对照是否用了同一端点与同一指标？"朴素 C1 更差"是否可能是实现选择（max over δ）导致，而非能力不存在？报告有没有把这一点写清楚？
6. **FA-05 的遮挡是否真的遮住全部表达**：`evaluate_section` 里 `visible` 是否同时排除隐藏点与组织外点；`multifield_ridge` 是否只在可见行上拟合；核心 Q 是否在不同半径下固定？
7. **统计单位**：所有汇总是否按单位/片聚合而不是按 spot 计数？有没有把重复中心当成独立样本？
8. **报告与产物一致性**：抽 5 个报告里的数字与 TSV/JSON 对齐（尤其 FA-04 的 0.776/0.578/0.79 与 FA-05 的 +0.245/−0.142）。
9. **结论是否越界**：REPORT_FA06 说"四条 AF 假设本轮均未获支持"，是否被所有产物支持？有没有哪条实际上获得了部分支持却被写成否定？

## 交付

写 `infra/field_analyzer_20261020/review/fa01_06_review.json`：

```json
{
  "reviewer": {"provider": "xai", "model": "grok-4.7", "thinking": "xhigh", "cwd": "/home/huyudi/013_spatial"},
  "commit_reviewed": "<git rev-parse --short HEAD>",
  "scope": {"read": [...], "ran": [...], "not_checked": [...]},
  "findings": [{"id": "G1", "node": "FA-0x", "severity": "BLOCKER|MAJOR|MINOR", "topic": "...", "evidence": "...", "required_fix": "..."}],
  "answers": {"1": "...", "2": "...", "...": "..."},
  "claims_that_do_not_hold": [],
  "allow_phase_conclusion": true,
  "limitations": "..."
}
```

`allow_phase_conclusion` 为 true 的条件：没有 BLOCKER；MAJOR 不改变"AF-01/02/03 未获支持"这一结论的范围。可读代码、跑测试与命令复算；不要修改仓库里除该 JSON 以外的任何文件。
