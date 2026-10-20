# FA-00 合同预审任务（交指定审核代理 xai/grok-4.7，high）

## 你的角色与边界

你是本项目的独立审核代理，不是实现者。本轮任务是**只读预审**：判定新采用的观测/掩膜合同与资产复用裁定是否足以支撑 FA-00 的实现，找出会改变科学解释的缺口。

硬边界：

- 不修改 `docs/plan.md`、`docs/roadmap.md`、`STATUS.md`、`TODO.md`、`docs/decisions.md`、`docs/ISSUES.md` 的任何内容；不改任何代码。
- 不读取人工 TLS 标签文件、不读取旧桥端点/标签派生产物、不运行新模型或真实数据拟合。
- 不下载任何外部数据、权重或依赖；不使用 GPU。
- 本轮结论只覆盖文档与接口规格，不得表述为"合同已被证明无缺陷"。

## 先读这些（顺序固定）

1. `AGENTS.md`（执行/审核分工与资源约束）
2. `docs/plan.md`（本轮 AF-01～AF-04 主线与第 8 节三条历史解释纠正）
3. `docs/roadmap.md` 顶部"当前执行口径"与"新节点状态"（FA-00～FA-07）
4. `docs/field_analyzer_killtest_20261020/05_COMMON_PROTOCOL.md`（**主要审核对象**）
5. `docs/field_analyzer_killtest_20261020/01_REPO_AUDIT_AND_REUSE.md` 与 `reuse_register.tsv`（复用裁定）
6. `docs/field_analyzer_killtest_20261020/07_ALGORITHM_ROUTES.md` §0 与 §2.3（共同输出、B 的联合求解能量式）
7. `docs/field_analyzer_killtest_20261020/08_REAL_DATA_AND_SCALE.md` §1 与 §6（真实资产边界）
8. `docs/ISSUES.md` I-053／I-054／I-055
9. 需要时对照源码：`scripts/field_discovery_20261018.py`（`read_matrix`/`feature_views`/`neighbor_graphs`）、`scripts/field_candidate_subfield_20261019.py`（`split_genes`/`cache_section`）、`scripts/field_biological_review_20261019.py`（`describe_field`/`band_autocorr`）。只看这些函数是否支持你的判断，不做全库审计。

## 必须回答的问题

按下列编号逐条给判定（PASS / GAP / BLOCKER + 一句理由）。GAP 必须写清"会影响什么解释"；BLOCKER 必须是会导致结果不能用于科学结论的问题。

1. **掩膜入口完整性**：05 §3 的"掩膜必须作用于整个信息入口"是否覆盖了 FA-00 计划实现的全部入口（读取、归一化、可测性/资格筛选、图特征、预平滑、PCA、查询点选择）？列出仍可能让隐藏值回流的路径。
2. **B 核验组隔离**：`exclude_all_witness_genes_from_all_input_layers` 与"跨多个程序须排除所有核验 B 的并集"是否在接口层可强制？旧 `cache_section` 同时保存 A/B 与全片统计是否会破坏该约束？
3. **坐标与轴语义**：05 §2 的 `ObservationBundle` 字段是否足以强制显式 barcode/轴/单位映射（禁止形状启发式转置）？缺口在哪？
4. **单位与分组**：同单位多片、同 base_id 的 mask/旋转重复不跨外层角色的要求在 FA-00 可核查的最小证据是什么？
5. **比较公平性**：05 §6 强参照清单与 07 §2 的 B 基线，是否足以避免"候选被优待"或"基线被刻意削弱"？多输出 LMC/GP 的最低配置（≥2 长度尺度 + 输出专属噪声）是否写死到可执行？
6. **毒化反例可执行性**：05 §3 的"替换/置换隐藏区与 B 后再重建全部输入，指纹与预测不变"在小样本上能否一次跑通？需要哪些最小断言？
7. **真实数据不可执行的降级路径**：08 §1 与 I-053 的降级顺序是否完备？"资产缺失"与"资产存在但不满足足量"是否被区分？
8. **B 的联合求解能量式**：05/07 给出的精度构造 `Q = τI + Σ LᵀWL` 与右端项 `MᵀR⁻¹y_O + τμ0 + ΣLᵀWb` 是否自洽、是否遗漏会让"已学背景被拉向零"的项？指出现式错误。
9. **FA-00 范围**：以 04 roadmap 的 FA-00 交付物（`reuse_resolution.tsv` + 信息入口反例 + 小几何样本）为准，判断当前范围是**太宽**、**合适**还是**太窄到无法支撑 FA-01/FA-03**，并给出最小增删。
10. **不可越界项**：是否有任何一条现有文字会把旧 TLS 端点、旧未校准 p 值、或"猜中即发现"的语义带回新分支？

## 交付

写到 `infra/field_analyzer_20261020/review/preflight_fa00.json`，严格 JSON：

```json
{
  "reviewer": {"provider": "xai", "model": "grok-4.7", "thinking": "high", "cwd": "/home/huyudi/013_spatial"},
  "scope": {"read": ["..."], "not_checked": ["..."]},
  "commit_reviewed": "<git rev-parse --short HEAD>",
  "items": [{"id": 1, "verdict": "PASS|GAP|BLOCKER", "reason": "...", "required_fix": null}],
  "allow_fa00_implementation": true,
  "blockers": [],
  "required_fixes_before_run": [],
  "recommended_minimal_scope_change": "...",
  "limitations": "..."
}
```

`allow_fa00_implementation` 为 true 的条件：没有 BLOCKER；GAP 可在实现中通过固定断言关闭。若存在会改变解释的泄漏/分组/坐标问题，置 false 并写清解除条件。

## 报告中必须自律

- 不承诺"绝对无缺陷"，只报你实际读到和实际推理出的内容。
- 明确区分"文档缺口"与"实现可能出错"；后者若无法只看文档判断，写 `NOT_ASSESSABLE_FROM_DOCS`。
- 若认为某项要求本身过度（例如要求先修好旧 p 值才开工），直接说它不应成为前置条件。
