# 修复后任务资格与方向性描述（2026-09-22）

决策：D-143，报告：`docs/DIRECTIONAL_ELIGIBILITY_20260922.md`。本地 CPU，无新增数据。

复现（输出目录必须不存在；现有产物先保留到另一个目录）：

```bash
LD_LIBRARY_PATH=/opt/anaconda3/lib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/spatial-repair-mpl PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=. python scripts/r16_directional_eligibility.py
```

`program_eligibility.tsv` 区分源覆盖与分析面板覆盖；`source_annotations.tsv` 保留原标签和缺失类别；`registry_scope.tsv` 标出本轮重放范围（不是重新审计所有来源）；`focus_eligibility.tsv` 列出所有 108 焦点及剔除原因；`focus_effects.tsv` 与 `patient_effects.tsv` 分开重复焦点与患者；`receipt.json` 保存输入/输出 hash 与未运行项。

方向差异是已知 TLS/TUM 锚点上的片内残差描述；同径向带不等于逐点精确距离匹配，不消除所有区室/形态混杂。无空间 p、无隐藏预测、无跨队列检验；不用于判定定位成立或否证所有空间场。两次前置失败及恢复见 execution_notes.json。
