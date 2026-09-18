# 013_spatial｜R-16 恢复推进施工包

审阅日期：2026-09-18。基准提交：`73996528634ffe0362f1a0341771658171f3b5e2`。

**本轮建议：保持“连续空间场及其可逆信息”的目标，不再扩充聚类算法；局部修正判读错误，把已有候选推进到独立验证与信息增量测试。**

这不是新的全项目版本，也不是要求重开 R-04 或 K 搜索。建议将解压后目录内的文件放入仓库 `docs/recovery_20260918/`，把 `00_AGENT_MASTER_PROMPT.md` 交给执行 agent。

## 阅读顺序

| 读者/用途 | 文件 |
|---|---|
| 用户先看：结论、问题、下一步 | `01_评估与解决方案.md` |
| Agent 直接启动 | `00_AGENT_MASTER_PROMPT.md` |
| 并行分工、输入输出、完成条件 | `02_施工任务与依赖.md`、`config/tasks.json` |
| 统计量、组成控制、外部验证、遮蔽实验 | `03_方法与证据契约.md` |
| 修订项目文档，防止再次漂移 | `04_文档变更与决策.md` |
| 结果卡、阶段判断、论文贡献边界 | `05_交付与阶段判断.md` |
| 可调整的运行默认值 | `config/recovery.json` |
| 具体代码问题与定位 | `evidence/审计证据表.md`、`evidence/sources.tsv` |

## 已做与未做

已完成仓库最新提交、关键代码、阶段报告和决策的审阅。已在独立环境运行本包的 **12 项微型测试**，验证符号匹配、统计量聚合、嵌套等高线等数学/软件问题。

**没有取得服务器上的原始切片或缓存，没有重跑真实患者数据，没有重训模型，没有运行仓库全量测试，也没有修改远端仓库。** `evidence/toy_audit.json` 只含人工合成示例，不能用作生物学结果。

## Agent 可用的第一条命令

在仓库根目录运行（包已放入上述位置）：

```bash
python docs/recovery_20260918/tools/inspect_local_artifacts.py \
  --repo . --out infra/r16/recovery_20260918/preflight.json
```

这个命令只检查提交、文件和 registry 元数据，不加载表达矩阵、不删除文件、不启动计算。真实分析由 agent 按施工文档接入现有代码实现。

本地微型测试复核：

```bash
python docs/recovery_20260918/diagnostics/audit_reproductions.py \
  --out infra/r16/recovery_20260918/toy_audit.json
python -m pytest -q docs/recovery_20260918/diagnostics
```

统计检验失效只阻止相应显著性/否证结论，不阻止其他合法的描述、候选整理和建模探索。新数据、付费资源、GPU 及原有 hold 项仍遵守仓库授权边界。
