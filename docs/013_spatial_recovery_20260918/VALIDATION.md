# 本包验证记录

日期：2026-09-18。范围：交付文件、人工示例及本包只读工具。

## 实际执行

- `python diagnostics/audit_reproductions.py --out evidence/toy_audit.json`：通过，7 组人工示例。
- `python -m pytest -q diagnostics`：**12 passed**（最后一次报告 0.28s；不代表真实项目运行时间）。
- `python -m compileall -q diagnostics tools`：通过。
- `inspect_local_artifacts.py --help`：通过。
- 临时合成目录上的 inventory CLI smoke：通过，2 行 registry 正确统计。
- 同一输出文件二次运行不带 `--overwrite`：按预期拒绝，原文件保持不变。
- JSON 配置可解析、任务 ID/依赖可解析；ZIP 与 SHA256 清单检查通过。

## 环境（只是测试环境，不是目标服务器要求）

Python 3.13.5；NumPy 2.3.5；SciPy 1.17.0；pytest 9.0.2。

## 未执行与限制

没有运行服务器上的原始表达/缓存数据、Lane A/B 真实重分析、PC loading 真实重匹配、外部队列验证、遮蔽实验或仓库全量回归测试。没有修改远端仓库。12 项测试不是整个 013_spatial 已经被修复的证据。

源代码与结果记录审阅基于提交 `73996528634ffe0362f1a0341771658171f3b5e2`。数据深度归一化的上游语义列为待核实，不在本包中假装已查明。

## 文件完整性

`MANIFEST.sha256` 覆盖包内除它自身以外的全部交付文件。解压后在目录内可用 `sha256sum -c MANIFEST.sha256` 核验。
