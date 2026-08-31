# 发布变更日志：R-01 增补 + R-02 验证谱系集成（2026-08-03）

本文件记录本次发布的审计结果、修复、再生与表面变化。当前操作契约见 `docs/PHASE_REPORT_R01_R02_2026-08-03.md`；科学动机不在此重述。

## 1. 发布范围

用户批准三条公开来源路径同时实施（D-037）并裁定人类标注者不区分资历、统一高置信度（D-036）后，完成：三个 validation 谱系的 R-01 增补注册与角色冻结、R-02 结构注册表的验证侧 GT 集成与门禁扩展、配套文档与清单更新。GSE274557 放弃索取；GSE274103 作者邮件由必需降为可选。

## 2. 包含的更新

- 数据获取：解包本地 `data/GEO/GSE175540/GSE175540_RAW.tar`（215 文件入 `raw/`）；下载 Zenodo 7760264 全部 15 个 zip（1.4 GB）并解压 14 个标注 CSV + 逐字提取 14 个坐标 CSV；下载 Zenodo 14620362 单 zip（2.1 GB）并选择性解压（高分辨率 tif 跳过），落盘 4.2 GB；三份 record/README JSON 快照存档。
- 决策：D-036（标注者资历规则）、D-037（三路径批准与角色治理）、D-038（实例定义统一 + USZ h5ad `obs/ground_truth` 标签列读取授权 + 失效条件）追加至 `docs/decisions.md`。

## 3. 审计结论（逐源）

- **GSE175540（KIRC）**：23 个标注文件两种表头（`TLS_2_cat` × 20、`TLS` × 3）；三个 `T_agg`-only 文件（GSM5924042/4045/4047）fail-closed 登记 NOT_AUDITABLE；GSM5924034 无沉积标注文件；`ffpe_c_21`/`frozen_b_13` 零 TLS 为已审计阴性（无行）。结果 18 个 AUDITABLE source、35 个 TLS 组件/18 患者。BioProject PRJNA732692 elink 无 BioSample，physical_specimen_id 采用 GSM。与 Atlas KIRC 行重复链接保持 study 级（GSM↔R_P crosswalk 未公开）；重算 35 组件 vs Atlas S4 的 30 条 R-ID 为定义粒度差异，以可重放重算为准。
- **TLS_VISIUM_USZ**：8 份 h5ad 的 `obs/ground_truth` 为 categorical；观测词汇 {TLS,INFL,TUM,NOR,UNASSIGNED,LN}（LN 33 个 spot，仅 LC3，超出记录描述的五类，TLS 范围外）；全部 spot in-tissue。结果 8 个 source、108 个 TLS 组件/8 患者。
- **ST_CRC_CMS**：14 个标注 CSV 四种表头变体（含 `Pathologist_KH`）；标注 spot 集与 in-tissue spot 集逐样本精确相等；39 项含拼写错误的词汇冻结进代码；仅 `tumor&stroma` 家族 4 标签入 TSB 范围，`IC aggregate*` 非 TLS fail-closed 排除；A798015 两个连续切片零 TSB（已审计阴性）。结果 12 个 source、551 个 TSB 组件/6 患者。
- **跨源重复检查**：Atlas S1/S2 全文扫描无 Valdeolivas/USZ 相关行；registry/repo 无 A416371/SN048 等样本重叠；USZ 与 Atlas LUAD 行（EGAS00001005021）确认为不同研究。

## 4. 代码与输出变化

- 新增 `scripts/r02_validation_gt.py`（三谱系记录构建 + 三个可执行 verifier + 冻结词汇）。
- `scripts/r02_build_registry.py`：`EXPECTED_TLS_COUNTS` 增 GSE175540=30、候选数 57→87、切分行 121→167；合并 validation 的 audits/instances/replay；policy 声明 4 个 verifier 与 `h5ad_obs_ground_truth_label_reads_allowed`；注册表 README 重写为四谱系表述。
- `scripts/r02_validate_gate.py`：候选计数与 lookup 扩展；verifier 调度泛化（`VERIFIER_FUNCTIONS`）；validation source/instance 独立重算比对（新错误码 `validation_gt_source_set_mismatch`、`validation_gt_source_mismatch`、`validation_instance_set_mismatch`、`validation_instance_mismatch`、`unknown_instance_source_class`）；replay 期望集合并 Heiser+validation（47→92）；conclusion_impact 重写。
- R-01 侧（本发布前段完成）：`r01_extract_validation_lineages.py`（新 extractor）、`r01_build_registry.py` 4 个 staging 合并、`r01_freeze_roles.py` ROLE_PLAN +3 行；逻辑单元 6→9、合格 physical rows 121→167。
- **修复**：`r02_validation_gt.py` 的 STCRC replay 行曾把相对路径传给 `Path.relative_to` 导致构建崩溃，已改为 root 绝对路径（首次实跑发现）。
- `tests/test_r02_structure_registry.py`：计数断言全部更新（候选 87、实例 1004、可审计 source 80、replay 92、切分 167/100/47、确认性 889=TLS 147+TSB 742）；新增 4 个 validation 篡改回归。`tests/test_r01_freeze_roles.py`、`tests/test_r01_extract_validation_lineages.py` 相应更新/新增。

## 5. 再生产物

`infra/sample-registry/`（physical_units 1002、role_freeze 9、source_assets 110、staging `validation_lineages_*.tsv` 4 表新增）与 `infra/structure-registry/` 全部输出（含 README）经两次确定性构建验证逐字节一致；`r01_gate.json`=`COMPLETE_WITH_EXCLUSIONS`、`r02_gate.json`=`PARTIAL_GT_READY`，门禁均为存储态重算一致。

## 6. 验证资产

- 全量测试：`python3 -m pytest tests/ -q` → 98 passed, 5 subtests passed（700.79 s）。
- 门禁计数（与预审计逐一相符）：候选 87（HTAN 44、KIRC 30、226997 4、274103 8、274557 1）；确认性 889（TLS 147、TSB 742，51 患者）；可审计 source 80（Heiser 42、KIRC 18、USZ 8、STCRC 12）；实例总数 1004；replay 92；切分 167 行/100 envelope/47 block group。
- 存储复查：2,982,261,653,504 bytes 可用（≥2 TB 约束满足，I-005 已更新）。

## 7. 规范表面变化

- 新增：`scripts/r02_validation_gt.py`；`infra/structure-registry/provenance/` 下三份 notes（Meylan/USZ/Valdeolivas）；`infra/sample-registry/staging/validation_lineages_*.tsv`；`data/GEO/GSE175540/raw/`、`data/other_sources/zenodo_st_crc_cms/`、`data/other_sources/zenodo_usz_tls_visium/`；`docs/PHASE_REPORT_R01_R02_2026-08-03.md` 与本文件。
- 更新：`docs/roadmap.md`（R-01/R-02 完成证据与假设影响）、`docs/ISSUES.md`（I-001 解除进度、I-005 存储复查）、`STATUS.md`（四段重写）、`data/other_sources/sources_manifest.tsv`（+3 行）、`data/DATA_BRIEFING.md`（GSE175540 解包行 + 两个 Zenodo 行）。
- 未删除任何文档：五份维护文档均仍为现行版本，无重叠摘要需要归档。

## 8. 最终发布状态

R-01 `COMPLETE_WITH_EXCLUSIONS`（9 单元、167 合格行）；R-02 `PARTIAL_GT_READY`（889 确认性实例、80 可审计 source、零完整性错误、零 blocker）；TLS 与 TUMOR_STROMA_BOUNDARY 对 HTAN/GSE175540/TLS_VISIUM_USZ/ST_CRC_CMS 四谱系冻结 claim-bearing；R-04 队列间一致性判据恢复可执行。

## 9. 发布后已知事项

- GSE226997/GSE274103/GSE211956 保持 GT-less；其任何确认性使用会触发角色治理违规（跟踪 I-001）。
- 血管与坏死无公开 GT；GSE274103 索取邮件仅剩该独占价值，是否发送由用户决定。
- STOmics 19 个 TLS 候选文件仍在冻结核心外，provenance 闭合未再尝试（三路径已补足，优先级取消）。
- `tmp/kirc_tls_audit.py` 为一次性审计脚本，属 scratch 非规范表面；`docs/`、`data/` 下 `.ipynb_checkpoints/` 为 Jupyter 残留，与本项目产物无关，未清理（无风险，留待用户决定）。
- USZ 门禁每次运行对 8 份 h5ad（约 2 GB）重算 sha256，单次构建+门禁约需分钟级，全量测试约 12 分钟，属预期。
