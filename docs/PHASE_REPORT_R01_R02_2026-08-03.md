# 阶段报告：R-01 + R-02 checkpoint（2026-08-03）

本文件是 R-01（样本注册与角色冻结）与 R-02（结构 GT 与切分冻结）两个基础设施节点的交接文档，面向接手研究的读者：只靠本文件与所引用的注册表，即可理解范围、产物、接口、流程与当前冻结表面。本文是时间点快照，不是持续维护文档；项目的五份维护文档见 `AGENTS.md`（`STATUS.md`、`docs/plan.md`、`docs/roadmap.md`、`docs/decisions.md`、`docs/ISSUES.md`）。

## 1. 范围与状态

- **科学背景**：项目检验"目标结构未被当前切片直接看见时，周围分子生态能否在排除捷径后给出可信位置概率"（假设 H-01…H-07 见 `docs/plan.md`）。R-01/R-02 是基础设施，不直接检验任何假设，但决定后续哪些数据与结构能承担确认性证据。
- **R-01 状态**：`COMPLETE_WITH_EXCLUSIONS`（`infra/sample-registry/r01_gate.json`）。1002 条注册记录中 167 条合格 physical rows，折叠为 9 个冻结逻辑单元；零身份完整性错误。
- **R-02 状态**：`PARTIAL_GT_READY`（`infra/structure-registry/r02_gate.json`）。889 个确认性结构实例、80 个可审计 GT source、87 条 scoped TLS 候选；零完整性错误；TLS 与 TUMOR_STROMA_BOUNDARY 对 4 个谱系冻结为 claim-bearing。

## 2. 本阶段产出什么

两个 fail-closed 控制面：

1. **样本注册表**（`infra/sample-registry/`）：回答"每条物理记录是谁、承担什么角色"。跨 HEST/GEO/10x/HTAN/STOmicsDB 等来源去重，恢复 patient/block/section 身份（或在 E2 互证下登记 patient-linked physical specimen 例外），并按 outcome-blind 规则冻结角色。
2. **结构注册表**（`infra/structure-registry/`）：回答"哪些结构实例可以承载科学主张、在什么输入与切分规则下"。登记结构本体、GT 来源审计、实例级几何、外层切分与输入政策；机器门禁独立重算，不信任任何自报字段。

## 3. 冻结角色（R-01 v1）

| 逻辑单元 | 角色 | 合格 physical rows | GT 状态 |
|---|---|---|---|
| HTAN_VANDERBILT_CRC | training | 47（piece 级 h5ad） | 可审计 GT：TLS 4 + TSB 191 |
| GEO::GSE175540 | external_validation | 24（一人一样本） | 可审计 GT：TLS 35 |
| TLS_VISIUM_USZ | external_validation | 8（肾 3 + 肺 5） | 可审计 GT：TLS 108 |
| ST_CRC_CMS | internal_validation | 14（7 患者 × 2 连续切片） | 可审计 GT：TSB 551 |
| GEO::GSE226997 | internal_validation | 4（specimen-equivalent） | GT-less，不承担确认性主张 |
| GEO::GSE274103 | internal_validation | 5（specimen-equivalent） | GT-less（作者标注未沉积，索取可选） |
| GEO::GSE211956 | external_validation | 8（specimen-equivalent） | GT-less（QuPath 几何未沉积） |
| GEO::GSE274557 | discovery | 55（specimen-equivalent） | 不在确认路径上，已放弃索取 |
| TENX_V1_BREAST_CANCER_BLOCK_A | serial-section validation | 2（同一真实 block） | 结构 GT 未审计（角色专用） |

角色分配只使用身份、谱系、平台与预注册 capability（outcome-blind）；冻结后不得以 TLS 阳性率、结构数量或模型表现更换单元。

## 4. 规范输入（上游）

- `paper/tables/science.adz2742_tables_s1_to_s8.xlsx` — Atlas Table S1–S8（TLS 候选与计数的 source-of-truth）。
- `repo/` — Atlas 分析代码与 `data_meta/`（含 `ST_CRC_cohort_meta2.csv`，Heiser piece 级 crosswalk 证据资产之一）。
- `data/GEO/` — 19 个 GSE supplementary；**GSE175540 的 RAW.tar 已解包**（`data/GEO/GSE175540/raw/`，215 文件：24 样本 Space Ranger 输出 + 23 个作者 TLS 标注 CSV），其余 18 个 RAW.tar 未解包。
- `data/other_sources/htan/spatial_CRC_atlas_repo/` — Heiser GitHub 仓库，commit `64e585453514801a3b730f86aec90bbc0f0595da` 钉死（42 个逐 spot 病理标注 CSV + visium_sample_key）。
- HTAN OSF h5ad ×49（`data/other_sources/htan/`，osf.io/hftq2，piece 级坐标重放资产）。
- `data/other_sources/zenodo_st_crc_cms/` — Zenodo 7760264：14 个病理医生标注 CSV（`Pathology_SpotAnnotations/`）+ 14 个坐标 CSV（`tissue_positions/`，从样本 zip 逐字提取）+ record JSON + GitHub README 快照。
- `data/other_sources/zenodo_usz_tls_visium/` — Zenodo 14620362：`TLS_VISIUM_USZ/h5ad_preprocessed/{KC1-3,LC1-5}.h5ad`（zip 选择性解包，高分辨率 tif 未取）+ record JSON。
- `infra/bioinf-data-index/` — 本地 metadata 索引（GEO soft 等审计证据）。
- `data/other_sources/10x_genomics/` — 10x Breast Block A 两个 section（serial-section 角色）。

## 5. 规范输出（下游契约表面）

### infra/sample-registry/

| 文件 | 行数 | 内容 |
|---|---|---|
| `physical_units.tsv` | 1002 | 全部注册物理记录；合格行 `record_status=RESOLVED_INCLUDED_CANDIDATE`（167 条） |
| `role_freeze.tsv` | 9 | 逻辑单元角色冻结（`record_status=FROZEN`） |
| `source_assets.tsv` | 110 | metadata 资产内容寻址（sha256 + 校验状态） |
| `identity_evidence.tsv` | 8168 | 身份证据（E2/E3 分级，raw↔normalized 闭合） |
| `duplicate_groups.tsv` | 234 | 重复/保守 leakage 分组 |
| `r01_gate.json` | — | 机器门禁：`COMPLETE_WITH_EXCLUSIONS` |
| `staging/` | — | 各 extractor 的可重建中间表（含 `validation_lineages_*.tsv` 4 表） |

### infra/structure-registry/

| 文件 | 行数 | 内容 |
|---|---|---|
| `structure_ontology.tsv` | 4 | TLS / BLOOD_VESSEL / NECROSIS / TUMOR_STROMA_BOUNDARY；前两者 TLS、TSB 为 `FROZEN_CLAIM_BEARING`，血管/坏死 `NOT_FROZEN_NO_SCOPED_GT` |
| `tls_candidate_summary.tsv` | 87 | source-reported TLS ID 候选（非确认 inventory 行） |
| `gt_source_audit.tsv` | 107 | GT 来源审计；80 行 `AUDITABLE_GT` |
| `structure_instances.tsv` | 1004 | 结构实例；889 行 `CONFIRMATORY`（TLS 147 + TSB 742，51 患者） |
| `h5ad_replay_index.tsv` | 92 | 坐标重放索引（spot 计数 + 指纹钉死） |
| `outer_splits.tsv` | 167 | 外层切分；100 个 patient-wide envelope + 47 个真实 block group |
| `cross_section_links.tsv` | 1 | 10x Block A 两 section 的 identity-only 链接 |
| `input_policy.tsv` | 7 | 输入通道政策（见 §8） |
| `leakage_selection_audit.tsv` | 3 | outcome/预处理选择偏差审计 |
| `r02_policy.json` / `r02_gate.json` | — | 政策声明与机器门禁（`PARTIAL_GT_READY`） |
| `README.md` | — | 注册表自述（每次构建重新生成） |
| `provenance/` | 5 | 四份 provenance notes（Heiser/Meylan/USZ/Valdeolivas）+ Heiser 全文 XML 快照 |

**Join keys**：`physical_unit_id`（两注册表主键）、`patient_id` / `block_id` / `logical_unit_id` / `leakage_group_id`（切分与角色）、`gt_source_id`（审计↔实例）、`instance_id`、`identity_envelope_id` / `block_group_id`（外层组）。

## 6. 处理流程

### R-01（按依赖顺序）

1. `r01_extract_atlas.py` / `r01_extract_hest.py` / `r01_extract_htan.py` / `r01_extract_tenx_explicit.py` / `r01_extract_external_geo.py` / `r01_extract_validation_lineages.py` — 各来源 staging 提取（输出均在 `staging/`，可重建）。
2. `r01_inventory_metadata.py` / `r01_prepare_metadata_request.py` / `r01_summarize_units.py` — metadata 盘点与获批外部 metadata 请求。
3. `r01_cross_source_duplicates.py` — 跨来源重复与保守 leakage 分组。
4. `r01_build_registry.py` — 合并 staging 为注册表。
5. `r01_freeze_roles.py` — outcome-blind 角色冻结（ROLE_PLAN 为唯一角色来源）。
6. `r01_validate_gate.py` — 独立复核，产出 `r01_gate.json`。

### R-02

1. `scripts/r02_build_registry.py` — 从 Atlas xlsx 重建 TLS 候选（87 条，与 Table S2 计数互核）；调用 `r02_heiser_gt.py::build_heiser_records` 与 `r02_validation_gt.py::build_validation_records` 重算 GT source/实例/replay；构建外层切分（保守患者 envelope，D-031）、cross-section 链接、输入政策；写出全部注册表并内嵌调用门禁。
2. `scripts/r02_validate_gate.py` — 独立门禁：sha256 复核全部 GT 资产、可执行 verifier 逐 source 重算、Heiser/validation 两套记录与注册表逐字段比对、Table S2/S4 计数互核、STOmics 候选清点（19 个，core 外）、切分/泄漏完整性。任何不一致即 `HARD_BLOCKED_*`。
3. `scripts/r02_heiser_gt.py` / `scripts/r02_validation_gt.py` — 共享确定性 GT 逻辑（六边形网格连通组件、坐标重放、verifier）；构建与门禁调用同一实现。

### 测试

`tests/`：15 个测试文件共 98 项测试（`python3 -m pytest tests/ -q`，约 12 分钟；USZ h5ad 校验和使单次构建较慢）。`test_r02_structure_registry.py` 含篡改回归：伪造/删除/改写 source、实例、几何、校验和、replay、切分与链接均须被门禁拒绝。

## 7. 关键参数（改变行为者）

- `EXPECTED_TLS_COUNTS`（build 与 gate 双处一致）：HTAN 44、GSE175540 30、GSE226997 4、GSE274103 8、GSE274557 1；候选总数 87。
- 可执行 verifier 四类：`SPOT_BARCODE_PATHOLOGY_ANNOTATION_CSV`（Heiser，44 行）、`SPOT_BARCODE_TLS_ANNOTATION_CSV`（KIRC，21 行）、`H5AD_OBS_GROUND_TRUTH_LABELS`（USZ，8 行）、`SPOT_BARCODE_PATHOLOGY_CATEGORY_CSV`（STCRC，12 行）；与 `r02_policy.json` 声明不一致即门禁失败。
- 实例定义：六边形网格连通组件（邻居 (0,±2)、(±1,±1)，in-tissue 过滤），四谱系统一。
- 标签词汇（fail-closed，未知值即构建失败）：KIRC `TLS_2_cat`{TLS,NO_TLS,空}/`TLS`{T_agg,空}；USZ {TLS,INFL,TUM,NOR,UNASSIGNED,LN}；STCRC 39 项冻结词汇，仅 `tumor&stroma` 家族 4 标签入 TSB 范围。
- 排除规则：T_agg ≠ TLS（3 个 KIRC 文件 NOT_AUDITABLE）；`IC aggregate*` ≠ TLS；空标注单元 = unknown 绝不作阴性；零目标标签样本为已审计阴性（无审计行）。
- h5ad 授权读取列：D-033（`obs/_index`、`array_row`、`array_col`）+ D-038（仅 USZ 八份 h5ad 的 `obs/ground_truth`）；表达值（X/layers/obsm）与图像像素一律禁止。
- 外层切分：真实 block 优先；block 未知的 specimen-equivalent 只冻结患者级 envelope（`PATIENT_ENVELOPE_BLOCK_UNKNOWN`），不计 block-level 证据（D-031）。

## 8. 下游使用契约（R-04 起）

允许：
- 使用 `RESOLVED_INCLUDED_CANDIDATE` 的 167 条 physical rows 及其 `outer_splits.tsv` 冻结外层折。
- 对 4 个有 GT 的谱系（HTAN、GSE175540、TLS_VISIUM_USZ、ST_CRC_CMS）的 TLS 与 TSB 做确认性检验；R-04 队列间一致性、R-05 两来源一致性判据限此范围执行。
- 分子-only 输入（表达矩阵及其合法派生）在训练折内使用。

禁止（`input_policy.tsv` 机器强制）：
- 表达矩阵/派生表达、图像像素/形态派生进入 R-02 审计面；GT 定义标注及其传递闭包作为模型输入（`FORBIDDEN_AS_MODEL_INPUT`）。
- outcome/response metadata 用于选择、切分或建模；identity/site/treatment/stage/sample/file metadata 作为预测输入。
- GSE226997/GSE274103/GSE211956/GSE274557 承担任何确认性主张（GT-less 或 discovery）。
- 血管/坏死的确认性主张（无公开 GT）；STOmics 19 个候选文件作为 GT；Atlas 汇总计数/ID 作为几何 GT。

## 9. 当前操作约束

- metadata-only + 已批准的标注内容范围（D-032/D-037）；新增外部数据须先说明证据缺口、不可替代性、规模与解锁判据并经用户批准。
- 禁 GPU；CPU/内存可用；存储余量 ≥2 TB（2026-08-03 实测 2.98 TB）。
- 构建与门禁不访问网络；下载须 `unset` 全部代理变量。
- 结构与样本依据只来自数据自带 metadata；文件名不得静默补成标签；来源不明 fail closed。
- Agent 不做任何 git mutation。
- 系统 python3（/opt/anaconda3）带有 anndata/h5py/openpyxl；`.venv_hest` 无 anndata，勿用于运行注册表构建。

## 10. 阅读顺序

1. `STATUS.md` — 两分钟通俗现状。
2. `docs/plan.md` — 科学问题与判定标准。
3. `docs/roadmap.md` — 节点图、R-01/R-02 完成证据、R-04 起的判据。
4. 本文件 — 接口与产物细节。
5. `infra/sample-registry/README.md`、`infra/structure-registry/README.md` — 注册表自述。
6. `infra/structure-registry/provenance/` — 四个 GT 来源的逐源审计依据。
7. `docs/decisions.md`（D-030…D-038）— 技术选择与失效条件；`docs/ISSUES.md`（I-001 等 OPEN 项）。
8. 审计/修复/再生历史：`docs/RELEASE_CHANGELOG_2026-08-03.md`。
