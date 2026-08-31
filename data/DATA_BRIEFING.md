# 肿瘤空间转录组数据简报

生成日期：2026-07-27 ｜ 数据根目录：`/home/huyudi/013_spatial/data/` ｜ 总量约 **1.9 TB**

## 一、论文配套数据（Pan-Cancer Spatial Atlas of TLS, Science 2026, doi:10.1126/science.adz2742）

依据原文 Table S1/S2（`paper/tables/science.adz2742_tables_s1_to_s8.xlsx`）逐项核对收集，分析代码仓库在 `repo/`。

| 来源 | 内容 | 状态 | 大小 |
|---|---|---|---|
| GEO（19 个数据集） | 论文全部 GEO 队列的 supplementary files（Visium 等 processed 数据，RAW.tar 未解包） | ✅ 19/19 到齐，大小与官方清单核对一致 | 69 G |
| 10x Genomics | 论文使用的 4 个样本：B37（乳腺 FF）、CO1（结直肠 FFPE）、Ov1（卵巢 FFPE）、L6（肺神经内分泌癌 FFPE，CytAssist 11mm） | ✅ 已核对到具体官方数据集 | 含于下表 |
| Dryad | doi:10.5061/dryad.h70rxwdmj（GBM，样本 Br21；28 个 UKF Visium 样本的 10XVisium_2.zip） | ✅ 7.8 G（仅 Visium 部分，IMC/MALDI 等未取） | 7.8 G |
| HTAN / OSF | HTAN Vanderbilt CRC 队列 49 个处理后 h5ad（osf.io/hftq2，与 repo 元数据样本号一一对应） | ✅ | 11 G |
| HTAN / GitHub | Heiser et al. 2023 分析仓库 `Ken-Lau-Lab/spatial_CRC_atlas`（commit 64e5854 钉死）：42 个逐 spot 病理标注 CSV + visium_sample_key，为 R-02 GT 来源（D-032 范围） | ✅ 2026-07-31 收录 | 96 M |
| GEO GSE175540 | Meylan et al. 2022（Immunity）ccRCC，论文 Table S2 KIRC 队列；RAW.tar 已解包：24 样本 Space Ranger 输出 + 23 个作者沉积逐 spot TLS 标注 CSV + family soft 快照 | ✅ 2026-08-03 解包（D-037 路径 1，R-02 external validation GT） | 1.1 G（raw/） |
| EGA / GSA / GSA-Human / NODE | EGAS00001005021、EGAD00001008031、PRJCA027506、HRA000437、HRA000979、OEP001756 | ⏸ 受控，待手动申请 | — |

## 二、补充肿瘤空转数据（超出论文范围）

| 数据集 | 规模 | 内容 | 大小 |
|---|---|---|---|
| HEST-1k 筛选子集 | 492 样本 / 5109 文件 | 人类肿瘤（362）+ 健康（130）组织，排除脑/脊髓；乳腺 140、肠 84、心 62、前列腺 61、肾 39、肺 31 等；Visium 244、ST 163、Xenium 61、Visium HD 24；含 h5ad、WSI、cellvit/tissue 分割、Xenium transcripts | ~1.0 T |
| 10x 近 3 年肿瘤全集 | 73 个数据集 | 2023 年后 10x 官网肿瘤空间数据：Visium/Visium HD/Xenium(XR/Prime)/CytAssist，覆盖乳腺、肺、卵巢、结直肠、前列腺、胰腺、肾、肝、皮肤、宫颈、GBM 等 | 592 G |
| STOmicsDB | 78 个肿瘤数据集 / 4356 文件 | 华大平台为主的肿瘤空转 processed 数据（乳腺、前列腺、结直肠、GBM、PDAC、HCC、卵巢、黑色素瘤、NSCLC 等） | 160 G |
| SpatialDB | 2 个数据集 | 前列腺癌（PMID 29925878）、黑色素瘤（PMID 30154148） 1.7 G |
| Zenodo 7760264（ST_CRC_CMS） | 7 患者 × 2 连续切片 / 15 zip | Valdeolivas et al. 2024（npj Precis Oncol）CRC Visium：14 个病理医生逐 spot 类别标注 CSV + 14 样本坐标文件；R-02 internal validation GT（肿瘤—基质边界，D-037 路径 2） | 1.4 G |
| Zenodo 14620362（USZ TLS Visium） | 8 样本（肾 3 + 肺 5） | USZ 专家人工 H&E 标注（TLS/Immune/Tumor/Normal），标签存于 h5ad `obs/ground_truth`；zip 选择性解包（高分辨率 tif 未取）；R-02 external validation GT（TLS，D-037 路径 3） | 4.2 G |

## 三、HOLD / 待手动事项

- **CROST**（1013 个肿瘤样本候选）：服务端文件接口持续宕机 + 现有数据量已足，存档 hold；恢复命令见 `other_sources/crost/HOLD.txt`。
- **Synapse/HTAN 原始数据**：用户手动处理。
- **受控数据申请**：EGA（DAC 审批 + pyega3）、GSA/GSA-Human（NGDC 在线申请）、NODE OEP001756（向所有者申请授权）——流程与入口 URL 汇总在 `access_control_needed/access_control_datasets.tsv`。

## 四、关键索引文件

- 总览：`data/README.md`
- 受控数据申请表：`data/access_control_needed/access_control_datasets.tsv`
- 非 GEO 来源明细：`data/other_sources/sources_manifest.tsv`
- 补充资源链接表：`data/other_sources/supplementary_tumor_ST_resources.tsv`
- HEST 样本清单：`data/other_sources/hest1k/selected_samples.tsv`
- 10x 清单：`data/other_sources/10x_genomics/manifest_3yr.tsv`
- STOmicsDB 清单：`data/other_sources/stomicsdb/manifest.tsv`

## 五、使用注意

- GEO 的 `*_RAW.tar` 除 GSE175540（R-02 GT 需要，D-037，已解包）外均未解包，按需 `tar -xf`（磁盘余量约 3.0 T）。
- Dryad 的 Br21 与 zip 内 UKF 样本编号的对应关系需查论文补充材料确认。
- HEST 数据为 CC-BY-NC-SA-4.0，仅限科研用途；HEST 下载环境在 `.venv_hest/`（token 方式见该目录脚本）。
- 所有后续大批量下载均已绕开本机付费代理（国内站直连 / hf-mirror）；勿在代理环境变量下直接跑下载脚本。
