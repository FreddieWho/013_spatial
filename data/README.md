# Pan-Cancer Spatial Atlas of TLS — 数据下载总览

对应论文 "Pan-Cancer Spatial Atlas of Tertiary Lymphoid Structures"（分析代码：repo/，
来源 GitHub Coolgenome/Pan-Cancer_Spatial_Atlas_TLS）。

## 目录结构

- `GEO/` — 19 个 GEO 数据集的 supplementary files（均为 processed data，未含 SRA 原始 fastq）；
  多数为 `<GSE>_RAW.tar` 整包（内含各 GSM 的矩阵/图像文件）+ `filelist.txt` 清单，使用前按需 `tar -xf` 解包。
- `other_sources/`
  - `10x_genomics/` — 10x 官网公开肿瘤 Visium 数据集 10 个（BRCA×4、Colorectal、Ovarian×2、GBM、Cervical、Lung SCC）
  - `dryad/` — doi:10.5061/dryad.dv41ns28r（RCC 抗 PD-1 试验 TLS scRNA-seq，推断为论文 ref 118）
  - `htan/` — HTAN Vanderbilt (HTA11) CRC Visium 处理后 h5ad ×49（来自 OSF osf.io/hftq2）
  - `sources_manifest.tsv` — 非 GEO 来源明细 manifest
  - `supplementary_tumor_ST_resources.tsv` — 补充的公开肿瘤空转资源表（CROST、HEST-1k、SpatialDB、STOmicsDB、10x、HTAN、TCIA）
- `access_control_needed/access_control_datasets.tsv` — 需手动申请的受控数据（EGA/GSA-HRA/NODE 等）

## GEO 下载状态

| GSE | 大小 | 状态 |
|---|---|---|
| GSE307534 | 9.4G | OK（RAW.tar 含 GSM 级 tar.gz） |
| GSE274557 | 6.3G | OK |
| GSE235672 | 4.7G | OK（含 GBM.spatial.rds.gz） |
| GSE175540 | 1.1G | OK（含 bulk TPM csv.gz） |
| GSE202740 | 1.0G | OK |
| GSE242311 | 857M | OK（含 8 个 per-sample zip） |
| GSE111672 | 722M | OK（含 inDrop 矩阵独立文件） |
| GSE203612 | 708M | OK |
| GSE225857 | 608M | OK |
| GSE319536 | 339M | OK |
| GSE213699 | 338M | OK |
| GSE211956 | 285M | OK |
| GSE246011 | 198M | OK |
| GSE274103 | 178M | OK |
| GSE144240 | 162M | OK |
| GSE208253 | 154M | OK |
| GSE210041 | 123M | OK |
| GSE171351 | 116M | OK（含 combined_visium.h5ad.gz） |
| GSE226997 | — | 跳过：唯一文件为 44.2GB RAW.tar，见 access_control_needed 表 |

合计约 27 GB（GEO）+ 约 12 GB（other_sources）。

## 补充资源下载状态（2026-07-27 最终）

- HEST-1k：492 个人类肿瘤/健康样本（排除脑/脊髓），5109 文件全部到齐（约 940 GB），selected_samples.tsv。
- 10x 近 3 年肿瘤数据集：73 个全部到齐（592 GB），manifest_3yr.tsv。
- STOmicsDB：78 个肿瘤数据集（172 GB），manifest.tsv。
- SpatialDB：2 个癌症数据集（1.8 GB），manifest.tsv。
- CROST：HOLD（服务端宕机+数据量已足），恢复方法见 other_sources/crost/HOLD.txt。
- Synapse/HTAN 原始数据：HOLD，用户手动处理。

## 备注

- GEO 文件大小均已与 FTP 清单逐一核对一致；RAW.tar 抽查可读。
- 10x 数据集清单与 Dryad 条目为基于论文背景与 repo 元数据的推断（论文付费墙，Table S1 未能核实），
  若取得 Table S1 建议核对；详见 other_sources/sources_manifest.tsv notes。
- 受控数据（EGA、GSA/GSA-Human、NODE OEP001756、HTAN Synapse 原始数据、HEST-1k）
  汇总于 access_control_needed/access_control_datasets.tsv，需手动申请/下载。
