# External bioinformatics metadata index

This directory records the approved metadata-only external inputs used for
R-01 identity resolution. `index.tsv` is the checksum inventory,
`manifest.json` records provenance and scope incidents, and `summary.json`
provides machine-readable retention totals.

Only official metadata and an official patient-to-GSM identity crosswalk are
retained here. When explicitly requested, R-04 also records official Zenodo
source metadata JSON and hashed pointers to existing local reference files;
the latter are read in place with an allowlist and are never copied into this
repository. No image, expression matrix or biological result file is retained.
Local GEO raw archives outside this directory are used only as header-level
GSM locators.

## 2026-10-03 已批准身份元数据补查

新增4份官方元数据（Zenodo14620362、GSE175540 SOFT、PRJNA732692、关联论文书目），见 `raw/identity_20261003/`。26个当前样本均未获得明确患者对应，不解除患者级确认限制。浏览器验证响应仅作访问诊断，不计入有效metadata；失败请求与累计字节登记于retrieval.json。详见[补查报告](../../docs/TLS_PATIENT_METADATA_LOOKUP_20261003.md)。

2026-10-03用户提供`gse175540_sd.xlsx`已保留于`raw/user_gse175540_20261003/`。24个GSM对应23个表内ID，当前18片对应17个ID；c_2有Frozen/FFPE两片。原出版来源及ID的患者含义待图注确认，非正式患者交叉表。
