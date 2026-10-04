# 患者身份元数据补查

用户已明确批准≤10MB、仅现有两来源及直接关联身份元数据的补查。本次完成有限范围检索，**未取得26个已分析样本的明确患者对应；I-033未解除**。这表示所取得记录没有交叉表，不表示任何地方都不存在该表。

| 来源 | 取得材料与检查结果 |
|---|---|
| [USZ Zenodo14620362](https://zenodo.org/records/14620362) | 官方JSON列出KC1–3/LC1–5及3个肾、5个肺肿瘤，未给患者分组字段。公开附件仅2.09GB数据压缩包，本轮未下载。样本数不直接转为独立患者数。 |
| [GSE175540](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175540) | 官方FTP SOFT成功取得，包含24个样本；已分析18个GSM全部找到，疾病/组织/样本别名字段无明确患者对应。未下载表达附件。 |
| [PRJNA732692](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA732692) | Entrez取得直接关联BioProject XML，未提供样本到患者交叉表。 |
| [关联论文书目](https://pubmed.ncbi.nlm.nih.gov/35231421/) | Entrez XML可用，仅含书目与研究摘要，没有样本交叉表或PMC全文编号。Cell直接关联全文/补充入口返回403，补充资料内容未能核查；不能写成已穷尽论文全部附件。 |

GEO网页请求超出单响应2MB上限后丢弃，改取官方FTP SOFT；PubMed网页返回浏览器验证页，未当作有效metadata。四份有效元数据合计29956字节。脚本记录的响应体累计2035547字节（包含丢弃响应和浏览器验证页），低于10MB；浏览工具另有少量页面/错误文本读取。没有获取表达矩阵、影像或联系他人，没有使用GPU。

[逐样本核查表](../infra/bioinf-data-index/raw/identity_20261003/sample_identity_audit.tsv)保留26行、患者字段均空白。[检索凭证](../infra/bioinf-data-index/raw/identity_20261003/retrieval.json)记录URL、时间、大小、SHA和失败原因；[数据索引](../infra/bioinf-data-index/index.tsv)、manifest与summary已同步登记4份有效元数据。浏览器验证文件只作为技术诊断，未列入有效数据索引。

这一步服务H-01/H-05及H-08的独立重复要求，但没有新增生物学支持。原全量结果和冻结identity不改写，未重跑表达分析。I-031的空间检验适用性/功效、短窄范围及肿瘤界面真值的缺口也不因此改变。下一项真正有用的身份依据是作者提供的明确样本—匿名患者表，或可追溯的样本独立性确认；当前未获授权联系作者，且本轮没有发送消息。
