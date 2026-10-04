# 用户提供的GSE175540补充表：可补信息与边界

文件`tmp/gse175540_sd.xlsx`，130行记录、12列；原始文件已保留并登记SHA。来源标记为用户提供，出版链接及原表图注待确认；工作簿作者属性不作为官方来源认证。

## 能补什么

GEO官方Sample_title的24个样本别名均能与表内ID及Spatial transcriptomics制备方式一致匹配，对应23个不同ID。当前已分析18片全部匹配，对应17个不同ID。这里使用官方标题和表内制备字段联合核对，不靠本地文件名或样本顺序推患者。

最关键的一行是Sheet1第74行：BIONIKK / c_2 / pT1a / Frozen / FFPE。GSM5924030的官方别名是ffpe_c_2，GSM5924050是frozen_c_2；两者均在当前18片内，应视为同表ID的保守分组候选，不能默认独立。当前18片包括BIONIKK14片（13个ID）、ExhauCRF2片、IMM2片。

还可以补充表中cohort、pT和各检测的开展情况。x/X表示表中列出的检测记录，不能改解释成生物学阳性；空白也不是阴性。表未给出实际TLS成熟度分数、细胞丰度、空间边界或治疗响应。

## 对既有结果的影响

旧逐切片留出并非按此ID分组留出；需要在后续分组验证中一起留出c_2的两片。1000μm时两片都有可评正负查询，已有逐片预测不得升级为独立患者验证。跨来源留出时两片同属GSE，不能由此断言发生了跨来源训练/测试泄漏；但切片汇总仍可能重复加权同一表ID。未重跑模型，旧数值保留，不声称影响已定量修正。

表头只有ID，没有明确patient/donor定义，因此现在可填table_id及保守分组，正式patient_id暂不赋值。若原表图注明确每行代表一位患者，可将本轮GSE部分升级为17个患者分组，随后按患者分组重算适用验证。该表不涉及USZ8片，也不解除空间校准或桥真值缺口。

[24样本对应表](../infra/bioinf-data-index/raw/user_gse175540_20261003/gsm_table_id_crosswalk.tsv)、[完整表格转录](../infra/bioinf-data-index/raw/user_gse175540_20261003/table_metadata.tsv)、[来源与SHA记录](../infra/bioinf-data-index/raw/user_gse175540_20261003/audit.json)。这一步补充证据资格，未新增H-01/H-05/H-08的生物学确认；不改变科学假设。

**后续确认（D-184）：** 用户已明确确认表内ID为患者；当前GSE18片为17位患者，c_2合组。此前“ID单位待确认”为历史状态，现已解除；原发表出处仍未独立核验。患者交叉表见`infra/bioinf-data-index/raw/user_gse175540_20261003/gsm_patient_crosswalk_confirmed.tsv`。模型未按患者重跑，USZ身份缺口不变。
