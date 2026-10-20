# COLON MAP 身份解阻材料（草稿，发送/获取均需用户批准）

**目标：** 把 7 个样本 token（6723 / 7003 / 7319 / 7794 / 8270 / 8578 / 8899，共 36 张切片）
映射到患者号，将独立单元从 8 个扩大。这是当前最大的科学瓶颈。

## 已知锚点（仓库内可溯源）

- 两篇 Cell 论文：PMID 34910928（Cell 2021，PMC8941949）、PMID 38065082（Cell 2023，PMC10756562）。
- 相关但未确认关联：40473992（Nat Genet 2025，CRC Visium HD）。
- 本地已排除的通道：HEST 元数据无 donor 字段；CZ CELLxGENE 397 个 collection 无 COLON MAP；
  公开网页检索无标本-患者对照表（D-243 负结果）。

## 路径 A：补表自查（成本最低，先做）

PMC 开放获取全文的 supplementary materials 通常含 specimen 表。步骤：
1. 从 PMC8941949 / PMC10756562 的 supplementary 区下载全部附表（Excel/CSV）。
2. 在表内检索 7 个 token；若命中，提取 token→patient 映射。
3. 命中即记录来源（PMID + 表号 + 行），写入 sections.tsv 的 identity_note，
   patient_status 升为 known_within_study，按 SPLIT_FREEZE 失效条件重新生成 fold（版本 +1）。

**这需要外网获取，属于"补充外部数据"，执行前需用户批准。**

## 路径 B：作者邮件（路径 A 未果时）

草稿（英文，PMID 38065082 通讯作者）：

> Dear Dr. [通讯作者姓名],
>
> We are reanalyzing the Visium sections from your COLON MAP resource
> (Cell 2023, PMID 38065082; companion PMID 34910928) as distributed through HEST.
> The released metadata list seven sample identifiers — 6723, 7003, 7319, 7794,
> 8270, 8578, 8899 — but do not state which sections derive from the same donor.
> For patient-level statistics we need only the mapping from these sample IDs
> to anonymized donor IDs (e.g., as in your Supplementary Table Sx).
> Could you share this mapping, or point us to the table containing it?
> We will cite both papers and follow any data-use conditions.

发送前需填：通讯作者姓名/邮箱（从 PMC 页面通讯作者栏取）、Supplementary 表号占位符。

## 解阻后的连锁动作（预先登记，免得到时漏）

1. sections.tsv 更新 + fold 表重新生成（SPLIT_FREEZE 失效条件已声明）；
2. 独立单元数从 8 → 8+（36 片分属几名患者决定上限）；
3. 患者级统计功效重估，评估报告口径更新。
