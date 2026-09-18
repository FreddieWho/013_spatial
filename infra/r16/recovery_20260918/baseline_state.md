# T0 基线状态（2026-09-18，主 agent 实测）

## 可直接复用
- cache 47 片 npz：原始 counts 原样拷贝（builder 只做 mat@sel），行和=library size 可算。
- loadings npz 65M：T1 只重算 matching，不重训 PCA/LISA。
- laneA/B 六＋一 JSON：历史记录保留，否证外推挂 METHOD_LIMITED（I-022/I-023/I-024）。
- ST-CRC 14/14、USZ 8/8 readers＋manifests：T4 可用；已知轴已在两队列打分过（exposure 必须披露）。

## 已知未知
- 深度混杂是否成立：未知（T3 敏感性分析回答；I-025）。
- remap 真实接受率：未知（历史 info 被丢弃；T2 重测；I-022）。
- 符号修复后 NEW 数：未知（T1 回答；I-023）。

## 受影响的旧结论
- D-118 NEW=0、D-126 全线阴性：外推效力暂停，数字保留。
- D-127 的 191：候选身份保留，"需排除亚型/状态"改为分层评价（04-B）。
- R-04 关闭、六轴描述性结果、四臂执行事实、PRECAST/S 跳过：不受影响。

## 缺口（不编造）
- 无独立结构锚点的程序只能做场预测（T5 合同 A）；USZ 跨癌种仅验共享部分。
- ILC 已退役；新 marker 轴需单独报告。
