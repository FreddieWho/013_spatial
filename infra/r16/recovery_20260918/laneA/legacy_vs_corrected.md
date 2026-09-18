# T2 legacy vs corrected（诊断集，6患者x1片x10基因，占位null）

## 机械修复验证（通过）
- SYNTH_NEG（纯噪声）：pmax=0.12, pcon=0.29 —— 不显著，机械未制造假阳性。
- SYNTH_BUMP（区内+1.5z注入）：tcon=+1.53, pcon=0.005 —— 固定对比检出注入；
  tmax p=0.215 —— 峰高统计量对宽 bump 不敏感（符合预期：max 追逐单箱噪声）。
- 同算子闭环：观察与 199 联合 draws 同一 patient-first+统计量算子（A05 关闭机械部分）。

## 与 legacy 口径的对照（CXCL13~B）
- legacy（D-126 口径）：全队列47片中位曲线峰高 + section-draw混合null + 切片投票 agree。
- corrected（本诊断集，占位null）：t_contrast=+1.08（区内-区外），t_max=2.06。
- 口径差异：legacy 追单箱峰位，corrected 追固定区内外对比；两者都看到 CXCL13 区内高，
  但 corrected 不回答“峰在哪一跳”（按03 §3，峰位仅辅助）。

## remap 近似性质实测（A08/A09 首次实测数，6片x199 draws）
- accept_rate=1.00（6/6片全接受；诊断片几何较规则，全量47片待测）。
- frac_within_tol 中位 0.71-0.84；unique_src_frac 中位 0.72-0.84（多对一确认，A08实锤量化：约20-28%源点被复用）。
- moran_ratio（surrogate/obs）中位 0.63-8.96，大片稳定（0.6-1.9），小 Moran obs 片（0.01）比值爆炸——
  说明 remap 不保变异函数尺度：surrogate 的平滑度与原场不一致。
- 结论：旧旋转 remap 降为诊断/敏感性对照（02要求），主显著性证据须等变异函数匹配 surrogate。
  本文件所有 p 值仍系占位 null（bin置换），记 NOT_CALIBRATED，不得作空间显著性证据。

## 遗留警告
- All-NaN RuntimeWarning：小片部分距离箱空箱（ns=0），nanmedian 已正确忽略；后续正式版加显式有效箱记录。
- MUC2 tcon=+0.31 pcon=0.02：占位null下擦线显著，不作任何解读（null未校准）。
