# 文献中已讨论的空间场形态（research 笔记，2026-09-20）

三轮来源：①灵感论文 TLS pan-cancer atlas（Science 2026, Cho et al.，全文24页已通读）；
②外部文献（chemokine网络/COMMOT/反应扩散/SpatialDE/肿瘤边界）；
③对照本项目 M2 现状（线性各向同性邻域均值，只能看单调扩散）。

## 第一轮：灵感论文说了什么场形态

论文核心场形态就两种，成对出现：

1. **high-to-low（TLS→肿瘤）**：免疫激活、抗原呈递（IFN-α/γ、MHC-II）在 TLS 核心最高，
   随到 TLS 距离单调衰减。方法：Semla RadialDistance（TLS core<0，margin 0–50μm，
   超出 fallback 0–100μm）+ 距离 rescale 0–1 + Monocle natural spline（df=3）拟合 +
   曲线分类。8 癌种保守。
2. **low-to-high（TLS→肿瘤）**：MYC、G2/M、EMT、KRAS 在 TLS 旁最低，随距离单调上升。
   同一套距离框架，方向反过来。

关键细节：
- 距离只做到 0.75（标准化后），不是无限远；分析限 tumor-annotated spots。
- 正常组织里的 TLS：免疫通路保留 proximity 梯度，肿瘤增殖通路不保留——场形态依赖组织上下文。
- 成熟度维度：E→P→S，DC/FDC/GC 组织化递增，场形态随成熟度变（不是静态场）。
- 配套：HookNet-TLS + YOLOv8 从 H&E 直接检 TLS（3071 WSI，25088 TLS），TLS composite score
  分层预后。这是"场→临床"的闭环，本项目未涉及。

对本项目的直接映射：high-to-low/low-to-high 正好是 M2 现在**唯一能看**的形态
（单调线性扩散）。论文用 spline df=3 而我们用线性均值——连人家 2026 年的基线都不如。

## 第二轮：别的场形态（M2 看不见的）

1. **core-periphery zoning（TLS 内部）**：2026 年 8 月 Cancer Immunol Immunother
   （EcoTyper 71 cell states）：TLS 不是均匀团，active core + diffuse effector；
   同一 lineage 的不同 state 对 TLS 距离关系相反（有的富集、有的排斥）。
   形态含义：场在结构内部还有径向分带，不是单调出不去。M2 的 8 邻均值把内外抹平了。
2. **配体-受体定向流（COMMOT, Nat Methods 2023）**：collective optimal transport，
   带空间距离上限、多配体竞争、总量约束；输出是**有方向**的信号流（谁→谁），
   不是无向相关。还用 PDE 正演仿真验证。形态含义：场有方向（源→汇），M2 无向。
3. **反应扩散图灵斑（经典）**：自激活+侧抑制、不同扩散系数→自发斑图，无需预设源。
   形态含义：场可以没有"源位置"，M1/M2 全预设了"位置→表达"方向，反了。
4. **GP 空间场（SpatialDE 2018）**：表达=GP（空间核+非空间核），lengthscale 可学。
   形态含义：场的"尺度"本身是参数，不是固定的 8 邻。M2 的 k=8 是拍脑袋。
5. **肿瘤边界/侵袭前沿 niche**：EMT-like、hybrid E/M 富集在肿瘤-基质交界 0.5mm 内；
   不是从点源扩散，而是沿**边界线**的带状场。M2 的点邻域对线源天然不敏感。

## 第三轮：M2 升级候选（按性价比排）

1. **spline 距离曲线（抄灵感论文）**：对程序分做 natural spline df=3 ~ 到锚点距离，
   分类 high-to-low / low-to-high / flat。成本最低（OLS 变 GAM），直接对标 2026 Science。
2. **距离加权+方向分裂**：邻域按距离衰减加权；再按相对肿瘤边界内外分裂。
   抓"带状场"和"衰减场"，仍是线性回归，加两列特征的事。
3. **配体-受体定向（COMMOT-lite）**：只对 chemokine 家族（CXCL13/CXCR5/CCL19/CCL21/CCR7）
   跑轻量 OT，看方向是否从 TLS  core 指向外。成本中，需 CellChatDB。
4. **GP lengthscale**：每程序学一个空间尺度，回答"场有多大"。成本中（47片×通路数）。
5. **图灵/反应扩散正演**：最后考虑。无源场假设与本项目 M1/M2 框架冲突，先记 LEADS。

## 与本项目假设的接口

- 本文假设灵感（TLS 周围距离梯度）= 上述形态 #1，正是 M2 现在唯一能看的形态——
  说明当前阴性（M2≈0）可能只是眼睛不行，不是没有场。
- P-mhc2（MHC-II）正好落在论文 high-to-low 的核心通路里：L-009 输给 B 轴是"丰度打不过丰度"，
  但距离梯度形态（spline 曲线）还没测过——这是 mhc2 唯一剩下的翻盘路径，记 LEADS。
- 下一步建议：先上候选 #1（spline），一天工作量，直接回答"M2 看不见的是不是真没有"。
