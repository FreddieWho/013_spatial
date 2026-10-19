# 实测场资格与冻结预测恢复：描述性清单

只描述，不判定成功/失败，不设生物学通过线，不给出 p 值。
长度单位一律为该切片经验点间距（pitch）；分位阈值是操作阈值，不是生物真值。
三本账分开：测量资格（coverage/支持）、输入回声关系（与可见模块场）、
已开发预测的探索性再分析（实测 vs 预测）。三者的列名与 evidence_role 不同，不共用“成功”一词。

## 0. 语料与身份

- 主清单：来源留出 `infra/field_source_holdout_20261018` 的 1000 份冻结预测包（内部 72 片/65 单位、来源 53 片/47 单位/6 来源）；`train_indices.npz` 不进入评分。
- 实测清单覆盖全部 72 片 × 4 目标；来源身份不足的 19 片保留在内部描述，不进入来源分母。
- 研究内确认身份 14 片 / 8 个单位；其余为未知单位，不升级为独立患者。
- 状态计数（实测）：{"DESCRIBED": 2, "DESCRIBED_WITH_NOT_ESTIMABLE_METRICS": 286}
- 状态计数（预测）：{"DESCRIBED": 6, "DESCRIBED_WITH_NOT_ESTIMABLE_METRICS": 3994, "NOT_ESTIMABLE_CONSTANT_PREDICTION": 1000}

## 1. 分子检测覆盖（连续报告，不分档删除）

- B_AXIS：72 片；检出基因数中位 47；任意目标基因检出点占比 中位 0.999（0.612–1.000）；有效基因数中位 2.87（最小 1.19）；首位基因占比中位 0.535（最大 0.918）
- HALLMARK_OXIDATIVE_PHOSPHORYLATION：72 片；检出基因数中位 191；任意目标基因检出点占比 中位 1.000（0.919–1.000）；有效基因数中位 75.21（最小 4.20）；首位基因占比中位 0.055（最大 0.482）
- HALLMARK_INTERFERON_GAMMA_RESPONSE：72 片；检出基因数中位 195；任意目标基因检出点占比 中位 1.000（0.950–1.000）；有效基因数中位 24.66（最小 10.83）；首位基因占比中位 0.120（最大 0.269）
- HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION：72 片；检出基因数中位 199；任意目标基因检出点占比 中位 1.000（0.962–1.000）；有效基因数中位 32.79（最小 7.66）；首位基因占比中位 0.100（最大 0.299）

完整 72 片 × 4 目标逐片数值见 `measured_fields.tsv`；分层中位/极值见 `measured_summary.tsv`。

## 2. 恢复描述（单位宏平均，整臂；不跨信息预算判胜负）

| 臂 | 目标 | 视图 | 模型 | 模型类别 | 单位 | RMSE | 去均值RMSE | 逐点相关 | 局部图样相关 | q75 Jaccard | q75 打乱中位 | 质心位移(pitch) | 自相关差 |
|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| internal | B_AXIS | full15 | fixed_local_hgb | module_local_hgb | 66 | 0.1316 | 0.1006 | 0.5142 | 0.2613 | 0.3375 | 0.1427 | 1.7303 | 0.0385 |
| internal | B_AXIS | full15 | fixed_point_hgb | module_point_hgb | 66 | 0.1357 | 0.1023 | 0.5059 | 0.2597 | 0.3357 | 0.1427 | 1.5376 | -0.0171 |
| internal | B_AXIS | full15 | gene_ridge | raw_gene_ridge | 66 | 0.1035 | 0.0880 | 0.6329 | 0.3455 | 0.4236 | 0.1430 | 1.1699 | -0.0219 |
| internal | B_AXIS | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 66 | 0.0972 | 0.0802 | 0.6828 | 0.3515 | 0.4528 | 0.1427 | 1.1563 | 0.2167 |
| internal | B_AXIS | full15 | train_mean | constant_prior | 66 | 0.1691 | 0.1166 | NA | NA | NA | NA | NA | NA |
| internal | B_AXIS | small5 | fixed_local_hgb | module_local_hgb | 66 | 0.1512 | 0.1096 | 0.4140 | 0.2056 | 0.2840 | 0.1427 | 1.9839 | 0.0385 |
| internal | B_AXIS | small5 | fixed_point_hgb | module_point_hgb | 66 | 0.1498 | 0.1100 | 0.4238 | 0.2231 | 0.2894 | 0.1430 | 1.8748 | -0.0164 |
| internal | B_AXIS | small5 | gene_ridge | raw_gene_ridge | 66 | 0.1186 | 0.0980 | 0.5551 | 0.2757 | 0.3735 | 0.1433 | 1.3503 | -0.0349 |
| internal | B_AXIS | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 66 | 0.1105 | 0.0877 | 0.6195 | 0.2874 | 0.4092 | 0.1429 | 1.2932 | 0.1950 |
| internal | B_AXIS | small5 | train_mean | constant_prior | 66 | 0.1691 | 0.1166 | NA | NA | NA | NA | NA | NA |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | fixed_local_hgb | module_local_hgb | 66 | 0.1435 | 0.1085 | 0.6882 | 0.4488 | 0.4384 | 0.1432 | 1.3737 | -0.0018 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | fixed_point_hgb | module_point_hgb | 66 | 0.1490 | 0.1126 | 0.6820 | 0.4537 | 0.4322 | 0.1424 | 1.3882 | -0.0409 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | gene_ridge | raw_gene_ridge | 66 | 0.0803 | 0.0715 | 0.8498 | 0.6257 | 0.6086 | 0.1427 | 0.6459 | -0.0266 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 66 | 0.0781 | 0.0690 | 0.8667 | 0.6124 | 0.6189 | 0.1423 | 0.8339 | 0.2184 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | train_mean | constant_prior | 66 | 0.2324 | 0.1522 | NA | NA | NA | NA | NA | NA |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | fixed_local_hgb | module_local_hgb | 66 | 0.1609 | 0.1165 | 0.6617 | 0.4211 | 0.4130 | 0.1426 | 1.5818 | 0.0021 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | fixed_point_hgb | module_point_hgb | 66 | 0.1628 | 0.1206 | 0.6477 | 0.4321 | 0.3985 | 0.1431 | 1.5382 | -0.0427 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | gene_ridge | raw_gene_ridge | 66 | 0.1080 | 0.0900 | 0.7794 | 0.5176 | 0.5280 | 0.1427 | 0.8568 | -0.0407 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 66 | 0.0991 | 0.0799 | 0.8210 | 0.5182 | 0.5608 | 0.1430 | 1.0253 | 0.1953 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | train_mean | constant_prior | 66 | 0.2324 | 0.1522 | NA | NA | NA | NA | NA | NA |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | fixed_local_hgb | module_local_hgb | 66 | 0.0729 | 0.0561 | 0.7621 | 0.5056 | 0.4564 | 0.1426 | 0.7649 | 0.0405 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | fixed_point_hgb | module_point_hgb | 66 | 0.0783 | 0.0587 | 0.7495 | 0.4844 | 0.4451 | 0.1428 | 0.7692 | -0.0129 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | gene_ridge | raw_gene_ridge | 66 | 0.0647 | 0.0531 | 0.7964 | 0.4908 | 0.4971 | 0.1423 | 0.6519 | -0.0143 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 66 | 0.0620 | 0.0498 | 0.8195 | 0.4861 | 0.5231 | 0.1423 | 0.6993 | 0.2384 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | train_mean | constant_prior | 66 | 0.1456 | 0.0913 | NA | NA | NA | NA | NA | NA |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | fixed_local_hgb | module_local_hgb | 66 | 0.0921 | 0.0641 | 0.6756 | 0.4723 | 0.3873 | 0.1428 | 1.0223 | 0.0366 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | fixed_point_hgb | module_point_hgb | 66 | 0.0927 | 0.0643 | 0.6841 | 0.4577 | 0.3866 | 0.1423 | 0.9837 | -0.0091 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | gene_ridge | raw_gene_ridge | 66 | 0.0751 | 0.0605 | 0.7431 | 0.4293 | 0.4420 | 0.1427 | 0.7299 | -0.0286 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 66 | 0.0706 | 0.0549 | 0.7803 | 0.4316 | 0.4803 | 0.1430 | 0.7864 | 0.2140 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | train_mean | constant_prior | 66 | 0.1456 | 0.0913 | NA | NA | NA | NA | NA | NA |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | fixed_local_hgb | module_local_hgb | 66 | 0.0790 | 0.0563 | 0.9064 | 0.6734 | 0.5744 | 0.1429 | 0.4781 | 0.0441 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | fixed_point_hgb | module_point_hgb | 66 | 0.0788 | 0.0572 | 0.9047 | 0.6567 | 0.5776 | 0.1431 | 0.4599 | 0.0087 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | gene_ridge | raw_gene_ridge | 66 | 0.0599 | 0.0494 | 0.9259 | 0.6760 | 0.6449 | 0.1429 | 0.3134 | 0.0114 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 66 | 0.0621 | 0.0518 | 0.9208 | 0.6483 | 0.6334 | 0.1423 | 0.5251 | 0.2877 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | train_mean | constant_prior | 66 | 0.2226 | 0.1402 | NA | NA | NA | NA | NA | NA |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | fixed_local_hgb | module_local_hgb | 66 | 0.0854 | 0.0622 | 0.8872 | 0.6232 | 0.5485 | 0.1428 | 0.5105 | 0.0486 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | fixed_point_hgb | module_point_hgb | 66 | 0.0877 | 0.0643 | 0.8798 | 0.5946 | 0.5479 | 0.1428 | 0.4743 | -0.0043 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | gene_ridge | raw_gene_ridge | 66 | 0.0672 | 0.0554 | 0.9084 | 0.6248 | 0.6059 | 0.1426 | 0.3438 | -0.0006 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 66 | 0.0674 | 0.0554 | 0.9094 | 0.6074 | 0.6099 | 0.1424 | 0.5506 | 0.2701 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | train_mean | constant_prior | 66 | 0.2226 | 0.1402 | NA | NA | NA | NA | NA | NA |
| source | B_AXIS | full15 | fixed_local_hgb | module_local_hgb | 47 | 0.1347 | 0.0965 | 0.4125 | 0.2073 | 0.2939 | 0.1428 | 1.8746 | 0.0349 |
| source | B_AXIS | full15 | fixed_point_hgb | module_point_hgb | 47 | 0.1375 | 0.0968 | 0.4115 | 0.2043 | 0.2917 | 0.1427 | 1.7848 | -0.0166 |
| source | B_AXIS | full15 | gene_ridge | raw_gene_ridge | 47 | 0.1135 | 0.0867 | 0.5402 | 0.2550 | 0.3645 | 0.1427 | 1.3126 | -0.0280 |
| source | B_AXIS | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.1049 | 0.0760 | 0.6128 | 0.2670 | 0.4044 | 0.1427 | 1.1945 | 0.1709 |
| source | B_AXIS | full15 | train_mean | constant_prior | 47 | 0.1611 | 0.1013 | NA | NA | NA | NA | NA | NA |
| source | B_AXIS | small5 | fixed_local_hgb | module_local_hgb | 47 | 0.1435 | 0.1000 | 0.3337 | 0.1470 | 0.2549 | 0.1427 | 2.0591 | 0.0267 |
| source | B_AXIS | small5 | fixed_point_hgb | module_point_hgb | 47 | 0.1432 | 0.0996 | 0.3632 | 0.1845 | 0.2633 | 0.1432 | 1.9598 | -0.0134 |
| source | B_AXIS | small5 | gene_ridge | raw_gene_ridge | 47 | 0.1212 | 0.0909 | 0.4848 | 0.2131 | 0.3350 | 0.1433 | 1.3756 | -0.0321 |
| source | B_AXIS | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.1138 | 0.0808 | 0.5622 | 0.2268 | 0.3744 | 0.1425 | 1.3191 | 0.1626 |
| source | B_AXIS | small5 | train_mean | constant_prior | 47 | 0.1611 | 0.1013 | NA | NA | NA | NA | NA | NA |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | fixed_local_hgb | module_local_hgb | 47 | 0.1661 | 0.1119 | 0.6672 | 0.4087 | 0.4307 | 0.1429 | 1.6178 | -0.0025 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | fixed_point_hgb | module_point_hgb | 47 | 0.1455 | 0.1095 | 0.6538 | 0.4026 | 0.4153 | 0.1425 | 1.5164 | -0.0478 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | gene_ridge | raw_gene_ridge | 47 | 0.0876 | 0.0692 | 0.8304 | 0.5766 | 0.5975 | 0.1428 | 0.7924 | -0.0303 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.0835 | 0.0634 | 0.8604 | 0.5755 | 0.6222 | 0.1424 | 0.8994 | 0.1865 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | train_mean | constant_prior | 47 | 0.2329 | 0.1378 | NA | NA | NA | NA | NA | NA |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | fixed_local_hgb | module_local_hgb | 47 | 0.1905 | 0.1113 | 0.6528 | 0.3686 | 0.4191 | 0.1429 | 1.5962 | 0.0098 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | fixed_point_hgb | module_point_hgb | 47 | 0.1708 | 0.1143 | 0.6221 | 0.3846 | 0.3876 | 0.1431 | 1.5929 | -0.0520 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | gene_ridge | raw_gene_ridge | 47 | 0.1145 | 0.0823 | 0.7693 | 0.4948 | 0.5310 | 0.1428 | 0.9205 | -0.0408 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.1079 | 0.0723 | 0.8178 | 0.5033 | 0.5748 | 0.1431 | 1.0521 | 0.1684 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | train_mean | constant_prior | 47 | 0.2329 | 0.1378 | NA | NA | NA | NA | NA | NA |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | fixed_local_hgb | module_local_hgb | 47 | 0.0864 | 0.0606 | 0.6223 | 0.4546 | 0.3777 | 0.1427 | 1.1186 | 0.0354 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | fixed_point_hgb | module_point_hgb | 47 | 0.0884 | 0.0584 | 0.6665 | 0.3772 | 0.4063 | 0.1428 | 0.9369 | -0.0173 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | gene_ridge | raw_gene_ridge | 47 | 0.0703 | 0.0526 | 0.7354 | 0.3855 | 0.4575 | 0.1428 | 0.7741 | -0.0192 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.0664 | 0.0475 | 0.7876 | 0.3931 | 0.5071 | 0.1425 | 0.8217 | 0.1933 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | train_mean | constant_prior | 47 | 0.1633 | 0.0792 | NA | NA | NA | NA | NA | NA |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | fixed_local_hgb | module_local_hgb | 47 | 0.1175 | 0.0631 | 0.6043 | 0.4270 | 0.3637 | 0.1424 | 1.2572 | 0.0377 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | fixed_point_hgb | module_point_hgb | 47 | 0.1145 | 0.0623 | 0.6366 | 0.3591 | 0.3771 | 0.1426 | 1.0148 | -0.0132 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | gene_ridge | raw_gene_ridge | 47 | 0.0857 | 0.0595 | 0.6694 | 0.3355 | 0.3974 | 0.1424 | 0.8552 | -0.0290 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.0803 | 0.0519 | 0.7314 | 0.3458 | 0.4451 | 0.1425 | 0.9343 | 0.1748 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | train_mean | constant_prior | 47 | 0.1633 | 0.0792 | NA | NA | NA | NA | NA | NA |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | fixed_local_hgb | module_local_hgb | 47 | 0.0830 | 0.0653 | 0.8796 | 0.6352 | 0.5473 | 0.1426 | 0.6520 | 0.0402 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | fixed_point_hgb | module_point_hgb | 47 | 0.0826 | 0.0640 | 0.8839 | 0.6104 | 0.5694 | 0.1432 | 0.5948 | 0.0066 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | gene_ridge | raw_gene_ridge | 47 | 0.0692 | 0.0543 | 0.9107 | 0.6279 | 0.6298 | 0.1426 | 0.4308 | 0.0032 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.0696 | 0.0549 | 0.9122 | 0.6144 | 0.6279 | 0.1420 | 0.6333 | 0.2512 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | train_mean | constant_prior | 47 | 0.2065 | 0.1400 | NA | NA | NA | NA | NA | NA |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | fixed_local_hgb | module_local_hgb | 47 | 0.0915 | 0.0673 | 0.8705 | 0.6068 | 0.5543 | 0.1431 | 0.6452 | 0.0418 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | fixed_point_hgb | module_point_hgb | 47 | 0.0933 | 0.0683 | 0.8675 | 0.5771 | 0.5665 | 0.1428 | 0.5837 | -0.0014 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | gene_ridge | raw_gene_ridge | 47 | 0.0717 | 0.0576 | 0.9003 | 0.5990 | 0.6087 | 0.1427 | 0.4152 | -0.0021 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.0713 | 0.0570 | 0.9048 | 0.5899 | 0.6161 | 0.1426 | 0.6306 | 0.2418 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | train_mean | constant_prior | 47 | 0.2065 | 0.1400 | NA | NA | NA | NA | NA | NA |

打乱中位是保持取值边际、打乱空间排列的描述性扰动参考；不保持空间支撑与密度，未校准，
不构成零模型、不产生 p 值、不用于筛选关系。
`gene_ridge_smooth_half` 的自相关上升是构造结果（后处理平滑），不得作为空间关系或恢复证据；
本阶段不做同点/邻域基线对照，也不在信息预算之间判胜负。

## 3. 边界与未运行

- 分位阈值（q75/q80/q90）是操作阈值，取严格大于分位值；逐场报告 area_frac 与 at_threshold_frac，
  不同场的点集大小可以不同。
- 预测场是可见模块/基因的函数；预测与可见模块的相关部分由输入结构决定，不能读作独立生物学证据。
- 相关、自相关与尺度在任一方方差低于合同 sd_floor 时为 NOT_ESTIMABLE，不填 0、不填 1；RMSE 仍保留。
- 数学上可算不等于有足够生物学参照：检出覆盖、单基因主导、有效基因数、边界/孤立点支撑与未确认的物理尺度
  仍会限制解释；这些证据逐项并列，不合成单一分数，也不用统一百分比删除切片/目标/行。
- 不得把"可计算"格子升级为已具生物学评价资格；NOT_TESTABLE/NOT_ESTIMABLE 也不读作场不存在或恢复失败。
- 生物学可用性由主代理结合上述连续证据判断，本清单只提供描述性事实与不可评价状态。
- 未训练/未重训、未做独立 marker 隔离、未读取人工 TLS、未做正式显著性、未做物理尺度换算。
- 已知患者仅 8 个研究内单位；切片数不冒称患者数。

