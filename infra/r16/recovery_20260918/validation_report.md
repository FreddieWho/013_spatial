# validation_report（T6，2026-09-19）

实际运行 21 测试，全绿：
- 本包 diagnostics 12（合成数学/软件问题说明，非生产修复证据）
- tests/test_r16_recovery_pc_matching.py 3（旧signed分组必失败项 + abs合并 + 定向平均conflict计数）
- tests/test_r16_recovery_laneA_stats.py 3（patient-first汇总 + 切片投票反例 + 同算子null）
- tests/test_r16_recovery_mask_tasks.py 3（窗与GT无关 + 隐藏改值特征不变 + 阴性窗参与）

未测（诚实列出）：仓库全量回归（313测试，历史产物）；变异函数匹配主surrogate（T2P2未做）；
正式条件推断（残差重抽样未做，主证据已换独立预测）；新registry派生（历史registry未改写）。
