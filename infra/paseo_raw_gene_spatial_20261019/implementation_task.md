# 执行任务：同原始基因输入空间对照

指定 Pi/opencode-go/deepseek-v4.1-flash，max，cwd=/home/huyudi/013_spatial。读取当前STATUS、plan及 docs/FIELD_RAW_GENE_SPATIAL_CONTRACT_20261019.md。先遵守CodeGraph定位现有实现；你不是唯一工作者，不覆盖他人修改。

所有权：新建 scripts/field_raw_gene_spatial_20261019.py、scripts/run_field_raw_gene_spatial_20261019.sh、必要的定向检查及 infra/field_raw_gene_spatial_20261019/ 输出；本目录implementation_handoff.md归你。科学合同已由主代理写定，你可以提供疑点与具体数值实现附件，不直接修改根STATUS/plan/roadmap/decisions/ISSUES。既有冻结脚本/预测不改。

第一阶段仅实现、定向检查、小样本prepare/fit/evaluate/summarize贯通与成本估计，暂不跑全量。只采用科学合同的固定候选和对照，全部4目标与11划分、full15范围不静默缩减。用Git提交你的实现及小结果，交接准确文件范围和未提交状态，停止代码修改供Grok独立审核。

重点核查/测试：训练和留出单位/来源互斥、抽样/基因/归一化/标准化与既有Ridge一致、目标基因不入输入/分母、条码坐标图顺序、孤立点语义；3份复制X经块缩放后与同点Ridge等价；置换图确实破坏表达与邻域对应而保留图结构，输出条码不变；heldout目标变化不影响任何特征构建/预测；65/47单位语义及TENX13/14合并、部分指标缺失、不把零值当缺失。按实际风险做必要小检查，不默认重跑所有旧测试。

评价尽量复用Grok已审 Git5fee80f 的生物学描述模块，新的适配/聚合仍须审核。新增范围绝对误差应先逐片abs(delta)，再按完整单位合同聚合；不能abs(全局有符号均值)。输出逐片/单位/来源配对差值及完整分母。不得利用留出评分改参数。

检查现有图/数组实际形状与训练样本数后选择CPU求解器与线程，保留收敛状态/训练残差并核实数值精度。不得把不收敛解当完成。避免所有样本/重复邻域缓存无限复制；根据现有机器可用内存与磁盘预算执行，不用GPU。至少200GB空闲；不得访问人工TLS、历史人工评价、HookNet或外部数据。

交付implementation_handoff.md：实际provider/model/thinking/agent ID、Git范围、模型及精确参数、输入与样本范围、运行命令、定向检查/贯通结果、预计全量成本、疑点/NOT_RUN。若依赖或模型不可用报告阻塞，不换模型或静默降级。准备好后停止等待审核。
