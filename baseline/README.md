# 噻虫嗪光催化的跨批次预测与停机决策

本项目对 Xie 等人发表于 Chemical Engineering Journal 499 (2024) 156029 的实验进行计算再分析。原论文 DOI 为 https://doi.org/10.1016/j.cej.2024.156029 。数据源为出版商公开补充材料 https://ars.els-cdn.com/content/image/1-s2.0-S138589472407520X-mmc1.docx 。原始实验与实验设计归原作者。本包交付提取的数据、派生结果和本次分析代码。本文所称“达标”指原文设定的 3.558 mg/L。没有新增反应器实验。

## 复跑

在本目录使用 Python 3.13（实际运行版本 3.13.9）。使用随包两个 CSV 即可离线运行实验。首次安装库需要网络。

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run.py
```

可选，从出版商 Word 重新提取同一份数据，再运行分析。

```sh
python prepare.py
python run.py
```

如已有补充材料，使用 `python prepare.py --supplement 本地文件.docx`。本次源文件 SHA-256 为 `0d705ec7c35b8a81c9b85b6bdc803316126ed85aea0b5bfb8f19f0708a315bee`。图像为嵌入的原尺寸 JPEG，坐标校准依赖原图尺寸。终端将输出哈希和图表终点核对结果。核验现有结果可运行 `python verify.py`。

## 数据字典

* `data/table_s1.csv` 共32批。event 对应原文 E1–E32；catalyst 是表脚注中的催化剂类型（1 对应复合物、2 对应纯 CN）；pH 是表列 PV；c0 为 mg/L；dose_mg 为每批称取毫克数；light_level 为光强档位，作为类别处理。de_pct 是表列 180 min 降解百分数。ratio180=1−de_pct/100，ct180=c0×ratio180，单位 mg/L。
* `dose_as_printed` 保留投加量的原始文本。E25 的“50..2”在主分析中按 50.2 mg 解析，并单独分析剔除 E25 的结果。正文与补充表脚注对光强的物理单位写法不同，所以没有计算光能和电费。
* `data/figure_s8.csv` 有 160 个目标坐标（32 批×5 时刻），包含无法读取的空值。figure_ratio 表示图 S8 的 C(t)/C0 近似值；pixels 为候选像素数；status 指示可读性。60、90、120、150 分钟共 120 个可读点，8 个缺失。180 分钟有 28 个可核查读数，仅用于与表核对。建模时每批 180 分钟端点只使用表 S1。
* 读图方式为按各图坐标轴校准、分离图例颜色，在时间刻度附近定位符号像素，并用左右相邻线段外推符号中心。少数符号加宽范围的人工定位窗口写在 prepare.py 中。重叠且无法确认的点保留缺失。小数位反映像素换算，不能当作化验精度。未对曲线插值造标签。发表前宜让第二名读图者独立盲读或以仪器数据替代。
* 原文 329 行原始仪器记录未纳入。暗吸附后开灯为时间原点。当前可复算的 152 个标签包括 120 个时间点读数和 32 个表列终点换算值。

## 方法与文件对应

`analysis.py` 完成开灯前预测。折外指当前事件的浓度标签均不参与该折训练。主结果逐批留一，固定随机种子 20261126。终点稳健性把催化剂、光强档位、四舍五入的整数 pH、最近的名义投加量（25/50/100 mg）、名义浓度（5/10/20/50 mg/L）完全相同的批次一起留出，共 29 折。这是事先指定的数值组合规则，研究总体为探索性再分析，未做前瞻预注册。岭回归只以训练批次拟合标准化，惩罚系数固定为4；森林固定240棵树、最小叶样本3、全部候选特征。曲线森林按每批标签条数的倒数加权。五种模型见报告。随机森林直接拟合浓度比例，曲线版本显式使用分钟数；一级模型终点实验从各训练批终点反算有效速率，曲线实验用可读时点拟合一条对数线性曲线。两项任务各自重新拟合。

`early_sample.py` 在每批揭示一次60分钟图读浓度，再预测90至180分钟。四个对照看到相同的60分钟读数。未来预测误差全部排除60分钟点。过门槛的这次读数可以在60分钟直接触发停止，这些批次必须与以后由模型建议的停机分开解释。曲线内同批相关性是该决策任务所允许的信息，原始检测标签的来源仍是同一作者、同一研究。模拟测量偏移单独加在60分钟输入上，评估使用未移动的参考值。0.25和0.50 mg/L为探索性绝对裕量，没有概率覆盖承诺。

`diagnostics.py` 将相同折外预测作批次配对，导出图表、来源核查、分层指标和规则改变后同批运行时长的变化。固定折外差值的2万次有放回重采样只描述这一组误差的敏感度，未重训模型，不是外部验证置信区间。within_batch_shape_diagnostic 使用本批终点仅诊断一级曲线形状，不得作为新批次预测得分。

## 结果索引

| 需要回答的问题 | 文件（均在 results 目录） |
|---|---|
| 180 min 新批次预测误差及每批原始预测 | endpoint_metrics.csv、endpoint_predictions.csv、paired_endpoint_errors.csv |
| 相似配方留出、分层和剔除 E25 的敏感性 | endpoint_metrics.csv 的 nominal_recipe 行、endpoint_subgroups.csv、endpoint_excluding_e25.csv |
| 仅判断 180 min 达标会错哪些批次 | endpoint_decisions.csv、endpoint_decision_counts.csv |
| 60–180 min 冷启动曲线 | curve_predictions.csv、curve_metrics.csv |
| 最早候选停机及裕量 | stop_decisions.csv、stop_counts.csv、paired_margins.csv、fixed_time_reference.csv |
| 揭示一次中途测量后的四个对照 | early_sample_predictions.csv、early_sample_future_metrics.csv、early_sample_endpoint_metrics.csv、early_sample_decisions.csv、early_sample_counts.csv |
| 改变60分钟读数的影响 | early_sample_input_sensitivity.csv |
| 图表一致性和曲线形状诊断 | figure_table_check.csv、within_batch_shape_diagnostic.csv |
| 接近研究指定范围的实测工况和时间搜索 | supported_recipes.csv、supported_recipe_time_grid.csv |

预测结果中 observed_ct 为表值换算或图读近似值。指标中浓度 MAE 先在批内求平均再按批等权；RMSE 与 R² 按可读观测点汇总；端点每批只有一个点。模型名称 pooled、catalyst、c0_only、condition_k、forest 分别指共享速率、按催化剂分组的速率、只用初始浓度的速率回归、全工况速率回归和随机森林。early_kinetics 为60分钟单点估计速率；early_linear 为60分钟线性外推（最低截在零）；observed60_anchor_prior_k 将浓度校正到60分钟而保留原速率；observed60_then_prior 在60分钟可直接停止，后续沿用原工况预测。

决策表中 abstain 表示截至180分钟在候选时刻均未建议停机，unreadable 表示该选中时刻参考标签缺失；两者分别统计。borderline 标记读图值距离阈值不超过初始浓度1%的内部点。E20、E29 在初始已满足目标，另行标为0分钟，规则比较分母为30。avg_min_saved_verified_pass 仅在已建议且达标的批次中相对180分钟计算，各行入选集合可能不同。paired_margins 对同一批次逐项比较，更适合评价增加裕量的代价。

figures 目录提供 4 张重新绘制的 PNG 及矢量 PDF。所有模型和比较都是这次计算再分析的结果；算法参数未针对最终分数做搜索，单批失误和未建议停机的情况均已保留。原文 ALG 未按原划分复现，与本实验没有同条件分数可比。

## 中文论文初稿增补实验

本版保留原分析，增加单次检测的等信息对照、检测周转时间及符号中心读图敏感性。prepare.py 同时输出 data/figure_s8_marker60.csv，保存60分钟中心窄带的有色像素中位数；主分析仍用线段与符号联合读数。两者来自同一幅图，不代表独立仪器测量。仅替换输入，后续参考值保持一致。

turnaround.py 对两条策略都安排在第60分钟取一次样。assay_only 只对该样本已经达到目标的批次建议停灯。rate_update 还为样本未达标的批次按拟一级速率选择后续时间。报告延迟为0、30、60、90、120分钟，操作时刻从60、90、120、150、180分钟中选择；得到结果前不执行。未建议停灯的批次统一计入180分钟的观察窗口，不能据此推定可以排放。情景计算忽略检测费和180分钟之后的时间，单位为照射分钟。参考浓度不可读时仍保留建议和计划时长并单列人数。

* turnaround_anchor_inputs.csv：32批的两种60分钟输入。
* turnaround_decisions.csv：每批、每种读法、两种操作策略和五种延迟的逐项决策。reason 区分60分钟读数已达标和需要对之后时刻做预测。
* turnaround_counts.csv：建议数、可核验结果、弃权、全30批的时间账目、相对固定180分钟和相对仅检测策略的时间差。
* future_scores_by_population.csv：所有32批及60分钟仍未达标的22批各自的未来时间点评分，附替代读图结果。
* paired_future_errors.csv：对同批预测的绝对误差做配对；固定预测值后重抽批次20000次，不代表对模型重新训练或外部风险估计。
* manuscript_figures.py：供论文使用的检测延迟图及 E24 图读敏感性图，见 figures 中对应 PNG 和 PDF。

原有 summary 中只在成功停机样本上计算的时间，与新增表中按全部30批的时间有不同分母。run.py 包含全部新旧实验，并检查延迟、计数与主分析一致性。
