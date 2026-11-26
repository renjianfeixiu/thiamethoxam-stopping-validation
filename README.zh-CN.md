# 噻虫嗪光催化停机决策　数据与复现

新增的绝对浓度停机验证位于 external/。Chen 等的聚酯纤维工况和 Farooq 等的负载量工况用于检验固定规则。图源、单位换算、判定口径、重算方法与逐工况结果见 [external/README.md](external/README.md)。顶层重算命令已包含这项验证。

本文是公开实验的回顾性二次分析。baseline 对应原研究的 32 批重建、整批留出预测、同次检测对照及回报延迟。extension 对应同批取样时机比较与跨论文浓度轨迹检验。研究目标 3.558 mg/L 来自 Xie 等的 ECOSAR 评价，适用于本文情景。外部文献评估使用相对浓度。

## 离线重算

在新环境中运行以下命令。首次安装依赖需要网络。随包的 CSV 已可用于完整分析，无需购买原文。

```sh
python -m pip install -r requirements.txt
python reproduce.py
```

使用 Python 3.13、NumPy 2.3.5、pandas 2.3.3、scikit-learn 1.7.2。baseline/README.md 详述基础实验。重新运行生成的表会写入各自的 results/。固定种子重采样的范围仅描述本组批次差值的变动。图表小数位多于文中展示位数，不能视为仪器精度。

## 论文数字对应

| 研究问题 | 可复算文件 |
|---|---|
| 全 30 批在 60 分钟取样的直接建议及额外照射时间 | baseline/results/turnaround_counts.csv 和 turnaround_decisions.csv |
| 同次测量在四种预测方法中的作用 | baseline/results/future_scores_by_population.csv 和 paired_future_errors.csv |
| 取样时刻 30、60、90 分钟对共同批次的影响 | extension/results/timing_metrics.csv、timing_counts.csv 及逐批 timing_decisions.csv |
| 两次取样方案及固定逐批误差的重采样 | extension/results/paired_endpoint_contrasts.csv |
| 外部曲线按来源分组的误差与覆盖数 | extension/results/transfer_metrics.csv、transfer_forecasts.csv |
| 附加的保护裕量敏感性 | extension/results/guard_counts.csv、guard_decisions.csv、transfer_decisions.csv |

在 timing_counts.csv 选 scope=identical_anchor_cohort、read=figure_ratio、method=exp、delay=0 即得论文中 27 批停机比较。相同范围的终点误差取 timing_metrics.csv 的 target=180，分母为 29。forecast 包括恰在 180 分钟建议的批次，forecast_early 仅含早于 180 分钟的建议。未建议的批次在窗口内按 180 分钟记账，其达标状态依然未知。表内的时间减少值均是计划量，个别情景含错误建议，不可作为实际节能量。外部曲线的“hold”意为沿用测量值，“linear”按初始至取样的平均绝对降速作线性外推并截在零，“exp”按相同比例降速作一级外推。外部案例按可获得的图片及时间范围选择，不构成系统综述。可读数以表内 n_series 为准。对于无催化剂对照单独列出。

## 原始来源与坐标提取

- Xie 等，2024。[出版记录](https://doi.org/10.1016/j.cej.2024.156029) · [补充材料](https://ars.els-cdn.com/content/image/1-s2.0-S138589472407520X-mmc1.docx)。表 S1 提供终点，图 S8 提供中途浓度。正文图 6f/g 和补充图 S6 属于同一研究。
- Farooq 等，2025。[Frontiers in Materials](https://www.frontiersin.org/journals/materials/articles/10.3389/fmats.2025.1677416/full)。图 7b 含 Ni-PbS、PbS 和空白。图 7d 未重复纳入。
- Arfanis 等，2023。[Catalysts](https://www.mdpi.com/2073-4344/13/3/516)。仅图 1a 的噻虫嗪；低于图面可分辨范围的重叠点保留空值。
- Emmanuel 与 Dilek，2024。[American Journal of Chemical Engineering](https://www.sciencepublishinggroup.com/article/10.11648/j.ajche.20241203.14)。取图 5 和 6，纯 TiO₂ 和珍珠岩负载 TiO₂，各两种反应器。将图示 30 分钟定为开灯时刻，以该时刻读数重新归一化，原始时间和分母均见 extension/data/reactor_lamp_origin.csv。

如需从原始图重建，运行 `python extension/source_extract.py --supplement 原论文补充材料.docx`，再运行 `python extension/transfer_extract.py --paper 原论文正文.pdf --supplement 原论文补充材料.docx`。后者需额外安装 PyMuPDF 并联网获取公开图。图像校准依赖原始尺寸；图片哈希见 extension/SOURCE_SHA256.txt。各算法、遮挡点、读取窗口和具体读数均可追溯。原始论文和出版商图像请由原出处获取。
