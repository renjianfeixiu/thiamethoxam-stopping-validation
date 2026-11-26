# 外部停机规则验证

本分析回放两篇论文中含催化剂的全部目标面板曲线。Chen 等图 5(b) 的三种纤维材料为本次新增外部来源。Farooq 等图 7(d) 的四个负载量工况为本次新增面板。先前分析用过同一研究的图 7(b)，不能把两个面板合并充当独立重复。空白和裸载体不属于含活性催化剂的工况；每一条图示曲线作为一个工况。参考时刻保留到照射 180 min，每个工况取 30、60、90、120、150、180 min 六个图示点，共 42 个。主判停计算仅用 60 min 锚点和初始浓度。其他点用于对照，不给决策程序作为预测输入。

## 图源与量纲

- Chen Y, Lu W, Shen H, et al. Applied Catalysis B: Environmental, 258 (2019), 117960. [原论文](https://doi.org/10.1016/j.apcatb.2019.117960)；[原始图5](https://ars.els-cdn.com/content/image/1-s2.0-S0926337319307064-gr5_lrg.jpg)。取下方面板 (b)。公开转载同一实验的图注注明噻虫嗪 2 × 10⁻⁵ mol/L、pH 7，见[图12(f)及其图注](https://www.mdpi.com/2073-4360/13/15/2568)。按[PubChem 的摩尔质量](https://pubchem.ncbi.nlm.nih.gov/compound/Thiamethoxam) 291.72 g/mol，名义初始浓度为 5.8344 mg/L。另两条对照曲线为空白和裸 LMPET。
- Farooq R, Singh M, Singh J, Mehta A. Frontiers in Materials, 12 (2025), 1677416. [原论文（方法和图7）](https://doi.org/10.3389/fmats.2025.1677416)；[原始图7](https://www.frontiersin.org/files/Articles/1677416/xml-images/fmats-12-1677416-g007.webp)。研究方法报告以 5 ppm 噻虫嗪储备液取 10 mL 用于光照实验，未另报稀释倍数。按水溶液名义 5.00 mg/L 换算。面板 (d) 为带催化剂负载量实验，无催化剂空白另予排除。0.5 g/L 条件与原论文图 7(b) 的 Ni-PbS 条件可能重复，图 7(b) 不进入本组计数。
- Xie 等开发数据源、3.558 mg/L 的由来及图表重建见上级 README 和 baseline 目录。该值是原研究的情景目标。两项外部研究均由相应研究团队独立开展。数据是由图形数字化得到的平均/代表性工况曲线，不能从图中恢复单次平行实验；名义浓度折算也不等同于建议时刻重新进行的仪器测量。

## 运行

在整个复现包的顶层执行 `python reproduce.py`，将重算基线、扩展实验和外部验证。只跑新实验时，先执行 `python external/validate.py`，再执行 `python external/check.py`。第二步调用原有 baseline/turnaround.py，对 7 个工况 × 2 种锚点读法 × 3 个浓度标度 × 5 个检测回报延迟的 210 组结果逐项比对。210 为情景组合数。

从出版原图重新数字化是可选步骤。将两个图源分别按原始字节保存为 chen5.jpg 和 front_fig7.webp 到同一目录。运行 `python external/extract.py 原图目录 --overlays 本地核图目录`，脚本先核对 SHA-256，再根据轴刻度和标记颜色恢复读数，覆盖 data/curves.csv 并生成位置叠图；随后重复上段的两个检查命令。两个图文件的哈希见 extract.py 的 PANELS。为减少第三方出版图像再分发，仅附数值、脚本、来源链接和哈希。左侧线段法只需 60 min，其中某些其他时刻无法提取时保留 NaN。标记主体法以颜色掩膜的连续宽行定位中心。根据相邻入射线段的读数仅供输入敏感性复核。

## 可核对的结果

- data/curves.csv 给出工况、取样时刻、名义初始浓度、图面坐标、标记中心读数与左线读数。
- results/decisions.csv 保存两种策略的计划时刻、建议类别、参考浓度、时间差及超标状态。results/summary.csv 按研究来源、读图法、名义浓度标度与延迟分别汇总。
- 主结果筛选 read=ratio、c0_factor=1、delay_min=0。以同次检测直接判停作对照。未发建议按 180 min 计时。预测建议恰落 180 min 时单列，不能算作提前。以建议时刻的原图读数核验，直接判停、未建议及到上限仍超标与“预测提前误停”分别记录。
- Chen 的总计划时间为 420 对 420 min，3 个工况中有 1 个直接判停，没有提前预测建议。Farooq 为 720 对 480 min，4 个工况中有 3 个提前预测建议，1 个超标。该错误发生在 0.3 g/L、120 min，参考浓度约 3.707 mg/L。两个来源没有合并成一个误停率。±5% 标度是人为设定的敏感性情景。
