# One-assay stopping decisions in thiamethoxam photocatalysis

This repository contains the digitized concentration data, analysis scripts and retrospective results for a study of batch planning and one-assay stopping decisions. The source studies performed the reactor experiments. The analysis uses published measurements and plotted points. The target of 3.558 mg/L comes from the development study's ECOSAR scenario.

## Reproduce

Tested with Python 3.13.9. Install the pinned dependencies in an isolated environment, then run from the repository root.

```sh
python -m pip install -r requirements.txt
python reproduce.py
python publication/figures.py
```

The first script refits the models, replays the schedules, recalculates sensitivity analyses and runs the verification checks. The second draws the two manuscript figures with Matplotlib. Reruns overwrite the derived results and figures. The bundled CSV inputs support offline analysis after dependency installation.

| Question | Inputs and outputs |
| --- | --- |
| How accurately can the endpoint be forecast before irradiation? | `baseline/data/table_s1.csv`, `baseline/results/endpoint_predictions.csv`, `baseline/results/endpoint_metrics.csv` |
| How much does extrapolation add to the same assay? | `baseline/results/turnaround_counts.csv`, `baseline/results/turnaround_decisions.csv` |
| Does the digitization change a decision? | `baseline/results/turnaround_anchor_inputs.csv`, E24 in the two files above |
| Does assay timing or kinetic form matter? | `extension/results/timing_metrics.csv`, `extension/results/timing_counts.csv`, `extension/results/transfer_metrics.csv` |
| Does the rule transfer without fitting? | `external/data/curves.csv`, `external/results/decisions.csv`, `external/results/summary.csv` |

## Independent sources and unit of analysis

The 32 development batches come from [Xie et al., Table S1 and Fig. S8](https://doi.org/10.1016/j.cej.2024.156029). In the external threshold analysis, the decision uses only the nominal initial concentration and the 60 min reading. The later observations serve as references. The seven external observations are seven condition curves from two publications, with shared conditions within each panel. They do not represent seven independent labs or seven sets of instrument replicates.

- [Chen et al., Fig. 5(b)](https://doi.org/10.1016/j.apcatb.2019.117960) provides three polyester-supported catalysts. The nominal starting concentration is 5.8344 mg/L. [The reproduced caption in Fig. 12(f)](https://doi.org/10.3390/polym13152568) gives the concentration in mol/L.
- [Farooq et al., Fig. 7(d)](https://doi.org/10.3389/fmats.2025.1677416) provides four Ni-doped PbS loadings at a nominal starting concentration of 5.00 mg/L. Fig. 7(b) of the same study supports a separate relative-error analysis. Panels are not pooled.
- Other normalized-trajectory comparisons use [Arfanis et al.](https://doi.org/10.3390/catal13030516) and [Emmanuel and Dilek](https://doi.org/10.11648/j.ajche.20241203.14). Source extraction details, image links, calibration methods and hashes appear in the component directories. Publisher images are obtained from their original sources.

With immediate results, extrapolation saved 420 scheduled minutes across 30 initially above-target development batches relative to using the same assay alone. In the external panels it yielded no increment in Chen's three conditions and 240 scheduled minutes across Farooq's four. One of three early forecast stops in Farooq's panel exceeded the target. Its 60 min contribution is included in the 240 min. The values describe replayed schedules. Energy use, the consequences of retreatment and prospective operating safety were not measured.

The research computations and preparation of the plotting scripts used OpenAI Codex. Verification scripts compare decisions, check source extraction assumptions, and reproduce tabulated totals. These checks address computational consistency; source measurements remain those of the cited publications.

[中文复现说明](README.zh-CN.md) · [External digitization and reference coordinates](external/README.md)
