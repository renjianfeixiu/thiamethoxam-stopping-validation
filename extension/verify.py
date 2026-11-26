from pathlib import Path
import pandas as pd
import numpy as np

root=Path(__file__).resolve().parent
def result(name):
    return pd.read_csv(root/'results'/(name+'.csv'))
design=result('design').iloc[0]
assert tuple(design[['all_events','treatment_events','shared_events','shared_treatment']])==(32,30,29,27)
m=result('timing_metrics').query('scope=="identical_anchor_cohort" and method=="exp" and read=="figure_ratio" and target==180').sort_values('anchor')
assert (m.n_series==29).all() and np.allclose(m.mae_ct,[.901933,.569408,.360176],atol=.000001)
r=result('timing_counts').query('scope=="identical_anchor_cohort" and method=="exp" and read=="figure_ratio" and delay==0').sort_values('anchor')
assert (r.n==27).all() and list(r.forecast)==[14,8,5] and list(r.forecast_fail)==[5,0,0]
a=result('timing_counts').query('scope=="all" and anchor==60 and method=="exp" and read=="marker_ratio" and delay==0').iloc[0]
assert a.forecast==8 and a.forecast_fail==1
m=result('transfer_metrics').query('anchor==60 and method=="exp"')
assert dict(m[m.target==180].set_index('source').n_series)=={'Farooq Fig7b':2,'Farooq blank':1,'Reactor Fig5-6':4,'Xie Fig6':8,'Xie FigS6':2}
assert m[(m.source=='Arfanis Fig1a')&(m.target==90)].n_series.iloc[0]==1
print('Verified manuscript denominators, decision outcomes and external coverage.')
