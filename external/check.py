"""Compare all decisions against the archived pre-external implementation."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from validate import ROOT

def main():
    sys.path.insert(0,str(ROOT.parent/'baseline'))
    from turnaround import replay
    points=pd.read_csv(ROOT/'data/curves.csv')
    out=pd.read_csv(ROOT/'results/decisions.csv')
    assert len(points)==42 and not points.duplicated(['panel','condition','time_min']).any()
    assert len(out)==210
    for (panel,condition,factor),rr in out.groupby(['panel','condition','c0_factor']):
        g=points[(points.panel==panel)&(points.condition==condition)]
        ref=g.assign(event=1,time=g.time_min,c0=g.nominal_c0_mg_l*factor)[['event','time','c0','ratio']]
        anchor=g[g.time_min==60].copy().set_index(pd.Index([1]))[['ratio','incoming_ratio']]
        rr0=replay(ref,anchor)
        for item in rr.itertuples():
            pair=rr0[(rr0.reading==item.read)&(rr0.delay==item.delay_min)].set_index('strategy')
            model=pair.loc['rate_update']
            assert item.assay_only_time==pair.loc['assay_only'].horizon_lamp_min
            assert item.policy_time==model.horizon_lamp_min
            assert item.forecast_false==int(model.reason=='forecast_from_above' and model.outcome=='fail')
            assert item.reason=={'forecast_from_above':'forecast','no_stop':'abstain','sample_already_met_target':'direct'}[model.reason]
    q=out[(out.read=='ratio')&(out.c0_factor==1)&(out.delay_min==0)]
    c=q[q.panel=='Chen5b']; f=q[q.panel=='Farooq7d']
    assert (len(c),sum(c.reason=='direct'),sum(c.reason=='forecast'),c.minutes_saved.sum())==(3,1,0,0)
    assert (len(f),sum(f.reason=='forecast'),f.premature_forecast_false.sum(),f.minutes_saved.sum())==(4,3,1,240)
    assert np.isclose(f.loc[f.premature_forecast_false==1,'reference_at_policy_mg_l'].iloc[0],3.706915)
    print('42 digitized points, 210 decisions cross-checked against the archived replay.')
if __name__=='__main__':main()
