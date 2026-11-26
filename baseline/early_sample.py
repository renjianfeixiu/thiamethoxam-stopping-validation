"""Additional-information study: reveal a single t=60 chart measurement."""
from analysis import *
import numpy as np
import pandas as pd

def main():
    d,f=read_data()
    g=pd.read_csv(OUT/"curve_predictions.csv").query("model=='condition_k'").copy()
    sixty=f[f.time==60].set_index("event").figure_ratio
    assert len(sixty)==32 and sixty.notna().all()
    r60=g.event.map(sixty)
    # Both arms know exactly the same 60-min reading. Beyond 60 min, only
    # the updated arm uses that reading to set its future slope.
    one=g.copy()
    one["model"]="early_kinetics"
    one["pred_ratio"]=r60**(one.time/60)
    common=g.copy()
    common["model"]="observed60_then_prior"
    common.loc[common.time==60,"pred_ratio"]=r60[common.time==60]
    anchor=g.copy()
    anchor["model"]="observed60_anchor_prior_k"
    prior60=g[g.time==60].set_index("event").pred_ratio
    anchor["pred_ratio"]=r60*anchor.pred_ratio/anchor.event.map(prior60)
    linear=g.copy()
    linear["model"]="early_linear"
    linear["pred_ratio"]=np.maximum(0,1+(r60-1)*(linear.time/60))
    out=pd.concat([one,common,anchor,linear],ignore_index=True)
    out["pred_ct"]=out.c0*out.pred_ratio
    write(out,"early_sample_predictions")
    write(score_rows(out[out.time>60]),"early_sample_future_metrics")
    decisions_out=decisions(out,[0,.25,.5])
    write(decisions_out,"early_sample_decisions")
    write(decision_counts(decisions_out),"early_sample_counts")
    matched=out[(out.time==180)].copy()
    write(score_rows(matched),"early_sample_endpoint_metrics")
    # Explore +/- 0.01 in relative concentration of the 60-min read. Perturb the
    # model input only, then score future held-out measurements unchanged.
    sensitivity=[]
    for shift in [-.01,0,.01]:
        v=one.copy()
        v["pred_ratio"]=np.clip(v.event.map(sixty)+shift,.001,1.0)**(v.time/60)
        v["pred_ct"]=v.c0*v.pred_ratio
        counts=decision_counts(decisions(v,[0,.25,.5]))
        counts["shift_input_ratio"]=shift
        sensitivity.append(counts)
    write(pd.concat(sensitivity,ignore_index=True),"early_sample_input_sensitivity")

if __name__=="__main__": main()
