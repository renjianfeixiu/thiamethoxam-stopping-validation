"""Check source accounting, out-of-fold multiplicity and decision denominators."""
from analysis import *

def main():
    d,f=read_data()
    assert len(d)==d.event.nunique()==32 and d.ct180.le(LIMIT).sum()==18
    assert d.c0.le(LIMIT).sum()==2
    assert np.allclose(d.ct180,d.c0*(1-d.de_pct/100))
    assert not f.duplicated(["event","time"]).any()
    assert f[f.time<180].figure_ratio.notna().sum()==120
    cross=f[f.time==180].merge(d,on="event").dropna(subset=["figure_ratio"])
    assert len(cross)==28 and (cross.figure_ratio-cross.ratio180).abs().max()<.004
    ep=pd.read_csv(OUT/"endpoint_predictions.csv")
    assert len(ep)==32*len(MODELS)*2 and not ep.duplicated(["event","model","split"]).any()
    grid=pd.read_csv(OUT/"curve_predictions.csv")
    assert len(grid)==32*len(MODELS)*len(TIMES)
    assert (grid.fold_id==grid.event).all()
    counts=pd.read_csv(OUT/"stop_counts.csv")
    assert (counts.confirmed_pass+counts.confirmed_fail+counts.unreadable==counts.recommended).all()
    assert (counts.abstain+counts.recommended==30).all()
    z=counts[(counts.model=="condition_k")&(counts.margin==.5)].iloc[0]
    assert (z.recommended,z.confirmed_pass,z.confirmed_fail)==(15,14,1)
    early=pd.read_csv(OUT/"early_sample_counts.csv")
    item=early[(early.model=="early_kinetics")&(early.margin==0)].iloc[0]
    assert (item.recommended,item.confirmed_pass,item.confirmed_fail,item.unreadable)==(16,16,0,0)
    dec=pd.read_csv(OUT/"early_sample_decisions.csv")
    later=dec[(dec.model=="early_kinetics")&(dec.margin==0)&(dec.time>60)&(dec.c0>LIMIT)]
    assert len(later)==8 and later.outcome.eq("pass").all()
    assert pd.read_csv(OUT/"early_sample_future_metrics.csv").n_read.eq(120).all()
    print("Validation passed: source data, folds, tested predictions, outcome partitions.")
if __name__=="__main__": main()
