"""Replay delayed availability of a single assay against an assay-only control."""
from pathlib import Path
import numpy as np
import pandas as pd
from analysis import LIMIT, ROOT, OUT, TIMES, assemble, read_data, write, score_rows

DELAYS=(0,30,60,90,120)

def replay(curves, anchors):
    rows=[]
    for event,g in curves.groupby("event"):
        g=g.set_index("time")
        c0=g.c0.iloc[0]
        if c0<=LIMIT:
            continue
        for source,r60 in anchors.loc[event].items():
            for delay in DELAYS:
                available=TIMES[TIMES>=60+delay]
                for strategy in ("assay_only","rate_update"):
                    concentrations=c0*(r60**(available/60))
                    if c0*r60<=LIMIT:
                        stop=int(available[0])
                        reason="sample_already_met_target"
                    elif strategy=="rate_update" and np.any(concentrations<=LIMIT):
                        stop=int(available[concentrations<=LIMIT][0])
                        reason="forecast_from_above"
                    else:
                        stop=None
                        reason="no_stop"
                    measured=c0*g.loc[stop,"ratio"] if stop is not None else np.nan
                    if stop is None:
                        status="no_stop"
                    elif pd.isna(measured):
                        status="missing"
                    else:
                        status="pass" if measured<=LIMIT else "fail"
                    rows.append(dict(event=event,reading=source,strategy=strategy,delay=delay,
                                     r60=r60,stop_min=stop,reason=reason,outcome=status,
                                     reference_ct=measured,horizon_lamp_min=180 if stop is None else stop))
    return pd.DataFrame(rows)

def summarize(rows):
    records=[]
    for (source,strategy,delay),g in rows.groupby(["reading","strategy","delay"]):
        early=g[g.reason=="forecast_from_above"]
        records.append(dict(reading=source,strategy=strategy,delay=delay,n=len(g),
            sample_already_met_target=sum(g.reason=="sample_already_met_target"),
            forecast_n=len(early),forecast_pass=sum(early.outcome=="pass"),
            forecast_fail=sum(early.outcome=="fail"),forecast_missing=sum(early.outcome=="missing"),
            proposed=g.stop_min.notna().sum(),pass_n=sum(g.outcome=="pass"),
            fail_n=sum(g.outcome=="fail"),missing_n=sum(g.outcome=="missing"),
            no_stop=sum(g.outcome=="no_stop"),capped_minutes=g.horizon_lamp_min.sum(),
            avg_minutes_vs_180=(180-g.horizon_lamp_min).mean()))
    summary=pd.DataFrame(records)
    control=summary[summary.strategy=="assay_only"].set_index(["reading","delay"]).capped_minutes
    summary["avg_minutes_vs_assay_only"]=[(control.loc[(r.reading,r.delay)]-r.capped_minutes)/r.n
                                           for r in summary.itertuples()]
    return summary

def main():
    d,f=read_data()
    curves=assemble(d,f)
    anchors=f[f.time==60].set_index("event")[["figure_ratio"]].rename(columns={"figure_ratio":"line_marker"})
    alt=pd.read_csv(ROOT/"data/figure_s8_marker60.csv").set_index("event")
    anchors["marker_only"]=alt.marker_ratio
    assert anchors.shape==(32,2) and anchors.notna().all().all()
    write(anchors.reset_index(),"turnaround_anchor_inputs")
    rows=replay(curves,anchors)
    summary=summarize(rows)
    write(rows,"turnaround_decisions")
    write(summary,"turnaround_counts")
    # Score only readings strictly after the assay. For each variant, also compute
    # error in the subset above target at t=60 under the primary reference reading.
    models=pd.read_csv(OUT/"early_sample_predictions.csv")
    base=models[models.model=="early_kinetics"].copy()
    base["model"]="early_marker_only"
    base["pred_ratio"]=base.event.map(anchors.marker_only)**(base.time/60)
    base["pred_ct"]=base.c0*base.pred_ratio
    pred=pd.concat([models,base],ignore_index=True)
    future=pred[pred.time>60].copy()
    cohort=d[d.c0>LIMIT].set_index("event").join(anchors)
    ids=cohort.index[cohort.c0*cohort.line_marker>LIMIT]
    assert len(ids)==22
    scores=[]
    pairs=[]
    rng=np.random.default_rng(20261126)
    for name,v in [("all_32",future),("still_above_at_60_n22",future[future.event.isin(ids)])]:
        result=score_rows(v)
        result["population"]=name
        scores.append(result)
        batch_errors=v.assign(err=(v.pred_ct-v.observed_ct).abs()).groupby(["event","model"]).err.mean().unstack()
        for name_other in ("observed60_then_prior","observed60_anchor_prior_k","early_linear"):
            dif=(batch_errors[name_other]-batch_errors.early_kinetics).dropna()
            res=rng.choice(dif.to_numpy(),(20000,len(dif)),replace=True).mean(axis=1)
            pairs.append(dict(population=name,compared_to=name_other,n_events=len(dif),n_lower=(dif>0).sum(),
                              n_equal=(dif==0).sum(),mean_mae_reduction=dif.mean(),
                              bootstrap_q025=np.quantile(res,.025),bootstrap_q975=np.quantile(res,.975)))
    write(pd.concat(scores,ignore_index=True),"future_scores_by_population")
    write(pd.DataFrame(pairs),"paired_future_errors")
    # Checks independent of aggregated tables.
    assert rows.groupby(["reading","strategy","delay"]).size().eq(30).all()
    available=rows.dropna(subset="stop_min")
    assert available.stop_min.ge(60+available.delay).all()
    assert rows.pivot(index=["reading","event","strategy"],columns="delay",values="horizon_lamp_min").diff(axis=1).iloc[:,1:].ge(0).all().all()
    assert (summary.proposed==summary.pass_n+summary.fail_n+summary.missing_n).all()
    primary=summary[(summary.reading=="line_marker")&(summary.strategy=="rate_update")&(summary.delay==0)].iloc[0]
    assert (primary.proposed,primary.forecast_n,primary.forecast_pass)==(16,8,8)
    print(summary.to_string(index=False))
if __name__=="__main__": main()
