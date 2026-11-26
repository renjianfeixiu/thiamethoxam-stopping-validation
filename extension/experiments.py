"""Reproduce paired timing trials, a second-assay ablation and cross-source tests."""
from pathlib import Path
import math
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'
OUT.mkdir(exist_ok=True)
GRID=np.array([30,60,90,120,150,180])
LIMIT=3.558
SEED=20261126

def save(rows,name):
    frame=pd.DataFrame(rows)
    frame.to_csv(OUT/(name+'.csv'),index=False,float_format='%.8f')
    return frame

def kinetic(r,a,t):
    return np.clip(r,0,1)**(np.array(t)/a)

def predicted(g,a,method,t):
    r=g.loc[a];t=np.asarray(t)
    if method=='exp':return kinetic(r,a,t)
    if method=='hold':return np.full(len(t),r)
    if method=='linear':return np.maximum(0,1+(r-1)*t/a)
    if not np.isfinite(g.loc[30]) or not np.isfinite(r):return np.full(len(t),np.nan)
    first,second=np.log(np.maximum([g.loc[30],r],1e-8))
    if method=='two_recent':
        slope=max(0,(first-second)/30)
    elif method=='two_origin':
        slope=max(0,-(30*first+60*second)/(30**2+60**2))
    else:raise ValueError(method)
    return r*np.exp(-slope*(t-a))

def summary(dec,key):
    rows=[]
    for k,g in dec.groupby(key,dropna=False,sort=False):
        k=k if isinstance(k,tuple) else (k,)
        nforecast=g[g.reason=='forecast']
        rows.append({**dict(zip(key,k)),'n':len(g),'issued':g.stop.notna().sum(),
                     'direct':(g.reason=='direct').sum(),'forecast':len(nforecast),
                     'forecast_early':sum(nforecast.stop<g.horizon.iloc[0]),
                     'forecast_pass':sum(nforecast.outcome=='pass'),'forecast_fail':sum(nforecast.outcome=='fail'),
                     'forecast_unreadable':sum(nforecast.outcome=='unreadable'),
                     'missing_anchor':sum(g.reason=='missing_anchor'),'abstain':sum(g.reason=='abstain'),
                     'pass':sum(g.outcome=='pass'),'fail':sum(g.outcome=='fail'),
                     'unreadable':sum(g.outcome=='unreadable'),
                     'mean_scheduled_min_saved_vs_horizon':(g.horizon-g.clock).mean(),
                     'mean_scheduled_min_saved_vs_same_assay_only':(g.assay_only_clock-g.clock).mean()})
    return pd.DataFrame(rows)

def decide(reference,a,input_values,method,cutoff,delay=0,guard=0):
    horizon=reference.index.max()
    times=reference.index[(reference.index>=a+delay)].to_numpy()
    sample=input_values.loc[a]
    stop=np.nan;reason='abstain'
    if pd.isna(sample) or (method.startswith('two_') and pd.isna(input_values.loc[30])):reason='missing_anchor'
    elif len(times):
        if sample<=cutoff:stop=times[0];reason='direct'
        else:
            forecast=predicted(input_values,a,method,times)+guard
            idx=np.flatnonzero((times>a)&(forecast<=cutoff))
            if len(idx):stop=times[idx[0]];reason='forecast'
    actual=reference.loc[stop] if np.isfinite(stop) else np.nan
    outcome=reason if not np.isfinite(stop) else ('unreadable' if pd.isna(actual) else ('pass' if actual<=cutoff else 'fail'))
    assay_only_clock=times[0] if len(times) and sample<=cutoff else horizon
    return dict(stop=stop,reference_ratio=actual,reason=reason,outcome=outcome,horizon=horizon,
                clock=horizon if pd.isna(stop) else stop,assay_only_clock=assay_only_clock)

def err_summary(df,keys):
    rows=[]
    for key,g in df.groupby(keys,sort=False):
        valid=g.dropna(subset=['error'])
        by=valid.groupby('curve').error.agg(lambda v:v.abs().mean())
        byct=valid.groupby('curve').ct_error.agg(lambda v:v.abs().mean())
        rows.append({**dict(zip(keys,key if isinstance(key,tuple) else (key,))),'n_series':len(by),
                     'n_read':len(valid),'mae_ratio':by.mean(),'mae_ct':byct.mean(),
                     'rmse_ratio':np.sqrt((valid.error**2).mean())})
    return pd.DataFrame(rows)

def main():
    conditions=pd.read_csv(ROOT/'data/table_s1.csv').set_index('event')
    s8=pd.read_csv(ROOT/'data/s8_grid.csv')
    ref=s8.pivot(index='event',columns='time',values='figure_ratio')
    ref[180]=conditions.ratio180
    eligible=conditions.index[conditions.c0>LIMIT]
    shared=ref.index[ref[[30,60,90]].notna().all(axis=1)]
    coverage=[]; dec=[]; errs=[];pairs=[]
    for read in ('figure_ratio','marker_ratio','left_ratio'):
        inp=s8.pivot(index='event',columns='time',values=read)
        for ev,truth in ref.iterrows():
            c0=conditions.loc[ev,'c0']
            for a,method in [(30,'exp'),(60,'exp'),(90,'exp'),(60,'two_origin'),(60,'two_recent')]:
                inp_row=inp.loc[ev]
                future=GRID[GRID>a]
                pred=predicted(inp_row,a,method,future)
                for t,yhat in zip(future,pred):
                    e=yhat-truth.loc[t]
                    errs.append(dict(read=read,curve=ev,anchor=a,method=method,target=t,
                                     shared=ev in shared,pred_ratio=yhat,observed_ratio=truth[t],error=e,ct_error=c0*e))
                if ev in eligible:
                    for delay in (0,30,60):
                        d=decide(truth,a,inp_row,method,LIMIT/c0,delay=delay)
                        dec.append(dict(read=read,curve=ev,anchor=a,method=method,delay=delay,shared=ev in shared,**d))
    err=save(errs,'timing_forecasts'); decisions=save(dec,'timing_decisions')
    for scope,sub in [('all',err),('identical_anchor_cohort',err[err.shared])]:
        for target in [120,150,180]:
            g=err_summary(sub[sub.target==target],['read','anchor','method'])
            g['target']=target;g['scope']=scope;coverage.append(g)
    save(pd.concat(coverage),'timing_metrics')
    all_counts=[]
    for scope,sub in [('all',decisions),('identical_anchor_cohort',decisions[decisions.shared])]:
        agg=summary(sub,['read','anchor','method','delay'])
        agg['scope']=scope;all_counts.append(agg)
    save(pd.concat(all_counts),'timing_counts')
    base=err[(err.read=='figure_ratio')&(err.target==180)&err.shared]
    for (a,m),(b,n) in [((30,'exp'),(60,'exp')),((60,'exp'),(90,'exp')),((60,'exp'),(60,'two_origin')),((60,'exp'),(60,'two_recent'))]:
        e1=base[(base.anchor==a)&(base.method==m)].set_index('curve').ct_error.abs()
        e2=base[(base.anchor==b)&(base.method==n)].set_index('curve').ct_error.abs()
        delta=(e1-e2).dropna(); rng=np.random.default_rng(SEED)
        q=np.quantile(rng.choice(delta.to_numpy(),(20000,len(delta)),replace=True).mean(axis=1),[.025,.975])
        pairs.append(dict(first=f'{m}_{a}',second=f'{n}_{b}',paired_n=len(delta),
                          second_lower=(delta>0).sum(),mean_reduction_mg_L=delta.mean(),
                          descriptive_bootstrap_low=q[0],descriptive_bootstrap_high=q[1]))
    save(pairs,'paired_endpoint_contrasts')
    # Safety margin calibrated solely from OTHER complete S8 batch curves.
    future=GRID[GRID>60]
    complete=ref.index[ref[[60,*future]].notna().all(axis=1)]
    scores={ev:max(0,np.max(ref.loc[ev,future].to_numpy()-kinetic(ref.loc[ev,60],60,future))) for ev in complete}
    def quantile(exclude=None):
        vals=sorted(v for k,v in scores.items() if k!=exclude)
        rank=math.ceil((len(vals)+1)*.90)
        return vals[rank-1] if rank<=len(vals) else np.inf
    save([dict(curve=ev,max_positive_ratio_residual=s) for ev,s in scores.items()],'calibration_scores')
    guarded=[]; audit=[]
    for ev in eligible:
        truth=ref.loc[ev];q=quantile(ev)
        for arm,value in [('unguarded',0),('other_batch_90',q)]:
            out=decide(truth,60,truth,'exp',LIMIT/conditions.loc[ev,'c0'],guard=value)
            guarded.append(dict(curve=ev,guard_rule=arm,buffer_ratio=value,n_cal=len(complete)-(ev in complete),**out))
        if ev in complete:
            audit.append(dict(curve=ev,loo_buffer=q,max_excess=scores[ev],covered=(scores[ev]<=q)))
    gg=save(guarded,'guard_decisions')
    save(summary(gg,['guard_rule']),'guard_counts')
    save(audit,'guard_coverage_complete_treatment')
    # Apply the frozen ratio rule at a=60 to OTHER figures without fitting on them.
    ext=pd.read_csv(ROOT/'data/transfer_curves.csv')
    ext['source']=ext.panel.map({'6f':'Xie Fig6','6g':'Xie Fig6','s6':'Xie FigS6',
                                 'F7b':'Farooq Fig7b','M1a':'Arfanis Fig1a','R5':'Reactor Fig5-6','R6':'Reactor Fig5-6'})
    ext.loc[ext.series=='no catalyst','source']='Farooq blank'
    errors=[]; decisions=[]; exact_refs=[]
    for (source,panel,series),g in ext.groupby(['source','panel','series']):
        truth=g.set_index('time').ratio.sort_index()
        # Do not extrapolate the guard past 180 min. Endpoints at 210 retained for error diagnostics only.
        horizon=min(180,truth.index.max())
        decision_ref=truth.loc[[t for t in truth.index if 60<=t<=horizon and t%30==0]]
        for a in (30,60,90):
            for method in ('exp','linear','hold'):
                times=truth.index[truth.index>a]
                if len(times)==0:continue
                pred=predicted(truth,a,method,times)
                for t,p in zip(times,pred):
                    e=p-truth.loc[t]
                    errors.append(dict(source=source,panel=panel,curve=series,anchor=a,method=method,
                                 target=t,error=e,ct_error=np.nan,pred_ratio=p,reference_ratio=truth.loc[t]))
        for cutoff in (.3,.5,.7):
            for rule,buffer in [('none',0),('S8_frozen',quantile())]:
                z=decide(decision_ref,60,truth,'exp',cutoff,guard=buffer)
                decisions.append(dict(source=source,panel=panel,curve=series,relative_target=cutoff,
                                      rule=rule,buffer_ratio=buffer,**z))
        # Error of the selected, unmodified one-point rule on the common tasks.
    ep=save(errors,'transfer_forecasts')
    es=err_summary(ep,['source','anchor','method','target'])
    save(es.drop(columns='mae_ct'),'transfer_metrics')
    dd=save(decisions,'transfer_decisions')
    save(summary(dd,['source','relative_target','rule']),'transfer_decision_counts')
    assert not err.duplicated(['read','curve','anchor','method','target']).any()
    assert not ep.duplicated(['source','curve','anchor','method','target']).any()
    for f in (gg,dd):
        assert not ((f.stop.notna()) & (f.stop>f.horizon)).any()
    print('all events',len(ref),'all treatment',len(eligible),'shared',len(shared),
          'shared treatment',len(set(eligible)&set(shared)),'calibration',len(scores),'frozen q',quantile())
    print(save([dict(all_events=len(ref),treatment_events=len(eligible),shared_events=len(shared),
          shared_treatment=len(set(eligible)&set(shared)),complete_calibration_n=len(scores),frozen_q=quantile())],'design').to_string(index=False))
if __name__=='__main__':main()
