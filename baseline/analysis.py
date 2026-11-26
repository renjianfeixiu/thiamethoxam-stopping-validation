"""Run batch-held-out endpoint prediction and retrospective stopping decisions."""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent
OUT = ROOT/"results"
OUT.mkdir(exist_ok=True)
LIMIT=3.558
TIMES=np.array([60,90,120,150,180])
SEED=20261126
MODELS=["pooled","catalyst","c0_only","condition_k","forest"]

def read_data():
    d=pd.read_csv(ROOT/"data/table_s1.csv")
    f=pd.read_csv(ROOT/"data/figure_s8.csv")
    return d,f

def make_features(d,kinetic=False):
    acidity=d.pH.to_numpy()-7
    dose=np.log(d.dose_mg.to_numpy()/50)
    columns=[np.log(d.c0),acidity, dose, d.catalyst.eq("CN/Bi2O3"),
             d.light_level.eq(331), d.light_level.eq(520)]
    if kinetic:
        columns.extend([acidity**2,dose**2])
    return np.column_stack(columns)

def recipe(row):
    def near(v,vals): return min(vals,key=lambda z:abs(z-v))
    return (row.catalyst, round(row.pH), near(row.c0,[5,10,20,50]),
            near(row.dose_mg,[25,50,100]),row.light_level)

def assemble(d,f):
    chart=f[f.time<180][["event","time","figure_ratio"]].rename(columns={"figure_ratio":"ratio"})
    end=d[["event","ratio180"]].assign(time=180).rename(columns={"ratio180":"ratio"})
    return pd.concat([chart,end],ignore_index=True).merge(d,on="event").sort_values(["event","time"])

def predict(train,test,model,task):
    # Rows in train contain labels only for fitting events. Each batch has
    # equal total weight in the forest, irrespective of visible chart markers.
    train=train.dropna(subset=["ratio"])
    n=train.groupby("event").size()
    if model=="forest":
        a,b=make_features(train),make_features(test)
        if task=="curve":
            a=np.column_stack([a,train.time])
            b=np.column_stack([b,test.time])
        rf=RandomForestRegressor(n_estimators=240,min_samples_leaf=3,max_features=1.0,
                                 random_state=SEED,n_jobs=1)
        rf.fit(a,train.ratio,sample_weight=train.event.map(1/n))
        return rf.predict(b)
    # Zero-intercept log-linear kinetic fit separately in each training batch.
    s=train.assign(num=-train.time*np.log(train.ratio),den=train.time**2)
    g=s.groupby("event")[["num","den"]].sum()
    k=g.num/g.den
    batches=train.drop_duplicates("event").set_index("event").loc[k.index]
    if model in ("condition_k","c0_only"):
        ridge=make_pipeline(StandardScaler(),Ridge(alpha=4.0))
        x=make_features(batches,True)
        xt=make_features(test,True)
        if model=="c0_only": x,xt=x[:,:1],xt[:,:1]
        ridge.fit(x,np.log(k))
        predicted=np.exp(ridge.predict(xt))
    elif model=="catalyst":
        predicted=test.catalyst.map(k.groupby(batches.catalyst).median()).to_numpy()
    else:
        predicted=np.full(len(test),k.median())
    return np.exp(-predicted*test.time.to_numpy())

def held_out(data,task,split,models=MODELS):
    result=[]
    for group,te in data.groupby(split):
        tr=data[data[split]!=group]
        assert set(te.event).isdisjoint(tr.event)
        for model in models:
            sub=te[["event","time","c0","ratio"]].copy()
            sub["pred_ratio"]=predict(tr,te,model,task)
            sub["pred_ct"]=sub.c0*sub.pred_ratio
            sub["observed_ct"]=sub.c0*sub.ratio
            sub["model"]=model
            sub["split"]=split
            sub["fold_id"]=group
            result.append(sub)
    return pd.concat(result,ignore_index=True)

def score_rows(pred):
    rows=[]
    for (split,model),g in pred.dropna(subset=["ratio"]).groupby(["split","model"]):
        err=g.pred_ct-g.observed_ct
        re=g.pred_ratio-g.ratio
        per_batch=err.abs().groupby(g.event).mean()
        r2=lambda e,y: 1-(e**2).sum()/((y-y.mean())**2).sum()
        rows.append(dict(split=split,model=model,n_events=g.event.nunique(),n_read=len(g),
                         mae_ct=per_batch.mean(),rmse_ct=np.sqrt((err**2).mean()),r2_ct=r2(err,g.observed_ct),
                         mae_ratio=re.abs().groupby(g.event).mean().mean(),r2_ratio=r2(re,g.ratio)))
    return pd.DataFrame(rows)

def decisions(pred,margins):
    records=[]
    for (split,model,event),g in pred.groupby(["split","model","event"]):
        g=g.sort_values("time")
        c0=g.c0.iloc[0]
        for margin in margins:
            choices=g[g.pred_ct<=LIMIT-margin]
            if c0<=LIMIT:
                t,predct,actual=0,c0,c0
            elif choices.empty:
                t,predct,actual=np.nan,np.nan,np.nan
            else:
                first=choices.iloc[0]
                t,predct,actual=first.time,first.pred_ct,first.observed_ct
            # Chart-reading sensitivity near threshold; 180-min endpoint as tabulated.
            obs="pass" if actual<=LIMIT else "fail"
            if np.isnan(t): obs="abstain"
            elif np.isnan(actual): obs="unreadable"
            borderline=np.isfinite(t) and 0<t<180 and abs(actual-LIMIT)<=.01*c0
            records.append(dict(split=split,model=model,event=event,margin=margin,c0=c0,time=t,
                                predicted_ct=predct,observed_ct=actual,outcome=obs,borderline=borderline))
    return pd.DataFrame(records)

def decision_counts(rows):
    out=[]
    for (split,model,margin),g in rows[rows.c0>LIMIT].groupby(["split","model","margin"]):
        valid=g[g.outcome.isin(["pass","fail"])]
        passed=g[g.outcome.eq("pass")]
        out.append(dict(split=split,model=model,margin=margin,n=len(g),
                        recommended=int(g.time.notna().sum()),confirmed_pass=len(passed),
                        confirmed_fail=sum(g.outcome.eq("fail")),unreadable=sum(g.outcome.eq("unreadable")),
                        abstain=sum(g.outcome.eq("abstain")),borderline=int(g.borderline.sum()),
                        avg_min_saved_verified_pass=(180-passed.time).mean(),
                        false_events=" ".join("E%d"%i for i in valid.loc[valid.outcome=="fail","event"])))
    return pd.DataFrame(out)

def write(df,name):
    df.to_csv(OUT/(name+".csv"),index=False,float_format="%.7f")

def main():
    d,f=read_data()
    curves=assemble(d,f)
    curves["nominal_recipe"]=curves.apply(lambda r: str(recipe(r)),axis=1)
    end=curves[curves.time==180].copy()
    endpoints=pd.concat([held_out(end,"endpoint",j) for j in ["event","nominal_recipe"]])
    write(endpoints,"endpoint_predictions")
    write(score_rows(endpoints),"endpoint_metrics")
    e=decisions(endpoints,[0])
    write(e,"endpoint_decisions")
    write(decision_counts(e),"endpoint_decision_counts")
    grid=held_out(curves,"curve","event")
    write(grid,"curve_predictions")
    write(score_rows(grid),"curve_metrics")
    s=decisions(grid,[0,.25,.50])
    write(s,"stop_decisions")
    write(decision_counts(s),"stop_counts")
    # Matched support audit: at 50 mL, 0.5-1 g/L targets are 25-50 mg.
    # 0.5 mg allowance explicitly accounts for weighed vs nominal amounts.
    supported=d[(d.pH.between(6,8)) & (abs(d.c0-5)<=1) & d.dose_mg.between(24.5,50.5)].copy()
    supported["material_usd_per_batch"]=supported.dose_mg/1000*supported.catalyst.map({"CN":.189,"CN/Bi2O3":.182})
    write(supported,"supported_recipes")
    # Retrospective time optimization with observed external light and fixed materials.
    # Each source event is evaluated using the model that did not train on its curve.
    proposals=grid[(grid.event.isin(supported.event)) & grid.model.eq("condition_k")].copy()
    proposals=proposals.merge(supported[["event","light_level","material_usd_per_batch"]],on="event")
    proposals["predicted_feasible"]=proposals.pred_ct<=LIMIT
    write(proposals,"supported_recipe_time_grid")
    # Fixed 180-min benchmark for the identical treatment-requiring population.
    base=end[end.c0>LIMIT]
    write(pd.DataFrame([dict(n=len(base),stop_min=180,pass_n=sum(base.ct180<=LIMIT),
                                 fail_n=sum(base.ct180>LIMIT))]),"fixed_time_reference")
    # Blocked ablation: exclude the single source dosage typo and retrain.
    ab=held_out(end[end.event!=25],"endpoint","event")
    write(score_rows(ab),"endpoint_excluding_e25")
    # Best-case model-free diagnostic: if the withheld batch's endpoint were given,
    # does interpolating that endpoint by first-order kinetics reproduce its interior?
    interior=curves[(curves.time<180)&curves.ratio.notna()].copy()
    interior["shape_pred"]=interior.ratio180**(interior.time/180)
    interior["abs_shape_error"]=(interior.shape_pred-interior.ratio).abs()
    write(interior[["event","time","ratio","shape_pred","abs_shape_error"]],"within_batch_shape_diagnostic")
    print("Independent nominal folds:",end.nominal_recipe.nunique())
    for filename in ["endpoint_metrics","endpoint_decision_counts","curve_metrics","stop_counts","supported_recipe_time_grid","fixed_time_reference"]:
        print(filename,pd.read_csv(OUT/(filename+".csv")).round(3).to_string(index=False))

if __name__=="__main__": main()
