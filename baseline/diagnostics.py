"""Summarize batch-paired errors and matched retrospective stopping times."""
from analysis import *
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator

FIG=ROOT/"figures"
FIG.mkdir(exist_ok=True)
names={"pooled":"Pooled k","catalyst":"Catalyst k","c0_only":"Initial concentration","condition_k":"Conditional kinetics","forest":"Random forest"}
teal="#17786f"
red="#b34438"
grey="#55636e"

def figsave(fig,name):
    fig.tight_layout()
    fig.savefig(FIG/(name+".png"),dpi=210,bbox_inches="tight")
    fig.savefig(FIG/(name+".pdf"),bbox_inches="tight")
    plt.close(fig)

def main():
    d,f=read_data()
    ep=pd.read_csv(OUT/"endpoint_predictions.csv")
    ev=ep[ep.split=="event"].copy()
    ev["absolute_error"]=abs(ev.pred_ct-ev.observed_ct)
    errors=ev.pivot(index="event",columns="model",values="absolute_error")
    rng=np.random.default_rng(SEED)
    diff=errors.forest-errors.condition_k
    res=rng.choice(diff.to_numpy(),size=(20000,len(diff)),replace=True).mean(1)
    write(pd.DataFrame([dict(mean_difference=diff.mean(),median_difference=diff.median(),
            conditional_smaller=int((diff>0).sum()),forest_smaller=int((diff<0).sum()),
            fixed_error_p025=np.quantile(res,.025),fixed_error_p975=np.quantile(res,.975))]),"paired_endpoint_errors")
    frame=ev.merge(d[["event","light_level","catalyst"]],on="event")
    frame["c0_band"]=np.where(frame.c0<7,"C0<7","C0>=7")
    summaries=[]
    for col in ["light_level","catalyst","c0_band"]:
        for (key,model),g in frame.groupby([col,"model"]):
            summaries.append(dict(group=col,value=key,model=model,n=len(g),mae_ct=g.absolute_error.mean()))
    write(pd.DataFrame(summaries),"endpoint_subgroups")
    st=pd.read_csv(OUT/"stop_decisions.csv")
    sk=st[(st.model=="condition_k") & (st.c0>LIMIT)]
    match=sk[sk.margin==0].merge(sk[sk.margin==.5],on="event",suffixes=("_0","_05"))
    match=match[match.time_0.notna()&match.time_05.notna()].copy()
    match["extra_minutes"]=match.time_05-match.time_0
    write(match[["event","time_0","time_05","outcome_0","outcome_05","extra_minutes"]],"paired_margins")
    # Source agreement only at 180 min; never treat the same endpoint as a second label.
    cross=f[f.time==180].merge(d[["event","ratio180"]],on="event")
    cross["absolute_difference"]=(cross.figure_ratio-cross.ratio180).abs()
    write(cross,"figure_table_check")
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,
                         "axes.labelcolor":"#273942","axes.edgecolor":"#9fa9ae",
                         "font.family":"DejaVu Sans","pdf.fonttype":42})
    fig,axs=plt.subplots(1,2,figsize=(8.4,3.35),sharex=True,sharey=True)
    for ax,m in zip(axs,["condition_k","forest"]):
        g=ev[ev.model==m]
        ax.plot([0,1],[0,1],color="#98a1a7",lw=1,ls="--")
        ax.scatter(g.ratio,g.pred_ratio,s=24,color=teal,alpha=.83,zorder=2)
        p=g[g.event==24].iloc[0]
        ax.scatter(p.ratio,p.pred_ratio,color=red,s=33,zorder=3)
        ax.annotate("E24",(p.ratio,p.pred_ratio),xytext=(4,-11),textcoords="offset points",color=red)
        ax.set(xlim=(0,1),ylim=(0,1),xlabel="Observed fraction remaining",title=names[m],aspect="equal")
    axs[0].set_ylabel("Held-out prediction")
    figsave(fig,"endpoint_predictions")
    fig,ax=plt.subplots(figsize=(7.6,3.15))
    for event,col,mark in [(12,teal,"o"),(24,red,"s"),(26,"#52668f","D")]:
        row=d[d.event==event].iloc[0]
        g=f[(f.event==event)&(f.time<180)].dropna(subset=["figure_ratio"])
        x=[0,*g.time.tolist(),180]; y=[row.c0,*(row.c0*g.figure_ratio).tolist(),row.ct180]
        ax.plot(x,y,color=col,lw=1.5,label="E%d"%event)
        ax.scatter(x[1:-1],y[1:-1],color="white",edgecolor=col,marker=mark,zorder=3)
        ax.scatter([0,180],[y[0],y[-1]],color=col,marker=mark,zorder=4)
    ax.axhline(LIMIT,color="#555",ls="--",lw=1,label="Study target 3.558 mg/L")
    ax.set(xlabel="Time after lamp-on (min)",ylabel="THM concentration (mg/L)",xticks=[0,60,90,120,150,180])
    ax.legend(frameon=False,ncol=2,fontsize=8,loc="best")
    figsave(fig,"nearby_batches")
    s=pd.read_csv(OUT/"stop_counts.csv");g=s[s.model=="condition_k"].sort_values("margin")
    fig,ax=plt.subplots(figsize=(7.4,2.6))
    left=np.zeros(len(g))
    for key,c,label in [("confirmed_pass",teal,"Met target"),("confirmed_fail",red,"Above target"),("abstain","#d8dfe0","No stop proposed")]:
        vals=g[key].to_numpy()
        ax.barh(range(len(g)),vals,left=left,color=c,height=.55,label=label)
        for i,v in enumerate(vals):
            if v: ax.text(left[i]+v/2,i,str(v),color="#243946" if key=="abstain" else "white",ha="center",va="center",size=9)
        left+=vals
    ax.set_yticks(range(len(g)),[f"{x:.2f} mg/L" for x in g.margin]); ax.invert_yaxis()
    ax.set_ylabel("Extra margin"); ax.set_xlabel("Treatment-requiring batches")
    ax.legend(frameon=False,fontsize=8,ncol=3,loc="lower center",bbox_to_anchor=(.5,1.0))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    figsave(fig,"stop_outcomes")
    early=pd.read_csv(OUT/"early_sample_predictions.csv")
    sample=early[early.event==24].copy()
    fig,ax=plt.subplots(figsize=(7.6,3.0))
    ref=sample[sample.model=="early_kinetics"].sort_values("time")
    for method,label,c,sty in [("observed60_then_prior","Continue prior forecast",grey,"--"),
                                    ("observed60_anchor_prior_k","Reset level; retain rate","#8c782e",":"),
                                    ("early_kinetics","Update rate from 60 min",teal,"-")]:
        g=sample[sample.model==method].sort_values("time")
        if method=="observed60_then_prior": g=g[g.time>60]
        ax.plot(g.time,g.pred_ct,label=label,color=c,ls=sty,lw=1.7)
    ax.scatter(ref.time,ref.observed_ct,marker="o",s=27,facecolor="white",edgecolor="#1a1a1a",zorder=4,label="Published reference")
    ax.axhline(LIMIT,color=red,ls="--",lw=.9,label="Study target")
    ax.set(xlabel="Time after lamp-on (min)",ylabel="THM concentration (mg/L)",xticks=TIMES)
    ax.legend(frameon=False,ncol=2,fontsize=8,loc="lower left")
    figsave(fig,"e24_one_sample")
    print("Matched times",len(match),match.extra_minutes.mean(),match.extra_minutes.median())
    print("Shape",pd.read_csv(OUT/"within_batch_shape_diagnostic.csv").groupby('event').abs_shape_error.mean().mean())
    print("Figure endpoint check",cross.absolute_difference.describe())
    print(pd.read_csv(OUT/"paired_endpoint_errors.csv").to_string(index=False))
if __name__=="__main__":
    main()
