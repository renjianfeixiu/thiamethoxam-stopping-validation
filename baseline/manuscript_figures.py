"""Additional paper figures from the batch-level calculations."""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from analysis import ROOT, OUT, LIMIT, read_data
from diagnostics import figsave, teal, red, grey

def main():
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"pdf.fonttype":42,
                         "axes.spines.top":False,"axes.spines.right":False})
    d,_=read_data()
    c0=d.loc[d.event==24,"c0"].item()
    r=pd.read_csv(OUT/"turnaround_anchor_inputs.csv").set_index("event").loc[24]
    xs=np.linspace(.885,.9,81)
    fig,ax=plt.subplots(figsize=(7.2,2.9))
    ax.plot(xs,c0*xs**3,color=grey,lw=1.5)
    ax.axhline(LIMIT,ls="--",color="#8c7c63",label="Stopping target 3.558 mg/L")
    for field,col,name,offset in [("marker_only",red,"Marker only",(12,-17)),
                 ("line_marker",teal,"Marker and line",(-84,14))]:
        x=r[field];y=c0*x**3
        ax.scatter([x],[y],color=col,s=32,zorder=3)
        ax.annotate(f"{name}: {y:.3f}",(x,y),xytext=offset,textcoords="offset points",color=col,fontsize=9)
    ax.set(xlabel="Fraction remaining at 60 min",ylabel="Predicted concentration\nat 180 min (mg/L)")
    ax.legend(frameon=False,fontsize=8,loc="lower right")
    figsave(fig,"e24_reading_threshold")
    s=pd.read_csv(OUT/"turnaround_counts.csv").query("reading=='line_marker'")
    fig,ax=plt.subplots(figsize=(7.2,2.9))
    for strategy,label,col,marker in [("rate_update","Assay + batch rate",teal,"o"),("assay_only","Assay alone",grey,"s")]:
        g=s[s.strategy==strategy].sort_values("delay")
        ax.plot(g.delay,g.avg_minutes_vs_180,marker=marker,lw=1.7,color=col,label=label)
    ax.set(xticks=[0,30,60,90,120],xlabel="Time from sampling to result (min)",
           ylabel="Mean lamp time saved\nwithin 180 min (min/batch)",ylim=(-2,52))
    ax.legend(frameon=False,fontsize=9)
    figsave(fig,"assay_turnaround")
if __name__=="__main__":main()
