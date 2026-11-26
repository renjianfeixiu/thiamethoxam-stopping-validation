"""Draw the manuscript figures from the published decision rules and extracted data."""
from pathlib import Path
import csv
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "publication" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})
ink, blue, teal = "#1e2d37", "#305873", "#146959"

def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=420, facecolor="white")
    plt.close(fig)

def workflow():
    fig, ax = plt.subplots(figsize=(7.2, 4.05))
    fig.subplots_adjust(left=.012, right=.988, bottom=.01, top=.99)
    ax.set(xlim=(0, 100), ylim=(0, 100))
    ax.axis("off")
    def node(x, y, w, h, title, body, shade, color):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.2,rounding_size=1.2",
                                    linewidth=.8, edgecolor=color, facecolor=shade))
        ax.text(x+2.2,y+h-3,title,va="top",color=color,fontsize=9.6,fontweight="bold")
        ax.text(x+2.2,y+h-10,body,va="top",color=ink,fontsize=8.7,linespacing=1.4)
    def arrow(a,b):
        ax.annotate("", xy=b, xytext=a, arrowprops={"arrowstyle":"-|>","color":"#536874","lw":1})
    node(1,69,43,29,"A   Before irradiation","Recipe and initial concentration\n→ forecast the 180 min concentration\n→ flag batches for closer monitoring","#f1f5f9",blue)
    node(53,69,46,29,"B   After 60 min of irradiation","Read the concentration in this batch.\nBoth policies receive the same result\nafter the same reporting delay.","#eef7f4",teal)
    arrow((44.5,84),(52.6,84))
    arrow((64,68.8),(58,61.8))
    arrow((86,68.8),(86,61.8))
    node(38,31,29,30,"Assay only","Target met: act when\nthe result is available.\nOtherwise, retain 180 min.","#f4f7f6",teal)
    node(70,31,29,30,"Assay + extrapolation","If still above target, use\nthe first predicted crossing\non the grid, or retain 180 min.","#f4f7f6",teal)
    node(1,31,33,30,"Validate the endpoint","32 whole-batch holdouts\nand 29 recipe-group holdouts.\nScore error at 180 min.","#f4f6f9",blue)
    arrow((16,68.8),(16,61.8))
    arrow((53,30.7),(53,25.7))
    arrow((85,30.7),(85,25.7))
    node(1,1,98,24,"Evaluate the stopping rule","Replay 30 initially above-target development batches. Transfer the same rule to 3 Chen and 4 Farooq\nconditions. Compare scheduled minutes and check reference concentrations at the proposed times.","#f8f8f6",ink)
    save(fig,"Fig1_workflow")

def sensitivity():
    def event_row(file):
        with file.open() as f:
            return next(row for row in csv.DictReader(f) if row["event"]=="24")
    source = event_row(ROOT / "baseline/data/table_s1.csv")
    anchor = event_row(ROOT / "baseline/results/turnaround_anchor_inputs.csv")
    c0 = float(source["c0"])
    x = np.linspace(.886, .900, 130)
    fig, ax = plt.subplots(figsize=(6.0,2.85))
    fig.subplots_adjust(left=.14,right=.965,bottom=.23,top=.945)
    ax.plot(x,c0*x**3,color=blue,lw=1.7)
    ax.axhline(3.558,color=ink,ls="--",lw=1)
    ax.text(.8999,3.554,"Target  3.558",va="top",ha="right",fontsize=8.2)
    for name,marker,label,offset in [("marker_only","s","Marker only",(7,-14)),
                                     ("line_marker","o","Marker + line",(-52,14))]:
        r=float(anchor[name])
        ax.scatter(r,c0*r**3,s=34,marker=marker,facecolors="white",edgecolors=blue,zorder=4)
        ax.annotate(label,(r,c0*r**3),xytext=offset,textcoords="offset points",fontsize=8.7)
    ax.set(xlim=(x.min(),x.max()),ylim=(3.475,3.66),xticks=[.886,.889,.892,.895,.898],
           ylabel="Predicted endpoint (mg/L)",
           xlabel="Fraction remaining at 60 min")
    ax.tick_params(labelsize=8.5)
    save(fig,"Fig2_digitization")

if __name__ == "__main__":
    workflow()
    sensitivity()
