"""Apply the unchanged stopping rule and then unmask the held-out observations."""
from pathlib import Path
import csv
from collections import defaultdict
from statistics import mean
import math

ROOT=Path(__file__).resolve().parent
LIMIT=3.558
GRID=(60,90,120,150,180)
def read_groups():
    groups=defaultdict(dict)
    with (ROOT/'data/curves.csv').open() as src:
        for row in csv.DictReader(src):
            key=row['panel'],row['condition'],float(row['nominal_c0_mg_l'])
            t=int(row['time_min'])
            if t in groups[key]:raise ValueError('Duplicate point')
            groups[key][t]=row
    return groups

def schedules(c0,anchor,delay):
    valid=[t for t in GRID if t>=60+delay]
    if c0*anchor<=LIMIT:return 'direct',valid[0],valid[0]
    for t in valid:
        if t>60 and c0*min(1,max(0,anchor))**(t/60)<=LIMIT:return 'forecast',t,180
    return 'abstain',180,180

def decisions(groups):
    for (panel,series,c0),points in groups.items():
        assert set(points)>={*GRID}
        for delay in (0,30,60,90,120):
            for read in ('ratio','incoming_ratio'):
                for factor in (.95,1.,1.05):
                    adjusted=c0*factor
                    r=float(points[60][read])
                    reason,t,base=schedules(adjusted,r,delay)
                    ref=adjusted*float(points[t]['ratio'])
                    comparator_ref=adjusted*float(points[base]['ratio'])
                    yield dict(panel=panel,condition=series,read=read,c0_factor=factor,delay_min=delay,
                        c0_mg_l=round(adjusted,6),anchor_mg_l=round(adjusted*r,6),reason=reason,
                        assay_only_time=base,policy_time=t,minutes_saved=base-t,
                        reference_at_policy_mg_l=round(ref,6),reference_at_assay_only_time_mg_l=round(comparator_ref,6),
                        forecast_false=int(reason=='forecast' and ref>LIMIT),
                        premature_forecast_false=int(reason=='forecast' and t<180 and ref>LIMIT),
                        reference_above_at_cap=int(t==180 and ref>LIMIT))
def write(name,rows):
    with (ROOT/'results'/name).open('w',newline='') as out:
        w=csv.DictWriter(out,rows[0]);w.writeheader();w.writerows(rows)
def main():
    (ROOT/'results').mkdir(exist_ok=True)
    rows=list(decisions(read_groups()));write('decisions.csv',rows)
    by=defaultdict(list)
    for r in rows:by[(r['panel'],r['read'],r['c0_factor'],r['delay_min'])].append(r)
    results=[]
    for (panel,read,factor,delay),rr in by.items():
        results.append(dict(panel=panel,read=read,c0_factor=factor,delay_min=delay,n=len(rr),
            direct=sum(r['reason']=='direct' for r in rr),forecast=sum(r['reason']=='forecast' for r in rr),
            early_forecast=sum(r['reason']=='forecast' and r['policy_time']<180 for r in rr),
            abstain=sum(r['reason']=='abstain' for r in rr),false_early=sum(r['premature_forecast_false'] for r in rr),
            false_forecast=sum(r['forecast_false'] for r in rr),total_comparator=sum(r['assay_only_time'] for r in rr),
            total_policy=sum(r['policy_time'] for r in rr),mean_saved_min=mean(r['minutes_saved'] for r in rr)))
    write('summary.csv',results)
    for r in rows:
        if (r['read'],r['c0_factor'],r['delay_min'])==('ratio',1.,0):print(r['panel'],r['condition'],r['anchor_mg_l'],r['reason'],r['policy_time'],r['reference_at_policy_mg_l'],r['minutes_saved'],r['premature_forecast_false'])
if __name__=='__main__':main()
