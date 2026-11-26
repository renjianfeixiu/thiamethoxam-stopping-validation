"""Extract Table S1 and a fixed grid of Figure S8 from the publisher's supplement."""
from pathlib import Path
from io import BytesIO
import argparse
import hashlib
import urllib.request
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parent
URL = "https://ars.els-cdn.com/content/image/1-s2.0-S138589472407520X-mmc1.docx"
W = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
# Pixel centers for t=0, t=180 and relative concentration 0, 1, calibrated from axis ticks.
AXES = [(267,1155.5,829.5,106.5), (1468.5,2361.5,830,103),
        (2681,3576,831,103), (264,1149.5,1821,1101),
        (1470,2358,1822.5,1100), (2682,3575.5,1825,1098.5)]
PANELS = [range(1,6),range(6,11),range(11,16),range(16,21),range(21,27),range(27,33)]
HUES = [(76,94),(14,26),(222,245),(135,152),(182,202)]
TIMES = (60,90,120,150,180)

def table(z):
    doc = ET.fromstring(z.read("word/document.xml"))
    source_table = doc.find("w:body/w:tbl", W)
    records = []
    for row in source_table.findall("w:tr", W)[1:-1]:
        cells = ["".join(t.text or "" for t in c.findall(".//w:t", W)) for c in row.findall("w:tc", W)]
        for offset in (0,7):
            event, ct, ph, c0, dose, light, de = cells[offset:offset+7]
            corrected = event=="E25" and dose=="50..2"
            rec = dict(event=int(event[1:]), catalyst={"1":"CN/Bi2O3","2":"CN"}[ct],
                       pH=float(ph), c0=float(c0), dose_mg=float("50.2" if corrected else dose),
                       dose_as_printed=dose, light_level=int(light), de_pct=float(de.rstrip("%")))
            rec["ratio180"] = 1-rec["de_pct"]/100
            rec["ct180"] = rec["ratio180"]*rec["c0"]
            records.append(rec)
    df = pd.DataFrame(records).sort_values("event")
    assert df.event.tolist()==list(range(1,33))
    return df

def figure(z):
    img = Image.open(BytesIO(z.read("word/media/image23.jpeg"))).convert("RGB")
    if img.size!=(3740,2057):
        raise ValueError("Image changed; recalibrate axes before digitizing.")
    hsv = np.asarray(img.convert("HSV"))
    rows = []
    for panel, events in enumerate(PANELS):
        x0,x180,y0,y1 = AXES[panel]
        for color,event in enumerate(events):
            for i,t in enumerate(TIMES):
                x = round(x0+(x180-x0)*t/180)
                # At 60 min the printed legend is on the lower left of some axes.
                lo_ratio = 0.40 if t==60 else 0.0
                top,bottom=round(y1-8), round(y0+(y1-y0)*lo_ratio)
                strip=hsv[top:bottom+1,x-5:x+6]
                h,s,v=strip.transpose(2,0,1)
                if color<5:
                    low,high=HUES[color]
                    mask=(h>=low)&(h<=high)&(s>=100)&(v>=75)
                else:
                    # Grey/black triangles of E26 and E32. Windows read from the
                    # printed graph, kept broad (±0.055) to isolate labels and axes.
                    approx={26:[.84,.76,.69,.64,.59],32:[.66,.50,.36,.27,.21]}[event][i]
                    coords=(np.arange(top,bottom+1)-y0)/(y1-y0)
                    mask=(v<85)&(s<80)&(np.abs(coords[:,None]-approx)<.055)
                # These four colored markers are visibly separate in the source;
                # a broad y-window excludes the printed legend or off-series pixels.
                guides={(3,60):.82,(5,120):.95,(24,90):.87,(31,150):.23}
                if (event,t) in guides:
                    coords=(np.arange(top,bottom+1)-y0)/(y1-y0)
                    mask &= np.abs(coords[:,None]-guides[event,t])<.04
                ys,xs=np.nonzero(mask)
                status="read"
                if len(ys)<15 or (len(ys) and np.ptp(ys)>35):
                    estimate=np.nan
                    status="overlap_or_insufficient_pixels"
                else:
                    estimate=round((top+np.median(ys)-y0)/(y1-y0),4)
                marker_ratio=estimate
                # Hollow colored markers have different edge contrasts. Estimate
                # their center from nearby line segments on either side where visible.
                if status=="read":
                    ymid=round(y0+(y1-y0)*estimate)
                    flank_estimates=[]
                    sides=[(-35,-17)] if t==180 else [(-35,-17),(17,35)]
                    for xa,xb in sides:
                        band=hsv[ymid-28:ymid+29,x+xa:x+xb+1]
                        bh,bs,bv=band.transpose(2,0,1)
                        if color<5:
                            bm=(bh>=low)&(bh<=high)&(bs>=85)&(bv>=70)
                        else:
                            bm=(bv<90)&(bs<80)
                        cols=np.flatnonzero(bm.any(axis=0))
                        if len(cols)<10: continue
                        centers=np.array([np.median(np.flatnonzero(bm[:,c])) for c in cols])
                        slope,intercept=np.polyfit(cols+xa,centers,1)
                        if np.max(np.abs(centers-(slope*(cols+xa)+intercept)))>3: continue
                        r=(ymid-28+intercept-y0)/(y1-y0)
                        if abs(r-estimate)<.024: flank_estimates.append(r)
                    if flank_estimates and np.ptp(flank_estimates)<.018:
                        estimate=round(float(np.mean(flank_estimates)),4)
                        status="line_and_marker"
                    else:
                        status="marker_only"
                rows.append(dict(event=event,time=t,figure_ratio=estimate,marker_ratio=marker_ratio,pixels=len(ys),status=status))
    return pd.DataFrame(rows)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--supplement",type=Path)
    args=p.parse_args()
    blob=args.supplement.read_bytes() if args.supplement else urllib.request.urlopen(URL,timeout=60).read()
    with zipfile.ZipFile(BytesIO(blob)) as z:
        d=table(z)
        f=figure(z)
    out=ROOT/"data"
    out.mkdir(exist_ok=True)
    d.to_csv(out/"table_s1.csv",index=False,float_format="%.6f")
    f.loc[f.time==60,["event","time","marker_ratio"]].to_csv(out/"figure_s8_marker60.csv",index=False)
    f.drop(columns="marker_ratio").to_csv(out/"figure_s8.csv",index=False)
    print("Source SHA-256",hashlib.sha256(blob).hexdigest())
    chk=f[f.time==180].merge(d,on="event")
    chk["abs_error"]=(chk.figure_ratio-chk.ratio180).abs()
    print(chk[["event","figure_ratio","ratio180","abs_error"]].to_string(index=False))
    print(f.pivot(index="event",columns="time",values="figure_ratio").round(3))
if __name__=="__main__":
    main()
