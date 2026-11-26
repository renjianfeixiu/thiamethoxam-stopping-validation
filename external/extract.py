"""Digitize two publisher figures. Supply original images in a directory (see README)."""
import argparse
import csv
import hashlib
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

BASE = Path(__file__).resolve().parent
# Coordinates on the documented crops; all y guides pick a series, not a value.
PANELS = {
    "Chen5b": {
        "file": "chen5.jpg", "sha256":"22fc6c0b310130b95fc06940cda7bb7f34d5e63b6e7078f638aff08aae19ef4b",
        "x0":157.5, "x180":932.1, "y0":640.4, "y1":75.4, "c0":2e-5*291.72*1000,
        "curves":[("g-C3N4@LMPET",(.76,.88),[91,113,136,158,177,193]),
                  ("TiO2@LMPET",(.63,.75),[104,158,196,232,265,295]),
                  ("g-C3N4-TiO2@LMPET",(.18,.38),[295,436,522,581,612,629])]},
    "Farooq7d": {
        "file":"front_fig7.webp", "sha256":"5864e8109ea3e64490081171904ab5b46c7597825ad9c6c4968165a79d716b49",
        "x0":147.2, "x180":698.6, "y0":607.8, "y1":82.7, "c0":5.,
        "curves":[("0.1 g/L",(.32,.48),[105,138,143,175,196,243]),
                  ("0.3 g/L",(.53,.64),[134,169,191,219,286,338]),
                  ("0.5 g/L",None,[188,220,282,340,398,468]),
                  ("0.7 g/L",(.94,1.04),[145,193,240,313,345,413])]}
}
TIMES=(30,60,90,120,150,180)

def pixels(hsv, hue):
    h,s,v = (hsv/255.).transpose(2,0,1)
    if hue is None:
        return (s<.11)&(v>.16)&(v<.55)
    low,high=hue
    in_range=(h>=low)&(h<=high)
    if high>1: in_range |= h<high-1
    return in_range&(s>.4)&(v>.27)

def extract(source, output, overlays=None):
    rows=[]
    for panel,p in PANELS.items():
        path=source/p["file"]
        if hashlib.sha256(path.read_bytes()).hexdigest()!=p['sha256']:
            raise ValueError(f'Unexpected source image: {path}')
        picture=Image.open(path).convert('RGB')
        if panel=='Chen5b':
            picture=picture.crop((0,int(picture.height*.49),picture.width,picture.height))
            picture.thumbnail((980,780))
        else:
            picture=picture.crop((820,698,1678,1409))
        hsv=np.asarray(picture.convert('HSV'))
        draw=ImageDraw.Draw(picture)
        for label,hue,guides in p['curves']:
            m=pixels(hsv,hue)
            for t,guess in zip(TIMES,guides):
                x=round(p['x0']+t*(p['x180']-p['x0'])/180)
                y,xs=np.nonzero(m[max(0,guess-20):guess+21,x-7:x+8])
                if len(y)<20: raise ValueError((panel,label,t,len(y)))
                # Center of the marker body; largely excludes thin vertical error bars.
                counts=np.bincount(y,minlength=41)
                body=np.flatnonzero(counts>=max(3,0.50*counts.max()))
                yc=guess-20+(body.min()+body.max())/2
                # Independent of later observations: extend short incoming line to the tick.
                flank=[]
                previous=guides[TIMES.index(t)-1] if t>30 else p['y1']
                step=(p['x180']-p['x0'])/6
                for dx in range(-30,-18):
                    expected_y=round(guess+dx*(guess-previous)/step)
                    start=expected_y-15
                    yy=np.flatnonzero(m[start:expected_y+16,x+dx])
                    if len(yy):flank.append((dx,start+np.median(yy)))
                yl=float('nan')
                if len(flank)>=6:
                    slope,intercept=np.polyfit(*np.array(flank).T,1)
                    if max(abs(f[1]-(slope*f[0]+intercept)) for f in flank)<=2.5:
                        yl=intercept
                if t==60 and not np.isfinite(yl):raise ValueError(('anchor line',panel,label))
                ratio=lambda y:(p['y0']-y)/(p['y0']-p['y1'])
                rows.append(dict(panel=panel,condition=label,time_min=t,nominal_c0_mg_l=f"{p['c0']:.4f}",marker_x=x,marker_y=f'{yc:.2f}',ratio=f'{ratio(yc):.6f}',incoming_ratio=f'{ratio(yl):.6f}'))
                draw.ellipse((x-12,yc-12,x+12,yc+12),outline='orange',width=2)
        if overlays:
            overlays.mkdir(exist_ok=True,parents=True)
            picture.save(overlays/(panel+'.png'))
    output.parent.mkdir(exist_ok=True,parents=True)
    with output.open('w',newline='') as f:
        w=csv.DictWriter(f,rows[0]);w.writeheader();w.writerows(rows)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('source',type=Path)
    parser.add_argument('--output',type=Path,default=BASE/'data'/'curves.csv')
    parser.add_argument('--overlays',type=Path)
    args=parser.parse_args();extract(args.source,args.output,args.overlays)
