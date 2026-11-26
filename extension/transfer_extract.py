"""Read normalized curves from one primary and three external studies; source details in the package README."""
from pathlib import Path
from io import BytesIO
import argparse
import hashlib
import urllib.request
import zipfile
import numpy as np
import pandas as pd
from PIL import Image,ImageDraw

BASE=Path(__file__).resolve().parent
SOURCE=BASE/'sources'
OUT=BASE/'data'
SOURCES={
  **{f'reactor{n}.png':f'https://article.sciencepublishinggroup.com/article/figure/2241250/2241250_Figure_{n}.png' for n in (5,6)},
  'mdpifig':'https://mdpi-res.com/catalysts/catalysts-13-00516/article_deploy/html/images/catalysts-13-00516-g001.png',
  'front_fig7.webp':'https://www.frontiersin.org/files/Articles/1677416/xml-images/fmats-12-1677416-g007.webp',
}
# x(0), horizontal pixels per minute, y(r=0), vertical pixels per unit r, image dimensions.
CAL={
 's6':('FigS6.jpg',865.,3948/180,4694.14,4417.14,(5091,4197),12),
 '6f':('main_115.png',1158.5,400.5/180,909.,326.,(2129,1445),2),
 '6g':('main_115.png',1697.,404.5/180,912.,329.,(2129,1445),2),
 'F7b':('front_fig7.webp',959.5,643/210,653.,582.6,(1678,1409),3),
 'M1a':('mdpifig',219.,1234/90,1130.5,1050.,(3193,1329),3)}
# Visually assigned loose windows select a trace at an actual labeled marker.
# Input and target values are calculated from colored pixels, not from these guides.
SERIES=[
 ('s6','tap water','green',[30,60,90,120,150,180],[.81,.68,.57,.47,.39,.33]),
 ('s6','secondary effluent','orange',[30,60,90,120,150,180],[.84,.71,.60,.51,.41,.35]),
 ('6f','CN 5 mg/L','green',[30,60,90,120,150,180],[.82,.68,.59,.49,.42,.34]),
 ('6f','CN 10 mg/L','orange',[30,60,90,120,150,180],[.90,.82,.74,.65,.58,.52]),
 ('6f','CN 20 mg/L','pink',[30,60,90,120,150,180],[.94,.86,.80,.78,.74,.69]),
 ('6f','CN 50 mg/L','blue',[30,60,90,120,150,180],[.94,.89,.87,.85,.83,.81]),
 ('6g','CN/Bi2O3 5 mg/L','green',[30,60,90,120,150,180],[.75,.57,.42,.30,.21,.15]),
 ('6g','CN/Bi2O3 10 mg/L','orange',[30,60,90,120,150,180],[.78,.55,.44,.33,.23,.17]),
 ('6g','CN/Bi2O3 20 mg/L','pink',[30,60,90,120,150,180],[.85,.71,.58,.47,.39,.32]),
 ('6g','CN/Bi2O3 50 mg/L','blue',[30,60,90,120,150,180],[.92,.84,.78,.73,.66,.63]),
 ('F7b','Ni-PbS','grey',[30,60,90,120,150,180,210],[.81,.75,.62,.51,.40,.27,.21]),
 ('F7b','PbS','red',[30,60,90,120,150,180,210],[.90,.83,.76,.65,.61,.50,.41]),
 ('F7b','no catalyst','blue',[30,60,90,120,150,180,210],[.99,.98,.97,.96,.95,.94,.93]),
 ('M1a','P25 20 ppm','grey',[15,30,45,60,75,90],[.50,.25,.16,.09,.05,.02]),
 ('M1a','P25 10 ppm','olive',[15,30,45,60,75,90],[.44,.18,.105,.07,.029,.012]),
 ('M1a','P25 5 ppm','blue',[15,30,45,60,75,90],[.25,.06,.022,.012,.006,.003]),
 ('M1a','P25 1 ppm','red',[15,30,45,60,75,90],[.20,.027,.012,.005,.003,.002])]
RANGES={'green':(72,98),'orange':(10,30),'pink':(221,246),'blue':(131,158),'olive':(49,65)}
def mask(ch,c):
    h,s,v=ch.transpose(2,0,1)
    if c=='grey':return (s<55)&(v<150)
    if c=='red':return ((h<12)|(h>245))&(s>65)&(v>85)
    lo,hi=RANGES[c];return (h>=lo)&(h<=hi)&(s>60)&(v>65)
def reactor_series():
    # Here the published zero is BEFORE dark equilibration. Use the marker at
    # printed t=30 min as the lamp-on denominator (article methods + UV-ON label).
    markers={
      5:([69,101,134,167,199,232,265],156.5,125.,
         [[58,73,107,128,135,140,150],[45,69,84,99,106,115,111]]),
      6:([70,102,133,165,196,227,258],161.,127.,
         [[50,67,98,105,116,126,128],[37,42,57,64,65,71,79]])}
    rows=[];details=[]
    for panel,(xx,yzero,scale,series) in markers.items():
        pic=Image.open(SOURCE/f'reactor{panel}.png').convert('RGB')
        grey=np.asarray(pic).mean(axis=2)
        dr=ImageDraw.Draw(pic)
        for j,guides in enumerate(series):
            label=('TiO2' if panel==5 else 'PST')+f' system {j+1}'
            denominator=None
            for i,(x,guide) in enumerate(zip(xx,guides)):
                h=grey[guide-4:guide+5,x-1:x+2]
                mask=(h>=119)&(h<=139) if j==0 else h<35
                yy,_=np.nonzero(mask)
                assert len(yy)>=5,(panel,label,i)
                center=guide-4+float(np.median(yy))
                raw=(yzero-center)/scale
                if i==0:denominator=raw
                ratio=raw/denominator
                details.append(dict(panel=panel,series=label,published_time=30+30*i,
                                    lamp_time=30*i,raw_ratio=raw,light_on_baseline=denominator,ratio=ratio))
                if i: rows.append(dict(panel=f'R{panel}',series=label,time=30*i,ratio=round(ratio,5),
                                   status='read_renormalized',mask_n=len(yy),marker_span=(np.ptp(yy)/scale)))
                dr.ellipse((x-4,center-4,x+4,center+4),outline='#c14e30')
        pic.resize((636,380)).save(BASE/'qa'/f'overlay_reactor{panel}.png')
    pd.DataFrame(details).to_csv(OUT/'reactor_lamp_origin.csv',index=False)
    return rows

def extract():
    imgs={k:Image.open(SOURCE/c[0]).convert('RGB') for k,c in CAL.items()}
    hsv={k:np.asarray(v.convert('HSV')) for k,v in imgs.items()}
    rows=[]
    for panel,series,color,times,guides in SERIES:
        filename,x0,dx,y0,dy,wh,strip=CAL[panel]
        assert imgs[panel].size==wh
        for t,guide in zip(times,guides):
            x=round(x0+t*dx);ymid=round(y0-guide*dy)
            radius=round(dy*(.025 if panel=='M1a' else .042))
            lower=ymid-radius; upper=ymid+radius
            narrow=hsv[panel][lower:upper+1,x-strip:x+strip+1]
            m=mask(narrow,color); yy,xx=np.nonzero(m); n=len(yy)
            r=np.nan;status='obscured';span=np.nan
            if n>=7:
                # Midrange is stable for symmetric error bars; trim tails to remove antialiasing.
                ylo,yhi=np.quantile(yy,[.10,.90]);center=lower+(ylo+yhi)/2
                r=(y0-center)/dy;span=(yhi-ylo)/dy
                status='read'
                if panel=='M1a' and guide<=.012:
                    r=np.nan;status='near_axis'; # no point substitution at the visual floor
                elif not 0<=r<=1.03: r=np.nan;status='obscured'
            # CN-20 and CN-50 at 30 min coincide to within marker width in panel 6f.
            # Hue-only extraction can pick an edge of the wrong visible center.
            if panel=='6f' and t==30 and ('20 mg' in series or '50 mg' in series):
                r=np.nan;status='overlapping_markers'
            # Large black triangle in F7b: the polyline passes through the marker center,
            # so use the extrapolated left flank, independent of later time points.
            if panel=='F7b' and status=='read':
                offset=np.arange(-25,-13)
                ys=[];ds=[]
                for xo in offset:
                    band=hsv[panel][lower:upper+1,x+xo:x+xo+1]
                    v=np.flatnonzero(mask(band,color)[:,0])
                    if len(v):ys.append(lower+np.median(v));ds.append(xo)
                if len(ys)>=8:
                    slope,z=np.polyfit(ds,ys,1)
                    if max(abs(np.array(ys)-(slope*np.array(ds)+z)))<2.7 and abs((y0-z)/dy-r)<.026:
                        r=(y0-z)/dy
            rows.append(dict(panel=panel,series=series,time=t,ratio=round(float(r),5),status=status,
                             mask_n=n,marker_span=round(float(span),5)))
            if np.isfinite(r):
                dr=ImageDraw.Draw(imgs[panel]);y=y0-r*dy;rr=max(5,round(dy*.016))
                dr.ellipse((x-rr,y-rr,x+rr,y+rr),outline='#4a4a4a',width=max(1,round(dy*.002)))
    for k,v in imgs.items():
        if k=='6g':v=v.crop((1630,535,2129,960))
        if k=='6f':v=v.crop((1080,535,1600,960))
        if k=='M1a':v=v.crop((0,0,1640,1329))
        if k=='F7b':v=v.crop((840,0,1678,702))
        v.thumbnail((800,650));v.save(BASE/'qa'/('overlay_'+k+'.jpg'))
    rows.extend(reactor_series())
    pd.DataFrame(rows).to_csv(OUT/'transfer_curves.csv',index=False)
def main():
    p=argparse.ArgumentParser()
    p.add_argument('--paper',type=Path,help='Xie et al. original publisher PDF')
    p.add_argument('--supplement',type=Path,help='Xie et al. publisher DOCX')
    args=p.parse_args()
    for folder in [SOURCE,OUT,BASE/'qa']:folder.mkdir(exist_ok=True)
    if args.paper:
        import pymupdf
        with pymupdf.open(args.paper) as doc:
            page=next(pg for pg in doc if 'Fig. 6. (a)' in pg.get_text())
            ref=next(k[0] for k in page.get_images() if (k[2],k[3])==(2129,1445))
            payload=doc.extract_image(ref)
            assert payload['ext']=='png';(SOURCE/'main_115.png').write_bytes(payload['image'])
    if args.supplement:
        with zipfile.ZipFile(args.supplement) as z:(SOURCE/'FigS6.jpg').write_bytes(z.read('word/media/image21.jpeg'))
    for file,url in SOURCES.items():
        if not (SOURCE/file).is_file():
            with urllib.request.urlopen(url,timeout=40) as resp:(SOURCE/file).write_bytes(resp.read())
    extract()
    for n in sorted(set(c[0] for c in CAL.values())|set(SOURCES)):
        print(n,hashlib.sha256((SOURCE/n).read_bytes()).hexdigest())
if __name__=='__main__':main()
