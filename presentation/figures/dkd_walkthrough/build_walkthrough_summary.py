#!/usr/bin/env python
"""Standalone, emailable version of the DKD walkthrough.

The committed dkd_walkthrough.html references figures by RELATIVE PATH (lightweight, in-repo). This
thin transform reads that file and inlines every <img> as a downscaled base64 JPEG, producing a single
self-contained file: dkd_walkthrough_summary.html. The *_summary.html suffix is REQUIRED — .gitignore
auto-excludes it (base64-embedded HTML trips secret scanning). Send via the file tool; do NOT commit.

Run:  conda run -n spatial python presentation/figures/dkd_walkthrough/build_walkthrough_summary.py
(run make_walkthrough.py first so the html + panels exist.)
"""
import os, io, re, base64
from PIL import Image
HERE=os.path.dirname(os.path.abspath(__file__))
SRC=os.path.join(HERE,"dkd_walkthrough.html"); DST=os.path.join(HERE,"dkd_walkthrough_summary.html")

def datauri(path,max_w=1500,q=84):
    im=Image.open(path).convert("RGB")
    if im.width>max_w: im=im.resize((max_w,round(im.height*max_w/im.width)),Image.LANCZOS)
    b=io.BytesIO(); im.save(b,format="JPEG",quality=q,optimize=True)
    return "data:image/jpeg;base64,"+base64.b64encode(b.getvalue()).decode()

html=open(SRC).read()
seen={}
def repl(m):
    rels=m.group(1); p=os.path.normpath(os.path.join(HERE,rels))
    if p not in seen:
        seen[p]=datauri(p); print(f"  embed {os.path.getsize(p)//1024:>5} KB  {rels}")
    return f'<img src="{seen[p]}" style="width:100%"'
html=re.sub(r'<img src="([^"]+)"', repl, html)

# add a small emailable banner so a forwarded file is self-explanatory
banner=('<div style="max-width:1020px;margin:0 auto;padding:10px 18px 0;color:#7a828f;font-size:12px">'
        'Self-contained copy — all figures embedded; forward as a single file. '
        'Pure science; single-section disease claims are hypotheses (see frame). Raw data read-only.</div>')
html=html.replace("<main>", banner+"<main>",1)

open(DST,"w").write(html)
print(f"\nwrote {DST} ({len(html)//1024} KB self-contained, {len(seen)} images embedded)")
