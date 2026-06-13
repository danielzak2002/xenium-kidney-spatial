#!/usr/bin/env python
"""Self-contained, emailable HTML summary of the per-sample composition-by-group layer.
Embeds the two dotplots as downscaled JPEGs; numbers pulled from composition_by_group.csv.
Named *_summary.html -> git-ignored (base64-embedded); send via the file tool, do not commit."""
import os, io, base64, numpy as np, pandas as pd
from PIL import Image
OUT="/Users/danie/ClaudeCode/pilot_analyses/xenium/analysis/dkd_xenium_disease"; FIG=f"{OUT}/figures"
def datauri(path,max_w=1600,q=85):
    im=Image.open(path).convert("RGB")
    if im.width>max_w: im=im.resize((max_w,round(im.height*max_w/im.width)),Image.LANCZOS)
    b=io.BytesIO(); im.save(b,format="JPEG",quality=q,optimize=True)
    return "data:image/jpeg;base64,"+base64.b64encode(b.getvalue()).decode()
def img(path,max_w=1600,q=85): return f'<img src="{datauri(path,max_w,q)}">'

df=pd.read_csv(f"{OUT}/composition_by_group.csv")
GORDER=["Control","DKD","IgAN","MN","AA amyloid","C3GN"]; NSING={"IgAN","MN","C3GN"}
def med(res,ct):
    d=df[(df.resolution==res)&(df.cell_type==ct)]
    return {g:float(np.median(d[d.group==g].fraction))*100 for g in GORDER}
LIN=["Epithelial","Immune","Stroma","Endothelial"]; IMM=["Myeloid","CD4 T","CD8 T","B","Plasma"]
def trow(name,vals,bold=None):
    cells="".join(f"<td class=n{' style=font-weight:700' if (bold and g==bold) else ''}>{vals[g]:.1f}</td>" for g in GORDER)
    return f"<tr><td>{name}</td>{cells}</tr>"
lin_rows="".join(trow(ct,med("coarse_lineage",ct)) for ct in LIN)
imm_rows="".join(trow(ct,med("immune_subtype",ct)) for ct in IMM)
hPT=med("epithelial_subtype","PT"); iPT=med("epithelial_subtype","iPT")
epi_rows=trow("healthy PT",hPT)+trow("injured iPT",iPT)
hdr="".join(f"<th class=n>{g}{'*' if g in NSING else ''}</th>" for g in GORDER)

HTML=f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Composition by disease group (per sample)</title><style>
body{{margin:0;font-family:-apple-system,'Helvetica Neue',Arial,sans-serif;background:#f4f5f7;color:#1c1c1c;line-height:1.55}}
header{{background:#15203a;color:#fff;padding:22px 30px}}
header h1{{margin:0;font-size:21px}} header p{{margin:8px 0 0;color:#b9c6e0;font-size:14px;max-width:920px}}
main{{max-width:1000px;margin:0 auto;padding:22px 18px 70px}}
.card{{background:#fff;border-radius:12px;box-shadow:0 1px 4px rgba(0,0,0,.1);margin:16px 0;padding:18px 22px}}
h2{{font-size:17px;margin:0 0 10px;border-left:4px solid #6A3D9A;padding-left:11px}}
h2 .s{{background:#6A3D9A;color:#fff;border-radius:5px;padding:1px 8px;font-size:12px;margin-right:7px}}
img{{width:100%;height:auto;border-radius:8px;border:1px solid #e6e8ec;display:block;margin:8px 0}}
p,li{{font-size:13.7px;color:#33384a}} .cap{{font-size:12.5px;color:#6a7180;margin:4px 2px 0}}
table{{border-collapse:collapse;width:100%;font-size:12.5px;margin:6px 0}}
th,td{{border-bottom:1px solid #eceef2;padding:4px 8px;text-align:left}} th{{color:#6a7180;font-weight:600}}
td.n,th.n{{text-align:right;font-variant-numeric:tabular-nums}}
.warn{{background:#fdeaea;border:1px solid #e3a6a6;border-radius:8px;padding:12px 16px;font-size:13px;color:#7a1a1a;margin:8px 0}}
.flag{{background:#fff7ed;border:1px solid #f1c899;border-radius:8px;padding:10px 14px;font-size:12.7px;color:#7a4a12;margin-top:8px}}
code{{background:#eef1f5;padding:1px 5px;border-radius:4px;font-size:12px}}
footer{{max-width:1000px;margin:0 auto;padding:0 18px 50px;color:#8a93a3;font-size:11.5px}}
</style></head><body>
<header><h1>Cell-type composition per sample, by disease group</h1>
<p>Top-level layer for the walkthrough. <b>Reuses the validated reannotation labels</b> (951,040 cells / 16 samples) —
no typing re-run. Per-sample fractions at coarse lineage and immune-subtype resolution, one dot per sample.
Dumoulin et&nbsp;al. 2026 DKD Xenium. Raw read-only.</p></header>
<main>

<div class="warn"><b>Descriptive · NO statistics.</b> IgAN, MN, C3GN are <b>single sections (n=1)</b>; AA amyloid n=2;
DKD n=8; Control n=3. Group-median bars are drawn <b>only where n&gt;1</b>; single dots (starred*) are one section.
Fractions are compositional (sum-to-one) — a CLR sensitivity row accompanies each panel and tracks the raw ordering.</div>

<div class="card"><h2><span class="s">COARSE</span>Lineage composition (epithelial / immune / stroma / endothelial)</h2>
{img(f"{FIG}/composition_by_group_coarse.png")}
<p>Control is strongly <b>epithelial-dominant (74%) with a near-bare immune compartment (2.8%)</b>. Every disease group
shows <b>epithelial loss with immune + stromal expansion</b>. The CLR (bottom row) tracks the raw fractions — the
ordering is not a closure artifact.</p>
<table><tr><th>% of section cells</th>{hdr}</tr>{lin_rows}</table></div>

<div class="card"><h2><span class="s">DRILL-DOWN</span>Immune-subtype composition</h2>
{img(f"{FIG}/composition_by_group_immune.png")}
<p>Single-section standouts: <b>MN</b> the most immune/plasma-skewed glomerular disease (immune 24%, plasma 1.8% —
echoing the B-lineage analysis); <b>C3GN</b> the most T-skewed (T 14%); <b>AA amyloid</b> the highest B-cell fraction (3.5%).
DKD's 8 sections fan widely; the IgAN/MN/C3GN single dots mostly sit within that spread except on their own standout axis.</p>
<table><tr><th>% of section cells</th>{hdr}</tr>{imm_rows}</table>
<div class="flag"><b>&#9888; NK and DC are not separately typed</b> in the validated labels, so they are not reported as their
own classes. The consistent epithelial shift is visible at sub-type level too: <b>healthy PT collapses as injured iPT rises</b>.</div>
<table><tr><th>epithelial % of section cells</th>{hdr}</tr>{epi_rows}</table></div>

<div class="card"><h2><span class="s">METHOD</span>How (one paragraph)</h2>
<p>Validated reannotation labels (<code>cells.parquet</code>, summary 01); no typing recomputed. Per-sample fraction of
total section cells at three resolutions — coarse lineage (<code>my_lineage</code>), immune subtype (<code>my_label</code>),
epithelial subtype (<code>my_coarse</code>); <code>Unresolved</code> (0.6%) excluded from denominators. Dotplots: one dot per
sample, x = group (Control→DKD→IgAN→MN→AA→C3GN), median bar only where n&gt;1; CLR small-multiple (0.5 pseudocount,
centred log-ratio over each resolution's parts) as a closure-aware sensitivity view — <b>not tested</b>. Folder:
<code>analysis/dkd_xenium_disease/</code> · <code>composition_by_group.py</code> + <code>REPORT_composition_by_group.md</code> + <code>composition_by_group.csv</code>.</p></div>

</main>
<footer>Self-contained summary &middot; pure science &middot; raw data read-only, git-ignored. <b>Descriptive; n=1 per non-DKD flagged; NO statistics.</b></footer>
</body></html>"""
path=f"{OUT}/composition_by_group_summary.html"
open(path,"w").write(HTML)
print(f"wrote {path} ({len(HTML)//1024} KB)")
