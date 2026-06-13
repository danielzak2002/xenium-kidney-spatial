#!/usr/bin/env python
"""Self-contained, emailable HTML summary of the B-lineage MECHANISTIC signatures analysis (extends 02).
Embeds the four new-layer figures as downscaled JPEGs; numbers pulled from the analysis CSVs.
Named *_summary.html -> git-ignored (base64-embedded); send via the file tool, do not commit."""
import os, io, base64, pandas as pd
from PIL import Image
OUT="/Users/danie/ClaudeCode/pilot_analyses/xenium/analysis/dkd_xenium_disease"; FIG=f"{OUT}/figures"
def datauri(path,max_w=1500,q=85):
    im=Image.open(path).convert("RGB")
    if im.width>max_w: im=im.resize((max_w,round(im.height*max_w/im.width)),Image.LANCZOS)
    b=io.BytesIO(); im.save(b,format="JPEG",quality=q,optimize=True)
    return "data:image/jpeg;base64,"+base64.b64encode(b.getvalue()).decode()
def img(path,max_w=1500,q=85): return f'<img src="{datauri(path,max_w,q)}">'

con=pd.read_csv(f"{OUT}/bcell_content.csv").set_index("sample")
iso=pd.read_csv(f"{OUT}/bcell_isotype.csv").set_index("sample")
tls=pd.read_csv(f"{OUT}/bcell_tls.csv").set_index("sample")
st =pd.read_csv(f"{OUT}/bcell_state.csv").set_index("sample")
def g(df,s,col):
    try: return df.loc[s,col]
    except Exception: return float("nan")

rows=[("1006","DKD B-rich"),("HK2695","DKD B-rich"),("1003","IgAN"),("1005","MN"),
      ("1004","AA amyloid"),("1007","C3GN")]
def trow(s,lbl):
    bp=g(con,s,"B_to_Plasma"); igg=g(iso,s,"pct_plasma_IgG"); iga=g(iso,s,"pct_plasma_IgA_subgate")
    nreg=g(tls,s,"n_region"); ccl=g(tls,s,"CCL19_pct"); cxcl=g(tls,s,"CXCL13_ncells")
    mzb=g(st,s,"P_MZB1"); bcma=g(st,s,"P_TNFRSF17")
    agg = "—" if (pd.isna(nreg) or nreg==0) else f"{ccl:.0f}% / {int(cxcl)}"
    hl=' style="background:#f3eef8;font-weight:600"' if lbl=="DKD B-rich" else (' style="background:#fff7ed"' if lbl in("IgAN","MN") else "")
    return (f"<tr{hl}><td>{s}</td><td>{lbl}</td><td class=n>{bp:.1f}</td>"
            f"<td class=n>{igg:.0f} ({iga:.1f})</td><td class=n>{agg}</td>"
            f"<td class=n>{mzb:.0f} / {bcma:.0f}</td></tr>")
tbody="".join(trow(s,l) for s,l in rows)

HTML=f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>B-lineage mechanism across nephropathies (anecdotal)</title><style>
body{{margin:0;font-family:-apple-system,'Helvetica Neue',Arial,sans-serif;background:#f4f5f7;color:#1c1c1c;line-height:1.55}}
header{{background:#15203a;color:#fff;padding:22px 30px}}
header h1{{margin:0;font-size:21px}} header p{{margin:8px 0 0;color:#b9c6e0;font-size:14px;max-width:920px}}
main{{max-width:980px;margin:0 auto;padding:22px 18px 70px}}
.card{{background:#fff;border-radius:12px;box-shadow:0 1px 4px rgba(0,0,0,.1);margin:16px 0;padding:18px 22px}}
h2{{font-size:17px;margin:0 0 10px;border-left:4px solid #6A3D9A;padding-left:11px}}
h2 .s{{background:#6A3D9A;color:#fff;border-radius:5px;padding:1px 8px;font-size:12px;margin-right:7px}}
img{{width:100%;height:auto;border-radius:8px;border:1px solid #e6e8ec;display:block;margin:8px 0}}
p,li{{font-size:13.7px;color:#33384a}} .cap{{font-size:12.5px;color:#6a7180;margin:4px 2px 0}}
table{{border-collapse:collapse;width:100%;font-size:12.5px;margin:6px 0}}
th,td{{border-bottom:1px solid #eceef2;padding:4px 8px;text-align:left}} th{{color:#6a7180;font-weight:600}}
td.n{{text-align:right;font-variant-numeric:tabular-nums}}
.warn{{background:#fdeaea;border:1px solid #e3a6a6;border-radius:8px;padding:12px 16px;font-size:13px;color:#7a1a1a;margin:8px 0}}
.flag{{background:#fff7ed;border:1px solid #f1c899;border-radius:8px;padding:10px 14px;font-size:12.7px;color:#7a4a12;margin-top:8px}}
.kpi{{display:flex;gap:14px;flex-wrap:wrap;margin:6px 0}}
.kpi div{{background:#f3eef8;border:1px solid #d9c9ec;border-radius:10px;padding:10px 16px;flex:1;min-width:170px}}
.kpi b{{font-size:16px;color:#6A3D9A;display:block}} .kpi span{{font-size:11.5px;color:#6a7180}}
code{{background:#eef1f5;padding:1px 5px;border-radius:4px;font-size:12px}}
footer{{max-width:980px;margin:0 auto;padding:0 18px 50px;color:#8a93a3;font-size:11.5px}}
</style></head><body>
<header><h1>B-lineage mechanism across nephropathies — anecdotal (extends the B-rich subgroup work)</h1>
<p>Layers B-lineage <b>mechanism</b> (B:Plasma split, Ig isotype, TLS organization, structural localization,
damage coupling, cell state) onto the B-lineage <b>burden</b> from summary 02. Validated reannotation labels.
Dumoulin et&nbsp;al. 2026 DKD Xenium. Raw read-only.</p></header>
<main>

<div class="warn"><b>n = 1 per non-DKD condition.</b> IgAN, MN, C3GN are single sections; AA amyloid n=2.
<b>No statistics, no p-values — nothing is tested.</b> Everything below is descriptive, single-section,
<b>hypothesis-generating</b>, framed only as <i>suggestive of a different disease process</i>. The three
sections look like three different B-lineage programs — an anecdote to test in a powered cohort, not a result.</div>

<div class="card"><h2><span class="s">RESULT</span>Three apparent B-lineage programs</h2>
<div class="kpi">
<div><b>DKD B-rich → organized lymphoid</b><span>B-predominant (B:Plasma 2–6); the only sections with follicular CXCL13 (HK2695, 303+ cells) over T-zone CCL19; B near glomeruli, far from injury.</span></div>
<div><b>MN → antibody-mediated-looking</b><span>only plasma-skewed glomerular disease (B:Plasma 0.73); IgG-dominant / IgA-absent plasma; mixed (T-zone) aggregates, not follicular.</span></div>
<div><b>IgAN → diffuse, unorganized</b><span>B-lineage present but NO aggregate / no TLS; peri-vascular not peri-glomerular; yet the most mature plasma (MZB1 88% / BCMA 37%).</span></div>
</div>
<div class="flag"><b>&#9888; Panel limits two headline axes.</b> (1) The <b>IgA-IgAN axis cannot be tested</b> — IGHA1 is
plasma-specific (88&times;) but detected in &lt;1% of plasma, and IGHM/IGKC/IGLC are absent from the Xenium panel;
only the <b>IgG</b> arm (IGHG1, 62% of plasma, 163&times;) is measurable. (2) The follicular-TLS marker <b>CXCL13 is
below the detection floor</b> (0.42%); only the T-zone chemokine <b>CCL19</b> is quantitative, with CXCL13/CXCR5/CCR7
as presence-only flags. iPT recall &asymp;0.64 → injured-tubule distances are conservative.</div></div>

<div class="card"><h2><span class="s">1 · 2</span>B:Plasma split &amp; IgG within plasma</h2>
{img(f"{FIG}/bcell_fig_A_split_isotype.png")}
<p>Splitting 02's combined B-lineage into <b>B vs Plasma</b>: DKD B-rich is B-skewed (TLS-like); <b>MN is the one
plasma-skewed glomerular disease</b> (B:Plasma 0.73) with <b>IgG-dominant, IgA-absent</b> plasma (69% IGHG1+) —
consistent with an IgG glomerular process. IgAN plasma is also IgG-detectable; its defining IgA is not measurable here.</p></div>

<div class="card"><h2><span class="s">3 · 6</span>TLS organization in the aggregate &amp; cell state</h2>
{img(f"{FIG}/bcell_fig_B_tls_state.png")}
<p><b>Only DKD B-rich (HK2695) carries substantial follicular CXCL13</b> on top of CCL19 — the most bona-fide TLS-like
structure. MN aggregates have the <b>highest T-zone CCL19 (27%) but little CXCL13</b> → mixed, not follicular.
IgAN has no aggregate region at all. Plasma maturity (MZB1/BCMA) is <b>highest in IgAN</b> and <b>lowest in AA amyloid</b>.</p></div>

<div class="card"><h2><span class="s">4 · 5</span>Localization &amp; damage&times;B coupling</h2>
{img(f"{FIG}/bcell_fig_C_localization_coupling.png")}
<p>B-lineage in <b>DKD B-rich and MN sits closer to glomeruli than to injured tubules</b> (B→injury 56–94&micro;m ≫
B→glom 26–34&micro;m), while <b>myeloid hugs injury everywhere</b> (~22–30&micro;m) — orthogonal to summary 06's injury→myeloid
axis. <b>IgAN B-cells are peri-vascular, not peri-glomerular</b> (0.49 vs 0.13).</p></div>

<div class="card"><h2><span class="s">CROPS</span>Glomerular-axis sections (IgAN / MN / DKD B-rich)</h2>
{img(f"{FIG}/bcell_fig_D_glom_crops.png")}
<p class="cap">B (blue) / Plasma (orange) / Myeloid (green) on grey tissue, glomeruli (pink), B-aggregate hulls (navy).
The per-participant 16-panel gallery is in summary 02 and is not reproduced here.</p></div>

<div class="card"><h2><span class="s">DATA</span>New-layer summary (focus rows)</h2>
<table><tr><th>sample</th><th>condition</th><th class=n>B:Plasma</th><th class=n>plasma IgG% (IgA sub)</th><th class=n>agg CCL19% / CXCL13+</th><th class=n>plasma MZB1/BCMA</th></tr>{tbody}</table>
<p class="cap">Hull composition pulled from 02: MN B11/Pl9/My12/CD4T25 vs B-rich DKD B36–38/Pl2/My4–5. Full per-sample CSVs in the analysis folder.</p></div>

<div class="card"><h2><span class="s">METHOD</span>How (one paragraph)</h2>
<p>Validated reannotation labels + coordinates; gated gene counts re-extracted per-cell from the DKD Xenium h5ad
(backed read, Xenium rows, aligned on (orig_ident, cell_id)). B-lineage burden, DBSCAN aggregates (eps=50&micro;m,
minPts=20) and hull composition are <b>pulled from 02</b>, not recomputed. New layers: B:Plasma from labels; Ig isotype
as within-plasma detection (gated); TLS markers among cells in the &le;50&micro;m aggregate region; localization by positional
nearest-anchor (&le;30&micro;m else interstitial); damage&times;B coupling as median distance to nearest injured-tubule vs
glomerulus; state as identity-conditioned detection. <b>No permutation, no test</b> — descriptive only. Folder:
<code>analysis/dkd_xenium_disease/</code> · <code>bcell_nephropathy_anecdotal.py</code> + <code>REPORT_bcell_anecdotal.md</code>.</p></div>

</main>
<footer>Self-contained summary &middot; pure science &middot; raw data read-only, git-ignored. <b>n=1 per non-DKD; descriptive; hypothesis-generating; NOT statistically tested.</b></footer>
</body></html>"""
path=f"{OUT}/bcell_anecdotal_summary.html"
open(path,"w").write(HTML)
print(f"wrote {path} ({len(HTML)//1024} KB)")
