#!/usr/bin/env python
"""Self-contained, emailable HTML summary of the DKD Xenium B-lineage disease analysis.
Embeds figures as downscaled JPEGs (single shareable file). Numbers pulled from the CSVs."""
import os, io, base64, pandas as pd
from PIL import Image
OUT="/Users/danie/ClaudeCode/pilot_analyses/xenium/analysis/dkd_xenium_disease"; FIG=f"{OUT}/figures"
def datauri(path,max_w=1400,q=85):
    im=Image.open(path).convert("RGB")
    if im.width>max_w: im=im.resize((max_w,round(im.height*max_w/im.width)),Image.LANCZOS)
    b=io.BytesIO(); im.save(b,format="JPEG",quality=q,optimize=True)
    return "data:image/jpeg;base64,"+base64.b64encode(b.getvalue()).decode()
def img(path,max_w=1400,q=85): return f'<img src="{datauri(path,max_w,q)}">'

sub=pd.read_csv(f"{OUT}/per_sample_substrate.csv"); sub["orig_ident"]=sub.orig_ident.astype(str)
tst=pd.read_csv(f"{OUT}/dkd_vs_control_test.csv").set_index("metric")
det=pd.read_csv(f"{OUT}/baff_detection_by_celltype.csv").set_index("cell_type")
pB=tst.loc["Blin_frac","p_value"]
order=["1006","HK2695","1013","1008","1012","1011","1001","1010","HK2753","HK3106","HK3626","1003","1005","1004","1009","1007"]
grp={"1006":"B-rich DKD","HK2695":"B-rich DKD","1003":"IgAN","1005":"MN"}
sub=sub.set_index("orig_ident")
def srow(s):
    r=sub.loc[s]; g=grp.get(s, "DKD" if r.Condition=="DKD" else ("Control" if r.Condition=="Control" else r.Condition))
    hl=' style="background:#f3eef8;font-weight:600"' if g=="B-rich DKD" else (' style="background:#fff7ed"' if g in("IgAN","MN") else "")
    return f"<tr{hl}><td>{s}</td><td>{g}</td><td class=n>{int(r.n_cells):,}</td><td class=n>{r.Blin_frac*100:.2f}</td><td class=n>{int(r.n_agg)}</td><td class=n>{r.agg_cells_per10k:.0f}</td><td class=n>{int(r.author_Bpredom_ME_20um)}</td></tr>"
subrows="".join(srow(s) for s in order)

HTML=f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>DKD Xenium — B-lineage disease analysis</title><style>
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
.kpi{{display:flex;gap:14px;flex-wrap:wrap;margin:6px 0}}
.kpi div{{background:#f3eef8;border:1px solid #d9c9ec;border-radius:10px;padding:10px 16px;flex:1;min-width:150px}}
.kpi b{{font-size:21px;color:#6A3D9A;display:block}} .kpi span{{font-size:11.5px;color:#6a7180}}
.go{{background:#eef6ee;border:1px solid #9ccc9c;border-radius:8px;padding:11px 15px;font-size:12.8px;color:#1a4d1a;margin-top:8px}}
.nogo{{background:#fdeaea;border:1px solid #e3a6a6;border-radius:8px;padding:11px 15px;font-size:12.8px;color:#7a1a1a;margin-top:8px}}
.flag{{background:#fff7ed;border:1px solid #f1c899;border-radius:8px;padding:10px 14px;font-size:12.6px;color:#7a4a12;margin-top:8px}}
code{{background:#eef1f5;padding:1px 5px;border-radius:4px;font-size:12px}}
footer{{max-width:980px;margin:0 auto;padding:0 18px 50px;color:#8a93a3;font-size:11.5px}}
</style></head><body>
<header><h1>DKD Xenium — B-lineage disease analysis: B-rich subgroup, IgAN/MN references, BAFF/APRIL</h1>
<p>Built on the validated independent reannotation. <b>B-lineage = B + Plasma.</b> Reproduce the paper's B-rich DKD subgroup within DKD, place IgAN/MN and other non-DKD as individual references, and re-test the BAFF/APRIL ligand axis under a correctly conditioned design. Dumoulin et&nbsp;al. 2026. Raw read-only.</p></header>
<main>

<div class="card"><h2><span class="s">RESULT</span>Headline numbers</h2>
<div class="kpi">
<div><b>{{1006, HK2695}}</b><span>B-rich DKD subgroup &mdash; 100% concordant (8/8) with authors' B-predom niche</span></div>
<div><b>BAFF: producer GO</b><span>myeloid 7.9&times; epithelial floor (tissue-wide); localized-stromal-niche claim retracted (spillover &mdash; see stress-test)</span></div>
<div><b>APRIL: NO-GO</b><span>myeloid below epithelial floor &mdash; ambient, not a producer/niche signal</span></div>
</div>
<p class="cap"><b>All comparisons descriptive / underpowered:</b> DKD n=8 vs Control n=3; non-DKD references n=1&ndash;2; one patient per sample. Architecture work inherits the reannotation typing caveats (CD8 recall 0.58, iPT/iTAL soft boundaries).</p>
</div>

<div class="card"><h2><span class="s">STEP 2</span>B-rich subgroup reproduced &amp; validated</h2>
{img(f"{FIG}/dkd_subgroup.png",1400,86)}
<p>The paper's "B cell-rich subgroup" is defined by B-<b>aggregates</b>/TLS, not raw B fraction (which mis-ranks the large HK2695 section). On B-aggregate burden the DKD distribution has a clean ~3&times; gap between the top-2 (1006, HK2695) and the rest &rarr; <b>B-rich DKD = {{1006, HK2695}}</b>. The authors' <code>B predom. Immune ME</code> niche is present in DKD <b>only</b> in these two (797 &amp; 1,576 cells; zero in the other six) &mdash; <b>100% concordant</b>. Reproduction, not discovery.</p></div>

<div class="card"><h2><span class="s">STEP 4</span>Per-participant B-lineage maps (all 16)</h2>
{img(f"{FIG}/b_lineage_gallery_16.png",1500,84)}
<p class="cap">B (blue) + Plasma (orange) on grey tissue, with B-lineage aggregate hulls (DBSCAN eps=50&micro;m/minPts=20). Ordered B-rich DKD &rarr; IgAN, MN (placed adjacent for contrast) &rarr; B-poor DKD &rarr; Control &rarr; other one-offs. B-rich DKD show compact B follicular aggregates; IgAN has none; MN has looser mixed aggregates; Controls are essentially bare.</p></div>

<div class="card"><h2><span class="s">STEP 3</span>DKD vs Control + non-DKD references</h2>
{img(f"{FIG}/dkd_vs_control.png",1400,86)}
<p>DKD &gt; Control in B-lineage fraction (median 1.32% vs 0.06%, Mann-Whitney p={pB}, <i>descriptive</i>); aggregate burden trends higher but n.s. (driven by the two B-rich samples). Individual non-DKD references overlaid as labelled points.</p>
<div class="flag"><b>IgAN (1003):</b> B-lineage 0.82%, <b>zero B-aggregates</b> &mdash; sits at the B-poor / no-TLS end, distinct from B-rich DKD. <b>MN (1005):</b> B-lineage 3.1% and it <b>does form aggregates</b>, but they are <b>plasma/myeloid-mixed</b> (in/around &le;50&micro;m: B 11% &middot; Plasma 9% &middot; Myeloid 12% &middot; CD4 T 25%) vs B-rich DKD which is strongly <b>B-dominated</b> (B 36&ndash;38% / Plasma 2% / Myeloid 4&ndash;5%). MN forms aggregates but of a different, non-B-follicular type. (n=1 each, exploratory.)</div></div>

<div class="card"><h2><span class="s">STEP 5</span>BAFF/APRIL conditioned re-assessment</h2>
{img(f"{FIG}/baff_april_panel.png",1500,86)}
<p>All five genes are on-panel and measured on Xenium. Detection was conditioned on <b>producers</b> (not a tissue-wide average, which washes out producer signal) and on <b>space</b> (near vs far from B-aggregates).</p>
<div class="go"><b>BAFF (TNFSF13B) &rarr; partial GO &mdash; overturns the prior global null at the producer level.</b> Myeloid detect <b>3.1% vs 0.39% epithelial floor = 7.9&times;</b> (global rate was 0.62% &mdash; the average hid it); B cells express BAFF-R (35.7%). A real <b>tissue-wide myeloid producer</b> signal.</div>
<div class="flag"><b>&#9888; Retracted sub-claim:</b> the earlier "peri-aggregate <b>stromal</b> BAFF (2.1&ndash;2.5&times;)" did <b>not</b> survive an ambient/spillover stress-test (<code>REPORT_baff_ambient.md</code>): the near/far BAFF rise is non-specific &mdash; non-producer epithelium (2.15&times;) and endothelium (2.86&times;) rise like stromal (2.44&times;), with flat transcript density &rarr; a <b>local spillover field, not stromal production</b>. Keep myeloid tissue-wide BAFF; drop the localized-niche reading. <b>Follow-up</b> (<code>REPORT_baff_receptor.md</code>): the myeloid producer is <b>anchored</b> (cell-intrinsic, 24&times; epithelium per-transcriptome, reproducible 16/16), and the receptors are <b>constitutive but NOT aggregate-concentrated</b> (BAFF-R "1.6&times;" matched the neutral control-gene inflation). Net: real tissue-wide myeloid BAFF, no localized niche.</div>
<div class="nogo"><b>APRIL (TNFSF13) &rarr; NO-GO &mdash; global null confirmed.</b> Myeloid 5.4% is <b>below</b> the epithelial floor (6.8%) &rarr; broad/ambient, not producer-specific. BCMA on plasma is real (17%) but not aggregate-concentrated.</div></div>

<div class="card"><h2><span class="s">DATA</span>Per-sample B-lineage substrate</h2>
<table><tr><th>sample</th><th>group</th><th class=n>cells</th><th class=n>B-lin %</th><th class=n>n agg</th><th class=n>agg cells/10k</th><th class=n>author B-predom ME</th></tr>{subrows}</table>
<p class="cap">B-lineage aggregates via DBSCAN (eps=50&micro;m, minPts=20). Author B-predom ME = cells in the authors' "B predom. Immune ME" niche (20&micro;m). Shaded: B-rich DKD (purple), IgAN/MN (amber).</p></div>

<div class="card"><h2><span class="s">METHOD</span>How (one paragraph)</h2>
<p>Validated reannotation labels (B-lineage = B + Plasma). Per sample: B-lineage fraction + DBSCAN aggregate burden (eps=50&micro;m/minPts=20), normalised per-10k and per-mm&sup2;. Within-DKD split on aggregate burden; validated vs the authors' B-predom niche. Non-DKD kept as individual references (no pooling). BAFF/APRIL detection conditioned on producer cell types and on distance to B-aggregates; receptors for context. Folder: <code>analysis/dkd_xenium_disease/</code> (commit ab02b20); reproduce with <code>substrate.py &rarr; baff_april.py &rarr; figures.py</code>.</p></div>

</main>
<footer>Self-contained summary &middot; figures from analysis/dkd_xenium_disease/figures &middot; pure science &middot; raw data read-only, git-ignored. Descriptive/underpowered (n=8 vs 3; references n=1&ndash;2; one patient per sample).</footer>
</body></html>"""
path=f"{OUT}/dkd_disease_summary.html"
open(path,"w").write(HTML)
print(f"wrote {path} ({len(HTML)//1024} KB)")
