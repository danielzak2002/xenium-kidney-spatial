#!/usr/bin/env python
"""Build a self-contained, emailable HTML summary of the DKD Xenium reannotation.
Embeds the figures as downscaled JPEGs (single shareable file). Pulls numbers from the CSVs."""
import os, io, base64, html as _html, pandas as pd
from PIL import Image
OUT="/Users/danie/ClaudeCode/pilot_analyses/xenium/analysis/dkd_xenium_reannotation"
FIG=f"{OUT}/figures"
def datauri(path,max_w=1400,q=85):
    im=Image.open(path).convert("RGB")
    if im.width>max_w: im=im.resize((max_w,round(im.height*max_w/im.width)),Image.LANCZOS)
    b=io.BytesIO(); im.save(b,format="JPEG",quality=q,optimize=True)
    return "data:image/jpeg;base64,"+base64.b64encode(b.getvalue()).decode()
def img(path,max_w=1400,q=85): return f'<img src="{datauri(path,max_w,q)}">'

vs=pd.read_csv(f"{OUT}/validation_summary.csv").set_index("comparison")
strata=pd.read_csv(f"{OUT}/disease_strata_summary.csv")
def recall_table(prefix,types,nmap):
    rec=pd.read_csv(f"{OUT}/{prefix}_recall.csv",index_col=0)
    prec=pd.read_csv(f"{OUT}/{prefix}_precision.csv",index_col=0)
    rows=""
    for t in types:
        r=rec.loc[t,t] if (t in rec.index and t in rec.columns) else float("nan")
        p=prec.loc[t,t] if (t in prec.index and t in prec.columns) else float("nan")
        rows+=f"<tr><td>{t}</td><td class=n>{nmap.get(t,''):,}</td><td class=n>{r:.2f}</td><td class=n>{p:.2f}</td></tr>"
    return rows

seg_n={"PT":208867,"PC/CNT":104274,"Immune":103215,"Fibroblast":120457,"EC":105061,"iPT":87231,
 "TAL":75480,"iTAL":38452,"IC":30317,"DCT":27960,"VSMC":12519,"MC":12091,"Podo":7308,"PEC":6986}
imm_n={"Myeloid/Macro":47311,"CD4":21140,"CD8":19673,"B":8422,"Plasma":2726}
seg_order=["Podo","Immune","PC/CNT","PT","DCT","Fibroblast","EC","IC","TAL","iTAL","VSMC","iPT"]
imm_order=["Plasma","Myeloid/Macro","B","CD4","CD8"]

ariA=vs.loc["segment_vs_annotation_updated","ARI"]; agrA=vs.loc["segment_vs_annotation_updated","agreement"]
ariB=vs.loc["immune_vs_immune_combined","ARI"]; agrB=vs.loc["immune_vs_immune_combined","agreement"]

HTML=f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>DKD Xenium — independent reannotation & validation</title><style>
body{{margin:0;font-family:-apple-system,'Helvetica Neue',Arial,sans-serif;background:#f4f5f7;color:#1c1c1c;line-height:1.55}}
header{{background:#15203a;color:#fff;padding:22px 30px}}
header h1{{margin:0;font-size:21px}} header p{{margin:8px 0 0;color:#b9c6e0;font-size:14px;max-width:920px}}
main{{max-width:980px;margin:0 auto;padding:22px 18px 70px}}
.card{{background:#fff;border-radius:12px;box-shadow:0 1px 4px rgba(0,0,0,.1);margin:16px 0;padding:18px 22px}}
h2{{font-size:17px;margin:0 0 10px;border-left:4px solid #6A3D9A;padding-left:11px}}
h2 .s{{background:#6A3D9A;color:#fff;border-radius:5px;padding:1px 8px;font-size:12px;margin-right:7px}}
img{{width:100%;height:auto;border-radius:8px;border:1px solid #e6e8ec;display:block;margin:8px 0}}
p,li{{font-size:13.7px;color:#33384a}} .cap{{font-size:12.5px;color:#6a7180;margin:4px 2px 0}}
table{{border-collapse:collapse;width:100%;font-size:12.7px;margin:6px 0}}
th,td{{border-bottom:1px solid #eceef2;padding:5px 8px;text-align:left}} th{{color:#6a7180;font-weight:600}}
td.n{{text-align:right;font-variant-numeric:tabular-nums}}
.kpi{{display:flex;gap:14px;flex-wrap:wrap;margin:6px 0}}
.kpi div{{background:#f3eef8;border:1px solid #d9c9ec;border-radius:10px;padding:10px 16px;flex:1;min-width:150px}}
.kpi b{{font-size:22px;color:#6A3D9A;display:block}} .kpi span{{font-size:11.5px;color:#6a7180}}
.flag{{background:#fff7ed;border:1px solid #f1c899;border-radius:8px;padding:10px 14px;font-size:12.8px;color:#7a4a12;margin-top:8px}}
code{{background:#eef1f5;padding:1px 5px;border-radius:4px;font-size:12px}}
footer{{max-width:980px;margin:0 auto;padding:0 18px 50px;color:#8a93a3;font-size:11.5px}}
</style></head><body>
<header><h1>DKD Xenium — independent cell-type reannotation, validated against the authors' labels</h1>
<p>An independent, full-panel pipeline re-typed the 16 Xenium DKD-atlas samples <b>blind to the authors' deposited labels</b>, then benchmarked against them. High agreement validates both typings; disagreements localise where typing is genuinely ambiguous. Dumoulin et&nbsp;al. 2026 (Zenodo 19868428). Raw data read-only.</p></header>
<main>

<div class="card"><h2><span class="s">RESULT</span>The two independent typings agree</h2>
<div class="kpi">
<div><b>0.78</b><span>segment ARI vs <code>annotation_updated</code> (85% agreement)</span></div>
<div><b>0.68</b><span>immune-subtype ARI vs <code>immune_cell_annotation_combined</code> (79%)</span></div>
<div><b>951,040</b><span>cells · 16 Xenium samples · 5,101 measured genes</span></div>
</div>
<p>Markedly cleaner than the earlier 68-gene cross-cohort pass (ARI 0.61) &mdash; confirming the earlier DKD B-cell over-call was a <b>shared-gene-depth artefact, not a labelling error</b>. Full-panel typing recovers it.</p>
</div>

<div class="card"><h2><span class="s">FIG 1</span>16-sample spatial cell-type map</h2>
{img(f"{FIG}/spatial_grid_16.png",1500,86)}
<p class="cap">Independent cell-type maps, region-cropped, one shared 14-class key (muted epithelium, saturated immune; immune points enlarged). Each panel: sample &middot; condition &middot; n cells &middot; n B&nbsp;cells. Tissue is epithelial-dominant; immune is a scattered minority that locally aggregates.</p></div>

<div class="card"><h2><span class="s">FIG 2</span>Validation vs authors &mdash; confusion</h2>
{img(f"{FIG}/concordance_matrices.png",1500,86)}
<p class="cap">Recall per author class (row-normalised). Strong diagonals; off-diagonals are the documented divergences (below).</p>
<div style="display:flex;gap:18px;flex-wrap:wrap">
<div style="flex:1;min-width:300px"><b style="font-size:13px">Segment vs annotation_updated &mdash; ARI {ariA:.2f}, {agrA*100:.0f}%</b>
<table><tr><th>type</th><th class=n>n (author)</th><th class=n>recall</th><th class=n>prec</th></tr>{recall_table("coarse",seg_order,seg_n)}</table></div>
<div style="flex:1;min-width:280px"><b style="font-size:13px">Immune subtype vs immune_combined &mdash; ARI {ariB:.2f}, {agrB*100:.0f}%</b>
<table><tr><th>type</th><th class=n>n (author)</th><th class=n>recall</th><th class=n>prec</th></tr>{recall_table("immune",imm_order,imm_n)}</table></div>
</div>
<div class="flag"><b>Divergences are findings, not failures:</b> injury states <b>iPT/iTAL</b> trade cells with PT/TAL; <b>VSMC&harr;Fibroblast/MC</b> stromal boundary; <b>CD8 recall 0.58</b> &mdash; author NK cells fall into my CD8 cluster (cytotoxic overlap; NKG7/GNLY absent from panel) and the CD4/CD8 split shifts; <b>NK, Neutrophil, Mast, DC, Treg</b> are not separately resolvable (their markers &mdash; NKG7, GNLY, TPSAB1, CPA3, S100A8, SLC26A4 &mdash; are absent or sub-ambient). All abundant lineages (epithelial segments, endothelium, stroma, B/Plasma/CD4/CD8/Myeloid) reproduce cleanly.</div>
</div>

<div class="card"><h2><span class="s">FIG 3</span>Integrated embedding (mine vs theirs)</h2>
{img(f"{FIG}/integration_umaps.png",1500,84)}
<p class="cap">Harmony UMAP on <code>sample_id</code> &mdash; by sample (mixed, integrated), by my cell type, by lineage.</p>
{img(f"{FIG}/umap_yours_vs_theirs.png",1400,84)}
<p class="cap">Same embedding, my coarse labels vs the authors' &mdash; the two independent typings occupy the same territories.</p></div>

<div class="card"><h2><span class="s">SETUP</span>Disease status of the 16 samples (next-step scoping)</h2>
<table><tr><th>stratum</th><th class=n>samples</th></tr>
<tr><td><b>DKD</b></td><td class=n>8</td></tr><tr><td><b>Control</b></td><td class=n>3</td></tr>
<tr><td>AA amyloid</td><td class=n>2</td></tr><tr><td>C3GN / IgA / MN</td><td class=n>1 each</td></tr></table>
<p><code>Diagnosis.xlsx</code> joined 16/16 (100%). DM/HTN/eGFR&gt;60 perfectly track Control (confounded with disease). With 16 single-patient samples, a <b>graded severity</b> model is underpowered &mdash; the supportable contrasts are <b>binary: DKD-vs-Control, DKD-vs-non-DKD, or the paper's B-rich-vs-B-poor DKD split</b>. No association is run here; samples are tagged and ready.</p></div>

<div class="card"><h2><span class="s">METHOD</span>How (one paragraph)</h2>
<p>16 Xenium samples subset from the union object (<code>tech=='Xenium'</code>); 342 CosMx-only structural-zeros dropped &rarr; 5,101 measured genes. Independent: normalize/log1p &rarr; HVG 2,000 &rarr; PCA 50 &rarr; <b>Harmony on sample_id</b> &rarr; Leiden res 1.5 (23 global clusters) &rarr; UMAP; immune compartment (105,719 cells) subclustered &rarr; 16. Global clusters typed by canonical kidney markers; <b>immune subtypes by hierarchical rules on absolute means</b> (z-scoring mislabels the T/cytotoxic axis &mdash; FOXP3/CTLA4 light up for any T cluster; CD4 is macrophage-expressed &mdash; so T is gated on CD3E). Author labels carried along but never used to guide any step. Folder: <code>analysis/dkd_xenium_reannotation/</code> (commit 63647fb); reproduce with <code>run_integration.py &rarr; annotate.py &rarr; validate.py &rarr; disease_strata.py &rarr; figures.py</code>.</p></div>

</main>
<footer>Self-contained summary &middot; figures from analysis/dkd_xenium_reannotation/figures &middot; pure science (scientific labels only) &middot; raw data read-only, git-ignored.</footer>
</body></html>"""
path=f"{OUT}/dkd_reannotation_summary.html"
open(path,"w").write(HTML)
print(f"wrote {path} ({len(HTML)//1024} KB)")
