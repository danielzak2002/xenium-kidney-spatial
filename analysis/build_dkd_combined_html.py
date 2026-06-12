#!/usr/bin/env python
"""Combined emailable HTML for the full DKD Xenium analysis arc: independent reannotation (validation)
-> B-rich subgroup reproduction + non-DKD references -> the BAFF/APRIL story (conditioned re-test,
ambient stress-test, receptor control + producer anchor). Self-contained (figures embedded as
downscaled JPEGs). Numbers pulled from committed CSVs."""
import os, io, base64, pandas as pd
from PIL import Image
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
RE=f"{REPO}/analysis/dkd_xenium_reannotation"; DIS=f"{REPO}/analysis/dkd_xenium_disease"
def img(path,max_w=1400,q=84):
    im=Image.open(path).convert("RGB")
    if im.width>max_w: im=im.resize((max_w,round(im.height*max_w/im.width)),Image.LANCZOS)
    b=io.BytesIO(); im.save(b,format="JPEG",quality=q,optimize=True)
    return f'<img src="data:image/jpeg;base64,{base64.b64encode(b.getvalue()).decode()}">'
val=pd.read_csv(f"{RE}/validation_summary.csv").set_index("comparison")
tst=pd.read_csv(f"{DIS}/dkd_vs_control_test.csv").set_index("metric")
rep=pd.read_csv(f"{DIS}/myeloid_baff_per_sample.csv")
ariS=val.loc["segment_vs_annotation_updated","ARI"]; agrS=val.loc["segment_vs_annotation_updated","agreement"]
ariI=val.loc["immune_vs_immune_combined","ARI"]; agrI=val.loc["immune_vs_immune_combined","agreement"]
pB=tst.loc["Blin_frac","p_value"]; medBAFF=rep.BAFF_det_pct.median()

def card(title,sidx,body): return f'<div class="card"><h2><span class="s">{sidx}</span> {title}</h2>{body}</div>'

HTML=f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>DKD Xenium — full analysis</title><style>
body{{margin:0;font-family:-apple-system,'Helvetica Neue',Arial,sans-serif;background:#f4f5f7;color:#1c1c1c;line-height:1.55}}
header{{background:#15203a;color:#fff;padding:24px 30px}} header h1{{margin:0;font-size:22px}} header p{{margin:9px 0 0;color:#b9c6e0;font-size:14px;max-width:940px}}
nav{{background:#22305a;color:#cdd8ef;font-size:12.5px;padding:9px 30px}} nav a{{color:#cdd8ef;text-decoration:none;margin-right:16px}}
main{{max-width:1000px;margin:0 auto;padding:22px 18px 70px}}
.sec{{font-size:13px;color:#6A3D9A;font-weight:700;letter-spacing:.04em;text-transform:uppercase;margin:30px 4px 6px;border-bottom:2px solid #6A3D9A;padding-bottom:4px}}
.card{{background:#fff;border-radius:12px;box-shadow:0 1px 4px rgba(0,0,0,.1);margin:14px 0;padding:18px 22px}}
h2{{font-size:16px;margin:0 0 10px;border-left:4px solid #6A3D9A;padding-left:11px}}
h2 .s{{background:#15203a;color:#fff;border-radius:5px;padding:1px 8px;font-size:12px;margin-right:7px}}
img{{width:100%;height:auto;border-radius:8px;border:1px solid #e6e8ec;display:block;margin:8px 0}}
p,li{{font-size:13.7px;color:#33384a}} .cap{{font-size:12.5px;color:#6a7180;margin-top:4px}}
.kpi{{display:flex;gap:12px;flex-wrap:wrap;margin:6px 0}}
.kpi div{{background:#f3eef8;border:1px solid #d9c9ec;border-radius:10px;padding:9px 15px;flex:1;min-width:150px}}
.kpi b{{font-size:19px;color:#6A3D9A;display:block}} .kpi span{{font-size:11px;color:#6a7180}}
.go{{background:#eef6ee;border:1px solid #9ccc9c;border-radius:8px;padding:10px 14px;color:#1a4d1a;font-size:12.6px;margin-top:8px}}
.no{{background:#fdeaea;border:1px solid #e3a6a6;border-radius:8px;padding:10px 14px;color:#7a1a1a;font-size:12.6px;margin-top:8px}}
.flag{{background:#fff7ed;border:1px solid #f1c899;border-radius:8px;padding:10px 14px;color:#7a4a12;font-size:12.6px;margin-top:8px}}
.synth{{background:#eef2f9;border:1px solid #aebbd6;border-radius:10px;padding:14px 18px;color:#16315e;font-size:14px}}
code{{background:#eef1f5;padding:1px 5px;border-radius:4px;font-size:12px}}
footer{{max-width:1000px;margin:0 auto;padding:0 18px 50px;color:#8a93a3;font-size:11.5px}}
</style></head><body>
<header><h1>DKD Xenium — independent reannotation, B-rich subgroup, and the BAFF axis</h1>
<p>The full analysis arc on the 16 Xenium samples of the diabetic-kidney-disease atlas (Dumoulin et&nbsp;al. 2026): an independent re-typing validated against the authors' labels, reproduction of the B-cell-rich subgroup, individual non-DKD references (IgAN/MN), and a stress-tested BAFF/APRIL ligand&ndash;receptor assessment. Raw data read-only; pure science. <b>All disease comparisons are descriptive / underpowered</b> (DKD n=8 vs Control n=3; references n=1&ndash;2; one patient per sample).</p></header>
<nav><b>Contents:</b> &nbsp; <a href="#a">1 · Reannotation &amp; validation</a> <a href="#b">2 · B-rich subgroup + references</a> <a href="#c">3 · BAFF/APRIL axis</a> <a href="#d">Bottom line</a></nav>
<main>

<div class="card"><div class="synth"><b>One-paragraph story.</b> An independent full-panel pipeline reproduces the authors' DKD typing (segment ARI {ariS}, immune ARI {ariI}), validating both. From those labels the paper's <b>B-rich DKD subgroup reproduces exactly</b> (= {{1006, HK2695}}, 100% concordant with the authors' B-predominant niche); <b>IgAN</b> lacks B-aggregates and <b>MN</b> forms a plasma/myeloid-mixed (non-B-follicular) aggregate. The BAFF axis was then stress-tested: <b>myeloid BAFF is a real, anchored, tissue-wide producer</b>, but <b>there is no spatially localized B-aggregate/TLS niche</b> &mdash; neither peri-aggregate stromal ligand nor in-aggregate receptor concentration survives transcript-density controls.</div></div>

<div class="sec" id="a">1 · Independent reannotation &amp; validation</div>
{card("The two independent typings agree","RESULT",f'''
<div class="kpi"><div><b>{ariS} / {agrS*100:.0f}%</b><span>segment ARI / agreement vs annotation_updated</span></div>
<div><b>{ariI} / {agrI*100:.0f}%</b><span>immune-subtype ARI / agreement</span></div>
<div><b>951,040</b><span>cells · 16 samples · 5,101 measured genes</span></div></div>
<p>Independent normalize &rarr; Harmony(sample) &rarr; Leiden &rarr; marker-based typing, blind to author labels; immune subtypes by hierarchical rules on absolute means (z-scoring mislabels the T/cytotoxic axis). High recall on every abundant lineage validates both pipelines; divergences (iPT/iTAL injury states, CD8&harr;NK, rare types) are interpretable and panel-limited.</p>''')}
{card("16-sample spatial cell-type map","FIG",img(f"{RE}/figures/spatial_grid_16.png",1400,82)+'<p class="cap">Independent cell-type maps, region-cropped, shared key (muted epithelium, saturated immune). Tissue is epithelial-dominant; immune is a scattered minority.</p>')}
{card("Confusion vs authors","FIG",img(f"{RE}/figures/concordance_matrices.png",1400,84)+'<p class="cap">Recall per author class. Strong diagonals; off-diagonals are the documented divergences (injury states, CD8&harr;NK, panel-limited rare immune types).</p>')}

<div class="sec" id="b">2 · B-rich subgroup reproduction &amp; non-DKD references</div>
{card("B-rich subgroup reproduced, 100% concordant","FIG",img(f"{DIS}/figures/dkd_subgroup.png",1400,86)+f'''<p>Split on B-aggregate burden (the paper's subgroup is TLS-based, not raw B fraction). Clean ~3&times; gap &rarr; <b>B-rich DKD = {{1006, HK2695}}</b>, which are the only DKD samples carrying the authors' <code>B predom. Immune ME</code> niche &mdash; <b>100% concordant (8/8)</b>.</p>''')}
{card("Per-participant B-lineage maps (all 16)","FIG",img(f"{DIS}/figures/b_lineage_gallery_16.png",1400,82)+'<p class="cap">B (blue) + Plasma (orange) with aggregate hulls (DBSCAN eps=50&micro;m/minPts=20). Ordered B-rich DKD &rarr; IgAN, MN (adjacent for contrast) &rarr; B-poor DKD &rarr; Control &rarr; one-offs. B-rich DKD show compact follicular aggregates; IgAN none; MN looser/mixed; Controls bare.</p>')}
{card("DKD vs Control + individual references","FIG",img(f"{DIS}/figures/dkd_vs_control.png",1400,86)+f'''<p>DKD &gt; Control in B-lineage fraction (median 1.32% vs 0.06%, Mann-Whitney p={pB}, <i>descriptive</i>).</p>
<div class="flag"><b>IgAN (1003):</b> zero B-aggregates &rarr; distinct from B-rich DKD. <b>MN (1005):</b> forms aggregates but <b>plasma/myeloid-mixed</b> (B 11% / Plasma 9% / Myeloid 12%) vs B-rich DKD B-dominated (B 36&ndash;38%) &mdash; aggregates of a different, non-B-follicular type. (n=1 each, exploratory.)</div>''')}

<div class="sec" id="c">3 · The BAFF/APRIL axis — conditioned, stress-tested, anchored</div>
{card("Step 1 — conditioned re-test (producers + space)","BAFF",img(f"{DIS}/figures/baff_april_panel.png",1400,84)+'''<div class="go"><b>BAFF producer signal is real:</b> myeloid 3.1% vs 0.39% epithelial floor = 7.9&times; (global rate 0.62% had washed it out). B cells express BAFF-R.</div><div class="no"><b>APRIL &rarr; NO-GO:</b> myeloid 5.4% below the epithelial floor 6.8% &mdash; ambient/broad, not producer-specific.</div>''')}
{card("Step 2 — peri-aggregate stromal BAFF: ambient/spillover (RETRACTED)","BAFF",img(f"{DIS}/figures/baff_ambient_control.png",1400,86)+'''<div class="no"><b>The "stromal BAFF near aggregates" claim does not survive.</b> The near/far BAFF rise is non-specific &mdash; non-producer epithelium (2.15&times;) and endothelium (2.86&times;) rise like stromal (2.44&times;), with flat transcript density &rarr; a local spillover field, not stromal production.</div>''')}
{card("Step 3a — receptors are NOT aggregate-concentrated","BAFF",img(f"{DIS}/figures/baff_receptor_control.png",1400,86)+'''<div class="no"><b>No receptor niche.</b> BAFF-R in/out ratio (1.37&times;) matches the neutral control-gene inflation (1.28&times;); BCMA below its control; TACI null; inconsistent across sections. The earlier "BAFF-R aggregate-elevated 1.6&times;" was transcript-count inflation. B/plasma express receptors <b>constitutively</b> (BAFF-R 27&ndash;37% of B, BCMA 12&ndash;17% of plasma), just not concentrated at aggregates.</div>''')}
{card("Step 3b — myeloid BAFF producer is anchored","BAFF",img(f"{DIS}/figures/baff_myeloid_anchor.png",1400,86)+f'''<div class="go"><b>Reliable producer.</b> Cell-intrinsic (per-transcriptome ~24&times; epithelium &mdash; and myeloid have <i>fewer</i> counts, so raw 7.9&times; understates it), bona fide myeloid identity (BAFF&plus; richer in CD68/CD14/C1QA/AIF1), and detected in all 16/16 samples (median {medBAFF:.1f}%).</div>''')}

<div class="sec" id="d">Bottom line</div>
<div class="card"><div class="synth"><b>What holds:</b> (1) the independent reannotation reproduces the authors' typing (ARI {ariS}/{ariI}); (2) the B-rich DKD subgroup reproduces exactly ({{1006, HK2695}}, 100% vs authors); (3) IgAN and MN are distinct from B-rich DKD (no aggregates / mixed aggregates); (4) <b>myeloid produce BAFF tissue-wide (anchored)</b> and B/plasma constitutively express the receptors &mdash; but <b>no localized B-aggregate BAFF niche</b> (neither ligand nor receptor concentration survives density/control-gene rigor). The global BAFF null was wrong about the producer, right about localization.<br><br>
<b>Caveats:</b> DKD n=8 vs Control n=3 and references n=1&ndash;2 are descriptive/underpowered; one patient per sample; architecture work inherits the reannotation typing caveats (CD8 recall 0.58, iPT/iTAL soft boundaries); receptors/ligand are sparsely detected.</div></div>

</main>
<footer>Self-contained combined summary &middot; figures from analysis/dkd_xenium_reannotation + analysis/dkd_xenium_disease &middot; commits 63647fb / ab02b20 / a76142e / 8c840ed &middot; pure science &middot; raw data read-only, git-ignored.</footer>
</body></html>"""
path=f"{REPO}/analysis/dkd_combined_summary.html"
open(path,"w").write(HTML)
print(f"wrote {path} ({len(HTML)//1024} KB)")
