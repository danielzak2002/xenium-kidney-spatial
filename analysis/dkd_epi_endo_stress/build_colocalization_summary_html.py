#!/usr/bin/env python
"""Self-contained emailable HTML for the injury<->immune co-localization addendum.
Named *_summary.html so .gitignore auto-excludes it (send, don't commit). Numbers from coloc_*.csv."""
import os, io, base64, pandas as pd
from PIL import Image
OUT="/Users/danie/ClaudeCode/pilot_analyses/xenium/analysis/dkd_epi_endo_stress"
def img(path,max_w=1400,q=85):
    im=Image.open(path).convert("RGB")
    if im.width>max_w: im=im.resize((max_w,round(im.height*max_w/im.width)),Image.LANCZOS)
    b=io.BytesIO(); im.save(b,format="JPEG",quality=q,optimize=True)
    return f'<img src="data:image/jpeg;base64,{base64.b64encode(b.getvalue()).decode()}">'
A=pd.read_csv(f"{OUT}/coloc_panelA_partials.csv"); A2=pd.read_csv(f"{OUT}/coloc_panelA_clr.csv")
C2=pd.read_csv(f"{OUT}/coloc_panelC2_immune_injury.csv")
def arow(r):
    hl=' style="background:#eef6ee"' if "NON-B" in r.test or "myeloid" in r.test else ""
    return f"<tr{hl}><td>{r.test}</td><td class=n>{r.rho} [{r.lo}, {r.hi}]</td></tr>"
arows="".join(arow(r) for _,r in A.iterrows())
k_my=int((C2.Myeloid_p<0.05).sum())
HTML=f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Injury ↔ immune co-localization</title><style>
body{{margin:0;font-family:-apple-system,'Helvetica Neue',Arial,sans-serif;background:#f4f5f7;color:#1c1c1c;line-height:1.55}}
header{{background:#15203a;color:#fff;padding:22px 30px}} header h1{{margin:0;font-size:20px}} header p{{margin:8px 0 0;color:#b9c6e0;font-size:14px;max-width:930px}}
main{{max-width:980px;margin:0 auto;padding:22px 18px 70px}}
.card{{background:#fff;border-radius:12px;box-shadow:0 1px 4px rgba(0,0,0,.1);margin:16px 0;padding:18px 22px}}
h2{{font-size:16px;margin:0 0 10px;border-left:4px solid #2ca02c;padding-left:11px}}
h2 .s{{background:#2ca02c;color:#fff;border-radius:5px;padding:1px 8px;font-size:12px;margin-right:7px}}
img{{width:100%;height:auto;border-radius:8px;border:1px solid #e6e8ec;display:block;margin:8px 0}}
p,li{{font-size:13.7px;color:#33384a}} .cap{{font-size:12.5px;color:#6a7180}}
table{{border-collapse:collapse;width:100%;font-size:12.7px;margin:6px 0}} th,td{{border-bottom:1px solid #eceef2;padding:5px 8px;text-align:left}} th{{color:#6a7180}} td.n{{text-align:right;font-variant-numeric:tabular-nums}}
.kpi{{display:flex;gap:12px;flex-wrap:wrap;margin:6px 0}}
.kpi div{{background:#eef6ee;border:1px solid #9ccc9c;border-radius:10px;padding:9px 15px;flex:1;min-width:160px}}
.kpi b{{font-size:18px;color:#1a4d1a;display:block}} .kpi span{{font-size:11px;color:#5a6}}
.no{{background:#fdeaea;border:1px solid #e3a6a6;border-radius:8px;padding:11px 15px;color:#7a1a1a;font-size:12.8px;margin-top:8px}}
.synth{{background:#eef2f9;border:1px solid #aebbd6;border-radius:10px;padding:14px 18px;color:#16315e;font-size:14px}}
code{{background:#eef1f5;padding:1px 5px;border-radius:4px;font-size:12px}}
footer{{max-width:980px;margin:0 auto;padding:0 18px 50px;color:#8a93a3;font-size:11.5px}}
</style></head><body>
<header><h1>Injury ↔ immune co-localization across nephropathies</h1>
<p>Door-closing addendum to the DKD consolidation: does the B-rich subgroup track tubular injury? <b>No — injury co-localizes with general, myeloid-led immune infiltration, not B-lineage specifically.</b> Compositional / descriptive; unit = sample/section; honest nulls first-class; no per-gene DE. Builds on the reannotation-label reconcile (B-lineage=B+Plasma aggregates, eps=50/minPts=20).</p></header>
<main>
<div class="card"><div class="synth"><b>Verdict.</b> Tubular injury (iPT) <b>co-localizes with myeloid infiltration</b> (cross-sample ρ 0.82) and injured epithelial states sit near immune aggregates — but <b>B-lineage adds nothing beyond non-B immune</b> (de-circularized partial 0.13, CI spans 0; CLR-robust), and the <b>near-injury infiltrate is myeloid-skewed / B-depleted</b>. The B-rich subgroup's significance is architectural/compositional, <b>not a tubular-injury role</b>.</div>
<div class="kpi">
<div><b>ρ 0.82 [0.46, 0.95]</b><span>iPT ~ myeloid fraction (the surviving positive)</span></div>
<div><b>0.13 [−0.37, 0.62]</b><span>iPT ~ B-lineage | non-B immune → collapses</span></div>
<div><b>{k_my} of 16</b><span>sections: myeloid enriched near injured tubule (perm p&lt;.05)</span></div>
</div></div>

<div class="card"><h2><span class="s">A</span>Cross-sample specificity — myeloid leads, B not specific</h2>
{img(f"{OUT}/coloc_fig_A_specificity.png",1400,86)}
<table><tr><th>test (raw fractions, n=16)</th><th class=n>Spearman ρ [95% CI]</th></tr>{arows}</table>
<p class="cap"><b>De-circularized:</b> total-immune contains B, so the original partial had B controlling for itself; against non-B immune the iPT~B partial is 0.13 (CI spans 0). <b>CLR closure-robust:</b> clr(iPT)~clr(Myeloid) {A2.iloc[1].rho} ≥ clr(iPT)~clr(B) {A2.iloc[0].rho}; non-B partial {A2.iloc[3].rho}. Closure biases positive part-correlations downward → the iPT~immune signals are, if anything, understated.</p></div>

<div class="card"><h2><span class="s">B</span>Cross-nephropathy — not disease-specific, not B-specific</h2>
{img(f"{OUT}/coloc_fig_B1_disease.png",1400,86)}
{img(f"{OUT}/coloc_fig_B2_spatial.png",1400,86)}
<p>iPT-fraction rises with both total-immune and myeloid; the three <b>Controls cluster in the low-injury/low-immune corner</b>. Near B-aggregates, iPT is enriched in 6/11 sections (incl. one control HK3626), myeloid in 4/11; <b>flat in MN, no aggregates in IgAN</b>, depleted in the B-richest 1006. ⚠ = spillover section.</p>
{img(f"{OUT}/coloc_fig_B3_gallery.png",1500,84)}
<p class="cap">One representative section per disease: injured tubule (iPT/iTAL) + immune (myeloid/B/plasma) with B-aggregate hulls.</p></div>

<div class="card"><h2><span class="s">C</span>Which cells participate</h2>
{img(f"{OUT}/coloc_fig_C1_epi_subtype.png",1400,86)}
<p><b>C1 — epithelium shifts to injured states near aggregates:</b> healthy PT depleted, iPT/iTAL enriched (beyond stroma/endo bystanders); injured gate passes 6/11 sections (fails in the B-richest 1006).</p>
{img(f"{OUT}/coloc_fig_C2_immune_injury.png",1400,86)}
<div class="no"><b>C2 — the near-injury infiltrate is MYELOID-skewed, B-depleted.</b> Within ≤30µm of an injured tubule cell, myeloid is enriched in {k_my}/16 sections (myeloid = 36–71% of the near-injury immune compartment), while B-lineage is depleted. <b>Injury attracts myeloid, not B.</b></div></div>

<div class="card"><h2><span class="s">Caveats</span></h2><p>n=16 samples / 11 aggregate-bearing; <b>n=1 per non-DKD nephropathy</b> (IgAN/MN/C3GN) — single sections, descriptive; iPT recall ~0.64 (label-fraction; reconcile showed program-score agrees); no donor column; spillover sections flagged; centroids + expression only (no morphology); associational, not causal. Stage-2 per-gene DE deliberately not run. Folder <code>analysis/dkd_epi_endo_stress/</code>; <code>colocalization.py</code>, <code>REPORT_colocalization.md</code>.</p></div>
</main>
<footer>Self-contained · door-closing exploratory addendum · no B-lineage→kidney-damage claim · raw data read-only, git-ignored.</footer>
</body></html>"""
open(f"{OUT}/colocalization_summary.html","w").write(HTML)
print(f"wrote colocalization_summary.html ({len(HTML)//1024} KB)")
