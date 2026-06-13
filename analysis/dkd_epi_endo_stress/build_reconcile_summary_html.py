#!/usr/bin/env python
"""Self-contained, emailable HTML for the injured-PT <-> B-lineage reconcile-and-extend.
Embeds figures as downscaled JPEGs; numbers pulled from the committed reconcile_*.csv.
Named *_summary.html so .gitignore auto-excludes it (send, don't commit)."""
import os, io, base64, pandas as pd
from PIL import Image
OUT="/Users/danie/ClaudeCode/pilot_analyses/xenium/analysis/dkd_epi_endo_stress"
def img(path,max_w=1400,q=86):
    im=Image.open(path).convert("RGB")
    if im.width>max_w: im=im.resize((max_w,round(im.height*max_w/im.width)),Image.LANCZOS)
    b=io.BytesIO(); im.save(b,format="JPEG",quality=q,optimize=True)
    return f'<img src="data:image/jpeg;base64,{base64.b64encode(b.getvalue()).decode()}">'
C=pd.read_csv(f"{OUT}/reconcile_correlations.csv")
S=pd.read_csv(f"{OUT}/reconcile_spatial_per_section.csv")
G=pd.read_csv(f"{OUT}/reconcile_usability_gate.csv")
def cor(x,y): r=C[(C.x==x)&(C.y==y)].iloc[0]; return f"{r.rho} [{r.ci_lo}, {r.ci_hi}]"
k_prog=int((S.prog_dz>0).sum()); n=len(S); k_exc=int((S.prog_minus_ctrl>0).sum())
k_perm=int((S.iPTfrac_perm_p<0.05).sum()); amb=int((S.nonPT_amb_dz>0).sum())
def srow(r):
    sig="✓" if r.iPTfrac_perm_p<0.05 else ""
    exc=' style="background:#eef6ee"' if r.prog_minus_ctrl>0 else ' style="background:#fdeaea"'
    return (f"<tr><td>{r.section}</td><td>{r.Condition}</td><td class=n>{r.prog_dz:+.2f}</td>"
            f"<td class=n>{r.ctrl_floor_dz:+.2f}</td><td class=n{exc[6:]}>{r.prog_minus_ctrl:+.2f}</td>"
            f"<td class=n>{r.nonPT_amb_dz:+.2f}</td><td class=n>{r.iPTfrac_near:.2f}/{r.iPTfrac_far:.2f}</td>"
            f"<td class=n>{r.iPTfrac_perm_p:.3f} {sig}</td></tr>")
gaterows="".join(f"<tr><td>{r.gene}</td><td class=n>{r.PT_detect:.3f}</td><td class=n>{r.ratio}×</td>"
                 f"<td>{'✓' if r.usable else '✗ (&lt;3% PT)'}</td></tr>" for _,r in G.iterrows())

HTML=f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>injured-PT ↔ B-lineage — reconcile</title><style>
body{{margin:0;font-family:-apple-system,'Helvetica Neue',Arial,sans-serif;background:#f4f5f7;color:#1c1c1c;line-height:1.55}}
header{{background:#15203a;color:#fff;padding:22px 30px}} header h1{{margin:0;font-size:20px}} header p{{margin:8px 0 0;color:#b9c6e0;font-size:14px;max-width:920px}}
main{{max-width:980px;margin:0 auto;padding:22px 18px 70px}}
.card{{background:#fff;border-radius:12px;box-shadow:0 1px 4px rgba(0,0,0,.1);margin:16px 0;padding:18px 22px}}
h2{{font-size:16px;margin:0 0 10px;border-left:4px solid #6A3D9A;padding-left:11px}}
h2 .s{{background:#6A3D9A;color:#fff;border-radius:5px;padding:1px 8px;font-size:12px;margin-right:7px}}
img{{width:100%;height:auto;border-radius:8px;border:1px solid #e6e8ec;display:block;margin:8px 0}}
p,li{{font-size:13.7px;color:#33384a}} .cap{{font-size:12.5px;color:#6a7180}}
table{{border-collapse:collapse;width:100%;font-size:12.6px;margin:6px 0}} th,td{{border-bottom:1px solid #eceef2;padding:5px 8px;text-align:left}} th{{color:#6a7180}} td.n{{text-align:right;font-variant-numeric:tabular-nums}}
.kpi{{display:flex;gap:12px;flex-wrap:wrap;margin:6px 0}}
.kpi div{{background:#f3eef8;border:1px solid #d9c9ec;border-radius:10px;padding:9px 15px;flex:1;min-width:160px}}
.kpi b{{font-size:18px;color:#6A3D9A;display:block}} .kpi span{{font-size:11px;color:#6a7180}}
.no{{background:#fdeaea;border:1px solid #e3a6a6;border-radius:8px;padding:11px 15px;color:#7a1a1a;font-size:12.8px;margin-top:8px}}
.mid{{background:#fff7ed;border:1px solid #f1c899;border-radius:8px;padding:11px 15px;color:#7a4a12;font-size:12.8px;margin-top:8px}}
.synth{{background:#eef2f9;border:1px solid #aebbd6;border-radius:10px;padding:14px 18px;color:#16315e;font-size:14px}}
code{{background:#eef1f5;padding:1px 5px;border-radius:4px;font-size:12px}}
footer{{max-width:980px;margin:0 auto;padding:0 18px 50px;color:#8a93a3;font-size:11.5px}}
</style></head><body>
<header><h1>Injured-PT ↔ B-lineage — reconcile &amp; extend under post-BAFF rigor</h1>
<p>Re-examines the prior tubular-injury-near-B-aggregates result (author labels + B-only aggregates: injPT near&gt;far Δz <b>+0.13, 6/9, p=0.038</b>) using the <b>validated reannotation labels</b> (incl. iPT) + <b>B-lineage (B+Plasma)</b> aggregates, with the same control-gene / specificity rigor that re-scoped the BAFF axis. Unit of replication = <b>sample</b>; descriptive throughout; read-only.</p></header>
<main>
<div class="card"><div class="synth"><b>Verdict.</b> The injured-PT ↔ B-lineage association is <b>real as co-localization but neither B-lineage-specific nor a program-intensity gradient.</b> Tubular injury tracks <b>general, myeloid-led immune infiltration</b>; injured-PT-<i>labeled</i> cells sit near immune aggregates compositionally, but the per-cell injury-<i>program intensity</i> near aggregates does not exceed neutral-gene transcript inflation. <b>The prior +0.13 spatial reading is retracted.</b></div>
<div class="kpi">
<div><b>0.07 [−0.46, 0.60]</b><span>partial iPT~B-lineage controlling for total-immune → collapses</span></div>
<div><b>1 of 9</b><span>sections where injPT near-rise exceeds the control-gene floor</span></div>
<div><b>6 of 9</b><span>sections with iPT cell-fraction enriched near aggregates (permutation)</span></div>
</div></div>

<div class="card"><h2><span class="s">PRIMARY</span>Cross-sample correlation + specificity (n = 16)</h2>
{img(f"{OUT}/reconcile_fig_primary.png",1400,86)}
<table><tr><th>x</th><th>y</th><th class=n>Spearman ρ [95% CI]</th></tr>
<tr><td>iPT frac</td><td>B-lineage frac</td><td class=n>{cor('iPT_frac_epi','Blin_frac')}</td></tr>
<tr><td>iPT frac</td><td>B-aggregate count</td><td class=n>{cor('iPT_frac_epi','n_agg')}</td></tr>
<tr style="background:#f3eef8"><td>iPT frac</td><td><b>total-immune frac</b></td><td class=n>{cor('iPT_frac_epi','total_imm_frac')}</td></tr>
<tr style="background:#f3eef8"><td>iPT frac</td><td><b>myeloid frac</b></td><td class=n>{cor('iPT_frac_epi','myeloid_frac')}</td></tr>
<tr style="background:#eef6ee"><td>iPT frac</td><td><b>B-lineage frac | total-immune (partial)</b></td><td class=n>{cor('iPT_frac_epi','Blin_frac | total_imm')}</td></tr>
</table>
<div class="no"><b>NOT B-lineage-specific.</b> iPT correlates with B-lineage (ρ 0.67) but <b>more strongly with total-immune (0.79) and myeloid (0.82)</b>, and the partial correlation collapses to <b>0.07</b> → report as "injury tracks general (myeloid-led) infiltration." The hard iPT-label fraction and the continuous injPT program score agree.</div></div>

<div class="card"><h2><span class="s">SECONDARY</span>Spatial near/far across the 9 aggregate-bearing sections</h2>
{img(f"{OUT}/reconcile_fig_spatial.png",1500,86)}
<table><tr><th>section</th><th>cond</th><th class=n>injPT Δz</th><th class=n>ctrl floor</th><th class=n>injPT−ctrl</th><th class=n>non-PT amb</th><th class=n>iPT-frac n/f</th><th class=n>perm p</th></tr>
{''.join(srow(r) for _,r in S.iterrows())}</table>
<div class="mid"><b>Program-intensity rise does not clear the control floor.</b> injPT near&gt;far {k_prog}/{n}, but exceeds the neutral control-gene floor only <b>{k_exc}/{n}</b> (control genes inflate +0.35 near aggregates vs injPT +0.07 pooled) → generic transcript inflation, not an injPT-specific gradient. What survives is <b>compositional</b>: the iPT cell-fraction is enriched near aggregates in <b>{k_perm}/{n}</b> sections (permutation) — but that co-localization is non-B-specific (PRIMARY). Non-PT cells also rise in {amb}/{n} (spillover). B-rich {{1006, HK2695}} disagree (1006 depleted).</div></div>

<div class="card"><h2><span class="s">GATE</span>Usability (injPT genes, PT detection vs immune floor)</h2>
<table><tr><th>gene</th><th class=n>PT detect</th><th class=n>vs floor</th><th>usable</th></tr>{gaterows}</table>
<p class="cap">4 of 5 pass; ITGB6 drops at the 3% PT-detection gate. Scored on VCAM1/HAVCR1/PROM1/SPP1. injPT genes, DBSCAN eps=50/minPts=20 inherited verbatim from the prior analysis.</p></div>

<div class="card"><h2><span class="s">METHOD</span>How &amp; caveats</h2>
<p>Validated reannotation labels (iPT/PT/B/Plasma) + coords; B-lineage = B+Plasma; DBSCAN aggregates (eps=50µm/minPts=20); CP-median log1p module score; near ≤50µm / far &gt;200µm to nearest aggregate member; within-section per-compartment z; control-gene panel (TPT1/PPIA/YWHAZ/TMSB10/UBB) = inflation floor; within-section permutation null for the iPT cell-fraction (positions+count fixed). <b>Underpowered</b> (n=16 samples / 9 aggregate-bearing sections; no patient column → donor clustering uncontrolled; iPT recall ~0.64 → both label-fraction and program-score reported; centroids + expression only, no morphology). Associational, not causal. Folder <code>analysis/dkd_epi_endo_stress/</code>; <code>reconcile_extend.py</code>, <code>REPORT_reconcile.md</code>.</p></div>
</main>
<footer>Self-contained · honest attenuation/retraction · reconciles prior +0.13/6-of-9 · raw data read-only, git-ignored.</footer>
</body></html>"""
open(f"{OUT}/reconcile_summary.html","w").write(HTML)
print(f"wrote reconcile_summary.html ({len(HTML)//1024} KB)")
