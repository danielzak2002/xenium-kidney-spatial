#!/usr/bin/env python
"""Self-contained, emailable HTML for the BAFF receptor nCount-control + myeloid producer anchor."""
import os, io, base64, pandas as pd
from PIL import Image
OUT="/Users/danie/ClaudeCode/pilot_analyses/xenium/analysis/dkd_xenium_disease"; FIG=f"{OUT}/figures"
def img(path,max_w=1500,q=88):
    im=Image.open(path).convert("RGB")
    if im.width>max_w: im=im.resize((max_w,round(im.height*max_w/im.width)),Image.LANCZOS)
    b=io.BytesIO(); im.save(b,format="JPEG",quality=q,optimize=True)
    return f'<img src="data:image/jpeg;base64,{base64.b64encode(b.getvalue()).decode()}">'
A=pd.read_csv(f"{OUT}/receptor_aggregate_control.csv"); P=A[A.sample_set=="POOLED"].set_index("receptor")
rep=pd.read_csv(f"{OUT}/myeloid_baff_per_sample.csv"); b1=pd.read_csv(f"{OUT}/myeloid_anchor_ncount.csv").iloc[0]
def rrow(r,tgt):
    p=P.loc[r]; return f"<tr><td>{r} <span style='color:#888'>({tgt})</span></td><td class=n>{p.det_in_pct}</td><td class=n>{p.det_out_pct}</td><td class=n>{p.det_ratio}&times;</td><td class=n>{p.normExpr_ratio}&times;</td><td class=n>{p.ctrl_gene_ratio}&times;</td><td>{p.verdict}</td></tr>"
rrows=rrow("BAFF-R","B")+rrow("TACI","B+Plasma")+rrow("BCMA","Plasma")
HTML=f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>BAFF receptors + myeloid anchor</title><style>
body{{margin:0;font-family:-apple-system,'Helvetica Neue',Arial,sans-serif;background:#f4f5f7;color:#1c1c1c;line-height:1.55}}
header{{background:#15203a;color:#fff;padding:22px 30px}} header h1{{margin:0;font-size:20px}} header p{{margin:8px 0 0;color:#b9c6e0;font-size:14px;max-width:910px}}
main{{max-width:960px;margin:0 auto;padding:22px 18px 70px}}
.card{{background:#fff;border-radius:12px;box-shadow:0 1px 4px rgba(0,0,0,.1);margin:16px 0;padding:18px 22px}}
h2{{font-size:16px;margin:0 0 10px;border-left:4px solid #6A3D9A;padding-left:11px}}
img{{width:100%;height:auto;border-radius:8px;border:1px solid #e6e8ec;display:block;margin:8px 0}}
p,li{{font-size:13.7px;color:#33384a}} .cap{{font-size:12.5px;color:#6a7180}}
table{{border-collapse:collapse;width:100%;font-size:12.6px;margin:6px 0}} th,td{{border-bottom:1px solid #eceef2;padding:5px 8px;text-align:left}} th{{color:#6a7180}} td.n{{text-align:right;font-variant-numeric:tabular-nums}}
.synth{{background:#eef2f9;border:1px solid #aebbd6;border-radius:10px;padding:14px 18px;color:#16315e;font-size:14px}}
.go{{background:#eef6ee;border:1px solid #9ccc9c;border-radius:8px;padding:11px 15px;color:#1a4d1a;font-size:12.8px;margin-top:8px}}
.no{{background:#fdeaea;border:1px solid #e3a6a6;border-radius:8px;padding:11px 15px;color:#7a1a1a;font-size:12.8px;margin-top:8px}}
code{{background:#eef1f5;padding:1px 5px;border-radius:4px;font-size:12px}}
footer{{max-width:960px;margin:0 auto;padding:0 18px 50px;color:#8a93a3;font-size:11.5px}}
</style></head><body>
<header><h1>BAFF axis &mdash; receptor nCount-control + myeloid producer anchor</h1>
<p>Closes the BAFF story with the same rigor that retracted the stromal-ligand claim: test the receptors (BAFF-R, BCMA, TACI) for genuine aggregate upregulation, and anchor the myeloid producer. Neutral control genes (TPT1/PPIA/YWHAZ/TMSB10/UBB) measure generic transcript inflation in aggregate cells; nCount is the density proxy. APRIL is absent &rarr; BAFF is the relevant ligand throughout.</p></header>
<main>
<div class="card"><div class="synth"><b>Synthesis:</b> myeloid produce BAFF <b>tissue-wide</b> (robustly anchored), and B/plasma cells express the receptors <b>constitutively</b> &mdash; but there is <b>no aggregate-specific concentration of either ligand or receptor</b>. The axis is real but <b>not a spatially localized aggregate/TLS niche</b>. The global null was wrong about the producer, right about localization.</div></div>

<div class="card"><h2>PART A &mdash; receptors are NOT aggregate-upregulated beyond control</h2>
{img(f"{FIG}/baff_receptor_control.png",1500,88)}
<table><tr><th>receptor (target)</th><th class=n>det in %</th><th class=n>det out %</th><th class=n>in/out</th><th class=n>norm-expr</th><th class=n>control-gene</th><th>verdict</th></tr>{rrows}</table>
<p>Aggregate-resident B/plasma cells detect modestly more of <i>everything</i>; the receptor "rise" <b>matches the neutral control-gene inflation</b> (BAFF-R 1.29&times; vs control 1.28&times;; BCMA 1.30&times; vs 1.34&times;; TACI flat) and is inconsistent across sections. <b>The earlier "BAFF-R aggregate-elevated 1.6&times;" does not survive.</b> What's real: B cells constitutively express BAFF-R (27&ndash;37%) and plasma BCMA (12&ndash;17%) &mdash; just not concentrated at aggregates.</p>
<div class="no"><b>No localized receptor niche.</b></div></div>

<div class="card"><h2>PART B &mdash; myeloid BAFF is a reliable anchored producer</h2>
{img(f"{FIG}/baff_myeloid_anchor.png",1500,88)}
<div class="go"><b>(1) cell-intrinsic, not a count artifact:</b> myeloid BAFF per-transcriptome <b>{b1.normExpr_fold}&times;</b> epithelium &mdash; and myeloid have <i>fewer</i> counts (median {int(b1.myeloid_nCount_med)} vs {int(b1.epi_nCount_med)}), so the raw 7.9&times; understates it (nCount-matched ~24&times;). <b>(2) bona fide myeloid:</b> BAFF&plus; cells are <i>richer</i> in CD68/CD14/C1QA/AIF1 than BAFF&minus; (activated-macrophage phenotype). <b>(3) reproducible:</b> detected in all 16/16 samples (median {rep.BAFF_det_pct.median():.1f}%, range {rep.BAFF_det_pct.min():.1f}&ndash;{rep.BAFF_det_pct.max():.1f}%), not driven by 1&ndash;2 sections.</div></div>

<div class="card"><h2>Method</h2><p>Reuses the <code>dkd_xenium_disease</code> object (validated reannotation labels + per-cell receptor/ligand/nCount + control-gene panel + myeloid markers). Aggregate "inside" = DBSCAN B-lineage member (eps=50&micro;m/minPts=20); tested in the 4 aggregate-bearing sections (1006, HK2695, 1007, 1009). Folder <code>analysis/dkd_xenium_disease/</code>; <code>baff_receptor_anchor.py</code>, <code>REPORT_baff_receptor.md</code>. Caveat: receptors/ligand sparsely detected; pooled across sections.</p></div>
</main>
<footer>Self-contained &middot; closes the BAFF story (producer real &amp; anchored; no localized niche) &middot; raw data read-only, git-ignored.</footer>
</body></html>"""
open(f"{OUT}/baff_receptor_summary.html","w").write(HTML)
print(f"wrote baff_receptor_summary.html ({len(HTML)//1024} KB)")
