#!/usr/bin/env python
"""Self-contained, emailable HTML summary of the squidpy neighbourhood-enrichment screen.
Embeds the two heatmap figures as downscaled JPEGs; numbers pulled from the classified CSV.
Named *_summary.html -> git-ignored (base64-embedded); send via the file tool, do not commit."""
import os, io, base64, pandas as pd
from PIL import Image
OUT="/Users/danie/ClaudeCode/pilot_analyses/xenium/analysis/dkd_xenium_disease"; FIG=f"{OUT}/figures"
def datauri(path,max_w=1700,q=85):
    im=Image.open(path).convert("RGB")
    if im.width>max_w: im=im.resize((max_w,round(im.height*max_w/im.width)),Image.LANCZOS)
    b=io.BytesIO(); im.save(b,format="JPEG",quality=q,optimize=True)
    return "data:image/jpeg;base64,"+base64.b64encode(b.getvalue()).decode()
def img(path,max_w=1700,q=85): return f'<img src="{datauri(path,max_w,q)}">'

cl=pd.read_csv(f"{OUT}/nhood_offtarget_classified.csv")
vc=cl[cl["class"]!="no reproducible DKD sign"]["class"].value_counts()
def n(k): return int(vc.get(k,0))
surv=cl[cl.disease_specific_survivor]
clsrows="".join(f"<tr><td>{k}</td><td class=n>{v}</td></tr>" for k,v in vc.items())

HTML=f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Unbiased neighbourhood-enrichment screen (squidpy)</title><style>
body{{margin:0;font-family:-apple-system,'Helvetica Neue',Arial,sans-serif;background:#f4f5f7;color:#1c1c1c;line-height:1.55}}
header{{background:#15203a;color:#fff;padding:22px 30px}} header h1{{margin:0;font-size:21px}}
header p{{margin:8px 0 0;color:#b9c6e0;font-size:14px;max-width:940px}}
main{{max-width:1000px;margin:0 auto;padding:22px 18px 70px}}
.card{{background:#fff;border-radius:12px;box-shadow:0 1px 4px rgba(0,0,0,.1);margin:16px 0;padding:18px 22px}}
h2{{font-size:17px;margin:0 0 10px;border-left:4px solid #2ca02c;padding-left:11px}}
h2 .s{{background:#2ca02c;color:#fff;border-radius:5px;padding:1px 8px;font-size:12px;margin-right:7px}}
img{{width:100%;height:auto;border-radius:8px;border:1px solid #e6e8ec;display:block;margin:8px 0}}
p,li{{font-size:13.7px;color:#33384a}} .cap{{font-size:12.5px;color:#6a7180;margin:4px 2px 0}}
table{{border-collapse:collapse;width:100%;font-size:12.7px;margin:6px 0}}
th,td{{border-bottom:1px solid #eceef2;padding:5px 9px;text-align:left}} th{{color:#6a7180}} td.n{{text-align:right;font-variant-numeric:tabular-nums}}
.null{{background:#eef6ee;border:1px solid #9ccc9c;border-radius:8px;padding:12px 16px;font-size:13.2px;color:#1a4d1a;margin:8px 0}}
.flag{{background:#fff7ed;border:1px solid #f1c899;border-radius:8px;padding:10px 14px;font-size:12.7px;color:#7a4a12;margin-top:8px}}
code{{background:#eef1f5;padding:1px 5px;border-radius:4px;font-size:12px}}
footer{{max-width:1000px;margin:0 auto;padding:0 18px 50px;color:#8a93a3;font-size:11.5px}}
</style></head><body>
<header><h1>Unbiased all-pairs neighbourhood-enrichment screen (squidpy)</h1>
<p>The hypothesis-FREE co-localization check on the DKD Xenium atlas. Per section: Delaunay graph pruned
to ≤50 µm, then <code>nhood_enrichment</code> over 21 validated subtype labels. Read by SIGN +
reproducibility, never z-magnitude. Dumoulin et&nbsp;al. 2026. Raw read-only.</p></header>
<main>

<div class="card"><h2><span class="s">RESULT</span>An honest null — no new disease-specific immune niche</h2>
<div class="null"><b>The screen recovers tissue ARCHITECTURE and reproduces the known biology, but
surfaces no new disease-specific off-target immune co-localization.</b> Of 201 off-target pairs, 146
carry a reproducible DKD sign — classified below. The Control→DKD shift map is mostly empty.</div>
<table><tr><th>class of reproducible off-target structure</th><th class=n>n pairs</th></tr>{clsrows}</table>
<div class="flag"><b>What the screen measures:</b> cell-TOUCH adjacency (do two types sit as graph
neighbours more than a label-shuffle expects) — a stricter, different question than section-level
co-abundance. Squidpy's z is computed against an abundance-preserving label-permutation null; we add an
explicit label-shuffle null and a spillover cross-check (analysis 06 flags).</div></div>

<div class="card"><h2><span class="s">FIG 1</span>Control vs DKD sign-reproducibility — the shift is mostly empty</h2>
{img(f"{FIG}/nhood_dkd_vs_control_sign.png")}
<p>Signed fraction of sections enriched (red) − avoided (blue), NOT z-magnitude. Control (left) and DKD
(middle) share the same architecture: immune block mutually enriched; immune × epithelium avoided
(compartmentalization); nephron/vascular segments enriched. <b>Right (DKD − Control) is mostly pale</b> —
the few non-white cells are immune rows where controls are too immune-sparse to register (power, not biology).</p></div>

<div class="card"><h2><span class="s">FIG 2</span>Every section shows the same block structure</h2>
{img(f"{FIG}/nhood_per_section_heatmaps.png")}
<p>Per-section z heatmaps (one per group + both DKD subgroups). The same architecture recurs everywhere —
SIGN is read, magnitude is not compared across sections. *spillover-flagged section.</p></div>

<div class="card"><h2><span class="s">DETAIL</span>What's real but not new, and the known biology reproducing</h2>
<ul>
<li><b>Constitutive (real, but in controls too → not disease-specific):</b> Myeloid×Fibroblast, Myeloid×CD4&nbsp;T,
Myeloid×CD8&nbsp;T, CD4&nbsp;T×CD8&nbsp;T — all enriched <b>8/8 DKD AND 3/3 control</b>.</li>
<li><b>Known reproduces:</b> B-aggregate composition (B×Plasma +7/8, B×CD4&nbsp;T +8/8, Plasma×Myeloid +8/8).</li>
<li><b>Injury×myeloid refined:</b> at cell-touch it is <i>avoidance</i> (iTAL×Myeloid −8/8) — NOT a
contradiction of 06's section-level ρ&nbsp;0.82, but its refinement: myeloid sits in the interstitium
<i>near</i> injured tubules, not in direct epithelial contact.</li>
<li><b>The 2 "disease-specific" survivors</b> (iPT×Fibroblast, iPT×EC_glom) are <b>non-immune avoidance
pairs dominated by iPT</b> (recall ≈ 0.64; ~2× abundance shift control→DKD) — most likely typing/abundance,
not a niche. No immune off-target niche survives.</li>
</ul></div>

<div class="card"><h2><span class="s">METHOD</span>How (one paragraph)</h2>
<p>squidpy 1.8 per physical section: <code>spatial_neighbors</code> (Delaunay, generic coords) pruned to ≤50 µm
(graph degree ≈5.9, within-tissue); <code>nhood_enrichment</code> (n_perms=1000) over 21 validated subtype
labels; z→sign (\|z\|>2); types &lt;20 cells masked. Unit = section, z never compared across sections —
SIGN + cross-section reproducibility instead. Abundance control = per-section label-shuffle null; spillover
control = analysis-06 flags (clean vs flagged DKD). Survivors followed up with <code>co_occurrence</code>.
Built from validated reannotation labels + coordinates only (8.7 GB h5ad never opened). Folder:
<code>analysis/dkd_xenium_disease/</code> · <code>nhood_enrichment_screen.py</code> + <code>REPORT_nhood_enrichment.md</code>.</p></div>

</main>
<footer>Self-contained summary · pure science · raw data read-only, git-ignored. <b>Honest null: unbiased
screen recovers architecture + known biology, no new disease-specific immune niche. Unit=section; SIGN not z.</b></footer>
</body></html>"""
path=f"{OUT}/nhood_enrichment_summary.html"
open(path,"w").write(HTML)
print(f"wrote {path} ({len(HTML)//1024} KB)")
