#!/usr/bin/env python
"""Self-contained, emailable HTML for the peri-aggregate stromal BAFF ambient/spillover stress-test."""
import os, io, base64, pandas as pd
from PIL import Image
OUT="/Users/danie/ClaudeCode/pilot_analyses/xenium/analysis/dkd_xenium_disease"; FIG=f"{OUT}/figures"
def img(path,max_w=1500,q=88):
    im=Image.open(path).convert("RGB")
    if im.width>max_w: im=im.resize((max_w,round(im.height*max_w/im.width)),Image.LANCZOS)
    b=io.BytesIO(); im.save(b,format="JPEG",quality=q,optimize=True)
    return f'<img src="data:image/jpeg;base64,{base64.b64encode(b.getvalue()).decode()}">'
s2=pd.read_csv(f"{OUT}/baff_ambient_step2.csv"); s3=pd.read_csv(f"{OUT}/baff_ambient_step3.csv")
def s3row(r):
    hl=' style="background:#fdeaea"' if r.cell_type in("Stromal","Epithelial (non-producer)","Endothelial") else ""
    return f"<tr{hl}><td>{r.cell_type}</td><td class=n>{r.baff_near_pct}</td><td class=n>{r.baff_far_pct}</td><td class=n><b>{r.baff_near_over_far}&times;</b></td><td class=n>{r.nCount_near_over_far}</td></tr>"
s3rows="".join(s3row(r) for _,r in s3.iterrows())
HTML=f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>BAFF peri-aggregate stress-test</title><style>
body{{margin:0;font-family:-apple-system,'Helvetica Neue',Arial,sans-serif;background:#f4f5f7;color:#1c1c1c;line-height:1.55}}
header{{background:#15203a;color:#fff;padding:22px 30px}} header h1{{margin:0;font-size:20px}} header p{{margin:8px 0 0;color:#b9c6e0;font-size:14px;max-width:900px}}
main{{max-width:960px;margin:0 auto;padding:22px 18px 70px}}
.card{{background:#fff;border-radius:12px;box-shadow:0 1px 4px rgba(0,0,0,.1);margin:16px 0;padding:18px 22px}}
h2{{font-size:16px;margin:0 0 10px;border-left:4px solid #b8860b;padding-left:11px}}
img{{width:100%;height:auto;border-radius:8px;border:1px solid #e6e8ec;display:block;margin:8px 0}}
p,li{{font-size:13.7px;color:#33384a}} .cap{{font-size:12.5px;color:#6a7180}}
table{{border-collapse:collapse;width:100%;font-size:12.7px;margin:6px 0}} th,td{{border-bottom:1px solid #eceef2;padding:5px 8px;text-align:left}} th{{color:#6a7180}} td.n{{text-align:right;font-variant-numeric:tabular-nums}}
.verdict{{background:#fdeaea;border:1px solid #e3a6a6;border-radius:10px;padding:14px 18px;color:#7a1a1a;font-size:14px}}
.keep{{background:#eef6ee;border:1px solid #9ccc9c;border-radius:8px;padding:11px 15px;color:#1a4d1a;font-size:12.8px;margin-top:8px}}
code{{background:#eef1f5;padding:1px 5px;border-radius:4px;font-size:12px}}
footer{{max-width:960px;margin:0 auto;padding:0 18px 50px;color:#8a93a3;font-size:11.5px}}
</style></head><body>
<header><h1>Peri-aggregate stromal BAFF &mdash; ambient/spillover stress-test</h1>
<p>The disease analysis flagged stromal BAFF "enriched ~2.1&ndash;2.5&times; near B-aggregates" with a TLS-niche reading. B-aggregates are dense regions, so a local spillover field would raise BAFF in <i>any</i> nearby cell. This tests whether the signal is stromal-<b>specific production</b> or shared <b>ambient/spillover</b>. (No Xenium negative-control probes in the object &rarr; epithelium = local non-producer ambient floor; nCount = density proxy.)</p></header>
<main>
<div class="card"><div class="verdict"><b>VERDICT: ambient/spillover &mdash; NOT localized stromal production.</b> The localized peri-aggregate stromal-BAFF claim is <b>retracted</b>. What survives: <b>myeloid tissue-wide BAFF</b> (the robust producer, 7.9&times; the epithelial floor) and the tissue-wide stroma&gt;epithelium baseline. B-aggregates do sit in a BAFF-elevated field, but it is spillover, not stromal production.</div></div>

<div class="card"><h2>The three controls</h2>{img(f"{FIG}/baff_ambient_control.png",1600,88)}
<p class="cap"><b>(1)</b> Stromal &gt; epithelial in BOTH near and far &rarr; a tissue-wide baseline, not a peri-aggregate effect. <b>(2)</b> Stromal BAFF rises ~2.4&times; near aggregates &mdash; looks like localization in isolation. <b>(3) Decisive:</b> the near/far rise is non-specific.</p></div>

<div class="card"><h2>Why it fails (STEP 3 &mdash; specificity)</h2>
<table><tr><th>cell type</th><th class=n>BAFF near %</th><th class=n>BAFF far %</th><th class=n>near/far</th><th class=n>nCount near/far</th></tr>{s3rows}</table>
<p>The near-aggregate BAFF rise is shared by <b>non-producer epithelium (2.15&times;) and endothelium (2.86&times;)</b> at the same magnitude as stromal (2.44&times;), while <b>transcript density is flat (~0.98&times;)</b> and the genuine tissue-wide producer <b>myeloid is flat (1.14&times;)</b>. A stromal-specific producer would rise <i>more</i> than bystander non-producers; instead all parenchymal/stromal cells in the neighborhood rise together &mdash; the signature of a local ambient/spillover field.</p>
<div class="keep"><b>Kept:</b> myeloid tissue-wide BAFF (real, dominant producer). <b>Open item:</b> the receptor-side BAFF-R aggregate-concentration (1.6&times;) was not re-tested and needs the same non-producer/near-far control before any "B cells receive BAFF in the niche" claim. Caveat: BAFF is rare (small positive counts; pooled over 4 sections).</div></div>

<div class="card"><h2>Method</h2><p>Reuses the <code>dkd_xenium_disease</code> object (validated reannotation labels + per-cell BAFF/TNFSF13B + nCount). Aggregate-bearing sections: B-rich DKD 1006, HK2695 + C3GN 1007 + AA amyloid 1009. Distance to nearest B-aggregate member (DBSCAN eps=50/minPts=20); near &le;50&micro;m, far &gt;200&micro;m. Folder <code>analysis/dkd_xenium_disease/</code>; <code>baff_ambient_control.py</code>, <code>REPORT_baff_ambient.md</code>.</p></div>
</main>
<footer>Self-contained &middot; honest negative result &middot; raw data read-only, git-ignored.</footer>
</body></html>"""
open(f"{OUT}/baff_ambient_summary.html","w").write(HTML)
print(f"wrote baff_ambient_summary.html ({len(HTML)//1024} KB)")
