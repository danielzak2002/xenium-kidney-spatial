#!/usr/bin/env python
"""Refresh the collected emailable summaries.
Run:  conda run -n spatial python refresh_summaries.py
Regenerates every analysis's self-contained *_summary.html from its build script, copies them into
summaries/ under stable descriptive names, and rebuilds summaries/index.html from the list below.

GENERAL PRACTICE: after any analysis that adds or changes a summary, (1) add a row to SUMMARIES,
(2) run this script. The summaries/ folder is git-ignored (base64-embedded HTML stays out of the
public repo); this script + the per-analysis build_*summary*.py are the committed source of truth."""
import os, sys, shutil, subprocess, html as _html
REPO=os.path.dirname(os.path.abspath(__file__)); SUM=os.path.join(REPO,"summaries"); os.makedirs(SUM,exist_ok=True)

# build_script (rel to repo) | source html it writes (rel) | dest name in summaries/ | title | one-line desc | headline | start?
SUMMARIES=[
 ("analysis/build_dkd_combined_html.py","analysis/dkd_combined_summary.html",
  "00_DKD_combined_overview_START_HERE.html","DKD combined overview",
  "The full arc in one file: independent reannotation → B-rich subgroup + references → the BAFF axis, with a table of contents.",
  "▶ START HERE — the single-file story.",True),
 ("analysis/dkd_xenium_reannotation/build_summary_html.py","analysis/dkd_xenium_reannotation/dkd_reannotation_summary.html",
  "01_DKD_reannotation_validated_vs_authors.html","Independent reannotation, validated vs authors",
  "Blind full-panel re-typing of 16 Xenium samples, benchmarked against the authors' labels.",
  "Segment ARI 0.78 / immune ARI 0.68 — both typings validated.",False),
 ("analysis/dkd_xenium_disease/build_summary_html.py","analysis/dkd_xenium_disease/dkd_disease_summary.html",
  "02_DKD_Brich_subgroup_and_references.html","B-rich subgroup + cross-nephropathy references",
  "Within-DKD B-rich/B-poor split; IgAN, MN, AA-amyloid, C3GN as individual references; per-participant B-lineage gallery.",
  "B-rich = {1006, HK2695}, 100% concordant with the authors' B-predominant niche.",False),
 ("analysis/dkd_xenium_disease/build_baff_ambient_html.py","analysis/dkd_xenium_disease/baff_ambient_summary.html",
  "03_BAFF_stromal_ambient_stresstest.html","BAFF — peri-aggregate stromal stress-test",
  "Tests whether 'stromal BAFF near aggregates' is real production or ambient/spillover.",
  "Retracted — non-specific spillover field, not stromal production.",False),
 ("analysis/dkd_xenium_disease/build_baff_receptor_html.py","analysis/dkd_xenium_disease/baff_receptor_summary.html",
  "04_BAFF_receptors_and_myeloid_anchor.html","BAFF — receptor control + myeloid producer anchor",
  "Receptors (BAFF-R/BCMA/TACI) under nCount control; anchoring the myeloid BAFF producer.",
  "Myeloid BAFF anchored (~24× epithelium, 16/16); no aggregate-specific receptor/ligand niche.",False),
 ("analysis/dkd_epi_endo_stress/build_reconcile_summary_html.py","analysis/dkd_epi_endo_stress/reconcile_summary.html",
  "05_injuredPT_vs_Blineage_reconcile.html","Injured-PT ↔ B-lineage — reconcile under post-BAFF rigor",
  "Re-examines the prior tubular-injury-near-B-aggregates result with validated labels + control-gene/specificity rigor.",
  "Real co-localization, but not B-specific and not a program-intensity gradient (prior +0.13 retracted).",False),
 ("analysis/dkd_epi_endo_stress/build_colocalization_summary_html.py","analysis/dkd_epi_endo_stress/colocalization_summary.html",
  "06_injury_immune_colocalization.html","Injury ↔ immune co-localization across nephropathies",
  "Which cells participate: de-circularized specificity, cross-disease, epithelial + immune composition.",
  "Injury co-localizes with myeloid (ρ 0.82), not B-lineage; near-injury infiltrate is myeloid-skewed.",False),
 ("analysis/dkd_xenium_disease/build_bcell_anecdotal_html.py","analysis/dkd_xenium_disease/bcell_anecdotal_summary.html",
  "07_Blineage_mechanism_across_nephropathies_anecdotal.html","B-lineage mechanism across nephropathies (anecdotal)",
  "Extends 02 with B:Plasma split, Ig isotype, TLS organization, localization, damage coupling, cell state — n=1 per non-DKD, descriptive only.",
  "Three apparent programs: DKD B-rich organized-lymphoid · MN antibody-mediated/plasma-skewed · IgAN diffuse/unorganized (n=1, not tested).",False),
 ("analysis/dkd_xenium_disease/build_composition_html.py","analysis/dkd_xenium_disease/composition_by_group_summary.html",
  "08_composition_by_disease_group.html","Cell-type composition per sample by disease group",
  "Walkthrough top-level layer: per-sample fractions (coarse lineage + immune subtypes) on the validated labels; dotplots with CLR sensitivity.",
  "Control epithelial-dominant (74%) / near-bare immune (2.8%); every disease loses epithelium, gains immune+stroma; descriptive, n=1 flagged.",False),
]

def main():
    ok=[]
    for build,src,dst,*_ in SUMMARIES:
        bp=os.path.join(REPO,build)
        if not os.path.exists(bp): print(f"  [skip] missing build script: {build}"); continue
        r=subprocess.run([sys.executable,bp],capture_output=True,text=True)
        sp=os.path.join(REPO,src)
        if r.returncode!=0 or not os.path.exists(sp):
            print(f"  [FAIL] {build}\n{r.stderr[-500:]}"); continue
        shutil.copyfile(sp,os.path.join(SUM,dst)); ok.append(dst); print(f"  [ok] {dst}")
    # rebuild index
    cards=[]
    for _,_,dst,title,desc,head,start in SUMMARIES:
        if dst not in ok: continue
        cls="card start" if start else "card"
        cards.append(f'<a class="{cls}" href="{dst}"><span class="t">{_html.escape(title)}</span>'
                     f'<span class="f">{dst}</span><div class="d">{_html.escape(desc)}</div>'
                     f'<div class="v">{_html.escape(head)}</div></a>')
    idx=f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>DKD Xenium — summary index</title><style>
body{{margin:0;font-family:-apple-system,'Helvetica Neue',Arial,sans-serif;background:#f4f5f7;color:#1c1c1c;line-height:1.5}}
header{{background:#15203a;color:#fff;padding:22px 30px}} header h1{{margin:0;font-size:20px}} header p{{margin:8px 0 0;color:#b9c6e0;font-size:13.5px;max-width:860px}}
main{{max-width:920px;margin:0 auto;padding:20px 18px 70px}}
a.card{{display:block;text-decoration:none;background:#fff;border-radius:11px;box-shadow:0 1px 4px rgba(0,0,0,.1);margin:13px 0;padding:15px 20px;border-left:5px solid #6A3D9A}}
a.card.start{{border-left-color:#2ca02c;background:#f3faf3}}
.t{{font-size:15px;font-weight:700;color:#15203a}} .f{{font-size:11.5px;color:#8a93a3;font-family:ui-monospace,Menlo,monospace;margin-left:8px}}
.d{{font-size:13px;color:#33384a;margin:5px 0 0}} .v{{font-size:12.3px;color:#6A3D9A;margin-top:4px}}
footer{{max-width:920px;margin:0 auto;padding:0 18px 50px;color:#8a93a3;font-size:11.5px}}
</style></head><body>
<header><h1>DKD Xenium — emailable summary index</h1>
<p>Self-contained HTML summaries of the diabetic-kidney-disease spatial analysis arc. Each file is standalone — open or forward any one. Pure science; descriptive/underpowered where noted. Regenerate with <code>conda run -n spatial python refresh_summaries.py</code>.</p></header>
<main>{''.join(cards)}</main>
<footer>Folder: <code>summaries/</code> (git-ignored). Built by <code>refresh_summaries.py</code> from each analysis's <code>build_*summary*.py</code>. Pure science · raw data read-only.</footer>
</body></html>"""
    open(os.path.join(SUM,"index.html"),"w").write(idx)
    print(f"\nrefreshed {len(ok)}/{len(SUMMARIES)} summaries + index.html -> summaries/")
if __name__=="__main__": main()
