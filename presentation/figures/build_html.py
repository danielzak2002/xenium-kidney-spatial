#!/usr/bin/env python
"""build_html.py — self-contained scrollable HTML gallery of the presentation figures
(talk order, captions). Figures referenced relatively (HTML lives beside the PNGs/SVGs)."""
import os, html, base64
import figstyle as fs
OUT=fs.OUT

def _datauri(fname, max_w=None):
    p=os.path.join(OUT,fname)
    if not os.path.exists(p): return None
    if fname.endswith(".svg"):
        return "data:image/svg+xml;base64,"+base64.b64encode(open(p,"rb").read()).decode()
    if max_w:  # downscale PNG for an emailable self-contained file (screen-res, not print)
        import io
        from PIL import Image
        im=Image.open(p).convert("RGB")
        if im.width>max_w: im=im.resize((max_w,round(im.height*max_w/im.width)),Image.LANCZOS)
        buf=io.BytesIO(); im.save(buf,format="JPEG",quality=88,optimize=True)
        return "data:image/jpeg;base64,"+base64.b64encode(buf.getvalue()).decode()
    return "data:image/png;base64,"+base64.b64encode(open(p,"rb").read()).decode()

CAP={
 "A1":("Spatial cell-type maps — ccRCC, cLN, DKD","Cells at spatial coordinates colored by harmonized lineage (shared palette). Orients the audience in each tissue; immune cells drawn on top of parenchyma. cLN (CosMx) cropped to its densest core.",False),
 "A3":("CD4/CD8 subtype: measured on Xenium, imputed on CosMx","Discriminating-marker (CD8A/CD4) detection in CD8- vs CD4-labelled cells, faceted by platform, AUROC annotated. On CosMx the CD4 AUROC ≈ 0.50 — the subtype split is reference-imputed, not measured. Dotted line = ambient floor.",False),
 "B2":("cLN: ambient CD3 mis-assigned to epithelium","CD3-family detection/level by compartment vs the per-cell negmean floor. ~35% of epithelial cells are CD3+, only ~2.3× T-vs-epithelial separation, ~7.4× above ambient — bounds (does not fix) the contamination.",False),
 "B3":("Usability gate — markers clearing vs failing the ambient floor","Detection in expected cells vs ambient floor for every candidate marker. Green = usable, red = sub-ambient/fail. BAFF/APRIL ligands and the fibroEMT program fail; injured-PT and endothelial-activation markers pass.",False),
 "B1":("★ Aggregate marker overlays — the impact figure","One representative ccRCC and one DKD B/plasma aggregate (~190 µm field). Blue underlay = B-core; marker+ cells overlaid (MS4A1, MZB1, FOXP3, CD8A, GZMB) plus a composite B/Treg/CD8 panel. Treg-around & cytotoxic-excluded in ccRCC vs cytotoxic-mixed-in in DKD — visible, not just numeric.",True),
 "C1":("DBSCAN B-cell density aggregates over tissue","B-cell aggregates delineated per section (DBSCAN ε=50 µm, minPts=20) for ccRCC and DKD. Light = dispersed B cells, dark = aggregated B cells.",False),
 "C2":("Per-aggregate Treg vs effector-CD8 enrichment","Per-aggregate log₂ enrichment (inside vs section background) for Treg-like and effector-CD8, ccRCC vs DKD. Tick = median. ccRCC excludes cytotoxic CD8; DKD co-enriches it.",False),
 "C3":("★ Burden-corrected differential — THE headline","Δlog₂ = log₂(Treg) − log₂(effector-CD8) per aggregate, count-pooled with bootstrap 95% CI. ccRCC +2.60 [+2.32,+2.84] (~6× Treg bias) vs DKD +0.24 [−0.31,+0.72] (~1×, no bias); CIs non-overlapping → the immunoregulatory bias is tumor-specific. Immune to the cytotoxic-burden confound.",True),
 "C4":("Count-pooled radial profile — no Treg collar","Treg-like and effector-CD8 enrichment in concentric rings from the B-core (count-pooled). Treg is flat across rings; a mild cytotoxic-core gradient is the only spatial structure — the earlier 'margin collar' was a mean-of-log2 artifact.",False),
 "C5":("Comparative non-immune × immune neighborhood enrichment","Per-section neighborhood z across kidney contexts. Absolute z is dominated by immune↔parenchyma geometry (immune cells aggregate, parenchyma 'avoids'); the biology is in the differentials. NOTE: the ccRCC stroma–immune inversion is provisional (not tile-verified).",False),
 "D1":("Conserved B/plasma scaffold, context-specific wiring","Schematic: one shared B/plasma core; three surrounds — ccRCC (Treg ring, cytotoxic excluded), cLN (+ myeloid / plasma–myeloid niche), DKD (cytotoxic mixed-in + injured tubule). Programmatic draft — may want manual vector polish.",False),
 "D2":("Platform capability matrix","What each panel can and cannot establish: T-lineage (both measured), CD4/CD8 subtype (Xenium measured / CosMx imputed), BAFF/APRIL ligands (sub-ambient both), receptors BCMA/BAFF-R/TACI (specific both), ambient/segmentation (cLN ~35% epithelial CD3+), IF anchors (decisive).",False),
 "A2":("cLN CosMx immune typing benchmark","InSituType recall and precision vs author labels per immune type (n annotated shown).",False),
 "B4":("DKD injured-PT program is elevated near B-aggregates","Stress-program score near vs far from B-aggregates (matched within section/cell-type) and the distance gradient. Injured-PT: Δz +0.13, 6/9 sections, p=0.038, monotonic gradient. Endothelial-activation and hypoxia are flat.",False),
 "Q1":("Data quality per dataset — the substrate is usable","Median transcripts/cell, genes/cell, panel depth, and ambient signal (cLN negmean). Panel depth differs explicitly (RCC 405 · PRCC 377 · cLN 957 · DKD 5 443) — read downstream sensitivity against this. DKD/Xenium neg-probes were dropped from the release, so ambient is not computable there.",False),
 "T1":("★ Marker dot-plot — the cell types are defined by the right markers","Canonical lineage markers × harmonized cell types; dot size = detection rate, color = scaled mean expression. B (MS4A1/CD79A), plasma (MZB1/DERL3), T (CD3D/E), myeloid (CD68/LYZ), endothelial (PECAM1/VWF), epithelial (EPCAM/KRT8), plus Treg (FOXP3) and cytotoxic (CD8A/GZMB) where resolved. The marker→type diagonal IS the typing veracity.",True),
 "T2":("UMAP by lineage with key-marker overlays","Stored embeddings colored by harmonized lineage, with MS4A1/CD3D/CD68/EPCAM overlays. The type structure is real, not imposed — marker hotspots coincide with the lineage islands.",False),
 "A3raw":("Raw counts behind the AUROC","Per-cell CD8A and CD4 count distributions in CD4- vs CD8-labelled cells, faceted Xenium vs CosMx, ambient/PT floor marked. The distributions visibly separate on Xenium and overlap on CosMx — the AUROC merely summarizes what you can see.",False),
 "C1raw":("Aggregates are real structures, not an algorithm artifact","Representative ccRCC and DKD sections: all cells grey, B cells highlighted, DBSCAN convex hulls overlaid; with an aggregate-size histogram and a nearest-neighbour test. Observed B–B spacing is far tighter than a permuted null — B cells genuinely cluster.",False),
 "C2raw":("Composition behind the enrichment — countable cells","Raw cells: inside-vs-background composition (ccRCC), per-aggregate Treg vs effector-CD8 counts (ccRCC above the diagonal = Treg>CD8; DKD on it), and per-aggregate Treg-enrichment consistency. The log₂ enrichment and the differential trace to real, countable cells.",False),
 "V1_qc":("ccRCC — data quality, two Xenium sections","Median transcripts/cell, genes/cell, panel size, and neg-control fraction for the discovery (ccRCC) and replication (pRCC) Xenium sections. Clean ambient (~0). The substrate is usable.",False),
 "V1_typing":("ccRCC — native cell typing","Native-label dot-plot (the study's own vocabulary, not harmonized) — each type is defined by its canonical markers — plus per-section UMAPs. The structure is visible in both sections.",False),
 "V1_gallery":("ccRCC — marker & population gallery (both sections)","Spatial maps of MS4A1/MZB1/FOXP3/CD8A/CD68 and the native B/plasma/Treg/CD8/myeloid populations in situ, ccRCC above pRCC. The biology reproduces across BOTH Xenium sections — the within-study veracity.",False),
 "V1_pattern":("ccRCC — observed immunoregulatory pattern","Per-aggregate Treg vs effector-CD8 counts (Treg>CD8 in nearly all aggregates) and the inside-vs-background composition: Treg enriched, cytotoxic-CD8 excluded. The native ccRCC signature, before any harmonization.",False),
 "V2_qc":("cLN — data quality across 14 slides","Median transcripts/genes per cell, cells per slide (4 control + 10 SLE), CosMx ambient (~0.10). Caveat stated up front: ~35% epithelial CD3 contamination → T-lineage is unreliable here and excluded from T claims.",False),
 "V2_typing":("cLN — native cell typing","Native author-label dot-plot + UMAP (CosMx 957-plex). See A2 for the InSituType benchmark vs author labels.",False),
 "V2_gallery":("cLN — marker & population gallery across slides","B/plasma/myeloid markers and populations in situ across 3 control + 3 SLE slides — the strongest reproducibility evidence in this study (a representative condition-spanning subset of the 14 slides).",False),
 "V2_niche":("cLN — observed plasma–myeloid niche","The reproducible cLN finding: plasma aggregates recruit myeloid cells, across slides and SLE classes (control-vs-SLE contrast). Plasma–myeloid neighborhood z and per-slide plasma-aggregate counts.",False),
 "V3_qc":("DKD — data quality, both platforms","Median transcripts/genes per cell and section counts for CosMx (n=48) and Xenium (n=16). Neg-probes were dropped from the release, so ambient is not computable.",False),
 "V3_typing":("DKD — native cell typing, both platforms","Native dot-plots for CosMx and Xenium (each platform's own labels) + UMAP. Subtype (CD4/CD8) is only trustworthy on Xenium — see A3/A3raw.",False),
 "V3_gallery":("DKD — marker & population gallery, both platforms","B/plasma/myeloid markers and populations in situ across representative CosMx and Xenium samples (can't show all 48 CosMx — a representative panel). The structural B/plasma finding is present on both platforms.",False),
 "V3_aggregates":("DKD — B/plasma aggregate composition on BOTH platforms","DBSCAN B-aggregate composition recovers on CosMx (lineage, new analysis) and Xenium — the structural finding is cross-platform (mirror of the ccRCC two-section replication). Subtype-resolved Treg/effector-CD8 only on Xenium.",False),
 "N_map":("Reconciling three native vocabularies → one common scheme","The native→common mapping table, shown as a deliberate step AFTER each study stands on its own. cLN T-lineage pooled and flagged unreliable; ccRCC epithelium is malignant (kept separate).",False),
}
BLOCKS=[("V1 · ccRCC vignette (Xenium ×2: discovery + replication)",["V1_qc","V1_typing","V1_gallery","V1_pattern"]),
        ("V2 · cLN vignette (CosMx, 14 slides)",["V2_qc","V2_typing","A2","V2_gallery","V2_niche","B2"]),
        ("V3 · DKD vignette (CosMx + Xenium)",["V3_qc","V3_typing","A3","A3raw","V3_gallery","V3_aggregates","C4","B4","B3"]),
        ("N · Normalization & harmonization (the deliberate bridge)",["N_map","T1","A1","T2"]),
        ("X · Cross-context — the earned comparison (across disease settings)",["B1","C3","C5","D1"]),
        ("P · Platforms & caveats",["D2"])]
EXTRA=["Q1","C1","C1raw","C2","C2raw"]
DS={k:v for k,v in fs.DATASET.items()}

def card(fid, embed):
    title,cap,hero=CAP[fid]
    star=' <span class="hero">★ hero</span>' if hero else ''
    if embed:  # embed downscaled JPEG once (no anchor -> no duplicate data-URI); emailable
        imgblock=f'<img src="{_datauri(f"{fid}.png", max_w=1500)}" alt="{fid}">'; links=""
    else:      # reference full 300-DPI PNG, click-through to vector SVG
        imgblock=f'<a href="{fid}.svg" target="_blank" title="open vector SVG"><img src="{fid}.png" alt="{fid}"></a>'
        links=f' <span class="links">[<a href="{fid}.png" target="_blank">PNG</a> · <a href="{fid}.svg" target="_blank">SVG</a>]</span>'
    return f'''<div class="card{' herocard' if hero else ''}" id="{fid}">
  <div class="meta"><span class="badge">{fid}</span><span class="ttl">{html.escape(title)}</span>{star}</div>
  {imgblock}
  <div class="cap">{html.escape(cap)}{links}</div>
</div>'''

def build(embed):
  sections=[]
  for name,ids in BLOCKS:
    cards="\n".join(card(i,embed) for i in ids)
    sections.append(f'<section><h2>{html.escape(name)}</h2>\n{cards}\n</section>')
  extra_cards="\n".join(card(i,embed) for i in EXTRA)
  sections.append(f'<section><h2>Supporting — global QC overview & cross-study aggregate evidence</h2>\n{extra_cards}\n</section>')
  contact=_datauri("CONTACT_SHEET.png", max_w=1800) if embed else "CONTACT_SHEET.png"
  foot=(f'<a href="{contact}" target="_blank">contact sheet</a>' if embed
        else '<a href="CONTACT_SHEET.png" target="_blank">contact sheet</a> · <a href="INDEX.md" target="_blank">index</a>')
  return f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Spatial-kidney — presentation figures</title>
<style>
:root{{--rcc:{DS['RCC']};--cln:{DS['cLN']};--dkd:{DS['DKD']};}}
*{{box-sizing:border-box}}
body{{margin:0;font-family:-apple-system,'Helvetica Neue',Arial,sans-serif;color:#1c1c1c;background:#f4f5f7;line-height:1.5}}
header{{position:sticky;top:0;z-index:10;background:#1f2430;color:#fff;padding:16px 28px;box-shadow:0 2px 8px rgba(0,0,0,.25)}}
header h1{{margin:0;font-size:20px}}
header p{{margin:6px 0 0;font-size:13px;color:#c9d2e0;max-width:1000px}}
nav{{margin-top:10px;font-size:12px}}
nav a{{color:#8fb7ff;text-decoration:none;margin-right:12px}}
main{{max-width:1180px;margin:0 auto;padding:24px 20px 80px}}
section{{margin:34px 0}}
h2{{font-size:16px;color:#3a4252;border-left:4px solid var(--rcc);padding-left:10px;margin:0 0 14px}}
.card{{background:#fff;border-radius:12px;box-shadow:0 1px 4px rgba(0,0,0,.10);padding:16px 18px;margin:16px 0}}
.herocard{{box-shadow:0 0 0 2px #f0b429,0 2px 10px rgba(0,0,0,.15)}}
.meta{{display:flex;align-items:center;gap:10px;margin-bottom:10px;flex-wrap:wrap}}
.badge{{background:#1f2430;color:#fff;font-weight:700;font-size:13px;border-radius:6px;padding:2px 9px;letter-spacing:.5px}}
.ttl{{font-weight:700;font-size:16px}}
.hero{{background:#f0b429;color:#1f2430;font-weight:700;font-size:11px;border-radius:5px;padding:2px 7px}}
.card img{{width:100%;height:auto;border-radius:8px;border:1px solid #e6e8ec;background:#fff;display:block}}
.cap{{font-size:13.5px;color:#3a4252;margin-top:10px}}
.links{{color:#8a93a3;font-size:12px;white-space:nowrap}}
.links a{{color:#5a78b0;text-decoration:none}}
.takehome{{background:#fff7e6;border:1px solid #f0b429;border-radius:10px;padding:12px 16px;font-size:13.5px;margin:18px 0}}
.caveats{{background:#fff;border-radius:10px;padding:14px 18px;font-size:13px;color:#555;border:1px solid #e6e8ec}}
.caveats li{{margin:4px 0}}
footer{{max-width:1180px;margin:0 auto;padding:0 20px 60px;color:#8a93a3;font-size:12px}}
</style></head><body>
<header>
  <h1>Spatial transcriptomics of kidney B/plasma niches — presentation figures</h1>
  <p>Each study stands on its own first — QC → native cell typing → observed pattern → reproducibility across many sections — BEFORE any cross-dataset normalization. Harmonization (N) and the cross-context comparison (X) come only after, as deliberate steps. The earned conclusion: a conserved B/plasma scaffold acquires context-specific immune wiring — immunoregulatory (Treg-in / cytotoxic-out) in tumor, absent in non-malignant disease.</p>
  <nav>{' '.join(f'<a href="#{i}">{i}</a>' for _,ids in BLOCKS for i in ids)} {' '.join(f'<a href="#{i}">{i}</a>' for i in EXTRA)}</nav>
</header>
<main>
  <div class="takehome"><b>How to read these:</b> ★ = hero (B1 the impact figure, C3 the headline). Click any figure to open it full size. Order follows the talk.</div>
  {''.join(sections)}
  <div class="caveats"><b>Caveats baked into the figures</b>
  <ul>
    <li>cLN T-lineage is unreliable (~35% epithelial CD3 contamination) — excluded from conserved-T claims.</li>
    <li>The ccRCC stroma–immune inversion (C5) is provisional — not tile-verified.</li>
    <li>B-aggregate niches are associational / colocalization, not communication or causation.</li>
    <li>No patient column for DKD (donor clustering uncontrolled); ccRCC epithelium is malignant (interpreted separately).</li>
    <li>D1 is a programmatic schematic — a starting point for manual vector polish.</li>
  </ul></div>
</main>
<footer>Generated from presentation/figures/ · PNG@300 + SVG per figure · pure science (scientific labels only).
&nbsp;{foot}</footer>
</body></html>'''

h_rel=build(embed=False)
open(os.path.join(OUT,"gallery.html"),"w").write(h_rel)
print(f"wrote gallery.html ({len(h_rel)//1024} KB, references {len(CAP)} figures relatively)")
h_emb=build(embed=True)
open(os.path.join(OUT,"gallery_standalone.html"),"w").write(h_emb)
print(f"wrote gallery_standalone.html ({len(h_emb)//1024//1024} MB, {len(CAP)} figures embedded as base64 — single portable file)")
