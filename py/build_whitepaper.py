#!/usr/bin/env python
"""
build_whitepaper.py — assemble the self-contained HTML whitepaper from the corrected
narrative markdown. Figures are base64-embedded (portable single file); tables are
rendered from the committed CSVs; the cLN per-slide set is a scrollable gallery; the
per-dataset technical reports are linked beneath. No fabricated numbers — every figure
and table is a committed artifact.

  conda run -n spatial python py/build_whitepaper.py   (or any python3)
"""
import os, re, base64, html

ROOT = os.environ.get("XENIUM_ROOT", os.getcwd())
FIG = os.path.join(ROOT, "outputs/figures/whitepaper"); TAB = os.path.join(ROOT, "outputs/tables")
MD = os.path.join(ROOT, "spatial_kidney_report_narrative.md")
OUT = os.path.join(ROOT, "outputs/whitepaper.html")

def b64(fname):
    p = os.path.join(FIG, fname)
    with open(p, "rb") as fh:
        return base64.b64encode(fh.read()).decode()

def figure(fname, caption=""):
    if not os.path.exists(os.path.join(FIG, fname)):
        return f'<div class="missing">MISSING FIGURE: {html.escape(fname)}</div>'
    cap = f'<figcaption>{caption}</figcaption>' if caption else ''
    return f'<figure><img alt="{html.escape(fname)}" src="data:image/png;base64,{b64(fname)}">{cap}</figure>'

def figures(fnames, caption=""):
    return '<div class="figrow">' + "".join(figure(f) for f in fnames) + '</div>' + \
           (f'<p class="cap">{caption}</p>' if caption else '')

def gallery(fnames, caption=""):
    items = "".join(f'<div class="gitem">{figure(f)}</div>' for f in fnames)
    return (f'<p class="cap"><b>Per-slide gallery</b> — {caption} (scroll horizontally)</p>'
            f'<div class="gallery">{items}</div>')

def csv_table(fname, caption="", maxrows=40):
    import csv
    p = os.path.join(TAB, fname)
    if not os.path.exists(p):
        return f'<div class="missing">MISSING TABLE: {html.escape(fname)}</div>'
    with open(p) as fh:
        rows = list(csv.reader(fh))
    head, body = rows[0], rows[1:maxrows + 1]
    def fmt(x):
        try:
            f = float(x)
            return f"{f:.3g}" if (f != int(f)) else str(int(f))
        except (ValueError, TypeError):
            return html.escape(x)
    th = "".join(f"<th>{html.escape(c)}</th>" for c in head)
    trs = "".join("<tr>" + "".join(f"<td>{fmt(c)}</td>" for c in r) + "</tr>" for r in body)
    cap = f'<p class="cap">{caption}</p>' if caption else ''
    return f'<div class="tblwrap"><table><thead><tr>{th}</tr></thead><tbody>{trs}</tbody></table></div>{cap}'

# ---- marker -> HTML mapping -------------------------------------------------
def slot(inner):
    t = inner.lower()
    M = lambda *ks: all(k in t for k in ks)
    if M("counts-per-cell"):
        return figures(["qcA_kidney_RCC_protein_dist.png", "qcA_kidney_preview_PRCC_dist.png",
                        "qcA_cln_cosmx_dist.png"],
                       "Per-cell QC distributions (counts, genes, ambient, cell area) for each dataset.")
    if M("genes-per-cell") or M("negative-probe", "distribution") or M("cell-area distribution"):
        return ""  # shown in the per-dataset distribution panels above
    if M("per-dataset qc flag summary"):
        return csv_table("qcA_flag_summary.csv", "Per-dataset QC flag summary (flag-don't-filter; only zero-count cells dropped).")
    if M("comparative qc panel"):
        return figure("qcA_comparative.png", "Comparative QC: counts, genes, and unit-matched neg-control fraction (cLN ≈ 65× RCC-Xenium).")
    if M("clustering umap"):
        return figure("qcB_umap_lineage.png", "Clustering UMAP per dataset, colored by major lineage.")
    if M("harmony integration"):
        return figure("qcB_harmony_cln.png", "cLN Harmony integration: slides mix after integration; lineage structure preserved.")
    if M("assignment-confidence") or M("posterior distribution"):
        return figure("qcB_insitutype_posterior.png", "InSituType per-cell assignment confidence (81% of cells ≥ 0.90).")
    if M("immune recall before"):
        return figure("qcB_insitutype_recall_precision.png", "Immune recall AND precision: two-stage vs InSituType.")
    if M("spatial neighbor graph"):
        return figure("qcB_spatial_graph_crop.png", "Example squidpy Delaunay neighbor graph on an RCC tissue crop.")
    if M("three-way benchmark"):
        return figure("qcC_threeway_benchmark.png", "Three labelings vs author annotation: clustering 82% / marker 74% / reference-transfer 50%.")
    if M("confusion of reference-transfer"):
        return figure("qcC_reftransfer_confusion.png", "Reference-transfer vs author: immune dumped into Myeloid; epithelial leakage; no B/Plasma/T/NK calls.")
    if M("ambient contamination illustration"):
        return figure("qcC_ambient_cd3e.png", "Ambient: CD3E detected in 21% of epithelial/tubular cells, not only T cells.")
    if M("where author-immune cells land"):
        return ""  # captured by the reference-transfer confusion above
    if M("insitutype per-type recall/precision"):  # the TABLE
        return csv_table("cln_cosmx_immune_benchmark_unified.csv", "InSituType per-type recall and precision vs two-stage (single source of truth).")
    if M("per-type recall/precision, two-stage"):
        return ""  # the recall/precision figure is shown in §4.2
    if M("t-cell confusion matrix"):
        return figure("qcC_tcell_confusion.png", "T-cell loss: author CD8 28% kept; CD4/Treg 11% (never as CD4); rest to de-novo/epithelial.")
    if M("gene-usability gate"):
        return figure("qcC_gene_usability_gate.png", "Gene-usability gate: B/plasma/CD68 markers pass; BAFF ligand and chemokines fail (≤2× ambient).")
    if M("rcc neighborhood-enrichment z-score heatmap"):
        return figure("qcD_nhood_heatmap.png", "RCC neighborhood-enrichment z: B×Treg +28 (enriched), B×effector-CD8 −117 (excluded), boxed.")
    if M("rcc delineated aggregates"):
        return figure("qcD_aggregates_overview.png", "37 delineated RCC B-cell aggregates, cells colored by type.") + \
               figure("qcF_rcc_immune_aggregate_region.png", "RCC immunoregulatory aggregate (per-region, E-format): B–Treg core, effector-CD8 excluded.")
    if M("representative aggregate crop"):
        return figure("qcD_aggregate_markers.png", "Representative aggregate (425 B / 316 Treg / 232 CD8 / 9 mregDC) beside MS4A1/FOXP3/LAMP3/CD8A transcript overlays.")
    if M("mregdc", "foci"):
        return figure("qcD_mregdc_ccr7_foci.png", "Representative mregDC–CCR7⁺T focus (1 of ~38 across the section).")
    if M("big vs preview replication summary"):
        return csv_table("preview_phaseB_replication_summary.csv", "RCC (BIG) vs PRCC-preview replication summary.")
    if M("side-by-side big vs preview"):
        return figure("qcD_big_vs_preview_comp.png", "B-aggregate composition: Treg enriched and effector-CD8 excluded replicate; mregDC/plasma underpowered in preview.")
    if M("across-slide summary"):
        return figure("qcE_across_slide_summary.png", "cLN per-slide points: aggregate-based metrics (right) are disease-only; nhood z (left) is supporting only.")
    if M("per-slide plasma count"):
        return csv_table("cln_phaseB_per_slide.csv", "cLN per-slide plasma–myeloid metrics (all 14 slides). nhood z is the per-section (cross-core-pruned) value.")
    if M("rcc vs cln plasma"):
        return csv_table("plasma_myeloid_rcc_vs_cln.csv", "RCC vs cLN plasma–myeloid, identical metrics and 50 µm scale.")
    if M("side-by-side representative plasma foci"):
        return figure("qcF_contrast_plasma_foci.png", "Contrast (annotated): RCC myeloid not enriched in plasma aggregates (log2 −0.05); cLN enriched (+1.64).") + \
               figure("qcF_rcc_plasma_myeloid_region.png", "RCC plasma–myeloid (per-region, E-format): plasma aggregate, myeloid at background (not enriched).")
    if M("myeloid-state log2 enrichment"):
        return figure("qcF_myeloid_state_in_vs_out.png", "Myeloid-state log2 in/out of plasma niche, per gene, per-slide points (weak C1QB lean).") + \
               figure("qcE_cln_marker_crop.png", "cLN plasma–myeloid focus with gate-PASSING transcripts (MZB1=plasma, CD68/C1QB=myeloid).")
    if M("usability gate result"):
        return figure("qcF_usability_gate_recruitment.png", "Usability gate: structural macrophage markers pass; all secreted/inducible recruitment genes fail.")
    # cLN per-slide breakout (three markers): emit gallery once, drop the others
    if M("plasma + myeloid spatial scatter"):
        return gallery(["qcE_slide_SP21_213_R1080_S1.png", "qcE_slide_SP19_1139_R1080_S3.png",
                        "qcE_slide_SP20_642_R1080_S3.png", "qcE_slide_SP20_10838_R1080_S1.png",
                        "qcE_slide_SP19_4061_R1087_S1.png", "qcE_slide_SP18_8471_R1087_S2.png",
                        "qcE_slide_SMI0016C_SP17SP19.png"],
                       "five plasma-bearing slides (faceted per tissue core) + one control")
    if M("delineated plasma aggregates with myeloid") or M("representative high-magnification crop"):
        return ""
    return ""  # unknown marker -> drop

MARK = re.compile(r"\[(?:FIGURE|TABLE)[^\]]*\]")

# ---- minimal, safe markdown -> HTML (headers/lists/bold/italic/code/hr) -----
def inline(s):
    s = html.escape(s)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", s)
    s = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<em>\1</em>", s)
    return s

def md_to_html(md):
    out = []
    blocks = re.split(r"\n\s*\n", md)
    for blk in blocks:
        # collect + strip figure/table markers in this block
        figs = "".join(slot(m.group(0)[1:-1]) for m in MARK.finditer(blk))
        blk = MARK.sub("", blk)
        lines = [l for l in blk.split("\n")]
        stripped = blk.strip()
        if not stripped and not figs:
            continue
        if stripped.startswith("### "):
            out.append(f"<h3>{inline(stripped[4:].strip())}</h3>")
        elif stripped.startswith("## "):
            out.append(f"<h2>{inline(stripped[3:].strip())}</h2>")
        elif stripped.startswith("# "):
            out.append(f"<h1>{inline(stripped[2:].strip())}</h1>")
        elif re.match(r"^-{3,}$", stripped):
            out.append("<hr>")
        elif all(l.strip().startswith("- ") or not l.strip() for l in lines) and stripped.startswith("- "):
            items = "".join(f"<li>{inline(l.strip()[2:])}</li>" for l in lines if l.strip().startswith("- "))
            out.append(f"<ul>{items}</ul>")
        elif stripped:
            out.append(f"<p>{inline(' '.join(l.strip() for l in lines if l.strip()))}</p>")
        if figs:
            out.append(figs)
    return "\n".join(out)

with open(MD) as fh:
    body_html = md_to_html(fh.read())

REPORTS = [("RCC (Xenium) — Phase-A report", "kidney_RCC_protein_report.html"),
           ("PRCC preview (Xenium) — Phase-A report", "kidney_preview_PRCC_report.html"),
           ("cLN (CosMx) — Phase-A report", "cln_cosmx_report.html")]
report_links = "".join(
    f'<li><a href="{fn}">{html.escape(name)}</a></li>' for name, fn in REPORTS
    if os.path.exists(os.path.join(ROOT, "outputs", fn)))

CSS = """
:root{--ink:#1a1a1a;--muted:#555;--line:#e2e2e2;--accent:#4c72b0}
*{box-sizing:border-box}
body{font:16px/1.65 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:var(--ink);
 max-width:980px;margin:0 auto;padding:32px 22px 80px}
h1{font-size:1.9em;line-height:1.25;border-bottom:3px solid var(--accent);padding-bottom:.3em}
h2{font-size:1.45em;margin-top:2em;border-bottom:1px solid var(--line);padding-bottom:.2em}
h3{font-size:1.18em;margin-top:1.6em;color:#222}
p{margin:.7em 0}code{background:#f3f3f3;padding:.1em .35em;border-radius:4px;font-size:.9em}
hr{border:0;border-top:1px solid var(--line);margin:2em 0}
ul{margin:.6em 0 .6em 1.2em}li{margin:.25em 0}
figure{margin:1.2em 0;text-align:center}
img{max-width:100%;height:auto;border:1px solid var(--line);border-radius:6px;background:#fff}
figcaption,.cap{font-size:.86em;color:var(--muted);margin:.4em auto 0;max-width:90%}
.figrow{display:flex;flex-wrap:wrap;gap:10px;justify-content:center}
.figrow figure{flex:1 1 300px;margin:.4em 0}
.gallery{display:flex;overflow-x:auto;gap:14px;padding:10px 4px 16px;border:1px solid var(--line);
 border-radius:8px;background:#fafafa;scroll-snap-type:x mandatory}
.gitem{flex:0 0 86%;scroll-snap-align:start}.gitem img{border:none}
.tblwrap{overflow-x:auto;margin:1em 0}
table{border-collapse:collapse;font-size:.84em;width:100%}
th,td{border:1px solid var(--line);padding:5px 9px;text-align:center}
th{background:var(--accent);color:#fff;font-weight:600}
tr:nth-child(even) td{background:#f7f7f7}
.missing{color:#c44e52;font-weight:bold;border:2px dashed #c44e52;padding:8px;border-radius:6px}
.reports{margin-top:2.5em;padding:18px 22px;background:#f4f6fa;border:1px solid var(--line);border-radius:8px}
.foot{margin-top:3em;color:var(--muted);font-size:.85em;border-top:1px solid var(--line);padding-top:1em}
"""

doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Spatial Transcriptomics of B-Lineage Immune Organization in Human Kidney</title>
<style>{CSS}</style></head><body>
{body_html}
<div class="reports"><h2 style="margin-top:0;border:none">Per-dataset technical reports</h2>
<p>Detailed Phase-A QC / clustering / annotation reports for each dataset:</p>
<ul>{report_links}</ul></div>
<p class="foot">Self-contained technical whitepaper. Figures are embedded; all figures and tables are
committed artifacts under <code>outputs/figures/whitepaper/</code> and <code>outputs/tables/</code>.
Data: 10x Genomics (Xenium) and the published childhood lupus-nephritis CosMx atlas, CC BY 4.0. Code: MIT.</p>
</body></html>"""

with open(OUT, "w") as fh:
    fh.write(doc)
n_missing = doc.count("MISSING")
print(f"wrote {OUT}  ({len(doc)//1024} KB; {doc.count('<figure>')} figures, "
      f"{doc.count('<table>')} tables, {n_missing} MISSING)")
