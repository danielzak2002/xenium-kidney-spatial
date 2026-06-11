#!/usr/bin/env python
"""
figstyle.py — shared slide-grade style for the spatial-kidney presentation figures.
Every figure imports this. Pure science: scientific labels only. Colorblind-safe palettes
reused from the analysis scripts (dataset identity + lineage). 300 DPI, PNG + SVG.
"""
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
OUT=os.path.join(REPO,"presentation/figures"); os.makedirs(OUT,exist_ok=True)

# ---- locked DATASET identity palette (reused: RCC blue / PRCC green / cLN red / DKD purple) ----
DATASET={"RCC":"#1F78B4","PRCC":"#33A02C","cLN":"#E31A1C","DKD":"#6A3D9A"}
DATASET_LONG={"RCC":"ccRCC (Xenium)","PRCC":"pRCC (Xenium)","cLN":"cLN (CosMx)","DKD":"DKD (Xenium)"}

# ---- lineage palette (common harmonized scheme; colorblind-safe) ----
LINEAGE={
 "B":"#1f77b4","Plasma":"#ff7f0e","Myeloid":"#9467bd","T_lineage":"#2ca02c","NK":"#17becf",
 "Tubular_epi":"#c7c7c7","Endothelial":"#8c564b","Stroma":"#e377c2","Podocyte":"#7f7f7f",
 "Malignant_epi":"#b2182b"}
LINEAGE_LABEL={"B":"B cell","Plasma":"Plasma","Myeloid":"Myeloid","T_lineage":"T lineage",
 "NK":"NK","Tubular_epi":"Tubular epithelium","Endothelial":"Endothelial","Stroma":"Stroma",
 "Podocyte":"Podocyte","Malignant_epi":"Malignant epithelium"}
# malignant epithelium gets a consistent DISTINCT marker style wherever it appears
MALIGNANT_STYLE=dict(marker="*", edgecolor="black", linewidth=0.2)

# ---- subtype/marker accents (used in aggregate overlays / swarms) ----
ACCENT={"Treg":"#d62728","CD8eff":"#2ca02c","cytotoxic":"#2ca02c","B":"#1f77b4",
        "Plasma":"#ff7f0e","mregDC":"#9467bd","injured":"#8c510a"}
# per-marker sequential colormaps for intensity renders (consistent across B1)
MARKER_CMAP={"MS4A1":"Blues","MZB1":"Oranges","FOXP3":"Reds","CD8A":"Greens",
             "GZMB":"YlGn","CD3E":"Purples"}
MARKER_LINEAGE={"MS4A1":"B cell","MZB1":"Plasma","FOXP3":"Treg","CD8A":"Cytotoxic CD8",
                "GZMB":"Cytotoxic (GZMB)","CD3E":"T lineage"}

# ---- slide-grade rcParams ----
def apply():
    plt.rcParams.update({
        "figure.dpi":120, "savefig.dpi":300, "figure.facecolor":"white", "savefig.facecolor":"white",
        "font.family":"sans-serif",
        "font.sans-serif":["Helvetica Neue","Helvetica","Arial","DejaVu Sans"],
        "axes.titlesize":18, "axes.labelsize":14, "xtick.labelsize":12, "ytick.labelsize":12,
        "legend.fontsize":11, "axes.titleweight":"bold",
        "axes.spines.top":False, "axes.spines.right":False,
        "axes.linewidth":1.0, "axes.edgecolor":"#444444",
        "figure.titlesize":19, "figure.titleweight":"bold",
        "svg.fonttype":"none",
    })
apply()

# ---- helpers ----
def save_fig(fig, fid, tight=True):
    if tight:
        try: fig.tight_layout()
        except Exception: pass
    png=os.path.join(OUT,f"{fid}.png"); svg=os.path.join(OUT,f"{fid}.svg")
    fig.savefig(png, dpi=300, bbox_inches="tight"); fig.savefig(svg, bbox_inches="tight")
    plt.close(fig); print(f"  [ok] {fid}  -> {fid}.png + {fid}.svg")
    return fid

def zeroline(ax, value=0.0, orient="v", label="no bias / chance"):
    """Every enrichment plot gets a zero/chance reference line."""
    if orient=="v": ax.axvline(value, color="#888888", ls="--", lw=1.2, zorder=0)
    else: ax.axhline(value, color="#888888", ls="--", lw=1.2, zorder=0)
    return ax

def lineage_legend(ax, keys, loc="center left", bbox=(1.01,0.5), title="Lineage", ncol=1):
    handles=[]
    for k in keys:
        mk=MALIGNANT_STYLE.get("marker","o") if k=="Malignant_epi" else "o"
        handles.append(Line2D([],[],marker=mk,ls="",mfc=LINEAGE[k],
                      mec="black" if k=="Malignant_epi" else "none", ms=9,
                      label=LINEAGE_LABEL.get(k,k)))
    ax.legend(handles=handles, loc=loc, bbox_to_anchor=bbox, title=title, frameon=False,
              ncol=ncol, handletextpad=0.4)

# ---- harmonization maps (reused from analysis/interaction_map) ----
RCC_LINEAGE={
 "Tumor_RCC":"Malignant_epi","Proximal_tubule":"Tubular_epi","Intercalated":"Tubular_epi",
 "Endothelial":"Endothelial","Stroma_mural":"Stroma","Fibroblast":"Stroma",
 "Naive B cells":"B","Switched memory B cells":"B","Plasmablasts":"Plasma",
 "Effector memory CD8 T cells":"T_lineage","CD8_T":"T_lineage","T regulatory cells":"T_lineage",
 "CCR7+ T (naive/CM)":"T_lineage","Natural killer cells":"NK",
 "Intermediate monocytes":"Myeloid","Classical monocytes":"Myeloid","myeloid/DC":"Myeloid",
 "Myeloid dendritic cells":"Myeloid","mregDC":"Myeloid","Plasmacytoid dendritic cells":"Myeloid"}
CLN_LINEAGE={
 "PCT":"Tubular_epi","Connecting.tubule":"Tubular_epi","Thick.ascending.limb.of.Loop.of.Henle":"Tubular_epi",
 "Type.A.intercalated.cell":"Tubular_epi","Type.B.intercalated.cell":"Tubular_epi","Principal.cell":"Tubular_epi",
 "Pelvic.epithelium":"Tubular_epi","Epithelial.progenitor.cell":"Tubular_epi","Indistinct.intercalated.cell":"Tubular_epi",
 "Transitional.urothelium":"Tubular_epi","Parietal.epithelium":"Tubular_epi",
 "Ascending.vasa.recta.endothelium":"Endothelial","Peritubular.capillary.endothelium.1":"Endothelial",
 "Peritubular.capillary.endothelium.2":"Endothelial","Glomerular.endothelium":"Endothelial",
 "Descending.vasa.recta.endothelium":"Endothelial",
 "Fibroblast":"Stroma","Myofibroblast":"Stroma","Vascular.pericyte":"Stroma","Mesangial.cell":"Stroma",
 "Podocyte":"Podocyte","B-cell":"B","plasmablast":"Plasma",
 "macrophage":"Myeloid","monocyte":"Myeloid","mDC":"Myeloid","pDC":"Myeloid",
 "Treg":"T_lineage","T CD8 memory":"T_lineage","T CD8 naive":"T_lineage","T CD4 memory":"T_lineage","T CD4 naive":"T_lineage",
 "NK":"NK"}
DKD_NONIMM={"PT":"Tubular_epi","TAL":"Tubular_epi","PC":"Tubular_epi","CNT":"Tubular_epi","DCT":"Tubular_epi",
 "IC A":"Tubular_epi","IC B":"Tubular_epi","DTL_ATL":"Tubular_epi","PEC":"Tubular_epi","iPT":"Tubular_epi","iTAL":"Tubular_epi",
 "EC_Peritub":"Endothelial","EC_glom":"Endothelial","EC_DVR":"Endothelial","EC_Lymph":"Endothelial",
 "Fibroblast":"Stroma","VSMC":"Stroma","MC1":"Stroma","Podo":"Podocyte"}
DKD_IMMUNE={"B":"B","Plasma":"Plasma","Macro":"Myeloid","cDC":"Myeloid","pDC":"Myeloid",
 "CD4+":"T_lineage","CD8+":"T_lineage","NK":"NK"}
LINEAGE_ORDER=["B","Plasma","T_lineage","Myeloid","NK","Tubular_epi","Endothelial","Stroma","Podocyte","Malignant_epi"]

if __name__=="__main__":
    print("figstyle loaded. OUT =", OUT)
    print("datasets:", list(DATASET), "| lineages:", LINEAGE_ORDER)
