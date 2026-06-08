#!/usr/bin/env python
"""
whitepaper_groupB_tools.py — Group B (Section 4, tools) whitepaper figures.

Inputs (committed/derived artifacts only):
  outputs/objects/wp_umap_<label>.csv            UMAP coords + label (R/13, from _05 objects)
  outputs/objects/wp_harmony_cln.csv             cLN PCA-vs-Harmony UMAP + slide + label (R/13)
  outputs/tables/cln_cosmx_insitutype_benchmark.csv   recall two-stage vs InSituType (R/10)
  outputs/objects/kidney_RCC_protein.h5ad        for the spatial-graph crop

Outputs -> outputs/figures/whitepaper/:
  qcB_umap_lineage.png        B1: per-dataset UMAP colored by major lineage
  qcB_harmony_cln.png         B2: cLN Harmony before/after (slide mixing + biology preserved)
  qcB_insitutype_posterior_PLACEHOLDER.png  B3: placeholder (posterior not persisted)
  qcB_insitutype_recall.png   B4: immune recall two-stage vs InSituType
  qcB_spatial_graph_crop.png  B5: squidpy Delaunay neighbor graph on a tissue crop

  conda run -n spatial python py/whitepaper_groupB_tools.py
"""
import os, re, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, squidpy as sq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.environ.get("XENIUM_ROOT", os.getcwd())
OBJ  = os.path.join(ROOT, "outputs/objects"); TAB = os.path.join(ROOT, "outputs/tables")
FIG  = os.path.join(ROOT, "outputs/figures/whitepaper"); os.makedirs(FIG, exist_ok=True)

DS = [("kidney_RCC_protein", "RCC (Xenium)"), ("kidney_preview_PRCC", "PRCC preview (Xenium)"),
      ("cln_cosmx", "cLN (CosMx)")]

# ---- coarse lineage mapper (handles all three label vocabularies) ----------
LIN_COL = {"Epithelial": "#8c9fca", "Endothelial": "#e377c2", "Stroma/fibroblast": "#bcbd22",
           "Tumor": "#9467bd", "Myeloid": "#2ca02c", "T cell": "#d62728", "B cell": "#1f77b4",
           "Plasma": "#ff7f0e", "NK": "#17becf", "Proliferating": "#7f7f7f",
           "Unassigned (de novo)": "#cfcfcf", "Other/low-q": "#e6e6e6"}
def lineage(x):
    s = str(x); l = s.lower()
    if len(s) == 1 and s.isalpha(): return "Unassigned (de novo)"   # cLN de-novo a-e
    if "tumor" in l or "tumour" in l: return "Tumor"
    if re.search(r"tubule|epitheli|podocyte|intercalated|principal|urotheli|henle|nephron|"
                 r"pelvic|connecting|progenitor|proximal", l): return "Epithelial"
    if re.search(r"endotheli|vasa.recta|capillary|glomerul", l): return "Endothelial"
    if re.search(r"fibroblast|myofibroblast|stroma|mural|mesangial", l): return "Stroma/fibroblast"
    if re.search(r"plasma", l): return "Plasma"
    if re.search(r"\bb[- ]?cell|switched memory b|naive b|ms4a1", l): return "B cell"
    if re.search(r"natural killer|\bnk\b", l): return "NK"
    if re.search(r"cd8|cd4|treg|regulatory t|t reg|ccr7\+ t|\bt cell|t cd|effector memory", l): return "T cell"
    if re.search(r"monocyte|macrophage|myeloid|\bdc\b|mdc|mregdc|dendritic|mast|neutrophil|pdc|plasmacytoid", l): return "Myeloid"
    if "proliferating" in l: return "Proliferating"
    if re.search(r"lowq|doublet|infiltrate/doublet|mtrnr2l", l): return "Other/low-q"
    return "Other/low-q"

def scatter_by_lineage(ax, x, y, labels, s=1.5, title=""):
    lin = np.array([lineage(v) for v in labels])
    for name, col in LIN_COL.items():
        m = lin == name
        if m.any():
            ax.scatter(x[m], y[m], s=s, c=col, linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=11)
    return lin

# ============================================================================
# B1 — per-dataset UMAP colored by major lineage
# ============================================================================
fig, axes = plt.subplots(1, 3, figsize=(19, 6.4))
present = set()
for ax, (lab, name) in zip(axes, DS):
    d = pd.read_csv(os.path.join(OBJ, f"wp_umap_{lab}.csv"))
    lin = scatter_by_lineage(ax, d.umap1.values, d.umap2.values, d.label.values, s=1.2,
                             title=f"{name}\n(n={len(d):,}, {d.label.nunique()} types)")
    present.update(np.unique(lin))
ax.set_xlabel("UMAP1")
handles = [plt.Line2D([], [], marker="o", ls="", mfc=LIN_COL[k], mec="none", ms=8, label=k)
           for k in LIN_COL if k in present]
fig.legend(handles=handles, loc="lower center", ncol=len(handles), fontsize=9,
           bbox_to_anchor=(0.5, -0.02))
fig.suptitle("Clustering UMAP per dataset, colored by major lineage "
             "(Seurat Leiden; unsupervised structure recovers the tissue compartments)", fontsize=13)
fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig(os.path.join(FIG, "qcB_umap_lineage.png"), dpi=150); plt.close(fig)
print("wrote qcB_umap_lineage.png")

# ============================================================================
# B2 — cLN Harmony before/after (slide mixing + biology preserved)
# ============================================================================
h = pd.read_csv(os.path.join(OBJ, "wp_harmony_cln.csv"))
slides = sorted(h.slide.unique()); scol = {s: plt.cm.tab20(i % 20) for i, s in enumerate(slides)}
fig, ax = plt.subplots(2, 2, figsize=(13, 12.5))
# top row: colored by SLIDE
for col, (xc, yc, ttl) in zip([0, 1], [("pca_umap1", "pca_umap2", "BEFORE integration (PCA)"),
                                        ("harm_umap1", "harm_umap2", "AFTER Harmony")]):
    a = ax[0, col]
    for s in slides:
        m = h.slide == s
        a.scatter(h[xc][m], h[yc][m], s=1.0, color=scol[s], linewidths=0, rasterized=True)
    a.set_xticks([]); a.set_yticks([]); a.set_title(ttl + " — by slide", fontsize=11)
# bottom row: colored by LINEAGE
for col, (xc, yc) in zip([0, 1], [("pca_umap1", "pca_umap2"), ("harm_umap1", "harm_umap2")]):
    scatter_by_lineage(ax[1, col], h[xc].values, h[yc].values, h.label.values, s=1.0,
                       title=("BEFORE — by lineage" if col == 0 else "AFTER — by lineage"))
sl_handles = [plt.Line2D([], [], marker="o", ls="", mfc=scol[s], mec="none", ms=6, label=s) for s in slides]
ax[0, 1].legend(handles=sl_handles, fontsize=6, ncol=1, loc="center left", bbox_to_anchor=(1.0, 0.5), title="slide")
lin_present = set(np.unique([lineage(v) for v in h.label.values]))
lin_handles = [plt.Line2D([], [], marker="o", ls="", mfc=LIN_COL[k], mec="none", ms=6, label=k)
               for k in LIN_COL if k in lin_present]
ax[1, 1].legend(handles=lin_handles, fontsize=6, ncol=1, loc="center left", bbox_to_anchor=(1.0, 0.5), title="lineage")
fig.suptitle("cLN Harmony integration (14 slides, n=60,000 shown): slides mix after integration, "
             "lineage structure preserved", fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcB_harmony_cln.png"), dpi=150); plt.close(fig)
print("wrote qcB_harmony_cln.png")

# ============================================================================
# B3 — InSituType posterior/confidence: PLACEHOLDER (not persisted)
# ============================================================================
fig, ax = plt.subplots(figsize=(9, 5)); ax.axis("off")
ax.add_patch(plt.Rectangle((0.02, 0.02), 0.96, 0.96, fill=False, ec="#c44e52", lw=2, ls="--"))
ax.text(0.5, 0.62, "PLACEHOLDER", ha="center", fontsize=22, color="#c44e52", weight="bold")
ax.text(0.5, 0.45, "InSituType per-cell posterior / assignment-confidence distribution",
        ha="center", fontsize=12)
ax.text(0.5, 0.30, "Not persisted: R/10 saved the InSituType labels but not res$prob.\n"
        "Regenerating requires re-running InSituType on 532k cells, which would re-fit the\n"
        "de-novo clusters and desync the committed labels/benchmark. Logged as a known gap.",
        ha="center", fontsize=9, color="#444")
fig.savefig(os.path.join(FIG, "qcB_insitutype_posterior_PLACEHOLDER.png"), dpi=150); plt.close(fig)
print("wrote qcB_insitutype_posterior_PLACEHOLDER.png  (LOGGED PLACEHOLDER)")

# ============================================================================
# B4 — immune recall: two-stage (before) vs InSituType (after)
# ============================================================================
b = pd.read_csv(os.path.join(TAB, "cln_cosmx_insitutype_benchmark.csv"))
b = b.sort_values("recall_IST", ascending=False)
x = np.arange(len(b)); w = 0.4
fig, ax = plt.subplots(figsize=(12, 5.2))
ax.bar(x - w/2, b.recall_09.fillna(0), w, color="#bbbbbb", label="two-stage (before)")
ax.bar(x + w/2, b.recall_IST.fillna(0), w, color="#c44e52", label="InSituType (after)")
ax.set_xticks(x); ax.set_xticklabels(b.immune_type, rotation=30, ha="right")
ax.set_ylabel("recall vs author labels"); ax.set_ylim(0, 1)
ax.set_title("cLN immune recall — two-stage co-expression vs InSituType "
             "(contamination-aware recovers B / DC / plasma; CD4 T stays 0 = panel limit)", fontsize=10)
ax.legend()
for xi, (r0, r1) in enumerate(zip(b.recall_09.fillna(0), b.recall_IST.fillna(0))):
    ax.text(xi + w/2, r1 + 0.02, f"{r1:.2f}", ha="center", fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcB_insitutype_recall.png"), dpi=150); plt.close(fig)
print("wrote qcB_insitutype_recall.png")

# ============================================================================
# B5 — squidpy Delaunay spatial neighbor graph on a tissue crop
# ============================================================================
a = ad.read_h5ad(os.path.join(OBJ, "kidney_RCC_protein.h5ad"))
xy = np.asarray(a.obsm["spatial"], float)
# pick an immune-rich crop: densest 300um window in B + Treg
key = "phase_b_label"; lab = a.obs[key].astype(str).values
foc = np.isin(lab, ["Naive B cells", "T regulatory cells"])
cx, cy = np.median(xy[foc, 0]), np.median(xy[foc, 1]); W = 200.0
m = (np.abs(xy[:, 0] - cx) < W) & (np.abs(xy[:, 1] - cy) < W)
crop = a[m].copy(); crop.obsm["spatial"] = np.asarray(crop.obsm["spatial"], float)
sq.gr.spatial_neighbors(crop, coord_type="generic", delaunay=True)
cxy = crop.obsm["spatial"]; clab = crop.obs[key].astype(str).values
G = crop.obsp["spatial_connectivities"].tocoo()
fig, ax = plt.subplots(figsize=(9, 9))
for i, j in zip(G.row, G.col):
    if i < j:
        ax.plot(cxy[[i, j], 0], cxy[[i, j], 1], color="#cfcfcf", lw=0.3, zorder=1)
clin = np.array([lineage(v) for v in clab])
for name, c in LIN_COL.items():
    mm = clin == name
    if mm.any():
        ax.scatter(cxy[mm, 0], cxy[mm, 1], s=22, c=c, linewidths=0.2, edgecolor="k", zorder=2, label=name)
ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
ax.set_title(f"squidpy Delaunay neighbor graph — RCC tissue crop ({crop.n_obs} cells, 400×400 µm)\n"
             "nodes = cells (by lineage), edges = spatial adjacency", fontsize=10)
ax.legend(fontsize=7, loc="upper left", markerscale=1.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcB_spatial_graph_crop.png"), dpi=150); plt.close(fig)
print("wrote qcB_spatial_graph_crop.png")
print("== Group B done ==")
