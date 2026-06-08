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
           "Unassigned (de novo)": "#8c564b", "Other/low-q": "#dddddd"}
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
# B3 — InSituType posterior / assignment-confidence distribution
# (REAL: R/15 reproduced the committed labels bit-for-bit (532,392/532,392) under the
#  same seed, so res$prob is recovered with zero desync.)
# ============================================================================
pp = os.path.join(OBJ, "wp_insitutype_posterior.csv")
if os.path.exists(pp):
    post = pd.read_csv(pp)
    post["lin"] = [lineage(v) for v in post["insitutype"].astype(str)]
    fig, (axp, axl) = plt.subplots(1, 2, figsize=(15, 5.2))
    axp.hist(post["posterior"], bins=60, color="#c44e52", alpha=0.85)
    med = post["posterior"].median(); hi = float((post["posterior"] >= 0.9).mean())
    axp.axvline(med, color="k", ls="--", lw=1, label=f"median {med:.3f}")
    axp.set_xlabel("InSituType posterior (max assignment probability)"); axp.set_ylabel("cells")
    axp.set_title(f"Assignment confidence — {hi*100:.0f}% of cells ≥ 0.90", fontsize=11); axp.legend()
    # by lineage: high-confidence FRACTION (median is ~1.0 for every lineage, so it is
    # uninformative; the ≥0.9 fraction is what differentiates the low-confidence tail).
    order_l = [k for k in LIN_COL if (post["lin"] == k).sum() > 100]
    frac = [float((post.loc[post.lin == k, "posterior"] >= 0.9).mean()) for k in order_l]
    axl.bar(range(len(order_l)), frac, color=[LIN_COL[k] for k in order_l])
    axl.set_xticks(range(len(order_l))); axl.set_xticklabels(order_l, rotation=30, ha="right")
    axl.set_ylabel("fraction of cells ≥ 0.90 posterior"); axl.set_ylim(0, 1.02)
    axl.set_title("High-confidence fraction by lineage (median ≈ 1.0 for all; the tail varies)", fontsize=11)
    fig.suptitle("InSituType per-cell assignment confidence (full cLN cohort, 532k cells)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(os.path.join(FIG, "qcB_insitutype_posterior.png"), dpi=150); plt.close(fig)
    print(f"wrote qcB_insitutype_posterior.png  (REAL; {hi*100:.0f}% cells >=0.90)")
else:
    fig, ax = plt.subplots(figsize=(9, 5)); ax.axis("off")
    ax.add_patch(plt.Rectangle((0.02, 0.02), 0.96, 0.96, fill=False, ec="#c44e52", lw=2, ls="--"))
    ax.text(0.5, 0.55, "PLACEHOLDER — InSituType posterior", ha="center", fontsize=18, color="#c44e52", weight="bold")
    ax.text(0.5, 0.38, "Run R/15 (exact-match gated) to populate; placeholder until then.", ha="center", fontsize=10)
    fig.savefig(os.path.join(FIG, "qcB_insitutype_posterior_PLACEHOLDER.png"), dpi=150); plt.close(fig)
    print("wrote qcB_insitutype_posterior_PLACEHOLDER.png  (LOGGED PLACEHOLDER)")

# ============================================================================
# B4 — immune RECALL *and* PRECISION: two-stage (R/09) vs InSituType
# (single source of truth: R/14 unified benchmark from the committed _05 object)
# ============================================================================
b = pd.read_csv(os.path.join(TAB, "cln_cosmx_immune_benchmark_unified.csv"))
b = b.sort_values("recall_insitutype", ascending=False).reset_index(drop=True)
x = np.arange(len(b)); w = 0.4
C_TS, C_IST = "#9e9e9e", "#c44e52"
fig, axes = plt.subplots(1, 2, figsize=(16, 5.4), sharex=True)
for ax, metric, ttl in [(axes[0], "recall", "Recall (vs author labels)"),
                        (axes[1], "precision", "Precision (vs author labels)")]:
    ts = b[f"{metric}_twostage"].astype(float); ist = b[f"{metric}_insitutype"].astype(float)
    ax.bar(x - w/2, ts.fillna(0), w, color=C_TS, label="two-stage co-expression (R/09)")
    ax.bar(x + w/2, ist.fillna(0), w, color=C_IST, label="InSituType")
    ax.set_xticks(x); ax.set_xticklabels(b.immune_type, rotation=30, ha="right")
    ax.set_ylabel(metric); ax.set_ylim(0, 1); ax.set_title(ttl, fontsize=11); ax.legend(fontsize=8)
# annotate two-stage over-calling (n predicted) on the precision panel, where huge
for xi, n in enumerate(b["n_pred_twostage"]):
    if n > 5000:
        axes[1].annotate(f"n={n//1000}k", (xi - w/2, 0.02), ha="center", fontsize=6.5, color="#555", rotation=90, va="bottom")
fig.suptitle("cLN immune typing — recall AND precision: two-stage (R/09) vs InSituType", fontsize=13)
fig.text(0.5, 0.005, "InSituType recovers B & DC that two-stage missed, with far higher precision; "
         "two-stage's monocyte/NK/CD8 recall is gross over-calling (n predicted ≫ truth → precision ≈ 0). "
         "CD4 T = 0 for both (panel limit).", ha="center", va="bottom", fontsize=8.5, color="#444")
fig.tight_layout(rect=(0, 0.045, 1, 0.95)); fig.savefig(os.path.join(FIG, "qcB_insitutype_recall_precision.png"), dpi=150); plt.close(fig)
print("wrote qcB_insitutype_recall_precision.png")

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
