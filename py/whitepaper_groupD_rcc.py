#!/usr/bin/env python
"""
whitepaper_groupD_rcc.py — Group D (Section 6.1, RCC spatial) whitepaper figures.

From the committed RCC h5ad (phase_b_label + coords + counts) and the committed Phase-B
composition tables. Recomputes the spatial graph / nhood / B-aggregate DBSCAN exactly as
Phase B (same params: Delaunay; eps=50um, min=20) so figures are self-contained.

Outputs -> outputs/figures/whitepaper/:
  qcD_nhood_heatmap.png        D1: neighborhood-enrichment z heatmap
  qcD_aggregates_overview.png  D2: delineated B-aggregates colored by cell type
  qcD_aggregate_markers.png    D3: representative aggregate crop + MS4A1/FOXP3/LAMP3/CD8A
  qcD_mregdc_ccr7_foci.png     D4: representative mregDC-CCR7+T focus crop
  qcD_big_vs_preview_comp.png  D5: BIG vs PRCC-preview aggregate composition

  conda run -n spatial python py/whitepaper_groupD_rcc.py
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, squidpy as sq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree

ROOT = os.environ.get("XENIUM_ROOT", os.getcwd())
OBJ = os.path.join(ROOT, "outputs/objects"); TAB = os.path.join(ROOT, "outputs/tables")
FIG = os.path.join(ROOT, "outputs/figures/whitepaper"); os.makedirs(FIG, exist_ok=True)
KEY = "phase_b_label"
B, TREG, PLASMA, MREG, CCR7 = ("Naive B cells", "T regulatory cells", "Plasmablasts",
                               "mregDC", "CCR7+ T (naive/CM)")
CD8 = ["Effector memory CD8 T cells", "CD8_T"]
FOCUS = [B, TREG] + CD8 + [PLASMA, MREG, CCR7]
COL = {B: "#1f77b4", TREG: "#d62728", "Effector memory CD8 T cells": "#2ca02c",
       "CD8_T": "#2ca02c", PLASMA: "#ff7f0e", MREG: "#9467bd", CCR7: "#17becf"}
LEG = [(B, "#1f77b4", "Naive B"), (TREG, "#d62728", "Treg"),
       ("eff-CD8", "#2ca02c", "effector CD8"), (PLASMA, "#ff7f0e", "Plasmablast"),
       (MREG, "#9467bd", "mregDC"), (CCR7, "#17becf", "CCR7+ T")]

a = ad.read_h5ad(os.path.join(OBJ, "kidney_RCC_protein.h5ad"))
xy = np.asarray(a.obsm["spatial"], float); lab = a.obs[KEY].astype(str).values
a.obs[KEY] = a.obs[KEY].astype("category")
print(f"loaded {a.n_obs} cells")
sq.gr.spatial_neighbors(a, coord_type="generic", delaunay=True)

# ============================================================================
# D1 — neighborhood-enrichment z heatmap
# ============================================================================
sq.gr.nhood_enrichment(a, cluster_key=KEY, seed=0, show_progress_bar=False)
cats = list(a.obs[KEY].cat.categories)
Z = pd.DataFrame(a.uns[f"{KEY}_nhood_enrichment"]["zscore"], index=cats, columns=cats)
# order: immune/aggregate types first for readability
imm_first = [c for c in [B, TREG, "Effector memory CD8 T cells", "CD8_T", PLASMA, MREG, CCR7] if c in cats]
rest = [c for c in cats if c not in imm_first]
order = imm_first + rest
Zo = Z.loc[order, order]
vmax = np.nanpercentile(np.abs(Zo.values), 98)
fig, ax = plt.subplots(figsize=(11, 9.5))
im = ax.imshow(Zo.values, cmap="coolwarm", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=90, fontsize=7)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=7)
ax.axhline(len(imm_first) - 0.5, color="k", lw=0.6); ax.axvline(len(imm_first) - 0.5, color="k", lw=0.6)
# box the load-bearing B x Treg (+) and B x eff-CD8 (-)
bi, ti = order.index(B), order.index(TREG); ci = order.index("Effector memory CD8 T cells")
for (r, c2, col) in [(bi, ti, "#000"), (bi, ci, "#000")]:
    ax.add_patch(plt.Rectangle((c2-0.5, r-0.5), 1, 1, fill=False, ec=col, lw=2))
fig.colorbar(im, ax=ax, fraction=0.046, label="neighborhood-enrichment z")
ax.set_title("RCC neighborhood-enrichment (z): B×Treg enriched (+), B×effector-CD8 excluded (−)\n"
             "immune/aggregate types upper-left block; boxes = the load-bearing pairs", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcD_nhood_heatmap.png"), dpi=150); plt.close(fig)
print(f"wrote qcD_nhood_heatmap.png  (B×Treg z={Z.loc[B,TREG]:.0f}, B×eff-CD8 z={Z.loc[B,'Effector memory CD8 T cells']:.0f})")

# ============================================================================
# B-aggregate DBSCAN (same params as Phase B pass 2)
# ============================================================================
isB = lab == B; Bidx = np.where(isB)[0]
cl = DBSCAN(eps=50.0, min_samples=20).fit(xy[Bidx]).labels_
agg_ids = [c for c in np.unique(cl) if c != -1]
tree = cKDTree(xy)
aggs = []
for c in agg_ids:
    members = Bidx[cl == c]; cen = xy[members].mean(0)
    aggs.append(dict(c=c, n=len(members), x=cen[0], y=cen[1], members=members))
aggs = sorted(aggs, key=lambda d: d["n"], reverse=True)
print(f"  {len(aggs)} B-aggregates (largest {aggs[0]['n']} B cells)")

# ============================================================================
# D2 — delineated aggregates overview, colored by cell type
# ============================================================================
fig, ax = plt.subplots(figsize=(13, 7))
foc = np.isin(lab, FOCUS)
ax.scatter(xy[~foc, 0], xy[~foc, 1], s=0.5, c="#e8e8e8", linewidths=0, rasterized=True)
for L, c, _n in LEG:
    sel = np.isin(lab, CD8) if L == "eff-CD8" else (lab == L)
    ax.scatter(xy[sel, 0], xy[sel, 1], s=1.5, c=c, linewidths=0, rasterized=True)
for d in aggs:
    ax.add_patch(plt.Circle((d["x"], d["y"]), 90, fill=False, ec="#222", lw=0.8, alpha=0.7))
ax.set_aspect("equal"); ax.axis("off")
handles = [plt.Line2D([], [], marker="o", ls="", mfc=c, mec="none", ms=7, label=n) for _, c, n in LEG]
handles.append(plt.Line2D([], [], marker="o", ls="", mfc="none", mec="#222", ms=10, label="delineated B-aggregate"))
ax.legend(handles=handles, fontsize=8, loc="upper right", markerscale=1.2)
ax.set_title(f"RCC: {len(aggs)} delineated B-cell aggregates (circles), cells colored by type", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcD_aggregates_overview.png"), dpi=160); plt.close(fig)
print("wrote qcD_aggregates_overview.png")

# ============================================================================
# D3 — representative aggregate crop + marker overlays (the derivation)
# ============================================================================
def gene_vec(g):
    j = a.var_names.get_loc(g)
    return np.asarray(a.X[:, j].todense()).ravel() if hasattr(a.X, "todense") else np.asarray(a.X[:, j]).ravel()
# pick a WELL-MIXED representative aggregate (not the largest, which is B-saturated):
# among mid-sized aggregates, maximize presence of the rarest focus types so all four
# markers (B/Treg/mregDC/CD8) are visible in the crop.
W = 180.0
def crop_counts(d):
    mm = (np.abs(xy[:, 0] - d["x"]) < W) & (np.abs(xy[:, 1] - d["y"]) < W)
    sl_ = lab[mm]
    return dict(B=int((sl_ == B).sum()), Treg=int((sl_ == TREG).sum()),
                CD8=int(np.isin(sl_, CD8).sum()), mreg=int((sl_ == MREG).sum()))
cand = [d for d in aggs if 150 <= d["n"] <= 900]
scored = [(min(c["Treg"], c["CD8"]) + 3 * c["mreg"], d) for d in cand for c in [crop_counts(d)]]
rep = max(scored, key=lambda t: t[0])[1] if scored else aggs[0]
cx, cy = rep["x"], rep["y"]
m = (np.abs(xy[:, 0] - cx) < W) & (np.abs(xy[:, 1] - cy) < W)
sx, sl = xy[m], lab[m]
print(f"  D3 representative aggregate: {rep['n']} B cells; crop {crop_counts(rep)}")
markers = [("MS4A1", "B cells"), ("FOXP3", "Treg"), ("LAMP3", "mregDC"), ("CD8A", "CD8 T")]
fig, ax = plt.subplots(1, 5, figsize=(24, 5))
# panel 0: cell type
ax[0].scatter(sx[~np.isin(sl, FOCUS), 0], sx[~np.isin(sl, FOCUS), 1], s=10, c="#dddddd", linewidths=0)
for L, c, n in LEG:
    sel = np.isin(sl, CD8) if L == "eff-CD8" else (sl == L)
    if sel.any(): ax[0].scatter(sx[sel, 0], sx[sel, 1], s=20, c=c, linewidths=0, label=n)
ax[0].set_title("cell type", fontsize=11); ax[0].legend(fontsize=7, loc="upper left")
for k, (g, who) in enumerate(markers, start=1):
    gv = gene_vec(g)[m]
    pos = gv > 0
    ax[k].scatter(sx[~pos, 0], sx[~pos, 1], s=8, c="#e3e3e3", linewidths=0)   # gray = marker-negative
    sc = ax[k].scatter(sx[pos, 0], sx[pos, 1], s=34, c=gv[pos], cmap="viridis", linewidths=0.2,
                       edgecolor="k", vmin=1, vmax=max(2, np.percentile(gv[pos], 95)) if pos.any() else 2)
    ax[k].set_title(f"{g}  ({who}) — {pos.sum()} pos cells", fontsize=11)
    fig.colorbar(sc, ax=ax[k], fraction=0.046, label="counts")
for x_ in ax:
    x_.set_aspect("equal"); x_.set_xticks([]); x_.set_yticks([])
fig.suptitle(f"Representative RCC B-aggregate ({rep['n']} B cells) — cell-type calls beside the "
             "transcript markers that back them (MS4A1=B, FOXP3=Treg, LAMP3=mregDC, CD8A=CD8)", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcD_aggregate_markers.png"), dpi=150); plt.close(fig)
print("wrote qcD_aggregate_markers.png")

# ============================================================================
# D4 — representative mregDC-CCR7+T focus crop
# ============================================================================
isc24 = np.isin(lab, [MREG, CCR7]); c24idx = np.where(isc24)[0]
cl2 = DBSCAN(eps=50.0, min_samples=10).fit(xy[c24idx]).labels_
foci = [(c, np.sum(cl2 == c)) for c in np.unique(cl2) if c != -1]
foci = sorted(foci, key=lambda t: t[1], reverse=True)
fc = c24idx[cl2 == foci[0][0]]; fx, fy = xy[fc].mean(0)
W2 = 120.0
m2 = (np.abs(xy[:, 0] - fx) < W2) & (np.abs(xy[:, 1] - fy) < W2)
sx2, sl2 = xy[m2], lab[m2]
fig, ax = plt.subplots(figsize=(7.5, 7))
ax.scatter(sx2[~np.isin(sl2, FOCUS), 0], sx2[~np.isin(sl2, FOCUS), 1], s=18, c="#dddddd", linewidths=0)
for L, c, n in LEG:
    sel = np.isin(sl2, CD8) if L == "eff-CD8" else (sl2 == L)
    if sel.any(): ax.scatter(sx2[sel, 0], sx2[sel, 1], s=40, c=c, linewidths=0.3, edgecolor="k", label=n)
ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.legend(fontsize=8, loc="upper left")
ax.set_title(f"Representative mregDC–CCR7+T focus (1 of {len(foci)} across the section)", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcD_mregdc_ccr7_foci.png"), dpi=150); plt.close(fig)
print("wrote qcD_mregdc_ccr7_foci.png")

# ============================================================================
# D5 — BIG vs PRCC-preview aggregate composition (log2 enrichment)
# ============================================================================
big = pd.read_csv(os.path.join(TAB, "rcc_phaseB2_aggregate_composition.csv")).set_index("cell_type")["log2_enrichment"]
rep_df = pd.read_csv(os.path.join(TAB, "preview_phaseB_replication_summary.csv"))
def prev_l2(name):
    s = rep_df[rep_df.effect == name]["preview"].iloc[0]
    return float(str(s).split("log2")[1].split("(")[0])
types = [("Treg", "Treg in B-agg", "Treg"), ("eff-CD8", "eff-CD8 in B-agg", "effector-CD8"),
         ("mregDC", "mregDC in B-agg", "mregDC"), ("Plasma", "Plasma in B-agg", "Plasma")]
big_v = [float(big[t[0]]) for t in types]
prev_v = [prev_l2(t[1]) for t in types]
x = np.arange(len(types)); w = 0.38
fig, ax = plt.subplots(figsize=(9, 5.4))
ax.bar(x - w/2, big_v, w, color="#4c72b0", label="RCC (BIG, Xenium)")
ax.bar(x + w/2, prev_v, w, color="#55a868", label="PRCC preview (Xenium)")
ax.axhline(0, color="#999", lw=0.7)
ax.set_xticks(x); ax.set_xticklabels([t[2] for t in types])
ax.set_ylabel("log2 enrichment inside B-aggregates")
ax.set_title("B-aggregate composition replicates: Treg enriched, effector-CD8 excluded\n"
             "(RCC vs independent PRCC preview; mregDC/Plasma underpowered in preview)", fontsize=10)
for xi, (bv, pv) in enumerate(zip(big_v, prev_v)):
    ax.text(xi - w/2, bv + (0.05 if bv >= 0 else -0.12), f"{bv:+.2f}", ha="center", fontsize=8)
    ax.text(xi + w/2, pv + (0.05 if pv >= 0 else -0.12), f"{pv:+.2f}", ha="center", fontsize=8)
ax.legend()
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcD_big_vs_preview_comp.png"), dpi=150); plt.close(fig)
print("wrote qcD_big_vs_preview_comp.png")
print("== Group D done ==")
