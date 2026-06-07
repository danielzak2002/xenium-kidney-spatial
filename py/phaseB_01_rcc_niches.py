#!/usr/bin/env python
"""
phaseB_01_rcc_niches.py — first squidpy pass on the RCC Xenium export.
Bounded goal: does TLS-like B/T/mregDC co-organization EXIST and where do BCMA+
plasma sit. NOT formal TLS segmentation / niche-domain clustering (that's pass 2).

Run in the `spatial` conda env:
  conda run -n spatial python py/phaseB_01_rcc_niches.py

RCC is a SINGLE section -> one spatial graph on the whole object. Reads effect SIZE
(z magnitude / enrichment), not significance (everything is "significant" at 465k).
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, squidpy as sq, scanpy as sc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu

ROOT = os.environ.get("XENIUM_ROOT", os.getcwd())
H5  = os.path.join(ROOT, "outputs/objects/kidney_RCC_protein.h5ad")
FIG = os.path.join(ROOT, "outputs/figures"); TAB = os.path.join(ROOT, "outputs/tables")
KEY = "phase_b_label"
B   = "Naive B cells"; MREG = "mregDC"; PLASMA = "Plasmablasts"
TCELLS = ["Effector memory CD8 T cells", "T regulatory cells", "CD8_T", "CCR7+ T (naive/CM)"]

a = ad.read_h5ad(H5)
assert a.uns["spatial_units"] == "um", a.uns["spatial_units"]
a.obs[KEY] = a.obs[KEY].astype("category")
a.obsm["spatial"] = np.asarray(a.obsm["spatial"], float)
print(f"loaded {a.n_obs} cells, {a.n_vars} genes; key={KEY} ({a.obs[KEY].nunique()} types)")

# ---- spatial graph (single section) ----------------------------------------
sq.gr.spatial_neighbors(a, coord_type="generic", delaunay=True)

# ---- 1. neighborhood enrichment (effect size = z) --------------------------
sq.gr.nhood_enrichment(a, cluster_key=KEY, seed=0, show_progress_bar=False)
cats = list(a.obs[KEY].cat.categories)
Z = pd.DataFrame(a.uns[f"{KEY}_nhood_enrichment"]["zscore"], index=cats, columns=cats)
fig, ax = plt.subplots(figsize=(11, 9))
sq.pl.nhood_enrichment(a, cluster_key=KEY, method="average", ax=ax, cmap="coolwarm")
fig.tight_layout(); fig.savefig(f"{FIG}/rcc_phaseB_nhood_heatmap.png", dpi=170); plt.close(fig)

def zpair(i, j):
    return float(Z.loc[i, j]) if (i in Z.index and j in Z.columns) else np.nan
tls = []
for t in TCELLS:
    tls.append(("B x " + t, zpair(B, t)))
tls += [("B x mregDC", zpair(B, MREG)), ("mregDC x CCR7+ T", zpair(MREG, "CCR7+ T (naive/CM)")),
        ("B x Plasmablasts", zpair(B, PLASMA)), ("mregDC x EM-CD8", zpair(MREG, "Effector memory CD8 T cells"))]
tls_df = pd.DataFrame(tls, columns=["pair", "nhood_zscore"]).sort_values("nhood_zscore", ascending=False)
print("\nTLS-pair neighborhood enrichment (z):"); print(tls_df.to_string(index=False))

# ---- 2. co-occurrence over distance (subsample for memory) ------------------
key_cl = [B, MREG, PLASMA] + TCELLS
sub = a[a.obs[KEY].isin(key_cl)].copy()
if sub.n_obs > 60000:
    idx = np.random.default_rng(0).choice(sub.n_obs, 60000, replace=False); sub = sub[idx].copy()
sub.obs[KEY] = sub.obs[KEY].cat.remove_unused_categories()
sq.gr.co_occurrence(sub, cluster_key=KEY, interval=np.linspace(0, 200, 25), show_progress_bar=False)
for anchor, fn in [(B, "rcc_phaseB_cooc_B.png"), (MREG, "rcc_phaseB_cooc_mregDC.png")]:
    if anchor in list(sub.obs[KEY].cat.categories):
        fig = sq.pl.co_occurrence(sub, cluster_key=KEY, clusters=anchor, figsize=(7, 5))
        plt.tight_layout(); plt.savefig(f"{FIG}/{fn}", dpi=160); plt.close("all")

# ---- 3. visualize: full + crops to candidate aggregates --------------------
xy = a.obsm["spatial"]; lab = a.obs[KEY].values
immune_focus = [B, MREG, PLASMA] + TCELLS
pal = {c: ("#cccccc" if c not in immune_focus else None) for c in cats}
fig = sq.pl.spatial_scatter(a, color=KEY, shape=None, size=2, figsize=(12, 10),
                            groups=immune_focus, return_ax=False)
plt.savefig(f"{FIG}/rcc_phaseB_spatial_immune.png", dpi=160, bbox_inches="tight"); plt.close("all")

# candidate aggregates = grid cells densest in B (TLS seeds), then crop +/-300um
isB = lab == B
if isB.sum() > 0:
    bxy = xy[isB]; gs = 150.0
    cell = np.floor(bxy / gs).astype(int)
    df = pd.DataFrame(cell, columns=["gx", "gy"]); top = df.value_counts().head(3)
    crops = [(gx * gs + gs / 2, gy * gs + gs / 2) for (gx, gy) in top.index]
    fig, axes = plt.subplots(1, len(crops), figsize=(6 * len(crops), 6))
    axes = np.atleast_1d(axes)
    for ax, (cx, cy) in zip(axes, crops):
        m = (np.abs(xy[:, 0] - cx) < 300) & (np.abs(xy[:, 1] - cy) < 300)
        sub_lab = lab[m]
        for c in immune_focus + ["_other"]:
            sel = (sub_lab == c) if c != "_other" else ~np.isin(sub_lab, immune_focus)
            ax.scatter(xy[m][sel, 0], xy[m][sel, 1], s=6,
                       c=("#dddddd" if c == "_other" else None), label=(None if c == "_other" else c))
        ax.set_title(f"crop @({cx:.0f},{cy:.0f})"); ax.set_aspect("equal"); ax.axis("off")
    axes[0].legend(markerscale=2, fontsize=7, loc="upper left", bbox_to_anchor=(0, 1))
    fig.tight_layout(); fig.savefig(f"{FIG}/rcc_phaseB_tls_crops.png", dpi=170); plt.close(fig)

# ---- 4. BCMA+ vs BCMA- plasma ----------------------------------------------
g = "TNFRSF17"; gi = a.var_names.get_loc(g)
cnt = np.asarray(a.X[:, gi].todense()).ravel() if hasattr(a.X, "todense") else a.X[:, gi]
isPL = (lab == PLASMA)
bcma_counts = cnt[isPL]
frac2 = float((bcma_counts >= 2).mean()); frac1 = float((bcma_counts >= 1).mean())
print(f"\nBCMA(TNFRSF17) in plasma (n={isPL.sum()}): det>=1={frac1:.2f}, det>=2={frac2:.2f}, "
      f"median(+cells)={np.median(bcma_counts[bcma_counts>0]) if (bcma_counts>0).any() else 0:.0f}")
bcma_pos = isPL & (cnt >= 2)
# distance to nearest B cell, BCMA+ vs BCMA- plasma
tree = cKDTree(xy[isB])
d_pos, _ = tree.query(xy[bcma_pos]); d_neg, _ = tree.query(xy[isPL & (cnt < 2)])
U, p = mannwhitneyu(d_pos, d_neg, alternative="less")
print(f"dist-to-nearest-B (um): BCMA+ median={np.median(d_pos):.0f}  BCMA- median={np.median(d_neg):.0f}  "
      f"(Mann-Whitney BCMA+<BCMA- p={p:.1e})")
fig, ax = plt.subplots(figsize=(7, 5))
for d, lb in [(d_pos, f"BCMA+ plasma (n={bcma_pos.sum()})"), (d_neg, f"BCMA- plasma (n={(isPL&(cnt<2)).sum()})")]:
    xs = np.sort(d); ax.plot(xs, np.linspace(0, 1, len(xs)), label=lb)
ax.set_xlim(0, 400); ax.set_xlabel("distance to nearest B cell (um)"); ax.set_ylabel("ECDF")
ax.set_title("BCMA+ vs BCMA- plasma: proximity to B cells"); ax.legend()
fig.tight_layout(); fig.savefig(f"{FIG}/rcc_phaseB_bcma_plasma.png", dpi=160); plt.close(fig)
# nhood enrichment with plasma split
a.obs["tls_key"] = a.obs[KEY].astype(str)
a.obs.loc[bcma_pos, "tls_key"] = "Plasma_BCMA+"
a.obs.loc[isPL & (cnt < 2), "tls_key"] = "Plasma_BCMA-"
a.obs["tls_key"] = a.obs["tls_key"].astype("category")
sq.gr.nhood_enrichment(a, cluster_key="tls_key", seed=0, show_progress_bar=False)
Zt = pd.DataFrame(a.uns["tls_key_nhood_enrichment"]["zscore"],
                  index=list(a.obs["tls_key"].cat.categories), columns=list(a.obs["tls_key"].cat.categories))

# ---- summary table ----------------------------------------------------------
rows = list(tls_df.itertuples(index=False, name=None))
rows += [("BCMA+plasma x B", float(Zt.loc["Plasma_BCMA+", B]) if "Plasma_BCMA+" in Zt.index else np.nan),
         ("BCMA-plasma x B", float(Zt.loc["Plasma_BCMA-", B]) if "Plasma_BCMA-" in Zt.index else np.nan),
         ("BCMA+ dist-to-B (um, median)", float(np.median(d_pos))),
         ("BCMA- dist-to-B (um, median)", float(np.median(d_neg))),
         ("plasma BCMA+ fraction (>=2)", frac2)]
pd.DataFrame(rows, columns=["metric", "value"]).to_csv(f"{TAB}/rcc_phaseB_summary.csv", index=False)
print("\nwrote figures + rcc_phaseB_summary.csv\n== phaseB_01 done ==")
