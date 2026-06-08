#!/usr/bin/env python
"""
whitepaper_recompute.py — two pre-assembly reconciliations:

(1) cLN per-slide nhood z on the PER-SECTION graph: Delaunay then prune cross-core edges
    (>50 um) so the graph no longer spans the mm gaps between tissue cores. These pruned
    values become PRIMARY in cln_phaseB_per_slide.csv (nhood_z_plasma_myeloid); the
    unpruned value is kept as nhood_z_unpruned for the sensitivity check.

(2) RCC B-aggregate composition at eps=50 um (37 aggregates — matches the committed
    Group-D overview), updating rcc_phaseB2_aggregate_composition.csv so the figure,
    table, and narrative all agree on the aggregate count.

  conda run -n spatial python py/whitepaper_recompute.py
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, squidpy as sq
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon

ROOT = os.environ.get("XENIUM_ROOT", os.getcwd())
OBJ = os.path.join(ROOT, "outputs/objects"); TAB = os.path.join(ROOT, "outputs/tables")
THR = 0.05  # mm cross-core prune / um footprint

# ============================================================================
# (1) cLN per-section (pruned) nhood z
# ============================================================================
a = ad.read_h5ad(os.path.join(OBJ, "cln_cosmx.h5ad"))
xy_all = np.asarray(a.obsm["spatial"], float); lab_all = a.obs["phase_b_label"].astype(str).values
slide_all = a.obs["sample"].astype(str).values
PLASMA = ["plasmablast"]; MYELOID = ["macrophage", "myeloid/DC", "monocyte"]; MIN = 10

def nhood_z(slide, prune):
    m = slide_all == slide; sub = a[m].copy(); sub.obsm["spatial"] = xy_all[m]
    nl = np.full(m.sum(), "Other", object)
    nl[np.isin(lab_all[m], PLASMA)] = "Plasma"; nl[np.isin(lab_all[m], MYELOID)] = "Myeloid"
    if (nl == "Plasma").sum() < MIN or (nl == "Myeloid").sum() < MIN:
        return np.nan
    sq.gr.spatial_neighbors(sub, coord_type="generic", delaunay=True)
    if prune:
        keep = sub.obsp["spatial_distances"].copy(); keep.data = (keep.data <= THR).astype(float)
        sub.obsp["spatial_connectivities"] = sub.obsp["spatial_connectivities"].multiply(keep).tocsr()
    present = [c for c in ["Plasma", "Myeloid", "Other"] if (nl == c).sum() > 0]
    sub.obs["nl2"] = pd.Categorical(nl, categories=present)
    sq.gr.nhood_enrichment(sub, cluster_key="nl2", seed=0, show_progress_bar=False)
    Z = pd.DataFrame(sub.uns["nl2_nhood_enrichment"]["zscore"], index=present, columns=present)
    return float(Z.loc["Plasma", "Myeloid"]) if "Plasma" in present and "Myeloid" in present else np.nan

ps = pd.read_csv(os.path.join(TAB, "cln_phaseB_per_slide.csv"))
ps["nhood_z_unpruned"] = ps["nhood_z_plasma_myeloid"]
pruned = []
for s in ps.slide:
    zp = nhood_z(s, prune=True)
    pruned.append(round(zp, 2) if zp == zp else np.nan)
    print(f"  {s}: unpruned {ps.loc[ps.slide==s,'nhood_z_unpruned'].iloc[0]} -> per-section {pruned[-1]}")
ps["nhood_z_plasma_myeloid"] = pruned
ps.to_csv(os.path.join(TAB, "cln_phaseB_per_slide.csv"), index=False)
dmax = float(np.nanmax(np.abs(ps["nhood_z_plasma_myeloid"] - ps["nhood_z_unpruned"])))
print(f"updated cln_phaseB_per_slide.csv (per-section primary). max |Δz| = {dmax:.2f}")

# ============================================================================
# (2) RCC B-aggregate composition at eps=50 (37 aggregates)
# ============================================================================
r = ad.read_h5ad(os.path.join(OBJ, "kidney_RCC_protein.h5ad"))
rxy = np.asarray(r.obsm["spatial"], float); rlab = r.obs["phase_b_label"].astype(str).values
B = "Naive B cells"; TREG = "T regulatory cells"; PLAS = "Plasmablasts"; MREG = "mregDC"
CD8 = ["Effector memory CD8 T cells", "CD8_T"]
isB = rlab == B; Bidx = np.where(isB)[0]
cl = DBSCAN(eps=50.0, min_samples=20).fit(rxy[Bidx]).labels_
aggs = [c for c in np.unique(cl) if c != -1]
tree = cKDTree(rxy); R = 50.0
bg = {"Treg": float((rlab == TREG).mean()), "eff-CD8": float(np.isin(rlab, CD8).mean()),
      "mregDC": float((rlab == MREG).mean()), "Plasma": float((rlab == PLAS).mean())}
inside = {k: [] for k in bg}
for c in aggs:
    mem = Bidx[cl == c]
    reg = np.unique(np.concatenate([np.asarray(t, int) for t in tree.query_ball_point(rxy[mem], r=R)]))
    inside["Treg"].append(float((rlab[reg] == TREG).mean()))
    inside["eff-CD8"].append(float(np.isin(rlab[reg], CD8).mean()))
    inside["mregDC"].append(float((rlab[reg] == MREG).mean()))
    inside["Plasma"].append(float((rlab[reg] == PLAS).mean()))
rows = []
for k in ["Treg", "eff-CD8", "Plasma", "mregDC"]:
    ins = np.array(inside[k]); b = bg[k]; im = float(ins.mean())
    l2 = float(np.log2((im + 1e-6) / (b + 1e-6)))
    try: _, p = wilcoxon(ins - b, alternative="two-sided", zero_method="zsplit")
    except Exception: p = np.nan
    rows.append(dict(cell_type=k, background_frac=round(b, 4), inside_mean_frac=round(im, 4),
                     log2_enrichment=round(l2, 3), direction=("ENRICHED" if l2 > 0 else "EXCLUDED"),
                     wilcoxon_p=round(p, 4) if p == p else np.nan,
                     n_agg_above_bg=int((ins > b).sum()), n_agg=len(aggs)))
comp = pd.DataFrame(rows)
comp.to_csv(os.path.join(TAB, "rcc_phaseB2_aggregate_composition.csv"), index=False)
print(f"\nRCC B-aggregates at eps=50: {len(aggs)} aggregates")
print(comp.to_string(index=False))

# ============================================================================
# (3) rebuild cLN rows of the RCC-vs-cLN contrast table so n_plasma_aggregates is a
#     consistent TOTAL across rows (was: overall=sum, per-class=mean) and nhood z uses
#     the per-section (pruned) values from ps. RCC rows are kept unchanged.
# ============================================================================
ctab = os.path.join(TAB, "plasma_myeloid_rcc_vs_cln.csv")
comp_t = pd.read_csv(ctab); rcc_rows = comp_t[comp_t.context.str.startswith("RCC")].copy()
cols = ["context", "dataset", "n_plasma", "n_myeloid", "myeloid_frac", "plasma_myeloid_nhood_z",
        "n_plasma_aggregates", "plasma_agg_rate", "myeloid_log2enrich_in_plasma_aggs"]
pb = ps[ps.n_plasma >= 10]
def crow(dataset, sub, overall):
    return {"context": "cLN (CosMx kidney)", "dataset": dataset,
            "n_plasma": int(sub.n_plasma.sum()) if overall else np.nan,
            "n_myeloid": int(sub.n_myeloid.sum()) if overall else np.nan,
            "myeloid_frac": round(sub.myeloid_frac.mean(), 4),
            "plasma_myeloid_nhood_z": round(sub.nhood_z_plasma_myeloid.median() if overall
                                            else sub.nhood_z_plasma_myeloid.mean(), 3),
            "n_plasma_aggregates": int(sub.n_plasma_aggregates.sum()),               # TOTAL everywhere
            "plasma_agg_rate": round(sub.plasma_agg_rate.median() if overall else sub.plasma_agg_rate.mean(), 3),
            "myeloid_log2enrich_in_plasma_aggs": round(sub.myeloid_log2enrich_in_plasma_aggs.mean(), 3)}
cln_rows = [crow("cLN_overall (plasma slides, totals)", pb, True)]
for c in ["III", "IV", "IV+V"]:
    cln_rows.append(crow(f"cLN_{c} (class total/mean)", ps[ps.condition == c], False))
pd.concat([rcc_rows[cols], pd.DataFrame(cln_rows)[cols]], ignore_index=True).to_csv(ctab, index=False)
print("\nrebuilt plasma_myeloid_rcc_vs_cln.csv (n_plasma_aggregates consistent = totals; per-section z)")
print("\n== recompute done ==")
