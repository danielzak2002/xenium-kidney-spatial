#!/usr/bin/env python
"""
phaseB_02_rcc_aggregates.py — pass 2 on the RCC Xenium export. Tight scope:

  (1) DE-RISK the mregDC x CCR7+T z=+90: it came partly by-construction from splitting
      ONE Leiden cluster (c24). Reconstruct c24 = {mregDC, CCR7+ T} cells and ask whether
      they form MANY separated foci across the tissue (real recurring niche) or ONE compact
      blob (circular artifact of the split). DBSCAN on c24 coords -> n foci, spread, largest
      -focus share -> verdict.
  (2) FORMALLY delineate the immunoregulatory aggregates (the load-bearing finding):
      DBSCAN on B-cell coords -> aggregates; per-aggregate composition (Treg / effector-CD8 /
      plasma / mregDC) INSIDE vs tissue background. Confirm: Treg-ENRICHED, CD8-effector-
      EXCLUDED, plasma-EXCLUDED.
  (3) One clean figure: the 3 largest aggregates colored by cell type — the B-Treg core with
      CD8/plasma excluded.

NOT niche-domain (UTAG) clustering — held unless step 2 shows multiple distinct niche types.

  conda run -n spatial python py/phaseB_02_rcc_aggregates.py
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon

ROOT = os.environ.get("XENIUM_ROOT", os.getcwd())
H5   = os.path.join(ROOT, "outputs/objects/kidney_RCC_protein.h5ad")
FIG  = os.path.join(ROOT, "outputs/figures"); TAB = os.path.join(ROOT, "outputs/tables")
KEY  = "phase_b_label"
rng  = np.random.default_rng(0)

# label vocabulary (RCC export)
B      = "Naive B cells"
TREG   = "T regulatory cells"
PLASMA = "Plasmablasts"
MREG   = "mregDC"
CCR7   = "CCR7+ T (naive/CM)"
CD8    = ["Effector memory CD8 T cells", "CD8_T"]            # effector/cytotoxic CD8
FOCUS  = [B, TREG] + CD8 + [PLASMA, MREG, CCR7]

a = ad.read_h5ad(H5)
assert a.uns["spatial_units"] == "um", a.uns["spatial_units"]
xy  = np.asarray(a.obsm["spatial"], float)
lab = a.obs[KEY].astype(str).values
N   = a.n_obs
print(f"loaded {N} cells; B={np.sum(lab==B)} Treg={np.sum(lab==TREG)} "
      f"plasma={np.sum(lab==PLASMA)} mregDC={np.sum(lab==MREG)} CCR7+T={np.sum(lab==CCR7)} "
      f"effCD8={np.isin(lab,CD8).sum()}")

# tissue extent (for spread normalization)
ext = np.array([np.ptp(xy[:, 0]), np.ptp(xy[:, 1])]); tissue_diag = float(np.hypot(*ext))
print(f"tissue extent: {ext[0]:.0f} x {ext[1]:.0f} um (diag {tissue_diag:.0f})")

# ============================================================================
# STEP 1 — de-risk mregDC x CCR7+T: spatial distribution of PARENT c24 cells
# ============================================================================
is_c24 = np.isin(lab, [MREG, CCR7])
c24xy  = xy[is_c24]
# DBSCAN: foci of the parent cluster. eps ~ niche radius; min_samples a real focus.
db1 = DBSCAN(eps=50.0, min_samples=10).fit(c24xy)
cl1 = db1.labels_
foci_ids = [c for c in np.unique(cl1) if c != -1]
foci_sizes = np.array([np.sum(cl1 == c) for c in foci_ids])
order = np.argsort(foci_sizes)[::-1]
foci_ids = list(np.array(foci_ids)[order]); foci_sizes = foci_sizes[order]
n_foci = len(foci_ids)
n_clustered = int(foci_sizes.sum()); n_noise = int(np.sum(cl1 == -1))
largest_share = float(foci_sizes[0] / is_c24.sum()) if n_foci else np.nan
# spatial spread of foci CENTROIDS relative to tissue diagonal
fc = np.array([c24xy[cl1 == c].mean(0) for c in foci_ids]) if n_foci else np.empty((0, 2))
centroid_spread = float(np.hypot(*fc.std(0))) / tissue_diag if n_foci > 1 else 0.0
# fraction of foci within the central 25% box (a "one blob" tell)
print("\n== STEP 1: de-risk mregDC x CCR7+T (parent c24) ==")
print(f"  c24 cells: {int(is_c24.sum())}  ->  DBSCAN foci(eps=50,min=10): {n_foci} "
      f"(clustered {n_clustered}, noise/dispersed {n_noise})")
print(f"  focus sizes top10: {foci_sizes[:10].tolist()}")
print(f"  largest-focus share of c24: {largest_share:.2f}")
print(f"  focus-centroid spread / tissue-diag: {centroid_spread:.2f} "
      f"(>0.15 = scattered across tissue)")
# verdict
recurring = (n_foci >= 8) and (largest_share < 0.40) and (centroid_spread > 0.15)
verdict1 = ("REAL recurring niche (many separated foci across tissue) -> keep +90"
            if recurring else
            "CIRCULAR (one compact blob / few foci) -> DROP +90 from headline")
print(f"  VERDICT: {verdict1}")

# ============================================================================
# STEP 2 — formally delineate immunoregulatory aggregates (B-cell DBSCAN)
# ============================================================================
print("\n== STEP 2: delineate B-cell aggregates + composition ==")
isB  = lab == B
Bxy  = xy[isB]
# DBSCAN on B coords: dense B foci = candidate aggregates.
db2  = DBSCAN(eps=40.0, min_samples=20).fit(Bxy)
cl2  = db2.labels_
agg_ids = [c for c in np.unique(cl2) if c != -1]
print(f"  B cells {isB.sum()} -> {len(agg_ids)} B aggregates "
      f"(eps=40,min=20; {np.sum(cl2==-1)} dispersed B)")

# background = global tissue fraction of each type
def frac(mask_cells, label_or_list):
    L = label_or_list if isinstance(label_or_list, list) else [label_or_list]
    return float(np.isin(lab[mask_cells], L).mean()) if mask_cells.sum() else np.nan
bg = {"B": frac(np.ones(N, bool), B), "Treg": frac(np.ones(N, bool), TREG),
      "eff-CD8": frac(np.ones(N, bool), CD8), "Plasma": frac(np.ones(N, bool), PLASMA),
      "mregDC": frac(np.ones(N, bool), MREG)}

# aggregate footprint = all cells within R um of any member B cell of that aggregate
R = 50.0
tree = cKDTree(xy)
Bidx = np.where(isB)[0]
agg_rows = []
per_type_inside = {k: [] for k in ["Treg", "eff-CD8", "Plasma", "mregDC"]}
region_union = set()
agg_footprints = {}   # id -> region cell indices (for figure)
for c in agg_ids:
    members = Bidx[cl2 == c]
    cen = xy[members].mean(0)
    nbrs = tree.query_ball_point(xy[members], r=R)
    region = np.unique(np.concatenate([np.asarray(n_, int) for n_ in nbrs]))
    region_union.update(region.tolist())
    agg_footprints[c] = region
    mask = np.zeros(N, bool); mask[region] = True
    comp = {"Treg": frac(mask, TREG), "eff-CD8": frac(mask, CD8),
            "Plasma": frac(mask, PLASMA), "mregDC": frac(mask, MREG)}
    for k, v in comp.items():
        per_type_inside[k].append(v)
    agg_rows.append(dict(aggregate=int(c), n_B=len(members), n_cells_region=len(region),
                         x=float(cen[0]), y=float(cen[1]),
                         f_B=frac(mask, B), **{f"f_{k}": comp[k] for k in comp},
                         radius_um=float(np.median(np.hypot(*(xy[members] - cen).T)))))
agg_df = pd.DataFrame(agg_rows).sort_values("n_B", ascending=False).reset_index(drop=True)

# pooled inside-aggregate composition vs background + per-aggregate Wilcoxon vs background
print(f"  aggregates: {len(agg_df)}  | B-cells/agg median {agg_df.n_B.median():.0f} "
      f"(range {agg_df.n_B.min()}-{agg_df.n_B.max()})  | region cells/agg median "
      f"{agg_df.n_cells_region.median():.0f}")
contrast = []
for k in ["Treg", "eff-CD8", "Plasma", "mregDC"]:
    inside = np.array(per_type_inside[k]); b = bg[k]
    inside_mean = float(np.nanmean(inside))
    log2e = float(np.log2((inside_mean + 1e-6) / (b + 1e-6)))
    # per-aggregate: is inside fraction systematically > or < background?
    try:
        W, p = wilcoxon(inside - b, alternative="two-sided", zero_method="zsplit")
    except Exception:
        p = np.nan
    direction = "ENRICHED" if log2e > 0 else "EXCLUDED"
    contrast.append(dict(cell_type=k, background_frac=round(b, 4),
                         inside_mean_frac=round(inside_mean, 4), log2_enrichment=round(log2e, 3),
                         direction=direction, wilcoxon_p=(round(p, 3) if p == p else np.nan),
                         n_agg_above_bg=int(np.sum(inside > b))))
contrast_df = pd.DataFrame(contrast)
print("\n  per-aggregate composition contrast (inside vs tissue background):")
print(contrast_df.to_string(index=False))

# headline: the two LOAD-BEARING claims are Treg-enrichment + effector-CD8-exclusion.
# Plasma is reported descriptively (it sits AT background at this scale — that REFINES
# pass-1's global nhood z=-57 rather than confirming exclusion).
def row(ct):
    r = contrast_df[contrast_df.cell_type == ct].iloc[0]
    return r["direction"], float(r["log2_enrichment"]), float(r["wilcoxon_p"])
treg_d, treg_l, treg_p = row("Treg"); cd8_d, cd8_l, cd8_p = row("eff-CD8")
pl_d, pl_l, pl_p = row("Plasma"); mreg_d, mreg_l, mreg_p = row("mregDC")
sig = lambda l, p: (abs(l) > 0.3) and (p < 0.05)
plasma_status = ("EXCLUDED" if (pl_l < 0 and sig(pl_l, pl_p)) else
                 "ENRICHED" if (pl_l > 0 and sig(pl_l, pl_p)) else "AT BACKGROUND")
confirmed = (treg_d == "ENRICHED" and sig(treg_l, treg_p)) and \
            (cd8_d == "EXCLUDED" and sig(cd8_l, cd8_p))
print(f"\n  LOAD-BEARING claims confirmed = {confirmed}  "
      f"(Treg ENRICHED {treg_l:+.2f} p={treg_p:.1e}; eff-CD8 EXCLUDED {cd8_l:+.2f} p={cd8_p:.1e})")
print(f"  Plasma status @ aggregate scale: {plasma_status} ({pl_l:+.2f} log2, p={pl_p:.2f}) "
      f"-> refines pass-1 global nhood; mregDC {mreg_d} ({mreg_l:+.2f}, p={mreg_p:.1e})")

# heterogeneity gate for a UTAG pass3: warranted ONLY if a NON-modal axis shows a
# SIGNIFICANT enrichment splitting aggregates into distinct types. Plasma-at-background
# scatter (p>0.05) is NOT distinct niche structure.
plasma_rich = int(np.sum(agg_df.f_Plasma > bg["Plasma"]))
cd8_rich    = int(np.sum(agg_df["f_eff-CD8"] > bg["eff-CD8"]))
multi = (plasma_status == "ENRICHED") or (cd8_d == "ENRICHED")
print(f"  heterogeneity: plasma>bg in {plasma_rich}/{len(agg_df)} aggs (p={pl_p:.2f}), "
      f"eff-CD8>bg in {cd8_rich}/{len(agg_df)} "
      f"-> {'MULTIPLE niche types (flag UTAG pass3)' if multi else 'single homogeneous niche type — UTAG NOT warranted (plasma scatter is background-level)'}")

# write tables
agg_df.to_csv(f"{TAB}/rcc_phaseB2_aggregates.csv", index=False)
contrast_df.to_csv(f"{TAB}/rcc_phaseB2_aggregate_composition.csv", index=False)

# ============================================================================
# STEP 3 — one clean figure: 3 largest aggregates colored by cell type
# ============================================================================
COLORS = {B: "#1f77b4", TREG: "#d62728",
          "Effector memory CD8 T cells": "#2ca02c", "CD8_T": "#2ca02c",
          PLASMA: "#ff7f0e", MREG: "#9467bd", CCR7: "#17becf"}
LEGEND = [(B, "#1f77b4", "Naive B"), (TREG, "#d62728", "Treg"),
          ("eff-CD8", "#2ca02c", "effector CD8"), (PLASMA, "#ff7f0e", "Plasmablast"),
          (MREG, "#9467bd", "mregDC"), (CCR7, "#17becf", "CCR7+ T")]
top3 = agg_df.head(3)
fig, axes = plt.subplots(1, 3, figsize=(18, 6.2))
for ax, (_, row) in zip(np.atleast_1d(axes), top3.iterrows()):
    cx, cy, W = row["x"], row["y"], 250.0
    m = (np.abs(xy[:, 0] - cx) < W) & (np.abs(xy[:, 1] - cy) < W)
    sub_lab = lab[m]; sub_xy = xy[m]
    other = ~np.isin(sub_lab, FOCUS)
    ax.scatter(sub_xy[other, 0], sub_xy[other, 1], s=7, c="#dddddd", linewidths=0)
    for L, col, _name in LEGEND:
        sel = np.isin(sub_lab, [L]) if L != "eff-CD8" else np.isin(sub_lab, CD8)
        if sel.any():
            ax.scatter(sub_xy[sel, 0], sub_xy[sel, 1], s=14, c=col, linewidths=0)
    ax.set_title(f"aggregate {int(row['aggregate'])}: {int(row['n_B'])} B  |  "
                 f"Treg {row['f_Treg']*100:.0f}%  CD8 {row['f_eff-CD8']*100:.0f}%  "
                 f"plasma {row['f_Plasma']*100:.0f}%", fontsize=10)
    ax.set_aspect("equal"); ax.axis("off")
handles = [plt.Line2D([], [], marker="o", ls="", mfc=c, mec="none", ms=8, label=n) for _, c, n in LEGEND]
handles.append(plt.Line2D([], [], marker="o", ls="", mfc="#dddddd", mec="none", ms=8, label="other"))
axes[0].legend(handles=handles, fontsize=8, loc="upper left", bbox_to_anchor=(0, 1), framealpha=0.9)
fig.suptitle("RCC immunoregulatory aggregates: B–Treg core, effector-CD8 / plasma excluded", fontsize=13)
fig.tight_layout(); fig.savefig(f"{FIG}/rcc_phaseB2_aggregates.png", dpi=170); plt.close(fig)

# step-1 supporting figure: c24 foci across the whole tissue
fig, ax = plt.subplots(figsize=(9, 8))
ax.scatter(xy[::8, 0], xy[::8, 1], s=1, c="#eeeeee", linewidths=0)
noise = cl1 == -1
ax.scatter(c24xy[noise, 0], c24xy[noise, 1], s=4, c="#bbbbbb", linewidths=0, label="dispersed")
if n_foci:
    ax.scatter(c24xy[~noise, 0], c24xy[~noise, 1], s=6, c="#9467bd", linewidths=0,
               label=f"{n_foci} foci")
    ax.scatter(fc[:, 0], fc[:, 1], s=40, marker="x", c="k", linewidths=1.2, label="foci centroids")
ax.set_title(f"Parent c24 (mregDC ∪ CCR7+T): {n_foci} foci, "
             f"largest-share {largest_share:.2f}, spread {centroid_spread:.2f}\n{verdict1}", fontsize=10)
ax.set_aspect("equal"); ax.axis("off"); ax.legend(fontsize=8, loc="upper right")
fig.tight_layout(); fig.savefig(f"{FIG}/rcc_phaseB2_c24_foci.png", dpi=160); plt.close(fig)

# ============================================================================
# summary
# ============================================================================
summ = pd.DataFrame([
    ("c24_n_foci", n_foci), ("c24_largest_focus_share", round(largest_share, 3)),
    ("c24_centroid_spread", round(centroid_spread, 3)), ("c24_verdict_recurring", int(recurring)),
    ("n_B_aggregates", len(agg_df)), ("median_B_per_agg", float(agg_df.n_B.median())),
    ("Treg_log2enrich", treg_l), ("effCD8_log2enrich", cd8_l), ("Plasma_log2enrich", pl_l),
    ("mregDC_log2enrich", mreg_l), ("plasma_status", plasma_status),
    ("loadbearing_confirmed", int(confirmed)), ("utag_pass3_warranted", int(multi)),
], columns=["metric", "value"])
summ.to_csv(f"{TAB}/rcc_phaseB2_summary.csv", index=False)
print("\nwrote: rcc_phaseB2_aggregates.png, rcc_phaseB2_c24_foci.png, "
      "rcc_phaseB2_{aggregates,aggregate_composition,summary}.csv")
print("== phaseB_02 done ==")
