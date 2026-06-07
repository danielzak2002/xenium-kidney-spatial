#!/usr/bin/env python
"""
phaseB_05_rcc_plasma_myeloid.py — APPLES-TO-APPLES plasma analysis on RCC.

The cLN plasma-myeloid niche (plasma self-aggregate + co-localize with myeloid) was NEVER
measured on RCC; RCC plasma were only looked at vs the B aggregates + BCMA distance. So
"RCC plasma disperse" is not yet a fair comparison. Here we run the IDENTICAL cLN
plasma-myeloid analysis (same metrics, same physical scales) on RCC-BIG and RCC-preview,
each a single Xenium section -> one Delaunay graph (the section = one "slide"). BIG and
preview are kept SEPARATE (never pooled).

Identical to cln_phaseB (phaseB_03):
  myeloid = monocyte/macrophage/myeloid-DC lineage (excl. the distinct mregDC) — same
            compartment, mapped per label-vocabulary (RCC has no 'macrophage'/'monocyte'
            labels but 'Intermediate/Classical monocytes' + 'myeloid/DC' are that lineage).
  1. plasma x myeloid nhood_enrichment z.
  2. plasma self-aggregation: DBSCAN on PLASMA coords, eps=50 um (== cLN 0.05 mm), min=10.
  3. per-plasma-aggregate myeloid enrichment (inside vs section bg), log2 fold, R=50 um.

  conda run -n spatial python py/phaseB_05_rcc_plasma_myeloid.py
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, squidpy as sq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree

ROOT = os.environ.get("XENIUM_ROOT", os.getcwd())
FIG  = os.path.join(ROOT, "outputs/figures"); TAB = os.path.join(ROOT, "outputs/tables")
OBJ  = os.path.join(ROOT, "outputs/objects")
rng  = np.random.default_rng(0)

PLASMA = "Plasmablasts"
EPS, MIN_SAMP, R = 50.0, 10, 50.0     # um — matches cLN physical scale (0.05 mm)
DATASETS = [("RCC_BIG", "kidney_RCC_protein.h5ad"),
            ("RCC_preview", "kidney_preview_PRCC.h5ad")]

def myeloid_labels(all_labels):
    out = []
    for L in all_labels:
        ll = str(L).lower()
        if "mregdc" in ll:
            continue
        if "monocyte" in ll or "macrophage" in ll or ll == "myeloid/dc":
            out.append(L)
    return out

def run_dataset(name, fn):
    a = ad.read_h5ad(os.path.join(OBJ, fn))
    assert a.uns["spatial_units"] == "um", a.uns["spatial_units"]
    xy  = np.asarray(a.obsm["spatial"], float)
    lab = a.obs["phase_b_label"].astype(str).values
    myl = myeloid_labels(pd.unique(lab))
    niche = np.full(a.n_obs, "Other", dtype=object)
    niche[np.isin(lab, [PLASMA])] = "Plasma"
    niche[np.isin(lab, myl)]      = "Myeloid"
    n_pl = int(np.sum(niche == "Plasma")); n_my = int(np.sum(niche == "Myeloid"))
    f_my = n_my / a.n_obs
    print(f"\n== {name} ({a.n_obs} cells) ==")
    print(f"  myeloid labels: {myl}")
    print(f"  Plasma={n_pl}  Myeloid={n_my} ({f_my*100:.1f}% of section)")

    # 1. plasma x myeloid nhood z (single Delaunay graph)
    a.obs["nl"] = pd.Categorical(niche, categories=["Plasma", "Myeloid", "Other"])
    sq.gr.spatial_neighbors(a, coord_type="generic", delaunay=True)
    present = [c for c in ["Plasma", "Myeloid", "Other"] if (niche == c).sum() > 0]
    a.obs["nl2"] = pd.Categorical(niche, categories=present)
    sq.gr.nhood_enrichment(a, cluster_key="nl2", seed=0, show_progress_bar=False)
    Z = pd.DataFrame(a.uns["nl2_nhood_enrichment"]["zscore"], index=present, columns=present)
    z_pm = float(Z.loc["Plasma", "Myeloid"])

    # 2. plasma self-aggregation (DBSCAN, same physical scale as cLN)
    pl_idx = np.where(niche == "Plasma")[0]
    cl = DBSCAN(eps=EPS, min_samples=MIN_SAMP).fit(xy[pl_idx]).labels_
    agg_ids = [c for c in np.unique(cl) if c != -1]
    plasma_in_agg = int(np.sum(cl != -1))
    agg_rate = plasma_in_agg / n_pl if n_pl else np.nan

    # 3. per-aggregate myeloid enrichment (inside vs section bg)
    tree = cKDTree(xy); ins_fracs = []; sizes = []; agg_id_full = np.full(a.n_obs, -1, int)
    agg_rows = []
    for c in agg_ids:
        members = pl_idx[cl == c]; agg_id_full[members] = c
        nbrs = tree.query_ball_point(xy[members], r=R)
        region = np.unique(np.concatenate([np.asarray(n_, int) for n_ in nbrs]))
        f_in = float(np.mean(niche[region] == "Myeloid"))
        ins_fracs.append(f_in); sizes.append(len(members))
        cen = xy[members].mean(0)
        agg_rows.append(dict(dataset=name, aggregate=int(c), n_plasma=len(members),
                             n_region=len(region), myeloid_frac_inside=round(f_in, 4),
                             myeloid_bg=round(f_my, 4),
                             myeloid_log2enrich=round(float(np.log2((f_in+1e-6)/(f_my+1e-6))), 3),
                             x=float(cen[0]), y=float(cen[1])))
    myeloid_l2 = (float(np.log2((np.mean(ins_fracs)+1e-6)/(f_my+1e-6))) if ins_fracs else np.nan)
    print(f"  nhood z(Plasma x Myeloid) = {z_pm:.2f}")
    print(f"  plasma aggregates (eps={EPS},min={MIN_SAMP}) = {len(agg_ids)}  "
          f"(sizes med {int(np.median(sizes)) if sizes else 0}, max {max(sizes) if sizes else 0}); "
          f"aggregation rate = {agg_rate:.3f}")
    print(f"  myeloid log2-enrich in plasma aggs = {myeloid_l2:+.3f}" if myeloid_l2==myeloid_l2
          else "  myeloid log2-enrich = NA (no plasma aggregates)")

    # persist per-cell plasma+myeloid
    keep = np.where(np.isin(niche, ["Plasma", "Myeloid"]))[0]
    pc = pd.DataFrame({"cell_id": a.obs_names[keep], "dataset": name,
                       "phase_b_label": lab[keep], "niche_label": niche[keep].astype(str),
                       "plasma_agg_id": agg_id_full[keep], "x_um": xy[keep, 0], "y_um": xy[keep, 1]})
    pc.to_parquet(os.path.join(OBJ, f"{name}_plasma_myeloid_cells.parquet"), index=False)

    return dict(dataset=name, n_plasma=n_pl, n_myeloid=n_my, myeloid_frac=round(f_my, 4),
                plasma_myeloid_nhood_z=round(z_pm, 2), n_plasma_aggregates=len(agg_ids),
                plasma_agg_rate=round(agg_rate, 3) if agg_rate == agg_rate else np.nan,
                myeloid_log2enrich_in_plasma_aggs=round(myeloid_l2, 3) if myeloid_l2 == myeloid_l2 else np.nan), agg_rows

# ---- run RCC datasets (separate) -------------------------------------------
rows = []; all_agg = []
for name, fn in DATASETS:
    r, ar = run_dataset(name, fn); rows.append(r); all_agg += ar
rcc = pd.DataFrame(rows)
pd.DataFrame(all_agg).to_csv(f"{TAB}/rcc_plasma_aggregates.csv", index=False)

# ---- load cLN reference numbers (per-slide -> summary) ----------------------
cln = pd.read_csv(f"{TAB}/cln_phaseB_per_slide.csv")
cln_pb = cln[cln.n_plasma >= 10]                          # plasma-bearing slides only
def med_range(col):
    v = cln_pb[col].dropna()
    return (round(v.median(), 3), round(v.min(), 3), round(v.max(), 3)) if len(v) else (np.nan,)*3
# RCC-metric name -> cLN per-slide column name
CLN_COL = {"plasma_myeloid_nhood_z": "nhood_z_plasma_myeloid",
           "plasma_agg_rate": "plasma_agg_rate",
           "myeloid_log2enrich_in_plasma_aggs": "myeloid_log2enrich_in_plasma_aggs"}
cln_summary = {m: med_range(CLN_COL[m]) for m in CLN_COL}

# ---- side-by-side comparison table -----------------------------------------
comp = rcc.copy()
comp["context"] = "RCC (Xenium tumour)"
# cLN rows: overall (plasma-bearing slides) + by class
cln_by = pd.read_csv(f"{TAB}/cln_phaseB_by_class.csv")
def cln_class_mean(metric, cl):
    s = cln_by[(cln_by.metric == metric) & (cln_by.condition == cl)]
    return round(float(s["mean"].iloc[0]), 3) if len(s) and pd.notna(s["mean"].iloc[0]) else np.nan
cln_rows = []
cln_rows.append(dict(context="cLN (CosMx kidney)", dataset="cLN_overall(plasma slides)",
                     n_plasma=int(cln_pb.n_plasma.sum()), n_myeloid=int(cln_pb.n_myeloid.sum()),
                     myeloid_frac=round(cln_pb.myeloid_frac.mean(), 4),
                     plasma_myeloid_nhood_z=cln_summary["plasma_myeloid_nhood_z"][0],
                     n_plasma_aggregates=int(cln_pb.n_plasma_aggregates.sum()),
                     plasma_agg_rate=cln_summary["plasma_agg_rate"][0],
                     myeloid_log2enrich_in_plasma_aggs=cln_summary["myeloid_log2enrich_in_plasma_aggs"][0]))
for cl in ["III", "IV", "IV+V"]:
    cln_rows.append(dict(context="cLN (CosMx kidney)", dataset=f"cLN_{cl}", n_plasma=np.nan, n_myeloid=np.nan,
                         myeloid_frac=cln_class_mean("myeloid_frac", cl),
                         plasma_myeloid_nhood_z=cln_class_mean("nhood_z_plasma_myeloid", cl),
                         n_plasma_aggregates=cln_class_mean("n_plasma_aggregates", cl),
                         plasma_agg_rate=cln_class_mean("plasma_agg_rate", cl),
                         myeloid_log2enrich_in_plasma_aggs=cln_class_mean("myeloid_log2enrich_in_plasma_aggs", cl)))
comp = pd.concat([comp, pd.DataFrame(cln_rows)], ignore_index=True)
comp = comp[["context", "dataset", "n_plasma", "n_myeloid", "myeloid_frac",
             "plasma_myeloid_nhood_z", "n_plasma_aggregates", "plasma_agg_rate",
             "myeloid_log2enrich_in_plasma_aggs"]]
comp.to_csv(f"{TAB}/plasma_myeloid_rcc_vs_cln.csv", index=False)
print("\n== RCC vs cLN plasma-myeloid (apples-to-apples) ==")
print(comp.to_string(index=False))

# ---- verdict ----------------------------------------------------------------
print("\n== VERDICT ==")
cln_z_med, cln_z_lo, cln_z_hi = cln_summary["plasma_myeloid_nhood_z"]
cln_r_med, cln_r_lo, cln_r_hi = cln_summary["plasma_agg_rate"]
cln_l2_med, cln_l2_lo, cln_l2_hi = cln_summary["myeloid_log2enrich_in_plasma_aggs"]
print(f"  cLN (plasma slides): nhood z {cln_z_lo}..{cln_z_hi} (median {cln_z_med}); "
      f"agg rate {cln_r_lo}..{cln_r_hi}; myeloid-l2 {cln_l2_lo}..{cln_l2_hi}")
for _, r in rcc.iterrows():
    self_agg   = (r.plasma_agg_rate or 0) >= 0.05                 # plasma self-aggregate?
    myeloid_co = (r.plasma_myeloid_nhood_z > 3) and (r.myeloid_log2enrich_in_plasma_aggs > 0.2)
    if myeloid_co:
        v = "plasma co-localize AND recruit myeloid -> cLN niche SHARED (contrast collapses)"
    elif self_agg:
        v = ("plasma SELF-AGGREGATE (rate %.2f) but myeloid EXCLUDED (nhood z<0, log2~0) "
             "-> plasma-myeloid niche is cLN-SPECIFIC; contrast HOLDS (refined)" % r.plasma_agg_rate)
    else:
        v = "plasma neither aggregate nor recruit myeloid"
    print(f"  {r.dataset}: nhood z={r.plasma_myeloid_nhood_z}  agg rate={r.plasma_agg_rate}  "
          f"myeloid-l2={r.myeloid_log2enrich_in_plasma_aggs}\n     -> {v}")

# ---- figure: RCC values against cLN per-slide distribution ------------------
CCOL = {"control": "#4c72b0", "III": "#dd8452", "IV": "#c44e52", "IV+V": "#8172b3"}
metrics = [("plasma_myeloid_nhood_z", "nhood_z_plasma_myeloid", "plasma×myeloid nhood z"),
           ("plasma_agg_rate", "plasma_agg_rate", "plasma aggregation rate"),
           ("myeloid_log2enrich_in_plasma_aggs", "myeloid_log2enrich_in_plasma_aggs", "myeloid log2-enrich in plasma aggs")]
fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
for ax, (rcc_col, cln_col, ttl) in zip(axes, metrics):
    # cLN per-slide points (x=0), colored by class
    for _, s in cln_pb.iterrows():
        if pd.notna(s[cln_col]):
            ax.scatter(rng.normal(0, 0.05), s[cln_col], c=CCOL.get(s["condition"], "#888"),
                       s=55, edgecolor="k", linewidth=0.3, zorder=3)
    # RCC markers (x=1 BIG, x=2 preview)
    for j, (_, r) in enumerate(rcc.iterrows(), start=1):
        if pd.notna(r[rcc_col]):
            ax.scatter(j, r[rcc_col], marker="D", s=150, c="#222", zorder=4)
            ax.annotate(f"{r[rcc_col]:.2f}", (j, r[rcc_col]), textcoords="offset points",
                        xytext=(8, 0), fontsize=9)
    ax.axhline(0, color="#999", lw=0.6, ls="--")
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["cLN slides\n(by class)", "RCC BIG", "RCC preview"])
    ax.set_title(ttl, fontsize=11)
handles = [plt.Line2D([], [], marker="o", ls="", mfc=CCOL[c], mec="k", ms=8, label=f"cLN {c}")
           for c in ["control", "III", "IV", "IV+V"]]
handles.append(plt.Line2D([], [], marker="D", ls="", mfc="#222", mec="none", ms=10, label="RCC (single section)"))
axes[0].legend(handles=handles, fontsize=8, loc="upper left")
fig.suptitle("Apples-to-apples plasma–myeloid: RCC (Xenium) vs cLN (CosMx) — identical metrics & 50µm scale", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/plasma_myeloid_rcc_vs_cln.png", dpi=160); plt.close(fig)
print("\nwrote: plasma_myeloid_rcc_vs_cln.png + plasma_myeloid_rcc_vs_cln.csv + rcc_plasma_aggregates.csv")
print("== phaseB_05 done ==")
