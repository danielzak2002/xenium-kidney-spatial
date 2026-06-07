#!/usr/bin/env python
"""
phaseB_03_cln_niches.py — cLN CosMx spatial pass (squidpy). Scoped to the TRUSTED
compartments only: PLASMA (plasmablast) + MYELOID (macrophage / myeloid-DC / monocyte).
B (n unreliable) and T (not recovered) are EXCLUDED from every spatial claim here.

CosMx cohort hazards handled explicitly:
  * PER-SLIDE graphs — cells on different slides are NOT neighbours. Build the neighbour
    graph WITHIN each of the 14 slides, compute per slide, THEN summarise by LN class.
    Never one graph on the merged object.
  * mm coordinates — adjacency via Delaunay (scale-invariant); co-occurrence bins set
    EXPLICITLY in mm.
  * Patient / n confound — control 4 / III 4 / IV 5 / IV+V 1 slides, 8 patients (SLE8 x3
    = pseudo-replication). Per-slide values are primary; class summaries are DESCRIPTIVE /
    exploratory, NOT inferential. IV+V (n=1) carries no class-level conclusion.

Cross-context comparability: plasma metrics are emitted as NORMALISED values (nhood z,
log2 fold, aggregation rate) — never absolute mm distances — so they line up against the
RCC (um) plasma numbers in the later deep-dive. Per-cell trusted-compartment table is
persisted (parquet) so the deep-dive reloads without recomputing.

  conda run -n spatial python py/phaseB_03_cln_niches.py
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, squidpy as sq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree

ROOT = os.environ.get("XENIUM_ROOT", os.getcwd())
H5   = os.path.join(ROOT, "outputs/objects/cln_cosmx.h5ad")
FIG  = os.path.join(ROOT, "outputs/figures"); TAB = os.path.join(ROOT, "outputs/tables")
OBJ  = os.path.join(ROOT, "outputs/objects")
KEY  = "phase_b_label"
rng  = np.random.default_rng(0)

PLASMA_LBL  = ["plasmablast"]
MYELOID_LBL = ["macrophage", "myeloid/DC", "monocyte"]
CLASS_ORDER = ["control", "III", "IV", "IV+V"]          # disease severity axis
MIN_CELLS   = 10                                         # min per category for a z
EPS_MM      = 0.05                                       # DBSCAN eps / footprint radius (mm = 50um)
MIN_SAMP    = 10                                         # DBSCAN min_samples for a plasma aggregate

a = ad.read_h5ad(H5)
assert a.uns["spatial_units"] == "mm", a.uns["spatial_units"]
xy_all = np.asarray(a.obsm["spatial"], float)
lab    = a.obs[KEY].astype(str).values
# coarse 3-level niche label (Plasma / Myeloid / Other) for the spatial-graph stats
niche = np.full(a.n_obs, "Other", dtype=object)
niche[np.isin(lab, PLASMA_LBL)]  = "Plasma"
niche[np.isin(lab, MYELOID_LBL)] = "Myeloid"
a.obs["niche_label"] = pd.Categorical(niche, categories=["Plasma", "Myeloid", "Other"])
slides  = a.obs["sample"].astype(str).values
classes = a.obs["condition"].astype(str).values
slide_class = (pd.DataFrame({"sample": slides, "condition": classes})
               .drop_duplicates().set_index("sample")["condition"].to_dict())
slide_ids = list(pd.unique(slides))
print(f"loaded {a.n_obs} cells, {a.n_vars} genes, {len(slide_ids)} slides; "
      f"Plasma={np.sum(niche=='Plasma')} Myeloid={np.sum(niche=='Myeloid')}")
print("class -> n slides:", {c: sum(v == c for v in slide_class.values()) for c in CLASS_ORDER})

# co-occurrence distance bins, mm (cell-niche scale 10um..0.5mm)
COOC_BINS = np.linspace(0.01, 0.5, 20)

# ---------------------------------------------------------------------------
per_slide   = []           # one row per slide: abundances + nhood z + aggregation
cooc_curves = {}           # slide -> plasma->myeloid co-occurrence ratio vs distance
agg_rows    = []           # one row per plasma aggregate
percell     = []           # trusted-compartment per-cell records (for deep-dive reload)

for s in slide_ids:
    m = slides == s
    cls = slide_class[s]
    sub = a[m].copy()
    sub.obsm["spatial"] = np.asarray(sub.obsm["spatial"], float)
    nl = sub.obs["niche_label"].values
    xy = sub.obsm["spatial"]
    n_pl = int(np.sum(nl == "Plasma")); n_my = int(np.sum(nl == "Myeloid"))
    f_pl = n_pl / sub.n_obs; f_my = n_my / sub.n_obs

    # ---- 1. per-slide spatial graph (Delaunay, scale-invariant) ----
    sq.gr.spatial_neighbors(sub, coord_type="generic", delaunay=True)

    # ---- nhood enrichment plasma x myeloid (z) ----
    z_pm = np.nan
    if n_pl >= MIN_CELLS and n_my >= MIN_CELLS:
        # restrict categorical to present categories so squidpy doesn't choke on empties
        present = [c for c in ["Plasma", "Myeloid", "Other"]
                   if (nl == c).sum() > 0]
        sub.obs["nl2"] = pd.Categorical(nl, categories=present)
        sq.gr.nhood_enrichment(sub, cluster_key="nl2", seed=0, show_progress_bar=False)
        Z = pd.DataFrame(sub.uns["nl2_nhood_enrichment"]["zscore"], index=present, columns=present)
        if "Plasma" in Z.index and "Myeloid" in Z.columns:
            z_pm = float(Z.loc["Plasma", "Myeloid"])

    # ---- 2. co-occurrence over mm distance (keep all P+M, subsample Other) ----
    if n_pl >= MIN_CELLS and n_my >= MIN_CELLS:
        keep = np.where((nl == "Plasma") | (nl == "Myeloid"))[0]
        oth  = np.where(nl == "Other")[0]
        cap  = max(0, 25000 - len(keep))
        if len(oth) > cap:
            oth = rng.choice(oth, cap, replace=False)
        csub = sub[np.concatenate([keep, oth])].copy()
        csub.obs["nl2"] = csub.obs["niche_label"].cat.remove_unused_categories()
        try:
            sq.gr.co_occurrence(csub, cluster_key="nl2", interval=COOC_BINS, show_progress_bar=False)
            occ = csub.uns["nl2_co_occurrence"]["occ"]          # (cat, cat, n_bins-1)
            cats = list(csub.obs["nl2"].cat.categories)
            if "Plasma" in cats and "Myeloid" in cats:
                cooc_curves[s] = occ[cats.index("Plasma"), cats.index("Myeloid"), :]
        except Exception as e:
            print(f"  [{s}] co-occurrence skipped: {e}")

    # ---- 3. plasma aggregation (DBSCAN, mm) + myeloid enrichment inside ----
    n_agg = 0; agg_myeloid_l2 = np.nan; plasma_in_agg = 0
    pl_idx = np.where(nl == "Plasma")[0]
    agg_id_slide = np.full(sub.n_obs, -1, int)
    if len(pl_idx) >= MIN_SAMP:
        db = DBSCAN(eps=EPS_MM, min_samples=MIN_SAMP).fit(xy[pl_idx])
        cl = db.labels_
        tree = cKDTree(xy)
        bg_my = f_my                                            # slide background myeloid frac
        ins_fracs = []
        for c in [c for c in np.unique(cl) if c != -1]:
            members = pl_idx[cl == c]
            agg_id_slide[members] = c
            nbrs = tree.query_ball_point(xy[members], r=EPS_MM)
            region = np.unique(np.concatenate([np.asarray(n_, int) for n_ in nbrs]))
            f_in = float(np.mean(nl[region] == "Myeloid"))
            ins_fracs.append(f_in)
            cen = xy[members].mean(0)
            agg_rows.append(dict(slide=s, condition=cls, aggregate=int(c), n_plasma=len(members),
                                 n_cells_region=len(region), myeloid_frac_inside=round(f_in, 4),
                                 myeloid_bg=round(bg_my, 4),
                                 myeloid_log2enrich=round(float(np.log2((f_in + 1e-6) / (bg_my + 1e-6))), 3),
                                 x=float(cen[0]), y=float(cen[1])))
        n_agg = len(ins_fracs)
        plasma_in_agg = int(np.sum(cl != -1))
        if n_agg:
            agg_myeloid_l2 = float(np.log2((np.mean(ins_fracs) + 1e-6) / (bg_my + 1e-6)))

    plasma_agg_rate = (plasma_in_agg / n_pl) if n_pl else np.nan

    # ---- orthogonal IF check: plasma/myeloid in CD45+ immune (not ambient tubule) ----
    cd45 = sub.obs["Mean.CD45"].values; panck = sub.obs["Mean.PanCK"].values
    cd45_pl = float(np.nanmean(cd45[nl == "Plasma"])) if n_pl else np.nan
    cd45_my = float(np.nanmean(cd45[nl == "Myeloid"])) if n_my else np.nan
    cd45_ep = float(np.nanmean(cd45[nl == "Other"]))
    panck_pl = float(np.nanmean(panck[nl == "Plasma"])) if n_pl else np.nan

    per_slide.append(dict(slide=s, condition=cls, n_cells=sub.n_obs,
        n_plasma=n_pl, n_myeloid=n_my, plasma_frac=round(f_pl, 4), myeloid_frac=round(f_my, 4),
        nhood_z_plasma_myeloid=round(z_pm, 2) if z_pm == z_pm else np.nan,
        n_plasma_aggregates=n_agg, plasma_agg_rate=round(plasma_agg_rate, 3) if plasma_agg_rate == plasma_agg_rate else np.nan,
        myeloid_log2enrich_in_plasma_aggs=round(agg_myeloid_l2, 3) if agg_myeloid_l2 == agg_myeloid_l2 else np.nan,
        cd45_plasma=round(cd45_pl, 1) if cd45_pl == cd45_pl else np.nan,
        cd45_myeloid=round(cd45_my, 1) if cd45_my == cd45_my else np.nan,
        cd45_other=round(cd45_ep, 1), panck_plasma=round(panck_pl, 1) if panck_pl == panck_pl else np.nan))

    # persist trusted-compartment per-cell records
    tr = np.where((nl == "Plasma") | (nl == "Myeloid"))[0]
    for i in tr:
        percell.append((sub.obs_names[i], s, cls, lab[m][i], str(nl[i]),
                        int(agg_id_slide[i]), float(xy[i, 0]), float(xy[i, 1]),
                        float(cd45[i]), float(panck[i])))
    print(f"  [{s:22s} {cls:7s}] P={n_pl:4d} M={n_my:5d}  z(PxM)={z_pm:6.1f}  "
          f"agg={n_agg:2d} rate={plasma_agg_rate if plasma_agg_rate==plasma_agg_rate else 0:.2f} "
          f"myeloid-l2={agg_myeloid_l2 if agg_myeloid_l2==agg_myeloid_l2 else float('nan'):+.2f}")

ps = pd.DataFrame(per_slide)
ps.to_csv(f"{TAB}/cln_phaseB_per_slide.csv", index=False)
pd.DataFrame(agg_rows).to_csv(f"{TAB}/cln_phaseB_plasma_aggregates.csv", index=False)

# ---- per-cell persistence (parquet, git-ignored) for the deep-dive ----
pc = pd.DataFrame(percell, columns=["cell_id", "slide", "condition", "phase_b_label",
                                    "niche_label", "plasma_agg_id", "x_mm", "y_mm",
                                    "Mean_CD45", "Mean_PanCK"])
pc.to_parquet(f"{OBJ}/cln_plasma_myeloid_cells.parquet", index=False)
print(f"\npersisted {len(pc)} trusted-compartment cells -> cln_plasma_myeloid_cells.parquet")

# ============================================================================
# summarise by LN class (DESCRIPTIVE — underpowered, patient-confounded)
# ============================================================================
def by_class(col):
    g = ps.groupby("condition")[col]
    return g.agg(["mean", "std", "count"]).reindex(CLASS_ORDER)

print("\n== BY-CLASS SUMMARY (descriptive; control 4 / III 4 / IV 5 / IV+V 1 slides) ==")
clsum_rows = []
for col in ["plasma_frac", "myeloid_frac", "nhood_z_plasma_myeloid",
            "n_plasma_aggregates", "plasma_agg_rate", "myeloid_log2enrich_in_plasma_aggs"]:
    bc = by_class(col)
    for cl in CLASS_ORDER:
        clsum_rows.append(dict(metric=col, condition=cl, n_slides=int(bc.loc[cl, "count"]),
                               mean=round(float(bc.loc[cl, "mean"]), 3) if bc.loc[cl, "count"] else np.nan,
                               std=round(float(bc.loc[cl, "std"]), 3) if bc.loc[cl, "count"] > 1 else np.nan))
    print(f"  {col:38s} " + "  ".join(
        f"{cl}={bc.loc[cl,'mean']:.3f}(n{int(bc.loc[cl,'count'])})" if bc.loc[cl, "count"] else f"{cl}=NA"
        for cl in CLASS_ORDER))
pd.DataFrame(clsum_rows).to_csv(f"{TAB}/cln_phaseB_by_class.csv", index=False)

# co-occurrence aggregated by class
bins_mid = (COOC_BINS[:-1] + COOC_BINS[1:]) / 2
cooc_by_class = {}
for cl in CLASS_ORDER:
    curves = [cooc_curves[s] for s in slide_ids if slide_class[s] == cl and s in cooc_curves]
    if curves:
        cooc_by_class[cl] = np.nanmean(np.vstack(curves), 0)
cdf = pd.DataFrame({"dist_mm": bins_mid})
for cl, v in cooc_by_class.items():
    cdf[cl] = v
cdf.to_csv(f"{TAB}/cln_phaseB_cooccurrence_by_class.csv", index=False)

# ============================================================================
# FIGURES
# ============================================================================
CCOL = {"control": "#4c72b0", "III": "#dd8452", "IV": "#c44e52", "IV+V": "#8172b3"}
nslide = {c: sum(v == c for v in slide_class.values()) for c in CLASS_ORDER}
caveat = ("descriptive only — control n=%d, III n=%d, IV n=%d, IV+V n=%d (single slide); "
          "patients contribute multiple regions (pseudo-replication)") % tuple(nslide[c] for c in CLASS_ORDER)

# Fig 1: disease-axis — per-slide points by class for the 4 key metrics
metrics = [("nhood_z_plasma_myeloid", "plasma×myeloid nhood z"),
           ("plasma_agg_rate", "plasma aggregation rate"),
           ("myeloid_log2enrich_in_plasma_aggs", "myeloid log2-enrich in plasma aggs"),
           ("plasma_frac", "plasma fraction of slide")]
fig, axes = plt.subplots(1, 4, figsize=(20, 5))
for ax, (col, ttl) in zip(axes, metrics):
    for j, cl in enumerate(CLASS_ORDER):
        vals = ps.loc[ps.condition == cl, col].dropna().values
        ax.scatter(np.full(len(vals), j) + rng.normal(0, 0.05, len(vals)), vals,
                   c=CCOL[cl], s=70, edgecolor="k", linewidth=0.4, zorder=3)
        if len(vals):
            ax.hlines(np.mean(vals), j - 0.25, j + 0.25, color=CCOL[cl], lw=2.5, zorder=2)
    ax.set_xticks(range(4)); ax.set_xticklabels([f"{c}\n(n{nslide[c]})" for c in CLASS_ORDER])
    ax.set_title(ttl, fontsize=11); ax.axhline(0, color="#999", lw=0.6, ls="--")
fig.suptitle("cLN disease axis: plasma–myeloid spatial metrics per slide\n" + caveat, fontsize=11)
fig.tight_layout(); fig.savefig(f"{FIG}/cln_phaseB_disease_axis.png", dpi=160); plt.close(fig)

# Fig 2: co-occurrence by class
fig, ax = plt.subplots(figsize=(7.5, 5.2))
for cl in CLASS_ORDER:
    if cl in cooc_by_class:
        ax.plot(bins_mid, cooc_by_class[cl], color=CCOL[cl], lw=2, label=f"{cl} (n{nslide[cl]})")
ax.axhline(1, color="#999", lw=0.7, ls="--")
ax.set_xlabel("distance (mm)"); ax.set_ylabel("plasma→myeloid co-occurrence ratio")
ax.set_title("cLN plasma–myeloid co-occurrence vs distance, by class\n" + caveat, fontsize=9)
ax.legend(); fig.tight_layout(); fig.savefig(f"{FIG}/cln_phaseB_cooccurrence.png", dpi=160); plt.close(fig)

# Fig 3: representative slide — the one with most plasma aggregates
rep = ps.sort_values("n_plasma_aggregates", ascending=False).iloc[0]["slide"]
m = slides == rep; sub_lab = niche[m]; sub_xy = xy_all[m]
fig, ax = plt.subplots(figsize=(9, 8))
ax.scatter(sub_xy[sub_lab == "Other", 0], sub_xy[sub_lab == "Other", 1], s=2, c="#eeeeee", linewidths=0)
ax.scatter(sub_xy[sub_lab == "Myeloid", 0], sub_xy[sub_lab == "Myeloid", 1], s=8, c="#2ca02c", linewidths=0, label="myeloid")
ax.scatter(sub_xy[sub_lab == "Plasma", 0], sub_xy[sub_lab == "Plasma", 1], s=12, c="#ff7f0e", linewidths=0, label="plasma")
ax.set_title(f"representative slide {rep} ({slide_class[rep]}): plasma + myeloid", fontsize=11)
ax.set_aspect("equal"); ax.axis("off"); ax.legend(markerscale=2, fontsize=9)
fig.tight_layout(); fig.savefig(f"{FIG}/cln_phaseB_rep_slide.png", dpi=160); plt.close(fig)

print("\nwrote: cln_phaseB_{disease_axis,cooccurrence,rep_slide}.png + "
      "cln_phaseB_{per_slide,plasma_aggregates,by_class,cooccurrence_by_class}.csv")
print("== phaseB_03 done ==")
