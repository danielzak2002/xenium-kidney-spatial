#!/usr/bin/env python
"""
phaseB_04_preview_validation.py — preview (PRCC) as INDEPENDENT, DIRECTIONAL replication
of the BIG RCC niche finding. Same Delaunay-graph approach as BIG (single Xenium section,
um). This is a directional check on SIGN/DIRECTION, NOT a second powered analysis — nhood
z scales with n, so preview z << BIG z; we read sign, not magnitude.

BIG reference effects to test (from phaseB_01/02 on RCC, 465k cells):
  B x Treg            nhood z +28   (ENRICHED)   ; aggregate Treg    log2 +1.22 (ENRICHED)
  B x effector-CD8    nhood z -117  (EXCLUDED)   ; aggregate eff-CD8 log2 -1.08 (EXCLUDED)
  B x mregDC          nhood z modest +           ; aggregate mregDC  log2 +1.03 (ENRICHED)
  mregDC x CCR7+T     nhood z +90   (ENRICHED, the de-risked recurring niche)

Preview wrinkles (handled): mregDC was FOLDED into 'myeloid/DC' in the preview export
(predates the R/11 mregDC-keep fix) -> RESTORED from cell_type. B here = 'Switched memory
B cells' (PRCC), not BIG's 'Naive B cells' (different subtype) -> non-replication is framed
as possibly PRCC-subtype/panel-specific, NOT a null.

  conda run -n spatial python py/phaseB_04_preview_validation.py
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, squidpy as sq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon

ROOT = os.environ.get("XENIUM_ROOT", os.getcwd())
H5   = os.path.join(ROOT, "outputs/objects/kidney_preview_PRCC.h5ad")
FIG  = os.path.join(ROOT, "outputs/figures"); TAB = os.path.join(ROOT, "outputs/tables")
OBJ  = os.path.join(ROOT, "outputs/objects")
rng  = np.random.default_rng(0)

B      = "Switched memory B cells"          # PRCC B subtype (BIG used 'Naive B cells')
TREG   = "T regulatory cells"
CD8    = ["Effector memory CD8 T cells"]    # no 'CD8_T' label in preview
PLASMA = "Plasmablasts"
MREG   = "mregDC"
CCR7   = "CCR7+ T (naive/CM)"
FOCUS  = [B, TREG] + CD8 + [PLASMA, MREG, CCR7]
MIN_POW = 200                               # population floor to call an arm "powered"
EPS, MIN_SAMP, R = 50.0, 15, 50.0           # B-aggregate DBSCAN (um) + footprint radius

a = ad.read_h5ad(H5)
assert a.uns["spatial_units"] == "um", a.uns["spatial_units"]
xy  = np.asarray(a.obsm["spatial"], float)
# working label = phase_b_label with mregDC RESTORED from cell_type (folded in export)
lab = a.obs["phase_b_label"].astype(str).values.copy()
ct  = a.obs["cell_type"].astype(str).values
lab[ct == MREG] = MREG
a.obs["nlab"] = pd.Categorical(lab)
N = a.n_obs

# ---- 1. POWER report ----
counts = {k: int(np.isin(lab, k if isinstance(k, list) else [k]).sum())
          for k in [B, TREG, "eff-CD8", PLASMA, MREG, CCR7]}
counts["eff-CD8"] = int(np.isin(lab, CD8).sum())
print(f"loaded {N} cells, {a.n_vars} genes (PRCC, Xenium um)")
print("population counts (power):")
for k, v in counts.items():
    print(f"  {k:28s} {v:6d}  {'OK' if v >= MIN_POW else 'SPARSE -> underpowered'}")

# ============================================================================
# 2. nhood enrichment — test BIG pairs (sign/direction)
# ============================================================================
sq.gr.spatial_neighbors(a, coord_type="generic", delaunay=True)
sq.gr.nhood_enrichment(a, cluster_key="nlab", seed=0, show_progress_bar=False)
cats = list(a.obs["nlab"].cat.categories)
Z = pd.DataFrame(a.uns["nlab_nhood_enrichment"]["zscore"], index=cats, columns=cats)
def zpair(i, j):
    return float(Z.loc[i, j]) if (i in Z.index and j in Z.columns) else np.nan
# BIG reference: (label_i, label_j, BIG_sign, BIG_value_str, powered)
pairs = [
    ("B x Treg",          B,    TREG, +1, "+28",  counts[B] >= MIN_POW and counts[TREG] >= MIN_POW),
    ("B x effector-CD8",  B,    CD8[0], -1, "-117", counts[B] >= MIN_POW and counts["eff-CD8"] >= MIN_POW),
    ("B x mregDC",        B,    MREG, +1, "modest+", counts[B] >= MIN_POW and counts[MREG] >= MIN_POW),
    ("mregDC x CCR7+T",   MREG, CCR7, +1, "+90",  counts[MREG] >= MIN_POW and counts[CCR7] >= MIN_POW),
]
nh_rows = []
print("\n== nhood enrichment: preview vs BIG (sign/direction) ==")
for name, i, j, big_sign, big_str, powered in pairs:
    z = zpair(i, j)
    match = (np.sign(z) == big_sign) if z == z else None
    nh_rows.append(dict(effect=name, big_nhood_z=big_str, preview_nhood_z=round(z, 2) if z == z else np.nan,
                        big_direction=("ENRICHED" if big_sign > 0 else "EXCLUDED"),
                        preview_direction=("ENRICHED" if (z == z and z > 0) else "EXCLUDED" if z == z else "NA"),
                        powered=powered, replicates=("YES" if match else "NO" if match is False else "NA")))
    print(f"  {name:20s} BIG {big_str:>7s} -> preview z={z:7.2f}  "
          f"{'[powered]' if powered else '[SPARSE]'}  replicates={nh_rows[-1]['replicates']}")

# ============================================================================
# 3. B-aggregate delineation + per-aggregate enrichment (mirror BIG pass-2)
# ============================================================================
agg_summary = {}
percell = []
isB = lab == B
print(f"\n== B-aggregate delineation (B n={isB.sum()}) ==")
if isB.sum() >= MIN_POW:
    Bidx = np.where(isB)[0]
    cl = DBSCAN(eps=EPS, min_samples=MIN_SAMP).fit(xy[Bidx]).labels_
    agg_ids = [c for c in np.unique(cl) if c != -1]
    print(f"  DBSCAN(eps={EPS},min={MIN_SAMP}) -> {len(agg_ids)} B aggregates "
          f"({np.sum(cl==-1)} dispersed B)")
    tree = cKDTree(xy)
    bg = {"Treg": float((lab == TREG).mean()), "eff-CD8": float(np.isin(lab, CD8).mean()),
          "mregDC": float((lab == MREG).mean()), "Plasma": float((lab == PLASMA).mean())}
    inside = {k: [] for k in bg}
    agg_id_full = np.full(N, -1, int)
    agg_rows = []
    for c in agg_ids:
        members = Bidx[cl == c]; agg_id_full[members] = c
        nbrs = tree.query_ball_point(xy[members], r=R)
        region = np.unique(np.concatenate([np.asarray(n_, int) for n_ in nbrs]))
        comp = {"Treg": float((lab[region] == TREG).mean()),
                "eff-CD8": float(np.isin(lab[region], CD8).mean()),
                "mregDC": float((lab[region] == MREG).mean()),
                "Plasma": float((lab[region] == PLASMA).mean())}
        for k in bg: inside[k].append(comp[k])
        cen = xy[members].mean(0)
        agg_rows.append(dict(aggregate=int(c), n_B=len(members), n_region=len(region),
                             x=float(cen[0]), y=float(cen[1]),
                             **{f"f_{k}": round(comp[k], 4) for k in comp}))
    agg_df = pd.DataFrame(agg_rows).sort_values("n_B", ascending=False).reset_index(drop=True)
    agg_df.to_csv(f"{TAB}/preview_phaseB_b_aggregates.csv", index=False)
    print("\n  per-aggregate composition (inside vs background), BIG sign in []:")
    big_agg = {"Treg": (+1, 1.22), "eff-CD8": (-1, -1.08), "mregDC": (+1, 1.03), "Plasma": (0, 0.04)}
    for k in ["Treg", "eff-CD8", "mregDC", "Plasma"]:
        ins = np.array(inside[k]); b = bg[k]; im = float(np.nanmean(ins))
        l2 = float(np.log2((im + 1e-6) / (b + 1e-6)))
        try:
            _, p = wilcoxon(ins - b, alternative="two-sided", zero_method="zsplit")
        except Exception:
            p = np.nan
        bs, bl = big_agg[k]
        match = (np.sign(round(l2, 3)) == bs) if bs != 0 else (abs(l2) < 0.3)
        agg_summary[k] = dict(bg=round(b, 4), inside=round(im, 4), log2=round(l2, 3),
                              p=round(p, 3) if p == p else np.nan,
                              big_log2=bl, replicates=("YES" if match else "NO"))
        print(f"    {k:9s} bg={b:.4f} inside={im:.4f} log2={l2:+.3f} (BIG {bl:+.2f}) "
              f"p={p if p==p else float('nan'):.3f} n_above_bg={int(np.sum(ins>b))}/{len(ins)} "
              f"replicates={agg_summary[k]['replicates']}")
    # figure: 3 largest B aggregates colored by cell type
    COLORS = {B: "#1f77b4", TREG: "#d62728", CD8[0]: "#2ca02c",
              PLASMA: "#ff7f0e", MREG: "#9467bd", CCR7: "#17becf"}
    top = agg_df.head(min(3, len(agg_df)))
    fig, axes = plt.subplots(1, len(top), figsize=(6 * len(top), 6.2)); axes = np.atleast_1d(axes)
    for ax, (_, row) in zip(axes, top.iterrows()):
        cx, cy, W = row["x"], row["y"], 250.0
        m = (np.abs(xy[:, 0] - cx) < W) & (np.abs(xy[:, 1] - cy) < W)
        sl, sx = lab[m], xy[m]
        ax.scatter(sx[~np.isin(sl, FOCUS), 0], sx[~np.isin(sl, FOCUS), 1], s=7, c="#dddddd", linewidths=0)
        for L, col in COLORS.items():
            sel = np.isin(sl, CD8) if L == CD8[0] else (sl == L)
            if sel.any(): ax.scatter(sx[sel, 0], sx[sel, 1], s=16, c=col, linewidths=0, label=L)
        ax.set_title(f"agg {int(row['aggregate'])}: {int(row['n_B'])} B | Treg {row['f_'+'Treg']*100:.0f}% "
                     f"CD8 {row['f_eff-CD8']*100:.0f}%", fontsize=9)
        ax.set_aspect("equal"); ax.axis("off")
    axes[0].legend(fontsize=7, loc="upper left", markerscale=1.5)
    fig.suptitle("Preview (PRCC) B aggregates — directional check of BIG B-Treg/CD8 finding", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{FIG}/preview_phaseB_b_aggregates.png", dpi=170); plt.close(fig)
    # persist per-cell (focus compartments)
    for i in np.where(np.isin(lab, FOCUS))[0]:
        percell.append((a.obs_names[i], lab[i], int(agg_id_full[i]), float(xy[i, 0]), float(xy[i, 1])))
else:
    print("  B insufficient -> aggregate arm skipped")

pc = pd.DataFrame(percell, columns=["cell_id", "phase_b_label", "b_agg_id", "x_um", "y_um"])
pc.to_parquet(f"{OBJ}/preview_focus_cells.parquet", index=False)
print(f"\npersisted {len(pc)} focus cells -> preview_focus_cells.parquet")

# ============================================================================
# replication summary table
# ============================================================================
rep_rows = []
for r in nh_rows:
    rep_rows.append(dict(level="nhood", effect=r["effect"], big=r["big_nhood_z"],
                         preview=r["preview_nhood_z"], powered=r["powered"], replicates=r["replicates"]))
for k, v in agg_summary.items():
    rep_rows.append(dict(level="B-aggregate", effect=f"{k} in B-agg", big=f"log2 {v['big_log2']:+.2f}",
                         preview=f"log2 {v['log2']:+.3f} (p={v['p']})", powered=(counts[B] >= MIN_POW),
                         replicates=v["replicates"]))
rep = pd.DataFrame(rep_rows)
rep.to_csv(f"{TAB}/preview_phaseB_replication_summary.csv", index=False)
print("\n== REPLICATION SUMMARY ==")
print(rep.to_string(index=False))
n_yes = (rep.replicates == "YES").sum(); n_test = rep.replicates.isin(["YES", "NO"]).sum()
print(f"\n{n_yes}/{n_test} testable effects replicate in direction. "
      "Sparse/NA arms are PRCC-subtype/panel-limited, not nulls.")
print("== phaseB_04 done ==")
