#!/usr/bin/env python
"""
whitepaper_groupE_cln.py — Group E (Section 6.2, cLN plasma-myeloid niche) PER-SLIDE.

The collective cLN figures are REPLACED by individual per-slide sets: for each plasma-
bearing slide (+ one control) a 3-panel figure — (1) plasma+myeloid on tissue, (2)
delineated plasma aggregates with myeloid overlaid, (3) a representative plasma-myeloid
focus crop — then an across-slide POINTS summary and the per-slide table.

From the committed cLN h5ad (coords + phase_b_label + slide/condition) and the committed
per-slide metric table. Plasma DBSCAN recomputed per slide at the Phase-B physical scale
(eps=0.05 mm, min=10); mm coordinates throughout.

  conda run -n spatial python py/whitepaper_groupE_cln.py
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree

ROOT = os.environ.get("XENIUM_ROOT", os.getcwd())
OBJ = os.path.join(ROOT, "outputs/objects"); TAB = os.path.join(ROOT, "outputs/tables")
FIG = os.path.join(ROOT, "outputs/figures/whitepaper"); os.makedirs(FIG, exist_ok=True)
EPS, MIN_SAMP, R = 0.05, 10, 0.05      # mm — Phase-B plasma-aggregate scale

PLASMA = ["plasmablast"]; MYELOID = ["macrophage", "myeloid/DC", "monocyte"]
PCOL, MCOL, TCOL = "#ff7f0e", "#2ca02c", "#e8e8e8"
CCOL = {"control": "#4c72b0", "III": "#dd8452", "IV": "#c44e52", "IV+V": "#8172b3"}
# (slide, class) — 5 plasma-bearing + 1 control (most-myeloid, ~0 plasma → clean contrast)
SLIDES = [("SP21_213_R1080_S1", "III"), ("SP19_1139_R1080_S3", "IV"),
          ("SP20_642_R1080_S3", "IV+V"), ("SP20_10838_R1080_S1", "III"),
          ("SP18_8471_R1087_S2", "III"), ("SMI0016C_SP17SP19", "control")]

a = ad.read_h5ad(os.path.join(OBJ, "cln_cosmx.h5ad"))
xy_all = np.asarray(a.obsm["spatial"], float)
lab_all = a.obs["phase_b_label"].astype(str).values
slide_all = a.obs["sample"].astype(str).values
print(f"loaded {a.n_obs} cells")

def plasma_aggregates(xy, isP):
    Pidx = np.where(isP)[0]
    if Pidx.size < MIN_SAMP:
        return [], np.full(len(isP), -1, int)
    cl = DBSCAN(eps=EPS, min_samples=MIN_SAMP).fit(xy[Pidx]).labels_
    agg_id = np.full(len(isP), -1, int); aggs = []
    for c in [c for c in np.unique(cl) if c != -1]:
        members = Pidx[cl == c]; agg_id[members] = c
        aggs.append(dict(c=c, members=members, cen=xy[members].mean(0), n=len(members)))
    aggs = sorted(aggs, key=lambda d: d["n"], reverse=True)
    return aggs, agg_id

# ---- per-slide 3-panel figures ---------------------------------------------
for slide, cls in SLIDES:
    sm = slide_all == slide
    xy = xy_all[sm]; lab = lab_all[sm]
    isP = np.isin(lab, PLASMA); isM = np.isin(lab, MYELOID)
    aggs, agg_id = plasma_aggregates(xy, isP)
    nP, nM, nA = int(isP.sum()), int(isM.sum()), len(aggs)
    ext = max(np.ptp(xy[:, 0]), np.ptp(xy[:, 1]))
    fig, ax = plt.subplots(1, 3, figsize=(21, 7))
    # P1: plasma + myeloid on tissue
    ax[0].scatter(xy[:, 0], xy[:, 1], s=1, c=TCOL, linewidths=0, rasterized=True)
    ax[0].scatter(xy[isM, 0], xy[isM, 1], s=4, c=MCOL, linewidths=0, label=f"myeloid (n={nM})")
    ax[0].scatter(xy[isP, 0], xy[isP, 1], s=7, c=PCOL, linewidths=0, label=f"plasma (n={nP})")
    ax[0].set_title(f"{slide} [{cls}] — plasma + myeloid on tissue", fontsize=10)
    ax[0].legend(fontsize=8, loc="upper right", markerscale=2)
    # P2: delineated plasma aggregates + myeloid overlaid
    ax[1].scatter(xy[:, 0], xy[:, 1], s=1, c=TCOL, linewidths=0, rasterized=True)
    ax[1].scatter(xy[isM, 0], xy[isM, 1], s=4, c=MCOL, linewidths=0)
    inA = agg_id >= 0
    ax[1].scatter(xy[isP & ~inA, 0], xy[isP & ~inA, 1], s=5, c="#ffd9a8", linewidths=0)  # dispersed plasma
    ax[1].scatter(xy[inA, 0], xy[inA, 1], s=10, c=PCOL, linewidths=0)                     # aggregated plasma
    for d in aggs:
        ax[1].add_patch(plt.Circle(d["cen"], R, fill=False, ec="#222", lw=1.0, alpha=0.8))
    ax[1].set_title(f"delineated plasma aggregates (n={nA}) + myeloid", fontsize=10)
    # P3: representative focus crop
    if aggs:
        cx, cy = aggs[0]["cen"]; W = 0.15
        mm = (np.abs(xy[:, 0] - cx) < W) & (np.abs(xy[:, 1] - cy) < W)
        sx, sl_ = xy[mm], lab[mm]
        ax[2].scatter(sx[:, 0], sx[:, 1], s=10, c=TCOL, linewidths=0)
        ax[2].scatter(sx[np.isin(sl_, MYELOID), 0], sx[np.isin(sl_, MYELOID), 1], s=40, c=MCOL, linewidths=0.2, edgecolor="k", label="myeloid")
        ax[2].scatter(sx[np.isin(sl_, PLASMA), 0], sx[np.isin(sl_, PLASMA), 1], s=40, c=PCOL, linewidths=0.2, edgecolor="k", label="plasma")
        ax[2].set_title(f"representative plasma–myeloid focus ({aggs[0]['n']} plasma)", fontsize=10)
        ax[2].legend(fontsize=8, loc="upper right", markerscale=1.3)
    else:
        ax[2].scatter(xy[:, 0], xy[:, 1], s=1, c=TCOL, linewidths=0)
        ax[2].scatter(xy[isM, 0], xy[isM, 1], s=4, c=MCOL, linewidths=0)
        ax[2].scatter(xy[isP, 0], xy[isP, 1], s=12, c=PCOL, linewidths=0)
        ax[2].set_title("no plasma aggregates (control / plasma-sparse)", fontsize=10)
    for x_ in ax:
        x_.set_aspect("equal"); x_.set_xticks([]); x_.set_yticks([])
    fig.suptitle(f"cLN slide {slide} ({cls}): {nP} plasma, {nM} myeloid, {nA} plasma aggregates", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fn = f"qcE_slide_{slide}.png"; fig.savefig(os.path.join(FIG, fn), dpi=140); plt.close(fig)
    print(f"wrote {fn}  (P={nP} M={nM} aggs={nA})")

# ---- across-slide POINTS summary (not collective means) --------------------
ps = pd.read_csv(os.path.join(TAB, "cln_phaseB_per_slide.csv"))
pb = ps[ps.n_plasma >= 10].copy().sort_values("nhood_z_plasma_myeloid", ascending=False)
fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
for j, (col, ttl, ylab) in enumerate([
        ("nhood_z_plasma_myeloid", "plasma×myeloid neighborhood enrichment (per slide)", "nhood z"),
        ("myeloid_log2enrich_in_plasma_aggs", "myeloid enrichment inside plasma aggregates (per slide)", "log2 enrichment")]):
    d = pb.dropna(subset=[col])
    for i, (_, r) in enumerate(d.iterrows()):
        ax[j].scatter(i, r[col], s=90, c=CCOL.get(r.condition, "#888"), edgecolor="k", linewidth=0.4, zorder=3)
        ax[j].annotate(r.slide.split("_")[0], (i, r[col]), fontsize=6.5, rotation=45,
                       textcoords="offset points", xytext=(3, 4))
    ax[j].axhline(0, color="#999", lw=0.6, ls="--")
    ax[j].set_xticks([]); ax[j].set_ylabel(ylab); ax[j].set_title(ttl, fontsize=10)
handles = [plt.Line2D([], [], marker="o", ls="", mfc=CCOL[c], mec="k", ms=8, label=c) for c in ["control", "III", "IV", "IV+V"]]
ax[0].legend(handles=handles, fontsize=8, loc="upper right", title="ISN class")
fig.suptitle("cLN plasma–myeloid niche across slides — per-slide points (consistency + variability); "
             "descriptive only (patient-confounded, IV+V n=1)", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(os.path.join(FIG, "qcE_across_slide_summary.png"), dpi=150); plt.close(fig)
print("wrote qcE_across_slide_summary.png")

# ---- per-slide table --------------------------------------------------------
tcols = ["slide", "condition", "n_plasma", "n_myeloid", "nhood_z_plasma_myeloid",
         "n_plasma_aggregates", "plasma_agg_rate", "myeloid_log2enrich_in_plasma_aggs"]
tt = ps[tcols].copy().sort_values(["condition", "n_plasma"], ascending=[True, False])
tt.columns = ["slide", "class", "n plasma", "n myeloid", "plasma×myeloid z", "n aggs", "agg rate", "myeloid log2"]
fig, ax = plt.subplots(figsize=(13, 5)); ax.axis("off")
tb = ax.table(cellText=tt.round(3).fillna("—").values, colLabels=tt.columns, cellLoc="center", loc="center")
tb.auto_set_font_size(False); tb.set_fontsize(8); tb.scale(1, 1.3)
for j in range(len(tt.columns)):
    tb[0, j].set_facecolor("#4c72b0"); tb[0, j].set_text_props(color="white", weight="bold")
for i, cls in enumerate(tt["class"].values, start=1):
    for j in range(len(tt.columns)):
        tb[i, j].set_facecolor(CCOL.get(cls, "#fff") + "22" if False else "#f4f4f4" if i % 2 else "#ffffff")
ax.set_title("cLN per-slide plasma–myeloid metrics (all 14 slides)", fontsize=11, pad=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcE_per_slide_table.png"), dpi=150); plt.close(fig)
print("wrote qcE_per_slide_table.png")
print("== Group E done ==")
