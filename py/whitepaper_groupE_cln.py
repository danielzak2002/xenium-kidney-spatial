#!/usr/bin/env python
"""
whitepaper_groupE_cln.py — Group E (Section 6.2, cLN plasma-myeloid niche) PER-SLIDE,
FACETED PER TISSUE CORE (the CosMx sections are multi-core with mm gaps; full-slide
scatter wastes most of the panel on emptiness).

For each plasma-bearing slide (+ a control) a figure faceted by data-bearing core, each
tightly cropped, showing tissue + myeloid + plasma + delineated plasma aggregates, with
the per-core aggregate count annotated. Plus: a cLN marker-overlay crop (gate-PASSING
transcripts MZB1/CD68/C1QB only), an across-slide POINTS summary, and the per-slide table.

Plasma DBSCAN at the Phase-B scale (eps=0.05 mm, min=10); cores via DBSCAN(eps=0.1 mm).
Analysis is gap-safe (50 µm DBSCAN cannot bridge mm core gaps); faceting is plotting-only.

  conda run -n spatial python py/whitepaper_groupE_cln.py
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN

ROOT = os.environ.get("XENIUM_ROOT", os.getcwd())
OBJ = os.path.join(ROOT, "outputs/objects"); TAB = os.path.join(ROOT, "outputs/tables")
FIG = os.path.join(ROOT, "outputs/figures/whitepaper"); os.makedirs(FIG, exist_ok=True)
EPS, MIN_SAMP = 0.05, 10           # mm — plasma-aggregate scale
CORE_EPS, CORE_MIN = 0.1, 20       # mm — tissue-core detection
PLASMA = ["plasmablast"]; MYELOID = ["macrophage", "myeloid/DC", "monocyte"]
PCOL, MCOL, TCOL = "#ff7f0e", "#2ca02c", "#e8e8e8"
CCOL = {"control": "#4c72b0", "III": "#dd8452", "IV": "#c44e52", "IV+V": "#8172b3"}
# 6 plasma-bearing slides (incl SP19_4061, 2 aggs) + 1 control
SLIDES = [("SP21_213_R1080_S1", "III"), ("SP19_1139_R1080_S3", "IV"),
          ("SP20_642_R1080_S3", "IV+V"), ("SP20_10838_R1080_S1", "III"),
          ("SP19_4061_R1087_S1", "IV"), ("SP18_8471_R1087_S2", "III"),
          ("SMI0016C_SP17SP19", "control")]

a = ad.read_h5ad(os.path.join(OBJ, "cln_cosmx.h5ad"))
xy_all = np.asarray(a.obsm["spatial"], float)
lab_all = a.obs["phase_b_label"].astype(str).values
slide_all = a.obs["sample"].astype(str).values
def gene_vec(g):
    j = a.var_names.get_loc(g)
    return np.asarray(a.X[:, j].todense()).ravel() if hasattr(a.X, "todense") else np.asarray(a.X[:, j]).ravel()
print(f"loaded {a.n_obs} cells")

def plasma_aggs(xy, isP):
    Pidx = np.where(isP)[0]
    agg_id = np.full(len(isP), -1, int); cents = []
    if Pidx.size >= MIN_SAMP:
        cl = DBSCAN(eps=EPS, min_samples=MIN_SAMP).fit(xy[Pidx]).labels_
        for c in [c for c in np.unique(cl) if c != -1]:
            members = Pidx[cl == c]; agg_id[members] = c; cents.append(xy[members].mean(0))
    return agg_id, np.array(cents) if cents else np.empty((0, 2))

# ============================================================================
# per-slide, faceted by tissue core
# ============================================================================
for slide, cls in SLIDES:
    sm = slide_all == slide; xy = xy_all[sm]; lab = lab_all[sm]
    isP = np.isin(lab, PLASMA); isM = np.isin(lab, MYELOID)
    agg_id, cents = plasma_aggs(xy, isP)
    cores = DBSCAN(eps=CORE_EPS, min_samples=CORE_MIN).fit(xy).labels_
    # per-core stats
    stats = []
    for cc in [c for c in np.unique(cores) if c != -1]:
        cm = cores == cc
        nP = int((isP & cm).sum()); nM = int((isM & cm).sum())
        naggs = 0
        if len(cents):
            bb = xy[cm]; xmn, ymn = bb.min(0); xmx, ymx = bb.max(0)
            naggs = int(((cents[:, 0] >= xmn) & (cents[:, 0] <= xmx) &
                         (cents[:, 1] >= ymn) & (cents[:, 1] <= ymx)).sum())
        stats.append(dict(core=cc, mask=cm, nP=nP, nM=nM, n=int(cm.sum()), naggs=naggs))
    # choose cores: plasma-bearing first (by plasma count); if none, top by myeloid
    pcores = sorted([s for s in stats if s["nP"] > 0], key=lambda s: s["nP"], reverse=True)
    sel = pcores[:6] if pcores else sorted(stats, key=lambda s: s["nM"], reverse=True)[:4]
    ncol = 3 if len(sel) > 1 else 1; nrow = int(np.ceil(len(sel) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.5 * ncol, 4.5 * nrow), squeeze=False)
    for k, st in enumerate(sel):
        ax = axes[k // ncol][k % ncol]; cm = st["mask"]; pad = 0.1
        xmn, xmx = xy[cm, 0].min() - pad, xy[cm, 0].max() + pad
        ymn, ymx = xy[cm, 1].min() - pad, xy[cm, 1].max() + pad
        ax.scatter(xy[cm, 0], xy[cm, 1], s=2, c=TCOL, linewidths=0, rasterized=True)
        mm = cm & isM; pp = cm & isP
        ax.scatter(xy[mm, 0], xy[mm, 1], s=6, c=MCOL, linewidths=0)
        ax.scatter(xy[pp, 0], xy[pp, 1], s=9, c=PCOL, linewidths=0)
        if len(cents):
            for cen in cents:
                if xmn <= cen[0] <= xmx and ymn <= cen[1] <= ymx:
                    ax.add_patch(plt.Circle(cen, EPS, fill=False, ec="#222", lw=1.0, alpha=0.85))
        ax.set_xlim(xmn, xmx); ax.set_ylim(ymn, ymx); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"core: {st['nP']} plasma · {st['nM']} myeloid · {st['naggs']} aggs", fontsize=9)
    for k in range(len(sel), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    handles = [plt.Line2D([], [], marker="o", ls="", mfc=PCOL, mec="none", ms=8, label="plasma"),
               plt.Line2D([], [], marker="o", ls="", mfc=MCOL, mec="none", ms=8, label="myeloid"),
               plt.Line2D([], [], marker="o", ls="", mfc="none", mec="#222", ms=10, label="plasma aggregate")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=9, bbox_to_anchor=(0.5, -0.01))
    tot_aggs = len(cents); tot_p = int(isP.sum()); tot_m = int(isM.sum())
    note = ("myeloid present, plasma essentially absent, no aggregates" if (cls == "control" or tot_aggs == 0)
            else "plasma and myeloid intermix within disease cores")
    fig.suptitle(f"cLN {slide} [{cls}] — {tot_p} plasma, {tot_m} myeloid, {tot_aggs} aggregates "
                 f"(faceted by tissue core; showing {len(sel)} of {len(stats)} cores). {note}.", fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95), h_pad=2.4, w_pad=1.2)
    fig.savefig(os.path.join(FIG, f"qcE_slide_{slide}.png"), dpi=125); plt.close(fig)
    print(f"wrote qcE_slide_{slide}.png ({len(sel)}/{len(stats)} cores; P={tot_p} aggs={tot_aggs})")

# ============================================================================
# cLN marker-overlay crop: gate-PASSING transcripts (MZB1 plasma; CD68/C1QB myeloid)
# ============================================================================
slide = "SP21_213_R1080_S1"; sm = slide_all == slide; xy = xy_all[sm]; lab = lab_all[sm]
isP = np.isin(lab, PLASMA); agg_id, cents = plasma_aggs(xy, isP)
# crop around the densest plasma aggregate
sizes = [(np.sum(agg_id == c), c) for c in np.unique(agg_id) if c != -1]
big_c = max(sizes)[1]; cen = xy[agg_id == big_c].mean(0); W = 0.18
mm = (np.abs(xy[:, 0] - cen[0]) < W) & (np.abs(xy[:, 1] - cen[1]) < W)
sx, sl = xy[mm], lab[mm]
gidx = np.where(sm)[0][mm]   # global indices for gene lookup
fig, ax = plt.subplots(1, 4, figsize=(16, 4.4))
ax[0].scatter(sx[:, 0], sx[:, 1], s=12, c=TCOL, linewidths=0)
ax[0].scatter(sx[np.isin(sl, MYELOID), 0], sx[np.isin(sl, MYELOID), 1], s=34, c=MCOL, linewidths=0.2, edgecolor="k", label="myeloid")
ax[0].scatter(sx[np.isin(sl, PLASMA), 0], sx[np.isin(sl, PLASMA), 1], s=34, c=PCOL, linewidths=0.2, edgecolor="k", label="plasma")
ax[0].set_title("cell type", fontsize=11); ax[0].legend(fontsize=8, loc="upper right")
for k, (g, who) in enumerate([("MZB1", "plasma"), ("CD68", "myeloid"), ("C1QB", "myeloid")], start=1):
    gv = gene_vec(g)[gidx]; pos = gv > 0
    ax[k].scatter(sx[~pos, 0], sx[~pos, 1], s=8, c="#e3e3e3", linewidths=0)
    scat = ax[k].scatter(sx[pos, 0], sx[pos, 1], s=34, c=gv[pos], cmap="viridis", linewidths=0.2, edgecolor="k",
                         vmin=1, vmax=max(2, np.percentile(gv[pos], 95)) if pos.any() else 2)
    ax[k].set_title(f"{g}  ({who}, gate-PASS) — {int(pos.sum())} pos", fontsize=11)
    fig.colorbar(scat, ax=ax[k], fraction=0.046, label="counts")
for x_ in ax: x_.set_aspect("equal"); x_.set_xticks([]); x_.set_yticks([])
fig.suptitle("cLN plasma–myeloid focus (SP21_213) — cell-type calls beside gate-PASSING transcripts only "
             "(MZB1=plasma, CD68/C1QB=myeloid; chemokines excluded — failed the usability gate)", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(os.path.join(FIG, "qcE_cln_marker_crop.png"), dpi=120); plt.close(fig)
print("wrote qcE_cln_marker_crop.png")

# ============================================================================
# across-slide POINTS summary (full slide IDs) + per-slide table
# ============================================================================
ps = pd.read_csv(os.path.join(TAB, "cln_phaseB_per_slide.csv"))
pb = ps[ps.n_plasma >= 10].copy().sort_values("nhood_z_plasma_myeloid", ascending=False)
fig, ax = plt.subplots(1, 2, figsize=(14, 6.0))
for j, (col, ttl, ylab) in enumerate([
        ("nhood_z_plasma_myeloid", "plasma×myeloid nhood z (supporting; NOT disease-specific)", "nhood z"),
        ("myeloid_log2enrich_in_plasma_aggs", "myeloid log2-enrich inside plasma aggregates (disease-only)", "log2 enrichment")]):
    d = pb.dropna(subset=[col])
    for i, (_, r) in enumerate(d.iterrows()):
        ax[j].scatter(i, r[col], s=90, c=CCOL.get(r.condition, "#888"), edgecolor="k", linewidth=0.4, zorder=3)
    ax[j].set_xticks(range(len(d))); ax[j].set_xticklabels(d.slide, rotation=45, ha="right", fontsize=7)
    ax[j].axhline(0, color="#999", lw=0.6, ls="--"); ax[j].set_ylabel(ylab); ax[j].set_title(ttl, fontsize=9.5)
handles = [plt.Line2D([], [], marker="o", ls="", mfc=CCOL[c], mec="k", ms=8, label=c) for c in ["control", "III", "IV", "IV+V"]]
ax[0].legend(handles=handles, fontsize=8, loc="upper right", title="ISN class")
fig.suptitle("cLN across slides — per-slide points. Disease specificity = AGGREGATE FORMATION "
             "(controls form 0 aggregates → right panel disease-only);\nnhood z (left) is positive even in a "
             "control, so it is supporting only. Descriptive (patient-confounded, IV+V n=1).", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(os.path.join(FIG, "qcE_across_slide_summary.png"), dpi=130); plt.close(fig)
print("wrote qcE_across_slide_summary.png")

tcols = ["slide", "condition", "n_plasma", "n_myeloid", "nhood_z_plasma_myeloid",
         "n_plasma_aggregates", "plasma_agg_rate", "myeloid_log2enrich_in_plasma_aggs"]
tt = ps[tcols].copy().sort_values(["condition", "n_plasma"], ascending=[True, False])
tt.columns = ["slide", "class", "n plasma", "n myeloid", "plasma×myeloid z", "n aggs", "agg rate", "myeloid log2"]
fig, ax = plt.subplots(figsize=(13, 5)); ax.axis("off")
tb = ax.table(cellText=tt.round(3).fillna("—").values, colLabels=tt.columns, cellLoc="center", loc="center")
tb.auto_set_font_size(False); tb.set_fontsize(8); tb.scale(1, 1.3)
for jj in range(len(tt.columns)):
    tb[0, jj].set_facecolor("#4c72b0"); tb[0, jj].set_text_props(color="white", weight="bold")
for ii in range(1, len(tt) + 1):
    for jj in range(len(tt.columns)):
        tb[ii, jj].set_facecolor("#f4f4f4" if ii % 2 else "#ffffff")
ax.set_title("cLN per-slide plasma–myeloid metrics (all 14 slides)", fontsize=11, pad=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcE_per_slide_table.png"), dpi=150); plt.close(fig)
print("wrote qcE_per_slide_table.png")
print("== Group E done ==")
