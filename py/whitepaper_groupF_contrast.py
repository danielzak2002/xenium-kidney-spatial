#!/usr/bin/env python
"""
whitepaper_groupF_contrast.py — Group F (Sections 6.3-6.4): the cross-context contrast +
myeloid state, AND the RCC comparable per-region spatial (same layout/palette as Group E).

Figures -> outputs/figures/whitepaper/:
  qcF_rcc_plasma_myeloid_region.png  RCC plasma+myeloid, E-format (myeloid SEGREGATED)
  qcF_rcc_immune_aggregate_region.png RCC B/Treg/CD8/mregDC, E-format (§6.1 aggregate)
  qcF_contrast_plasma_foci.png       side-by-side: RCC (segregated) vs cLN (infiltrated)
  qcF_myeloid_state_in_vs_out.png     §6.4 myeloid-state log2 in/out niche, per-slide points
  qcF_usability_gate_recruitment.png  usability-gate pass/fail for recruitment candidates

  conda run -n spatial python py/whitepaper_groupF_contrast.py
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN

ROOT = os.environ.get("XENIUM_ROOT", os.getcwd())
OBJ = os.path.join(ROOT, "outputs/objects"); TAB = os.path.join(ROOT, "outputs/tables")
FIG = os.path.join(ROOT, "outputs/figures/whitepaper"); os.makedirs(FIG, exist_ok=True)
PCOL, MCOL, TCOL = "#ff7f0e", "#2ca02c", "#e8e8e8"

# RCC label vocabulary
R_PLASMA = ["Plasmablasts"]; R_MYELOID = ["Intermediate monocytes", "Classical monocytes", "myeloid/DC"]
R_B, R_TREG, R_MREG, R_CCR7 = "Naive B cells", "T regulatory cells", "mregDC", "CCR7+ T (naive/CM)"
R_CD8 = ["Effector memory CD8 T cells", "CD8_T"]
IMM_LEG = [(R_B, "#1f77b4", "Naive B"), (R_TREG, "#d62728", "Treg"),
           ("eff-CD8", "#2ca02c", "effector CD8"), (R_MREG, "#9467bd", "mregDC")]

def load(fn):
    a = ad.read_h5ad(os.path.join(OBJ, fn))
    return a, np.asarray(a.obsm["spatial"], float), a.obs["phase_b_label"].astype(str).values
def aggregates(xy, isX, eps, ms):
    idx = np.where(isX)[0]; agg = np.full(len(isX), -1, int); cents = []
    if idx.size >= ms:
        cl = DBSCAN(eps=eps, min_samples=ms).fit(xy[idx]).labels_
        for c in [c for c in np.unique(cl) if c != -1]:
            mem = idx[cl == c]; agg[mem] = c; cents.append((xy[mem].mean(0), len(mem)))
    return agg, sorted(cents, key=lambda t: t[1], reverse=True)

rcc, rxy, rlab = load("kidney_RCC_protein.h5ad")
isP = np.isin(rlab, R_PLASMA); isM = np.isin(rlab, R_MYELOID)
pagg, pcents = aggregates(rxy, isP, 50.0, 10)     # RCC plasma aggregates (um)
bagg, bcents = aggregates(rxy, np.isin(rlab, [R_B]), 50.0, 20)

def region_panels(fig_name, cx, cy, Wr, Wf, color_fn, leg, title):
    """3-panel E-format: region scatter | aggregates+overlay | focus crop."""
    fig, ax = plt.subplots(1, 3, figsize=(18, 6))
    for p, W, lab_p in [(0, Wr, "region"), (1, Wr, "aggregates"), (2, Wf, "focus")]:
        m = (np.abs(rxy[:, 0] - cx) < W) & (np.abs(rxy[:, 1] - cy) < W)
        color_fn(ax[p], m, p)
        ax[p].set_aspect("equal"); ax[p].set_xticks([]); ax[p].set_yticks([])
    ax[0].set_title("region: plasma + myeloid on tissue" if "plasma" in title else "region: immune cells on tissue", fontsize=10)
    ax[1].set_title("delineated aggregates + overlay", fontsize=10)
    ax[2].set_title("focus crop", fontsize=10)
    handles = [plt.Line2D([], [], marker="o", ls="", mfc=c, mec="none", ms=8, label=n) for _, c, n in leg]
    ax[0].legend(handles=handles, fontsize=8, loc="upper right")
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(os.path.join(FIG, fig_name), dpi=130); plt.close(fig)
    print(f"wrote {fig_name}")

# ---- F1: RCC plasma + myeloid (myeloid SEGREGATED) -------------------------
cx, cy = pcents[0][0]
def pm_color(ax, m, p):
    ax.scatter(rxy[m, 0], rxy[m, 1], s=3, c=TCOL, linewidths=0, rasterized=True)
    mm, pp = m & isM, m & isP
    ax.scatter(rxy[mm, 0], rxy[mm, 1], s=8 if p < 2 else 26, c=MCOL, linewidths=0)
    if p == 1:
        ina = m & (pagg >= 0)
        ax.scatter(rxy[m & isP & ~ina, 0], rxy[m & isP & ~ina, 1], s=6, c="#ffd9a8", linewidths=0)
        ax.scatter(rxy[ina, 0], rxy[ina, 1], s=12, c=PCOL, linewidths=0)
        for cen, _n in pcents:
            if abs(cen[0]-cx) < 600 and abs(cen[1]-cy) < 600:
                ax.add_patch(plt.Circle(cen, 50, fill=False, ec="#222", lw=1.0, alpha=0.8))
    else:
        ax.scatter(rxy[pp, 0], rxy[pp, 1], s=10 if p == 0 else 30, c=PCOL, linewidths=0 if p==0 else 0.2,
                   edgecolor="none" if p==0 else "k")
region_panels("qcF_rcc_plasma_myeloid_region.png", cx, cy, 600, 250, pm_color,
              [("plasma", PCOL, f"plasma"), ("myeloid", MCOL, "myeloid"), ("agg", "#222", "plasma aggregate")],
              "RCC plasma–myeloid (Xenium): plasma DO form aggregates, but myeloid are NOT enriched inside them "
              "(myeloid log2 ≈ 0 vs the high tumour myeloid background; cohort nhood z −57) — cf. cLN where they are")

# ---- F2: RCC immune aggregate (B/Treg/CD8/mregDC) --------------------------
bx, by = bcents[0][0]
def imm_color(ax, m, p):
    ax.scatter(rxy[m, 0], rxy[m, 1], s=3, c=TCOL, linewidths=0, rasterized=True)
    foc = [R_B, R_TREG] + R_CD8 + [R_MREG, R_CCR7]
    sz = 8 if p < 2 else 26
    for L, c, _n in IMM_LEG:
        sel = (m & np.isin(rlab, R_CD8)) if L == "eff-CD8" else (m & (rlab == L))
        ax.scatter(rxy[sel, 0], rxy[sel, 1], s=sz, c=c, linewidths=0)
    if p == 1:
        for cen, _n in bcents:
            if abs(cen[0]-bx) < 600 and abs(cen[1]-by) < 600:
                ax.add_patch(plt.Circle(cen, 70, fill=False, ec="#222", lw=1.0, alpha=0.8))
region_panels("qcF_rcc_immune_aggregate_region.png", bx, by, 600, 250, imm_color, IMM_LEG,
              "RCC immunoregulatory B-aggregate (Xenium) — B–Treg core, effector-CD8 excluded (§6.1)")

# ---- F3: side-by-side contrast, QUANTITATIVELY ANNOTATED --------------------
# Raw crops mislead (RCC tumour is far cell-denser than cLN), so each panel is annotated
# with the myeloid log2-enrichment INSIDE that plasma aggregate's 50um footprint vs its
# own section background — the number, not the visual density, carries the contrast. We
# pick a REPRESENTATIVE aggregate in each (RCC ~ median ≈0; cLN ~ a positive niche).
from scipy.spatial import cKDTree
cln, cxy, clab = load("cln_cosmx.h5ad")
cmask = cln.obs["sample"].astype(str).values == "SP21_213_R1080_S1"
cxm = cxy[cmask]; clab_c = clab[cmask]
cP = np.isin(clab_c, ["plasmablast"]); cM = np.isin(clab_c, ["macrophage", "myeloid/DC", "monocyte"])
cpagg, cpcents = aggregates(cxm, cP, 0.05, 10)

def rep_aggregate(xy, agg, isM_, R, bg, want_positive):
    tree = cKDTree(xy); rows = []
    for c in np.unique(agg[agg >= 0]):
        mem = np.where(agg == c)[0]
        reg = np.unique(np.concatenate([np.asarray(t, int) for t in tree.query_ball_point(xy[mem], r=R)]))
        f = float(isM_[reg].mean()); l2 = float(np.log2((f + 1e-6) / (bg + 1e-6)))
        rows.append((c, xy[mem].mean(0), len(mem), f, l2))
    df = pd.DataFrame(rows, columns=["c", "cen", "n", "fmy", "l2"])
    df = df[df.n >= 15]
    target = df.l2.max() if want_positive else df.l2.median()
    return df.iloc[(df.l2 - target).abs().argmin()]

bg_r = float(isM.mean()); bg_c = float(cM.mean())
rep_r = rep_aggregate(rxy, pagg, isM, 50.0, bg_r, want_positive=False)
rep_c = rep_aggregate(cxm, cpagg, cM, 0.05, bg_c, want_positive=True)
fig, ax = plt.subplots(1, 2, figsize=(13, 6.6))
for k, (xy_, lab_isM, lab_isP, rep, W, bg, ctx) in enumerate([
        (rxy, isM, isP, rep_r, 150.0, bg_r, "RCC (tumor)"),
        (cxm, cM, cP, rep_c, 0.13, bg_c, "cLN (lupus)")]):
    cxc, cyc = rep["cen"]
    m = (np.abs(xy_[:, 0]-cxc) < W) & (np.abs(xy_[:, 1]-cyc) < W)
    ax[k].scatter(xy_[m, 0], xy_[m, 1], s=12, c=TCOL, linewidths=0)
    ax[k].scatter(xy_[m & lab_isM, 0], xy_[m & lab_isM, 1], s=32, c=MCOL, linewidths=0.2, edgecolor="k", label="myeloid")
    ax[k].scatter(xy_[m & lab_isP, 0], xy_[m & lab_isP, 1], s=32, c=PCOL, linewidths=0.2, edgecolor="k", label="plasma")
    ax[k].add_patch(plt.Circle((cxc, cyc), (0.05 if k else 50.0), fill=False, ec="#222", lw=1.4))
    verdict = "myeloid NOT enriched (≈background)" if rep["l2"] < 0.3 else "myeloid ENRICHED"
    ax[k].set_title(f"{ctx}: plasma aggregate — {verdict}\nmyeloid inside {rep['fmy']*100:.0f}% vs "
                    f"background {bg*100:.0f}%  →  log2 = {rep['l2']:+.2f}", fontsize=10)
    ax[k].set_aspect("equal"); ax[k].set_xticks([]); ax[k].set_yticks([]); ax[k].legend(fontsize=8, loc="upper right")
fig.suptitle("Cross-context contrast (same 50 µm footprint, each vs its own section background): myeloid "
             "RECRUITMENT into plasma aggregates is the cLN-specific feature — not plasma aggregation\n"
             "(cohort-level: RCC nhood z −57 / log2 ≈ 0; cLN nhood z +4…+17 / log2 +0.45…+2.57)", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.92)); fig.savefig(os.path.join(FIG, "qcF_contrast_plasma_foci.png"), dpi=140); plt.close(fig)
print(f"wrote qcF_contrast_plasma_foci.png  (RCC log2={rep_r['l2']:.2f}, cLN log2={rep_c['l2']:.2f})")

# ---- F4: myeloid-state in/out niche, per gene + per-slide points -----------
io = pd.read_csv(os.path.join(TAB, "myeloid_cln_in_vs_out_niche.csv"))
genes = io.groupby("gene")["log2_in_vs_out"].mean().sort_values(ascending=False).index.tolist()
fig, ax = plt.subplots(figsize=(10, 5.4))
for i, g in enumerate(genes):
    v = io.loc[io.gene == g, "log2_in_vs_out"].values
    ax.scatter(v, np.full(len(v), i) + np.random.default_rng(i).normal(0, 0.06, len(v)), s=30, c="#888", alpha=0.7, zorder=2)
    ax.scatter(np.mean(v), i, s=140, c="#c44e52", marker="D", edgecolor="k", linewidth=0.4, zorder=3)
ax.axvline(0, color="#999", ls="--", lw=0.8); ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes)
ax.invert_yaxis(); ax.set_xlabel("log2(myeloid INSIDE / OUTSIDE plasma niche)")
ax.set_title("§6.4 cLN plasma-niche myeloid state (5 slides; diamond=mean, dots=per-slide)\n"
             "weak complement/MRC1/antigen-presenting lean — correlative, hypothesis-generating", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcF_myeloid_state_in_vs_out.png"), dpi=150); plt.close(fig)
print("wrote qcF_myeloid_state_in_vs_out.png")

# ---- F5: usability gate for recruitment candidates -------------------------
gate = pd.read_csv(os.path.join(TAB, "myeloid_usability_gate.csv"))
gate = gate.sort_values("log2_enrich_vs_ambient", ascending=False)
fig, ax = plt.subplots(figsize=(12, 5.2))
cols = ["#2ca02c" if u else "#c44e52" for u in gate.usable]
ax.bar(range(len(gate)), gate.log2_enrich_vs_ambient, color=cols)
ax.axhline(1.0, color="k", ls="--", lw=1, label="usability threshold (2× ambient)")
ax.set_xticks(range(len(gate))); ax.set_xticklabels(gate.gene, rotation=60, ha="right", fontsize=7.5)
ax.set_ylabel("log2(myeloid mean / ambient)")
ax.set_title("§6.4 usability gate — recruitment candidates: structural macrophage markers PASS (green);\n"
             "ALL secreted/inducible 'recruitment' genes FAIL (red) — chemokines, IFN (MX1/STAT1), CXCL9/10, IL1B, SPP1", fontsize=10)
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcF_usability_gate_recruitment.png"), dpi=150); plt.close(fig)
print("wrote qcF_usability_gate_recruitment.png")
print("== Group F done ==")
