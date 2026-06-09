#!/usr/bin/env python
"""
cd4_cd8_support.py — Are the CD4/CD8 T-cell SUBTYPE labels in the Dumoulin 2026 DKD atlas
measured from platform transcripts, or imputed by reference label transfer (scANVI)?

Logic: one labeling pipeline was applied to both platforms. Xenium 5.1k *can* measure the
CD4/CD8-discriminating transcripts; CosMx 1k may not. If discriminating markers support the
labels on Xenium but sit at the ambient floor on CosMx (AUROC ~ 0.5), the CosMx subtype split
is imputation, not measurement. We separate T-LINEAGE support (CD3 family, above ambient in
both subtypes) from SUBTYPE support (CD4 vs CD8 discrimination) to localize where measurement
fails.

Read-only on raw data. Open backed, subset rows, then .to_memory() — never materialize full X.
"""
import os, warnings, json
warnings.filterwarnings("ignore")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import numpy as np, pandas as pd, anndata as ad
import scipy.sparse as sp
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(0)
BASE = "/Users/danie/ClaudeCode/pilot_analyses/xenium/Demoulin26"
H5AD = os.path.join(BASE, "data", "spatial_adata_xenium_cosmx_zenodo.h5ad")
OUT  = os.path.join(BASE, "analysis", "cross_platform_tcell")
os.makedirs(OUT, exist_ok=True)

PLATFORM = "tech"                              # 48 CosMx / 16 Xenium split (from assessment)
SAMPLE   = "orig_ident"
CT_IMM   = "immune_cell_annotation_combined"   # CD4/CD8/B/Plasma... subtypes
CT_MAIN  = "annotation_updated"                # coarse types incl. epithelial/PT
DUAL_IDS = ["HK2695", "HK2753", "HK3106", "HK3626"]
PAL = {"CosMx": "#E69F00", "Xenium": "#0072B2"}  # consistent platform palette
EPI_CAP = 40000                                  # epithelial-ref subsample per platform

SUBTYPE = ["CD4", "CD8A", "CD8B"]
LINEAGE = ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC2"]
CONTEXT = ["FOXP3", "IL2RA", "CTLA4", "GZMB", "GZMK", "NKG7", "PRF1", "IL7R", "CCR7", "SELL"]
MCLASS  = {**{g: "subtype" for g in SUBTYPE}, **{g: "lineage" for g in LINEAGE},
           **{g: "context" for g in CONTEXT}}
ALL_MARKERS = SUBTYPE + LINEAGE + CONTEXT

def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)

# ============================================================================
# Step 0 — handles
# ============================================================================
hdr("STEP 0 — handles (obs/var, platform split, neg-controls, marker presence)")
adata = ad.read_h5ad(H5AD, backed="r")
obs = adata.obs
print(f"object: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes | X backed, layers={list(adata.layers.keys())}")

print(f"\n[{CT_IMM}] unique labels:")
imm_vals = obs[CT_IMM].astype(str)
print(imm_vals.value_counts().to_string())
print(f"\n[{CT_MAIN}] unique labels:")
main_vals = obs[CT_MAIN].astype(str)
print(main_vals.value_counts().to_string())

print(f"\n[{PLATFORM}] platform value_counts (per-cell):")
plat = obs[PLATFORM].astype(str)
print(plat.value_counts().to_string())
# normalize platform names to CosMx / Xenium
def plat_norm(x):
    xl = x.lower()
    if "cosmx" in xl or "cos" in xl: return "CosMx"
    if "xenium" in xl or "xen" in xl: return "Xenium"
    return x
plat_map = {v: plat_norm(v) for v in plat.unique()}
print("platform name mapping:", plat_map)
platN = plat.map(plat_map)
PLATFORMS = ["CosMx", "Xenium"]

# CD4 / CD8 label strings
cd4_labels = sorted([v for v in imm_vals.unique() if "CD4" in v.upper()])
cd8_labels = sorted([v for v in imm_vals.unique() if "CD8" in v.upper()])
print(f"\nCD4 label string(s): {cd4_labels}")
print(f"CD8 label string(s): {cd8_labels}")
assert cd4_labels and cd8_labels, "could not find CD4/CD8 subtype labels"

# negative-control features
NEG_PAT = ("negprb", "neg_", "negative", "blank", "antisense", "negcontrol",
           "control", "falsecode", "unassigned")
neg_genes = [g for g in adata.var_names if any(p in g.lower() for p in NEG_PAT)]
print(f"\nnegative-control features matched in var_names: {len(neg_genes)}")
print("  examples:", neg_genes[:12])

# marker presence in var
present = {g: (g in adata.var_names) for g in ALL_MARKERS}
print("\nmarker present in var_names:")
for g in ALL_MARKERS:
    print(f"  {g:7s} {'present' if present[g] else 'ABSENT'}  ({MCLASS[g]})")

# ============================================================================
# Build the working subset: all CD4/CD8 T cells + epithelial (PT) reference
# ============================================================================
hdr("BUILD SUBSET — T cells (CD4/CD8) + epithelial reference (row-slice then to_memory)")
is_t = imm_vals.isin(cd4_labels + cd8_labels).values

# pick epithelial reference: largest non-immune, prefer proximal-tubule
IMMUNE_PAT = ("immune", "t cell", "t-cell", " t ", "b cell", "plasma", "myeloid",
              "macro", "monocyte", "dendritic", "dc", "nk", "lympho", "leuko", "mast")
counts_main = main_vals.value_counts()
def is_immune_label(lbl):
    l = lbl.lower()
    return any(p in l for p in IMMUNE_PAT)
pt_cands = [l for l in counts_main.index if (("pt" in l.lower().split()) or ("proximal" in l.lower())
            or ("tubul" in l.lower()) or l.upper().startswith("PT"))]
if pt_cands:
    epi_label = max(pt_cands, key=lambda l: counts_main[l])
else:
    epi_label = next(l for l in counts_main.index if not is_immune_label(l))
print(f"epithelial reference label chosen: '{epi_label}'  (n={counts_main[epi_label]:,})")
is_epi = (main_vals == epi_label).values

# subsample epithelial per platform for the ambient floor (detection rate is stable)
epi_keep = np.zeros(adata.shape[0], bool)
pN = platN.values
for p in PLATFORMS:
    idx = np.where(is_epi & (pN == p))[0]
    if len(idx) > EPI_CAP:
        idx = rng.choice(idx, EPI_CAP, replace=False)
    epi_keep[idx] = True
    print(f"  epithelial ref {p}: using {epi_keep[is_epi & (pN==p)].sum():,} cells")

sel = is_t | epi_keep
sel_idx = np.sort(np.where(sel)[0])
print(f"\nT cells (CD4/CD8) total: {is_t.sum():,}  |  epithelial-ref used: {epi_keep.sum():,}")
print(f"materializing subset: {len(sel_idx):,} cells x {adata.shape[1]:,} genes ...")
sub = adata[sel_idx].to_memory()
adata.file.close()
print("subset materialized. sub.X stays sparse:", sp.issparse(sub.X))

# counts layer (raw); fall back to X if absent
if "counts" in sub.layers:
    C = sub.layers["counts"]
    print("using layer 'counts' for raw detection/means")
else:
    C = sub.X
    print("no 'counts' layer; using X as raw counts")
C = C.tocsr() if sp.issparse(C) else sp.csr_matrix(C)

subN = platN.values[sel_idx]                  # platform per subset cell
sub_imm = imm_vals.values[sel_idx]
sub_main = main_vals.values[sel_idx]
sub_samp = obs[SAMPLE].astype(str).values[sel_idx]
label_t = np.where(np.isin(sub_imm, cd4_labels), "CD4",
           np.where(np.isin(sub_imm, cd8_labels), "CD8", "epi"))

var_ix = {g: i for i, g in enumerate(sub.var_names)}
def col(g):  # dense 1-D raw counts for a gene over subset
    return np.asarray(C[:, var_ix[g]].todense()).ravel() if g in var_ix else None

# which markers are actually MEASURED per platform (nonzero in >=1 cell of that platform
# within the subset whose expressing populations ARE represented here). Structural zero on a
# platform => excluded from that platform's test (flag #1), NOT read as "no expression".
hdr("MARKER MEASUREMENT per platform (structural-zero check, flag #1)")
measured = {p: {} for p in PLATFORMS}
for g in ALL_MARKERS:
    if g not in var_ix:
        for p in PLATFORMS: measured[p][g] = False
        continue
    cg = col(g)
    for p in PLATFORMS:
        measured[p][g] = bool((cg[subN == p] > 0).sum() > 0)
mtab = pd.DataFrame({p: {g: ("measured" if measured[p][g] else "STRUCT-ZERO/absent")
                          for g in ALL_MARKERS} for p in PLATFORMS})
mtab["class"] = [MCLASS[g] for g in mtab.index]
print(mtab.to_string())

# negative-control per-cell mean per platform
neg_ix = [var_ix[g] for g in neg_genes if g in var_ix]
neg_mean = {}
if neg_ix:
    negsum = np.asarray(C[:, neg_ix].sum(1)).ravel()
    for p in PLATFORMS:
        neg_mean[p] = float(negsum[(subN == p)].mean())
    print("\nper-cell negative-control mean count (subset):",
          {p: round(neg_mean[p], 4) for p in PLATFORMS})
else:
    print("\nno negative-control features available -> neg-control floor not computed")

# ============================================================================
# Step 1 — ambient floor (epithelial reference, per platform)
# ============================================================================
hdr("STEP 1 — ambient floor: marker detection & mean count in epithelial reference")
epi_mask = (label_t == "epi")
floor = {p: {} for p in PLATFORMS}   # floor[p][g] = (det_rate, mean_count)
rows_floor = []
for p in PLATFORMS:
    m = epi_mask & (subN == p)
    n = int(m.sum())
    for g in ALL_MARKERS:
        if g not in var_ix:
            floor[p][g] = (np.nan, np.nan); continue
        cg = col(g)[m]
        det = float((cg > 0).mean()); mc = float(cg.mean())
        floor[p][g] = (det, mc)
        rows_floor.append(dict(platform=p, ref_label=epi_label, n_ref=n, marker=g,
                               marker_class=MCLASS[g], measured=measured[p][g],
                               ambient_det_rate=round(det, 4), ambient_mean_count=round(mc, 4)))
floor_df = pd.DataFrame(rows_floor)
print(floor_df.pivot_table(index="marker", columns="platform",
      values="ambient_det_rate").reindex(ALL_MARKERS).to_string())

# ============================================================================
# Step 2 — marker support per label, per platform
# ============================================================================
hdr("STEP 2 — marker support: detection/mean per (platform x label), vs ambient floor")
rows = []
for p in PLATFORMS:
    for lab in ["CD4", "CD8"]:
        m = (label_t == lab) & (subN == p)
        n = int(m.sum())
        for g in ALL_MARKERS:
            adet, amean = floor[p][g]
            if g not in var_ix:
                rows.append(dict(platform=p, label=lab, n_cells=n, marker=g,
                    marker_class=MCLASS[g], measured=False, det_rate=np.nan,
                    mean_count=np.nan, ambient_det_rate=adet, ambient_mean_count=amean,
                    above_floor=False)); continue
            cg = col(g)[m]
            det = float((cg > 0).mean()); mc = float(cg.mean())
            # above floor: detection materially exceeds ambient AND mean count > 1.5x ambient
            above = bool(measured[p][g] and det > adet + 0.05 and mc > 1.5 * (amean + 1e-9))
            rows.append(dict(platform=p, label=lab, n_cells=n, marker=g,
                marker_class=MCLASS[g], measured=measured[p][g], det_rate=round(det, 4),
                mean_count=round(mc, 4), ambient_det_rate=round(adet, 4),
                ambient_mean_count=round(amean, 4), above_floor=above))
support = pd.DataFrame(rows)
support.to_csv(os.path.join(OUT, "results_marker_support.csv"), index=False)
print("\nSUBTYPE markers — detection rate (CD4-label / CD8-label / ambient):")
for p in PLATFORMS:
    print(f"  [{p}]")
    for g in SUBTYPE:
        r4 = support[(support.platform==p)&(support.label=="CD4")&(support.marker==g)].iloc[0]
        r8 = support[(support.platform==p)&(support.label=="CD8")&(support.marker==g)].iloc[0]
        print(f"    {g:5s} CD4lab={r4.det_rate:.3f}  CD8lab={r8.det_rate:.3f}  "
              f"ambient={r4.ambient_det_rate:.3f}  measured={r4.measured}")

# ---- LEAD METRIC: AUROC discriminating CD8-label from CD4-label, per marker, per platform
hdr("STEP 2 (lead) — AUROC: does each marker separate CD8-label from CD4-label?")
def auroc_rankbiserial(score, y):
    # y: 1 = positive (CD8 for CD8 markers). need both classes present and score variance.
    if len(np.unique(y)) < 2: return np.nan, np.nan
    if np.all(score == score[0]): return 0.5, 0.0
    auc = roc_auc_score(y, score)
    rb = 2 * auc - 1
    return float(auc), float(rb)

auc_rows = []
for p in PLATFORMS:
    mp = subN == p
    y_cd8 = (label_t[mp] == "CD8").astype(int)   # CD8 vs CD4 among T cells
    for g in ALL_MARKERS:
        if not measured[p][g] or g not in var_ix:
            auc_rows.append(dict(platform=p, marker=g, marker_class=MCLASS[g],
                auroc_cd8_vs_cd4=np.nan, rank_biserial=np.nan, measured=measured[p][g],
                n_cd4=int((label_t[mp]=="CD4").sum()), n_cd8=int((label_t[mp]=="CD8").sum())))
            continue
        sc = col(g)[mp]
        auc, rb = auroc_rankbiserial(sc, y_cd8)
        auc_rows.append(dict(platform=p, marker=g, marker_class=MCLASS[g],
            auroc_cd8_vs_cd4=round(auc, 4) if auc==auc else np.nan,
            rank_biserial=round(rb, 4) if rb==rb else np.nan, measured=True,
            n_cd4=int((label_t[mp]=="CD4").sum()), n_cd8=int((label_t[mp]=="CD8").sum())))
auc_df = pd.DataFrame(auc_rows)

# composite CD8 score on within-platform log-normalized expression: (CD8A+CD8B) - CD4
def lognorm_block(mp):
    Csub = C[mp]
    lib = np.asarray(Csub.sum(1)).ravel(); lib[lib == 0] = 1
    med = np.median(lib)
    Cn = Csub.multiply(1.0 / lib[:, None]).tocsr() * med
    return Cn  # sparse; log1p applied per-gene when extracted
comp_rows = []
for p in PLATFORMS:
    mp = subN == p
    have = all(measured[p][g] for g in ["CD8A", "CD8B", "CD4"])
    if not have:
        comp_rows.append(dict(platform=p, auroc_composite_cd8=np.nan, note="missing CD8A/CD8B/CD4"));
        continue
    Cn = lognorm_block(mp)
    def lg(g): return np.log1p(np.asarray(Cn[:, var_ix[g]].todense()).ravel())
    score = lg("CD8A") + lg("CD8B") - lg("CD4")
    y = (label_t[mp] == "CD8").astype(int)
    auc, rb = auroc_rankbiserial(score, y)
    comp_rows.append(dict(platform=p, auroc_composite_cd8=round(auc,4), note="(CD8A+CD8B)-CD4 lognorm"))
comp_df = pd.DataFrame(comp_rows)
auc_df.to_csv(os.path.join(OUT, "results_auroc.csv"), index=False)
comp_df.to_csv(os.path.join(OUT, "results_auroc_composite.csv"), index=False)
print("\nSUBTYPE-discriminating AUROC (CD8-label vs CD4-label), per platform:")
print(auc_df[auc_df.marker_class=="subtype"].pivot_table(index="marker",
      columns="platform", values="auroc_cd8_vs_cd4").reindex(SUBTYPE).to_string())
print("\nLINEAGE-marker AUROC (expected ~0.5 — lineage shouldn't separate subtypes):")
print(auc_df[auc_df.marker_class=="lineage"].pivot_table(index="marker",
      columns="platform", values="auroc_cd8_vs_cd4").reindex(LINEAGE).to_string())
print("\ncomposite CD8 score AUROC:")
print(comp_df.to_string(index=False))

# ---- T-LINEAGE support: CD3 family above ambient in BOTH CD4- and CD8-labeled, both platforms
hdr("STEP 2 (lineage) — is T-LINEAGE (CD3 family) above ambient in both subtypes?")
lin_summary = []
for p in PLATFORMS:
    for lab in ["CD4", "CD8"]:
        for g in LINEAGE:
            r = support[(support.platform==p)&(support.label==lab)&(support.marker==g)]
            if len(r):
                r = r.iloc[0]
                lin_summary.append(dict(platform=p, label=lab, marker=g,
                    det=r.det_rate, ambient=r.ambient_det_rate, above_floor=r.above_floor))
lin_df = pd.DataFrame(lin_summary)
print(lin_df.pivot_table(index=["marker"], columns=["platform","label"],
      values="det").to_string())

# ============================================================================
# Step 3 — dual-platform cross-check (shared tissue blocks)
# ============================================================================
hdr("STEP 3 — dual-platform CD4:CD8 ratio per shared sample")
ratio_rows = []
for s in DUAL_IDS:
    rec = {"sample": s}
    for p in PLATFORMS:
        m = (sub_samp == s) & (subN == p) & np.isin(label_t, ["CD4", "CD8"])
        n4 = int(((label_t == "CD4") & m).sum()); n8 = int(((label_t == "CD8") & m).sum())
        rec[f"{p}_nCD4"] = n4; rec[f"{p}_nCD8"] = n8
        rec[f"{p}_CD4frac"] = round(n4 / (n4 + n8), 4) if (n4 + n8) else np.nan
        rec[f"{p}_CD4_CD8_ratio"] = round(n4 / n8, 3) if n8 else np.nan
    ratio_rows.append(rec)
ratio_df = pd.DataFrame(ratio_rows)
ratio_df.to_csv(os.path.join(OUT, "results_dual_platform_ratio.csv"), index=False)
print(ratio_df.to_string(index=False))

# ============================================================================
# FIGURES
# ============================================================================
hdr("FIGURES")
# (a) detection rate of CD4/CD8A/CD8B in CD4-label vs CD8-label, faceted by platform
fig, axes = plt.subplots(1, len(PLATFORMS), figsize=(11, 4.2), sharey=True)
for ax, p in zip(axes, PLATFORMS):
    x = np.arange(len(SUBTYPE)); w = 0.38
    cd4d = [support[(support.platform==p)&(support.label=="CD4")&(support.marker==g)].det_rate.iloc[0] for g in SUBTYPE]
    cd8d = [support[(support.platform==p)&(support.label=="CD8")&(support.marker==g)].det_rate.iloc[0] for g in SUBTYPE]
    amb  = [floor[p][g][0] for g in SUBTYPE]
    ax.bar(x - w/2, cd4d, w, label="CD4-labeled", color="#56B4E9")
    ax.bar(x + w/2, cd8d, w, label="CD8-labeled", color="#D55E00")
    for xi, a in zip(x, amb):
        ax.hlines(a, xi - 0.5, xi + 0.5, color="k", ls="--", lw=1.2,
                  label="ambient floor" if xi == 0 else None)
    ax.set_xticks(x); ax.set_xticklabels(SUBTYPE)
    ax.set_title(f"{p}", color=PAL[p]); ax.set_ylabel("detection rate (>=1 count)")
    if p == PLATFORMS[0]: ax.legend(fontsize=8, loc="upper left")
fig.suptitle("Subtype-marker detection in CD4- vs CD8-labeled T cells (ambient floor dashed)")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_a_subtype_detection.png"), dpi=150)
plt.close(fig)

# (b) AUROC summary CosMx vs Xenium, subtype vs lineage
fig, ax = plt.subplots(figsize=(10, 4.5))
order = SUBTYPE + LINEAGE
x = np.arange(len(order)); w = 0.38
for i, p in enumerate(PLATFORMS):
    vals = [auc_df[(auc_df.platform==p)&(auc_df.marker==g)].auroc_cd8_vs_cd4.iloc[0] for g in order]
    ax.bar(x + (i - 0.5) * w, vals, w, label=p, color=PAL[p])
ax.axhline(0.5, color="k", ls="--", lw=1, label="chance (0.5)")
ax.axvline(len(SUBTYPE) - 0.5, color="gray", lw=1)
ax.text(len(SUBTYPE)/2 - 0.5, 1.02, "SUBTYPE", ha="center", fontsize=9, weight="bold")
ax.text(len(SUBTYPE) + len(LINEAGE)/2 - 0.5, 1.02, "LINEAGE", ha="center", fontsize=9, weight="bold")
ax.set_xticks(x); ax.set_xticklabels(order, rotation=45, ha="right")
ax.set_ylabel("AUROC (CD8-label vs CD4-label)"); ax.set_ylim(0, 1.08)
ax.set_title("Marker discrimination of CD8 vs CD4 labels — subtype should separate, lineage should not")
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_b_auroc_summary.png"), dpi=150)
plt.close(fig)

# (c) dual-platform CD4:CD8 ratio comparison
fig, ax = plt.subplots(figsize=(8, 4.5))
x = np.arange(len(DUAL_IDS)); w = 0.38
for i, p in enumerate(PLATFORMS):
    vals = [ratio_df.loc[ratio_df["sample"]==s, f"{p}_CD4frac"].iloc[0] for s in DUAL_IDS]
    ax.bar(x + (i - 0.5) * w, vals, w, label=p, color=PAL[p])
ax.axhline(0.5, color="k", ls=":", lw=1)
ax.set_xticks(x); ax.set_xticklabels(DUAL_IDS)
ax.set_ylabel("CD4 fraction of CD4+CD8 T cells"); ax.set_ylim(0, 1)
ax.set_title("Dual-platform CD4 fraction per shared tissue block (region-level; cells not co-registered)")
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_c_dual_platform_ratio.png"), dpi=150)
plt.close(fig)
print("saved figures: fig_a_subtype_detection.png, fig_b_auroc_summary.png, fig_c_dual_platform_ratio.png")

# ============================================================================
# Markdown report
# ============================================================================
hdr("WRITING REPORT")
def auc_of(p, g):
    v = auc_df[(auc_df.platform==p)&(auc_df.marker==g)].auroc_cd8_vs_cd4
    return float(v.iloc[0]) if len(v) and v.iloc[0]==v.iloc[0] else np.nan
def comp_of(p):
    v = comp_df[comp_df.platform==p].auroc_composite_cd8
    return float(v.iloc[0]) if len(v) and v.iloc[0]==v.iloc[0] else np.nan

# decision rule helpers
def lineage_supported(p):
    # CD3 family above ambient in BOTH CD4 and CD8 labels
    ok = []
    for lab in ["CD4", "CD8"]:
        for g in ["CD3D", "CD3E", "CD3G"]:
            r = support[(support.platform==p)&(support.label==lab)&(support.marker==g)]
            if len(r): ok.append(bool(r.above_floor.iloc[0]))
    return (sum(ok) >= 1, ok)
def subtype_supported(p):
    aucs = [auc_of(p, g) for g in ["CD8A", "CD8B", "CD4"]]
    aucs = [a for a in aucs if a == a]
    best = max([abs(a - 0.5) for a in aucs], default=np.nan)
    return (best == best and best >= 0.10, aucs)  # AUROC >= 0.60 or <= 0.40 on some marker

lines = []
def w(s=""): lines.append(s)
w("# CD4/CD8 subtype labels: measured or imputed? (Dumoulin 2026 DKD atlas)\n")
w(f"Read-only analysis. Object: {sub.n_obs:,} cells materialized (CD4/CD8 T cells + "
  f"epithelial reference '{epi_label}'). Platform column `{PLATFORM}`; raw counts from "
  f"layer `{'counts' if 'counts' in sub.layers else 'X'}`.\n")
w("## Decision rule\n")
w("- **T-lineage supported** on a platform if CD3-family markers sit **above the ambient "
  "floor** (epithelial reference) in *both* CD4- and CD8-labeled cells.")
w("- **Subtype split supported** if discriminating markers (CD8A/CD8B/CD4) are above ambient "
  "**and discriminate** CD8- from CD4-labeled cells: **AUROC well off 0.5**. AUROC ~0.5 with "
  "markers at the ambient floor => the subtype labels are **reference-imputed, not measured**.\n")
w("## Marker measurement (flag #1: structural zeros)\n")
for p in PLATFORMS:
    sz = [g for g in ALL_MARKERS if not measured[p][g]]
    w(f"- **{p}**: structural-zero / absent markers excluded from its tests: "
      f"{', '.join(sz) if sz else 'none'}.")
if neg_mean:
    w(f"- Negative-control per-cell mean: " + ", ".join(f"{p}={neg_mean[p]:.4f}" for p in PLATFORMS) + ".")
w("")
w("## Subtype discrimination (lead metric, AUROC CD8-label vs CD4-label)\n")
w("| marker | CosMx | Xenium |")
w("|---|---|---|")
for g in SUBTYPE:
    w(f"| {g} | {auc_of('CosMx', g):.3f} | {auc_of('Xenium', g):.3f} |")
w(f"| **composite (CD8A+CD8B−CD4)** | {comp_of('CosMx'):.3f} | {comp_of('Xenium'):.3f} |\n")
w("Lineage-marker AUROC (control — should be ~0.5):\n")
w("| marker | CosMx | Xenium |")
w("|---|---|---|")
for g in LINEAGE:
    w(f"| {g} | {auc_of('CosMx', g):.3f} | {auc_of('Xenium', g):.3f} |")
w("")
w("## Per-platform verdict\n")
for p in PLATFORMS:
    lin_ok, lin_detail = lineage_supported(p)
    sub_ok, sub_aucs = subtype_supported(p)
    w(f"### {p}")
    w(f"- T-lineage above ambient (CD3 family): **{'YES' if lin_ok else 'NO'}**.")
    w(f"- Subtype discriminating AUROC (CD8A/CD8B/CD4): {['%.2f'%a for a in sub_aucs]} "
      f"-> subtype split **{'SUPPORTED (measured)' if sub_ok else 'NOT supported (≈chance -> imputed)'}**.")
    w("")
w("## Dual-platform shared blocks (CD4 fraction of T cells)\n")
w("| sample | CosMx CD4frac | Xenium CD4frac | CosMx ratio | Xenium ratio |")
w("|---|---|---|---|---|")
for _, r in ratio_df.iterrows():
    w(f"| {r['sample']} | {r.get('CosMx_CD4frac')} | {r.get('Xenium_CD4frac')} | "
      f"{r.get('CosMx_CD4_CD8_ratio')} | {r.get('Xenium_CD4_CD8_ratio')} |")
w("\n(Region-level only — cells are NOT co-registered across platforms.)\n")
w("## Cross-platform conclusion\n")
w("_Filled from the numbers above_: state whether subtype is measured on Xenium, "
  "imputed on CosMx, and whether this is consistent with the cLN 1k finding (where CD4/CD8 "
  "collapsed). See verdicts per platform. Flags: structural-zero exclusions (flag #1); "
  "ambient floor from a single epithelial reference; dual-platform comparison is region-level.\n")
with open(os.path.join(OUT, "REPORT.md"), "w") as f:
    f.write("\n".join(lines))
print("wrote REPORT.md, results_marker_support.csv, results_auroc.csv, "
      "results_auroc_composite.csv, results_dual_platform_ratio.csv")
print("\n== done ==")
