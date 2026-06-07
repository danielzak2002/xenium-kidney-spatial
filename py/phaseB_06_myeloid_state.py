#!/usr/bin/env python
"""
phaseB_06_myeloid_state.py — EXPLORATORY, hypothesis-generating myeloid state in the
plasma niche. Strictly correlative: a static snapshot cannot establish recruitment or
directionality. We ask which (usable) myeloid-state markers distinguish the plasma-niche
myeloid in cLN, and how that compares to RCC myeloid on the shared genes.

Pipeline:
  1. USABILITY GATE (the BAFF lesson) on CosMx: a candidate is kept only if it is both
     DETECTED in myeloid and ENRICHED in myeloid over the ambient (epithelial/stroma)
     floor. Ambient-flat genes (BAFF-like) are dropped before any biology.
  2. cLN in-vs-out-of-niche: myeloid INSIDE plasma aggregates (within 50um of an
     aggregate plasma cell) vs OTHER myeloid, same slide. Per-slide log2 fold, summarised.
  3. Cross-context (cautious): SHARED usable genes (CosMx ∩ Xenium), relative/normalised
     myeloid state, cLN plasma-niche myeloid vs RCC myeloid. Different-panel caveat stated.

Reuses persisted parquets (plasma_agg_id on plasma cells) + the h5ads for expression.
  conda run -n spatial python py/phaseB_06_myeloid_state.py
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

ROOT = os.environ.get("XENIUM_ROOT", os.getcwd())
FIG  = os.path.join(ROOT, "outputs/figures"); TAB = os.path.join(ROOT, "outputs/tables")
OBJ  = os.path.join(ROOT, "outputs/objects")

CAND = {
 "myeloid": ["C1QA","C1QB","C1QC","CD68","CD163","MRC1","MERTK","TREM2","APOE","SPP1","LYZ",
             "S100A8","S100A9","IL1B","CXCL9","CXCL10","HLA-DRA","CD74"],
 "ifn":     ["ISG15","MX1","IFI6","STAT1","IFI44L"],
 "chemokine": ["CCL2","CSF1","CXCL12"],
}
ALL_CAND = [g for v in CAND.values() for g in v]
G2GRP = {g: grp for grp, gs in CAND.items() for g in gs}
MYELOID_LBL_PAT = ("monocyte", "macrophage", "myeloid/dc")    # excl mregDC (handled below)
EPI_PAT = ("tubule","epitheli","podocyte","intercalated","principal","urotheli","henle",
           "vasa.recta","capillary","glomerul","endotheli","fibroblast","myofibroblast",
           "pelvic","stroma","connecting","progenitor")
R_MM = 0.05; R_UM = 50.0
DET_MIN, ENRICH_MIN = 0.05, 1.0     # gate: >=5% myeloid positive AND >=2x over ambient

def is_myeloid(labels):
    ll = pd.Series(labels).str.lower()
    return (ll.str.contains("|".join(MYELOID_LBL_PAT)) & ~ll.str.contains("mregdc")).values
def is_epi(labels):
    return pd.Series(labels).str.lower().str.contains("|".join(EPI_PAT)).values

def col_mean(a, gene, mask, layer="lognorm"):
    j = a.var_names.get_loc(gene)
    X = a.layers[layer] if layer in a.layers else a.X
    v = np.asarray(X[mask, j].todense()).ravel() if hasattr(X, "todense") else np.asarray(X[mask, j]).ravel()
    return v
def col_counts(a, gene, mask):
    j = a.var_names.get_loc(gene)
    v = np.asarray(a.X[mask, j].todense()).ravel() if hasattr(a.X, "todense") else np.asarray(a.X[mask, j]).ravel()
    return v

# ============================================================================
# 1. USABILITY GATE (cLN CosMx)
# ============================================================================
print("== loading cLN CosMx ==")
cln = ad.read_h5ad(os.path.join(OBJ, "cln_cosmx.h5ad"))
lab = cln.obs["phase_b_label"].astype(str).values
m_my, m_epi = is_myeloid(lab), is_epi(lab)
print(f"  myeloid={m_my.sum()}  ambient(epi/stroma)={m_epi.sum()}")
present = [g for g in ALL_CAND if g in cln.var_names]
rows = []
for g in present:
    cnt_my = col_counts(cln, g, m_my); det = float((cnt_my > 0).mean())
    mean_my = float(col_mean(cln, g, m_my).mean())
    mean_epi = float(col_mean(cln, g, m_epi).mean())
    enrich = float(np.log2((mean_my + 1e-6) / (mean_epi + 1e-6)))
    usable = (det >= DET_MIN) and (enrich >= ENRICH_MIN)
    rows.append(dict(gene=g, group=G2GRP[g], det_myeloid=round(det, 3),
                     mean_myeloid=round(mean_my, 3), mean_ambient=round(mean_epi, 3),
                     log2_enrich_vs_ambient=round(enrich, 2), usable=usable))
gate = pd.DataFrame(rows).sort_values(["usable", "log2_enrich_vs_ambient"], ascending=[False, False])
gate.to_csv(f"{TAB}/myeloid_usability_gate.csv", index=False)
usable_genes = gate.loc[gate.usable, "gene"].tolist()
print("\n== USABILITY GATE (cLN CosMx) ==")
print(gate.to_string(index=False))
print(f"\n  USABLE ({len(usable_genes)}): {usable_genes}")
print(f"  DROPPED: {gate.loc[~gate.usable,'gene'].tolist()}  (ambient-flat / not myeloid-enriched)")

# ============================================================================
# 2. cLN: myeloid INSIDE vs OUTSIDE plasma aggregates (per slide -> summary)
# ============================================================================
pc = pd.read_parquet(os.path.join(OBJ, "cln_plasma_myeloid_cells.parquet"))
# index expression by cell_id (align parquet <-> h5ad)
cln_idx = {c: i for i, c in enumerate(cln.obs_names)}
pc["h5_row"] = pc["cell_id"].map(cln_idx)
pc = pc[pc.h5_row.notna()].copy(); pc["h5_row"] = pc["h5_row"].astype(int)
print(f"\n== cLN in/out niche (usable genes; per slide) ==")
inout_rows = []
for s, sd in pc.groupby("slide"):
    agg_pl = sd[(sd.niche_label == "Plasma") & (sd.plasma_agg_id >= 0)]
    my = sd[sd.niche_label == "Myeloid"]
    if len(agg_pl) < 10 or len(my) < 30:
        continue
    tree = cKDTree(agg_pl[["x_mm", "y_mm"]].values)
    d, _ = tree.query(my[["x_mm", "y_mm"]].values)
    inside = d <= R_MM
    if inside.sum() < 10 or (~inside).sum() < 10:
        continue
    rin = my.h5_row.values[inside]; rout = my.h5_row.values[~inside]
    for g in usable_genes:
        vi = col_mean(cln, g, rin).mean(); vo = col_mean(cln, g, rout).mean()
        inout_rows.append(dict(slide=s, condition=sd.condition.iloc[0], gene=g,
                               n_in=int(inside.sum()), n_out=int((~inside).sum()),
                               mean_in=round(float(vi), 3), mean_out=round(float(vo), 3),
                               log2_in_vs_out=round(float(np.log2((vi + 1e-6) / (vo + 1e-6))), 3)))
inout = pd.DataFrame(inout_rows)
inout.to_csv(f"{TAB}/myeloid_cln_in_vs_out_niche.csv", index=False)
# summarise across slides
summ = (inout.groupby("gene")["log2_in_vs_out"]
        .agg(mean_log2="mean", sd_log2="std", n_slides="count",
             n_pos=lambda x: int((x > 0).sum()))
        .reset_index().sort_values("mean_log2", ascending=False))
print(summ.to_string(index=False))
n_slides_used = inout.slide.nunique()
print(f"  ({n_slides_used} slides with an evaluable in/out split)")

# ============================================================================
# 3. cross-context: SHARED usable genes, cLN niche-myeloid vs RCC myeloid
# ============================================================================
def myeloid_detection(fn, genes):
    a = ad.read_h5ad(os.path.join(OBJ, fn))
    lab = a.obs["phase_b_label"].astype(str).values
    mm = is_myeloid(lab)
    out = {}
    for g in genes:
        if g in a.var_names:
            out[g] = float((col_counts(a, g, mm) > 0).mean())
        else:
            out[g] = np.nan
    return out, int(mm.sum())

rcc_present = set(ad.read_h5ad(os.path.join(OBJ, "kidney_RCC_protein.h5ad"), backed="r").var_names)
shared = [g for g in usable_genes if g in rcc_present]
print(f"\n== cross-context (shared usable genes: {shared}) ==")
# cLN niche-myeloid detection (inside) vs cLN non-niche vs RCC/preview myeloid
cln_in_rows = pc[pc.niche_label == "Myeloid"]
# recompute inside flag pooled across slides for detection
inside_ids = []
for s, sd in pc.groupby("slide"):
    agg_pl = sd[(sd.niche_label == "Plasma") & (sd.plasma_agg_id >= 0)]
    my = sd[sd.niche_label == "Myeloid"]
    if len(agg_pl) < 10 or len(my) == 0:
        continue
    d, _ = cKDTree(agg_pl[["x_mm", "y_mm"]].values).query(my[["x_mm", "y_mm"]].values)
    inside_ids += my.h5_row.values[d <= R_MM].tolist()
inside_ids = np.array(sorted(set(inside_ids)), dtype=int)
det_rows = []
for g in shared:
    d_in  = float((col_counts(cln, g, inside_ids) > 0).mean()) if len(inside_ids) else np.nan
    det_rows.append(dict(gene=g, cLN_niche_myeloid=round(d_in, 3)))
ddf = pd.DataFrame(det_rows).set_index("gene")
rcc_det, n_rcc = myeloid_detection("kidney_RCC_protein.h5ad", shared)
prev_det, n_prev = myeloid_detection("kidney_preview_PRCC.h5ad", shared)
ddf["RCC_myeloid"] = pd.Series(rcc_det); ddf["preview_myeloid"] = pd.Series(prev_det)
ddf = ddf.reset_index()
ddf.to_csv(f"{TAB}/myeloid_crosscontext_shared.csv", index=False)
print(f"  detection fraction (frac>0) — cLN niche-myeloid (n={len(inside_ids)}) vs RCC myeloid "
      f"(n={n_rcc}) vs preview (n={n_prev}):")
print(ddf.to_string(index=False))
print("  CAVEAT: different panels/platforms (CosMx 957 vs Xenium 405/377) -> detection")
print("  fractions are NOT absolute-comparable; read relative/within-marker patterns only.")

# ============================================================================
# FIGURE
# ============================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5.6))
# (A) cLN in-vs-out niche: per-gene log2 fold, per-slide points + mean
order = summ.gene.tolist()
for i, g in enumerate(order):
    v = inout.loc[inout.gene == g, "log2_in_vs_out"].values
    ax1.scatter(v, np.full(len(v), i) + np.random.default_rng(i).normal(0, 0.06, len(v)),
                s=28, c="#888", alpha=0.7, zorder=2)
    ax1.scatter(np.mean(v), i, s=130, c="#c44e52", marker="D", zorder=3, edgecolor="k", linewidth=0.4)
ax1.axvline(0, color="#999", lw=0.8, ls="--")
ax1.set_yticks(range(len(order))); ax1.set_yticklabels(order)
ax1.set_xlabel("log2(myeloid INSIDE / OUTSIDE plasma niche)"); ax1.invert_yaxis()
ax1.set_title(f"cLN: plasma-niche vs other myeloid\n({n_slides_used} slides; diamond=mean, dots=per-slide)", fontsize=10)
# (B) cross-context detection of shared genes
gx = np.arange(len(shared)); w = 0.27
ax2.bar(gx - w, ddf["cLN_niche_myeloid"], w, label="cLN niche-myeloid", color="#dd8452")
ax2.bar(gx,     ddf["RCC_myeloid"],       w, label="RCC myeloid",       color="#4c72b0")
ax2.bar(gx + w, ddf["preview_myeloid"],   w, label="preview myeloid",   color="#8172b3")
ax2.set_xticks(gx); ax2.set_xticklabels(ddf.gene, rotation=30, ha="right")
ax2.set_ylabel("detection fraction (frac > 0)")
ax2.set_title("Shared myeloid markers (panel caveat: not absolute-comparable)", fontsize=10)
ax2.legend(fontsize=8)
fig.suptitle("Exploratory myeloid state in the plasma niche (correlative; hypothesis-generating)", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/myeloid_state_plasma_niche.png", dpi=160); plt.close(fig)

# ---- hypothesis statement (printed; correlative only) ----------------------
top_up = summ.loc[summ.mean_log2 > 0.2].sort_values("mean_log2", ascending=False)
top_up = top_up[top_up.n_pos >= np.ceil(top_up.n_slides / 2)]
markers = top_up.gene.tolist()
print("\n== HYPOTHESIS (correlative; NOT causal/directional) ==")
if markers:
    print(f"  Plasma-niche myeloid in cLN are enriched (inside>outside, consistent across slides) for: "
          f"{markers}. This RAISES THE HYPOTHESIS that myeloid in this state may associate with /")
    print( "  support the plasma niche. A static snapshot cannot establish recruitment or direction.")
else:
    print("  No usable marker consistently distinguishes niche myeloid (inside~outside) -> no state")
    print("  signal detectable on this panel; hypothesis not supported by available markers.")
print("\nwrote: myeloid_usability_gate.csv, myeloid_cln_in_vs_out_niche.csv, "
      "myeloid_crosscontext_shared.csv, myeloid_state_plasma_niche.png")
print("== phaseB_06 done ==")
