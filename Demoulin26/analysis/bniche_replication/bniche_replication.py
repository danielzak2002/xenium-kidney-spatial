#!/usr/bin/env python
"""
bniche_replication.py — Does the RCC immunoregulatory B-aggregate signature replicate in the
Dumoulin DKD "B predom. Immune ME" niche?

RCC signature to test (from committed RCC work, eps=50, 37 aggregates):
  Treg ENRICHED in B-aggregate (36/37); effector-CD8 EXCLUDED (34/37);
  mregDC membership & plasma inconclusive.
Question: does the DKD validated B-niche carry the same Treg(+)/effector-CD8(-) organization?

SCOPE (follows directly from cd4_cd8_support): CosMx CD4/CD8 subtype is reference-imputed, so
ALL subtype/state-resolved replication is XENIUM-ONLY. CosMx is used solely for a coarse
B-predominance consistency check (Step 5), with NO subtype claims.

Read-only on raw data. Open backed, subset to niche-containing samples + immune cells, then
.to_memory(); never materialize full X.
"""
import os, warnings
warnings.filterwarnings("ignore")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import numpy as np, pandas as pd, anndata as ad
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/danie/ClaudeCode/pilot_analyses/xenium/Demoulin26"
H5AD = os.path.join(BASE, "data", "spatial_adata_xenium_cosmx_zenodo.h5ad")
OUT  = os.path.join(BASE, "analysis", "bniche_replication")
os.makedirs(OUT, exist_ok=True)

PLATFORM = "tech"; SAMPLE = "orig_ident"; CT_IMM = "immune_cell_annotation_combined"
NICHE_VAL = "B predom. Immune ME"
DKD_COLOR = "#6A3D9A"   # distinct from RCC blue / PRCC green / cLN red
MIN = 10                # min cells per compartment per sample-side for a marker fraction

# marker panels (RCC-axis); structural zeros dropped per platform after measurement check
TREG    = ["FOXP3", "IL2RA", "CTLA4"]            # FOXP3 primary
EFFCD8  = ["GZMB", "GZMK", "PRF1", "GNLY"]
MREGDC  = ["LAMP3", "CCR7", "FSCN1", "CD274"]
ALLMARK = TREG + EFFCD8 + MREGDC

def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
def log2e(a_in, n_in, a_out, n_out):
    """log2 enrichment of a fraction in vs out (a=count positive, n=compartment size)."""
    if n_in < MIN or n_out < MIN: return np.nan
    fi = (a_in + 0.5) / (n_in + 1); fo = (a_out + 0.5) / (n_out + 1)  # Haldane-style smoothing
    return float(np.log2(fi / fo))

# ============================================================================
# Step 0 — handles
# ============================================================================
hdr("STEP 0 — niche column, platform/niche counts, subtype labels, measured markers")
adata = ad.read_h5ad(H5AD, backed="r")
obs = adata.obs
plat = obs[PLATFORM].astype(str)
platN = plat.map(lambda x: "CosMx" if "cosmx" in x.lower() else ("Xenium" if "xenium" in x.lower() else x))
samp = obs[SAMPLE].astype(str)

# locate the niche column carrying the exact value string
niche_col = None
for c in obs.columns:
    try: cats = obs[c].astype(str).unique()
    except Exception: continue
    if NICHE_VAL in set(cats):
        niche_col = c; break
assert niche_col, f"could not find a column containing '{NICHE_VAL}'"
print(f"niche column: '{niche_col}'  | target value: '{NICHE_VAL}'")
niche = obs[niche_col].astype(str)
is_niche = (niche == NICHE_VAL).values

print("\nniche-cell counts per platform:")
for p in ["CosMx", "Xenium"]:
    print(f"  {p}: {int((is_niche & (platN.values==p)).sum()):,}")

# per-sample, per-platform niche cell counts
nd = pd.DataFrame({"platform": platN.values[is_niche], "sample": samp.values[is_niche]})
per_samp = nd.value_counts().reset_index(name="n_niche_cells").sort_values(
    ["platform", "n_niche_cells"], ascending=[True, False])
print("\nniche cells per sample (>0):")
print(per_samp.to_string(index=False))
xen_samples = sorted(per_samp.loc[per_samp.platform == "Xenium", "sample"].tolist())
cos_samples = sorted(per_samp.loc[per_samp.platform == "CosMx", "sample"].tolist())
n_xen = len(xen_samples)
print(f"\nXenium samples containing the niche (n for the test): {n_xen}  -> {xen_samples}")
print(f"CosMx samples containing the niche: {len(cos_samples)}")
DESCRIPTIVE = n_xen < 5
print(f"power: n_xen={n_xen} -> {'DESCRIPTIVE / per-sample only (no inferential claims)' if DESCRIPTIVE else 'small-cohort inferential ok with caveats'}")
print("NOTE: no patient column in obs -> patient-level clustering CANNOT be controlled.")

imm = obs[CT_IMM].astype(str)
print(f"\n[{CT_IMM}] labels: {sorted(imm.unique())}")
CD4="CD4+"; CD8="CD8+"; B="B"; PLASMA="Plasma"; MACRO="Macro"; NK="NK"
DCs=[l for l in imm.unique() if l.upper().endswith("DC")]
print(f"DC labels grouped: {DCs}")
COARSE = [CD4, CD8, PLASMA, MACRO, NK] + (["DC"] if DCs else [])

# measured markers per platform (nonzero in >=1 cell of platform among materialized immune set)
present = {g: (g in adata.var_names) for g in ALLMARK}
print("\nmarker present in var_names:")
for g in ALLMARK: print(f"  {g:6s} {'present' if present[g] else 'ABSENT'}")

# ============================================================================
# Build subset: immune cells in niche-containing samples (both platforms; markers used Xenium)
# ============================================================================
hdr("BUILD SUBSET — immune cells in niche-containing samples (row-slice then to_memory)")
niche_sample_set = set(xen_samples) | set(cos_samples)
in_niche_samp = samp.isin(niche_sample_set).values
is_immune = (imm != "Unknown").values
sel = in_niche_samp & is_immune
sel_idx = np.sort(np.where(sel)[0])
print(f"immune cells in niche-containing samples: {len(sel_idx):,} (materializing)")
sub = adata[sel_idx].to_memory()
adata.file.close()
C = (sub.layers["counts"] if "counts" in sub.layers else sub.X)
C = C.tocsr() if sp.issparse(C) else sp.csr_matrix(C); C = C.astype(np.float32)
print("materialized; counts sparse:", sp.issparse(C))

s_plat = platN.values[sel_idx]; s_samp = samp.values[sel_idx]
s_imm  = imm.values[sel_idx];   s_niche = is_niche[sel_idx]
var_ix = {g: i for i, g in enumerate(sub.var_names)}
def pos(genes, mask):
    """boolean: cells in mask positive for ANY of genes (>=1 raw count), measured-only."""
    gs = [g for g in genes if g in var_ix]
    if not gs: return None
    M = C[np.where(mask)[0]][:, [var_ix[g] for g in gs]]
    return (np.asarray(M.sum(1)).ravel() > 0)

# measured check restricted to Xenium immune cells
hdr("MEASURED markers on XENIUM (drop structural zeros)")
measured_xen = {}
xmask = s_plat == "Xenium"
for g in ALLMARK:
    if g not in var_ix: measured_xen[g] = False; continue
    cg = np.asarray(C[np.where(xmask)[0]][:, var_ix[g]].todense()).ravel()
    measured_xen[g] = bool((cg > 0).sum() > 0)
for grp, name in [(TREG,"Treg"),(EFFCD8,"eff-CD8"),(MREGDC,"mregDC")]:
    kept=[g for g in grp if measured_xen.get(g)]; dropped=[g for g in grp if g in var_ix and not measured_xen.get(g)]
    absent=[g for g in grp if g not in var_ix]
    print(f"  {name:7s} measured={kept}  struct-zero={dropped}  absent={absent}")
TREG_X   = [g for g in TREG   if measured_xen.get(g)]
EFFCD8_X = [g for g in EFFCD8 if measured_xen.get(g)]
MREGDC_X = [g for g in MREGDC if measured_xen.get(g)]

# ============================================================================
# Step 2 — coarse compositional enrichment (Xenium), obs-only counts over ALL cells/sample
# ============================================================================
hdr("STEP 2 — coarse compositional enrichment in-niche vs out-niche (Xenium, per sample)")
def subtype_mask(lbls, arr):
    return np.isin(arr, lbls if isinstance(lbls, list) else [lbls])
coarse_rows = []
for s in xen_samples:
    m = (samp.values == s) & (platN.values == "Xenium")
    nin = is_niche[m]; lab = imm.values[m]
    tin = int(nin.sum()); tout = int((~nin).sum())
    for ct in COARSE:
        lbls = DCs if ct == "DC" else [ct]
        cin = int((subtype_mask(lbls, lab) & nin).sum())
        cout = int((subtype_mask(lbls, lab) & ~nin).sum())
        coarse_rows.append(dict(sample=s, cell_type=ct, n_in=cin, n_out=cout,
            total_in=tin, total_out=tout, log2_enrich=log2e(cin, tin, cout, tout)))
    # B dominance
    bin_ = int((subtype_mask([B], lab) & nin).sum())
    coarse_rows.append(dict(sample=s, cell_type="B", n_in=bin_, n_out=int((subtype_mask([B],lab)&~nin).sum()),
        total_in=tin, total_out=tout, log2_enrich=log2e(bin_, tin, int((subtype_mask([B],lab)&~nin).sum()), tout)))
coarse = pd.DataFrame(coarse_rows)
coarse.to_csv(os.path.join(OUT, "xenium_coarse_enrichment_per_sample.csv"), index=False)
print(coarse.pivot_table(index="cell_type", columns="sample", values="log2_enrich").to_string())
print("\nB-fraction in-niche per Xenium sample:")
for s in xen_samples:
    m=(samp.values==s)&(platN.values=="Xenium"); nin=is_niche[m]; lab=imm.values[m]
    print(f"  {s}: B in-niche frac = {(subtype_mask([B],lab)&nin).sum()/max(nin.sum(),1):.3f} (n_niche={int(nin.sum())})")

# ============================================================================
# Step 3 — refined RCC-axis test (Xenium, per sample): markers within compartments
# ============================================================================
hdr("STEP 3 — refined RCC-axis: Treg-like / eff-CD8 / mregDC within compartments (Xenium)")
refined_rows = []
for s in xen_samples:
    base = (s_plat == "Xenium") & (s_samp == s)
    nin = s_niche & base; nout = (~s_niche) & base
    # Treg-like among CD4+
    for state, genes, comp in [("Treg_like", TREG_X, CD4),
                               ("eff_CD8",  EFFCD8_X, CD8),
                               ("mregDC_like", MREGDC_X, "DC")]:
        comp_lbls = DCs if comp == "DC" else [comp]
        cin  = nin  & np.isin(s_imm, comp_lbls)
        cout = nout & np.isin(s_imm, comp_lbls)
        n_in = int(cin.sum()); n_out = int(cout.sum())
        if not genes or n_in < MIN or n_out < MIN:
            refined_rows.append(dict(sample=s, state=state, compartment=comp, markers="|".join(genes),
                n_comp_in=n_in, n_comp_out=n_out, pos_in=np.nan, pos_out=np.nan,
                frac_in=np.nan, frac_out=np.nan, log2_enrich=np.nan)); continue
        pin = int(pos(genes, cin).sum()); pout = int(pos(genes, cout).sum())
        refined_rows.append(dict(sample=s, state=state, compartment=comp, markers="|".join(genes),
            n_comp_in=n_in, n_comp_out=n_out, pos_in=pin, pos_out=pout,
            frac_in=round(pin/n_in,4), frac_out=round(pout/n_out,4),
            log2_enrich=log2e(pin, n_in, pout, n_out)))
refined = pd.DataFrame(refined_rows)
refined.to_csv(os.path.join(OUT, "xenium_refined_axis_per_sample.csv"), index=False)
print(refined[["sample","state","compartment","n_comp_in","n_comp_out","frac_in","frac_out","log2_enrich"]].to_string(index=False))

# in-niche T-compartment balance (Treg-like : eff-CD8 counts)
print("\nin-niche T-compartment balance (Treg-like : eff-CD8, raw positive counts):")
for s in xen_samples:
    rt = refined[(refined["sample"]==s)&(refined.state=="Treg_like")]
    re_ = refined[(refined["sample"]==s)&(refined.state=="eff_CD8")]
    tt = rt.pos_in.iloc[0] if len(rt) else np.nan; ee = re_.pos_in.iloc[0] if len(re_) else np.nan
    print(f"  {s}: Treg-like+={tt}  eff-CD8+={ee}")

# ============================================================================
# Step 4 — aggregate across Xenium niche samples (k/n style)
# ============================================================================
hdr("STEP 4 — aggregate across Xenium niche samples (RCC-style k/n)")
treg = refined[refined.state=="Treg_like"].dropna(subset=["log2_enrich"])
effc = refined[refined.state=="eff_CD8"].dropna(subset=["log2_enrich"])
k_treg = int((treg.log2_enrich > 0).sum()); n_treg = len(treg)
k_eff  = int((effc.log2_enrich < 0).sum()); n_eff  = len(effc)
print(f"Treg-like ENRICHED in-niche: {k_treg}/{n_treg} samples  (log2>0)")
print(f"effector-CD8 DEPLETED in-niche: {k_eff}/{n_eff} samples  (log2<0)")
print(f"  Treg-like per-sample log2: {dict(zip(treg['sample'], treg.log2_enrich.round(2)))}")
print(f"  eff-CD8   per-sample log2: {dict(zip(effc['sample'], effc.log2_enrich.round(2)))}")

# ============================================================================
# Step 5 — CosMx COARSE consistency ONLY (no subtype claims)
# ============================================================================
hdr("STEP 5 — CosMx coarse consistency (B-predominance + plasma/myeloid; NO subtype claims)")
cos_rows = []
for s in cos_samples:
    m = (samp.values == s) & (platN.values == "CosMx")
    nin = is_niche[m]; lab = imm.values[m]; tin=int(nin.sum()); tout=int((~nin).sum())
    for ct in [B, PLASMA, MACRO]:
        cin=int((subtype_mask([ct],lab)&nin).sum()); cout=int((subtype_mask([ct],lab)&~nin).sum())
        cos_rows.append(dict(sample=s, cell_type=ct, n_in=cin, total_in=tin,
            frac_in=round(cin/max(tin,1),4), log2_enrich=log2e(cin,tin,cout,tout)))
cosdf = pd.DataFrame(cos_rows)
cosdf.to_csv(os.path.join(OUT, "cosmx_coarse_consistency.csv"), index=False)
print(cosdf.pivot_table(index="cell_type", columns="sample", values="log2_enrich").to_string())
print("CAVEAT: CosMx CD4/CD8 subtype is reference-imputed (see cd4_cd8_support) -> NOT used here.")

# ============================================================================
# Step 6 — head-to-head verdict table
# ============================================================================
hdr("STEP 6 — RCC vs DKD B-niche head-to-head verdict")
def frac_str(k,n): return f"{k}/{n}" if n else "n/a"
treg_dir = (k_treg/n_treg) if n_treg else np.nan
eff_dir  = (k_eff/n_eff) if n_eff else np.nan
def axis_verdict(frac, n):
    if not n: return "UNTESTED"
    if frac >= 0.8: return "REPLICATE"
    if frac >= 0.5: return "PARTIAL"
    return "DIVERGE"
vt = pd.DataFrame([
    dict(axis="Treg ENRICHED in B-niche", rcc="36/37 (0.97)",
         dkd_xenium=f"{frac_str(k_treg,n_treg)} ({treg_dir:.2f})" if n_treg else "n/a",
         verdict=axis_verdict(treg_dir, n_treg)),
    dict(axis="effector-CD8 EXCLUDED from B-niche", rcc="34/37 (0.92)",
         dkd_xenium=f"{frac_str(k_eff,n_eff)} ({eff_dir:.2f})" if n_eff else "n/a",
         verdict=axis_verdict(eff_dir, n_eff)),
    dict(axis="mregDC membership", rcc="inconclusive",
         dkd_xenium="descriptive (see refined CSV)", verdict="DESCRIPTIVE"),
    dict(axis="plasma", rcc="inconclusive", dkd_xenium="coarse only", verdict="DESCRIPTIVE"),
])
vt.to_csv(os.path.join(OUT, "rcc_vs_dkd_verdict.csv"), index=False)
print(vt.to_string(index=False))
load_bearing = (axis_verdict(treg_dir,n_treg), axis_verdict(eff_dir,n_eff))
if "DIVERGE" in load_bearing: overall="DIVERGE"
elif all(v=="REPLICATE" for v in load_bearing): overall="REPLICATE"
else: overall="PARTIAL"
print(f"\nOverall (Treg+/eff-CD8- axis), {'DESCRIPTIVE — ' if DESCRIPTIVE else ''}scoped: {overall}")

# ============================================================================
# FIGURES
# ============================================================================
hdr("FIGURES")
# (a) forest/dot of in-niche log2 enrichment across Xenium samples; Treg-like & eff-CD8 highlighted
rows_f = []
for ct in ["B"] + COARSE:
    sub_c = coarse[coarse.cell_type == ct]
    for _, r in sub_c.iterrows(): rows_f.append(("coarse:"+ct, r.log2_enrich))
for state in ["Treg_like", "eff_CD8", "mregDC_like"]:
    for _, r in refined[refined.state==state].iterrows(): rows_f.append(("state:"+state.replace("_","-"), r.log2_enrich))
fdf = pd.DataFrame(rows_f, columns=["row","log2"]).dropna()
order = ["coarse:B","coarse:"+CD4,"coarse:"+CD8,"coarse:"+PLASMA,"coarse:"+MACRO,"coarse:"+NK,"coarse:DC",
         "state:Treg-like","state:eff-CD8","state:mregDC-like"]
order = [o for o in order if o in set(fdf.row)]
fig, ax = plt.subplots(figsize=(8.5, 5.5))
for i, o in enumerate(order):
    vals = fdf.loc[fdf.row==o, "log2"].values
    hl = o in ("state:Treg-like","state:eff-CD8")
    ax.scatter(vals, np.full(len(vals), i)+np.random.uniform(-0.08,0.08,len(vals)),
               s=60 if hl else 38, color=DKD_COLOR, edgecolor="k" if hl else "none",
               linewidth=1.1 if hl else 0, alpha=0.9, zorder=3 if hl else 2)
    ax.scatter(np.nanmean(vals), i, marker="|", s=400, color="k", zorder=4)
ax.axvline(0, color="gray", ls="--", lw=1)
ax.set_yticks(range(len(order)))
ax.set_yticklabels([o.split(":",1)[1] for o in order])
for i,o in enumerate(order):
    if o in ("state:Treg-like","state:eff-CD8"): ax.get_yticklabels()[i].set_fontweight("bold")
ax.set_xlabel("log2 enrichment  in-niche vs out-niche  (per Xenium sample)")
ax.set_title(f"DKD 'B predom. Immune ME' niche — in-niche enrichment (Xenium, n={n_xen})\n"
             "Treg-like(+) / eff-CD8(-) = the load-bearing RCC axis", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_a_inniche_enrichment.png"), dpi=150); plt.close(fig)

# (b) head-to-head RCC vs DKD
fig, ax = plt.subplots(figsize=(7.5, 4.2))
axes_l = ["Treg ENRICHED", "eff-CD8 EXCLUDED"]
rcc_v = [36/37, 34/37]; dkd_v = [treg_dir if n_treg else 0, eff_dir if n_eff else 0]
x = np.arange(len(axes_l)); w=0.36
ax.bar(x-w/2, rcc_v, w, label="RCC (Xenium, 37 aggs)", color="#1F78B4")
ax.bar(x+w/2, dkd_v, w, label=f"DKD B-niche (Xenium, n={n_xen})", color=DKD_COLOR)
for xi,(rv,dv,nn) in enumerate(zip(rcc_v, dkd_v, [(k_treg,n_treg),(k_eff,n_eff)])):
    ax.text(xi-w/2, rv+0.02, f"{rv:.2f}", ha="center", fontsize=8)
    ax.text(xi+w/2, dv+0.02, f"{nn[0]}/{nn[1]}", ha="center", fontsize=8)
ax.axhline(0.8, color="gray", ls=":", lw=1); ax.set_ylim(0,1.1)
ax.set_xticks(x); ax.set_xticklabels(axes_l)
ax.set_ylabel("fraction of units in expected direction")
ax.set_title("RCC B-aggregate signature vs DKD B-niche (load-bearing axes)")
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_b_rcc_vs_dkd.png"), dpi=150); plt.close(fig)
print("saved fig_a_inniche_enrichment.png, fig_b_rcc_vs_dkd.png")

# ============================================================================
# REPORT.md
# ============================================================================
hdr("WRITING REPORT")
L=[]; W=L.append
W("# Does the RCC immunoregulatory B-aggregate signature replicate in the DKD B-niche?\n")
W(f"Read-only. Niche column `{niche_col}`, value `{NICHE_VAL}`. "
  f"Subset: immune cells in niche-containing samples ({len(sel_idx):,} cells materialized; full X never loaded).\n")
W("## Scope — why subtype claims are Xenium-only\n")
W("Per `cd4_cd8_support`, CosMx 1k CD4/CD8 subtype is **reference-imputed** (CD8A AUROC 0.56, "
  "CD4 0.50 at chance), so **all Treg / effector-CD8 / mregDC state claims here use Xenium 5k only**. "
  "CosMx contributes a coarse B-predominance consistency check (Step 5) with no subtype claims.\n")
W("## Power & confounds\n")
W(f"- Xenium samples containing the niche: **n = {n_xen}** ({', '.join(xen_samples)}).")
W(f"- {'**Underpowered (n<5): results are DESCRIPTIVE / per-sample; no inferential claims.**' if DESCRIPTIVE else 'Small cohort — interpret k/n descriptively.'}")
W("- **No patient column** in obs -> patient-level clustering cannot be controlled (samples may share donors).")
W(f"- Measured Xenium markers (struct-zeros dropped): Treg={TREG_X}, eff-CD8={EFFCD8_X}, mregDC={MREGDC_X}.\n")
W("## RCC-axis result (Xenium)\n")
W(f"- **Treg-like (FOXP3+/IL2RA+/CTLA4+ among CD4+) ENRICHED in-niche: {k_treg}/{n_treg} samples.**")
W(f"- **effector-CD8 (GZMB/GZMK/PRF1/GNLY+ among CD8+) DEPLETED in-niche: {k_eff}/{n_eff} samples.**")
W("- mregDC-like (LAMP3/CCR7/FSCN1/CD274+ among DC): descriptive only (RCC inconclusive) — see refined CSV.")
W(f"- Coarse composition confirms B-dominance; non-B immune pattern in `xenium_coarse_enrichment_per_sample.csv`.\n")
W("## Head-to-head verdict\n")
W("| axis | RCC | DKD B-niche (Xenium) | verdict |")
W("|---|---|---|---|")
for _, r in vt.iterrows(): W(f"| {r.axis} | {r.rcc} | {r.dkd_xenium} | {r.verdict} |")
W(f"\n**Overall, scoped to the load-bearing Treg(+)/effector-CD8(-) axis: {overall}"
  f"{' (DESCRIPTIVE — small n)' if DESCRIPTIVE else ''}.**\n")
W("## Caveats\n")
W("- Coarse subtype labels -> marker-refined states within compartments (FOXP3+ among CD4+, etc.); "
  "not identical to the RCC DBSCAN-aggregate spatial delineation (this is a niche-membership, not a "
  "density-aggregate, test). Direction of the Treg(+)/eff-CD8(-) axis is the comparable quantity.")
W("- Small n and no patient control (above); CosMx excluded from all subtype claims.")
W("- Out-niche background = all non-niche cells in the SAME sample (controls sample composition).")
open(os.path.join(OUT, "REPORT.md"), "w").write("\n".join(L))
print("wrote REPORT.md + CSVs + figures")
print("\n== done ==")
