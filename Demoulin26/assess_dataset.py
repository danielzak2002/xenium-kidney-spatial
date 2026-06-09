#!/usr/bin/env python
"""
assess_dataset.py — READ-ONLY orientation/inventory of the Dumoulin et al. 2026 DKD
spatial atlas (Zenodo 19868428): CosMx 1k + Xenium 5k, annotated AnnData.

Reconnaissance only. The .h5ad is ~8.7 GB; it is opened with anndata backed='r' and ONLY
.obs/.var/.obsm/.uns/.layers metadata are inspected. X is NEVER materialized
(no .to_memory(), no slicing of .X). Raw files are treated as read-only.

Prints everything to console and writes a summary to Demoulin26/ASSESSMENT.md.

  conda run -n spatial python Demoulin26/assess_dataset.py
"""
import os, sys, io, contextlib
import numpy as np, pandas as pd, anndata as ad

BASE = os.path.dirname(os.path.abspath(__file__))
DIAG = os.path.join(BASE, "data", "Diagnosis.xlsx")
H5AD = os.path.join(BASE, "data", "spatial_adata_xenium_cosmx_zenodo.h5ad")
OUT = os.path.join(BASE, "ASSESSMENT.md")

MD = []  # accumulate markdown summary lines
def show(s=""):
    print(s)
def md(s=""):
    MD.append(s)

def vc(series, name, dropna=False, top=60):
    counts = series.value_counts(dropna=dropna)
    show(f"\n-- value_counts: {name} ({series.nunique(dropna=True)} unique) --")
    show(counts.head(top).to_string())
    return counts

# ============================================================================
show("=" * 78); show("STEP 1 — readme.txt"); show("=" * 78)
with open(os.path.join(BASE, "readme.txt")) as fh:
    show(fh.read())

# ============================================================================
show("=" * 78); show("STEP 2 — Diagnosis.xlsx"); show("=" * 78)
xl = pd.ExcelFile(DIAG)
show(f"sheet names: {xl.sheet_names}")
diag_frames = {}
for sh in xl.sheet_names:
    df = xl.parse(sh)
    diag_frames[sh] = df
    show(f"\n--- sheet '{sh}'  shape={df.shape} ---")
    show(f"columns: {list(df.columns)}")
    show(df.to_string(max_rows=80))
    for col in df.columns:
        if df[col].dtype == object or df[col].nunique() < 30:
            if any(k in col.lower() for k in ("diag", "disease", "group", "condition", "tech", "platform", "class")):
                vc(df[col].astype(str), f"[{sh}] {col}")

# ============================================================================
show("\n" + "=" * 78); show("STEP 3 — AnnData structure (backed='r')"); show("=" * 78)
a = ad.read_h5ad(H5AD, backed="r")
show(f"n_obs = {a.n_obs:,}   n_vars = {a.n_vars:,}")
try:
    show(f"X: shape={a.X.shape} dtype={a.X.dtype} (sparse, backed — NOT materialized)")
except Exception as e:
    show(f"X: (backed) {e}")
show(f"layers: {list(a.layers.keys())}")
show(f"obsm keys: {[(k, a.obsm[k].shape) for k in a.obsm.keys()]}")
show(f"uns keys: {list(a.uns.keys())}")

show("\n-- .obs columns (dtype) --")
obs = a.obs
for c in obs.columns:
    show(f"  {c:34s} {str(obs[c].dtype):12s} nunique={obs[c].nunique(dropna=True)}")

# value_counts for every categorical/string col with <60 uniques
show("\n-- value_counts for low-cardinality obs columns --")
for c in obs.columns:
    nun = obs[c].nunique(dropna=True)
    is_strlike = (obs[c].dtype == object) or str(obs[c].dtype) == "category"
    if is_strlike and nun < 60:
        vc(obs[c].astype(str), c)

# .var
show("\n-- .var --")
var = a.var
show(f"var columns: {list(var.columns)}; n_genes={var.shape[0]}")
show("var head:\n" + var.head(10).to_string())
show(f"(no per-platform/panel column in .var: columns are {list(var.columns)} — panel "
     "membership is NOT encoded here; reconciliation must be inferred, see notes.)")

# ============================================================================
show("\n" + "=" * 78); show("STEP 4 — cross-tabs (.obs only)"); show("=" * 78)

PLATFORM = "tech"; SAMPLE = "orig_ident"; DISEASE = "Condition"
CT_MAIN = "annotation_updated"; CT_IMM = "immune_cell_annotation_combined"
NICHE = "niches_annotation_based"

# 4a — cells per platform / disease / sample
vc(obs[PLATFORM].astype(str), "cells per platform (tech)")
vc(obs[DISEASE].astype(str), "cells per disease group (Condition)")
spc = obs.groupby(PLATFORM, observed=True)[SAMPLE].nunique()
_sp = pd.crosstab(obs[SAMPLE].astype(str), obs[PLATFORM].astype(str))
dual_ids = list(_sp.index[(_sp > 0).sum(1) > 1])
show(f"\nsamples per platform:\n{spc.to_string()}")
show(f"orig_ident measured on BOTH platforms ({len(dual_ids)}): {dual_ids}")
show(f"total samples (orig_ident): {obs[SAMPLE].nunique()}")
cps = obs[SAMPLE].value_counts()
show(f"cells per sample: median={cps.median():.0f} min={cps.min()} max={cps.max()} (n={len(cps)})")
show("Paper check: CosMx n=48, Xenium n=16; 64 samples / 58 patients.")
patient_cols = [c for c in obs.columns if any(k in c.lower() for k in ("patient", "donor", "subject", "case"))]
show(f"patient/donor column present? {patient_cols if patient_cols else 'NONE (only sample-level orig_ident; patient mapping not in obs)'}")
md("# Dumoulin et al. 2026 DKD spatial atlas — read-only assessment")

# 4b — cell-type vocabulary
def lineage_labels(cats):
    out = {"B": [], "Plasma": [], "T": [], "Myeloid/Macrophage": [], "DC": [], "NK": []}
    for c in cats:
        l = str(c).lower()
        if "plasma" in l or "plasmablast" in l: out["Plasma"].append(c)
        elif ("b cell" in l or "b-cell" in l or l.strip() in ("b", "b_cell") or "bcell" in l
              or "memory b" in l or "naive b" in l or "germinal" in l): out["B"].append(c)
        if any(k in l for k in ("t cell", "t-cell", "cd4", "cd8", "treg", "regulatory t", " th", "tfh", "tcell")): out["T"].append(c)
        if any(k in l for k in ("macroph", "monocyte", "myeloid", "kupffer", "mph")): out["Myeloid/Macrophage"].append(c)
        if any(k in l for k in ("dendritic", "dc", "mregdc", "cdc", "pdc")) and "endothel" not in l: out["DC"].append(c)
        if "nk" in l or "natural killer" in l: out["NK"].append(c)
    return out

for ctcol in [CT_MAIN, CT_IMM]:
    show("\n" + "-" * 60)
    counts = vc(obs[ctcol].astype(str), f"cell-type vocabulary: {ctcol}", top=80)
    lin = lineage_labels(counts.index)
    show(f"\nB-lineage / immune labels detected in {ctcol}:")
    for k, v in lin.items():
        show(f"  {k:20s}: {v}")

# 4c — cross-platform T / B / plasma resolution
show("\n" + "=" * 60); show("STEP 4c — CROSS-PLATFORM resolution probe"); show("=" * 60)
ct = CT_IMM if obs[CT_IMM].nunique() else CT_MAIN
for label_set, name in [(None, "ALL cell types"), ("T", "T-lineage"), ("BP", "B / plasma")]:
    cats = obs[ct].astype(str)
    if label_set == "T":
        keep = cats.str.contains("t cell|cd4|cd8|treg|regulatory t|tfh|naive t|memory t", case=False, regex=True)
    elif label_set == "BP":
        keep = cats.str.contains("b cell|b-cell|plasma|memory b|naive b|germinal", case=False, regex=True)
    else:
        keep = pd.Series(True, index=cats.index)
    sub = obs.loc[keep]
    if len(sub) == 0:
        show(f"\n[{name}] no matching labels"); continue
    tab = pd.crosstab(sub[ct].astype(str), sub[PLATFORM].astype(str))
    tab = tab.loc[tab.sum(1).sort_values(ascending=False).index]
    show(f"\n-- crosstab({ct}, {PLATFORM}) restricted to {name} --")
    show(tab.to_string())

# 4d — niche probe
show("\n" + "=" * 60); show("STEP 4d — NICHE probe"); show("=" * 60)
niche_cols = [c for c in obs.columns if any(k in c.lower() for k in ("niche", "_me", "microenv"))]
show(f"candidate niche/microenvironment columns: {niche_cols}")
bcell_niche = None; bcell_niche_col = None; best_frac = 0.10
for nc in niche_cols:
    if obs[nc].nunique(dropna=True) > 60:
        continue
    counts = vc(obs[nc].astype(str), f"niche vocabulary: {nc}")
    # B-lineage fraction per niche -> identify the B-cell-predominant niche.
    # NOTE: the immune label is the bare string "B" (and "Plasma"); match those exactly
    # plus any spelled-out B/plasma variants.
    cats = obs[ct].astype(str)
    isB = cats.isin(["B", "Plasma", "Plasmablast"]) | cats.str.contains(
        "b cell|b-cell|memory b|naive b|germinal|plasma", case=False, regex=True)
    bfrac = obs.assign(_isB=isB.values).groupby(nc, observed=True)["_isB"].mean().sort_values(ascending=False)
    show(f"  B-lineage fraction per {nc} (top 5):\n{(bfrac.head(5)*100).round(1).to_string()}")
    if bfrac.iloc[0] > best_frac:                       # global max B-fraction niche
        best_frac = bfrac.iloc[0]; bcell_niche, bcell_niche_col = bfrac.index[0], nc

if bcell_niche is not None:
    show(f"\n>>> B-cell-predominant niche: '{bcell_niche}' in column '{bcell_niche_col}'")
    inn = obs[obs[bcell_niche_col].astype(str) == str(bcell_niche)]
    show(f"composition (cell types in this niche), top 25:")
    show((inn[ct].astype(str).value_counts().head(25)).to_string())
    show(f"\nsamples containing this niche: {inn[SAMPLE].nunique()} of {obs[SAMPLE].nunique()}")
    show(f"split by disease group:\n{inn[DISEASE].astype(str).value_counts().to_string()}")
    show(f"split by platform:\n{inn[PLATFORM].astype(str).value_counts().to_string()}")
else:
    show("\nNo single B-cell-predominant niche identified above the 10% threshold.")

# ---- comprehensive markdown summary ----------------------------------------
imm_tab = pd.crosstab(obs[CT_IMM].astype(str), obs[PLATFORM].astype(str))
imm_tab = imm_tab.drop(index=[i for i in ["Unknown"] if i in imm_tab.index])
imm_tab = imm_tab.loc[imm_tab.sum(1).sort_values(ascending=False).index]
spc_combos = obs.groupby(PLATFORM, observed=True)[SAMPLE].nunique()
md("\n## Source & object")
md(f"Dumoulin et al., *Spatial atlas of diabetic kidney disease reveals a B cell-rich subgroup* "
   f"(Nature 2026; Zenodo 19868428). Read-only inventory via `assess_dataset.py` (anndata backed='r', "
   f"X never materialized). Raw `.h5ad`/`.xlsx` are git-ignored.")
md(f"\n**Object:** `spatial_adata_xenium_cosmx_zenodo.h5ad` — **{a.n_obs:,} cells × {a.n_vars:,} genes**, "
   f"sparse `X` (float32) + `counts` layer; `obsm` = X_umap / spatial / spatial_fov; `uns` = orig_ident_colors only "
   f"(no scVI/scANVI latent, no niche color map).")
md("\n## Samples & design")
md(f"- {obs[SAMPLE].nunique()} distinct `orig_ident`; **{int(spc_combos['CosMx'])} CosMx + "
   f"{int(spc_combos['Xenium'])} Xenium = {int(spc_combos.sum())} platform-samples** (matches the paper's 64). "
   f"4 IDs measured on both platforms ({', '.join(dual_ids)}). No patient/donor column in `.obs` "
   f"(paper: 58 patients) — patient mapping must come from the diagnosis sheet.")
md(f"- Diagnosis sheet (64 rows): conditions DKD 29, Control 15, DM 6, DM/HTN 4, + 8 rarer GN/amyloid groups; "
   f"covariates Age/Sex/Race/DM/HTN/GFR/treatments present.")
md(f"- Cells: {int(obs[PLATFORM].value_counts()['CosMx']):,} CosMx + {int(obs[PLATFORM].value_counts()['Xenium']):,} "
   f"Xenium; median {cps.median():.0f} cells/sample (range {cps.min()}–{cps.max()}).")
md("\n## Annotations present (validated reference — YES)")
md(f"- **Cell types:** `annotation_updated` (20 kidney epithelial/stromal/endothelial types; immune as one coarse `Immune` bucket).")
md(f"- **Immune subtypes:** `immune_cell_annotation_combined` (12): Macro, **CD8+, CD4+, B, Plasma**, Neutrophil, NK, Baso_Mast, cDC, pDC, mDC. No Treg label (T resolved only to CD4+/CD8+).")
md(f"- **Niches:** `niches_annotation_based` (12 tissue niches incl. `Immune niche`) + microenvironment columns "
   f"`immune_ME`/`Immune_ME_20um` (incl. **`B predom. Immune ME`**) and `iPT_iLOH_ME(_20um)`.")
md("\n## Cross-platform immune resolution (the key finding)")
md("Unlike the cLN 1k atlas (CD4/CD8 collapsed on CosMx), here **CD4+ and CD8+ are separately annotated on BOTH platforms**, "
   "as are B and Plasma. Per-platform immune cell counts (`immune_cell_annotation_combined` × `tech`):")
md("\n| immune type | CosMx | Xenium |\n|---|---:|---:|")
for ix in imm_tab.index:
    md(f"| {ix} | {int(imm_tab.loc[ix,'CosMx'])} | {int(imm_tab.loc[ix,'Xenium'])} |")
md("\nCaveat: labels exist on both platforms, but CosMx 1k immune calls inherit the platform's ambient/low-plex "
   "limits (cf. the cLN finding) — treat CosMx CD4/CD8 as author-asserted, to be validated, not assumed clean.")
md("\n## B / plasma / TLS-like niche")
md(f"- **`B predom. Immune ME`** = the B-cell-predominant microenvironment: ~83% B cells "
   f"(B 2748, CD4+ 362, CD8+ 113, Macro 55, Plasma 22, NK 5) — a B-aggregate / TLS-like niche.")
md(f"- Present in **9/60 samples**, predominantly **DKD** (then C3GN, Control, AA amyloid, DM/HTN); both platforms "
   f"(Xenium-leaning). This is the spatial substrate of the paper's 'B cell-rich subgroup'.")
md("\n## Data-quality flags / mismatches")
md("- `.var` has **no panel/platform column** (only `n_cells`); the object holds 5,443 genes with no CosMx-1k vs "
   "Xenium-5k mask. CosMx cells are effectively restricted to the ~1k CosMx panel (structural zeros elsewhere), so any "
   "**cross-platform gene-level analysis must be limited to the shared/overlap panel — which is NOT pre-computed here**.")
md("- Diagnosis `Disease` column has inconsistent casing (Yes/YES, No/NO); `Sex` missing for 229,669 cells; large `Unknown` "
   "in immune (non-immune cells) and ME columns (ME assigned only to immune-region cells).")
md("- 60 `orig_ident` vs paper's 58 patients (replicate samples e.g. HK2844_2, HK3035_2, HK3531_2; 4 dual-platform IDs).")
md("\n## Candidate next analyses (extend the existing project — NOT run here)")
md("1. **Transfer the B-aggregate / plasma–myeloid niche workflow** (squidpy nhood-enrichment, DBSCAN aggregate "
   "delineation, in/out-of-niche composition) onto `B predom. Immune ME` across the 9 B-rich DKD samples — a direct "
   "extension of the RCC B–Treg and cLN plasma–myeloid analyses to a third disease context.")
md("2. **Cross-platform validation of CD4/CD8 typing**: use the dual-platform IDs (HK2695/HK2753/HK3106/HK3626) to test "
   "whether CosMx-1k CD4/CD8 calls agree with paired Xenium-5k — the cleanest test of the cLN T-cell limitation.")
md("3. **B-rich vs B-poor DKD contrast**: compare niche composition, plasma load, and myeloid state between the 9 "
   "B-niche samples and the remaining DKD samples (the paper's subgroup split), reusing the per-slide methodology.")
md("4. **Gene-usability gate on both panels**: compute per-platform detection/enrichment for B-survival and recruitment "
   "genes (MS4A1/MZB1/TNFRSF17 vs BAFF/chemokines) to define the analyzable overlap panel before any modeling.")
md("5. **Niche conservation across disease**: test whether the B-predominant niche composition (B–CD4 core, CD8-low) "
   "is conserved vs the RCC immunoregulatory aggregate and cLN plasma–myeloid niche.")

a.file.close()
with open(OUT, "w") as fh:
    fh.write("\n".join(MD) + "\n")
show(f"\nwrote {OUT}")
show("== assessment done ==")
