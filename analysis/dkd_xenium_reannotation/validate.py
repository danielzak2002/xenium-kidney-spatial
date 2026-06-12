#!/usr/bin/env python
"""STEP 3 — validation vs author labels (confidence check). Per-cell confusion of MY labels vs
annotation_updated (segment + coarse Immune) and vs immune_cell_annotation_combined (immune
subtypes). ARI (label-invariant) + per-type recall/precision after harmonizing names to a shared
vocabulary. Divergences are reported as findings (e.g. Treg is mine-only -> folds into author CD4)."""
import os, numpy as np, pandas as pd
from sklearn.metrics import adjusted_rand_score
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
OUT=f"{REPO}/analysis/dkd_xenium_reannotation"
cells=pd.read_parquet(f"{OUT}/cells.parquet")

# ---- shared coarse vocabulary (segment + coarse immune) ----
def harm_coarse(x):
    # shared vocabulary: keep injury states (iPT/iTAL) distinct; merge genuinely fuzzy boundaries
    # (CNT+PC principal/connecting; IC A+B; all EC subtypes) since marker support / panel limits
    # make finer splits unreliable. MC1->MC. EC_Lymph and DTL_ATL kept where present.
    x=str(x)
    m={"PT":"PT","iPT":"iPT","TAL":"TAL","iTAL":"iTAL","DTL_ATL":"DTL_ATL","DCT":"DCT",
       "CNT":"PC/CNT","PC":"PC/CNT","IC A":"IC","IC B":"IC","Podo":"Podo","PEC":"PEC",
       "Fibroblast":"Fibroblast","VSMC":"VSMC","MC":"MC","MC1":"MC","Immune":"Immune"}
    if x.startswith("EC"): return "EC"
    return m.get(x,x)
cells["my_coarse"]=np.where(cells.is_immune_cluster,"Immune",cells.global_label).astype(str)
cells["my_coarse"]=cells.my_coarse.map(harm_coarse)
cells["auth_coarse"]=cells.author_annotation.map(harm_coarse)

def confusion(true,pred,fn_prefix):
    tab=pd.crosstab(true,pred)
    tab.to_csv(f"{OUT}/{fn_prefix}_confusion_counts.csv")
    rec=tab.div(tab.sum(1),axis=0)        # recall per author (true) class
    prec=tab.div(tab.sum(0),axis=1)       # precision per my (pred) class
    rec.to_csv(f"{OUT}/{fn_prefix}_recall.csv"); prec.to_csv(f"{OUT}/{fn_prefix}_precision.csv")
    ari=adjusted_rand_score(true,pred)
    # per-type recall/precision where names match
    shared=[t for t in tab.index if t in tab.columns]
    pt=pd.DataFrame({"type":shared,
        "n_author":[int(tab.loc[t].sum()) for t in shared],
        "recall":[round(float(rec.loc[t,t]),3) for t in shared],
        "precision":[round(float(prec.loc[t,t]),3) if t in prec.columns else np.nan for t in shared]})
    # agreement = fraction of cells where harmonized names equal
    agree=float((np.asarray(true)==np.asarray(pred)).mean())
    return ari,agree,pt,tab

# ===== Confusion A: segment-level vs annotation_updated =====
ariA,agrA,ptA,tabA=confusion(cells.auth_coarse.values,cells.my_coarse.values,"coarse")
print(f"[A segment] ARI={ariA:.3f}  agreement={agrA*100:.1f}%  n={len(cells):,}")
print(ptA.to_string(index=False))

# ===== Confusion B: immune subtypes vs immune_cell_annotation_combined =====
def harm_imm(x):
    x=str(x)
    m={"Macro":"Myeloid/Macro","Myeloid":"Myeloid/Macro","CD4+":"CD4","CD4 T":"CD4","Treg":"CD4",
       "CD8+":"CD8","CD8 T":"CD8","B":"B","Plasma":"Plasma","NK":"NK","Neutrophil":"Neutrophil",
       "cDC":"DC","pDC":"DC","mDC":"DC","DC":"DC","Baso_Mast":"Mast/Baso","Mast_Baso":"Mast/Baso"}
    return m.get(x,x)
# author-immune-labeled cells (the authors called an immune subtype here)
im=cells[cells.author_immune.isin(["Macro","CD8+","CD4+","B","Plasma","Neutrophil","NK","Baso_Mast","cDC","pDC"])].copy()
im["auth_imm"]=im.author_immune.map(harm_imm)
im["my_imm"]=np.where(im.my_immune_label.notna(),im.my_immune_label.map(harm_imm),"non-immune (missed)")
ariB,agrB,ptB,tabB=confusion(im.auth_imm.values,im.my_imm.values,"immune")
print(f"\n[B immune] ARI={ariB:.3f}  agreement={agrB*100:.1f}%  n={len(im):,} (author-immune cells)")
print(ptB.to_string(index=False))
# Treg divergence: where do MY Treg cells land in author labels?
treg=cells[cells.my_immune_label=="Treg"]
if len(treg):
    print(f"\nTreg divergence (mine-only): {len(treg):,} Treg cells; author_immune of those:")
    print(treg.author_immune.value_counts().head(6).to_string())

summ=pd.DataFrame([
  dict(comparison="segment_vs_annotation_updated",ARI=round(ariA,3),agreement=round(agrA,3),n_cells=len(cells)),
  dict(comparison="immune_vs_immune_combined",ARI=round(ariB,3),agreement=round(agrB,3),n_cells=len(im))])
summ.to_csv(f"{OUT}/validation_summary.csv",index=False)
if len(treg): treg.author_immune.value_counts().to_csv(f"{OUT}/treg_divergence.csv")
cells.to_parquet(f"{OUT}/cells.parquet")   # persist my_coarse / auth_coarse for the figures
print("\n",summ.to_string(index=False)); print("== validate done ==")
