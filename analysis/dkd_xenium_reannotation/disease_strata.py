#!/usr/bin/env python
"""STEP 4 — disease-status tagging (SETUP only, no association). Join Diagnosis.xlsx to the
16 Xenium samples by Sample ID = '<orig_ident>_Xenium'. Report join coverage, strata counts
(samples + cells), and covariate availability among these 16. Read-only raw."""
import os, pandas as pd, numpy as np
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
OUT=f"{REPO}/analysis/dkd_xenium_reannotation"
dx=pd.read_excel(f"{REPO}/Demoulin26/data/Diagnosis.xlsx")
pc=pd.read_csv(f"{OUT}/per_sample_counts.csv")          # 16 Xenium samples + n_cells (from STEP0)
pc["orig_ident"]=pc["orig_ident"].astype(str)
pc=pc.drop(columns=[c for c in ["Condition"] if c in pc.columns])  # Condition is authoritative from the sheet
dx["orig_ident"]=dx["Sample ID"].astype(str).str.replace("_Xenium","",regex=False).str.replace("_CosMx","",regex=False)
# restrict the sheet to the Xenium rows
dx_x=dx[dx["Sample ID"].astype(str).str.endswith("_Xenium")].copy()

j=pc.merge(dx_x,on="orig_ident",how="left",indicator=True)
cov=(j._merge=="both").mean()
print(f"join coverage: {(j._merge=='both').sum()}/{len(j)} samples matched ({cov*100:.0f}%)")
miss=j[j._merge!="both"]["orig_ident"].tolist()
if miss: print("UNMATCHED:",miss)

cols=["orig_ident","n_cells","Condition","Disease","DM","HTN","GFR","Age","Sex","Race",
      "Treatment: RAAS","Treatment: MRA","Treatment: SGLT2i"]
tab=j[cols].sort_values(["Condition","orig_ident"]).reset_index(drop=True)
tab.to_csv(f"{OUT}/disease_strata_per_sample.csv",index=False)
print("\n=== per-sample disease table ===")
print(tab.to_string(index=False))

def strata(col):
    g=j.groupby(col,observed=True).agg(n_samples=("orig_ident","nunique"),n_cells=("n_cells","sum")).sort_values("n_samples",ascending=False)
    return g
print("\n=== strata: Condition ===");  print(strata("Condition").to_string())
print("\n=== strata: Disease (Yes/No) ===");  print(strata("Disease").to_string())
print("\n=== strata: GFR band ===");  print(strata("GFR").to_string())
print("\n=== strata: DM ===");  print(strata("DM").to_string())
print("\n=== strata: HTN ===");  print(strata("HTN").to_string())

# scoping verdict
nb={k:int((j.Condition==k).sum()) for k in j.Condition.unique()}
dkd=nb.get("DKD",0); ctrl=nb.get("Control",0); other=len(j)-dkd-ctrl
summary=pd.DataFrame([dict(strata="DKD",n_samples=dkd),dict(strata="Control",n_samples=ctrl),
                      dict(strata="other glomerulopathy",n_samples=other)])
summary.to_csv(f"{OUT}/disease_strata_summary.csv",index=False)
print(f"\nSCOPING: DKD={dkd}, Control={ctrl}, other(GN/amyloid)={other}.")
print("16 samples cannot power a graded severity model; the supportable contrasts are binary:")
print(" (a) DKD vs Control, (b) DKD vs non-DKD, (c) B-rich vs B-poor DKD subgroup (per the paper).")
print("== disease_strata done ==")
