#!/usr/bin/env python
"""STEP 5 — BAFF/APRIL conditioned re-assessment. The prior global null washes out producer-specific
signal; here we CONDITION on producers (myeloid/macro primary; stroma/fibroblast) and on SPACE
(near vs far from B-lineage aggregates). Ligands TNFSF13B(BAFF)/TNFSF13(APRIL); receptors
TNFRSF17(BCMA)/TNFRSF13B(TACI)/TNFRSF13C(BAFF-R). Verdict GO/NO-GO, honest either way. Read-only."""
import os, numpy as np, pandas as pd
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
RE=f"{REPO}/analysis/dkd_xenium_reannotation"; OUT=f"{REPO}/analysis/dkd_xenium_disease"
EPS,MINPTS=50,20
LIG={"TNFSF13B":"BAFF","TNFSF13":"APRIL"}; REC={"TNFRSF17":"BCMA","TNFRSF13B":"TACI","TNFRSF13C":"BAFF-R"}
GENES=list(LIG)+list(REC)

cells=pd.read_parquet(f"{RE}/cells.parquet").reset_index(drop=True)
lig=pd.read_parquet(f"{OUT}/ligand_me_percell.parquet").reset_index(drop=True)
assert (cells.cell_id.astype(str).values==lig.cell_id.astype(str).values).all()
for g in GENES: cells[g]=lig[g].values
cells["is_Blin"]=cells.my_label.isin(["B","Plasma"])
brich=["1006","HK2695"]

# ---------- 1. detection by cell type (producers vs ambient floor = epithelium) ----------
def det(mask,g):
    s=cells.loc[mask,g]; return (s>0).mean()
groups={"Myeloid":cells.my_label=="Myeloid","Fibroblast":cells.my_label=="Fibroblast",
        "VSMC/MC":cells.my_label.isin(["VSMC","MC"]),"Endothelial":cells.my_lineage=="Endothelial",
        "Epithelial (floor)":cells.my_lineage=="Epithelial","B":cells.my_label=="B",
        "Plasma":cells.my_label=="Plasma","CD4 T":cells.my_label=="CD4 T","CD8 T":cells.my_label=="CD8 T"}
rows=[]
floor={g:det(groups["Epithelial (floor)"],g) for g in GENES}
for gname,m in groups.items():
    r={"cell_type":gname,"n":int(m.sum())}
    for g in GENES:
        d=det(m,g); r[f"{LIG.get(g,REC.get(g))}_det%"]=round(d*100,2); r[f"{LIG.get(g,REC.get(g))}_vsFloor"]=round(d/(floor[g]+1e-9),1)
    rows.append(r)
dett=pd.DataFrame(rows); dett.to_csv(f"{OUT}/baff_detection_by_celltype.csv",index=False)
pd.set_option("display.width",240)
print("=== detection by cell type (BAFF/APRIL = ligands; vsFloor = x epithelial ambient) ===")
print(dett[["cell_type","n","BAFF_det%","BAFF_vsFloor","APRIL_det%","APRIL_vsFloor","BCMA_det%","TACI_det%","BAFF-R_det%"]].to_string(index=False))

# ---------- 2. producer ligand detection NEAR vs FAR from B-lineage aggregates ----------
def agg_members(g):
    B=g[g.is_Blin];
    if len(B)<MINPTS: return np.empty((0,2))
    lab=DBSCAN(eps=EPS,min_samples=MINPTS).fit(B[["spatial_x","spatial_y"]].values).labels_
    return B[["spatial_x","spatial_y"]].values[lab!=-1]
def near_far(sample_ids,producer_mask,gene):
    near_pos=near_tot=far_pos=far_tot=0
    for sid in sample_ids:
        g=cells[cells.orig_ident==sid]; mem=agg_members(g)
        if len(mem)==0: continue
        prod=g[producer_mask.loc[g.index].values]
        if len(prod)==0: continue
        d,_=cKDTree(mem).query(prod[["spatial_x","spatial_y"]].values,k=1)
        pos=(prod[gene].values>0); near=d<=EPS
        near_pos+=int(pos[near].sum()); near_tot+=int(near.sum())
        far_pos+=int(pos[~near].sum()); far_tot+=int((~near).sum())
    return (near_pos,near_tot,near_pos/max(near_tot,1)),(far_pos,far_tot,far_pos/max(far_tot,1))
prodmask={"Myeloid":cells.my_label=="Myeloid","Stroma":cells.my_lineage=="Stroma"}
aggbear=cells.groupby("orig_ident",observed=True).apply(lambda g:len(agg_members(g))>0)
aggbear=list(aggbear[aggbear].index)
rows=[]
for setname,sids in [("B-rich DKD",brich),("all aggregate-bearing",aggbear)]:
    for pname,pm in prodmask.items():
        for g in LIG:
            (np_,nt,nr),(fp,ft,fr)=near_far(sids,pm,g)
            rows.append(dict(producer=pname,ligand=LIG[g],sample_set=setname,
                near_det=f"{np_}/{nt}",near_pct=round(nr*100,2),far_det=f"{fp}/{ft}",far_pct=round(fr*100,2),
                enrich=round(nr/(fr+1e-9),2)))
nf=pd.DataFrame(rows); nf.to_csv(f"{OUT}/baff_near_far.csv",index=False)
print("\n=== producer ligand detection NEAR (<=50um) vs FAR from B-aggregates ===")
print(nf.to_string(index=False))

# ---------- 3. receptors: aggregate-concentration (B-lineage near vs far own aggregates) ----------
rows=[]
for g in REC:
    (np_,nt,nr),(fp,ft,fr)=near_far(aggbear,cells.is_Blin,g)   # receptors on B-lineage cells
    # producer-side receptor on all cells near vs far
    rows.append(dict(receptor=REC[g],on="B-lineage",near_pct=round(nr*100,2),far_pct=round(fr*100,2),enrich=round(nr/(fr+1e-9),2)))
rt=pd.DataFrame(rows); rt.to_csv(f"{OUT}/baff_receptors.csv",index=False)
print("\n=== receptors on B-lineage, near vs far own aggregates ===")
print(rt.to_string(index=False))

# ---------- 4. verdict ----------
my_baff=nf[(nf.producer=="Myeloid")&(nf.ligand=="BAFF")&(nf.sample_set=="B-rich DKD")].iloc[0]
my_april=nf[(nf.producer=="Myeloid")&(nf.ligand=="APRIL")&(nf.sample_set=="B-rich DKD")].iloc[0]
print("\n=== VERDICT ===")
for lab,row in [("BAFF",my_baff),("APRIL",my_april)]:
    floorpct=dett.loc[dett.cell_type=="Epithelial (floor)",f"{lab}_det%"].iloc[0]
    prodpct=dett.loc[dett.cell_type=="Myeloid",f"{lab}_det%"].iloc[0]
    print(f"{lab}: myeloid {prodpct}% vs epithelial floor {floorpct}% ; near-agg {row.near_pct}% vs far {row.far_pct}% (enrich {row.enrich}x)")
print("== baff_april done ==")
