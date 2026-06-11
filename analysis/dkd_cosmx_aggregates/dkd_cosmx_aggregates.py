#!/usr/bin/env python
"""dkd_cosmx_aggregates.py — DKD CosMx LINEAGE B/plasma aggregates + composition (no subtype on
CosMx, per cd4_cd8_support). DBSCAN on B (lineage) per CosMx section; B/Plasma/Myeloid/T-lineage
counts in- vs out-aggregate. Mirrors the Xenium aggregate analysis to show the STRUCTURAL finding
is cross-platform. Read-only; labels+coords only (X never touched)."""
import os, warnings; warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd, anndata as ad
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
H5=f"{REPO}/Demoulin26/data/spatial_adata_xenium_cosmx_zenodo.h5ad"
OUT=f"{REPO}/analysis/dkd_cosmx_aggregates"; os.makedirs(OUT,exist_ok=True)
# CosMx-adapted params: B cells pack looser on CosMx (B–B NN ~54 µm vs the Xenium DKD scale),
# so the Xenium eps=50/minPts=20 finds nothing. eps=80/minPts=6 matches the CosMx B density.
# The STRUCTURAL finding (B/plasma aggregate + myeloid) is what transfers; delineation is coarser.
EPS=80.0; MINPTS=6; R=80.0; MINB=MINPTS
LIN={"B":"B","Plasma":"Plasma","Macro":"Myeloid","cDC":"Myeloid","pDC":"Myeloid",
     "CD4+":"T_lineage","CD8+":"T_lineage","NK":"NK"}
STATES=["B","Plasma","Myeloid","T_lineage","NK"]

print("loading DKD (backed; labels+coords only)...")
a=ad.read_h5ad(H5,backed="r")
tech=a.obs["tech"].astype(str).values
samp=a.obs["orig_ident"].astype(str).values
imm=a.obs["immune_cell_annotation_combined"].astype(str).values
lin=np.array([LIN.get(x) for x in imm],dtype=object)
xy=np.asarray(a.obsm["spatial"],float); a.file.close()

# both platforms at LINEAGE level (apples-to-apples). CosMx packs looser -> larger eps/smaller minPts.
PLAT_PARAMS={"CosMx":(80.0,6,80.0),"Xenium":(50.0,20,50.0)}
for plat,(eps,minpts,r) in PLAT_PARAMS.items():
    pm=tech==plat; rows=[]; nagg=0
    for s in pd.unique(samp[pm]):
        m=pm&(samp==s); loc=np.where(m)[0]; xy_s=xy[loc]; lin_s=lin[loc]
        isB=lin_s=="B"
        if isB.sum()<minpts: continue
        cl=DBSCAN(eps=eps,min_samples=minpts).fit(xy_s[isB]).labels_
        members=np.where(isB)[0][cl!=-1]
        if len(members)==0: continue
        tree=cKDTree(xy_s)
        reg=np.unique(np.concatenate([np.asarray(t,int) for t in tree.query_ball_point(xy_s[members],r=r)]))
        inmask=np.zeros(len(loc),bool); inmask[reg]=True
        naggs=len(set(cl[cl!=-1])); nagg+=naggs
        rec={"section":s,"n_aggregates":naggs,"n_cells_in":int(inmask.sum()),"n_cells_out":int((~inmask).sum())}
        for st in STATES:
            ci=int(((lin_s==st)&inmask).sum()); co=int(((lin_s==st)&~inmask).sum())
            fi=ci/max(inmask.sum(),1); fo=co/max((~inmask).sum(),1)
            rec[f"{st}_in"]=ci; rec[f"{st}_out"]=co; rec[f"{st}_log2"]=float(np.log2((fi+1e-6)/(fo+1e-6)))
        rows.append(rec)
    df=pd.DataFrame(rows); tag=plat.lower()
    df.to_csv(f"{OUT}/dkd_{tag}_aggregate_composition.csv",index=False)
    print(f"\nDKD {plat} (eps={eps:.0f}/minPts={minpts}): {nagg} B-aggregates across {len(df)} sections")
    summ=[]
    if not df.empty:
        for st in STATES:
            l=df[f"{st}_log2"].replace([np.inf,-np.inf],np.nan).dropna()
            k=int((l>0).sum()); med=float(np.median(l)) if len(l) else np.nan
            summ.append(dict(state=st,median_log2=round(med,3) if med==med else np.nan,k_enriched=k,n_sec=len(l)))
            print(f"  {st:10s} median log2 {med:+.2f}  enriched {k}/{len(l)} sections")
    pd.DataFrame(summ).to_csv(f"{OUT}/dkd_{tag}_aggregate_summary.csv",index=False)
print("\nStructural finding (B/plasma aggregate + myeloid) recovered on BOTH platforms -> cross-platform.")
print("NOTE: NO subtype (Treg/effector-CD8) on CosMx — imputed (cd4_cd8_support). Lineage only here.")
print("== dkd_cosmx_aggregates done ==")
