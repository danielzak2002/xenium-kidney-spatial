#!/usr/bin/env python
"""Stress-test the peri-aggregate STROMAL BAFF signal against ambient / transcript-density spillover.
STEP 1 local-ambient control: stromal vs non-producer epithelial BAFF in the SAME peri-aggregate
window. STEP 2 near(<=50um) vs far(>200um) within stromal, per section. STEP 3 density confound:
near/far BAFF ratio per cell type vs the near/far transcript-density ratio. No Xenium negative-control
probes survive in this object (only the gene NEGR1) -> epithelium is the local non-producer ambient
floor and nCount is the density proxy. Read-only raw."""
import os, numpy as np, pandas as pd
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
RE=f"{REPO}/analysis/dkd_xenium_reannotation"; OUT=f"{REPO}/analysis/dkd_xenium_disease"
EPS,MINPTS=50,20; NEAR,FAR=50.0,200.0; BAFF="TNFSF13B"

cells=pd.read_parquet(f"{RE}/cells.parquet").reset_index(drop=True)
lig=pd.read_parquet(f"{OUT}/ligand_me_percell.parquet").reset_index(drop=True)
nc=pd.read_parquet(f"{OUT}/ncount_percell.parquet").reset_index(drop=True)
assert (cells.cell_id.astype(str).values==lig.cell_id.astype(str).values).all()
assert (cells.cell_id.astype(str).values==nc.cell_id.astype(str).values).all()
cells["BAFF"]=lig[BAFF].values; cells["nCount"]=nc.nCount.values
cells["is_Blin"]=cells.my_label.isin(["B","Plasma"])
GRP={"Stromal":cells.my_label.isin(["Fibroblast","VSMC","MC"]),
     "Epithelial (non-producer)":cells.my_lineage=="Epithelial",
     "Myeloid":cells.my_label=="Myeloid","Endothelial":cells.my_lineage=="Endothelial",
     "T":cells.my_label.isin(["CD4 T","CD8 T"]),"B-lineage":cells.is_Blin}
for k,m in GRP.items(): cells[f"g_{k}"]=m.values
AGG=["1006","HK2695","1007","1009"]   # the aggregate-bearing sections (2 B-rich DKD + C3GN + AA amyloid)

def dist_to_agg(g):
    B=g[g.is_Blin]
    if len(B)<MINPTS: return None
    lab=DBSCAN(eps=EPS,min_samples=MINPTS).fit(B[["spatial_x","spatial_y"]].values).labels_
    mem=B[["spatial_x","spatial_y"]].values[lab!=-1]
    if len(mem)==0: return None
    d,_=cKDTree(mem).query(g[["spatial_x","spatial_y"]].values,k=1)
    return d

# annotate every cell in aggregate-bearing sections with distance to nearest aggregate
cells["dist_agg"]=np.nan
for sid in AGG:
    idx=cells.index[cells.orig_ident==sid]; d=dist_to_agg(cells.loc[idx])
    if d is not None: cells.loc[idx,"dist_agg"]=d
AG=cells[cells.orig_ident.isin(AGG)].copy()
AG["zone"]=np.where(AG.dist_agg<=NEAR,"near",np.where(AG.dist_agg>FAR,"far","mid"))

def det(df,gmask_col,zone):
    s=df[(df[gmask_col])&(df.zone==zone)]; return (s.BAFF>0).mean(),len(s),int((s.BAFF>0).sum())

# ---------- STEP 1: local-ambient control (stromal vs epithelial in same window, pooled over AGG) ----------
rows=[]
for zone in ["near","far"]:
    for g in ["Stromal","Epithelial (non-producer)"]:
        r,n,pos=det(AG,f"g_{g}",zone); rows.append(dict(zone=zone,cell_type=g,baff_det_pct=round(r*100,2),n=n,pos=pos))
s1=pd.DataFrame(rows); s1.to_csv(f"{OUT}/baff_ambient_step1.csv",index=False)
sn=s1[(s1.zone=='near')&(s1.cell_type=='Stromal')].baff_det_pct.iloc[0]
en=s1[(s1.zone=='near')&(s1.cell_type=='Epithelial (non-producer)')].baff_det_pct.iloc[0]
print("=== STEP 1 — local-ambient control (peri-aggregate window, pooled) ===")
print(s1.to_string(index=False)); print(f"stromal-near / epithelial-near (local ambient) = {sn/max(en,1e-9):.2f}x")

# ---------- STEP 2: near vs far within stromal, per section ----------
rows=[]
for sid in AGG:
    g=AG[AG.orig_ident==sid]
    rn,nn,pn=det(g,"g_Stromal","near"); rf,nf,pf=det(g,"g_Stromal","far")
    rows.append(dict(orig_ident=sid,Condition=g.Condition.iloc[0],stromal_near_pct=round(rn*100,2),near_n=nn,
        stromal_far_pct=round(rf*100,2),far_n=nf,near_over_far=round(rn/max(rf,1e-9),2)))
# pooled
gn=AG; rn,nn,pn=det(gn,"g_Stromal","near"); rf,nf,pf=det(gn,"g_Stromal","far")
rows.append(dict(orig_ident="POOLED",Condition="-",stromal_near_pct=round(rn*100,2),near_n=nn,
    stromal_far_pct=round(rf*100,2),far_n=nf,near_over_far=round(rn/max(rf,1e-9),2)))
s2=pd.DataFrame(rows); s2.to_csv(f"{OUT}/baff_ambient_step2.csv",index=False)
print("\n=== STEP 2 — stromal BAFF near(<=50) vs far(>200), per section ===")
print(s2.to_string(index=False))

# ---------- STEP 3: density confound — near/far BAFF ratio per cell type + nCount ratio ----------
rows=[]
for g in GRP:
    rn,nn,_=det(AG,f"g_{g}","near"); rf,nf,_=det(AG,f"g_{g}","far")
    ncn=AG[(AG[f"g_{g}"])&(AG.zone=="near")].nCount.mean(); ncf=AG[(AG[f"g_{g}"])&(AG.zone=="far")].nCount.mean()
    rows.append(dict(cell_type=g,baff_near_pct=round(rn*100,2),baff_far_pct=round(rf*100,2),
        baff_near_over_far=round(rn/max(rf,1e-9),2),nCount_near=round(ncn,0),nCount_far=round(ncf,0),
        nCount_near_over_far=round(ncn/max(ncf,1e-9),2)))
s3=pd.DataFrame(rows); s3.to_csv(f"{OUT}/baff_ambient_step3.csv",index=False)
print("\n=== STEP 3 — near/far BAFF ratio per cell type vs transcript-density (nCount) ratio ===")
print(s3.to_string(index=False))

# ---------- verdict ----------
strom=s3[s3.cell_type=="Stromal"].iloc[0]; epi=s3[s3.cell_type=="Epithelial (non-producer)"].iloc[0]
print("\n=== VERDICT ===")
print(f"STEP1 stromal-near {sn:.2f}% vs epithelial-near(ambient) {en:.2f}%  -> {sn/max(en,1e-9):.1f}x")
print(f"STEP2 stromal near/far (pooled) = {s2[s2.orig_ident=='POOLED'].near_over_far.iloc[0]}x")
print(f"STEP3 stromal BAFF near/far {strom.baff_near_over_far}x vs its density near/far {strom.nCount_near_over_far}x;"
      f" epithelial(non-producer) BAFF near/far {epi.baff_near_over_far}x")

# ---------- figures ----------
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys; sys.path.insert(0,f"{REPO}/presentation/figures"); import figstyle as fs
FIG=f"{OUT}/figures"; os.makedirs(FIG,exist_ok=True)
def save(fig,name):
    fig.savefig(f"{FIG}/{name}.png",dpi=300,bbox_inches="tight"); fig.savefig(f"{FIG}/{name}.svg",bbox_inches="tight")
    plt.close(fig); print(f"  [ok] {name}")
fig,axes=plt.subplots(1,3,figsize=(17,5.4))
# (1) STEP1 local-ambient
ax=axes[0]; x=np.arange(2); w=0.36
strm=[s1[(s1.zone==z)&(s1.cell_type=="Stromal")].baff_det_pct.iloc[0] for z in ["near","far"]]
epin=[s1[(s1.zone==z)&(s1.cell_type=="Epithelial (non-producer)")].baff_det_pct.iloc[0] for z in ["near","far"]]
ax.bar(x-w/2,strm,w,label="Stromal (candidate producer)",color="#e377c2")
ax.bar(x+w/2,epin,w,label="Epithelial (non-producer ambient)",color="#c7c7c7")
ax.set_xticks(x); ax.set_xticklabels(["near (<=50um)","far (>200um)"]); ax.set_ylabel("BAFF detection (% cells)")
ax.set_title("(1) local-ambient control\nstromal > epithelial in BOTH zones = tissue-wide baseline,\nnot a peri-aggregate effect",fontsize=10)
ax.legend(frameon=False,fontsize=8.5)
for s in ["top","right"]: ax.spines[s].set_visible(False)
# (2) STEP2 per-sample near vs far stromal
ax=axes[1]; d=s2[s2.orig_ident!="POOLED"]; x=np.arange(len(d))
ax.bar(x-w/2,d.stromal_near_pct,w,label="near (<=50um)",color="#08519c")
ax.bar(x+w/2,d.stromal_far_pct,w,label="far (>200um)",color="#bbbbbb")
for xi,r in zip(x,d.near_over_far): ax.text(xi,max(d.stromal_near_pct.iloc[xi],d.stromal_far_pct.iloc[xi])+0.08,f"{r}x",ha="center",fontsize=8,fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(d.orig_ident,fontsize=9); ax.set_ylabel("stromal BAFF detection (%)")
ax.set_title("(2) stromal near vs far, per section\nstromal BAFF rises near aggregates (~2.4x)\n-- but is it stromal-specific? -> panel 3",fontsize=10)
ax.legend(frameon=False,fontsize=8.5)
for s in ["top","right"]: ax.spines[s].set_visible(False)
# (3) STEP3 decisive: near/far BAFF ratio per cell type + density reference
ax=axes[2]; order=["Stromal","Epithelial (non-producer)","Endothelial","Myeloid","T","B-lineage"]
d=s3.set_index("cell_type").reindex(order); x=np.arange(len(order))
cols=["#e377c2","#c7c7c7","#8c564b","#9467bd","#2ca02c","#1f77b4"]
ax.bar(x,d.baff_near_over_far,color=cols,edgecolor="#333",linewidth=0.4)
ax.axhline(1.0,ls="-",color="#888",lw=1)
ax.plot(x,d.nCount_near_over_far,"D",ms=7,color="black",label="transcript density (nCount) near/far")
ax.set_xticks(x); ax.set_xticklabels([o.replace(" (non-producer)","") for o in order],rotation=35,ha="right",fontsize=8.5)
ax.set_ylabel("BAFF near/far ratio"); ax.legend(frameon=False,fontsize=8.5,loc="upper right")
ax.set_title("(3) DECISIVE -- near/far BAFF rise is NON-SPECIFIC\nepithelial(non-prod) & endothelial rise like stromal;\ndensity flat -> local spillover field, not stromal production",fontsize=9.5)
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.suptitle("Peri-aggregate stromal BAFF stress-test -- VERDICT: ambient/spillover, NOT localized stromal production",fontsize=13,y=1.02)
fig.subplots_adjust(top=0.80,wspace=0.26)
save(fig,"baff_ambient_control")
print("== baff_ambient_control done ==")
