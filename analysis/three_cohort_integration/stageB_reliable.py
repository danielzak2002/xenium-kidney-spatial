#!/usr/bin/env python
"""STAGE 0b+1b — reliable-in-all-3 gene set (above per-cohort ambient via non-expressing-lineage
reference), then RE-TYPE + RE-INTEGRATE using ONLY those genes. Saves reliable parquet + figs."""
import os, warnings; warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp, harmonypy, json
import matplotlib.pyplot as plt
import sys; sys.path.insert(0,"/Users/danie/ClaudeCode/pilot_analyses/xenium/presentation/figures")
import figstyle as fs
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"; OBJ=f"{REPO}/outputs/objects"; OUT=f"{REPO}/analysis/three_cohort_integration"
rng=np.random.default_rng(0); COH=["RCC_big","RCC_figshare","DKD"]; COH_COL={"RCC_big":"#1F78B4","RCC_figshare":"#6BAED6","DKD":"#6A3D9A"}
MAXTHR=0.025; RATIO=1.5   # reliable: best-lineage detection >= 2.5% AND >= 1.5x cross-lineage median
def harmony(a):
    ho=harmonypy.run_harmony(a.obsm["X_pca"],a.obs,["cohort"],max_iter_harmony=20)
    Z=np.asarray(ho.Z_corr); a.obsm["X_pca_harmony"]=Z if Z.shape[0]==a.n_obs else Z.T
LIN={"B":["MS4A1","CD79A"],"Plasma":["MZB1","TNFRSF17","DERL3"],"T":["CD3E","CD3G","CD3D","TRBC2","TRAC"],
 "Myeloid":["CD68","CD14","CD163","LYZ","AIF1","ITGAX","C1QA"],"NK":["GNLY","KLRD1","NKG7","KLRF1"],
 "Endothelial":["PECAM1","VWF","EGFL7","CLDN5","CD34"],"Epithelial":["EPCAM","CDH1","KRT8","KRT18","PAX8","CA9"],
 "Stroma":["PDGFRA","PDGFRB","ACTA2","COL1A1","DCN","LUM"]}
TREG=["FOXP3","IL2RA","CTLA4"]; CYTO=["CD8A","GZMK","GZMB"]

print("loading + prior lineage (for gene QC) ...")
A=ad.read_h5ad(f"{OBJ}/three_cohort_123.h5ad"); genes=list(map(str,A.var_names)); gi={g:i for i,g in enumerate(genes)}
A.obs["lineage0"]=pd.read_parquet(f"{OUT}/cells_labeled.parquet")["lineage"].values
A.layers["counts"]=A.X.copy()

# ---- STAGE 0b: reliable set ----
rel_by={}; maxd_by={}
for c in COH:
    sub=A[A.obs.cohort.values==c]; lins=pd.unique(sub.obs.lineage0)
    det=np.array([np.asarray((sub[sub.obs.lineage0.values==l].X>0).mean(0)).ravel() for l in lins])
    maxd=det.max(0); medd=np.median(det,0); rel_by[c]=(maxd>=MAXTHR)&(maxd>=RATIO*medd); maxd_by[c]=maxd
allrel=rel_by["RCC_big"]&rel_by["RCC_figshare"]&rel_by["DKD"]
RG=[genes[i] for i in range(len(genes)) if allrel[i]]; RSET=set(RG)
open(f"{OUT}/reliable_genes.txt","w").write("\n".join(RG)+"\n")
print(f"\nRELIABLE-in-all-3: {len(RG)} genes (of 123)")
drops=[]
for i,g in enumerate(genes):
    if not allrel[i]:
        fails=[c for c in COH if not rel_by[c][gi[g]]]
        drops.append(dict(gene=g,fails_in=";".join(fails),**{f"maxdet_{c}":round(float(maxd_by[c][gi[g]]),4) for c in COH}))
pd.DataFrame(drops).to_csv(f"{OUT}/reliable_drops.csv",index=False)
print(f"dropped {len(drops)} genes; e.g. DKD-driven drops:",[d['gene'] for d in drops if 'DKD' in d['fails_in']][:12])

# ---- STAGE 1b: re-type on reliable markers only ----
relLIN={k:[g for g in v if g in RSET] for k,v in LIN.items()}
print("\nlineage marker availability on RELIABLE set:")
for k,v in relLIN.items(): print(f"  {k:11s} {v if v else '*** NO reliable markers ***'}")
relTREG=[g for g in TREG if g in RSET]; relCYTO=[g for g in CYTO if g in RSET]
print(f"  Treg gate (reliable): {relTREG}  | cytotoxic gate (reliable): {relCYTO}")
sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
for k,gs in relLIN.items():
    sc.tl.score_genes(A,gs,score_name=f"r_{k}") if gs else A.obs.__setitem__(f"r_{k}",-9)
S=A.obs[[f"r_{k}" for k in LIN]].values; A.obs["lineage"]=np.array(list(LIN))[S.argmax(1)]
cnt=A.layers["counts"]
def anypos(gs):
    gs=[g for g in gs if g in RSET];
    return np.asarray(cnt[:,[gi[g] for g in gs]].sum(1)).ravel()>0 if gs else np.zeros(A.n_obs,bool)
isT=A.obs["lineage"].values=="T"; fp=anypos(relTREG); cy=anypos(relCYTO)
A.obs["is_B"]=A.obs.lineage.values=="B"; A.obs["is_endo"]=A.obs.lineage.values=="Endothelial"
A.obs["is_Treg"]=isT&fp; A.obs["is_cyto"]=isT&cy&(~fp)
print("\nreliable-set lineage counts:\n",A.obs.lineage.value_counts())
print(f"T {isT.sum()}: Treg {A.obs.is_Treg.sum()}, cytotoxic {A.obs.is_cyto.sum()}")
o=A.obs[["cohort","sample","lineage","is_B","is_endo","is_Treg","is_cyto"]].copy()
o["x"]=A.obsm["spatial"][:,0]; o["y"]=A.obsm["spatial"][:,1]
o.to_parquet(f"{OUT}/cells_labeled_reliable.parquet"); print("saved cells_labeled_reliable.parquet")

# ---- re-validate (Harmony on reliable-set PCA, balanced subsample) ----
idx=np.sort(np.concatenate([rng.choice(np.where(A.obs.cohort.values==c)[0],min(60000,(A.obs.cohort.values==c).sum()),replace=False) for c in COH]))
R=A[idx,[g for g in genes if g in RSET]].copy(); R.obs["cohort"]=R.obs["cohort"].astype(str)
sc.pp.scale(R,max_value=10); sc.tl.pca(R,n_comps=min(30,R.n_vars-1))
rcc=R[R.obs.cohort.isin(["RCC_big","RCC_figshare"])].copy(); rcc.obs["cohort"]=rcc.obs["cohort"].astype(str)
harmony(rcc); sc.pp.neighbors(rcc,use_rep="X_pca_harmony"); sc.tl.umap(rcc)
harmony(R); sc.pp.neighbors(R,use_rep="X_pca_harmony"); sc.tl.umap(R)
fig,axes=plt.subplots(1,3,figsize=(17,5.2))
for c in COH:
    m=R.obs.cohort.values==c; axes[0].scatter(R.obsm["X_umap"][m,0],R.obsm["X_umap"][m,1],s=2,c=COH_COL[c],linewidths=0,rasterized=True,label=c)
axes[0].set_title(f"All three — Harmony UMAP ({len(RG)}-gene reliable set)"); axes[0].axis("off"); axes[0].legend(markerscale=4,frameon=False,fontsize=9)
cmap=plt.get_cmap("tab10")
for i,k in enumerate(LIN):
    m=R.obs.lineage.values==k; axes[1].scatter(R.obsm["X_umap"][m,0],R.obsm["X_umap"][m,1],s=2,c=[cmap(i)],linewidths=0,rasterized=True,label=k)
axes[1].set_title("by lineage (reliable-marker typing)"); axes[1].axis("off"); axes[1].legend(markerscale=4,frameon=False,fontsize=8,ncol=2)
for c in ["RCC_big","RCC_figshare"]:
    m=rcc.obs.cohort.values==c; axes[2].scatter(rcc.obsm["X_umap"][m,0],rcc.obsm["X_umap"][m,1],s=2,c=COH_COL[c],linewidths=0,rasterized=True,label=c)
axes[2].set_title("Two RCC cohorts merge (batch baseline)"); axes[2].axis("off"); axes[2].legend(markerscale=4,frameon=False,fontsize=9)
fig.suptitle(f"Stage 1b · re-integration on the {len(RG)}-gene reliable set",fontsize=15); fs.save_fig(fig,"INTrel_umap")

# reliable-set per-cohort dot-plot
DOT=[g for g in ["MS4A1","CD79A","MZB1","CD3E","CD3G","FOXP3","CTLA4","CD8A","GZMK","CD68","C1QA","KLRD1","PECAM1","EGFL7","EPCAM","PDGFRB","ACTA2"] if g in RSET]
rows=list(LIN); fig,axes=plt.subplots(1,3,figsize=(18,5.4),sharey=True)
for ax,c in zip(axes,COH):
    sub=A[A.obs.cohort.values==c]; det=np.zeros((len(rows),len(DOT))); mn=np.zeros_like(det)
    for i,r in enumerate(rows):
        sl=sub[sub.obs.lineage.values==r]
        if sl.n_obs<20: continue
        for j,g in enumerate(DOT):
            v=np.asarray(sl.layers["counts"][:,gi[g]].todense()).ravel(); det[i,j]=(v>0).mean(); mn[i,j]=np.log1p(v).mean()
    z=(mn-mn.mean(0))/(mn.std(0)+1e-9)
    for i in range(len(rows)):
        for j in range(len(DOT)): ax.scatter(j,i,s=6+det[i,j]*230,c=[z[i,j]],cmap="RdBu_r",vmin=-1.2,vmax=2,edgecolor="#999",linewidth=0.2)
    ax.set_xticks(range(len(DOT))); ax.set_xticklabels(DOT,rotation=45,ha="right",fontsize=8)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows,fontsize=9); ax.invert_yaxis(); ax.set_title(c,fontsize=12)
fig.suptitle(f"Stage 1b validation · reliable-marker typing faithful in every cohort ({len(RG)}-gene set)",fontsize=14)
fs.save_fig(fig,"INTrel_dotplot",tight=False)
print("== stageB reliable done ==")
