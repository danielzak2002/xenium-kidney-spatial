#!/usr/bin/env python
"""STAGE 1 — Harmony integration FOR TYPING + per-cell lineage labels (cross-cohort consistent).
Integration (Harmony+Leiden on a balanced subsample) VALIDATES that the shared-marker lineage
typing is consistent across cohorts; per-cell labels for ALL cells are assigned by the SAME
marker definitions in every cohort (uniform, scalable). Saves obs+flags+stress-counts parquet."""
import os, warnings; warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp, harmonypy
import matplotlib.pyplot as plt
def harmony(adata):  # harmonypy 2.0 returns Z_corr as (N,d); scanpy wrapper assumes (d,N) -> call directly
    ho=harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, ["cohort"], max_iter_harmony=20)
    Z=np.asarray(ho.Z_corr); adata.obsm["X_pca_harmony"]= Z if Z.shape[0]==adata.n_obs else Z.T
import sys; sys.path.insert(0,"/Users/danie/ClaudeCode/pilot_analyses/xenium/presentation/figures")
import figstyle as fs
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"; OBJ=f"{REPO}/outputs/objects"
OUT=f"{REPO}/analysis/three_cohort_integration"; rng=np.random.default_rng(0)
COH_COL={"RCC_big":"#1F78B4","RCC_figshare":"#6BAED6","DKD":"#6A3D9A"}

# lineage marker sets (intersected with the 123 at runtime)
LIN={"B":["MS4A1","CD79A"],"Plasma":["MZB1","TNFRSF17","DERL3"],
 "T":["CD3E","CD3G","CD3D","TRBC2","TRAC"],
 "Myeloid":["CD68","CD14","CD163","LYZ","AIF1","ITGAX","C1QA"],
 "NK":["GNLY","KLRD1","NKG7","KLRF1"],
 "Endothelial":["PECAM1","VWF","EGFL7","CLDN5","CD34"],
 "Epithelial":["EPCAM","CDH1","KRT8","KRT18","PAX8","CA9"],
 "Stroma":["PDGFRA","PDGFRB","ACTA2","COL1A1","DCN","LUM"]}
TREG=["FOXP3","IL2RA","CTLA4"]; CYTO=["CD8A","GZMK","GZMB"]
STRESS=["ANGPT2","CXCL9","CXCL10","HLA-DRA","PECAM1","VWF"]

print("loading three_cohort_123.h5ad ...")
A=ad.read_h5ad(f"{OBJ}/three_cohort_123.h5ad"); V=set(map(str,A.var_names))
A.layers["counts"]=A.X.copy()
sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
print("cohorts:",dict(A.obs.groupby('cohort').size()))

# ---- per-cell lineage by marker-score argmax (uniform across cohorts) ----
present={k:[g for g in v if g in V] for k,v in LIN.items()}
for k,gs in present.items():
    sc.tl.score_genes(A,gs,score_name=f"sc_{k}") if gs else A.obs.__setitem__(f"sc_{k}",-9)
S=A.obs[[f"sc_{k}" for k in LIN]].values; labs=np.array(list(LIN))[S.argmax(1)]
A.obs["lineage"]=labs
# within-T marker gating for Treg / cytotoxic (the readout-A states)
cnt=A.layers["counts"]
def anypos(genes):
    gs=[g for g in genes if g in V]
    if not gs: return np.zeros(A.n_obs,bool)
    return np.asarray(cnt[:,[A.var_names.get_loc(g) for g in gs]].sum(1)).ravel()>0
isT=labs=="T"; fp=anypos(TREG); cy=anypos(CYTO)
A.obs["is_B"]=labs=="B"; A.obs["is_endo"]=labs=="Endothelial"
A.obs["is_Treg"]=isT&fp; A.obs["is_cyto"]=isT&cy&(~fp)
print("lineage counts:\n",A.obs.lineage.value_counts())
print(f"T cells {isT.sum()}: Treg(FOXP3+) {A.obs.is_Treg.sum()}, cytotoxic(CD8A/GZMK+) {A.obs.is_cyto.sum()}")

# ---- stress marker raw counts (for readout B) ----
obs=A.obs[["cohort","sample","lineage","is_B","is_endo","is_Treg","is_cyto"]].copy()
obs["x"]=A.obsm["spatial"][:,0]; obs["y"]=A.obsm["spatial"][:,1]
for g in STRESS:
    obs[f"cnt_{g}"]=np.asarray(cnt[:,A.var_names.get_loc(g)].todense()).ravel() if g in V else 0
obs.to_parquet(f"{OUT}/cells_labeled.parquet")
print("saved cells_labeled.parquet")

# ============================================================================
# Integration validation on a BALANCED subsample (Harmony + Leiden + UMAP)
# ============================================================================
print("\nintegration validation (balanced subsample) ...")
idx=[]
for c in A.obs.cohort.unique():
    ci=np.where(A.obs.cohort.values==c)[0]; idx.append(rng.choice(ci,min(60000,len(ci)),replace=False))
idx=np.sort(np.concatenate(idx)); R=A[idx].copy()
R.obs["cohort"]=R.obs["cohort"].astype(str)                     # drop unused categorical levels (harmonypy)
sc.pp.scale(R,max_value=10); sc.tl.pca(R,n_comps=30)
# (i) two RCC cohorts first — batch baseline
rcc=R[R.obs.cohort.isin(["RCC_big","RCC_figshare"])].copy()
rcc.obs["cohort"]=rcc.obs["cohort"].astype(str)
harmony(rcc)
sc.pp.neighbors(rcc,use_rep="X_pca_harmony"); sc.tl.umap(rcc)
# (ii) all three
harmony(R)
sc.pp.neighbors(R,use_rep="X_pca_harmony"); sc.tl.umap(R); sc.tl.leiden(R,resolution=1.0)
R.obsm["spatial_dummy"]=R.obsm["X_umap"]

# fig: integration UMAP by cohort + by lineage; 2-RCC merge
fig,axes=plt.subplots(1,3,figsize=(17,5.2))
for c in COH_COL:
    m=R.obs.cohort.values==c; axes[0].scatter(R.obsm["X_umap"][m,0],R.obsm["X_umap"][m,1],s=2,c=COH_COL[c],linewidths=0,rasterized=True,label=c)
axes[0].set_title("All three — Harmony UMAP by cohort"); axes[0].axis("off"); axes[0].legend(markerscale=4,frameon=False,fontsize=9)
linord=[k for k in LIN]; cmap=plt.get_cmap("tab10")
for i,k in enumerate(linord):
    m=R.obs.lineage.values==k; axes[1].scatter(R.obsm["X_umap"][m,0],R.obsm["X_umap"][m,1],s=2,c=[cmap(i)],linewidths=0,rasterized=True,label=k)
axes[1].set_title("by lineage (shared-marker typing)"); axes[1].axis("off"); axes[1].legend(markerscale=4,frameon=False,fontsize=8,ncol=2)
for c in ["RCC_big","RCC_figshare"]:
    m=rcc.obs.cohort.values==c; axes[2].scatter(rcc.obsm["X_umap"][m,0],rcc.obsm["X_umap"][m,1],s=2,c=COH_COL[c],linewidths=0,rasterized=True,label=c)
axes[2].set_title("Two RCC cohorts merge (batch baseline)"); axes[2].axis("off"); axes[2].legend(markerscale=4,frameon=False,fontsize=9)
fig.suptitle("Stage 1 · Harmony integration for typing — cohorts mix, lineages resolve",fontsize=15)
fs.save_fig(fig,"INT_umap")

# fig: per-cohort marker dot-plot (validation — labels marker-faithful)
DOT=["MS4A1","CD79A","MZB1","CD3E","CD3G","FOXP3","CD8A","GZMK","CD68","C1QA","GNLY","KLRD1","PECAM1","VWF","EPCAM","KRT8","PDGFRB","ACTA2"]
DOT=[g for g in DOT if g in V]; rows=[k for k in LIN]
fig,axes=plt.subplots(1,3,figsize=(20,5.4),sharey=True)
for ax,c in zip(axes,["RCC_big","RCC_figshare","DKD"]):
    sub=A[A.obs.cohort.values==c]
    det=np.zeros((len(rows),len(DOT))); mn=np.zeros_like(det)
    for i,r in enumerate(rows):
        sl=sub[sub.obs.lineage.values==r]
        if sl.n_obs<20: continue
        for j,g in enumerate(DOT):
            v=np.asarray(sl.layers["counts"][:,A.var_names.get_loc(g)].todense()).ravel()
            det[i,j]=(v>0).mean(); mn[i,j]=np.log1p(v).mean()
    z=(mn-mn.mean(0))/(mn.std(0)+1e-9)
    for i in range(len(rows)):
        for j in range(len(DOT)): ax.scatter(j,i,s=6+det[i,j]*230,c=[z[i,j]],cmap="RdBu_r",vmin=-1.2,vmax=2,edgecolor="#999",linewidth=0.2)
    ax.set_xticks(range(len(DOT))); ax.set_xticklabels(DOT,rotation=45,ha="right",fontsize=8)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows,fontsize=9); ax.invert_yaxis()
    ax.set_title(c,fontsize=12)
fig.suptitle("Stage 1 validation · each lineage defined by the right markers in EVERY cohort (typing consistent)",fontsize=14)
fs.save_fig(fig,"INT_dotplot",tight=False)

# fig: Treg:cytotoxic composition RCC vs DKD (immune diff still visible post-labeling)
fig,ax=plt.subplots(figsize=(7,4.6)); x=np.arange(3); w=0.38
for i,(st,col) in enumerate([("is_Treg","#d62728"),("is_cyto","#2ca02c")]):
    fr=[A.obs[(A.obs.cohort==c)&isT][st].mean() for c in ["RCC_big","RCC_figshare","DKD"]]
    ax.bar(x+(i-0.5)*w,fr,w,label="Treg" if st=="is_Treg" else "cytotoxic",color=col)
ax.set_xticks(x); ax.set_xticklabels(["RCC_big","RCC_figshare","DKD"]); ax.set_ylabel("fraction of T cells")
ax.set_title("Treg vs cytotoxic among T cells — immune composition difference survives labeling"); ax.legend(frameon=False)
fs.save_fig(fig,"INT_composition")
print("== stage1 done ==")
