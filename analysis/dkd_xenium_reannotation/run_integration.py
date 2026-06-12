#!/usr/bin/env python
"""STEP 0-1 — DKD Xenium-only re-integration, BLIND to author labels.
Subset 16 Xenium samples -> drop Xenium structural-zero genes -> normalize/log1p ->
HVG -> PCA -> Harmony(orig_ident) -> neighbors -> Leiden -> UMAP. Then subcluster the
immune compartment (resolves B/Plasma/CD4/CD8/Treg/Myeloid/NK/DC for validation vs the
authors' immune_cell_annotation_combined). Author labels are carried ALONG as columns but
NEVER used to guide clustering. Saves cells.parquet + per-cluster marker means + counts.
Read-only raw; outputs (parquet/h5ad) are git-ignored under analysis/."""
import os, warnings, time
warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad, scanpy as sc
sc.settings.verbosity=1
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
OBJ=f"{REPO}/Demoulin26/data/spatial_adata_xenium_cosmx_zenodo.h5ad"
OUT=f"{REPO}/analysis/dkd_xenium_reannotation"; os.makedirs(OUT,exist_ok=True)
SEED=0; np.random.seed(SEED)
AUTH=["annotation_updated","immune_cell_annotation_combined"]   # ground truth, set ASIDE
t0=time.time(); log=lambda m: print(f"[{time.time()-t0:6.0f}s] {m}",flush=True)

# ---------- STEP 0: subset to the 16 Xenium samples ----------
log("loading backed handle ...")
A=ad.read_h5ad(OBJ,backed="r")
mask=(A.obs["tech"].values=="Xenium")
log(f"materializing Xenium subset: {int(mask.sum()):,} cells ...")
adata=A[mask].to_memory(); del A
adata.X=adata.layers["counts"].copy()        # independent renormalization from raw counts
del adata.layers["counts"]
adata.obs["orig_ident"]=adata.obs["orig_ident"].astype(str)
# drop unused samples / empty categories
keep_cols=["cell_id","orig_ident","tech","nCount_RNA","nFeature_RNA","Condition","Age","Sex"]+AUTH
adata.obs=adata.obs[keep_cols].copy()
percell=adata.obs.groupby("orig_ident",observed=True).agg(
    n_cells=("cell_id","size"),median_counts=("nCount_RNA","median"),
    median_genes=("nFeature_RNA","median"),Condition=("Condition","first")).sort_index()
percell.to_csv(f"{OUT}/per_sample_counts.csv")
log(f"samples={percell.shape[0]} total={percell.n_cells.sum():,} genes(panel)={adata.n_vars}")

# ---------- drop Xenium structural-zero genes (CosMx-only in the union panel) ----------
gsum=np.asarray(adata.X.sum(0)).ravel(); measured=gsum>0
ndrop=int((~measured).sum())
adata=adata[:,measured].copy()
log(f"dropped {ndrop} structural-zero genes -> {adata.n_vars} measured on Xenium")

# ---------- STEP 1: independent integration ----------
sc.pp.normalize_total(adata,target_sum=1e4); sc.pp.log1p(adata)
adata.raw=adata                                # keep log-norm full panel for marker means
sc.pp.highly_variable_genes(adata,n_top_genes=2000,flavor="seurat")
log(f"HVG={int(adata.var.highly_variable.sum())}")
hv=adata[:,adata.var.highly_variable].copy()
sc.pp.scale(hv,max_value=10)
sc.tl.pca(hv,n_comps=50,svd_solver="arpack",random_state=SEED)
adata.obsm["X_pca"]=hv.obsm["X_pca"]; del hv
log("Harmony on orig_ident (harmonypy direct; 2.0 returns Z_corr as (N,d)) ...")
import harmonypy
ho=harmonypy.run_harmony(adata.obsm["X_pca"],adata.obs,["orig_ident"],max_iter_harmony=30,random_state=SEED)
Z=np.asarray(ho.Z_corr); adata.obsm["X_pca_harmony"]=Z if Z.shape[0]==adata.n_obs else Z.T
assert adata.obsm["X_pca_harmony"].shape==adata.obsm["X_pca"].shape, "harmony orientation"
np.save(f"{OUT}/X_pca_harmony.npy",adata.obsm["X_pca_harmony"])   # save embedding for cheap re-clustering
sc.pp.neighbors(adata,n_neighbors=15,use_rep="X_pca_harmony",random_state=SEED)
log("Leiden (global, res=1.5) ...")
sc.tl.leiden(adata,resolution=1.5,key_added="leiden",flavor="igraph",n_iterations=2,directed=False,random_state=SEED)
sc.tl.umap(adata,random_state=SEED)
adata.obs["leiden"]=adata.obs["leiden"].astype(str)
log(f"global clusters={adata.obs.leiden.nunique()}")

# ---------- per-(global cluster) mean log-norm expression (for marker scoring) ----------
def cluster_means(ad_obj,key):
    Xr=ad_obj.raw.X; genes=list(map(str,ad_obj.raw.var_names))
    cats=sorted(ad_obj.obs[key].unique(),key=lambda v:(len(v),v))
    rows={}
    for c in cats:
        m=(ad_obj.obs[key].values==c)
        rows[c]=np.asarray(Xr[m].mean(0)).ravel()
    return pd.DataFrame(rows,index=genes).T   # clusters x genes
cluster_means(adata,"leiden").to_csv(f"{OUT}/global_marker_means.csv")
log("saved global_marker_means.csv")

# ---------- immune compartment: provisional gate (markers only, NOT author labels) ----------
IMM_MARKERS=["PTPRC","CD3E","CD3D","CD2","CD8A","CD4","CD68","LYZ","ITGAM","C1QA",
             "MS4A1","CD79A","BANK1","MZB1","DERL3","NKG7","KLRD1","GNLY","FCER1A","LILRA4"]
imm=[g for g in IMM_MARKERS if g in set(map(str,adata.raw.var_names))]
sc.tl.score_genes(adata,imm,score_name="immune_score",use_raw=True)
gm=adata.obs.groupby("leiden",observed=True)["immune_score"].mean()
thr=gm.mean()+0.25*gm.std()                    # clusters clearly above-mean immune signal
imm_clusters=sorted(gm[gm>thr].index.tolist(),key=lambda v:(len(v),v))
adata.obs["is_immune_cluster"]=adata.obs["leiden"].isin(imm_clusters)
log(f"immune-gate clusters ({len(imm_clusters)}): {imm_clusters}  thr={thr:.3f}")

# ---------- subcluster immune cells (resolve subtypes) ----------
sub=adata[adata.obs.is_immune_cluster.values].copy()
log(f"immune subset: {sub.n_obs:,} cells")
sc.pp.neighbors(sub,n_neighbors=15,use_rep="X_pca_harmony",random_state=SEED)
sc.tl.leiden(sub,resolution=1.4,key_added="immune_leiden",flavor="igraph",n_iterations=2,directed=False,random_state=SEED)
sc.tl.umap(sub,random_state=SEED)
sub.obs["immune_leiden"]=sub.obs["immune_leiden"].astype(str)
log(f"immune clusters={sub.obs.immune_leiden.nunique()}")
cluster_means(sub,"immune_leiden").to_csv(f"{OUT}/immune_marker_means.csv")
log("saved immune_marker_means.csv")

# ---------- assemble per-cell table ----------
sp_xy=adata.obsm["spatial"]; um=adata.obsm["X_umap"]
cells=pd.DataFrame(dict(
    cell_id=adata.obs.cell_id.values, orig_ident=adata.obs.orig_ident.values,
    Condition=adata.obs.Condition.astype(str).values,
    leiden=adata.obs.leiden.values, is_immune_cluster=adata.obs.is_immune_cluster.values,
    immune_score=adata.obs.immune_score.values,
    umap_x=um[:,0], umap_y=um[:,1], spatial_x=sp_xy[:,0], spatial_y=sp_xy[:,1],
    author_annotation=adata.obs.annotation_updated.astype(str).values,
    author_immune=adata.obs.immune_cell_annotation_combined.astype(str).values),
    index=np.arange(adata.n_obs))
cells["immune_leiden"]=np.nan
iidx=np.where(adata.obs.is_immune_cluster.values)[0]
cells.loc[iidx,"immune_leiden"]=sub.obs["immune_leiden"].values
iu=sub.obsm["X_umap"]
cells["immune_umap_x"]=np.nan; cells["immune_umap_y"]=np.nan
cells.loc[iidx,"immune_umap_x"]=iu[:,0]; cells.loc[iidx,"immune_umap_y"]=iu[:,1]
cells.to_parquet(f"{OUT}/cells.parquet")
log(f"saved cells.parquet ({len(cells):,} rows)")
print("== run_integration done ==",flush=True)
