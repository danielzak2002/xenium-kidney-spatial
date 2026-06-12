#!/usr/bin/env python
"""STAGE C (revision) — tumor-only figshare re-integration + metrics + concordance.
1b re-integrate (Harmony, Leiden tuned ~20); 1c iLISI/cLISI pre vs post; 1d Readout A tumor-only +
difference test; 1e DKD author-concordance (ARI + confusion). Read-only raw."""
import os, warnings; warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp, harmonypy
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import adjusted_rand_score
from scipy.spatial import cKDTree
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"; INT=f"{REPO}/analysis/three_cohort_integration"; OBJ=f"{REPO}/outputs/objects"
rng=np.random.default_rng(0); COL=["RCC_big","RCC_figshare","DKD"]
def harmony(a):
    ho=harmonypy.run_harmony(a.obsm["X_pca"],a.obs,["cohort"],max_iter_harmony=20)
    Z=np.asarray(ho.Z_corr); a.obsm["X_pca_harmony"]=Z if Z.shape[0]==a.n_obs else Z.T
def lisi(emb,labels,k=90):                                  # inverse-Simpson LISI (vectorized)
    nn=NearestNeighbors(n_neighbors=k+1).fit(emb); _,idx=nn.kneighbors(emb)
    codes=pd.factorize(labels)[0]; nb=codes[idx[:,1:]]; ncat=int(codes.max()+1)
    p=np.stack([(nb==c).mean(1) for c in range(ncat)],1); return float((1.0/np.sum(p**2,1)).mean())

# ---- 1a result + tumor-only labels ----
cls=pd.read_csv(f"{INT}/figshare_tumor_classification.csv"); tumor=set(cls[cls.region=="tumor"]["sample"])
print("figshare tumor samples:",sorted(tumor))
df=pd.read_parquet(f"{INT}/cells_labeled_reliable.parquet")
keep=(df.cohort!="RCC_figshare")|(df["sample"].isin(tumor))
dft=df[keep].copy(); dft.to_parquet(f"{INT}/cells_labeled_reliable_tumoronly.parquet")
print("tumor-only cells per cohort:",dict(dft.groupby("cohort").size()))

# ---- 1b re-integrate (balanced subsample) + 1c metrics ----
A=ad.read_h5ad(f"{OBJ}/three_cohort_123.h5ad")
RG=set(open(f"{INT}/reliable_genes.txt").read().split())
lab=pd.read_parquet(f"{INT}/cells_labeled_reliable.parquet")
A.obs["lineage"]=lab["lineage"].values
# tumor-only mask on A (align by row order = same as df)
A.obs["_keep"]=keep.values
A=A[A.obs._keep.values][:, [g for g in map(str,A.var_names) if g in RG]].copy()
idx=np.sort(np.concatenate([rng.choice(np.where(A.obs.cohort.values==c)[0],min(50000,int((A.obs.cohort.values==c).sum())),replace=False) for c in COL]))
R=A[idx].copy(); R.obs["cohort"]=R.obs["cohort"].astype(str)
sc.pp.normalize_total(R,target_sum=1e4); sc.pp.log1p(R); sc.pp.scale(R,max_value=10); sc.tl.pca(R,n_comps=30)
sc.pp.neighbors(R,use_rep="X_pca"); sc.tl.umap(R); pre=R.obsm["X_umap"].copy()
harmony(R); sc.pp.neighbors(R,use_rep="X_pca_harmony"); sc.tl.umap(R); post=R.obsm["X_umap"].copy()
# tune Leiden resolution to ~15-25 clusters
res_used=None
for res in [0.05,0.08,0.12,0.2,0.3]:
    sc.tl.leiden(R,resolution=res,key_added="lt"); nc=R.obs['lt'].nunique()
    print(f"  leiden res={res}: {nc} clusters")
    if 15<=nc<=25: res_used=res; break
if res_used is None: res_used=0.12; sc.tl.leiden(R,resolution=res_used,key_added="lt")
print(f"chosen leiden resolution {res_used} -> {R.obs['lt'].nunique()} clusters")
# 1c metrics
nC=R.obs.cohort.nunique(); nL=R.obs.lineage.nunique()
m=dict(iLISI_pre=lisi(pre,R.obs.cohort.values),iLISI_post=lisi(post,R.obs.cohort.values),
       cLISI_pre=lisi(pre,R.obs.lineage.values),cLISI_post=lisi(post,R.obs.lineage.values),n_cohorts=nC,n_lineages=nL)
print("\n=== 1c integration metrics (LISI, inverse-Simpson) ===")
print(f"  iLISI (mixing; 1=segregated .. {nC}=perfectly mixed):  pre {m['iLISI_pre']:.2f} -> post {m['iLISI_post']:.2f}")
print(f"  cLISI (cell-type; 1=preserved .. {nL}=scrambled):       pre {m['cLISI_pre']:.2f} -> post {m['cLISI_post']:.2f}")
pd.DataFrame([m]).to_csv(f"{INT}/integration_metrics.csv",index=False)
fm=df[["is_B","is_Treg","is_cyto"]].reindex(R.obs_names)   # align immune flags by cell name
mod=pd.DataFrame(dict(cohort=R.obs.cohort.values,lineage=R.obs.lineage.values,leiden=R.obs['lt'].astype(str).values,
    is_B=fm["is_B"].values,is_Treg=fm["is_Treg"].values,is_cyto=fm["is_cyto"].values,
    pre_x=pre[:,0],pre_y=pre[:,1],post_x=post[:,0],post_y=post[:,1]))
mod.to_parquet(f"{INT}/model_tumoronly.parquet"); print("saved model_tumoronly.parquet")
pd.Series({"leiden_resolution":res_used,"n_clusters":int(R.obs['lt'].nunique())}).to_csv(f"{INT}/leiden_choice.csv")

# ---- 1d Readout A tumor-only + difference test ----
EPS=50.0;MINPTS=20;Rr=50.0
def aggs(cd):
    rows=[]
    for s,g in cd.groupby("sample"):
        xy=g[["x","y"]].values; isB=g.is_B.values
        if isB.sum()<MINPTS: continue
        cl=DBSCAN(eps=EPS,min_samples=MINPTS).fit(xy[isB]).labels_; mem0=np.where(isB)[0][cl!=-1]
        if len(mem0)==0: continue
        tree=cKDTree(xy); bgT=g.is_Treg.mean(); bgC=g.is_cyto.mean(); tr=g.is_Treg.values; cy=g.is_cyto.values
        for c in [c for c in np.unique(cl) if c!=-1]:
            mem=np.where(isB)[0][cl==c]; reg=np.unique(np.concatenate([np.asarray(t,int) for t in tree.query_ball_point(xy[mem],r=Rr)]))
            rows.append((len(reg),int(tr[reg].sum()),int(cy[reg].sum()),float(bgT),float(bgC),s))
    return pd.DataFrame(rows,columns=["n","Tin","Cin","bgT","bgC","section"])
def dl(a):
    ti,te,ci,ce=a.Tin.values,(a.n*a.bgT).values,a.Cin.values,(a.n*a.bgC).values
    return np.log2((ti.sum()+1e-9)/(te.sum()+1e-9))-np.log2((ci.sum()+1e-9)/(ce.sum()+1e-9))
def boot_arr(a,nb=5000,seed=0):
    r=np.random.default_rng(seed);n=len(a);ti,te,ci,ce=a.Tin.values,(a.n*a.bgT).values,a.Cin.values,(a.n*a.bgC).values
    return np.array([np.log2((ti[i].sum()+1e-9)/(te[i].sum()+1e-9))-np.log2((ci[i].sum()+1e-9)/(ce[i].sum()+1e-9)) for i in (r.integers(0,n,n) for _ in range(nb))])
AG={c:aggs(dft[dft.cohort==c]) for c in COL}; AG["RCC_pooled"]=pd.concat([AG["RCC_big"],AG["RCC_figshare"]],ignore_index=True)
print("\n=== 1d Readout A (tumor-only figshare) ===")
rowsA=[]
for c in ["RCC_big","RCC_figshare","RCC_pooled","DKD"]:
    a=AG[c]; d=dl(a); bs=boot_arr(a); lo,hi=np.percentile(bs,[2.5,97.5])
    rowsA.append(dict(cohort=c,n_agg=len(a),n_sec=a.section.nunique(),delta=round(d,3),lo=round(lo,3),hi=round(hi,3),fold=round(2**d,2)))
    print(f"  {c:13s}: {len(a)} aggs/{a.section.nunique()} sec | Δ {d:+.2f} [{lo:+.2f},{hi:+.2f}]")
pd.DataFrame(rowsA).to_csv(f"{INT}/readoutA_tumoronly.csv",index=False)
print("difference (cohort − DKD):")
rowsD=[]
for ca in ["RCC_pooled","RCC_big","RCC_figshare"]:
    a=AG[ca];b=AG["DKD"];r=np.random.default_rng(1);na=len(a);nb_=len(b)
    pt=dl(a)-dl(b)
    bs=[dl(a.iloc[r.integers(0,na,na)])-dl(b.iloc[r.integers(0,nb_,nb_)]) for _ in range(5000)]
    lo,hi=np.percentile(bs,[2.5,97.5]); excl=(lo>0)or(hi<0)
    rowsD.append(dict(contrast=f"{ca} − DKD",diff=round(pt,3),lo=round(lo,3),hi=round(hi,3),excludes_zero=excl))
    print(f"  {ca:13s} − DKD: {pt:+.2f} [{lo:+.2f},{hi:+.2f}]  excludes 0: {excl}")
pd.DataFrame(rowsD).to_csv(f"{INT}/readoutA_tumoronly_difference.csv",index=False)

# ---- 1e DKD author-concordance ----
print("\n=== 1e DKD author-concordance ===")
d2=ad.read_h5ad(f"{OBJ}/three_cohort_123.h5ad",backed="r")  # for DKD obs_names order
dkd_names=[n for n in d2.obs_names if str(n).startswith("DKD_")]; d2.file.close()
dem=ad.read_h5ad(f"{REPO}/Demoulin26/data/spatial_adata_xenium_cosmx_zenodo.h5ad",backed="r")
xen=dem.obs["tech"].astype(str).values=="Xenium"; gi=np.where(xen)[0]
ann=dem.obs["annotation_updated"].astype(str).values; imm=dem.obs["immune_cell_annotation_combined"].astype(str).values; dem.file.close()
our=df[df.cohort=="DKD"].copy(); jj=np.array([int(n.split("_")[1]) for n in our.index])
our["author_ann"]=ann[gi[jj]]; our["author_imm"]=imm[gi[jj]]
# map authors coarse lineage to ours
A2L={"PT":"Epithelial","iPT":"Epithelial","TAL":"Epithelial","iTAL":"Epithelial","DCT":"Epithelial","CNT":"Epithelial","PC":"Epithelial","IC A":"Epithelial","IC B":"Epithelial","DTL_ATL":"Epithelial","PEC":"Epithelial","Podo":"Podocyte","EC_Peritub":"Endothelial","EC_glom":"Endothelial","EC_DVR":"Endothelial","EC_Lymph":"Endothelial","Fibroblast":"Stroma","VSMC":"Stroma","MC1":"Stroma","Immune":"Immune"}
I2L={"B":"B","Plasma":"Plasma","Macro":"Myeloid","cDC":"Myeloid","pDC":"Myeloid","CD4+":"T","CD8+":"T","NK":"NK","Neutrophil":"Myeloid","Baso_Mast":"Myeloid"}
def authlin(a,i): return I2L.get(i, A2L.get(a,"other")) if i!="Unknown" else A2L.get(a,"other")
our["author_lineage"]=[authlin(a,i) for a,i in zip(our.author_ann,our.author_imm)]
common=[l for l in ["B","Plasma","T","Myeloid","NK","Endothelial","Epithelial","Stroma","Podocyte"]]
ari=adjusted_rand_score(our.author_lineage,our.lineage)
agree=float((our.author_lineage==our.lineage).mean())
print(f"  ARI(our lineage vs author lineage) = {ari:.3f} | overall agreement = {agree:.3f}")
conf=pd.crosstab(our.author_lineage,our.lineage).reindex(index=common,columns=common).fillna(0).astype(int)
conf.to_csv(f"{INT}/dkd_concordance_confusion.csv")
pd.Series({"ARI":round(ari,3),"agreement":round(agree,3),"n_cells":len(our)}).to_csv(f"{INT}/dkd_concordance_summary.csv")
print(conf.to_string())
print("== stageC done ==")
