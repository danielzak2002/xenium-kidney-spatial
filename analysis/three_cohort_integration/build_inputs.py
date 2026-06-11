#!/usr/bin/env python
"""STAGE 0 — build three-cohort inputs on the 123-gene shared space (raw counts + centroids).
figshare transcript tables -> per-cell counts (native cell_id), streamed one file at a time.
Read-only raw; writes outputs/objects/three_cohort_123.h5ad (git-ignored)."""
import os, json, gc, warnings, subprocess
warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, scanpy as sc
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
OUT=f"{REPO}/analysis/three_cohort_integration"; OBJ=f"{REPO}/outputs/objects"
GENES=[l.strip() for l in open(f"{REPO}/analysis/three_cohort_assessment/three_way_shared_genes.txt") if l.strip()]
GSET=set(GENES); GIDX={g:i for i,g in enumerate(GENES)}
print(f"123-gene shared space: {len(GENES)} genes")

# ---------- figshare: stream transcript table -> cell x gene (123) + centroids ----------
def build_figshare_sample(url, sid, chunk=3_000_000):
    gz=f"/tmp/fig_{sid}.csv.gz"
    print(f"  [{sid}] downloading {url} ...", flush=True)
    subprocess.run(["curl","-sL","--max-time","1200",url,"-o",gz],check=True)
    cnt_parts=[]; cen_parts=[]   # per-chunk groupby results
    for ch in pd.read_csv(gz, usecols=["cell_id","feature_name","x_location","y_location"],
                          chunksize=chunk, low_memory=False):
        ch=ch[ch.cell_id>0]                                   # assigned transcripts only
        cen=ch.groupby("cell_id").agg(xs=("x_location","sum"),ys=("y_location","sum"),n=("x_location","size"))
        cen_parts.append(cen)
        g=ch[ch.feature_name.isin(GSET)]
        if len(g): cnt_parts.append(g.groupby(["cell_id","feature_name"]).size().rename("c").reset_index())
    os.remove(gz)
    # centroids
    cen=pd.concat(cen_parts).groupby(level=0).sum(); cen["x"]=cen["xs"]/cen["n"]; cen["y"]=cen["ys"]/cen["n"]
    cells=cen.index.values
    cpos={c:i for i,c in enumerate(cells)}
    # counts -> sparse cells x 123
    cc=pd.concat(cnt_parts).groupby(["cell_id","feature_name"],as_index=False)["c"].sum()
    rows=np.array([cpos[c] for c in cc.cell_id]); cols=np.array([GIDX[g] for g in cc.feature_name])
    X=sp.coo_matrix((cc.c.values,(rows,cols)),shape=(len(cells),len(GENES))).tocsr()
    a=ad.AnnData(X=X.astype(np.float32))
    a.var_names=GENES; a.obs_names=[f"{sid}_{c}" for c in cells]
    a.obs["cohort"]="RCC_figshare"; a.obs["sample"]=sid
    a.obsm["spatial"]=cen[["x","y"]].values.astype(float)
    a.obs["n_counts"]=np.asarray(X.sum(1)).ravel()
    a=a[a.obs["n_counts"]>=5].copy()                          # drop near-empty cells
    print(f"  [{sid}] {a.n_obs} cells x {a.n_vars} genes (median counts {np.median(a.obs.n_counts):.0f})",flush=True)
    return a

def stage_figshare():
    fs=sorted(json.load(open("/tmp/figshare_25685961.json"))["files"], key=lambda x:x["size"])
    parts=[]
    for k,f in enumerate(fs,1):
        parts.append(build_figshare_sample(f["download_url"], f"figS{k:02d}"))
        gc.collect()
    A=ad.concat(parts, join="outer"); A.var_names=GENES
    return A

# ---------- BIG RCC ----------
def stage_big():
    print("BIG RCC: load cell_feature_matrix.h5 ...",flush=True)
    a=sc.read_10x_h5(f"{REPO}/kidney_10x/data/cell_feature_matrix.h5")
    a.var_names_make_unique()
    keep=[g for g in GENES if g in set(map(str,a.var_names))]
    a=a[:, keep].copy()
    # reorder/pad to full 123
    full=ad.AnnData(sp.csr_matrix((a.n_obs,len(GENES)),dtype=np.float32)); full.var_names=GENES; full.obs_names=a.obs_names
    full[:, keep]=a.X
    cells=pd.read_parquet(f"{REPO}/kidney_10x/data/cells.parquet")
    cmap=cells.set_index(cells.cell_id.astype(str))
    idx=full.obs_names.astype(str)
    full.obsm["spatial"]=cmap.reindex(idx)[["x_centroid","y_centroid"]].values.astype(float)
    full.obs["cohort"]="RCC_big"; full.obs["sample"]="RCC_big"
    full.obs["n_counts"]=np.asarray(full.X.sum(1)).ravel()
    full=full[(full.obs.n_counts>=5)&np.isfinite(full.obsm["spatial"]).all(1)].copy()
    print(f"BIG RCC: {full.n_obs} cells; genes present {len(keep)}/123",flush=True)
    return full

# ---------- DKD ----------
def stage_dkd():
    print("DKD: load Xenium subset ...",flush=True)
    a=ad.read_h5ad(f"{REPO}/Demoulin26/data/spatial_adata_xenium_cosmx_zenodo.h5ad", backed="r")
    xen=a.obs["tech"].astype(str).values=="Xenium"
    present=[g for g in GENES if g in set(map(str,a.var_names))]
    idx=np.where(xen)[0]
    sub=a[idx, present]; M=sub.to_memory()
    samp=a.obs["orig_ident"].astype(str).values[idx]; xy=np.asarray(a.obsm["spatial"],float)[idx]
    C=(M.layers["counts"] if "counts" in M.layers else M.X); C=C.tocsr() if sp.issparse(C) else sp.csr_matrix(C)
    a.file.close()
    # struct-zero check on the 123 present genes
    nz=np.asarray((C>0).sum(0)).ravel(); sz=[present[i] for i in range(len(present)) if nz[i]==0]
    keep=[present[i] for i in range(len(present)) if nz[i]>0]
    print(f"DKD: {len(present)}/123 in var; structural-zero on Xenium dropped: {sz}",flush=True)
    full=ad.AnnData(sp.csr_matrix((M.n_obs,len(GENES)),dtype=np.float32)); full.var_names=GENES; full.obs_names=[f"DKD_{i}" for i in range(M.n_obs)]
    ci={g:i for i,g in enumerate(present)}
    full[:, keep]=C[:, [ci[g] for g in keep]]
    full.obsm["spatial"]=xy; full.obs["cohort"]="DKD"; full.obs["sample"]=samp
    full.obs["n_counts"]=np.asarray(full.X.sum(1)).ravel()
    full=full[full.obs.n_counts>=5].copy()
    print(f"DKD: {full.n_obs} cells; struct-zeros={sz}",flush=True)
    return full, sz

if __name__=="__main__":
    fig=stage_figshare()
    big=stage_big()
    dkd,dkd_sz=stage_dkd()
    A=ad.concat([fig,big,dkd], join="outer", index_unique=None); A.var_names=GENES
    A.obs["cohort"]=A.obs["cohort"].astype(str); A.obs["sample"]=A.obs["sample"].astype(str)
    A.uns["dkd_structzero_dropped"]=dkd_sz
    os.makedirs(OBJ,exist_ok=True)
    A.write(f"{OBJ}/three_cohort_123.h5ad")
    print("\n=== combined ===")
    print(A.obs.groupby("cohort").size())
    print("saved outputs/objects/three_cohort_123.h5ad")
    print("== stage0 done ==")
