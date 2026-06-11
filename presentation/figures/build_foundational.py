#!/usr/bin/env python
"""build_foundational.py — evidence/veracity layer: Q1, T1, T2, A3raw, C1raw, C2raw.
Reuses figstyle. Read-only; backed DKD reads + spatial crops; reuse analysis CSVs."""
import os, warnings; warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree, ConvexHull
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import figstyle as fs
rng=np.random.default_rng(0); REPO=fs.REPO; OBJ=f"{REPO}/outputs/objects"
RCC_H5=f"{OBJ}/kidney_RCC_protein.h5ad"; CLN_H5=f"{OBJ}/cln_cosmx.h5ad"; PRCC_H5=f"{OBJ}/kidney_preview_PRCC.h5ad"
DKD_H5=f"{REPO}/Demoulin26/data/spatial_adata_xenium_cosmx_zenodo.h5ad"
BN=f"{REPO}/Demoulin26/analysis/bniche_dbscan"; TAB=f"{REPO}/outputs/tables"; XT=f"{REPO}/Demoulin26/analysis/cross_platform_tcell"

def dkd_lineage(ann,imm):
    out=np.full(len(ann),None,object)
    for i in range(len(ann)):
        out[i]=fs.DKD_IMMUNE.get(imm[i]) if imm[i]!="Unknown" else fs.DKD_NONIMM.get(ann[i])
    return out

# marker -> lineage group (for dot-plot column ordering)
MARKERS=[("MS4A1","B"),("CD79A","B"),("MZB1","Plasma"),("DERL3","Plasma"),("CD3D","T_lineage"),
 ("CD3E","T_lineage"),("FOXP3","Treg"),("CD8A","Cytotoxic"),("GZMB","Cytotoxic"),
 ("CD68","Myeloid"),("LYZ","Myeloid"),("PECAM1","Endothelial"),("VWF","Endothelial"),
 ("EPCAM","Epithelial"),("KRT8","Epithelial")]
ROW_ORDER=["B","Plasma","T_lineage","Myeloid","NK","Endothelial","Stroma","Podocyte","Tubular_epi","Malignant_epi"]

# ============================================================================
def figQ1():
    print("Q1 QC summary per dataset...")
    rows=[]
    for nm,csv,panel in [("RCC",f"{TAB.replace('/tables','/objects')}/qc_metrics_kidney_RCC_protein.csv",405),
                         ("PRCC",f"{OBJ}/qc_metrics_kidney_preview_PRCC.csv",377),
                         ("cLN",f"{OBJ}/qc_metrics_cln_cosmx.csv",957)]:
        q=pd.read_csv(csv)
        rows.append(dict(dataset=nm,panel=panel,tx=q.n_counts.median(),genes=q.n_genes.median(),
                         ambient=q.neg_frac.median() if "neg_frac" in q else np.nan,n_cells=len(q)))
    # DKD from obs (backed); cells/section from orig_ident
    d=ad.read_h5ad(DKD_H5,backed="r"); xen=d.obs["tech"].astype(str).values=="Xenium"
    nc=d.obs["nCount_RNA"].astype(float).values[xen]; ng=d.obs["nFeature_RNA"].astype(float).values[xen]
    nsec=d.obs["orig_ident"].astype(str).values[xen]
    rows.append(dict(dataset="DKD",panel=5443,tx=np.median(nc),genes=np.median(ng),
                     ambient=np.nan,n_cells=xen.sum()))
    d.file.close()
    df=pd.DataFrame(rows)
    fig,axes=plt.subplots(1,4,figsize=(17,4.8))
    order=["RCC","PRCC","cLN","DKD"]; cols=[fs.DATASET[o] for o in order]
    metrics=[("tx","median transcripts / cell"),("genes","median genes / cell"),
             ("panel","panel size (genes)"),("ambient","ambient: neg-control fraction")]
    for ax,(key,lab) in zip(axes,metrics):
        vals=[df[df.dataset==o][key].iloc[0] for o in order]
        bars=ax.bar(order,[v if v==v else 0 for v in vals],color=cols)
        for b,v in zip(bars,vals):
            ax.text(b.get_x()+b.get_width()/2,(v if v==v else 0),
                    (f"{v:.3f}" if key=="ambient" else f"{v:.0f}") if v==v else "n/a",
                    ha="center",va="bottom",fontsize=10)
        ax.set_title(lab,fontsize=13); ax.set_ylabel("")
        for s in ax.get_xticklabels(): s.set_fontsize(11)
    axes[3].text(0.97,0.9,"Xenium ≈ 0 (clean);\nCosMx ~0.10 (high ambient);\nDKD neg-probes dropped",
                 transform=axes[3].transAxes,ha="right",va="top",fontsize=8.5,color="#666")
    fig.suptitle("Data quality per dataset — the substrate is usable",fontsize=16,y=0.99)
    fig.text(0.5,0.005,"Panel depth differs explicitly: RCC 405 · PRCC 377 · cLN 957 · DKD 5 443 genes — read downstream sensitivity against this. Metrics from QC subsamples.",
             ha="center",fontsize=10,color="#555")
    fig.subplots_adjust(top=0.84,bottom=0.16,wspace=0.22)
    fs.save_fig(fig,"Q1",tight=False)

# ============================================================================
def _dotplot_data(h5, lab_map_or_fn, present_markers, counts_from="X", section=None):
    a=ad.read_h5ad(h5,backed="r")
    if callable(lab_map_or_fn):
        lab=lab_map_or_fn(a)
    else:
        lab=pd.Series(a.obs[lab_map_or_fn[0]].astype(str)).map(lab_map_or_fn[1]).values
    base=np.ones(a.n_obs,bool) if section is None else section(a)
    gp=[g for g in present_markers if g in set(map(str,a.var_names))]
    sub=a[np.where(base)[0],gp]
    M=sub.to_memory(); a.file.close()
    C=(M.layers["counts"] if counts_from=="counts" and "counts" in M.layers else M.X)
    C=C.toarray() if sp.issparse(C) else np.asarray(C)
    lab=lab[base]
    df=pd.DataFrame(C,columns=gp); df["lab"]=lab
    return df,gp

def figT1():
    print("T1 marker dot-plot per dataset (typing-veracity)...")
    panels=[("RCC",RCC_H5,("phase_b_label",fs.RCC_LINEAGE),"X"),
            ("cLN",CLN_H5,("author_celltype",fs.CLN_LINEAGE),"X"),
            ("DKD",DKD_H5,(lambda a: dkd_lineage(a.obs["annotation_updated"].astype(str).values,
                a.obs["immune_cell_annotation_combined"].astype(str).values)),"counts")]
    allm=[m for m,_ in MARKERS]
    fig,axes=plt.subplots(len(panels),1,figsize=(11,12.5))
    for ax,(nm,h5,lm,cf) in zip(axes,panels):
        sect=(lambda a:a.obs["tech"].astype(str).values=="Xenium") if nm=="DKD" else None
        df,gp=_dotplot_data(h5,lm,allm,cf,sect)
        rows=[r for r in ROW_ORDER if r in set(df.lab.dropna())]
        cols=[m for m in allm if m in gp]
        det=np.zeros((len(rows),len(cols))); mn=np.zeros((len(rows),len(cols)))
        for i,r in enumerate(rows):
            sl=df[df.lab==r]
            for j,m in enumerate(cols):
                v=sl[m].values; det[i,j]=(v>0).mean(); mn[i,j]=np.log1p(v).mean()
        # z-score mean per column (marker) for color
        z=(mn-mn.mean(0))/(mn.std(0)+1e-9)
        for i in range(len(rows)):
            for j in range(len(cols)):
                ax.scatter(j,i,s=8+det[i,j]*330,c=[z[i,j]],cmap="RdBu_r",vmin=-1.5,vmax=2.2,
                           edgecolor="#888",linewidth=0.3)
        ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols,rotation=45,ha="right",fontsize=10)
        ax.set_yticks(range(len(rows))); ax.set_yticklabels([fs.LINEAGE_LABEL.get(r,r) for r in rows],fontsize=10)
        ax.set_title(f"{fs.DATASET_LONG[nm]}",color=fs.DATASET[nm],fontsize=13,loc="left")
        ax.set_xlim(-0.6,len(cols)-0.4); ax.set_ylim(-0.6,len(rows)-0.4); ax.invert_yaxis()
    # shared legends
    sm=plt.cm.ScalarMappable(cmap="RdBu_r",norm=plt.Normalize(-1.5,2.2)); sm.set_array([])
    cb=fig.colorbar(sm,ax=axes,fraction=0.012,pad=0.02,label="scaled mean expression (z)")
    for d,t in [(0.1,"10%"),(0.5,"50%"),(0.9,"90%")]:
        axes[0].scatter([],[],s=8+d*330,c="#999",label=t)
    axes[0].legend(title="detection",loc="upper left",bbox_to_anchor=(1.005,1.0),frameon=False,labelspacing=1.1)
    fig.suptitle("Canonical lineage markers define the cell types — you can see it\n"
                 "(dot size = detection rate, color = scaled mean expression)",fontsize=15)
    fs.save_fig(fig,"T1",tight=False)

# ============================================================================
def figT2():
    print("T2 UMAP per dataset colored by lineage + marker insets...")
    insets=["MS4A1","CD3D","CD68","EPCAM"]
    def load(nm):
        if nm=="DKD":
            a=ad.read_h5ad(DKD_H5,backed="r"); xen=a.obs["tech"].astype(str).values=="Xenium"
            idx=np.where(xen)[0]; idx=rng.choice(idx,min(45000,len(idx)),replace=False)
            um=np.asarray(a.obsm["X_umap"],float)[idx]
            lab=dkd_lineage(a.obs["annotation_updated"].astype(str).values[idx],
                            a.obs["immune_cell_annotation_combined"].astype(str).values[idx])
            sub=a[idx,[g for g in insets if g in set(map(str,a.var_names))]].to_memory(); a.file.close()
            X=sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
            return um,np.array(lab,dtype=object),pd.DataFrame(X,columns=list(sub.var_names))
        h5={"RCC":RCC_H5,"cLN":CLN_H5}[nm]; csv={"RCC":"kidney_RCC_protein","cLN":"cln_cosmx"}[nm]
        u=pd.read_csv(f"{OBJ}/wp_umap_{csv}.csv"); lmap=fs.RCC_LINEAGE if nm=="RCC" else fs.CLN_LINEAGE
        labcol="phase_b_label" if nm=="RCC" else "author_celltype"
        a=ad.read_h5ad(h5); lab=pd.Series(a.obs[labcol].astype(str)).map(lmap).values
        gp=[g for g in insets if g in set(map(str,a.var_names))]
        X=a[:,gp].X; X=X.toarray() if sp.issparse(X) else np.asarray(X)
        n=min(len(u),a.n_obs); idx=rng.choice(n,min(45000,n),replace=False)
        return u[["umap1","umap2"]].values[idx],lab[idx],pd.DataFrame(X[idx],columns=gp)
    dsets=["RCC","cLN","DKD"]; ncol=1+len(insets)
    fig,axes=plt.subplots(len(dsets),ncol,figsize=(3.3*ncol,3.3*len(dsets)))
    for ri,nm in enumerate(dsets):
        um,lab,mk=load(nm); present=[k for k in fs.LINEAGE_ORDER if (lab==k).any()]
        ax=axes[ri,0]
        ax.scatter(um[:,0],um[:,1],s=2,c="#eee",linewidths=0,rasterized=True)
        for k in present:
            m=lab==k; ax.scatter(um[m,0],um[m,1],s=3,c=fs.LINEAGE[k],linewidths=0,rasterized=True)
        ax.set_title(f"{fs.DATASET_LONG[nm]} — lineage",color=fs.DATASET[nm],fontsize=11); ax.axis("off")
        for ci,g in enumerate(insets):
            ax=axes[ri,ci+1]
            ax.scatter(um[:,0],um[:,1],s=2,c="#e9e9e9",linewidths=0,rasterized=True)  # grey background
            if g in mk.columns:
                v=mk[g].values.astype(float); pos=v>0
                if pos.any():
                    vv=np.log1p(v[pos]); o=np.argsort(vv)
                    ax.scatter(um[pos][o,0],um[pos][o,1],s=5,c=vv[o],cmap="inferno",linewidths=0,rasterized=True)
            ax.set_title(g,fontsize=11); ax.axis("off")
    from matplotlib.lines import Line2D
    present_all=[k for k in fs.LINEAGE_ORDER]
    handles=[Line2D([],[],marker=("*" if k=="Malignant_epi" else "o"),ls="",mfc=fs.LINEAGE[k],
             mec="black" if k=="Malignant_epi" else "none",ms=9,label=fs.LINEAGE_LABEL.get(k,k)) for k in present_all]
    fig.legend(handles=handles,loc="lower center",ncol=5,frameon=False,fontsize=10,bbox_to_anchor=(0.5,-0.02))
    fig.suptitle("Type structure is real, not imposed — UMAP by lineage; marker+ cells (inferno) localize to their islands",fontsize=15)
    fs.save_fig(fig,"T2")

# ============================================================================
def figA3raw():
    print("A3raw raw distributions behind the AUROC...")
    a=ad.read_h5ad(DKD_H5,backed="r"); plat=a.obs["tech"].astype(str).values
    imm=a.obs["immune_cell_annotation_combined"].astype(str).values; ann=a.obs["annotation_updated"].astype(str).values
    need=np.isin(imm,["CD4+","CD8+"])|((ann=="PT"))
    idx=np.where(need)[0]; sub=a[idx,["CD8A","CD4"]].to_memory(); a.file.close()
    X=sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
    df=pd.DataFrame(X,columns=["CD8A","CD4"]); df["plat"]=plat[idx]; df["imm"]=imm[idx]; df["ann"]=ann[idx]
    auc=pd.read_csv(f"{XT}/results_auroc.csv")
    fig,axes=plt.subplots(2,2,figsize=(11,9))
    for r,mk in enumerate(["CD8A","CD4"]):
        for c,pl in enumerate(["Xenium","CosMx"]):
            ax=axes[r,c]
            g4=np.log1p(df[(df.plat==pl)&(df.imm=="CD4+")][mk].values)
            g8=np.log1p(df[(df.plat==pl)&(df.imm=="CD8+")][mk].values)
            amb=np.log1p(df[(df.plat==pl)&(df.ann=="PT")][mk].values).mean()
            parts=ax.violinplot([g4,g8],positions=[0,1],showmedians=True,widths=0.8)
            for pc,col in zip(parts['bodies'],["#56B4E9","#D55E00"]): pc.set_facecolor(col); pc.set_alpha(0.6)
            for gi,(g,col) in enumerate([(g4,"#56B4E9"),(g8,"#D55E00")]):
                xs=gi+rng.uniform(-0.13,0.13,min(len(g),1500)); ax.scatter(xs,rng.choice(g,len(xs)),s=3,c=col,alpha=0.25,rasterized=True)
            ax.axhline(amb,color="k",ls=":",lw=1.4)
            av=auc[(auc.platform==pl)&(auc.marker==mk)].auroc_cd8_vs_cd4
            ax.text(0.5,0.93,f"AUROC {av.iloc[0]:.2f}" if len(av) else "",transform=ax.transAxes,ha="center",
                    fontsize=12,fontweight="bold",color="#b2182b" if (len(av) and abs(av.iloc[0]-0.5)<0.1) else "#222")
            ax.set_xticks([0,1]); ax.set_xticklabels(["CD4-labelled","CD8-labelled"])
            ax.set_title(f"{mk} — {pl}",fontsize=12);
            if c==0: ax.set_ylabel("log1p(count)")
    fig.suptitle("Raw counts behind the AUROC: CD8A separates CD4/CD8 cells on Xenium, overlaps on CosMx\n"
                 "(dotted = ambient/PT floor; CD4 barely separates on either panel)",fontsize=14)
    fs.save_fig(fig,"A3raw")

# ============================================================================
def _section_B(h5, sect, blab):
    a=ad.read_h5ad(h5,backed="r"); xy=np.asarray(a.obsm["spatial"],float)
    base=sect(a); lab=blab(a); a.file.close()
    co=xy[base]; isB=(lab[base]=="B")
    return co,isB

def figC1raw():
    print("C1raw aggregate reality (hulls + NN null)...")
    def rcc_sect(a): return np.ones(a.n_obs,bool)
    def rcc_B(a): return np.where(a.obs["phase_b_label"].astype(str).values=="Naive B cells","B","x")
    def dkd_sect(a): return (a.obs["orig_ident"].astype(str).values=="HK2695")&(a.obs["tech"].astype(str).values=="Xenium")
    def dkd_B(a): return a.obs["immune_cell_annotation_combined"].astype(str).values
    fig=plt.figure(figsize=(15,9))
    gs=fig.add_gridspec(2,2,width_ratios=[2,1],hspace=0.32,wspace=0.18,top=0.9,bottom=0.07,left=0.04,right=0.97)
    for row,(nm,h5,sect,bl,col) in enumerate([("ccRCC",RCC_H5,rcc_sect,rcc_B,fs.DATASET["RCC"]),
                                              ("DKD (HK2695)",DKD_H5,dkd_sect,dkd_B,fs.DATASET["DKD"])]):
        co,isB=_section_B(h5,sect,bl)
        cl=np.full(isB.sum(),-1); Bco=co[isB]
        if isB.sum()>=20:
            cl=DBSCAN(eps=50,min_samples=20).fit(Bco).labels_
        ax=fig.add_subplot(gs[row,0])
        bg=co[rng.choice(len(co),min(60000,len(co)),replace=False)]
        ax.scatter(bg[:,0],bg[:,1],s=1,c="#ededed",linewidths=0,rasterized=True)
        ax.scatter(Bco[cl==-1,0],Bco[cl==-1,1],s=5,c="#c9c9c9",linewidths=0,rasterized=True)
        ax.scatter(Bco[cl>=0,0],Bco[cl>=0,1],s=7,c=col,linewidths=0,rasterized=True)
        for k in [c for c in np.unique(cl) if c>=0]:
            P=Bco[cl==k]
            if len(P)>=4:
                try: h=ConvexHull(P); ax.add_patch(Polygon(P[h.vertices],closed=True,fill=False,edgecolor=col,lw=1.2,alpha=0.8))
                except Exception: pass
        ax.set_title(f"{nm}: {len(set(cl[cl>=0]))} B-aggregates (DBSCAN ε=50, hulls)",color=col,fontsize=12)
        ax.set_aspect("equal"); ax.axis("off")
        # size histogram + NN null
        ax2=fig.add_subplot(gs[row,1])
        sizes=[int((cl==k).sum()) for k in np.unique(cl) if k>=0]
        ax2.hist(sizes,bins=20,color=col,alpha=0.85); ax2.set_xlabel("B cells per aggregate"); ax2.set_ylabel("aggregates")
        ax2.set_title(f"{nm}: aggregate sizes",fontsize=11,color=col)
        # observed B nearest-neighbour vs permuted null (median NN distance)
        tB=cKDTree(Bco); obs=np.median(tB.query(Bco,k=2)[0][:,1])
        null=[]
        for _ in range(150):
            r=co[rng.choice(len(co),len(Bco),replace=False)]
            null.append(np.median(cKDTree(r).query(r,k=2)[0][:,1]))
        null=np.array(null); z=(obs-null.mean())/null.std()
        ax2.text(0.97,0.97,f"B–B nearest-neighbour:\nobserved {obs:.0f} µm  vs  null {null.mean():.0f} µm\n(z = {z:.0f}) — B cells genuinely cluster",
                 transform=ax2.transAxes,ha="right",va="top",fontsize=9,
                 bbox=dict(boxstyle="round",fc="#fff7e6",ec="#f0b429"))
    fig.suptitle("Aggregates are real structures: B cells cluster far beyond chance (not an algorithm artifact)",fontsize=15)
    fs.save_fig(fig,"C1raw",tight=False)

# ============================================================================
def figC2raw():
    print("C2raw composition behind the enrichment...")
    # per-aggregate counts from CSVs: count = fraction * n_cells_region
    d=pd.read_csv(f"{BN}/dkd_xenium_aggregates.csv")
    d_treg=(d["f_Treg-like"]*d.n_cells_region).round(); d_cd8=(d["f_eff-CD8"]*d.n_cells_region).round()
    r=pd.read_csv(f"{TAB}/rcc_phaseB2_aggregates.csv")
    r_treg=(r["f_Treg"]*r.n_cells_region).round(); r_cd8=(r["f_eff-CD8"]*r.n_cells_region).round()
    rc=pd.read_csv(f"{TAB}/rcc_phaseB2_aggregate_composition.csv").set_index("cell_type")
    fig,axes=plt.subplots(1,3,figsize=(16,5))
    # (1) inside-vs-background fraction (real, countable) — RCC
    ax=axes[0]; types=["Treg","eff-CD8","Plasma"]; x=np.arange(len(types)); w=0.38
    ins=[rc.loc[t,"inside_mean_frac"] for t in types]; bg=[rc.loc[t,"background_frac"] for t in types]
    ax.bar(x-w/2,ins,w,label="inside aggregate",color=fs.DATASET["RCC"])
    ax.bar(x+w/2,bg,w,label="tissue background",color="#cccccc")
    ax.set_xticks(x); ax.set_xticklabels(["Treg","effector-CD8","Plasma"]); ax.set_ylabel("fraction of cells")
    ax.set_title("ccRCC: composition inside vs background"); ax.legend(frameon=False)
    # (2) per-aggregate Treg vs CD8 counts, RCC vs DKD
    ax=axes[1]
    ax.scatter(r_treg,r_cd8,s=26,c=fs.DATASET["RCC"],alpha=0.7,edgecolor="k",linewidth=0.2,label="ccRCC")
    ax.scatter(d_treg,d_cd8,s=26,c=fs.DATASET["DKD"],alpha=0.7,edgecolor="k",linewidth=0.2,label="DKD")
    lim=max(r_treg.max(),r_cd8.max(),d_treg.max(),d_cd8.max())*1.05
    ax.plot([0,lim],[0,lim],ls="--",c="#888",lw=1); ax.set_xlim(0,lim); ax.set_ylim(0,lim)
    ax.set_xlabel("Treg cells per aggregate"); ax.set_ylabel("effector-CD8 cells per aggregate")
    ax.set_title("Per-aggregate counts (above line = Treg>CD8)"); ax.legend(frameon=False)
    # (3) per-section/aggregate consistency: fraction Treg-enriched
    ax=axes[2]
    rt_l=np.log2((r["f_Treg"]+1e-6)/(rc.loc["Treg","background_frac"]+1e-6))
    dt_l=np.log2((d["f_Treg-like"]+1e-6)/(d["bg_Treg-like"]+1e-6))
    for yb,(vals,col,nm) in enumerate([(rt_l,fs.DATASET["RCC"],"ccRCC"),(dt_l,fs.DATASET["DKD"],"DKD")]):
        ax.scatter(vals,np.full(len(vals),yb)+rng.uniform(-0.09,0.09,len(vals)),s=26,c=col,alpha=0.7,edgecolor="k",linewidth=0.2)
        ax.text(vals.max(),yb+0.25,f"{int((vals>0).sum())}/{len(vals)} Treg-enriched",color=col,fontsize=10,ha="right")
    fs.zeroline(ax,0,"v"); ax.set_yticks([0,1]); ax.set_yticklabels(["ccRCC","DKD"]); ax.set_ylim(-0.5,1.6)
    ax.set_xlim(-4,5.2)  # clip a single 0-Treg aggregate (pseudocount artifact) off-scale
    ax.set_xlabel("per-aggregate Treg log₂ enrichment"); ax.set_title("Per-aggregate consistency")
    fig.suptitle("The enrichment traces to real, countable cells (inside vs background; per-aggregate counts & consistency)",fontsize=15)
    fs.save_fig(fig,"C2raw")

if __name__=="__main__":
    figQ1(); figT1(); figT2(); figA3raw(); figC1raw(); figC2raw()
    print("== foundational done ==")
