#!/usr/bin/env python
"""build_vignettes.py — per-study vignette figures (native vocabulary, all data) +
normalization/mapping section. Reuses figstyle. Read-only; backed DKD; per-section; downsample
for display only (saved files full-res 300 DPI)."""
import os, warnings; warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from sklearn.cluster import DBSCAN
from scipy import ndimage
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figstyle as fs

def region_cores(xy, nbins=55, min_cells=400, margin=0.06, max_cores=2):
    """Detect spatially separated tissue cores on a (multi-core CosMx) slide via dense-bin
    connected components; return tight bounding boxes (drop empty space). Unit-agnostic."""
    if len(xy)<min_cells: return [(xy[:,0].min(),xy[:,0].max(),xy[:,1].min(),xy[:,1].max(),len(xy))]
    H,xe,ye=np.histogram2d(xy[:,0],xy[:,1],bins=nbins)
    thr=max(1,np.percentile(H[H>0],45))
    lab,n=ndimage.label(H>=thr); cores=[]
    for k in range(1,n+1):
        b=np.argwhere(lab==k); x0,x1=xe[b[:,0].min()],xe[b[:,0].max()+1]; y0,y1=ye[b[:,1].min()],ye[b[:,1].max()+1]
        mx=(x1-x0)*margin; my=(y1-y0)*margin
        sel=(xy[:,0]>=x0-mx)&(xy[:,0]<=x1+mx)&(xy[:,1]>=y0-my)&(xy[:,1]<=y1+my)
        if sel.sum()>=min_cells: cores.append((x0-mx,x1+mx,y0-my,y1+my,int(sel.sum())))
    cores.sort(key=lambda c:-c[4])
    return cores[:max_cores] if cores else [(xy[:,0].min(),xy[:,0].max(),xy[:,1].min(),xy[:,1].max(),len(xy))]
def crop_mask(xy,bb): return (xy[:,0]>=bb[0])&(xy[:,0]<=bb[1])&(xy[:,1]>=bb[2])&(xy[:,1]<=bb[3])
rng=np.random.default_rng(0); REPO=fs.REPO; OBJ=f"{REPO}/outputs/objects"
RCC_H5=f"{OBJ}/kidney_RCC_protein.h5ad"; PRCC_H5=f"{OBJ}/kidney_preview_PRCC.h5ad"
CLN_H5=f"{OBJ}/cln_cosmx.h5ad"; DKD_H5=f"{REPO}/Demoulin26/data/spatial_adata_xenium_cosmx_zenodo.h5ad"

MARKERS=["MS4A1","CD79A","MZB1","CD3D","CD3E","FOXP3","CD8A","GZMB","CD68","LYZ","PECAM1","VWF","EPCAM","KRT8"]
GAL_MARK=["MS4A1","MZB1","FOXP3","CD8A","CD68"]   # markers for the spatial gallery
def hdr(s): print("\n=== "+s+" ===")

# ---------- generic spatial panel ----------
def sp_marker(ax,xy,v,title,cmap="inferno",bgs=2):
    ax.scatter(xy[:,0],xy[:,1],s=bgs,c="#ececec",linewidths=0,rasterized=True)
    pos=v>0
    if pos.any():
        vv=np.log1p(v[pos]); o=np.argsort(vv)
        ax.scatter(xy[pos][o,0],xy[pos][o,1],s=7,c=vv[o],cmap=cmap,linewidths=0,rasterized=True)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title,fontsize=11)
def sp_pop(ax,xy,mask,title,col,bgs=2):
    ax.scatter(xy[:,0],xy[:,1],s=bgs,c="#ececec",linewidths=0,rasterized=True)
    ax.scatter(xy[mask,0],xy[mask,1],s=9,c=col,linewidths=0,rasterized=True)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title,fontsize=11)

# ---------- generic native dot-plot ----------
def native_dotplot(ax,df,labels,markers,row_order,title,color):
    rows=[r for r in row_order if (labels==r).sum()>=30]
    cols=[m for m in markers if m in df.columns]
    mn=np.zeros((len(rows),len(cols))); det=np.zeros_like(mn)
    for i,r in enumerate(rows):
        sl=df[labels==r]
        for j,m in enumerate(cols):
            v=sl[m].values; det[i,j]=(v>0).mean(); mn[i,j]=np.log1p(v).mean()
    z=(mn-mn.mean(0))/(mn.std(0)+1e-9)
    for i in range(len(rows)):
        for j in range(len(cols)):
            ax.scatter(j,i,s=6+det[i,j]*230,c=[z[i,j]],cmap="RdBu_r",vmin=-1.5,vmax=2.2,edgecolor="#999",linewidth=0.25)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols,rotation=45,ha="right",fontsize=9)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows,fontsize=8)
    ax.set_xlim(-0.6,len(cols)-0.4); ax.set_ylim(-0.6,len(rows)-0.4); ax.invert_yaxis()
    ax.set_title(title,color=color,fontsize=12,loc="left")
    return rows,cols

def _load(h5,markers,label_col,section=None,sub=None):
    a=ad.read_h5ad(h5,backed="r"); xy=np.asarray(a.obsm["spatial"],float)
    base=np.ones(a.n_obs,bool) if section is None else section(a)
    lab=a.obs[label_col].astype(str).values
    idx=np.where(base)[0]
    if sub and len(idx)>sub: idx=rng.choice(idx,sub,replace=False); idx.sort()
    gp=[g for g in markers if g in set(map(str,a.var_names))]
    M=a[idx,gp].to_memory(); a.file.close()
    X=M.X.toarray() if sp.issparse(M.X) else np.asarray(M.X)
    return xy[idx],pd.DataFrame(X,columns=gp),lab[idx]

# qc panel reused
def qc_panel(ax,order,vals,colors,lab,fmt="{:.0f}",ref=None):
    bars=ax.bar(order,[v if v==v else 0 for v in vals],color=colors)
    for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,(v if v==v else 0),
        fmt.format(v) if v==v else "n/a",ha="center",va="bottom",fontsize=9)
    ax.set_title(lab,fontsize=12)
    if ref is not None: ax.axhline(ref,color="#888",ls="--",lw=1)

# ============================================================================
# V1 — RCC (Xenium x2: ccRCC + pRCC)
# ============================================================================
RCC_POPS={"B cell":(["Naive B cells","Switched memory B cells"],"#1f77b4"),
 "Plasma":(["Plasmablasts"],"#ff7f0e"),"Treg":(["T regulatory cells"],"#d62728"),
 "Effector CD8":(["Effector memory CD8 T cells","CD8_T"],"#2ca02c"),
 "Myeloid":(["Intermediate monocytes","Classical monocytes","myeloid/DC"],"#9467bd")}
RCC_ROWS=["Naive B cells","Plasmablasts","T regulatory cells","Effector memory CD8 T cells","CD8_T",
 "CCR7+ T (naive/CM)","Intermediate monocytes","Classical monocytes","myeloid/DC","mregDC",
 "Endothelial","Stroma_mural","Fibroblast","Proximal_tubule","Tumor_RCC","Mast"]

def V1_qc():
    print("V1_qc..."); fig,axes=plt.subplots(1,4,figsize=(15,4.4))
    order=["ccRCC","pRCC"]; cols=[fs.DATASET["RCC"],fs.DATASET["PRCC"]]
    q1=pd.read_csv(f"{OBJ}/qc_metrics_kidney_RCC_protein.csv"); q2=pd.read_csv(f"{OBJ}/qc_metrics_kidney_preview_PRCC.csv")
    qc_panel(axes[0],order,[q1.n_counts.median(),q2.n_counts.median()],cols,"median transcripts / cell")
    qc_panel(axes[1],order,[q1.n_genes.median(),q2.n_genes.median()],cols,"median genes / cell")
    qc_panel(axes[2],order,[405,377],cols,"panel size (genes)")
    qc_panel(axes[3],order,[q1.neg_frac.median(),q2.neg_frac.median()],cols,"neg-control fraction","{:.3f}")
    fig.suptitle("V1 · ccRCC — data quality, two Xenium sections (discovery + replication)",fontsize=15)
    fs.save_fig(fig,"V1_qc")

def V1_typing():
    print("V1_typing...")
    xr,dr,lr=_load(RCC_H5,MARKERS,"phase_b_label",sub=120000)
    fig=plt.figure(figsize=(16,6)); gs=fig.add_gridspec(1,3,width_ratios=[1.4,1,1])
    ax=fig.add_subplot(gs[0]); native_dotplot(ax,dr,lr,MARKERS,RCC_ROWS,"ccRCC native typing (markers → types)",fs.DATASET["RCC"])
    # umap RCC + PRCC colored by native lineage (tab20)
    for k,(nm,csv,col) in enumerate([("ccRCC","kidney_RCC_protein",fs.DATASET["RCC"]),("pRCC","kidney_preview_PRCC",fs.DATASET["PRCC"])]):
        ax=fig.add_subplot(gs[k+1]); u=pd.read_csv(f"{OBJ}/wp_umap_{csv}.csv")
        labs=u["label"].astype(str).values; cats=pd.Series(labs).value_counts().index[:18]
        cmap=plt.get_cmap("tab20")
        ax.scatter(u.umap1,u.umap2,s=2,c="#eee",linewidths=0,rasterized=True)
        for ci,c in enumerate(cats):
            m=labs==c; ax.scatter(u.umap1[m],u.umap2[m],s=3,c=[cmap(ci%20)],linewidths=0,rasterized=True)
        ax.set_title(f"{nm} — UMAP (native clusters)",color=col,fontsize=12); ax.axis("off")
    fig.suptitle("V1 · ccRCC native cell typing — defined by canonical markers, structure visible in UMAP",fontsize=15)
    fs.save_fig(fig,"V1_typing",tight=False)

def V1_gallery():
    print("V1_gallery (markers + populations, both sections)...")
    secs=[("ccRCC",RCC_H5,None,fs.DATASET["RCC"]),("pRCC",PRCC_H5,None,fs.DATASET["PRCC"])]
    ncol=len(GAL_MARK)+len(RCC_POPS)
    fig,axes=plt.subplots(len(secs),ncol,figsize=(2.5*ncol,5.4))
    for ri,(nm,h5,sect,col) in enumerate(secs):
        xy,d,lab=_load(h5,GAL_MARK,"phase_b_label",sect,sub=120000)
        for ci,g in enumerate(GAL_MARK):
            sp_marker(axes[ri,ci],xy,d[g].values.astype(float) if g in d else np.zeros(len(xy)),g if ri==0 else "")
        for ci,(pn,(labs,pcol)) in enumerate(RCC_POPS.items()):
            sp_pop(axes[ri,len(GAL_MARK)+ci],xy,np.isin(lab,labs),pn if ri==0 else "",pcol)
        axes[ri,0].text(-0.12,0.5,nm,transform=axes[ri,0].transAxes,rotation=90,va="center",ha="center",
                        fontsize=13,fontweight="bold",color=col)
    fig.suptitle("V1 · ccRCC marker & population gallery — biology reproduces across BOTH Xenium sections\n"
                 "(left: marker intensity; right: native populations in situ)",fontsize=14)
    fs.save_fig(fig,"V1_gallery")

def V1_pattern():
    print("V1_pattern (Treg-around / CD8-excluded, both sections)...")
    rcc=pd.read_csv(f"{REPO}/outputs/tables/rcc_phaseB2_aggregates.csv")
    comp=pd.read_csv(f"{REPO}/outputs/tables/rcc_phaseB2_aggregate_composition.csv").set_index("cell_type")
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    # per-aggregate Treg vs CD8 counts (ccRCC) + composition inside vs bg
    rt=(rcc["f_Treg"]*rcc.n_cells_region).round(); rc8=(rcc["f_eff-CD8"]*rcc.n_cells_region).round()
    ax=axes[0]; ax.scatter(rt,rc8,s=30,c=fs.DATASET["RCC"],alpha=0.75,edgecolor="k",linewidth=0.2)
    lim=max(rt.max(),rc8.max())*1.05; ax.plot([0,lim],[0,lim],"--",c="#888",lw=1); ax.set_xlim(0,lim); ax.set_ylim(0,lim)
    ax.set_xlabel("Treg cells / aggregate"); ax.set_ylabel("effector-CD8 cells / aggregate")
    ax.set_title(f"ccRCC: Treg>CD8 in {int((rt>rc8).sum())}/{len(rt)} aggregates")
    ax=axes[1]; types=["Treg","eff-CD8","Plasma"]; x=np.arange(len(types)); w=0.38
    ax.bar(x-w/2,[comp.loc[t,"inside_mean_frac"] for t in types],w,label="inside aggregate",color=fs.DATASET["RCC"])
    ax.bar(x+w/2,[comp.loc[t,"background_frac"] for t in types],w,label="background",color="#ccc")
    ax.set_xticks(x); ax.set_xticklabels(["Treg","effector-CD8","Plasma"]); ax.set_ylabel("fraction of cells")
    ax.set_title("Treg enriched, cytotoxic-CD8 excluded"); ax.legend(frameon=False)
    fig.suptitle("V1 · ccRCC observed pattern: immunoregulatory aggregate (Treg-around / cytotoxic-excluded)",fontsize=14)
    fs.save_fig(fig,"V1_pattern")

def _crop_markers(h5,markers,cx,cy,W,label_col,bval,sect=None):
    a=ad.read_h5ad(h5,backed="r"); xy=np.asarray(a.obsm["spatial"],float)
    base=np.ones(a.n_obs,bool) if sect is None else sect(a); lab=a.obs[label_col].astype(str).values
    m=base&(np.abs(xy[:,0]-cx)<W)&(np.abs(xy[:,1]-cy)<W); idx=np.where(m)[0]
    gp=[g for g in markers if g in set(map(str,a.var_names))]; M=a[idx,gp].to_memory(); a.file.close()
    X=M.X.toarray() if sp.issparse(M.X) else np.asarray(M.X)
    return xy[idx],pd.DataFrame(X,columns=gp),(lab[idx]==bval)

def V1_closeup():
    print("V1_closeup (ccRCC aggregate marker overlay, B1-style)...")
    markers=["MS4A1","MZB1","FOXP3","CD8A","GZMB"]; W=95.0
    ragg=pd.read_csv(f"{REPO}/outputs/tables/rcc_phaseB2_aggregates.csv").sort_values("n_B",ascending=False).iloc[0]
    xy,mk,B=_crop_markers(RCC_H5,markers,ragg.x,ragg.y,W,"phase_b_label","Naive B cells")
    ncol=len(markers)+1; fig,axes=plt.subplots(1,ncol,figsize=(3.1*ncol,3.7))
    for ci,g in enumerate(markers):
        ax=axes[ci]; ax.scatter(xy[:,0],xy[:,1],s=6,c="#ededed",linewidths=0,rasterized=True)
        ax.scatter(xy[B,0],xy[B,1],s=11,c="#bcd6ee",linewidths=0,rasterized=True)
        if g in mk.columns:
            v=mk[g].values.astype(float); pos=v>0; vmax=max(np.percentile(v[pos],98),1) if pos.any() else 1
            ax.scatter(xy[pos,0],xy[pos,1],s=36,c=np.clip(v[pos],0,vmax),cmap=fs.MARKER_CMAP.get(g,"magma"),
                       vmin=0,vmax=vmax,edgecolor="white",linewidth=0.2,rasterized=True)
        ax.set_aspect("equal"); ax.axis("off"); ax.set_title(f"{g}\n{fs.MARKER_LINEAGE.get(g,'')}",fontsize=12)
    ax=axes[ncol-1]; ax.scatter(xy[:,0],xy[:,1],s=6,c="#ededed",linewidths=0,rasterized=True)
    ax.scatter(xy[B,0],xy[B,1],s=14,c="#1f77b4",alpha=0.8,linewidths=0,rasterized=True)
    if "FOXP3" in mk: fp=mk["FOXP3"].values>0; ax.scatter(xy[fp,0],xy[fp,1],s=42,c="#d62728",edgecolor="white",linewidth=0.3,rasterized=True)
    if "CD8A" in mk: c8=mk["CD8A"].values>0; ax.scatter(xy[c8,0],xy[c8,1],s=42,c="#2ca02c",marker="^",edgecolor="white",linewidth=0.3,rasterized=True)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title("Composite\nB / Treg / CD8",fontsize=12)
    ax.legend(handles=[Line2D([],[],marker="o",ls="",mfc="#1f77b4",mec="none",ms=8,label="B cell"),
        Line2D([],[],marker="o",ls="",mfc="#d62728",mec="white",ms=8,label="FOXP3+ (Treg)"),
        Line2D([],[],marker="^",ls="",mfc="#2ca02c",mec="white",ms=8,label="CD8A+ (cytotoxic)")],
        loc="upper center",bbox_to_anchor=(0.5,-0.02),frameon=False,fontsize=9)
    fig.suptitle("V1 · ccRCC aggregate close-up — Treg-around & cytotoxic-excluded, SEEN in one representative aggregate\n"
                 "(blue underlay = B-core; ~190 µm field)",fontsize=14)
    fs.save_fig(fig,"V1_closeup")

# ============================================================================
# V2 — cLN (CosMx, 14 slides)
# ============================================================================
CLN_POPS={"B cell":(["B-cell"],"#1f77b4"),"Plasma":(["plasmablast"],"#ff7f0e"),
 "Myeloid":(["macrophage","monocyte","mDC","pDC"],"#9467bd"),"T (low-conf)":(["T cell (low-conf)","T CD8 memory"],"#2ca02c")}
CLN_ROWS=["B-cell","plasmablast","macrophage","monocyte","mDC","NK","Treg","T CD8 memory",
 "Fibroblast","Myofibroblast","Glomerular.endothelium","Peritubular.capillary.endothelium.1",
 "Podocyte","PCT","Connecting.tubule","Thick.ascending.limb.of.Loop.of.Henle","Principal.cell"]

def V2_qc():
    print("V2_qc..."); q=pd.read_csv(f"{OBJ}/qc_metrics_cln_cosmx.csv")
    a=ad.read_h5ad(CLN_H5,backed="r"); cond=a.obs["condition"].astype(str).values; samp=a.obs["sample"].astype(str).values; a.file.close()
    cps=pd.Series(samp).value_counts()
    fig,axes=plt.subplots(1,3,figsize=(14,4.4))
    axes[0].bar(["transcripts","genes"],[q.n_counts.median(),q.n_genes.median()],color=[fs.DATASET["cLN"],"#f08"])
    axes[0].set_title("median per cell (CosMx 957-plex)")
    axes[1].bar(cps.index,cps.values,color=fs.DATASET["cLN"]); axes[1].set_xticklabels(cps.index,rotation=90,fontsize=7)
    axes[1].set_title(f"cells per slide (14 slides, {len(samp):,} cells)"); axes[1].set_ylabel("cells")
    axes[2].bar(["neg-control\nfraction"],[q.neg_frac.median()],color=fs.DATASET["cLN"]); axes[2].set_ylim(0,0.15)
    axes[2].text(0,q.neg_frac.median(),f"{q.neg_frac.median():.3f}",ha="center",va="bottom")
    axes[2].set_title("ambient (CosMx ~0.10)")
    fig.suptitle("V2 · cLN (childhood lupus nephritis) — CosMx 957-plex, 14 slides (4 control + 10 SLE)\n"
                 "CAVEAT up front: ~35% epithelial CD3 contamination → T-lineage UNRELIABLE here",fontsize=13)
    fs.save_fig(fig,"V2_qc")

def V2_typing():
    print("V2_typing...")
    xc,dc,lc=_load(CLN_H5,MARKERS,"author_celltype",sub=120000)
    fig=plt.figure(figsize=(15,6)); gs=fig.add_gridspec(1,2,width_ratios=[1.5,1])
    ax=fig.add_subplot(gs[0]); native_dotplot(ax,dc,lc,MARKERS,CLN_ROWS,"cLN native typing (author labels)",fs.DATASET["cLN"])
    ax=fig.add_subplot(gs[1]); u=pd.read_csv(f"{OBJ}/wp_umap_cln_cosmx.csv")
    labs=u["label"].astype(str).values; cats=pd.Series(labs).value_counts().index[:18]; cmap=plt.get_cmap("tab20")
    ax.scatter(u.umap1,u.umap2,s=2,c="#eee",linewidths=0,rasterized=True)
    for ci,c in enumerate(cats):
        m=labs==c; ax.scatter(u.umap1[m],u.umap2[m],s=3,c=[cmap(ci%20)],linewidths=0,rasterized=True)
    ax.set_title("cLN — UMAP (native clusters)",color=fs.DATASET["cLN"],fontsize=12); ax.axis("off")
    fig.suptitle("V2 · cLN native cell typing (see A2 for InSituType benchmark vs author labels)",fontsize=15)
    fs.save_fig(fig,"V2_typing",tight=False)

def _cln_pick():
    a=ad.read_h5ad(CLN_H5,backed="r"); samp=a.obs["sample"].astype(str).values
    cond=a.obs["condition"].astype(str).values; lab=a.obs["author_celltype"].astype(str).values; a.file.close()
    bcount={s:int(((np.isin(lab,["B-cell","plasmablast"]))&(samp==s)).sum()) for s in pd.unique(samp)}
    ctrl=[s for s in pd.unique(samp) if cond[samp==s][0]=="control"]
    sle=[s for s in pd.unique(samp) if cond[samp==s][0]!="control"]
    pick=sorted(ctrl,key=lambda s:-bcount[s])[:3]+sorted(sle,key=lambda s:-bcount[s])[:3]
    return pick,{s:cond[samp==s][0] for s in pick}

def V2_gallery():
    print("V2_gallery (region-cropped cores, condition-spanning slides)...")
    pick,cmap=_cln_pick(); gmark=["MS4A1","MZB1","CD68"]; ncol=len(gmark)+len(CLN_POPS)
    fig,axes=plt.subplots(len(pick),ncol,figsize=(2.4*ncol,2.4*len(pick)))
    for ri,s in enumerate(pick):
        def sect(a,s=s): return a.obs["sample"].astype(str).values==s
        xy,d,lb=_load(CLN_H5,gmark,"author_celltype",sect)          # full slide
        bb=region_cores(xy)[0]; cm=crop_mask(xy,bb)                  # crop to densest core
        xy,d,lb=xy[cm],d[cm].reset_index(drop=True),lb[cm]
        cnd=cmap[s]
        for ci,g in enumerate(gmark): sp_marker(axes[ri,ci],xy,d[g].values.astype(float) if g in d else np.zeros(len(xy)),g if ri==0 else "")
        for ci,(pn,(labs,pcol)) in enumerate(CLN_POPS.items()): sp_pop(axes[ri,len(gmark)+ci],xy,np.isin(lb,labs),pn if ri==0 else "",pcol)
        axes[ri,0].text(-0.18,0.5,f"{s}\n({cnd})",transform=axes[ri,0].transAxes,rotation=90,va="center",ha="center",
                        fontsize=9,fontweight="bold",color="#444" if cnd=="control" else fs.DATASET["cLN"])
    fig.suptitle("V2 · cLN marker & population gallery — cropped to the dense tissue CORE per slide (3 control + 3 SLE)\n"
                 "cells fill each frame; reproducibility across the cohort",fontsize=14)
    fs.save_fig(fig,"V2_gallery")

def V2_niche():
    print("V2_niche (region-cropped spatial plasma–myeloid, control vs SLE)...")
    pick,cmap=_cln_pick()
    pops={"plasmablast":"#ff7f0e","macrophage":"#9467bd","monocyte":"#9467bd","mDC":"#9467bd"}
    fig,axes=plt.subplots(len(pick),3,figsize=(11,2.4*len(pick)))
    for ri,s in enumerate(pick):
        def sect(a,s=s): return a.obs["sample"].astype(str).values==s
        xy,_,lb=_load(CLN_H5,["MS4A1"],"author_celltype",sect)
        bb=region_cores(xy)[0]; cm=crop_mask(xy,bb); xy,lb=xy[cm],lb[cm]
        plas=lb=="plasmablast"; mye=np.isin(lb,["macrophage","monocyte","mDC","pDC"])
        sp_pop(axes[ri,0],xy,plas,"Plasma" if ri==0 else "","#ff7f0e")
        sp_pop(axes[ri,1],xy,mye,"Myeloid" if ri==0 else "","#9467bd")
        ax=axes[ri,2]; ax.scatter(xy[:,0],xy[:,1],s=2,c="#ececec",linewidths=0,rasterized=True)
        ax.scatter(xy[mye,0],xy[mye,1],s=9,c="#9467bd",linewidths=0,rasterized=True)
        ax.scatter(xy[plas,0],xy[plas,1],s=11,c="#ff7f0e",linewidths=0,rasterized=True)
        ax.set_aspect("equal"); ax.axis("off")
        if ri==0: ax.set_title("Plasma + Myeloid (overlay)",fontsize=11)
        cnd=cmap[s]
        axes[ri,0].text(-0.2,0.5,f"{s}\n({cnd})",transform=axes[ri,0].transAxes,rotation=90,va="center",ha="center",
                        fontsize=9,fontweight="bold",color="#444" if cnd=="control" else fs.DATASET["cLN"])
    fig.suptitle("V2 · cLN plasma–myeloid niche in situ, per tissue core (control vs SLE)\n"
                 "plasma aggregates recruit myeloid cells — the reproducible cLN finding (region-cropped)",fontsize=13)
    fs.save_fig(fig,"V2_niche")

# ============================================================================
# V3 — DKD (CosMx + Xenium)
# ============================================================================
def dkd_native(a,plat):
    ann=a.obs["annotation_updated"].astype(str).values; imm=a.obs["immune_cell_annotation_combined"].astype(str).values
    lab=np.where(imm!="Unknown",imm,ann); return lab
DKD_ROWS=["B","Plasma","Macro","cDC","pDC","NK","CD4+","CD8+","Fibroblast","EC_Peritub","EC_glom",
 "Podo","PT","iPT","TAL","DCT","CNT","PC","IC A"]
DKD_POPS={"B cell":(["B"],"#1f77b4"),"Plasma":(["Plasma"],"#ff7f0e"),"Myeloid":(["Macro","cDC","pDC"],"#9467bd"),
 "T lineage":(["CD4+","CD8+"],"#2ca02c")}

def V3_qc():
    print("V3_qc...")
    a=ad.read_h5ad(DKD_H5,backed="r"); tech=a.obs["tech"].astype(str).values
    nc=a.obs["nCount_RNA"].astype(float).values; ng=a.obs["nFeature_RNA"].astype(float).values
    samp=a.obs["orig_ident"].astype(str).values; a.file.close()
    fig,axes=plt.subplots(1,3,figsize=(14,4.4)); order=["CosMx","Xenium"]; cols=["#8a6db5","#6A3D9A"]
    qc_panel(axes[0],order,[np.median(nc[tech=='CosMx']),np.median(nc[tech=='Xenium'])],cols,"median transcripts / cell")
    qc_panel(axes[1],order,[np.median(ng[tech=='CosMx']),np.median(ng[tech=='Xenium'])],cols,"median genes / cell")
    nsec=[pd.Series(samp[tech=='CosMx']).nunique(),pd.Series(samp[tech=='Xenium']).nunique()]
    qc_panel(axes[2],order,nsec,cols,"number of sections")
    fig.suptitle("V3 · DKD atlas — both platforms (CosMx n=48 + Xenium n=16). "
                 "Neg-probes dropped from release → ambient not computable",fontsize=13)
    fs.save_fig(fig,"V3_qc")

def V3_typing():
    print("V3_typing (both platforms)...")
    fig,axes=plt.subplots(1,3,figsize=(17,6))
    for k,pl in enumerate(["CosMx","Xenium"]):
        xy2,d2,lab2=_load_native(DKD_H5,MARKERS,pl,sub=120000)
        native_dotplot(axes[k],d2,lab2,MARKERS,DKD_ROWS,f"DKD {pl} native typing",("#8a6db5" if pl=="CosMx" else "#6A3D9A"))
    # umap (Xenium) by native
    a=ad.read_h5ad(DKD_H5,backed="r"); xen=a.obs["tech"].astype(str).values=="Xenium"
    idx=np.where(xen)[0]; idx=rng.choice(idx,45000,replace=False)
    um=np.asarray(a.obsm["X_umap"],float)[idx]; lab=dkd_native(a,"Xenium")[idx]; a.file.close()
    ax=axes[2]; cats=pd.Series(lab).value_counts().index[:18]; cmap=plt.get_cmap("tab20")
    ax.scatter(um[:,0],um[:,1],s=2,c="#eee",linewidths=0,rasterized=True)
    for ci,c in enumerate(cats):
        m=lab==c; ax.scatter(um[m,0],um[m,1],s=3,c=[cmap(ci%20)],linewidths=0,rasterized=True)
    ax.set_title("DKD Xenium — UMAP (native)",color="#6A3D9A",fontsize=12); ax.axis("off")
    fig.suptitle("V3 · DKD native cell typing on BOTH platforms (subtype only trustworthy on Xenium — see A3/A3raw)",fontsize=14)
    fs.save_fig(fig,"V3_typing",tight=False)

def _load_native(h5,markers,plat,sub=None):
    a=ad.read_h5ad(h5,backed="r"); xy=np.asarray(a.obsm["spatial"],float)
    base=a.obs["tech"].astype(str).values==plat; lab=dkd_native(a,plat)
    idx=np.where(base)[0]
    if sub and len(idx)>sub: idx=rng.choice(idx,sub,replace=False); idx.sort()
    gp=[g for g in markers if g in set(map(str,a.var_names))]
    M=a[idx,gp].to_memory(); a.file.close()
    X=M.X.toarray() if sp.issparse(M.X) else np.asarray(M.X)
    return xy[idx],pd.DataFrame(X,columns=gp),lab[idx]

def V3_gallery():
    print("V3_gallery (markers + pops per platform, representative samples)...")
    reps={"CosMx":["HK3631","HK3035","HK3474"],"Xenium":["HK2695","1006","HK3626"]}
    gmark=["MS4A1","MZB1","CD68"]; ncol=len(gmark)+len(DKD_POPS)
    fig,axes=plt.subplots(6,ncol,figsize=(2.4*ncol,2.4*6))
    ri=0
    for pl,sl in reps.items():
        for s in sl:
            def sect(a,pl=pl,s=s): return (a.obs["tech"].astype(str).values==pl)&(a.obs["orig_ident"].astype(str).values==s)
            xy,d,imm=_pop_load(DKD_H5,gmark,sect,sub=45000)
            for ci,g in enumerate(gmark): sp_marker(axes[ri,ci],xy,d[g].values.astype(float) if g in d else np.zeros(len(xy)),g if ri==0 else "")
            for ci,(pn,(labs,pcol)) in enumerate(DKD_POPS.items()): sp_pop(axes[ri,len(gmark)+ci],xy,np.isin(imm,labs),pn if ri==0 else "",pcol)
            axes[ri,0].text(-0.16,0.5,f"{pl}\n{s}",transform=axes[ri,0].transAxes,rotation=90,va="center",ha="center",
                            fontsize=9,fontweight="bold",color="#8a6db5" if pl=="CosMx" else "#6A3D9A")
            ri+=1
    fig.suptitle("V3 · DKD marker & population gallery — both platforms, representative samples (structural finding cross-platform)",fontsize=14)
    fs.save_fig(fig,"V3_gallery")

def _pop_load(h5,markers,sect,sub=None):
    a=ad.read_h5ad(h5,backed="r"); xy=np.asarray(a.obsm["spatial"],float); base=sect(a)
    imm=a.obs["immune_cell_annotation_combined"].astype(str).values; idx=np.where(base)[0]
    if sub and len(idx)>sub: idx=rng.choice(idx,sub,replace=False); idx.sort()
    gp=[g for g in markers if g in set(map(str,a.var_names))]; M=a[idx,gp].to_memory(); a.file.close()
    X=M.X.toarray() if sp.issparse(M.X) else np.asarray(M.X)
    return xy[idx],pd.DataFrame(X,columns=gp),imm[idx]

def V3_aggregates():
    print("V3_aggregates (lineage composition both platforms, apples-to-apples)...")
    base=f"{REPO}/analysis/dkd_cosmx_aggregates"
    cos=pd.read_csv(f"{base}/dkd_cosmx_aggregate_summary.csv").set_index("state")
    xen=pd.read_csv(f"{base}/dkd_xenium_aggregate_summary.csv").set_index("state")
    fig,axes=plt.subplots(1,2,figsize=(13,5),sharey=True)
    states=["B","Plasma","Myeloid","T_lineage","NK"]; x=np.arange(len(states))
    for ax,(df,nm,col,eps) in zip(axes,[(cos,"CosMx","#8a6db5",80),(xen,"Xenium","#6A3D9A",50)]):
        vals=[df.loc[s,"median_log2"] if s in df.index else np.nan for s in states]
        ks=[f"{int(df.loc[s,'k_enriched'])}/{int(df.loc[s,'n_sec'])}" if s in df.index and df.loc[s,'n_sec']>0 else "" for s in states]
        bars=ax.bar(x,[v if v==v else 0 for v in vals],color=col)
        for b,v,k in zip(bars,vals,ks):
            if v==v: ax.text(b.get_x()+b.get_width()/2,v,k,ha="center",va="bottom" if v>=0 else "top",fontsize=9)
        fs.zeroline(ax,0,"h"); ax.set_xticks(x); ax.set_xticklabels(states,rotation=30,ha="right")
        ax.set_ylim(-3,5.2)  # clip CosMx NK single-section pseudocount artifact
        ax.set_title(f"DKD {nm} (lineage; ε={eps})",color=col); ax.set_ylabel("median log₂ enrichment in-aggregate")
    fig.suptitle("V3 · DKD B/plasma aggregate composition recovers on BOTH platforms — structural finding is cross-platform\n"
                 "(lineage only; k/N sections labelled; subtype-resolved Treg/effector-CD8 only on Xenium — see C2/C3)",fontsize=13)
    fs.save_fig(fig,"V3_aggregates")

# ============================================================================
# N — mapping / harmonization explainer
# ============================================================================
def N_map():
    print("N_map (native -> common mapping table)...")
    rowmap=[("B / plasma","Naive B cells; Plasmablasts","B-cell; plasmablast","B; Plasma","B, Plasma"),
     ("T lineage","Eff-CD8; Treg; CD8_T; CCR7+T","Treg; T CD8 mem (low-conf)","CD4+; CD8+","T_lineage (pooled)"),
     ("Myeloid","Inter/Classical mono; myeloid-DC; mregDC","macrophage; monocyte; mDC; pDC","Macro; cDC; pDC","Myeloid"),
     ("NK","Natural killer cells","NK","NK","NK"),
     ("Endothelial","Endothelial","*.endothelium (5 types)","EC_Peritub/glom/DVR/Lymph","Endothelial"),
     ("Stroma","Stroma_mural; Fibroblast","Fibroblast; Myofibroblast; pericyte; mesangial","Fibroblast; VSMC; MC1","Stroma"),
     ("Podocyte","—","Podocyte","Podo","Podocyte"),
     ("Tubular epi.","Proximal_tubule; Intercalated","PCT; CT; TAL; IC; PC ...","PT; TAL; DCT; CNT; IC; PC ...","Tubular_epi"),
     ("Malignant epi.","Tumor_RCC","—","—","Malignant_epi (RCC only)")]
    cols=["Common lineage","ccRCC native","cLN native","DKD native","→ harmonized"]
    fig,ax=plt.subplots(figsize=(15,6)); ax.axis("off")
    cw=[0.13,0.26,0.26,0.23,0.12]; x0=np.concatenate([[0],np.cumsum(cw)]); nrow=len(rowmap); yh=1.0/(nrow+1)
    for j,c in enumerate(cols):
        ax.add_patch(plt.Rectangle((x0[j],1-yh),cw[j],yh,facecolor="#333")); ax.text(x0[j]+cw[j]/2,1-yh/2,c,ha="center",va="center",color="white",fontweight="bold",fontsize=11)
    for i,row in enumerate(rowmap):
        y=1-(i+2)*yh
        for j,val in enumerate(row):
            fc="#eef3f8" if j in (0,4) else "white"
            ax.add_patch(plt.Rectangle((x0[j],y),cw[j],yh,facecolor=fc,edgecolor="#ddd"))
            ax.text(x0[j]+cw[j]/2,y+yh/2,val,ha="center",va="center",fontsize=8.5,fontweight="bold" if j in(0,4) else "normal")
    ax.set_xlim(0,1); ax.set_ylim(1-(nrow+1)*yh,1.0)
    ax.set_title("N · Reconciling three native vocabularies → one common lineage scheme\n"
                 "(done deliberately AFTER each study stands on its own; cLN T pooled & flagged unreliable; RCC epithelium malignant)",fontsize=14)
    fs.save_fig(fig,"N_map",tight=False)

if __name__=="__main__":
    hdr("V1 RCC"); V1_qc(); V1_typing(); V1_gallery(); V1_pattern()
    hdr("V2 cLN"); V2_qc(); V2_typing(); V2_gallery(); V2_niche()
    hdr("V3 DKD"); V3_qc(); V3_typing(); V3_gallery(); V3_aggregates()
    hdr("N mapping"); N_map()
    print("\n== vignettes done ==")
