#!/usr/bin/env python
"""build_hero.py — hero six presentation figures (A1, A3, B1, C3, D1, D2). Imports figstyle."""
import os, warnings; warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
import figstyle as fs
rng=np.random.default_rng(0)
REPO=fs.REPO
RCC_H5=f"{REPO}/outputs/objects/kidney_RCC_protein.h5ad"
CLN_H5=f"{REPO}/outputs/objects/cln_cosmx.h5ad"
DKD_H5=f"{REPO}/Demoulin26/data/spatial_adata_xenium_cosmx_zenodo.h5ad"
XT=f"{REPO}/Demoulin26/analysis/cross_platform_tcell"
BN=f"{REPO}/Demoulin26/analysis/bniche_dbscan"

IMM_SET={"B","Plasma","T_lineage","Myeloid","NK"}
def scatter_lineage(ax, xy, lab, title, present_order, s=2.0):
    # parenchyma first, immune ON TOP (bigger) so the interesting cells pop; malignant starred
    ax.scatter(xy[:,0],xy[:,1],s=s*0.5,c="#f0f0f0",linewidths=0,rasterized=True)
    draw=[k for k in present_order if k not in IMM_SET]+[k for k in present_order if k in IMM_SET]
    for k in draw:
        m=lab==k
        if not m.any(): continue
        if k=="Malignant_epi":
            ax.scatter(xy[m,0],xy[m,1],s=s*2.0,c=fs.LINEAGE[k],marker="*",
                       edgecolor="black",linewidth=0.15,rasterized=True)
        else:
            ax.scatter(xy[m,0],xy[m,1],s=(s*2.2 if k in IMM_SET else s),
                       c=fs.LINEAGE[k],linewidths=0,rasterized=True)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title)

def densest_crop(xy, half):
    """crop to a (2*half)-wide window around the densest region (for multi-core CosMx slides)."""
    H,xe,ye=np.histogram2d(xy[:,0],xy[:,1],bins=24)
    i,j=np.unravel_index(np.argmax(H),H.shape)
    cx=(xe[i]+xe[i+1])/2; cy=(ye[j]+ye[j+1])/2
    return (np.abs(xy[:,0]-cx)<half)&(np.abs(xy[:,1]-cy)<half)

# ============================================================================
def figA1():
    print("A1 spatial cell-type maps...")
    panels=[]
    # RCC (single section)
    r=ad.read_h5ad(RCC_H5)
    lab=pd.Series(r.obs["phase_b_label"].astype(str)).map(fs.RCC_LINEAGE).values
    xy=np.asarray(r.obsm["spatial"],float); m=pd.notna(lab)
    panels.append(("RCC",xy[m],lab[m].astype(str))); del r
    # cLN (largest section)
    c=ad.read_h5ad(CLN_H5); samp=c.obs["sample"].astype(str).values
    big=pd.Series(samp).value_counts().index[0]
    lab=pd.Series(c.obs["author_celltype"].astype(str)).map(fs.CLN_LINEAGE).values
    xy=np.asarray(c.obsm["spatial"],float); m=pd.notna(lab)&(samp==big)
    cxy=xy[m]; clab=lab[m].astype(str); crop=densest_crop(cxy,0.85)  # mm units; crop to densest core
    panels.append(("cLN",cxy[crop],clab[crop])); del c
    # DKD (HK2695 Xenium section)
    d=ad.read_h5ad(DKD_H5,backed="r"); xen=d.obs["tech"].astype(str).values=="Xenium"
    samp=d.obs["orig_ident"].astype(str).values; sel_s="HK2695"
    ann=d.obs["annotation_updated"].astype(str).values; imm=d.obs["immune_cell_annotation_combined"].astype(str).values
    xy=np.asarray(d.obsm["spatial"],float); d.file.close()
    lab=np.array([fs.DKD_IMMUNE.get(i) if i!="Unknown" else fs.DKD_NONIMM.get(a)
                  for a,i in zip(ann,imm)],dtype=object)
    m=xen&(samp==sel_s)&np.array([x is not None for x in lab])
    panels.append(("DKD",xy[m],lab[m].astype(str)))
    present=[k for k in fs.LINEAGE_ORDER]
    fig,axes=plt.subplots(1,3,figsize=(18,6.6))
    for ax,(nm,xy,lab) in zip(axes,panels):
        if len(lab)>90000:
            idx=rng.choice(len(lab),90000,replace=False); xy=xy[idx]; lab=lab[idx]
        po=[k for k in present if (lab==k).any()]
        scatter_lineage(ax,xy,lab,fs.DATASET_LONG[nm],po,s=2.4)
        ax.set_title(fs.DATASET_LONG[nm],color=fs.DATASET[nm])
    fs.lineage_legend(axes[2],[k for k in present],bbox=(1.02,0.5))
    fig.suptitle("Spatial cell-type maps across three kidney contexts")
    fs.save_fig(fig,"A1")

# ============================================================================
def figA3():
    print("A3 measured-vs-imputed...")
    sup=pd.read_csv(f"{XT}/results_marker_support.csv"); auc=pd.read_csv(f"{XT}/results_auroc.csv")
    plats=["Xenium","CosMx"]; markers=["CD8A","CD4"]
    fig,axes=plt.subplots(1,2,figsize=(12,5),sharey=True)
    for ax,plat in zip(axes,plats):
        x=np.arange(len(markers)); w=0.36
        for i,(grp,col) in enumerate([("CD4","#56B4E9"),("CD8","#D55E00")]):
            vals=[sup[(sup.platform==plat)&(sup.label==grp)&(sup.marker==mk)].det_rate.iloc[0] for mk in markers]
            ax.bar(x+(i-0.5)*w,vals,w,label=f"{grp}-labeled",color=col)
        amb=[sup[(sup.platform==plat)&(sup.label=="CD4")&(sup.marker==mk)].ambient_det_rate.iloc[0] for mk in markers]
        for xi,a in zip(x,amb): ax.hlines(a,xi-0.45,xi+0.45,color="k",ls=":",lw=1.4)
        for xi,mk in zip(x,markers):
            a=auc[(auc.platform==plat)&(auc.marker==mk)].auroc_cd8_vs_cd4
            if len(a): ax.text(xi,ax.get_ylim()[1]*0.92 if ax.get_ylim()[1]>0 else 0.4,
                               f"AUROC\n{a.iloc[0]:.2f}",ha="center",fontsize=11,
                               color="#222",fontweight="bold")
        ax.set_xticks(x); ax.set_xticklabels(markers); ax.set_title(plat)
        ax.set_ylabel("detection rate (≥1 count)") if plat=="Xenium" else None
    axes[0].legend(frameon=False,loc="upper right")
    axes[1].text(0.5,0.78,"CD4 AUROC ≈ 0.50\n-> subtype is IMPUTED",transform=axes[1].transAxes,
                 ha="center",color="#b2182b",fontsize=12,fontweight="bold")
    fig.suptitle("CD4/CD8 subtype is measured on Xenium, imputed on CosMx\n"
                 "(discriminating-marker detection in CD4- vs CD8-labelled cells; dotted = ambient floor)",
                 fontsize=15)
    fs.save_fig(fig,"A3")

# ============================================================================
def _markers_in_crop(h5, markers, cx, cy, W, label_col, bval, section_mask_fn=None):
    a=ad.read_h5ad(h5,backed="r"); xy=np.asarray(a.obsm["spatial"],float)
    base=np.ones(a.n_obs,bool) if section_mask_fn is None else section_mask_fn(a)
    lab=a.obs[label_col].astype(str).values
    m=base&(np.abs(xy[:,0]-cx)<W)&(np.abs(xy[:,1]-cy)<W)
    idx=np.where(m)[0]; present=[g for g in markers if g in set(map(str,a.var_names))]
    sub=a[idx,present].to_memory(); a.file.close()
    X=sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
    return xy[idx], pd.DataFrame(X,columns=present), (lab[idx]==bval)

def figB1():
    print("B1 aggregate marker overlays (impact figure)...")
    markers=["MS4A1","MZB1","FOXP3","CD8A","GZMB"]; W=95.0
    ragg=pd.read_csv(f"{REPO}/outputs/tables/rcc_phaseB2_aggregates.csv").sort_values("n_B",ascending=False).iloc[0]
    rxy,rmk,rB=_markers_in_crop(RCC_H5,markers,ragg.x,ragg.y,W,"phase_b_label","Naive B cells")
    dagg=pd.read_csv(f"{BN}/dkd_xenium_aggregates.csv").sort_values("n_B",ascending=False).iloc[0]
    dsec=str(dagg["sample"])
    def dkd_sec(a): return (a.obs["orig_ident"].astype(str).values==dsec)&(a.obs["tech"].astype(str).values=="Xenium")
    dxy,dmk,dB=_markers_in_crop(DKD_H5,markers,dagg.x,dagg.y,W,"immune_cell_annotation_combined","B",dkd_sec)
    ncol=len(markers)+1
    fig,axes=plt.subplots(2,ncol,figsize=(3.4*ncol,7.4))
    for row,(xy,mk,Bm,tag,dcol) in enumerate([(rxy,rmk,rB,"ccRCC aggregate",fs.DATASET["RCC"]),
                                              (dxy,dmk,dB,"DKD aggregate",fs.DATASET["DKD"])]):
        for col,g in enumerate(markers):
            ax=axes[row,col]
            ax.scatter(xy[:,0],xy[:,1],s=5,c="#ededed",linewidths=0,rasterized=True)
            ax.scatter(xy[Bm,0],xy[Bm,1],s=10,c="#bcd6ee",linewidths=0,rasterized=True)  # B-core reference
            if g in mk.columns:
                v=mk[g].values.astype(float); pos=v>0
                vmax=max(np.percentile(v[pos],98),1) if pos.any() else 1
                ax.scatter(xy[pos,0],xy[pos,1],s=34,c=np.clip(v[pos],0,vmax),cmap=fs.MARKER_CMAP.get(g,"magma"),
                           vmin=0,vmax=vmax,linewidths=0.2,edgecolor="white",rasterized=True)
            ax.set_aspect("equal"); ax.axis("off")
            if row==0: ax.set_title(f"{g}\n{fs.MARKER_LINEAGE.get(g,'')}",fontsize=13)
        # composite column: B (blue) + Treg/FOXP3+ (red) + cytotoxic/CD8A+ (green)
        ax=axes[row,ncol-1]
        ax.scatter(xy[:,0],xy[:,1],s=5,c="#ededed",linewidths=0,rasterized=True)
        ax.scatter(xy[Bm,0],xy[Bm,1],s=14,c="#1f77b4",linewidths=0,alpha=0.8,rasterized=True)
        if "FOXP3" in mk: fp=mk["FOXP3"].values>0; ax.scatter(xy[fp,0],xy[fp,1],s=40,c="#d62728",edgecolor="white",linewidth=0.3,rasterized=True)
        if "CD8A" in mk: c8=mk["CD8A"].values>0; ax.scatter(xy[c8,0],xy[c8,1],s=40,c="#2ca02c",edgecolor="white",linewidth=0.3,marker="^",rasterized=True)
        ax.set_aspect("equal"); ax.axis("off")
        if row==0: ax.set_title("Composite\nB / Treg / CD8",fontsize=13)
        axes[row,0].text(-0.10,0.5,tag,transform=axes[row,0].transAxes,rotation=90,
                         va="center",ha="center",fontsize=15,fontweight="bold",color=dcol)
    # composite legend
    from matplotlib.lines import Line2D
    axes[1,ncol-1].legend(handles=[Line2D([],[],marker="o",ls="",mfc="#1f77b4",mec="none",ms=8,label="B cell"),
        Line2D([],[],marker="o",ls="",mfc="#d62728",mec="white",ms=8,label="FOXP3+ (Treg)"),
        Line2D([],[],marker="^",ls="",mfc="#2ca02c",mec="white",ms=8,label="CD8A+ (cytotoxic)")],
        loc="upper center",bbox_to_anchor=(0.5,-0.02),frameon=False,fontsize=9,ncol=1)
    fig.suptitle("B/plasma aggregates: Treg-around & cytotoxic-excluded (ccRCC) vs cytotoxic-mixed-in (DKD)\n"
                 "blue underlay = B-core; marker+ cells overlaid (single representative aggregate per context, ~190 µm field)",fontsize=15)
    fs.save_fig(fig,"B1")

# ============================================================================
def figC3():
    print("C3 burden-corrected differential (HEADLINE)...")
    d=pd.read_csv(f"{BN}/differential_treg_vs_cd8.csv")
    fig,ax=plt.subplots(figsize=(9,4.6))
    ymap={"RCC (tumor)":1,"DKD (kidney)":0}
    for _,r in d.iterrows():
        nm="RCC" if "RCC" in r.cohort else "DKD"; y=ymap[r.cohort]; col=fs.DATASET[nm]
        ax.errorbar(r.delta_log2,y,xerr=[[r.delta_log2-r.ci_lo],[r.ci_hi-r.delta_log2]],
                    fmt="o",ms=15,color=col,capsize=8,lw=3,markeredgecolor="black",zorder=3)
        ax.text(r.delta_log2,y+0.22,f"Δ = {r.delta_log2:+.2f}  (~{r.fold_bias:.0f}× Treg bias)",
                ha="center",fontsize=13,color=col,fontweight="bold")
    fs.zeroline(ax,0,"v",label="no bias")
    ax.text(0.02,-0.46,"no Treg-vs-cytotoxic bias",fontsize=10,color="#666",rotation=0)
    ax.text(0.5,1.62,"CIs non-overlapping -> immunoregulatory bias is TUMOR-SPECIFIC",
            transform=ax.transAxes,ha="center",fontsize=13,fontweight="bold",color="#222")
    ax.set_yticks([0,1]); ax.set_yticklabels(["DKD (kidney)","ccRCC (tumor)"]); ax.set_ylim(-0.55,1.5)
    ax.set_xlabel("differential enrichment  Δlog₂ = log₂(Treg) − log₂(effector-CD8)\n(count-pooled, bootstrap 95% CI; burden-corrected)")
    ax.set_title("Treg-over-cytotoxic bias in B-aggregates: tumor vs non-malignant kidney")
    fs.save_fig(fig,"C3")

# ============================================================================
def figD1():
    print("D1 conserved-scaffold / context-wiring schematic...")
    fig,axes=plt.subplots(1,3,figsize=(16.5,6))
    contexts=[("RCC","ccRCC: Treg ring, cytotoxic EXCLUDED",
               [("Treg","#d62728","ring"),("CD8","#2ca02c","excluded")]),
              ("cLN","cLN: + myeloid, plasma–myeloid niche",
               [("Myeloid","#9467bd","ring"),("Plasma","#ff7f0e","mixed")]),
              ("DKD","DKD: cytotoxic MIXED-IN + injured tubule",
               [("CD8","#2ca02c","mixed"),("Treg","#d62728","mixed"),("injured PT","#8c510a","adjacent")])]
    for ax,(nm,title,surround) in zip(axes,contexts):
        ax.set_xlim(-1.4,1.4); ax.set_ylim(-1.4,1.5); ax.set_aspect("equal"); ax.axis("off")
        # shared B/plasma core
        ax.add_patch(Circle((0,0),0.55,facecolor="#1f77b4",alpha=0.25,edgecolor="#1f77b4",lw=2))
        ax.text(0,0,"B / plasma\ncore",ha="center",va="center",fontsize=11,fontweight="bold",color="#10406b")
        for name,col,mode in surround:
            if mode=="ring":
                th=np.linspace(0,2*np.pi,22)
                ax.scatter(0.85*np.cos(th),0.85*np.sin(th),s=70,c=col,edgecolor="white",lw=0.5,zorder=4)
                ax.text(0,1.05,name,ha="center",color=col,fontsize=12,fontweight="bold")
            elif mode=="excluded":
                th=np.linspace(0,2*np.pi,16)
                ax.scatter(1.25*np.cos(th),1.25*np.sin(th),s=45,c=col,alpha=0.55,zorder=2)
                ax.text(0,-1.28,f"{name} excluded",ha="center",color=col,fontsize=12,fontweight="bold")
                ax.add_patch(Circle((0,0),0.95,facecolor="none",edgecolor=col,ls=":",lw=1.8))
            elif mode=="mixed":
                pts=rng.uniform(-0.5,0.5,(14,2))
                ax.scatter(pts[:,0],pts[:,1],s=42,c=col,alpha=0.75,zorder=5,edgecolor="white",lw=0.4)
            elif mode=="adjacent":
                ax.add_patch(Rectangle((-1.3,-1.45),2.6,0.32,facecolor=col,alpha=0.5))
                ax.text(0,-1.29,name,ha="center",color="white",fontsize=11,fontweight="bold")
        ax.set_title(title,fontsize=13,color=fs.DATASET[nm])
    fig.suptitle("Conserved B/plasma scaffold, context-specific immune wiring",fontsize=18)
    fs.save_fig(fig,"D1")
    print("    NOTE: D1 is a programmatic draft — may want manual vector polish (Illustrator/Inkscape).")

# ============================================================================
def figD2():
    print("D2 platform capability matrix...")
    rows=[("T-lineage (CD3 family)","measured","measured","Anchored by CD3D/E/G on both panels"),
          ("CD4/CD8 subtype","measured","IMPUTED (AUROC≈0.5)","Xenium 5k resolves; CosMx 1k reference-imputed"),
          ("Low-abundance ligands\nBAFF / APRIL","sub-ambient","sub-ambient","Below ambient floor on both panels"),
          ("Receptors\nBCMA / BAFF-R / TACI","specific","specific","343–756× ambient; cell-type specific"),
          ("Ambient / segmentation","present","present (~35% epi CD3+)","cLN: ambient CD3 mis-assigned to epithelium"),
          ("IF anchors (PanCK/CD45)","—","decisive","Anchored lineages recoverable; unanchored unfalsifiable")]
    cols=["Capability","Xenium 5k","CosMx 1k","Note"]
    def cellcolor(t):
        t=t.lower()
        if "measured" in t or "specific" in t or "decisive" in t: return "#cdeccd"
        if "imput" in t or "sub-ambient" in t or "35%" in t: return "#f6cccc"
        if "present" in t: return "#fdebcd"
        return "white"
    fig,ax=plt.subplots(figsize=(13.5,5.4)); ax.axis("off")
    ncol=len(cols); nrow=len(rows)
    cw=[0.27,0.16,0.22,0.35]; x0=np.concatenate([[0],np.cumsum(cw)])
    yh=1.0/(nrow+1)
    for j,c in enumerate(cols):
        ax.add_patch(Rectangle((x0[j],1-yh),cw[j],yh,facecolor="#333333"))
        ax.text(x0[j]+cw[j]/2,1-yh/2,c,ha="center",va="center",color="white",fontweight="bold",fontsize=12)
    for i,row in enumerate(rows):
        y=1-(i+2)*yh
        for j,val in enumerate(row):
            fcol=cellcolor(val) if j in (1,2) else ("#f3f3f3" if j==0 else "white")
            ax.add_patch(Rectangle((x0[j],y),cw[j],yh,facecolor=fcol,edgecolor="#cccccc"))
            ax.text(x0[j]+cw[j]/2,y+yh/2,val,ha="center",va="center",fontsize=10.5,
                    fontweight="bold" if j==0 else "normal",wrap=True)
    ax.set_xlim(0,1); ax.set_ylim(1-(nrow+1)*yh-0.02,1.01)
    ax.set_title("Platform capability matrix: what each panel can and cannot establish",fontsize=16,pad=12)
    fs.save_fig(fig,"D2",tight=False)

if __name__=="__main__":
    figA1(); figA3(); figB1(); figC3(); figD1(); figD2()
    print("== hero six done ==")
