#!/usr/bin/env python
"""build_support.py — supporting figures A2,B2,B3,B4,C1,C2,C4,C5. Imports figstyle."""
import os, warnings; warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd, anndata as ad
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import figstyle as fs
REPO=fs.REPO
RCC_H5=f"{REPO}/outputs/objects/kidney_RCC_protein.h5ad"; DKD_H5=f"{REPO}/Demoulin26/data/spatial_adata_xenium_cosmx_zenodo.h5ad"
TAB=f"{REPO}/outputs/tables"; CD3=f"{REPO}/analysis/cln_cd3_contamination"; BAFF=f"{REPO}/analysis/dkd_baff_april"
STRESS=f"{REPO}/analysis/dkd_epi_endo_stress"; BN=f"{REPO}/Demoulin26/analysis/bniche_dbscan"; IM=f"{REPO}/analysis/interaction_map"

# ============================================================================
def figA2():
    print("A2 cLN InSituType recall-precision...")
    b=pd.read_csv(f"{TAB}/cln_cosmx_immune_benchmark_unified.csv")
    b=b.sort_values("n_author",ascending=False)
    fig,ax=plt.subplots(figsize=(9,5)); x=np.arange(len(b)); w=0.38
    ax.bar(x-w/2,b.recall_insitutype,w,label="recall",color="#4C72B0")
    ax.bar(x+w/2,b.precision_insitutype,w,label="precision",color="#DD8452")
    for xi,n in zip(x,b.n_author): ax.text(xi,1.02,f"n={int(n)}",ha="center",fontsize=9,color="#666")
    ax.set_xticks(x); ax.set_xticklabels(b.immune_type,rotation=30,ha="right")
    ax.set_ylabel("InSituType vs author labels"); ax.set_ylim(0,1.12)
    ax.axhline(0.5,color="#888",ls="--",lw=1)
    ax.set_title("cLN CosMx immune typing benchmark (InSituType vs author labels)")
    ax.legend(frameon=False,loc="lower right")
    fs.save_fig(fig,"A2")

# ============================================================================
def figB2():
    print("B2 cLN CD3-in-epithelium contamination...")
    s=pd.read_csv(f"{CD3}/cln_cd3_contamination_summary.csv").set_index("compartment")
    order=[c for c in ["author_epithelial","IF_epithelial(PanCK+CD45-)","author_Treg","IF_immune(CD45+PanCK-)"] if c in s.index]
    lab=[o.split("(")[0].replace("author_","").replace("_"," ") for o in order]
    fig,axes=plt.subplots(1,2,figsize=(13,4.8))
    x=np.arange(len(order)); w=0.25
    for i,g in enumerate(["CD3D","CD3E","CD3G"]):
        axes[0].bar(x+(i-1)*w,[s.loc[o,f"{g}_detect"] for o in order],w,label=g)
    axes[0].set_xticks(x); axes[0].set_xticklabels(lab,fontsize=10); axes[0].set_ylabel("detection rate")
    axes[0].set_title("CD3-family detection by compartment"); axes[0].legend(frameon=False,fontsize=10)
    for i,g in enumerate(["CD3D","CD3E","CD3G"]):
        axes[1].bar(x+(i-1)*w,[s.loc[o,f"{g}_mean"] for o in order],w,label=g)
    axes[1].plot(x,[s.loc[o,"negmean"] for o in order],"k--o",lw=1.6,label="negmean (ambient)")
    axes[1].set_xticks(x); axes[1].set_xticklabels(lab,fontsize=10); axes[1].set_ylabel("mean count")
    axes[1].set_title("CD3-family mean count vs negmean floor"); axes[1].legend(frameon=False,fontsize=10)
    epi_fp=s.loc["author_epithelial","any_CD3_detect"]*100
    fig.suptitle(f"cLN ambient CD3 mis-assigned to epithelium: ~{epi_fp:.0f}% epithelial CD3+, "
                 "~2.3× T-vs-epithelial, ~7.4× above negmean",fontsize=14)
    fs.save_fig(fig,"B2")

# ============================================================================
def figB3():
    print("B3 usability-gate strip...")
    g1=pd.read_csv(f"{BAFF}/usability_gate.csv"); g1=g1.assign(
        name=g1.protein,detect=g1.expect_detect,floor=g1.floor_detect,
        ok=g1.verdict.eq("usable"),group="BAFF/APRIL axis")
    g2=pd.read_csv(f"{STRESS}/usability_gate.csv"); g2=g2.assign(
        name=g2.program+":"+g2.gene,detect=g2.target_detect,floor=g2.floor_detect,
        ok=g2.usable,group=g2.program)
    df=pd.concat([g1[["name","detect","floor","ok","group"]],g2[["name","detect","floor","ok","group"]]],ignore_index=True)
    df=df.dropna(subset=["detect"]).sort_values(["ok","detect"])
    fig,ax=plt.subplots(figsize=(9,7.5)); y=np.arange(len(df))
    ax.hlines(y,df.floor,df.detect,color="#bbbbbb",lw=2,zorder=1)
    ax.scatter(df.floor,y,s=40,c="#999999",label="ambient floor",zorder=2,marker="|")
    ax.scatter(df.detect,y,s=70,c=["#2a9d3a" if o else "#c0392b" for o in df.ok],zorder=3,
               edgecolor="black",linewidth=0.3)
    ax.axvline(0.03,color="#888",ls="--",lw=1.2)
    ax.text(0.031,len(df)-0.5,"3% gate",fontsize=10,color="#555")
    ax.set_yticks(y); ax.set_yticklabels(df.name,fontsize=9); ax.set_xlabel("detection rate in expected cells")
    ax.set_title("Usability gate: markers clearing vs failing the ambient floor\n"
                 "(green = usable, red = sub-ambient/fail; BAFF/APRIL & fibroEMT fail)")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([],[],marker="o",ls="",mfc="#2a9d3a",mec="k",ms=9,label="usable"),
                       Line2D([],[],marker="o",ls="",mfc="#c0392b",mec="k",ms=9,label="fail"),
                       Line2D([],[],marker="|",ls="",mfc="#999",mec="#999",ms=11,label="ambient floor")],
              frameon=False,loc="lower right")
    fs.save_fig(fig,"B3")

# ============================================================================
def figB4():
    print("B4 injured-PT near vs far + gradient...")
    nf=pd.read_csv(f"{STRESS}/near_vs_far_summary.csv"); gr=pd.read_csv(f"{STRESS}/distance_gradient.csv")
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    progs=list(nf.program); cols={"injPT":fs.DATASET["DKD"],"endoAct":"#E69F00","hypoxia":"#d62728"}
    x=np.arange(len(progs))
    axes[0].bar(x,nf.pooled_delta,color=[cols.get(p,"#777") for p in progs])
    fs.zeroline(axes[0],0,"h")
    for xi,(_,r) in zip(x,nf.iterrows()):
        axes[0].text(xi,r.pooled_delta+0.01,f"{int(r.k_sections)}/{int(r.n_sections)}",ha="center",fontsize=10)
    axes[0].set_xticks(x); axes[0].set_xticklabels(progs); axes[0].set_ylabel("pooled Δz (near − far)")
    axes[0].set_title("Stress program: near vs far B-aggregate\n(k/N sections positive labelled)")
    blab=["0-50","50-100","100-200","200-500",">500"]; xb=np.arange(len(blab))
    for p in nf.program:
        vals=[gr.loc[(gr.program==p)&(gr.bin==b),"mean_score"].iloc[0] for b in blab]
        axes[1].plot(xb,vals,marker="o",lw=2,color=cols.get(p,"#777"),label=p)
    fs.zeroline(axes[1],0,"h")
    axes[1].set_xticks(xb); axes[1].set_xticklabels(blab); axes[1].set_xlabel("distance to nearest B-aggregate (µm)")
    axes[1].set_ylabel("mean program z-score"); axes[1].set_title("Distance gradient"); axes[1].legend(frameon=False)
    fig.suptitle("DKD injured-PT program is elevated near B-aggregates (Δz +0.13, 6/9 sections, p=0.038)",fontsize=14)
    fs.save_fig(fig,"B4")

# ============================================================================
def _bcrop_section(h5, sect_fn, label_fn, eps=50, mp=20):
    a=ad.read_h5ad(h5,backed="r"); xy=np.asarray(a.obsm["spatial"],float)
    base=sect_fn(a); lab=label_fn(a); a.file.close()
    isB=base&(lab=="B");
    cl=np.full(int(base.sum()),-2)
    Bsub=isB[base]; coords=xy[base]
    if Bsub.sum()>=mp:
        d=DBSCAN(eps=eps,min_samples=mp).fit(coords[Bsub]).labels_
        cl[np.where(Bsub)[0]]=d
    return coords, Bsub, cl

def figC1():
    print("C1 DBSCAN aggregate delineation maps...")
    fig,axes=plt.subplots(1,2,figsize=(15,7.2))
    # RCC
    def rcc_sect(a): return np.ones(a.n_obs,bool)
    def rcc_lab(a): return pd.Series(a.obs["phase_b_label"].astype(str)).map(
        lambda v:"B" if v=="Naive B cells" else "x").values
    co,Bs,cl=_bcrop_section(RCC_H5,rcc_sect,rcc_lab)
    ax=axes[0]; ax.scatter(co[::6,0],co[::6,1],s=1,c="#eeeeee",linewidths=0,rasterized=True)
    ax.scatter(co[Bs][cl[Bs]==-1,0],co[Bs][cl[Bs]==-1,1],s=4,c="#9ecae1",linewidths=0,rasterized=True)
    ax.scatter(co[Bs][cl[Bs]>=0,0],co[Bs][cl[Bs]>=0,1],s=8,c=fs.DATASET["RCC"],linewidths=0,rasterized=True)
    ax.set_title(f"ccRCC: {len(set(cl[cl>=0]))} B-aggregates",color=fs.DATASET["RCC"]); ax.set_aspect("equal"); ax.axis("off")
    # DKD HK2695
    def dkd_sect(a): return (a.obs["orig_ident"].astype(str).values=="HK2695")&(a.obs["tech"].astype(str).values=="Xenium")
    def dkd_lab(a): return a.obs["immune_cell_annotation_combined"].astype(str).values
    co,Bs,cl=_bcrop_section(DKD_H5,dkd_sect,dkd_lab)
    ax=axes[1]; ax.scatter(co[::6,0],co[::6,1],s=1,c="#eeeeee",linewidths=0,rasterized=True)
    ax.scatter(co[Bs][cl[Bs]==-1,0],co[Bs][cl[Bs]==-1,1],s=4,c="#c0a5d8",linewidths=0,rasterized=True)
    ax.scatter(co[Bs][cl[Bs]>=0,0],co[Bs][cl[Bs]>=0,1],s=8,c=fs.DATASET["DKD"],linewidths=0,rasterized=True)
    ax.set_title(f"DKD (HK2695): {len(set(cl[cl>=0]))} B-aggregates",color=fs.DATASET["DKD"]); ax.set_aspect("equal"); ax.axis("off")
    fig.suptitle("B-cell density aggregates (DBSCAN ε=50 µm, minPts=20): light = dispersed B, dark = aggregated B",fontsize=14)
    fs.save_fig(fig,"C1")

# ============================================================================
def figC2():
    print("C2 per-aggregate Treg vs eff-CD8 swarm...")
    d=pd.read_csv(f"{BN}/dkd_xenium_aggregates.csv")
    dt=np.log2((d["f_Treg-like"]+1e-6)/(d["bg_Treg-like"]+1e-6)); de=np.log2((d["f_eff-CD8"]+1e-6)/(d["bg_eff-CD8"]+1e-6))
    r=pd.read_csv(f"{TAB}/rcc_phaseB2_aggregates.csv"); rc=pd.read_csv(f"{TAB}/rcc_phaseB2_aggregate_composition.csv").set_index("cell_type")
    rt=np.log2((r["f_Treg"]+1e-6)/(rc.loc["Treg","background_frac"]+1e-6)); re=np.log2((r["f_eff-CD8"]+1e-6)/(rc.loc["eff-CD8","background_frac"]+1e-6))
    fig,ax=plt.subplots(figsize=(9,5.2))
    data=[("Treg-like",0,rt,dt),("effector-CD8",1,re,de)]
    for name,yb,rv,dv in data:
        ax.scatter(rv,np.full(len(rv),yb+0.18)+np.random.uniform(-0.06,0.06,len(rv)),s=30,
                   color=fs.DATASET["RCC"],alpha=0.7,edgecolor="k",linewidth=0.2,label="ccRCC" if yb==0 else None)
        ax.scatter(dv,np.full(len(dv),yb-0.18)+np.random.uniform(-0.06,0.06,len(dv)),s=30,
                   color=fs.DATASET["DKD"],alpha=0.7,edgecolor="k",linewidth=0.2,label="DKD" if yb==0 else None)
        ax.scatter(np.median(rv),yb+0.18,marker="|",s=400,c="k"); ax.scatter(np.median(dv),yb-0.18,marker="|",s=400,c="k")
    fs.zeroline(ax,0,"v")
    ax.set_yticks([0,1]); ax.set_yticklabels(["Treg-like","effector-CD8"])
    ax.set_xlabel("per-aggregate log₂ enrichment (inside vs section background)")
    ax.set_title("Per-aggregate Treg vs cytotoxic-CD8 enrichment: ccRCC vs DKD")
    ax.legend(frameon=False,loc="lower right")
    fs.save_fig(fig,"C2")

# ============================================================================
def figC4():
    print("C4 count-pooled radial profile...")
    rad=pd.read_csv(f"{BN}/dkd_xenium_radial.csv"); rings=["core 0-50","margin 50-100","outer 100-150"]
    fig,ax=plt.subplots(figsize=(8.5,5)); xb=np.arange(len(rings))
    for state,col,mk in [("Treg-like",fs.DATASET["DKD"],"o"),("eff-CD8","#E69F00","s")]:
        ys=[]
        for rn in rings:
            sub=rad[(rad.state==state)&(rad.ring==rn)]
            ins=sub.state_cells.sum(); exp=(sub.n_cells*sub.bg).sum()
            ys.append(np.log2((ins+1e-9)/(exp+1e-9)) if exp>0 else np.nan)
        ax.plot(xb,ys,marker=mk,lw=2.5,ms=10,color=col,label=state)
    fs.zeroline(ax,0,"h")
    ax.set_xticks(xb); ax.set_xticklabels(rings); ax.set_xlabel("ring (distance from B-core centroid, µm)")
    ax.set_ylabel("count-pooled log₂ enrichment")
    ax.set_title("DKD radial profile (count-pooled): Treg flat across rings,\nmild cytotoxic-core gradient — no Treg collar")
    ax.legend(frameon=False)
    fs.save_fig(fig,"C4")

# ============================================================================
def figC5():
    print("C5 comparative interaction heatmaps...")
    order_ni=["Tubular_epi","Endothelial","Stroma","Podocyte"]; imm=["B","Plasma","Myeloid","T_lineage","NK"]
    dsets=["RCC","cLN","DKD"]
    fig,axes=plt.subplots(1,len(dsets),figsize=(5*len(dsets),4.8),sharey=True)
    for ax,dn in zip(axes,dsets):
        z=pd.read_csv(f"{IM}/zmatrix_{dn}.csv")
        mat=np.full((len(order_ni),len(imm)),np.nan)
        for i,r in enumerate(order_ni):
            for j,c in enumerate(imm):
                s=z[(z.non_immune==r)&(z.immune==c)]
                if len(s) and s.n_sec.iloc[0]>0: mat[i,j]=s.mean_z.iloc[0]
        im=ax.imshow(mat,cmap="RdBu_r",vmin=-8,vmax=8,aspect="auto")
        ax.set_xticks(range(len(imm))); ax.set_xticklabels([fs.LINEAGE_LABEL.get(c,c) for c in imm],rotation=45,ha="right",fontsize=10)
        ax.set_yticks(range(len(order_ni))); ax.set_yticklabels([fs.LINEAGE_LABEL.get(r,r) for r in order_ni],fontsize=10)
        for i in range(len(order_ni)):
            for j in range(len(imm)):
                if mat[i,j]==mat[i,j]: ax.text(j,i,f"{mat[i,j]:.0f}",ha="center",va="center",fontsize=8,
                                               color="white" if abs(mat[i,j])>4 else "black")
        ax.set_title(fs.DATASET_LONG[dn],color=fs.DATASET[dn],fontsize=13)
    fig.colorbar(im,ax=axes,fraction=0.02,label="mean neighborhood z")
    fig.suptitle("Comparative non-immune × immune spatial neighborhood enrichment\n"
                 "INSET: absolute z = geometry (immune cells aggregate, parenchyma 'avoids'); biology = differentials",
                 fontsize=14)
    # inset note
    axes[0].text(-0.02,-0.32,"Note: ccRCC stroma–immune inversion is PROVISIONAL (not tile-verified) — interpret with caution.",
                 transform=axes[0].transAxes,fontsize=9,color="#b2182b")
    fs.save_fig(fig,"C5",tight=False)

if __name__=="__main__":
    figA2(); figB2(); figB3(); figB4(); figC1(); figC2(); figC4(); figC5()
    print("== supporting done ==")
