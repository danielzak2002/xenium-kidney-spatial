#!/usr/bin/env python
"""STEP 5 — presentable deliverables. (1) 4x4 spatial cell-type grid of the 16 Xenium samples,
shared color key, sample id + disease status, region-cropped. (2) integration UMAPs (by sample /
by my cell type / by lineage). (3) concordance: yours-vs-theirs UMAP + confusion heatmaps.
Reuses figstyle.py rcParams; saves PNG@300 + SVG into this analysis folder."""
import os, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys; sys.path.insert(0,"/Users/danie/ClaudeCode/pilot_analyses/xenium/presentation/figures")
import figstyle as fs   # applies slide-grade rcParams on import
from sklearn.metrics import adjusted_rand_score
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
OUT=f"{REPO}/analysis/dkd_xenium_reannotation"; FIG=f"{OUT}/figures"; os.makedirs(FIG,exist_ok=True)
def save(fig,name):
    fig.savefig(f"{FIG}/{name}.png",dpi=300,bbox_inches="tight"); fig.savefig(f"{FIG}/{name}.svg",bbox_inches="tight")
    plt.close(fig); print(f"  [ok] {name}")

cells=pd.read_parquet(f"{OUT}/cells.parquet")
strata=pd.read_csv(f"{OUT}/disease_strata_per_sample.csv").set_index("orig_ident")
cells["orig_ident"]=cells.orig_ident.astype(str)

# ---- display cell-type scheme: muted epithelium, saturated immune (DKD B-cell story) ----
DISPLAY={"PT":"PT","iPT":"PT","TAL":"TAL","iTAL":"TAL","DTL_ATL":"DTL/ATL","DCT":"DCT/CNT","CNT":"DCT/CNT",
 "PC":"PC/IC","IC A":"PC/IC","IC B":"PC/IC","Podo":"Glomerular","PEC":"Glomerular",
 "EC_glom":"Endothelial","EC_Peritub":"Endothelial","EC_DVR":"Endothelial","EC_Lymph":"Endothelial",
 "Fibroblast":"Fibroblast","VSMC":"VSMC/MC","MC":"VSMC/MC",
 "B":"B cell","Plasma":"Plasma","CD4 T":"T cell","CD8 T":"T cell","Treg":"T cell",
 "Myeloid":"Myeloid","NK":"NK/DC","DC":"NK/DC","Mast_Baso":"NK/DC","Neutrophil":"NK/DC"}
DCOL={"PT":"#E8D9B5","TAL":"#C9B27C","DCT/CNT":"#A8863B","PC/IC":"#8C6D1F","DTL/ATL":"#D6C2A8",
 "Glomerular":"#B5651D","Endothelial":"#8c564b","Fibroblast":"#e377c2","VSMC/MC":"#b89cc4",
 "B cell":"#1f77b4","Plasma":"#ff7f0e","T cell":"#2ca02c","Myeloid":"#9467bd","NK/DC":"#17becf"}
DORDER=["PT","TAL","DCT/CNT","PC/IC","DTL/ATL","Glomerular","Endothelial","Fibroblast","VSMC/MC",
        "B cell","Plasma","T cell","Myeloid","NK/DC"]
IMMUNE_DISP={"B cell","Plasma","T cell","Myeloid","NK/DC"}
cells["disp"]=cells.my_label.map(DISPLAY).fillna("PT")

# ============ FIG 1 — 4x4 spatial cell-type grid ============
def spatial_grid():
    order=strata.sort_values(["Condition","orig_ident"]).index.tolist()
    fig,axes=plt.subplots(4,4,figsize=(17,17.6))
    for ax,sid in zip(axes.ravel(),order):
        g=cells[cells.orig_ident==sid]; x=g.spatial_x.values; y=g.spatial_y.values
        # region crop: drop stray 0.3-99.7 pct whitespace
        xlo,xhi=np.percentile(x,[0.3,99.7]); ylo,yhi=np.percentile(y,[0.3,99.7])
        epi=~g.disp.isin(IMMUNE_DISP)
        ax.scatter(x[epi],y[epi],s=0.6,c=g.disp[epi].map(DCOL),linewidths=0,rasterized=True)
        im=g.disp.isin(IMMUNE_DISP)            # immune on top, slightly larger
        ax.scatter(x[im],y[im],s=3.0,c=g.disp[im].map(DCOL),linewidths=0,rasterized=True)
        ax.set_xlim(xlo,xhi); ax.set_ylim(yhi,ylo); ax.set_aspect("equal"); ax.axis("off")
        cond=strata.loc[sid,"Condition"]; n=int(strata.loc[sid,"n_cells"])
        nB=int((g.disp=="B cell").sum())
        ax.set_title(f"{sid} · {cond}\n{n:,} cells · {nB:,} B",fontsize=11,fontweight="bold",
                     color=fs.DATASET["DKD"] if cond=="DKD" else ("#555" if cond=="Control" else "#888"))
    handles=[Line2D([],[],marker="o",ls="",mfc=DCOL[k],mec="none",ms=9,label=k) for k in DORDER]
    fig.legend(handles=handles,loc="lower center",ncol=7,frameon=False,fontsize=12,bbox_to_anchor=(0.5,-0.012))
    fig.suptitle("DKD Xenium — independent cell-type maps, 16 samples (region-cropped, shared key; immune enlarged)",fontsize=18,y=0.995)
    fig.subplots_adjust(top=0.95,bottom=0.05,wspace=0.04,hspace=0.12)
    save(fig,"spatial_grid_16")

# ============ FIG 2 — integration UMAPs ============
def umaps():
    fig,axes=plt.subplots(1,3,figsize=(21,7))
    # by sample
    ax=axes[0]; sids=sorted(cells.orig_ident.unique()); cmap=plt.cm.tab20(np.linspace(0,1,len(sids)))
    for c,sid in zip(cmap,sids):
        m=cells.orig_ident==sid; ax.scatter(cells.umap_x[m],cells.umap_y[m],s=0.4,color=c,linewidths=0,rasterized=True)
    ax.set_title("by sample (16 — Harmony-integrated)"); ax.axis("off")
    ax.legend(handles=[Line2D([],[],marker="o",ls="",mfc=c,mec="none",ms=6,label=s) for c,s in zip(cmap,sids)],
              loc="center left",bbox_to_anchor=(1.0,0.5),frameon=False,fontsize=7,ncol=1)
    # by my cell type
    ax=axes[1]
    for k in DORDER:
        m=cells.disp==k; ax.scatter(cells.umap_x[m],cells.umap_y[m],s=0.4,color=DCOL[k],linewidths=0,rasterized=True)
    ax.set_title("by independent cell type"); ax.axis("off")
    ax.legend(handles=[Line2D([],[],marker="o",ls="",mfc=DCOL[k],mec="none",ms=6,label=k) for k in DORDER],
              loc="center left",bbox_to_anchor=(1.0,0.5),frameon=False,fontsize=8)
    # by lineage
    ax=axes[2]; LC={"Epithelial":"#c7c7c7","Endothelial":"#8c564b","Stroma":"#e377c2","Immune":"#d62728"}
    for k,c in LC.items():
        m=cells.my_lineage==k; ax.scatter(cells.umap_x[m],cells.umap_y[m],s=0.4,color=c,linewidths=0,rasterized=True)
    ax.set_title("by lineage (immune highlighted)"); ax.axis("off")
    ax.legend(handles=[Line2D([],[],marker="o",ls="",mfc=c,mec="none",ms=7,label=k) for k,c in LC.items()],
              loc="center left",bbox_to_anchor=(1.0,0.5),frameon=False,fontsize=9)
    fig.suptitle("DKD Xenium integration UMAP (Harmony on sample_id; 951,040 cells)",fontsize=17)
    save(fig,"integration_umaps")

# ============ FIG 3 — yours-vs-theirs UMAP ============
def umap_concordance():
    sub=cells.sample(min(200000,len(cells)),random_state=0)
    fig,axes=plt.subplots(1,2,figsize=(15,7.2))
    for ax,(col,ttl) in zip(axes,[("my_coarse","MY labels (independent)"),("auth_coarse","AUTHOR annotation_updated")]):
        cats=sorted(pd.unique(cells[col].dropna()))
        cm=plt.cm.tab20(np.linspace(0,1,len(cats))); cd=dict(zip(cats,cm))
        for k in cats:
            m=sub[col]==k; ax.scatter(sub.umap_x[m],sub.umap_y[m],s=0.5,color=cd[k],linewidths=0,rasterized=True)
        ax.set_title(ttl); ax.axis("off")
        if col=="auth_coarse":
            ax.legend(handles=[Line2D([],[],marker="o",ls="",mfc=cd[k],mec="none",ms=6,label=k) for k in cats],
                      loc="center left",bbox_to_anchor=(1.0,0.5),frameon=False,fontsize=8,ncol=1)
    fig.suptitle("Same integrated UMAP, my labels vs authors' — spatial concordance of the two independent typings",fontsize=15)
    save(fig,"umap_yours_vs_theirs")

# ============ FIG 4 — concordance heatmaps ============
def concordance_heatmaps():
    vs=pd.read_csv(f"{OUT}/validation_summary.csv").set_index("comparison")
    fig,axes=plt.subplots(1,2,figsize=(17,7.6))
    for ax,(fn,ttl,key) in zip(axes,[("coarse","segment vs annotation_updated","segment_vs_annotation_updated"),
                                     ("immune","immune subtype vs immune_cell_annotation_combined","immune_vs_immune_combined")]):
        tab=pd.read_csv(f"{OUT}/{fn}_confusion_counts.csv",index_col=0)
        rec=tab.div(tab.sum(1),axis=0)  # row=author(true), recall
        im=ax.imshow(rec.values,cmap="Purples",vmin=0,vmax=1,aspect="auto")
        ax.set_xticks(range(rec.shape[1])); ax.set_xticklabels(rec.columns,rotation=45,ha="right",fontsize=9)
        ax.set_yticks(range(rec.shape[0])); ax.set_yticklabels(rec.index,fontsize=9)
        ax.set_xlabel("MY label"); ax.set_ylabel("author label")
        for i in range(rec.shape[0]):
            for j in range(rec.shape[1]):
                v=rec.values[i,j]
                if v>0.05: ax.text(j,i,f"{v*100:.0f}",ha="center",va="center",fontsize=7,color="white" if v>0.5 else "#333")
        ari=vs.loc[key,"ARI"]; agr=vs.loc[key,"agreement"]
        ax.set_title(f"{ttl}\nARI={ari:.2f} · agreement={agr*100:.0f}%",fontsize=11)
    fig.colorbar(im,ax=axes,fraction=0.025,label="row fraction (recall per author class)")
    fig.suptitle("Validation vs authors — two independent typings agree (recall per author class); off-diagonals are the documented divergences",fontsize=14)
    save(fig,"concordance_matrices")

if __name__=="__main__":
    spatial_grid(); umaps(); umap_concordance(); concordance_heatmaps()
    print("== figures done ==")
