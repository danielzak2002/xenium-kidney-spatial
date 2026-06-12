#!/usr/bin/env python
"""STEP 4 figures + STEP 2/3/5 plots. Per-participant B-lineage gallery (all 16, grouped by
condition, IgAN+MN adjacent to B-rich DKD); within-DKD subgroup plot; DKD-vs-Control descriptive
boxplots with one-offs overlaid; BAFF/APRIL producer-near-B panel. Reuses figstyle. Read-only."""
import os, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull
from sklearn.cluster import DBSCAN
import sys; sys.path.insert(0,"/Users/danie/ClaudeCode/pilot_analyses/xenium/presentation/figures")
import figstyle as fs
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
RE=f"{REPO}/analysis/dkd_xenium_reannotation"; OUT=f"{REPO}/analysis/dkd_xenium_disease"
FIG=f"{OUT}/figures"; os.makedirs(FIG,exist_ok=True); EPS,MINPTS=50,20
def save(fig,name,svg=True):
    fig.savefig(f"{FIG}/{name}.png",dpi=300,bbox_inches="tight")
    if svg: fig.savefig(f"{FIG}/{name}.svg",bbox_inches="tight")
    plt.close(fig); print(f"  [ok] {name}")

cells=pd.read_parquet(f"{RE}/cells.parquet").reset_index(drop=True)
cells["is_B"]=cells.my_label=="B"; cells["is_Plasma"]=cells.my_label=="Plasma"
cells["is_Blin"]=cells.is_B|cells.is_Plasma
sub=pd.read_csv(f"{OUT}/per_sample_substrate.csv"); sub["orig_ident"]=sub.orig_ident.astype(str)
dkd=pd.read_csv(f"{OUT}/dkd_subgroup_split.csv"); dkd["orig_ident"]=dkd.orig_ident.astype(str)
FR=dict(zip(sub.orig_ident,sub.Blin_frac)); COND=dict(zip(sub.orig_ident,sub.Condition))
BRICH=set(dkd[dkd.our_subgroup=="B-rich"].orig_ident)
BCOL="#1f77b4"; PCOL="#ff7f0e"

def aggregates(g):
    B=g[g.is_Blin]
    if len(B)<MINPTS: return B,np.array([])
    lab=DBSCAN(eps=EPS,min_samples=MINPTS).fit(B[["spatial_x","spatial_y"]].values).labels_
    return B,lab

# ============ STEP 4 — per-participant B-lineage gallery (16) ============
# order: B-rich DKD | IgAN, MN (adjacent for contrast) | B-poor DKD | Control | other one-offs
def grouptag(s):
    c=COND[s]
    if c=="DKD": return "B-rich DKD" if s in BRICH else "B-poor DKD"
    if c=="IgA": return "IgAN";
    if c=="MN": return "MN"
    if c=="Control": return "Control"
    return c
ORDER=( list(dkd[dkd.our_subgroup=="B-rich"].orig_ident)             # B-rich DKD
      + [s for s in sub.orig_ident if COND[s]=="IgA"]                # IgAN  (adjacent to B-rich)
      + [s for s in sub.orig_ident if COND[s]=="MN"]                 # MN    (adjacent to B-rich)
      + list(dkd[dkd.our_subgroup=="B-poor"].orig_ident)            # B-poor DKD
      + [s for s in sub.orig_ident if COND[s]=="Control"]           # Control
      + [s for s in sub.orig_ident if COND[s] in ("AA amyloid","C3GN")])  # other one-offs
def gallery():
    fig,axes=plt.subplots(4,4,figsize=(17,17.4))
    for ax,sid in zip(axes.ravel(),ORDER):
        g=cells[cells.orig_ident==sid]; xy=g[["spatial_x","spatial_y"]].values
        xlo,xhi=np.percentile(xy[:,0],[0.3,99.7]); ylo,yhi=np.percentile(xy[:,1],[0.3,99.7])
        ax.scatter(xy[:,0],xy[:,1],s=0.5,c="#e9e9e9",linewidths=0,rasterized=True)
        B,lab=aggregates(g); Bxy=B[["spatial_x","spatial_y"]].values
        ax.scatter(Bxy[B.is_B.values,0],Bxy[B.is_B.values,1],s=5,c=BCOL,linewidths=0,rasterized=True)
        ax.scatter(Bxy[B.is_Plasma.values,0],Bxy[B.is_Plasma.values,1],s=7,c=PCOL,linewidths=0,rasterized=True)
        for k in [k for k in set(lab) if k!=-1]:
            P=Bxy[lab==k]
            if len(P)>=4:
                try: h=ConvexHull(P); ax.add_patch(Polygon(P[h.vertices],closed=True,fill=False,edgecolor="#08306b",lw=1.3,alpha=0.9))
                except Exception: pass
        ax.set_xlim(xlo,xhi); ax.set_ylim(yhi,ylo); ax.set_aspect("equal"); ax.axis("off")
        tag=grouptag(sid); col=fs.DATASET["DKD"] if "DKD" in tag else ("#b8860b" if tag in("IgAN","MN") else ("#555" if tag=="Control" else "#888"))
        ax.set_title(f"{sid} · {tag}\nB-lineage {FR[sid]*100:.2f}% · {int((lab!=-1).sum())} agg cells",
                     fontsize=10.5,color=col,fontweight="bold")
    h=[Line2D([],[],marker="o",ls="",mfc=BCOL,mec="none",ms=9,label="B cell"),
       Line2D([],[],marker="o",ls="",mfc=PCOL,mec="none",ms=9,label="Plasma"),
       Line2D([],[],marker="s",ls="",mfc="none",mec="#08306b",ms=11,label="B-lineage aggregate (DBSCAN eps=50,minPts=20)"),
       Line2D([],[],marker="o",ls="",mfc="#e9e9e9",mec="none",ms=9,label="other cells")]
    fig.legend(handles=h,loc="lower center",ncol=4,frameon=False,fontsize=12,bbox_to_anchor=(0.5,-0.008))
    fig.suptitle("Per-participant B-lineage (B+Plasma) maps — B-rich DKD & IgAN/MN adjacent for contrast; Controls lack aggregates",fontsize=17,y=0.997)
    fig.subplots_adjust(top=0.955,bottom=0.045,wspace=0.04,hspace=0.13)
    save(fig,"b_lineage_gallery_16",svg=False)

# ============ STEP 2 — within-DKD subgroup plot ============
def subgroup():
    d=dkd.sort_values("agg_cells_per10k",ascending=False)
    fig,axes=plt.subplots(1,2,figsize=(14,5.4))
    ax=axes[0]; x=np.arange(len(d)); cols=[fs.DATASET["DKD"] if s=="B-rich" else "#c9b8de" for s in d.our_subgroup]
    ax.bar(x,d.agg_cells_per10k,color=cols,edgecolor="#333",linewidth=0.4)
    ax.axhline(75,ls="--",color="#888",lw=1.2); ax.text(len(d)-0.5,80,"B-rich cut (75)",fontsize=9,color="#555",ha="right")
    for xi,(_,r) in zip(x,d.iterrows()):
        if r.author_Bpredom_ME_20um>0: ax.text(xi,r.agg_cells_per10k+4,f"ME {int(r.author_Bpredom_ME_20um)}",ha="center",fontsize=7.5,color="#6A3D9A")
    ax.set_xticks(x); ax.set_xticklabels(d.orig_ident,rotation=45,ha="right",fontsize=9)
    ax.set_ylabel("B-aggregate cells per 10k"); ax.set_title("within-DKD B-aggregate burden\n(purple = authors' B predom. ME count)",fontsize=11)
    for s in ["top","right"]: ax.spines[s].set_visible(False)
    ax=axes[1]
    for sg,c in [("B-rich",fs.DATASET["DKD"]),("B-poor","#c9b8de")]:
        m=dkd.our_subgroup==sg; ax.scatter(dkd.Blin_frac[m]*100,dkd.agg_cells_per10k[m],s=90,c=c,edgecolor="black",label=sg,zorder=3)
    for _,r in dkd.iterrows(): ax.annotate(r.orig_ident,(r.Blin_frac*100,r.agg_cells_per10k),fontsize=7.5,xytext=(4,3),textcoords="offset points")
    ax.set_xlabel("B-lineage fraction (%)"); ax.set_ylabel("B-aggregate cells per 10k")
    ax.legend(frameon=False,fontsize=10,title="our split"); ax.set_title("B-rich = {1006, HK2695}\n100% concordant with authors' B-predom niche",fontsize=11)
    for s in ["top","right"]: ax.spines[s].set_visible(False)
    fig.suptitle("STEP 2 · within-DKD B-rich/B-poor reproduction (validated vs authors)",fontsize=14)
    save(fig,"dkd_subgroup")

# ============ STEP 3 — DKD vs Control boxplots with one-offs ============
def dkd_control():
    tst=pd.read_csv(f"{OUT}/dkd_vs_control_test.csv").set_index("metric")
    oneoff={"IgA":"#33A02C","MN":"#E31A1C","AA amyloid":"#FF7F00","C3GN":"#6A3D9A"}
    fig,axes=plt.subplots(1,2,figsize=(13,5.6))
    for ax,(metric,lab) in zip(axes,[("Blin_frac","B-lineage fraction"),("agg_cells_per10k","B-aggregate cells per 10k")]):
        groups=[sub[sub.Condition=="DKD"][metric].values,sub[sub.Condition=="Control"][metric].values]
        sc=100 if metric=="Blin_frac" else 1
        bp=ax.boxplot([g*sc for g in groups],positions=[0,1],widths=0.5,patch_artist=True,showfliers=False)
        for patch,c in zip(bp["boxes"],[fs.DATASET["DKD"],"#888"]): patch.set_facecolor(c); patch.set_alpha(0.35)
        for xi,g in zip([0,1],groups): ax.scatter(np.full(len(g),xi)+np.random.uniform(-0.08,0.08,len(g)),g*sc,s=30,c="#333",zorder=3)
        # one-offs overlaid as labeled points at x=2
        for _,r in sub[~sub.Condition.isin(["DKD","Control"])].iterrows():
            ax.scatter(2+np.random.uniform(-0.06,0.06),r[metric]*sc,s=70,c=oneoff.get(r.Condition,"#000"),edgecolor="black",zorder=4)
            ax.annotate(f"{r.orig_ident}",(2,r[metric]*sc),fontsize=7,xytext=(6,0),textcoords="offset points",va="center")
        ax.set_xticks([0,1,2]); ax.set_xticklabels(["DKD\n(n=8)","Control\n(n=3)","one-offs\n(n=1-2)"],fontsize=9)
        p=tst.loc[metric,"p_value"]; ax.set_title(f"{lab}\nMann-Whitney p={p} (descriptive, underpowered)",fontsize=10.5)
        ax.set_ylabel(lab);
        for s in ["top","right"]: ax.spines[s].set_visible(False)
    h=[Line2D([],[],marker="o",ls="",mfc=c,mec="black",ms=9,label=k) for k,c in oneoff.items()]
    axes[1].legend(handles=h,frameon=False,fontsize=8,title="one-offs",loc="upper right")
    fig.suptitle("STEP 3 · DKD vs Control (descriptive) with individual non-DKD references overlaid",fontsize=13)
    save(fig,"dkd_vs_control")

# ============ STEP 5 — BAFF/APRIL producer-near-B panel ============
def baff_panel():
    det=pd.read_csv(f"{OUT}/baff_detection_by_celltype.csv")
    nf=pd.read_csv(f"{OUT}/baff_near_far.csv"); rt=pd.read_csv(f"{OUT}/baff_receptors.csv")
    fig,axes=plt.subplots(1,3,figsize=(17,5.4))
    # (a) detection by cell type
    ax=axes[0]; order=["Myeloid","Fibroblast","VSMC/MC","Endothelial","Epithelial (floor)","B","Plasma","CD4 T","CD8 T"]
    d=det.set_index("cell_type").reindex(order); x=np.arange(len(order)); w=0.4
    ax.bar(x-w/2,d["BAFF_det%"],w,label="BAFF",color="#1f77b4"); ax.bar(x+w/2,d["APRIL_det%"],w,label="APRIL",color="#ff7f0e")
    ax.axhline(d.loc["Epithelial (floor)","BAFF_det%"],ls=":",color="#1f77b4",lw=1)
    ax.set_xticks(x); ax.set_xticklabels(order,rotation=45,ha="right",fontsize=8.5); ax.set_ylabel("% cells detected")
    ax.set_title("(a) ligand detection by cell type\nBAFF myeloid 7.9x epithelial floor; APRIL ambient/broad",fontsize=10)
    ax.legend(frameon=False,fontsize=9)
    for s in ["top","right"]: ax.spines[s].set_visible(False)
    # (b) near/far producer ligand (B-rich DKD)
    ax=axes[1]; g=nf[nf.sample_set=="B-rich DKD"].copy(); g["k"]=g.producer+" "+g.ligand
    x=np.arange(len(g))
    ax.bar(x-0.2,g.near_pct,0.4,label="near agg (<=50um)",color="#08519c")
    ax.bar(x+0.2,g.far_pct,0.4,label="far",color="#bbbbbb")
    for xi,e in zip(x,g.enrich): ax.text(xi,max(g.near_pct.iloc[0],g.far_pct.iloc[0])+0.2 if False else g[["near_pct","far_pct"]].iloc[xi].max()+0.15,f"{e}x",ha="center",fontsize=8,fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(g.k,rotation=30,ha="right",fontsize=8.5); ax.set_ylabel("% producer cells detecting ligand")
    ax.set_title("(b) producer ligand near vs far B-aggregates\nstromal BAFF peri-aggregate (2.1x); myeloid tissue-wide",fontsize=10)
    ax.legend(frameon=False,fontsize=8.5)
    for s in ["top","right"]: ax.spines[s].set_visible(False)
    # (c) receptors on B-lineage near vs far
    ax=axes[2]; x=np.arange(len(rt))
    ax.bar(x-0.2,rt.near_pct,0.4,label="near agg",color="#6A3D9A"); ax.bar(x+0.2,rt.far_pct,0.4,label="far",color="#cbb8e0")
    for xi,e in zip(x,rt.enrich): ax.text(xi,rt[["near_pct","far_pct"]].iloc[xi].max()+0.6,f"{e}x",ha="center",fontsize=8,fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(rt.receptor,fontsize=9); ax.set_ylabel("% B-lineage cells detecting receptor")
    ax.set_title("(c) receptors on B-lineage\nBAFF-R high (36%) & aggregate-elevated (1.6x)",fontsize=10)
    ax.legend(frameon=False,fontsize=8.5)
    for s in ["top","right"]: ax.spines[s].set_visible(False)
    fig.suptitle("STEP 5 · BAFF/APRIL conditioned re-assessment — GO for BAFF/BAFF-R axis, NO-GO for APRIL (ambient)",fontsize=13)
    save(fig,"baff_april_panel")

if __name__=="__main__":
    gallery(); subgroup(); dkd_control(); baff_panel()
    print("== figures done ==")
