#!/usr/bin/env python
"""Integration walkthrough — robust/honest visual narrative for the 68-gene three-cohort
integration. Generates new panels (P1,P2,P4,P5,P6,P7,P8), reuses INTrel_umap/INTrel_dotplot,
assembles a self-contained integration_walkthrough.html. figstyle; region-cropped; read-only."""
import os, warnings, base64, io, html as _html
warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree, ConvexHull
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, FancyArrow
import sys; sys.path.insert(0,"/Users/danie/ClaudeCode/pilot_analyses/xenium/presentation/figures")
import figstyle as fs
from PIL import Image
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"; INT=f"{REPO}/analysis/three_cohort_integration"
PF=f"{REPO}/presentation/figures"; OUT=f"{PF}/integration_walkthrough"; os.makedirs(OUT,exist_ok=True)
COL={"RCC_big":"#1F78B4","RCC_figshare":"#6BAED6","DKD":"#6A3D9A"}
rng=np.random.default_rng(0)
def save(fig,name,tight=True):
    if tight:
        try: fig.tight_layout()
        except Exception: pass
    fig.savefig(f"{OUT}/{name}.png",dpi=300,bbox_inches="tight"); fig.savefig(f"{OUT}/{name}.svg",bbox_inches="tight")
    plt.close(fig); print(f"  [ok] {name}")

# ---------------- P1 — meet the cohorts ----------------
COHORTS=[("RCC_big","kidney_10x (10x Xenium)","clear-cell RCC · 1 section · 380-gene custom panel","455,954"),
         ("RCC_figshare","figshare 25685961 (ProSeg)","clear-cell RCC · 4 patients tumor+adjacent (10 samples)","1,474,861"),
         ("DKD","Demoulin et al. Xenium","diabetic kidney disease · 16 sections · 5k panel","444,675")]
def P1():
    fig,ax=plt.subplots(figsize=(14,3.6)); ax.axis("off"); ax.set_xlim(0,3); ax.set_ylim(0,1)
    for i,(c,src,desc,n) in enumerate(COHORTS):
        ax.add_patch(Rectangle((i+0.04,0.08),0.92,0.84,facecolor="white",edgecolor=COL[c],lw=2.5))
        ax.add_patch(Rectangle((i+0.04,0.78),0.92,0.14,facecolor=COL[c]))
        ax.text(i+0.5,0.85,c,ha="center",va="center",color="white",fontsize=15,fontweight="bold")
        ax.text(i+0.5,0.6,n+" cells",ha="center",fontsize=18,fontweight="bold",color=COL[c])
        ax.text(i+0.5,0.42,src,ha="center",fontsize=10,color="#333")
        ax.text(i+0.5,0.26,desc,ha="center",fontsize=9.5,color="#555",wrap=True)
    fig.suptitle("Three Xenium kidney cohorts — harmonized to ask one question",fontsize=15,y=1.02)
    save(fig,"P1_cohorts",tight=False)

# ---------------- P2 — gene funnel ----------------
def P2():
    d=pd.read_csv(f"{INT}/reliable_drops.csv")
    sz=list(d[d.maxdet_DKD==0].gene); amb=list(d[d.maxdet_DKD>0].gene)
    fig,ax=plt.subplots(figsize=(13,4.6)); ax.axis("off"); ax.set_xlim(0,10); ax.set_ylim(0,5)
    steps=[(123,"123","name-intersection\n(figshare ∩ BIG ∩ DKD)","#8aa9c9"),
           (104,"104","measured on DKD Xenium\n(−19 CosMx-only struct-zeros)","#5a87b8"),
           (68,"68","reliably above ambient\nin ALL three (−36 sub-ambient)","#08519c")]
    xs=[0.6,4.0,7.4]
    for k,(n,big,lab,col) in enumerate(steps):
        w=0.4+2.4*(n/123); ax.add_patch(Rectangle((xs[k],2.5-w/2*0.9),2.4,w*0.9,facecolor=col,alpha=0.9))
        ax.text(xs[k]+1.2,2.5,big,ha="center",va="center",color="white",fontsize=26,fontweight="bold")
        ax.text(xs[k]+1.2,0.7,lab,ha="center",va="center",fontsize=9.5,color="#222")
        if k<2: ax.annotate("",xy=(xs[k+1]-0.05,2.5),xytext=(xs[k]+2.45,2.5),arrowprops=dict(arrowstyle="-|>",lw=2,color="#444"))
    ax.text(5,4.55,"Honest reduction: cross-cohort reliability costs depth (mostly DKD-driven sparsity)",ha="center",fontsize=12,fontweight="bold")
    ax.text(2.2,4.05,"dropped (struct-zero on DKD): "+", ".join(sz[:6])+" …",ha="center",fontsize=8.5,color="#a33")
    ax.text(7.0,4.05,"dropped (sub-ambient): ANGPT2, CXCL9, CXCL10, IL2RA, CCR7 …",ha="center",fontsize=8.5,color="#a33")
    save(fig,"P2_funnel",tight=False)

# ---------------- P4 — lineage availability ----------------
def P4():
    rows=[("B","MS4A1, CD79A","robust"),("Plasma","MZB1, TNFRSF17","robust"),("T","CD3E, CD3G","robust"),
     ("Myeloid","CD68, CD14, CD163, AIF1, ITGAX","robust"),("Endothelial","PECAM1, EGFL7","robust (VWF dropped on DKD)"),
     ("Epithelial","EPCAM, CDH1","robust"),("Stroma","PDGFRA, PDGFRB, ACTA2","robust"),
     ("NK","KLRD1 only","WEAKENED (GNLY, NKG7 dropped on DKD)"),
     ("Treg gate","FOXP3, CTLA4","reduced (IL2RA dropped on DKD)"),("cytotoxic gate","CD8A, GZMK","robust")]
    fig,ax=plt.subplots(figsize=(12,5)); ax.axis("off"); n=len(rows); yh=1.0/(n+1)
    cw=[0.18,0.45,0.37]; x0=[0,0.18,0.63]
    for j,h in enumerate(["lineage / gate","reliable-marker definition","status"]):
        ax.add_patch(Rectangle((x0[j],1-yh),cw[j],yh,facecolor="#333")); ax.text(x0[j]+cw[j]/2,1-yh/2,h,ha="center",va="center",color="white",fontweight="bold",fontsize=11)
    for i,(a,b,c) in enumerate(rows):
        y=1-(i+2)*yh; weak=("WEAK" in c.upper()) or ("reduced" in c)
        for j,val in enumerate([a,b,c]):
            fc="#fdeaea" if (weak and j==2) else ("#eef3f8" if j==0 else "white")
            ax.add_patch(Rectangle((x0[j],y),cw[j],yh,facecolor=fc,edgecolor="#ddd"))
            ax.text(x0[j]+cw[j]/2,y+yh/2,val,ha="center",va="center",fontsize=9,
                    fontweight="bold" if (j==0 or weak and j==2) else "normal",color="#a33" if (weak and j==2) else "#222")
    ax.set_xlim(0,1); ax.set_ylim(1-(n+1)*yh,1.0)
    fig.suptitle("Reliability filter: typing stays marker-faithful, but NK & Treg are weakened (stated, not hidden)",fontsize=13)
    save(fig,"P4_lineage",tight=False)

# ---------------- P5 — B-aggregate object (region-cropped) ----------------
def hull(ax,P,col):
    if len(P)>=4:
        try: h=ConvexHull(P); ax.add_patch(Polygon(P[h.vertices],closed=True,fill=False,edgecolor=col,lw=1.4,alpha=0.85))
        except Exception: pass
def P5():
    df=pd.read_parquet(f"{INT}/cells_labeled_reliable.parquet")
    reps={"RCC_big":"RCC_big","RCC_figshare":"figS10","DKD":"HK2695"}
    fig,axes=plt.subplots(1,3,figsize=(16,5.6))
    for ax,(c,s) in zip(axes,reps.items()):
        g=df[(df.cohort==c)&(df["sample"]==s)]; xy=g[["x","y"]].values; isB=g.is_B.values
        cl=DBSCAN(eps=50,min_samples=20).fit(xy[isB]).labels_
        Bc=xy[isB]; members=Bc[cl!=-1]
        # crop to the largest aggregate's neighbourhood
        if (cl!=-1).any():
            sizes={k:(cl==k).sum() for k in set(cl) if k!=-1}; big=max(sizes,key=sizes.get); cen=Bc[cl==big].mean(0)
        else: cen=xy.mean(0)
        W=350.0; m=(np.abs(xy[:,0]-cen[0])<W)&(np.abs(xy[:,1]-cen[1])<W)
        ax.scatter(xy[m,0],xy[m,1],s=3,c="#ececec",linewidths=0,rasterized=True)
        bm=m&isB; ax.scatter(xy[bm,0],xy[bm,1],s=10,c="#1f77b4",linewidths=0,rasterized=True,label="B cell")
        tm=m&g.is_Treg.values; cm=m&g.is_cyto.values
        ax.scatter(xy[tm,0],xy[tm,1],s=26,c="#d62728",edgecolor="white",linewidth=0.3,rasterized=True,label="Treg")
        ax.scatter(xy[cm,0],xy[cm,1],s=26,c="#2ca02c",marker="^",edgecolor="white",linewidth=0.3,rasterized=True,label="cytotoxic")
        for k in [k for k in set(cl) if k!=-1]:
            P=Bc[cl==k]
            if np.all(np.abs(P.mean(0)-cen)<W): hull(ax,P,"#08519c")
        ax.set_aspect("equal"); ax.axis("off"); ax.set_title(f"{c}  ({s})",color=COL[c],fontsize=12)
        if c=="RCC_big": ax.legend(loc="upper left",fontsize=8,frameon=False,markerscale=1.2)
    fig.suptitle("The B-aggregate object per cohort — DBSCAN B-aggregates (blue, hulls) with Treg/cytotoxic in & around (region-cropped ~700 µm)",fontsize=13)
    save(fig,"P5_spatial",tight=False)

# ---------------- P6 — replication vs separation ----------------
def P6():
    ra=pd.read_csv(f"{INT}/readoutA_reliable.csv").set_index("cohort")
    rd=pd.read_csv(f"{INT}/readoutA_difference.csv")
    nat=pd.read_csv(f"{REPO}/Demoulin26/analysis/bniche_dbscan/differential_treg_vs_cd8.csv")
    fig,axes=plt.subplots(1,2,figsize=(15,5.2))
    # LEFT — integrated replication + difference
    ax=axes[0]; order=["RCC_big","RCC_figshare","RCC_pooled","DKD"]; yp={c:i for i,c in enumerate(order[::-1])}
    cc={"RCC_big":"#1F78B4","RCC_figshare":"#6BAED6","RCC_pooled":"#08519c","DKD":"#6A3D9A"}
    for c in order:
        r=ra.loc[c]; y=yp[c]; ax.errorbar(r.delta,y,xerr=[[r.delta-r.lo],[r.hi-r.delta]],fmt="o",ms=11,color=cc[c],capsize=6,lw=2.2,markeredgecolor="black",zorder=3)
        ax.text(r.delta,y+0.2,f"{r.delta:+.2f}",ha="center",fontsize=10,color=cc[c],fontweight="bold")
    fs.zeroline(ax,0,"v"); ax.set_yticks(list(yp.values())); ax.set_yticklabels(list(yp.keys())); ax.set_ylim(-0.6,len(order)-0.4)
    ax.set_xlabel("Δlog₂ (Treg − cytotoxic) per aggregate"); ax.set_xlim(-1.6,2.2)
    ax.set_title("INTEGRATED — replication (reliable 68-gene set)",fontsize=12,color="#08519c")
    txt="difference (cohort − DKD), 95% CI:\n"+"\n".join(f"  {r.contrast}: {r['diff']:+.2f} [{r.lo:+.2f},{r.hi:+.2f}] {'excl 0' if r.excludes_zero else 'INCL 0'}" for _,r in rd.iterrows())
    ax.text(0.02,0.02,txt,transform=ax.transAxes,fontsize=8.5,va="bottom",family="monospace",bbox=dict(boxstyle="round",fc="#f3f6fa",ec="#bbb"))
    # RIGHT — native separation
    ax=axes[1]; nm={"RCC (tumor)":1,"DKD (kidney)":0}; nc={"RCC (tumor)":"#1F78B4","DKD (kidney)":"#6A3D9A"}
    for _,r in nat.iterrows():
        y=nm[r.cohort]; ax.errorbar(r.delta_log2,y,xerr=[[r.delta_log2-r.ci_lo],[r.ci_hi-r.delta_log2]],fmt="o",ms=14,color=nc[r.cohort],capsize=7,lw=3,markeredgecolor="black",zorder=3)
        ax.text(r.delta_log2,y+0.2,f"{r.delta_log2:+.2f}  (~{r.fold_bias:.0f}×)",ha="center",fontsize=11,color=nc[r.cohort],fontweight="bold")
    fs.zeroline(ax,0,"v"); ax.set_yticks([0,1]); ax.set_yticklabels(["DKD (kidney)","RCC (tumor)"]); ax.set_ylim(-0.5,1.6)
    ax.set_xlabel("Δlog₂ (Treg − cytotoxic), native labels"); ax.set_title("NATIVE-DEPTH — separation (full panels, native subtypes)",fontsize=12,color="#1F78B4")
    ax.text(0.5,0.06,"CIs non-overlapping ->\nmagnitude/separation lives here",transform=ax.transAxes,ha="center",fontsize=10,fontweight="bold",color="#222")
    fig.suptitle("Readout A — two complementary pillars (NOT the same claim): direction REPLICATES (left), separation is at NATIVE depth (right)",fontsize=13)
    save(fig,"P6_pillars")

# ---------------- P7 — Readout B panel limitation ----------------
def P7():
    g=pd.read_csv(f"{INT}/readoutB_usability.csv"); piv=g.pivot(index="gene",columns="cohort",values="endo_detect")
    genes=["ANGPT2","CXCL9","CXCL10","HLA-DRA"]; cohorts=["RCC_big","RCC_figshare","DKD"]
    fig,ax=plt.subplots(figsize=(9,4.4)); x=np.arange(len(genes)); w=0.26
    for i,c in enumerate(cohorts):
        ax.bar(x+(i-1)*w,[piv.loc[gn,c] for gn in genes],w,label=c,color=COL[c])
    ax.axhline(0.02,color="#888",ls="--",lw=1.2); ax.text(len(genes)-0.5,0.022,"2% gate",fontsize=9,color="#555")
    ax.set_xticks(x); ax.set_xticklabels(genes); ax.set_ylabel("detection in endothelial cells")
    ax.set_title("Readout B is a PANEL LIMITATION, not a result:\nANGPT2 / CXCL9 / CXCL10 are sub-ambient on DKD Xenium -> not measurable cross-cohort",fontsize=11)
    ax.legend(frameon=False,fontsize=9)
    save(fig,"P7_readoutB")

# ---------------- P8 — synthesis ----------------
def P8():
    fig,ax=plt.subplots(figsize=(13,5)); ax.axis("off"); ax.set_xlim(0,2); ax.set_ylim(0,1)
    ax.add_patch(Rectangle((0.03,0.1),0.94,0.8,facecolor="#eef6ee",edgecolor="#2ca02c",lw=2))
    ax.add_patch(Rectangle((1.03,0.1),0.94,0.8,facecolor="#eef2f9",edgecolor="#1F78B4",lw=2))
    ax.text(0.5,0.83,"WHAT GENERALIZES",ha="center",fontsize=14,fontweight="bold",color="#2ca02c")
    ax.text(0.5,0.62,"• Conserved B/plasma aggregate scaffold\n   (present in all three cohorts, both platforms)\n\n• Treg-favoring DIRECTION replicates across\n   TWO independent RCC cohorts (different\n   patients, panels, labs): +0.44 to +0.55",
            ha="center",va="center",fontsize=10.5,color="#1a4d1a")
    ax.text(1.5,0.83,"ESTABLISHED ONLY AT NATIVE DEPTH",ha="center",fontsize=13,fontweight="bold",color="#1F78B4")
    ax.text(1.5,0.62,"• MAGNITUDE / SEPARATION of RCC vs DKD\n   (Δ≈+2.6 vs +0.24, non-overlapping CIs)\n\n• On the reliable 68-gene set the integrated\n   difference is positive but its 95% CI\n   INCLUDES zero — underpowered, not absent",
            ha="center",va="center",fontsize=10.5,color="#103a63")
    ax.text(1.0,0.05,"Caveats: RCC_big single-section (wide CI) · sparse DKD aggregates + thin immune detection · reliability filter weakens Treg/NK · native segmentation (ProSeg RCC-only sensitivity pending) · association, not causation",
            ha="center",fontsize=8.5,color="#555")
    fig.suptitle("Synthesis — two pillars: replicated direction (integration) + separation (native depth)",fontsize=14,y=0.99)
    save(fig,"P8_synthesis",tight=False)

# ---------------- assemble HTML ----------------
def datauri(path,max_w=1500):
    im=Image.open(path).convert("RGB")
    if im.width>max_w: im=im.resize((max_w,round(im.height*max_w/im.width)),Image.LANCZOS)
    b=io.BytesIO(); im.save(b,format="JPEG",quality=88,optimize=True)
    return "data:image/jpeg;base64,"+base64.b64encode(b.getvalue()).decode()
def img(path,max_w=1500): return f'<img src="{datauri(path,max_w)}">'

# ================= integration model (pre/post UMAP + Leiden), cached =================
MODEL=f"{OUT}/model.parquet"
def build_model():
    if os.path.exists(MODEL): return pd.read_parquet(MODEL)
    import anndata as ad, scanpy as sc
    print("  building integration model (pre/post UMAP + Leiden) ...")
    A=ad.read_h5ad(f"{REPO}/outputs/objects/three_cohort_123.h5ad")
    lab=pd.read_parquet(f"{INT}/cells_labeled_reliable.parquet")
    for k in ["lineage","is_B","is_Treg","is_cyto"]: A.obs[k]=lab[k].values
    RG=set(open(f"{INT}/reliable_genes.txt").read().split())
    A=A[:,[g for g in map(str,A.var_names) if g in RG]].copy()
    idx=np.sort(np.concatenate([rng.choice(np.where(A.obs.cohort.values==c)[0],min(50000,int((A.obs.cohort.values==c).sum())),replace=False) for c in COL]))
    R=A[idx].copy(); R.obs["cohort"]=R.obs["cohort"].astype(str)
    sc.pp.normalize_total(R,target_sum=1e4); sc.pp.log1p(R); sc.pp.scale(R,max_value=10); sc.tl.pca(R,n_comps=30)
    sc.pp.neighbors(R,use_rep="X_pca"); sc.tl.umap(R); pre=R.obsm["X_umap"].copy()
    import harmonypy; ho=harmonypy.run_harmony(R.obsm["X_pca"],R.obs,["cohort"],max_iter_harmony=20)
    Z=np.asarray(ho.Z_corr); R.obsm["X_pca_harmony"]=Z if Z.shape[0]==R.n_obs else Z.T
    sc.pp.neighbors(R,use_rep="X_pca_harmony"); sc.tl.umap(R); sc.tl.leiden(R,resolution=0.3,flavor="igraph",n_iterations=2); post=R.obsm["X_umap"]
    m=pd.DataFrame(dict(cohort=R.obs.cohort.values,lineage=R.obs.lineage.values,leiden=R.obs.leiden.astype(str).values,
        is_B=R.obs.is_B.values,is_Treg=R.obs.is_Treg.values,is_cyto=R.obs.is_cyto.values,
        pre_x=pre[:,0],pre_y=pre[:,1],post_x=post[:,0],post_y=post[:,1]))
    m.to_parquet(MODEL); return m

LINORD=["B","Plasma","T","Myeloid","NK","Endothelial","Epithelial","Stroma"]
def _knn_purity(emb,coh,k=30):
    from sklearn.neighbors import NearestNeighbors
    nn=NearestNeighbors(n_neighbors=k+1).fit(emb); _,idx=nn.kneighbors(emb)
    codes=pd.factorize(coh)[0]; nb=codes[idx[:,1:]]
    return float((nb==codes[:,None]).mean())   # 1.00 = fully cohort-segregated; ~0.33 = perfectly mixed (3 balanced cohorts)
def P3pre(m):
    coh=m.cohort.values
    pp=_knn_purity(m[["pre_x","pre_y"]].values,coh); qp=_knn_purity(m[["post_x","post_y"]].values,coh)
    fig,axes=plt.subplots(1,2,figsize=(13,5.6))
    for ax,(col,ttl,sc) in zip(axes,[("pre","BEFORE Harmony (68-gene PCA)",pp),("post","AFTER Harmony",qp)]):
        for c in COL:
            mk=coh==c; ax.scatter(m[f"{col}_x"][mk],m[f"{col}_y"][mk],s=2,c=COL[c],linewidths=0,rasterized=True,label=c)
        ax.set_title(f"{ttl}\ncohort kNN purity = {sc:.2f}",fontsize=11); ax.axis("off")
    axes[0].legend(markerscale=4,frameon=False,fontsize=9,loc="upper right")
    fig.text(0.5,0.015,f"cohort kNN purity: 1.00 = fully batch-segregated · 0.33 = perfectly mixed (3 balanced cohorts).  "
             f"Harmony lowers it {pp:.2f} to {qp:.2f} — neighbourhoods become cohort-mixed.",ha="center",fontsize=10,color="#333")
    fig.suptitle("P3-pre · integration was necessary AND did something (pre- vs post-Harmony, same 68-gene features)",fontsize=13)
    fig.subplots_adjust(bottom=0.1); save(fig,"P3pre_umap",tight=False)

def P3assoc(m):
    full=pd.crosstab(m.leiden,m.lineage).reindex(columns=[l for l in LINORD if l in m.lineage.unique()]).fillna(0)
    TOPN=22; ntot=len(full)
    ct=full.loc[full.sum(1).sort_values(ascending=False).index][:TOPN] # top clusters by size
    cov=ct.values.sum()/full.values.sum()
    row=ct.div(ct.sum(1),axis=0)                                       # row-normalize
    fig,ax=plt.subplots(figsize=(9,8.5))
    im=ax.imshow(row.values,cmap="Blues",vmin=0,vmax=1,aspect="auto")
    ax.set_xticks(range(row.shape[1])); ax.set_xticklabels(row.columns,rotation=40,ha="right",fontsize=10)
    ax.set_yticks(range(row.shape[0])); ax.set_yticklabels([f"c{i}" for i in row.index],fontsize=8)
    MK={"B":"MS4A1","Plasma":"MZB1","T":"CD3E","Myeloid":"CD68","NK":"KLRD1","Endothelial":"PECAM1","Epithelial":"EPCAM","Stroma":"PDGFRB"}
    for i in range(row.shape[0]):
        j=int(np.argmax(row.values[i])); dom=row.columns[j]
        ax.text(j,i,f"{row.values[i,j]*100:.0f}",ha="center",va="center",fontsize=7,color="white" if row.values[i,j]>0.5 else "#333")
        ax.text(row.shape[1]-0.35,i,f"-> {dom} ({MK.get(dom,'')})",va="center",fontsize=8.5,color="#222")
    fig.colorbar(im,ax=ax,fraction=0.035,label="fraction of cluster")
    ax.set_title(f"P3-assoc · Leiden cluster -> dominant cell-type ({TOPN} largest of {ntot} clusters, {cov*100:.0f}% of cells)\n"
                 "the 68-gene embedding over-clusters, but clusters still map cleanly to coherent populations",fontsize=10.5)
    save(fig,"P3assoc_heatmap")

def B1_umap(m):
    fig,ax=plt.subplots(figsize=(7.5,6.5))
    ax.scatter(m.post_x,m.post_y,s=2,c="#e8e8e8",linewidths=0,rasterized=True)
    for lab,col,mask in [("B/Plasma","#1f77b4",m.lineage.isin(["B","Plasma"]).values),
                         ("Treg","#d62728",m.is_Treg.values),("cytotoxic","#2ca02c",m.is_cyto.values)]:
        ax.scatter(m.post_x[mask],m.post_y[mask],s=6,c=col,linewidths=0,rasterized=True,label=lab)
    ax.axis("off"); ax.legend(markerscale=3,frameon=False,fontsize=11,loc="upper right")
    ax.set_title("B1 · from the embedding — B/plasma, Treg, cytotoxic in the integrated UMAP\n(the same populations we now follow into tissue)",fontsize=12)
    save(fig,"P5b1_umap")

# ---------------- aggregates helper ----------------
REPS={"RCC_big":"RCC_big","RCC_figshare":"figS10","DKD":"HK2695"}
def get_aggs(df,c,s,topn=3):
    g=df[(df.cohort==c)&(df["sample"]==s)]; xy=g[["x","y"]].values; isB=g.is_B.values
    cl=DBSCAN(eps=50,min_samples=20).fit(xy[isB]).labels_; Bc=xy[isB]
    sizes=[(k,int((cl==k).sum())) for k in set(cl) if k!=-1]; sizes.sort(key=lambda t:-t[1])
    tops=[k for k,_ in sizes[:topn]]; cents=[Bc[cl==k].mean(0) for k in tops]
    return g,xy,isB,cl,Bc,tops,cents

def P5_context(df):
    fig,axes=plt.subplots(2,3,figsize=(16,9.4))
    LC={"B":"#1f77b4","Plasma":"#ff7f0e","T":"#2ca02c","Myeloid":"#9467bd","NK":"#17becf","Endothelial":"#8c564b","Epithelial":"#cccccc","Stroma":"#e377c2"}
    for j,(c,s) in enumerate(REPS.items()):
        g,xy,isB,cl,Bc,tops,cents=get_aggs(df,c,s)
        disp=xy[rng.choice(len(xy),min(70000,len(xy)),replace=False)]
        ax=axes[0,j]; ax.scatter(disp[:,0],disp[:,1],s=1,c="#e9e9e9",linewidths=0,rasterized=True)
        ax.scatter(Bc[cl!=-1,0],Bc[cl!=-1,1],s=4,c="#1f77b4",linewidths=0,rasterized=True)
        for n,cen in enumerate(cents,1):
            ax.add_patch(Rectangle((cen[0]-350,cen[1]-350),700,700,fill=False,edgecolor="#d62728",lw=1.6))
            ax.text(cen[0],cen[1]+430,str(n),color="#d62728",fontsize=12,fontweight="bold",ha="center")
        ax.set_aspect("equal"); ax.axis("off"); ax.set_title(f"{c} ({s}) — whole section, aggregates boxed",color=COL[c],fontsize=11)
        ax=axes[1,j]
        for l,col in LC.items():
            mm=g.lineage.values==l; ax.scatter(xy[mm,0],xy[mm,1],s=1.5,c=col,linewidths=0,rasterized=True)
        ax.set_aspect("equal"); ax.axis("off"); ax.set_title("same section — by lineage",fontsize=10)
    handles=[plt.Line2D([],[],marker="o",ls="",mfc=LC[l],mec="none",ms=8,label=l) for l in LINORD if l in LC]
    axes[1,2].legend(handles=handles,loc="center left",bbox_to_anchor=(1.0,0.5),frameon=False,fontsize=9)
    fig.suptitle("P5 context (overview) · where the zooms come from — whole sections with aggregates boxed, and the same sections by lineage",fontsize=14)
    save(fig,"P5_context",tight=False)

def hull2(ax,P,col):
    if len(P)>=4:
        try: h=ConvexHull(P); ax.add_patch(Polygon(P[h.vertices],closed=True,fill=False,edgecolor=col,lw=1.6,alpha=0.9))
        except Exception: pass
def P5_gallery(df):
    fig,axes=plt.subplots(3,3,figsize=(13,13))
    for r,(c,s) in enumerate(REPS.items()):
        g,xy,isB,cl,Bc,tops,cents=get_aggs(df,c,s)
        for n,(k,cen) in enumerate(zip(tops,cents)):
            ax=axes[r,n]; W=350.0; m=(np.abs(xy[:,0]-cen[0])<W)&(np.abs(xy[:,1]-cen[1])<W)
            ax.scatter(xy[m,0],xy[m,1],s=4,c="#ececec",linewidths=0,rasterized=True)
            ax.scatter(xy[m&isB,0],xy[m&isB,1],s=9,c="#1f77b4",linewidths=0,rasterized=True)
            ax.scatter(xy[m&g.is_Treg.values,0],xy[m&g.is_Treg.values,1],s=30,c="#d62728",edgecolor="white",linewidth=0.3,rasterized=True)
            ax.scatter(xy[m&g.is_cyto.values,0],xy[m&g.is_cyto.values,1],s=30,c="#2ca02c",marker="^",edgecolor="white",linewidth=0.3,rasterized=True)
            hull2(ax,Bc[cl==k],"#08519c")
            ax.set_aspect("equal"); ax.axis("off"); ax.set_title(f"{c} #{n+1}",color=COL[c],fontsize=11)
        for n in range(len(tops),3): axes[r,n].axis("off")
    fig.legend(handles=[plt.Line2D([],[],marker="o",ls="",mfc="#1f77b4",mec="none",ms=9,label="B cell"),
        plt.Line2D([],[],marker="o",ls="",mfc="#d62728",mec="white",ms=9,label="Treg"),
        plt.Line2D([],[],marker="^",ls="",mfc="#2ca02c",mec="white",ms=9,label="cytotoxic")],
        loc="lower center",ncol=3,frameon=False,fontsize=11,bbox_to_anchor=(0.5,0.0))
    fig.suptitle("P5 detail (gallery) · 9 B-aggregates across all three cohorts — Treg/cytotoxic in & around the B-core (region-cropped)",fontsize=14)
    save(fig,"P5_gallery",tight=False)

def S_scoring(df):
    # build the metric on ONE representative aggregate (RCC_figshare #1)
    c,s="RCC_figshare","figS10"; g,xy,isB,cl,Bc,tops,cents=get_aggs(df,c,s,topn=1)
    k=tops[0]; cen=cents[0]; mem=Bc[cl==k]; tree=cKDTree(xy)
    reg=np.unique(np.concatenate([np.asarray(t,int) for t in tree.query_ball_point(mem,r=50)]))
    nreg=len(reg); tin=int(g.is_Treg.values[reg].sum()); cin=int(g.is_cyto.values[reg].sum())
    bgT=float(g.is_Treg.mean()); bgC=float(g.is_cyto.mean())
    fTin=tin/nreg; fCin=cin/nreg; eT=fTin/(bgT+1e-9); eC=fCin/(bgC+1e-9); dl=np.log2(eT+1e-9)-np.log2(eC+1e-9)
    ra=pd.read_csv(f"{INT}/readoutA_reliable.csv").set_index("cohort")
    fig,axes=plt.subplots(1,3,figsize=(16,5.4)); W=300.0
    ax=axes[0]; m=(np.abs(xy[:,0]-cen[0])<W)&(np.abs(xy[:,1]-cen[1])<W)
    ax.scatter(xy[m,0],xy[m,1],s=5,c="#ececec",linewidths=0,rasterized=True)
    ax.scatter(xy[m&isB,0],xy[m&isB,1],s=11,c="#1f77b4",linewidths=0,rasterized=True)
    hull2(ax,mem,"#08519c"); ax.scatter(xy[m&g.is_Treg.values,0],xy[m&g.is_Treg.values,1],s=34,c="#d62728",edgecolor="white",linewidth=0.3)
    ax.scatter(xy[m&g.is_cyto.values,0],xy[m&g.is_cyto.values,1],s=34,c="#2ca02c",marker="^",edgecolor="white",linewidth=0.3)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title("S1–S2 · the aggregate (B hull) +\nTreg/cytotoxic counted within 50 µm",fontsize=11)
    ax=axes[1]; ax.axis("off")
    txt=(f"S2  inside the aggregate ({nreg} cells):\n     Treg (FOXP3/CTLA4+)  = {tin}   ->  {fTin*100:.1f}%\n     cytotoxic (CD8A/GZMK+) = {cin}   ->  {fCin*100:.1f}%\n"
         f"     section background:  Treg {bgT*100:.2f}% · cyto {bgC*100:.2f}%\n\n"
         f"S3  enrichment = frac_inside / frac_background\n     Treg      : {fTin*100:.1f}% / {bgT*100:.2f}% = {eT:.2f}×\n     cytotoxic : {fCin*100:.1f}% / {bgC*100:.2f}% = {eC:.2f}×\n\n"
         f"S4  Δlog₂ = log₂({eT:.2f}) − log₂({eC:.2f}) = {dl:+.2f}\n     (this single aggregate)")
    ax.text(0.0,0.98,txt,va="top",ha="left",fontsize=11,family="monospace",bbox=dict(boxstyle="round",fc="#f3f6fa",ec="#bbb"))
    ax=axes[2]; order=["RCC_big","RCC_figshare","RCC_pooled","DKD"]; yp={cc:i for i,cc in enumerate(order[::-1])}
    cc={"RCC_big":"#1F78B4","RCC_figshare":"#6BAED6","RCC_pooled":"#08519c","DKD":"#6A3D9A"}
    for cohx in order:
        rr=ra.loc[cohx]; y=yp[cohx]; ax.errorbar(rr.delta,y,xerr=[[rr.delta-rr.lo],[rr.hi-rr.delta]],fmt="o",ms=10,color=cc[cohx],capsize=5,lw=2,markeredgecolor="black")
        ax.text(rr.delta,y+0.22,f"{rr.delta:+.2f}",ha="center",fontsize=9,color=cc[cohx],fontweight="bold")
    fs.zeroline(ax,0,"v"); ax.set_yticks(list(yp.values())); ax.set_yticklabels(list(yp.keys()),fontsize=9); ax.set_ylim(-0.6,3.4)
    ax.set_xlabel("Δlog₂ (Treg − cytotoxic)"); ax.set_title("S5 · pool aggregates -> bootstrap ->\nper-cohort Δlog₂ + 95% CI  (-> P6)",fontsize=11)
    fig.suptitle("P5->P6 · object -> number: how one aggregate becomes the per-cohort differential (real values)",fontsize=14)
    save(fig,"P5_scoring")

def Methods():
    rows=[("load / restrict / normalize","anndata 0.12, scanpy 1.11","read_h5ad; pp.normalize_total; pp.log1p; pp.scale"),
     ("dimensionality reduction","scanpy (scikit-learn)","tl.pca (68-gene reliable features → 30 PCs)"),
     ("batch integration (typing)","harmonypy 2.0","run_harmony(X_pca, cohort) — called directly (2.0 returns Z_corr as (N,d); scanpy wrapper assumed (d,N))"),
     ("graph / clustering / embedding","scanpy + leidenalg 0.12 + umap 0.5","pp.neighbors; tl.leiden; tl.umap (pre- and post-Harmony)"),
     ("cell-type labelling","scanpy","tl.score_genes per lineage on reliable markers (uniform, all cohorts)"),
     ("B-aggregate detection","scikit-learn 1.9","cluster.DBSCAN(eps=50 µm, min_samples=20) on B-cell coordinates, per section"),
     ("aggregate footprint / hulls","scipy 1.17","spatial.cKDTree (cells within 50 µm); spatial.ConvexHull (display hulls)"),
     ("enrichment + bootstrap CIs","numpy 2.4","count-pooled Δlog₂; bootstrap by resampling aggregates within cohort (5000×)"),
     ("native-depth separation (ref.)","squidpy 1.8","gr.spatial_neighbors / gr.nhood_enrichment in the native-label bniche_dbscan analysis (P6 right)"),
     ("figures","matplotlib 3.10 + figstyle","slide-grade panels; PNG@300 + SVG; self-contained HTML")]
    fig,ax=plt.subplots(figsize=(15,5.6)); ax.axis("off"); n=len(rows); yh=1.0/(n+1)
    cw=[0.24,0.22,0.54]; x0=[0,0.24,0.46]
    for j,h in enumerate(["step","package(s)","function — what it does"]):
        ax.add_patch(Rectangle((x0[j],1-yh),cw[j],yh,facecolor="#15203a")); ax.text(x0[j]+0.008,1-yh/2,h,va="center",color="white",fontweight="bold",fontsize=10)
    for i,(a,b,d) in enumerate(rows):
        y=1-(i+2)*yh
        for j,val in enumerate([a,b,d]):
            ax.add_patch(Rectangle((x0[j],y),cw[j],yh,facecolor="#eef3f8" if j==0 else "white",edgecolor="#ddd"))
            ax.text(x0[j]+0.008,y+yh/2,val,va="center",fontsize=8.5,fontweight="bold" if j==0 else "normal")
    ax.set_xlim(0,1); ax.set_ylim(1-(n+1)*yh,1.0)
    fig.suptitle("Methods · per-step packages & functions (reproducibility appendix)",fontsize=14)
    save(fig,"Methods",tight=False)

def build_html():
 PANELS=[
 ("P0",None,"Does the RCC immunoregulatory B-aggregate bias generalize across kidney disease?",
   "Three Xenium cohorts (two ccRCC, one DKD), harmonized on a reliably-detected gene set, asked whether the Treg-favoring B-aggregate bias reproduces."),
 ("P1",f"{OUT}/P1_cohorts.png","Meet the cohorts","Three Xenium kidney datasets — all imaged on Xenium, integrated to ask one question."),
 ("P2",f"{OUT}/P2_funnel.png","The gene funnel (honest reduction)","123 name-intersection → 104 measured on DKD Xenium → 68 reliably above ambient in all three. The cost is stated up front; most drops are DKD-driven sparsity."),
 ("P3",f"{PF}/INTrel_umap.png","Integration is sound — embedding","Harmony UMAP on the 68-gene set: the three cohorts MIX (left), lineages RESOLVE (middle), and the two RCC cohorts merge as a batch baseline (right)."),
 ("P3-pre",f"{OUT}/P3pre_umap.png","Was integration necessary? (pre vs post)","UMAP on the SAME 68-gene features before vs after Harmony. Before, cohorts segregate by batch; after, they mix — integration was necessary and did something real."),
 ("P3-assoc",f"{OUT}/P3assoc_heatmap.png","Clustering recovers cell types (without a clean UMAP)","Leiden cluster × cell-type contingency (row-normalized): each cluster maps to a dominant lineage with its defining marker. The 68-gene UMAP is noisy by design, but clustering still recovers coherent populations — shown directly."),
 ("P3-type",f"{PF}/INTrel_dotplot.png","Integration is sound — typing","Per-cohort marker dot-plot: each lineage is defined by the right markers in EVERY cohort — typing is consistent, not blended."),
 ("P4",f"{OUT}/P4_lineage.png","Lineage availability (what reliability cost)","Typing stays marker-faithful, but NK collapses to KLRD1-only and the Treg gate loses IL2RA — stated, not hidden."),
 ("B1",f"{OUT}/P5b1_umap.png","From the embedding — the populations we follow","B/plasma, Treg and cytotoxic cells highlighted in the integrated UMAP — the same populations we now trace from embedding into tissue."),
 ("P5-context",f"{OUT}/P5_context.png","Back to the slides — overview (where the zooms come from)","Whole tissue sections (all cells, light grey) with B-aggregate locations BOXED, and the same sections colored by lineage. Section-scale immune organization and the origin of every zoom."),
 ("P5-detail",f"{OUT}/P5_gallery.png","To the aggregates — detail gallery (9 examples)","Nine B-aggregates spanning all three cohorts, region-cropped: DBSCAN B-core (blue, hull) with Treg (red) and cytotoxic (green) in & around it. The reproducible spatial object behind Readout A."),
 ("P5→P6",f"{OUT}/P5_scoring.png","Object → number (the scoring, transparently)","On one representative aggregate: count Treg/cytotoxic inside vs background (S2) → enrichment (S3) → Δlog₂ for that aggregate (S4) → pool + bootstrap → per-cohort Δlog₂ with CI (S5). Real values, fully traceable into P6."),
 ("P6",f"{OUT}/P6_pillars.png","Readout A — replication vs separation (the core)","Two complementary pillars, NOT the same claim. LEFT (integrated, reliable set): the Treg-favoring direction REPLICATES across two independent RCC cohorts, but the direct cohort−DKD difference CIs INCLUDE zero (underpowered at this depth — CIs shown, not hidden). RIGHT (native depth): the magnitude/separation (Δ≈+2.6 vs +0.24, non-overlapping CIs)."),
 ("P7",f"{OUT}/P7_readoutB.png","Readout B — a panel limitation","ANGPT2/CXCL9/CXCL10 are sub-ambient on the DKD Xenium panel → the vascular/inflammatory axis is not measurable cross-cohort; the within-RCC signal is null. Reported as a limitation, not a result."),
 ("P8",f"{OUT}/P8_synthesis.png","Synthesis","What generalizes (conserved scaffold + replicated Treg direction) vs what is established only at native depth (the magnitude/separation)."),
 ("Methods",f"{OUT}/Methods.png","Methods & packages (appendix)","Per-step packages and the exact functions — scanpy (IO/normalize/PCA/neighbors/Leiden/UMAP), harmonypy (Harmony, called directly), scikit-learn (DBSCAN), scipy (cKDTree/ConvexHull), numpy (bootstrap), squidpy (native-depth nhood), matplotlib + figstyle (figures)."),
 ]
 cards=[]
 for pid,path,title,cap in PANELS:
    body=f'<div class="big">{_html.escape(cap)}</div>' if pid=="P0" else (img(path) if path else "")
    cards.append(f'<section class="card{" title" if pid=="P0" else ""}"><h2><span class="pid">{pid.replace("P3a","P3").replace("P3b","P3")}</span> {_html.escape(title)}</h2>{body}{"" if pid=="P0" else f"<p>{_html.escape(cap)}</p>"}</section>')
 HTML=f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Integration walkthrough — three-cohort kidney</title><style>
body{{margin:0;font-family:-apple-system,'Helvetica Neue',Arial,sans-serif;background:#f4f5f7;color:#1c1c1c;line-height:1.5}}
header{{background:#15203a;color:#fff;padding:20px 28px;position:sticky;top:0;z-index:5;box-shadow:0 2px 8px rgba(0,0,0,.3)}}
header h1{{margin:0;font-size:19px}} header p{{margin:6px 0 0;color:#aebfdc;font-size:13px;max-width:1000px}}
main{{max-width:1080px;margin:0 auto;padding:24px 18px 80px}}
.card{{background:#fff;border-radius:12px;box-shadow:0 1px 4px rgba(0,0,0,.1);margin:18px 0;padding:16px 20px}}
.card.title{{background:#15203a;color:#fff}} .card.title h2{{color:#fff;border:none}} .big{{font-size:16px;color:#dce4f5;margin-top:8px}}
h2{{font-size:16px;margin:0 0 10px;border-left:4px solid #1F78B4;padding-left:10px}}
.pid{{background:#15203a;color:#fff;border-radius:5px;padding:1px 8px;font-size:13px;margin-right:6px}}
.card.title .pid{{background:#fff;color:#15203a}}
img{{width:100%;height:auto;border-radius:8px;border:1px solid #e6e8ec;display:block;margin:6px 0}}
p{{font-size:13.5px;color:#3a4252;margin:10px 2px 2px}}
footer{{max-width:1080px;margin:0 auto;padding:0 18px 60px;color:#8a93a3;font-size:12px}}
</style></head><body>
<header><h1>Does the RCC immunoregulatory B-aggregate bias generalize across kidney disease?</h1>
<p>A robust, honest walkthrough of the 68-gene three-cohort Xenium integration: the integration is sound; the Treg-favoring direction replicates across two independent RCC cohorts; the magnitude/separation is established at native depth. Marginal results are shown WITH their uncertainty.</p></header>
<main>{''.join(cards)}
<div class="card"><h2>Read it as two pillars</h2><p><b>Replication</b> (integration, reliable set): direction reproduces across two independent RCC cohorts — but the direct difference is underpowered (CIs include 0). <b>Separation</b> (native depth): magnitude/non-overlap. These are complementary, not interchangeable; the walkthrough never lets one masquerade as the other.</p></div>
</main><footer>Self-contained · figures from analysis/three_cohort_integration + presentation/figures · pure science (scientific labels only).</footer>
</body></html>"""
 open(f"{OUT}/integration_walkthrough.html","w").write(HTML)
 print(f"wrote integration_walkthrough.html ({len(HTML)//1024} KB)")

if __name__=="__main__":
    m=build_model(); df=pd.read_parquet(f"{INT}/cells_labeled_reliable.parquet")
    P1(); P2(); P4()
    P3pre(m); P3assoc(m); B1_umap(m)
    P5_context(df); P5_gallery(df); S_scoring(df)
    P6(); P7(); P8(); Methods()
    build_html()   # assemble last so it embeds fresh figures
    print("== walkthrough done ==")
