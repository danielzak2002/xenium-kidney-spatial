#!/usr/bin/env python
"""Explanatory squidpy panels — make the neighbourhood-enrichment heatmap legible to a bench biologist.

PEDAGOGICAL / DESCRIPTIVE ONLY (no new tests). Builds the z-matrix up step by step on ONE representative
DKD section (1006), and visualises the injury×myeloid SCALE reconciliation (touch-avoidance vs
interstitial co-abundance) in one figure. Reuses the per-section graph recipe + validated labels from
nhood_enrichment_screen.py. SIGN/scale is read, never z-as-importance. Memory-safe (coords+labels only).

Panels:
 A graph→counts→z build-up: (1) Delaunay graph pruned ≤50µm on a region crop ("neighbours = who touches
   whom", mean degree ~5.9); (2) sq.gr.interaction_matrix raw A–B neighbour COUNTS (the numerator);
   (3) the z-matrix as the normalised, compact summary of (2).
 B permutation null worked for TWO pairs — Myeloid×Fibroblast (enriched) and iTAL×Myeloid (avoided):
   label-shuffle expected-count histograms with the observed count marked → z = (obs−mean)/sd.
 C co_occurrence(distance) for iTAL×Myeloid: rises at interstitial radius (06's section co-abundance)
   while the touch-z is negative (panels A/B) → scale reconciliation; B×Plasma as a tight-positive contrast.
 D (optional) Ripley's L for B vs the CSR envelope — "is this type clustered at all" single-type warm-up.
"""
import os, warnings
warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd, anndata as ad, squidpy as sq, scipy.sparse as sp
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
RE=f"{REPO}/analysis/dkd_xenium_reannotation"; DIS=f"{REPO}/analysis/dkd_xenium_disease"; FIG=f"{DIS}/figures"
os.makedirs(FIG,exist_ok=True)
SECTION="1006"; EDGE=50.0; NPERM=1000; SEED=0; rng=np.random.default_rng(SEED)
def savep(fig,n): fig.savefig(f"{FIG}/{n}.png",dpi=165,bbox_inches="tight"); plt.close(fig); print("  [panel]",n)

CATS=["B","Plasma","Myeloid","CD4 T","CD8 T","PT","iPT","TAL","iTAL","DCT","CNT","PC","IC A",
      "Podo","MC","PEC","EC_glom","Fibroblast","VSMC","EC_Peritub","EC_DVR"]
PAL={"B":"#1f77b4","Plasma":"#ff7f0e","Myeloid":"#2ca02c","CD4 T":"#7d3ac1","CD8 T":"#e7298a",
     "PT":"#8a8a8a","iPT":"#c79a5b","TAL":"#8a8a8a","iTAL":"#c79a5b","DCT":"#8a8a8a","CNT":"#8a8a8a",
     "PC":"#8a8a8a","IC A":"#8a8a8a","Podo":"#e7559a","MC":"#e7559a","PEC":"#e7559a","EC_glom":"#e7559a",
     "Fibroblast":"#3fb6a8","VSMC":"#c8a91f","EC_Peritub":"#8c6d5c","EC_DVR":"#8c6d5c"}

# ---------------- build the section graph (same recipe as the screen) ----------------
c=pd.read_parquet(f"{RE}/cells.parquet"); c["sample"]=c.orig_ident.astype(str)
d=c[(c["sample"]==SECTION)&(c.my_label!="Unresolved")].copy()
XY=d[["spatial_x","spatial_y"]].values; lab=d.my_label.values; n=len(d)
a=ad.AnnData(np.zeros((n,1),np.float32),obs=pd.DataFrame({"ct":pd.Categorical(lab,categories=CATS)}))
a.obsm["spatial"]=XY
sq.gr.spatial_neighbors(a,coord_type="generic",delaunay=True)
Dco=a.obsp["spatial_distances"].tocoo(); k=Dco.data<=EDGE
conn=sp.coo_matrix((np.ones(k.sum()),(Dco.row[k],Dco.col[k])),shape=(n,n)).tocsr()
a.obsp["spatial_connectivities"]=conn
degree=conn.getnnz()/n
# unique undirected edge endpoints (i<j)
cc=conn.tocoo(); em=cc.row<cc.col; ER,EC=cc.row[em],cc.col[em]
print(f"section {SECTION}: n={n}, mean degree={degree:.2f}, undirected edges={len(ER)}")

# ---------------- squidpy interaction_matrix (raw counts) + nhood z ----------------
sq.gr.interaction_matrix(a,cluster_key="ct",normalized=False)
IM=np.asarray(a.uns["ct_interactions"],float)
sq.gr.nhood_enrichment(a,cluster_key="ct",n_perms=NPERM,seed=SEED,show_progress_bar=False)
Z=a.uns["ct_nhood_enrichment"]["zscore"].astype(float)

def tickcolor(ax):
    ax.set_xticks(range(len(CATS))); ax.set_yticks(range(len(CATS)))
    ax.set_xticklabels(CATS,rotation=90,fontsize=6); ax.set_yticklabels(CATS,fontsize=6)
    for tl,t in zip(ax.get_xticklabels(),CATS): tl.set_color(PAL[t])
    for tl,t in zip(ax.get_yticklabels(),CATS): tl.set_color(PAL[t])

# ============================ PANEL A: graph -> counts -> z ============================
fig=plt.figure(figsize=(20,6.6)); gs=fig.add_gridspec(1,3,width_ratios=[1.05,1,1])
# (1) graph crop
axg=fig.add_subplot(gs[0,0])
ctr=XY[lab!="PT"].mean(0) if (lab!="PT").any() else XY.mean(0)
# center on a B-lineage-dense spot for an interesting neighbourhood
bl=np.isin(lab,["B","Plasma"]);
if bl.sum()>20:
    from scipy.spatial import cKDTree
    bp=XY[bl]; t=cKDTree(bp); dens=np.array([len(t.query_ball_point(p,r=75)) for p in bp]); ctr=bp[int(np.argmax(dens))]
PAD=170; x0,x1,y0,y1=ctr[0]-PAD,ctr[0]+PAD,ctr[1]-PAD,ctr[1]+PAD
inw=(XY[:,0]>=x0)&(XY[:,0]<=x1)&(XY[:,1]>=y0)&(XY[:,1]<=y1); idxw=set(np.where(inw)[0])
em2=[(i,j) for i,j in zip(ER,EC) if i in idxw and j in idxw]
segs=[[(XY[i,0],XY[i,1]),(XY[j,0],XY[j,1])] for i,j in em2]
axg.add_collection(LineCollection(segs,colors="#cfcfcf",linewidths=0.5,zorder=1))
for t in CATS:
    m=inw&(lab==t); axg.scatter(XY[m,0],XY[m,1],s=26,c=PAL[t],linewidths=0.2,edgecolor="white",zorder=2)
axg.set_aspect("equal"); axg.axis("off")
axg.set_title(f"1 · The graph — Delaunay, pruned ≤{EDGE:.0f}µm\n“neighbours = who touches whom”  (mean degree {degree:.1f})",fontsize=11)
# (2) raw interaction counts
axc=fig.add_subplot(gs[0,1]); M=np.log10(IM+1)
im=axc.imshow(M,cmap="magma"); tickcolor(axc)
axc.set_title("2 · Raw adjacency counts (interaction_matrix)\nlog₁₀(observed A–B neighbour pairs + 1) — the NUMERATOR",fontsize=10.5)
fig.colorbar(im,ax=axc,fraction=0.046,pad=0.04,label="log₁₀ count")
# (3) z matrix
axz=fig.add_subplot(gs[0,2]); imz=axz.imshow(Z,cmap="RdBu_r",vmin=-8,vmax=8); tickcolor(axz)
axz.set_title("3 · Neighbourhood-enrichment z\nraw counts NORMALISED vs label-shuffle (panel B)",fontsize=10.5)
fig.colorbar(imz,ax=axz,fraction=0.046,pad=0.04,label="z (within-section)")
fig.suptitle(f"From graph to z — worked on representative DKD section {SECTION}. SIGN/scale is read, NOT z-as-importance.",fontsize=13)
savep(fig,"nhood_explain_buildup")

# ============================ PANEL B: permutation null for two pairs ============================
def paircount(L,A,B):
    La=L[ER]; Lb=L[EC]
    return int(np.sum((La==A)&(Lb==B))+np.sum((La==B)&(Lb==A)))
PAIRS=[("Myeloid","Fibroblast","#2ca02c","ENRICHED — myeloid sits with stroma"),
       ("iTAL","Myeloid","#c79a5b","AVOIDED — injured tubule not touched by myeloid")]
fig,axes=plt.subplots(1,2,figsize=(14,5))
for ax,(A,B,col,tag) in zip(axes,PAIRS):
    obs=paircount(lab,A,B)
    null=np.array([paircount(rng.permutation(lab),A,B) for _ in range(NPERM)])
    mu,sd=null.mean(),null.std(); z=(obs-mu)/sd
    ax.hist(null,bins=34,color="#bdbdbd",edgecolor="white",label="label-shuffle expected (n=1000)")
    ax.axvline(obs,color=col,lw=2.6,label=f"observed = {obs:,}")
    ax.axvline(mu,color="#444",ls="--",lw=1.2,label=f"null mean = {mu:,.0f}")
    ax.set_title(f"{A} × {B}\n{tag}\nz = (obs − mean)/sd = ({obs:,} − {mu:,.0f})/{sd:,.0f} = {z:+.1f}",fontsize=11)
    ax.set_xlabel("number of A–B neighbour pairs"); ax.set_ylabel("permutations"); ax.legend(fontsize=8.5)
fig.suptitle("The normalisation, worked twice: z = how many SDs the OBSERVED touch-count sits above/below random labelling (section "+SECTION+")",fontsize=12)
fig.tight_layout(rect=[0,0,1,0.84]); savep(fig,"nhood_explain_permnull")

# ============================ PANEL C: co_occurrence scale reconciliation ============================
SUBN=18000
ds=d if n<=SUBN else d.sample(SUBN,random_state=SEED)
ac=ad.AnnData(np.zeros((len(ds),1),np.float32),obs=pd.DataFrame({"ct":pd.Categorical(ds.my_label.values,categories=CATS)}))
ac.obsm["spatial"]=ds[["spatial_x","spatial_y"]].values
iv=np.linspace(0,200,26)
sq.gr.co_occurrence(ac,cluster_key="ct",interval=iv,show_progress_bar=False)
occ=ac.uns["ct_co_occurrence"]["occ"]; ivc=(iv[:-1]+iv[1:])/2
def curve(A,B): return occ[CATS.index(A),CATS.index(B),:]
fig,ax=plt.subplots(1,2,figsize=(15,5.2))
# left: injury x myeloid reconciliation
a0=ax[0]
a0.axhline(1,color="#999",ls=":",lw=1)
a0.plot(ivc,curve("iTAL","Myeloid"),"-o",ms=4,color="#c79a5b",label="P(Myeloid | iTAL) / P(Myeloid)")
a0.axvspan(0,EDGE,color="#e7298a",alpha=0.08); a0.text(EDGE/2,a0.get_ylim()[1]*0.96,"touch scale\n(panels A/B: z<0)",fontsize=8,ha="center",va="top",color="#a11")
a0.set_xlabel("radius from an iTAL cell (µm)"); a0.set_ylabel("co-occurrence ratio")
a0.set_title("Injury × myeloid — the SCALE reconciliation\ntouch-z NEGATIVE, yet co-occurrence RISES at interstitial radius",fontsize=10.5); a0.legend(fontsize=9)
# right: tight-positive contrast B x Plasma
a1=ax[1]; a1.axhline(1,color="#999",ls=":",lw=1)
a1.plot(ivc,curve("B","Plasma"),"-o",ms=4,color="#1f77b4",label="P(Plasma | B) / P(Plasma)")
a1.plot(ivc,curve("iTAL","Myeloid"),"-o",ms=3,color="#c79a5b",alpha=0.5,label="iTAL×Myeloid (for scale)")
a1.axvspan(0,EDGE,color="#e7298a",alpha=0.08)
a1.set_xlabel("radius (µm)"); a1.set_ylabel("co-occurrence ratio")
a1.set_title("Contrast: B × Plasma is a TIGHT positive\nhigh AT the touch scale and decays — a real adjacency niche",fontsize=10.5); a1.legend(fontsize=9)
fig.suptitle(f"Distance-resolved co_occurrence reconciles the scales (section {SECTION}): cell-TOUCH adjacency (squidpy z) ≠ section-level co-abundance (analysis 06, ρ0.82).",fontsize=11.5)
fig.tight_layout(rect=[0,0,1,0.86]); savep(fig,"nhood_explain_cooccurrence")

# ============================ PANEL D (optional): Ripley's L for B vs CSR ============================
try:
    sq.gr.ripley(ac,cluster_key="ct",mode="L",n_simulations=20,seed=SEED)
    rip=ac.uns["ct_ripley_L"]; stat=rip["L_stat"]; sims=rip["sims_stat"]
    RMAX=400.0   # focus on the aggregate scale; beyond ~1.5mm L saturates (finite-window artefact, not dispersion)
    bstat=stat[stat["ct"]=="B"]; bstat=bstat[bstat["bins"]<=RMAX]
    sm=sims[sims["bins"]<=RMAX]
    fig,ax=plt.subplots(figsize=(7.2,5.2))
    lo=sm.groupby("bins")["stats"].min(); hi=sm.groupby("bins")["stats"].max(); md=sm.groupby("bins")["stats"].median()
    ax.fill_between(lo.index,lo.values,hi.values,color="#cdcdcd",alpha=0.8,label="CSR envelope (20 sims, min–max)")
    ax.plot(md.index,md.values,color="#888",lw=1,ls="--",label="CSR median (random)")
    ax.plot(bstat["bins"],bstat["stats"],color=PAL["B"],lw=2.6,label="B cells (observed L)")
    ax.set_xlim(0,RMAX); ax.set_xlabel("radius (µm)"); ax.set_ylabel("Ripley's L"); ax.legend(fontsize=9,loc="upper left")
    ax.set_title(f"Single-type warm-up: B-cell Ripley's L vs CSR (section {SECTION})\nat the aggregate scale B's L sits FAR above the random envelope → B cells ARE clustered",fontsize=10.3)
    savep(fig,"nhood_explain_ripley")
except Exception as e:
    print("  [ripley skipped]",repr(e))

print("== nhood_enrichment_explain done ==")
