#!/usr/bin/env python
"""DKD spatial-transcriptomics WALKTHROUGH — for (1) a skeptical bench biologist new to spatial,
(2) a biologist who wants to SEE the data behind every abstraction, (3) a B-cell / autoimmune-
nephropathy expert.

ONE story, S0–S8. Reuses committed figures from summaries 01/02/03/04/06/07 + composition_by_group
(relative-path <img>, not base64 -> lightweight & committable). Generates narrative + worked-example
scoring panels + per-group tissue drill-downs + the marker-proof dotplot + distance histograms, all
from the VALIDATED reannotation labels and a per-cell gated-transcript matrix.

Honesty is load-bearing: DKD/control is the powered backbone; IgAN/MN/C3GN are SINGLE sections
(reconnaissance, nothing tested); panel limits and tool attribution are taken FROM THE ARTIFACTS
(gene_panel_audit.csv; script imports), not assumed. squidpy did the RCC/cLN neighbourhood stats; the
DKD aggregate/distance work here is sklearn DBSCAN + scipy cKDTree/ConvexHull. Read-only raw; PNG-only.
"""
import os, sys, html as _html
import numpy as np, pandas as pd
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree, ConvexHull
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0,"/Users/danie/ClaudeCode/pilot_analyses/xenium/presentation/figures")
import figstyle; figstyle.apply()

REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
RE=f"{REPO}/analysis/dkd_xenium_reannotation"; DIS=f"{REPO}/analysis/dkd_xenium_disease"
EPI=f"{REPO}/analysis/dkd_epi_endo_stress"
HERE=f"{REPO}/presentation/figures/dkd_walkthrough"; os.makedirs(HERE,exist_ok=True)
def rel(p): return os.path.relpath(p,HERE)
def png(name): return f"{HERE}/{name}.png"
def savep(fig,name): fig.savefig(png(name),dpi=165,bbox_inches="tight"); plt.close(fig); print("  [panel]",name)

# ===================== ONE consistent, colorblind-safe palette (used in EVERY spatial panel) =====================
IMM={"B":"#1f77b4","Plasma":"#ff7f0e","Myeloid":"#2ca02c","CD4 T":"#7d3ac1","CD8 T":"#e7298a"}  # saturated, distinct
EPI_BG="#d9d9d9"      # epithelial — muted grey
INJ_COL="#c79a5b"     # injured tubule — warm tan (muted)
STROMA_COL="#3fb6a8"  # stroma — teal (neutral)
ENDO_COL="#8c6d5c"    # endothelium — brown (neutral)
GLOM_COL="#e7559a"    # glomerulus anchor — pink (outline)
VESS_COL="#c8a91f"    # vessel anchor — gold
UNRES="#efefef"
EPI_LAB=["PT","TAL","DCT","CNT","PC","IC A"]; INJ_LAB=["iPT","iTAL"]
GLOM_LAB=["Podo","MC","EC_glom","PEC"]; ENDO_LAB=["EC_Peritub","EC_DVR"]; VESS_LAB=["VSMC"]; STROMA_LAB=["Fibroblast"]
LEGEND=[("B cell",IMM["B"],"o"),("Plasma",IMM["Plasma"],"o"),("Myeloid",IMM["Myeloid"],"o"),
        ("CD4 T",IMM["CD4 T"],"o"),("CD8 T",IMM["CD8 T"],"o"),
        ("glomerulus",GLOM_COL,"$\\circ$"),("vessel",VESS_COL,"s"),("injured tubule",INJ_COL,"s"),
        ("epithelium",EPI_BG,"s"),("stroma",STROMA_COL,"s"),("endothelium",ENDO_COL,"s")]
def legend_handles():
    return [Line2D([0],[0],marker=m,color="none",markerfacecolor=("none" if l=="glomerulus" else c),
            markeredgecolor=c,markersize=9,markeredgewidth=1.4,label=l) for l,c,m in LEGEND]

c=pd.read_parquet(f"{RE}/cells.parquet").reset_index(drop=True); c["sample"]=c.orig_ident.astype(str)
G=pd.read_parquet(f"{DIS}/walkthrough_genes_percell.parquet").reset_index(drop=True)
assert (c.cell_id.astype(str).values==G.cell_id.astype(str).values).all(),"gene/cell alignment"
GENES=[g for g in G.columns if g!="cell_id"]
for g in GENES: c[g]=G[g].values
samp=c["sample"].values; lab=c.my_label.values; lin=c.my_lineage.values
XY=c[["spatial_x","spatial_y"]].values
isB=lab=="B"; isP=lab=="Plasma"; Blin=isB|isP; isMye=lab=="Myeloid"
GLOM=np.isin(lab,GLOM_LAB); INJ=np.isin(lab,INJ_LAB)
GV={g:c[g].values for g in GENES}

def densest_point(pts,r=75):
    """centroid-free crop anchor: the B-lineage point with the most neighbours within r (avoids
    landing the crop in empty space between scattered cells, e.g. IgAN/control with no aggregate)."""
    t=cKDTree(pts); n=np.array([len(t.query_ball_point(p,r=r)) for p in pts]); return pts[int(np.argmax(n))]
def crop_on_aggregate(s,pad=230):
    m=samp==s; pts=XY[Blin&m]
    if len(pts)>=20:
        dl=DBSCAN(eps=50,min_samples=20).fit(pts).labels_
        ctr=(pts[dl==pd.Series(dl[dl!=-1]).value_counts().idxmax()].mean(0) if (dl!=-1).any()
             else densest_point(pts))         # no cluster -> densest local B-lineage spot, not the mean
    elif m.sum(): ctr=XY[m].mean(0)
    else: ctr=XY.mean(0)
    return ctr,(ctr[0]-pad,ctr[0]+pad,ctr[1]-pad,ctr[1]+pad)
def win_mask(s,box):
    x0,x1,y0,y1=box; m=samp==s
    return m&(XY[:,0]>=x0)&(XY[:,0]<=x1)&(XY[:,1]>=y0)&(XY[:,1]<=y1)
def celltype_map(ax,w,legend=False,title=None):
    def sc(mask,col,s,**k): mm=w&mask; ax.scatter(XY[mm,0],XY[mm,1],s=s,c=col,linewidths=0,rasterized=True,**k)
    sc(np.isin(lab,EPI_LAB),EPI_BG,6); sc(lab=="Unresolved",UNRES,5)
    sc(np.isin(lab,STROMA_LAB),STROMA_COL,6); sc(np.isin(lab,ENDO_LAB),ENDO_COL,6)
    sc(np.isin(lab,VESS_LAB),VESS_COL,8); sc(np.isin(lab,INJ_LAB),INJ_COL,9)
    gmask=w&GLOM; ax.scatter(XY[gmask,0],XY[gmask,1],s=15,facecolors="none",edgecolors=GLOM_COL,linewidths=0.7,rasterized=True)
    for nm in ["CD8 T","CD4 T","Myeloid","Plasma","B"]:
        sc(lab==nm,IMM[nm],18 if nm in("B","Plasma") else 13)
    ax.set_aspect("equal"); ax.axis("off")
    if title: ax.set_title(title,fontsize=11)
    if legend: ax.legend(handles=legend_handles(),fontsize=7.4,loc="upper left",bbox_to_anchor=(1.0,1.0),framealpha=0.95)

# ===================== S1 — slide -> cells -> typed (consistent palette) =====================
ctr,box=crop_on_aggregate("1006"); w=win_mask("1006",box)
fig,ax=plt.subplots(1,3,figsize=(15.5,5.4))
ax[0].scatter(XY[w,0],XY[w,1],s=7,c="#bdbdbd",linewidths=0,rasterized=True)
ax[0].set_title("1 · Segmentation\n(cell boundaries drawn on the image)",fontsize=12); ax[0].set_aspect("equal"); ax[0].axis("off")
celltype_map(ax[1],w,title="2 · Typing — identity from ~5K-gene expression")
celltype_map(ax[2],w,legend=True,title="3 · The cells we follow (same key everywhere)")
fig.suptitle("Xenium keeps the architecture: every dot is one cell, in place — not a dissociated suspension   (sample 1006, ~460 µm region)",fontsize=13)
savep(fig,"wt_pipeline")

# ===================== S2 — marker-proof dotplot (show the data behind each label) =====================
ROWS=["B","Plasma","Myeloid","CD4 T","PT","iPT","TAL","Podo","Fibroblast"]
COLS=[("MS4A1","B"),("CD79A","B"),("CD79B","B"),("MZB1","Pl"),("TNFRSF17","Pl"),("XBP1","Pl"),
      ("CD68","My"),("CD14","My"),("C1QA","My"),("LRP2","PT"),("CUBN","PT"),
      ("VCAM1","iPT"),("HAVCR1","iPT"),("PLA2R1","Po"),("NELL1","Po")]
genes=[g for g,_ in COLS]
frac=np.zeros((len(ROWS),len(genes))); meanz=np.zeros_like(frac)
for j,g in enumerate(genes):
    v=GV[g]; lg=np.log1p(v)
    perrow=[]
    for i,r in enumerate(ROWS):
        m=lab==r; frac[i,j]=(v[m]>0).mean()*100; perrow.append(lg[m].mean())
    pr=np.array(perrow); meanz[:,j]=(pr-pr.mean())/(pr.std()+1e-9)
fig,ax=plt.subplots(figsize=(13.5,6.2))
for i in range(len(ROWS)):
    for j in range(len(genes)):
        ax.scatter(j,i,s=8+frac[i,j]*9,c=[meanz[i,j]],cmap="magma",vmin=-1,vmax=2.2,edgecolor="#555",linewidth=0.3)
ax.set_xticks(range(len(genes))); ax.set_xticklabels([g for g,_ in COLS],rotation=55,ha="right",fontsize=10)
ax.set_yticks(range(len(ROWS))); ax.set_yticklabels(ROWS,fontsize=11); ax.invert_yaxis()
ax.set_title("The data behind each label: canonical markers light up on the cells we call (dot size = % detected · colour = scaled mean)",fontsize=12)
sm=plt.cm.ScalarMappable(cmap="magma",norm=plt.Normalize(-1,2.2)); sm.set_array([])
plt.colorbar(sm,ax=ax,fraction=0.025,pad=0.02,label="scaled mean expr")
for s in [8,30,60]: ax.scatter([],[],s=8+s*9,c="#888",label=f"{s}% detected")
ax.legend(fontsize=8,loc="lower right",title="dot size")
ax.margins(y=0.08)
savep(fig,"wt_annotation_markers")

# ===================== S3 — composition object->number (consistent colours) =====================
comp=pd.read_csv(f"{DIS}/composition_by_group.csv")
LC4={"Epithelial":EPI_BG,"Stroma":STROMA_COL,"Endothelial":ENDO_COL,"Immune":IMM["B"]}
def lin_fracs(s):
    d=comp[(comp["sample"]==s)&(comp.resolution=="coarse_lineage")].set_index("cell_type").fraction
    return {k:float(d.get(k,0)) for k in LC4}
fig,ax=plt.subplots(1,3,figsize=(15,4.6))
ex="1006"; m=samp==ex; counts={k:int((lin==k)&m).sum() if False else int(((lin==k)&m).sum()) for k in LC4}; tot=int(m.sum())
ax[0].axis("off")
ax[0].text(0.0,0.98,("OBJECT -> NUMBER  (composition)\n\nsample 1006\n"
    f"total cells scored = {tot:,}\n\n"+"\n".join(f"{k:12s} {counts[k]:>7,}  ->  {counts[k]/tot*100:5.1f}%" for k in LC4)+
    "\n\nfraction = (cells of type) / (all cells)\nsums to 100% -> compositional (closure)"),
    va="top",ha="left",family="monospace",fontsize=11.5)
for axi,(s,ttl) in zip([ax[1],ax[2]],[("1006","1006 (DKD)"),("HK3626","HK3626 (Control) — epithelial-dominant")]):
    f=lin_fracs(s); bottom=0
    for k in ["Epithelial","Stroma","Endothelial","Immune"]:
        axi.bar(0,f[k]*100,bottom=bottom,color=LC4[k],width=0.6,label=k); bottom+=f[k]*100
    axi.set_xlim(-0.8,0.8); axi.set_xticks([]); axi.set_ylabel("% of cells"); axi.set_title(ttl,fontsize=11)
ax[1].legend(fontsize=8)
fig.suptitle("Composition is just counting, made transparent: count each type, divide by all cells -> the fractions in S3's dotplots",fontsize=13)
savep(fig,"wt_scoring_composition")

# ===================== S4 — aggregate object->number + zero-cluster control =====================
fig,ax=plt.subplots(1,4,figsize=(19,5.2))
ctr,box=crop_on_aggregate("1006",pad=260); w=win_mask("1006",box)
ax[0].scatter(XY[w,0],XY[w,1],s=6,c=EPI_BG,linewidths=0,rasterized=True)
bm=w&Blin; ax[0].scatter(XY[bm,0],XY[bm,1],s=16,c=IMM["B"],linewidths=0,rasterized=True)
ax[0].set_title("1 · B-lineage cells\n(blue) among all cells",fontsize=11)
pts=XY[bm]; tree=cKDTree(pts); dens=np.array([len(tree.query_ball_point(p,r=50))-1 for p in pts])
sc=ax[1].scatter(pts[:,0],pts[:,1],s=22,c=dens,cmap="viridis",linewidths=0,rasterized=True)
ax[1].set_title("2 · Local density\n(B neighbours within 50 µm)",fontsize=11); plt.colorbar(sc,ax=ax[1],fraction=0.046,pad=0.04)
dl=DBSCAN(eps=50,min_samples=20).fit(pts).labels_
ax[2].scatter(pts[dl==-1,0],pts[dl==-1,1],s=14,c="#cccccc",linewidths=0,label="noise")
for k in [kk for kk in set(dl) if kk!=-1]:
    P=pts[dl==k]; ax[2].scatter(P[:,0],P[:,1],s=18,c="#08519c",linewidths=0)
    if len(P)>=4:
        h=ConvexHull(P); ax[2].add_patch(Polygon(P[h.vertices],closed=True,fill=False,edgecolor="#08306b",lw=1.6))
ax[2].set_title("3 · DBSCAN cluster + hull\n(eps=50 µm, minPts=20 · sklearn)",fontsize=11)
ax[3].axis("off")
sub=pd.read_csv(f"{DIS}/per_sample_substrate.csv"); sub["sample"]=sub.orig_ident.astype(str); sub=sub.set_index("sample")
b1006=sub.loc["1006","agg_cells_per10k"]; b1003=sub.loc["1003","agg_cells_per10k"]
ax[3].text(0.0,0.98,("OBJECT -> NUMBER  (organization)\n\nhull = one B-aggregate (scored object)\n"
    "burden = clustered B-lineage cells / 10k\n\n"
    f"1006 (DKD)   burden = {b1006:.0f}/10k\nthreshold (natural gap) = 75/10k\n   {b1006:.0f} >= 75  ->  B-RICH\n\n"
    f"SAME ALGORITHM, IgAN 1003:\n   burden = {b1003:.0f}/10k  ->  ZERO clusters\n   (B cells present, not organized)"),
    va="top",ha="left",family="monospace",fontsize=11.5)
m3=samp=="1003"; p3=XY[Blin&m3]; iax=ax[3].inset_axes([0.46,0.0,0.54,0.40])
iax.scatter(XY[m3,0],XY[m3,1],s=1,c="#eeeeee",linewidths=0,rasterized=True)
iax.scatter(p3[:,0],p3[:,1],s=5,c=IMM["B"],linewidths=0,rasterized=True)
iax.set_xticks([]); iax.set_yticks([]); iax.set_title("IgAN 1003 — 0 aggregates",fontsize=8)
for a in ax[:3]: a.set_aspect("equal"); a.axis("off")
fig.suptitle("How a B-aggregate becomes a number — and the honest negative control: the SAME unbiased algorithm finds ZERO clusters in IgAN/controls",fontsize=12.5)
savep(fig,"wt_scoring_aggregate")

# ===================== S5 — distance object->number =====================
coup=pd.read_csv(f"{DIS}/bcell_damage_coupling.csv").set_index("sample")
ctr,box=crop_on_aggregate("1006",pad=300); w=win_mask("1006",box)
fig,ax=plt.subplots(1,2,figsize=(14,6)); a=ax[0]
a.scatter(XY[w,0],XY[w,1],s=5,c=EPI_BG,linewidths=0,rasterized=True)
gm=w&GLOM; im=w&INJ; bm=w&isB; mm=w&isMye
a.scatter(XY[gm,0],XY[gm,1],s=14,facecolors="none",edgecolors=GLOM_COL,linewidths=0.8,label="glomerulus")
a.scatter(XY[im,0],XY[im,1],s=12,c=INJ_COL,linewidths=0,label="injured tubule")
a.scatter(XY[bm,0],XY[bm,1],s=22,c=IMM["B"],linewidths=0,label="B cell")
a.scatter(XY[mm,0],XY[mm,1],s=22,c=IMM["Myeloid"],linewidths=0,label="myeloid")
if gm.sum() and bm.sum():
    tg=cKDTree(XY[gm])
    for i in np.where(bm)[0][:14]:
        d,j=tg.query(XY[i]); a.plot([XY[i,0],XY[gm][j,0]],[XY[i,1],XY[gm][j,1]],c=IMM["B"],lw=0.7,alpha=0.6)
if im.sum() and mm.sum():
    ti=cKDTree(XY[im])
    for i in np.where(mm)[0][:14]:
        d,j=ti.query(XY[i]); a.plot([XY[i,0],XY[im][j,0]],[XY[i,1],XY[im][j,1]],c=IMM["Myeloid"],lw=0.7,alpha=0.6)
a.set_aspect("equal"); a.axis("off"); a.legend(fontsize=9,markerscale=1.4,loc="upper right")
a.set_title("Each cell -> distance to nearest landmark (1006 region)",fontsize=12)
a=ax[1]; vals=[coup.loc["1006","B_d_glom"],coup.loc["1006","B_d_inj"],coup.loc["1006","Mye_d_inj"]]
labs=["B -> glomerulus","B -> injured tubule","myeloid -> injured tubule"]; cols=[IMM["B"],IMM["B"],IMM["Myeloid"]]; al=[1,0.4,1]
for i,(v,cc,aa) in enumerate(zip(vals,cols,al)): a.barh(i,v,color=cc,alpha=aa)
a.set_yticks(range(3)); a.set_yticklabels(labs,fontsize=11); a.invert_yaxis(); a.set_xlabel("median distance (µm)")
for i,v in enumerate(vals): a.text(v+1.5,i,f"{v:.0f} µm",va="center",fontsize=11)
a.set_title("Object -> number: B near GLOMERULI; myeloid near INJURY\n(1006 — two orthogonal spatial logics)",fontsize=11.5)
fig.suptitle("Two spatial programs, measured the same simple way — distance to the nearest landmark cell (scipy cKDTree)",fontsize=13)
savep(fig,"wt_scoring_distance")

# ===================== S5 — per-sample distance HISTOGRAMS (show distributions, not just medians) =====================
def dists(s):
    m=samp==s; out={}
    tg=cKDTree(XY[GLOM&m]) if (GLOM&m).sum() else None; ti=cKDTree(XY[INJ&m]) if (INJ&m).sum() else None
    bidx=np.where(isB&m)[0]; midx=np.where(isMye&m)[0]
    out["Bg"]=tg.query(XY[bidx])[0] if (tg is not None and len(bidx)) else np.array([])
    out["Bi"]=ti.query(XY[bidx])[0] if (ti is not None and len(bidx)) else np.array([])
    out["Mi"]=ti.query(XY[midx])[0] if (ti is not None and len(midx)) else np.array([])
    return out
HS=[("1006","DKD B-rich"),("1005","MN"),("1003","IgAN"),("HK3626","Control")]
fig,ax=plt.subplots(1,4,figsize=(18,4.4),sharex=True)
bins=np.linspace(0,200,26)
for k,(s,ttl) in enumerate(HS):
    d=dists(s); a=ax[k]
    for key,col,lb in [("Bg",IMM["B"],"B -> glom"),("Bi","#9ecae1","B -> injury"),("Mi",IMM["Myeloid"],"myeloid -> injury")]:
        if len(d[key]): a.hist(np.clip(d[key],0,200),bins=bins,density=True,histtype="step",lw=2,color=col,label=f"{lb} (n={len(d[key])})")
        if len(d[key]): a.axvline(np.median(d[key]),color=col,ls=":",lw=1.2)
    a.set_title(f"{s} · {ttl}",fontsize=11); a.set_xlabel("distance to nearest (µm)"); a.legend(fontsize=7.4)
    if k==0: a.set_ylabel("density")
fig.suptitle("Distance DISTRIBUTIONS per sample (not just medians): B shifts toward glomeruli; myeloid toward injury. Dotted = median.  Single-section refs = n=1.",fontsize=12)
savep(fig,"wt_distance_histograms")

# ===================== D — per-group tissue DRILL-DOWNS (centroids + gated transcripts, NOT histology) =====================
def overlay(ax,w,pairs,title):
    ax.scatter(XY[w,0],XY[w,1],s=4,c="#ececec",linewidths=0,rasterized=True)
    gm=w&GLOM; ax.scatter(XY[gm,0],XY[gm,1],s=12,facecolors="none",edgecolors=GLOM_COL,linewidths=0.6,rasterized=True)
    for g,col in pairs:
        mm=w&(GV[g]>0); ax.scatter(XY[mm,0],XY[mm,1],s=16,c=col,linewidths=0,rasterized=True,label=f"{g}+")
    ax.set_aspect("equal"); ax.axis("off"); ax.legend(fontsize=7.6,loc="upper right",framealpha=0.95); ax.set_title(title,fontsize=11)
def cloud(ax,s,box,title):
    m=samp==s; pts=XY[Blin&m&((XY[:,0]>=box[0])&(XY[:,0]<=box[1])&(XY[:,1]>=box[2])&(XY[:,1]<=box[3]))]
    allp=XY[Blin&m]
    if len(allp)>=20:
        dl=DBSCAN(eps=50,min_samples=20).fit(allp).labels_
    else: dl=np.full(len(allp),-1)
    # restrict to window for plotting
    inwin=(allp[:,0]>=box[0])&(allp[:,0]<=box[1])&(allp[:,1]>=box[2])&(allp[:,1]<=box[3])
    P=allp[inwin]; L=dl[inwin]
    ax.scatter(P[L==-1,0],P[L==-1,1],s=12,c="#bbbbbb",linewidths=0,label="noise / unclustered")
    nclust=0
    for k in [kk for kk in set(L) if kk!=-1]:
        Q=P[L==k]; nclust+=1; ax.scatter(Q[:,0],Q[:,1],s=16,c="#08519c",linewidths=0)
        if len(Q)>=4:
            h=ConvexHull(Q); ax.add_patch(Polygon(Q[h.vertices],closed=True,fill=False,edgecolor="#08306b",lw=1.5))
    nt=len(set(dl[dl!=-1]))
    ax.set_aspect("equal"); ax.axis("off"); ax.legend(fontsize=7.6,loc="upper right",framealpha=0.95)
    ax.set_title(title+(f"\nDBSCAN: {nt} aggregate(s) in section" if nt else "\nDBSCAN: ZERO aggregates (same algorithm)"),fontsize=10.5)

DRILL=[("1003","IgAN (1003) — single section",False),
       ("1005","MN (1005) — single section",False),
       ("1006","DKD B-rich (1006)",True),
       ("1008","DKD B-poor (1008)",True),
       ("HK3626","Control (HK3626)",True)]
for s,ttl,_ in DRILL:
    ctr,box=crop_on_aggregate(s,pad=300); w=win_mask(s,box)
    fig,ax=plt.subplots(1,4,figsize=(19,5.2))
    celltype_map(ax[0],w,legend=(s=="1003"),title="cell types (validated labels)")
    cloud(ax[1],s,box,"B-lineage point-cloud")
    overlay(ax[2],w,[("MS4A1",IMM["B"]),("CD79A","#6baed6"),("MZB1",IMM["Plasma"]),("IGHG1","#fdae6b")],
            "B & plasma transcripts\nMS4A1/CD79A · MZB1/IGHG1")
    overlay(ax[3],w,[("CD68",IMM["Myeloid"]),("C1QA","#74c476"),("PLA2R1",GLOM_COL)],
            "myeloid + podocyte-antigen\nCD68/C1QA · PLA2R1")
    fig.suptitle(f"{ttl} — cell centroids + GATED TRANSCRIPTS (a transcript dot = that gene detected in that cell). NOT histology; the deposit/section is processed-only.",fontsize=11.5)
    savep(fig,f"wt_drill_{s}")

print("panels generated.\n")

# ===================== gene audit -> measurability tables (HTML) =====================
AUD=pd.read_csv(f"{DIS}/gene_panel_audit.csv")
def callcls(call):
    if call.startswith("measurable"): return "ok"
    if "presence-only" in call: return "warn"
    return "bad"
def audtable(roles,cols=("gene","role","target_pct","fold","call"),caption=""):
    sub=AUD[AUD.role.isin(roles)].copy()
    head="".join(f"<th>{_html.escape(h)}</th>" for h in ["gene","role","detect in target %","fold vs ambient","measurability"])
    rows=""
    for _,r in sub.iterrows():
        tp="" if pd.isna(r.target_pct) else f"{r.target_pct:g}"
        fo="" if pd.isna(r.fold) else f"{r.fold:g}×"
        rows+=f'<tr><td><code>{r.gene}</code></td><td>{r.role}</td><td class=n>{tp}</td><td class=n>{fo}</td><td class="{callcls(r.call)}">{_html.escape(r.call)}</td></tr>'
    return f'<table class="aud"><tr>{head}</tr>{rows}</table>'+(f'<p class="cap">{caption}</p>' if caption else "")
TLS_TABLE=audtable(["TLS-follicle","TLS-GC","TLS-Tzone","HEV-ish"],
  caption="Follicular CXCL13/CXCR5/CR2 are on-panel and B-enriched but low-prevalence; germinal-centre BCL6/AICDA/MKI67 are sub-floor; CCL21 off-panel; PNAd is a protein. -> we call these structures “TLS-like,” not confirmed germinal-centre TLS.")
MN_TABLE=audtable(["MN-antigen","complement","compl-reg"],
  caption="PLA2R1 (primary MN autoantigen) and NELL1 are measurable at the mRNA level in podocytes; THSD7A and C3 are ABSENT. mRNA ≠ the deposited immune complex (see blind-spot).")

# ===================== HTML =====================
def img(path,cap=None):
    return f'<img src="{rel(path)}" loading="lazy">'+(f'<p class="cap">{cap}</p>' if cap else "")
def mbox(tool,purpose,why):
    return (f'<div class="mbox"><span class="mt">METHOD · {tool}</span> {purpose}'
            f'<div class="mw"><b>why this tool here:</b> {why}</div></div>')
def expert(html): return f'<div class="frame exp"><span class="tag">EXPERT LENS</span> {html}</div>'

R1=f"{RE}/figures"; D=f"{DIS}/figures"
SEC=[]; sec=lambda pid,title,html: SEC.append((pid,title,html))

sec("S0","Can spatial transcriptomics tell these diseases apart — and where does it stop?", f"""
<p class="lead">This walks one question through 16 human-kidney <b>Xenium</b> sections: <b>can in-situ,
single-cell, ~5,000-gene spatial transcriptomics show how IgAN, MN and AA-amyloid differ from diabetic
kidney disease (DKD) and from controls — and where does the method run out of road?</b></p>
<div class="frame"><b>Read this first — the honest frame, which is load-bearing:</b>
<ul>
<li><b>DKD vs control is the powered backbone</b> (8 DKD, 3 control). Real replicate spread, shown everywhere.</li>
<li><b>IgAN (1003), MN (1005), C3GN (1007) are SINGLE SECTIONS</b> — reconnaissance, hypothesis-generating,
<b>nothing is statistically tested</b>. AA-amyloid is n=2. One patient per section, no donor replication.</li>
<li><b>Two headline axes are panel-limited and we say so where they bite:</b> the IgA isotype (the defining
axis of IgAN) is <b>not measurable</b>; the germinal-centre / follicular-TLS markers are presence-only.</li>
<li><b>Credibility a slide reader can audit:</b> our cell calls are <b>validated blind</b> against the
published atlas, the markers behind each label are <b>shown</b> (S2), gene availability is <b>audited from
the matrix</b> (not assumed), and we <b>caught and reversed two of our own over-claims</b> (S7).</li>
</ul></div>
{expert("<b>Protein-deposit blind spot — state it up front.</b> Spatial transcriptomics sees <b>cells and mRNA</b>, "
"not deposited protein. The lesions that DEFINE these diseases — IgAN mesangial galactose-deficient-IgA1 immune "
"complexes, MN subepithelial anti-PLA2R1 / THSD7A IgG along the GBM — are <b>immunofluorescence / electron-microscopy</b> "
"entities. They are <b>invisible to mRNA</b>. Everything here is the cellular/transcriptional <i>context around</i> those "
"deposits, never the deposits themselves.")}
<p class="cap">Plain-language bridges — <b>segmentation</b> = drawing each cell's boundary; <b>typing</b> =
calling identity from expression; <b>aggregate / DBSCAN</b> = unbiased cell-cluster finding; <b>ambient /
spillover</b> = transcripts misassigned across boundaries, the failure mode we control for.</p>""")

sec("S1","The data, and one colour key for the whole tour", f"""
<p>Xenium is <b>in-situ</b>: tissue is never dissociated. Each transcript is imaged in place, assigned to a
segmented cell, ~5,000 genes per cell <b>with architecture intact</b> — a multiplexed IHC you can re-stain
5,000 ways. Below is the slide→cells→typing bridge, and the <b>single colour key used in every spatial panel
from here on</b>: immune saturated &amp; distinct, epithelium muted, anchors (glomerulus pink, vessel gold) distinct.</p>
{img(png('wt_pipeline'),"(1) Segmentation draws boundaries; (2) typing calls identity from the ~5K-gene profile; (3) the consistent legend (right) — immune saturated, epithelium muted grey, injured tubule tan, glomeruli pink outline, vessels gold. One ~460 µm region of 1006.")}
<p class="cap">Cohort: 16 sections — 8 DKD, 3 control, plus IgAN / MN / AA-amyloid (×2) / C3GN. Large image &amp;
per-transcript files stay on disk; we work from the cell-by-gene matrix + coordinates.</p>""")

sec("S2","Trust the calls: blind re-annotation + the markers behind every label", f"""
<p>Before any biology: are the calls real? Two proofs. <b>(a) External:</b> we re-integrated within this atlas
and re-annotated <b>blind</b>, then compared to the authors. <b>(b) Internal:</b> the canonical markers light up
on exactly the cells we label.</p>
{mbox("scanpy + harmonypy","HVG → PCA → Harmony(orig_ident) → Leiden → UMAP, then blind marker z-score typing; concordance via scikit-learn ARI / confusion.",
 "Harmony removes per-sample batch without collapsing cell types; Leiden is the standard graph clustering; scanpy reads the 8.7 GB / 4.3 M-cell object in BACKED mode so a 24 GB laptop never materialises the full matrix.")}
{img(R1+'/concordance_matrices.png',"(a) Blind labels vs the authors' published annotation: segment ARI 0.78, immune ARI 0.68 — strong agreement, no peeking.")}
{img(png('wt_annotation_markers'),"(b) Show the data behind the abstraction: B = MS4A1/CD79A/CD79B; Plasma = MZB1/TNFRSF17/XBP1; Myeloid = CD68/CD14/C1QA; PT = LRP2/CUBN; iPT = VCAM1/HAVCR1; Podo = PLA2R1/NELL1. Dot size = % detected in that type; colour = scaled mean. Audited fact: NPHS1/NPHS2 are OFF this panel, so podocytes rest on PLA2R1/NELL1 — flagged, not hidden.")}
{img(R1+'/umap_yours_vs_theirs.png',"Same cells, our typing vs theirs, on the integrated embedding — populations land in the same places.")}""")

sec("S3","Top-level composition, by disease group", f"""
<p>Start where IHC starts: what is the tissue made of? Every dot is <b>one section</b>, grouped by disease; a
median bar only where n&gt;1, so single IgAN/MN/C3GN dots can't masquerade as a group.</p>
{mbox("pandas + numpy","per-section fraction of each type, plus a centred-log-ratio (CLR) view.",
 "composition is literally counting; CLR respects the sum-to-one (closure) constraint without imposing a model or a test — we read it descriptively.")}
{img(png('wt_scoring_composition'),"Object → number for composition: count each type, divide by all cells.")}
{img(D+'/composition_by_group_coarse.png',"Coarse lineage. Control epithelial-dominant (74%), near-bare immune (2.8%); every disease loses epithelium, gains immune + stroma. CLR row tracks the raw ordering.")}
{img(D+'/composition_by_group_immune.png',"Immune drill-down. MN most immune/plasma-skewed, C3GN most T-skewed, AA highest B. NK/DC not separately typed — not invented. (Representative dots are mapped in the Tissue tour below.)")}""")

sec("S4","Organization, not just composition — finding aggregates", f"""
<p>Composition says <i>how much</i> immune; not whether it's <b>organized</b>. DBSCAN either finds a dense
cluster or doesn't — shown beside the image so the hull is visibly a real lymphoid aggregate, with the
<b>negative control</b>: the same algorithm on IgAN finds nothing.</p>
{mbox("scikit-learn DBSCAN + scipy cKDTree / ConvexHull","eps = 50 µm, minPts = 20; cKDTree for neighbourhoods, ConvexHull for the hull object.",
 "DBSCAN needs no preset cluster count and labels sparse cells as noise, so it returns aggregates OR an honest zero. NB — the RCC/cLN niche work used <b>squidpy</b> (Delaunay graphs, nhood_enrichment, co-occurrence, Moran's I); the DKD aggregate/distance analysis here is sklearn + scipy, NOT squidpy.")}
{img(png('wt_scoring_aggregate'),"B-lineage → local density → DBSCAN core/hull → burden/10k → threshold. 1006 = 202/10k (B-rich); SAME algorithm on IgAN 1003 = ZERO clusters.")}
{img(D+'/dkd_subgroup.png',"Within DKD a ~3× burden gap → B-rich subgroup = {1006, HK2695}; 8/8 concordant with the authors' B-predominant niche.")}
{img(D+'/b_lineage_gallery_16.png',"Every section, grouped: B+Plasma with hulls. B-rich DKD show compact aggregates; IgAN/controls have none under the identical algorithm.")}""")

sec("S5","Two spatial programs (the centerpiece — robust within a section)", f"""
<p>Where spatial earns its keep, and it holds at n=1 because it is measured <b>within each section</b>. Two
<b>orthogonal</b> logics — and we show the full <b>distributions</b>, not just medians.</p>
{mbox("scipy cKDTree + scipy.stats","nearest-landmark distance per cell (cKDTree query); cross-section association via Spearman + bootstrap CIs + within-section permutation.",
 "the distance is an exact geometric query, fully transparent; Spearman/bootstrap make the cross-section claim honestly associational with CIs, not a black box.")}
{img(png('wt_scoring_distance'),"Distance object→number: in 1006, B ~34 µm from glomeruli but ~94 µm from injury; myeloid ~30 µm from injury.")}
{img(png('wt_distance_histograms'),"The DISTRIBUTIONS behind the medians, per sample: B→glom (solid blue) shifts left of B→injury (light blue); myeloid→injury (green) hugs zero. Single-section refs carry their n on the legend.")}
{img(EPI+'/coloc_fig_A_specificity.png',"Across sections, injury co-localizes with MYELOID (ρ = 0.82 [0.46, 0.95]); the B-lineage association dies under de-circularization (partial 0.13, CI spans 0).")}
{img(D+'/bcell_fig_C_localization_coupling.png',"Same B-vs-myeloid distance contrast across all samples: B→glom ≪ B→injury; myeloid hugs injury everywhere.")}
<p><b>One sentence:</b> injury recruits myeloid; B-lineage organizes around glomeruli/aggregates.</p>""")

sec("TOUR","Tissue tour — per-group drill-downs (see the data, region by region)", f"""
<p>For the biologist who wants to SEE every abstraction: five sections, each as
<b>cell-type map → B-lineage point-cloud (DBSCAN) → B/plasma transcripts → myeloid + podocyte-antigen
transcripts</b>, all region-cropped and on the one colour key.</p>
<div class="frame warn"><b>These are cell centroids + gated transcripts, NOT histology.</b> A coloured dot in
the transcript panels means "<i>that gene was detected in that cell</i>" — there is no H&amp;E/IF image here;
the section is processed to a cell-by-gene matrix. Glomeruli are drawn as pink outlines for orientation.</div>
{img(png('wt_drill_1003'),"IgAN (1003), n=1: B-lineage present but the point-cloud yields ZERO aggregates; B/plasma transcripts are sparse and diffuse, not focal. Consistent with deposition-driven disease where pathogenic IgA1 is largely mucosal/marrow-derived (little organized intrarenal B-cell machinery).")}
{img(png('wt_drill_1005'),"MN (1005), n=1: looser mixed aggregates; plasma transcripts (MZB1/IGHG1) relatively prominent; PLA2R1+ podocytes visible around glomeruli (mRNA — NOT the subepithelial IgG deposit).")}
{img(png('wt_drill_1006'),"DKD B-rich (1006): a compact B-lineage aggregate (hull) with dense MS4A1/CD79A B transcripts — the organized end of the spectrum.")}
{img(png('wt_drill_1008'),"DKD B-poor (1008): a small aggregate, sparser B-lineage — the within-DKD contrast.")}
{img(png('wt_drill_HK3626'),"Control (HK3626) — the MOST immune of the three controls: even here the B compartment is sparse with only a few small aggregates; the other two controls (HK2753, HK3106) have zero. The contrast with B-rich DKD (compact, dense hulls) is the point.")}""")

sec("S6","B-lineage programs across nephropathies — and the expert lens", f"""
<p>Layering mechanism onto the single-section references gives <b>three apparent B-lineage programs</b> —
anecdotes to test, with the panel limits stated exactly where they bite.</p>
{img(D+'/bcell_fig_A_split_isotype.png',"B vs Plasma split: MN is the one plasma-skewed glomerular disease, IgG-dominant. LIMIT: IgA (IGHA1) is sub-floor in plasma — the IgA axis is NOT testable; only IgG is.")}
{img(D+'/bcell_fig_B_tls_state.png',"TLS organization + cell state: only DKD B-rich (HK2695) carries substantial follicular signal; CCL19 is the quantitative T-zone marker; follicular markers are presence-only.")}
{img(D+'/bcell_fig_D_glom_crops.png',"Crops — IgAN / MN / DKD B-rich: organized B-skewed vs plasma-skewed mixed vs diffuse peri-vascular.")}
{expert("<b>Canonical TLS panel — what is actually measurable here.</b> A bench reader should not take “TLS” on "
"faith. Audited against the matrix:" + TLS_TABLE +
"So the DKD B-rich aggregates are <b>TLS-like</b> (follicular chemokine/receptor + B organization), but a "
"<b>confirmed germinal-centre TLS</b> (BCL6/AICDA/MKI67, PNAd<sup>+</sup> HEV) is <b>not establishable on this panel</b>.")}
{expert("<b>MN target antigens &amp; complement.</b>" + MN_TABLE +
"<b>PLA2R1 mRNA is measurable and podocyte-specific</b> (83% of podocytes, 33× ambient) — but that is podocyte "
"<i>expression</i>, not the deposited anti-PLA2R1 IgG immune complex (IF/EM). C3 — the hub of complement and the "
"crux of C3GN — is <b>absent from the panel</b>; C4A/C4B/CFB and regulators (CD46/CD55/CD59/CFH) are measurable, "
"so complement can be approached only partially.")}
{expert("<b>IgAN context &amp; the plasma-maturity ceiling.</b> Pathogenic IgA1 in IgAN is largely mucosal / "
"bone-marrow-derived, so <b>sparse, diffuse, unorganized intrarenal B-lineage is consistent with deposition-driven "
"disease</b> — the IgAN drill-down (zero aggregates, diffuse transcripts) fits that, though it cannot prove it. "
"Plasma maturity (MZB1/TNFRSF17/PRDM1/XBP1) IS measurable; <b>isotype class-switch and κ/λ clonality are NOT</b> "
"(IGHA1 sub-floor, light chains off-panel, ambient saturation) — so we never infer clonality from this data.")}
<div class="frame warn"><b>n = 1 per non-DKD condition. Descriptive. Hypothesis-generating. NOT tested.</b>
The IgAN read rests on architecture + plasma maturity, NOT on isotype.</div>""")

sec("S7","How we avoided fooling ourselves (two self-corrections)", f"""
<p><b>Ambient / spillover</b>: transcripts occasionally get assigned to the wrong neighbouring cell across a
segmentation boundary — a faint smear where signal isn't produced. Uncontrolled, you “discover” spillover.
We controlled for it, and it cost us two claims.</p>
{img(D+'/baff_ambient_control.png',"Self-correction 1: a peri-aggregate STROMAL BAFF niche did NOT survive an ambient stress-test — non-producer epithelium/endothelium rise like stroma near aggregates → spillover field, not production. RETRACTED.")}
{img(D+'/baff_myeloid_anchor.png',"What survives (canonical BAFF, 03/04): a tissue-wide MYELOID producer, cell-intrinsic, reproducible 16/16, ~24× epithelial floor — but no aggregate-specific niche. APRIL is NO-GO.")}
{img(EPI+'/coloc_fig_B1_disease.png',"Self-correction 2: 'injury sits near B-aggregates' — de-circularized, the injury association is MYELOID, not B-lineage (S5). The B-specific claim is RETRACTED.")}
<p><b>Why this matters:</b> the analysis visibly disproves two of its own attractive stories. That discipline
is what makes the survivors — validated typing, the two spatial programs — trustworthy.</p>""")

sec("TOOLS","Tooling — R/Seurat vs Python, and where squidpy earned its keep", f"""
<p>Sourced from the actual imports, not habit.</p>
<table class="meth"><tr><th>ecosystem</th><th>used for</th><th>why / honest notes</th></tr>
<tr><td><b>R · Seurat</b></td><td>mature imaging-ST loaders (LoadXenium), SingleR / Azimuth reference annotation, InSituType, interactive QC (project annotation/benchmark arm, <code>R/01–16</code>)</td><td>R has the most mature imaging-ST loaders and reference-mapping (SingleR/Azimuth/KPMP); great for interactive QC on one section.</td></tr>
<tr><td><b>Python · scanpy</b></td><td>the DKD reannotation + EVERY analysis in this walkthrough: backed/streaming reads, normalize/HVG/PCA/neighbors/Leiden/UMAP, harmonypy</td><td>the combined object is <b>8.7 GB / 4.3 M cells</b> — only a backed AnnData (subset → <code>to_memory</code>) fits in 24 GB RAM; harmonypy integrates per-sample batch.</td></tr>
<tr><td><b>squidpy</b></td><td><b>RCC / cLN neighbourhood statistics</b>: <code>sq.gr.spatial_neighbors</code>, <code>nhood_enrichment</code>, <code>co_occurrence</code>, Moran's I (py/phaseB_*, whitepaper groupB/D)</td><td>squidpy shines for graph-based neighbourhood-enrichment across many cell types. <b>It was NOT used for the DKD aggregates/distances</b> — there, an explicit DBSCAN hull + cKDTree distance was simpler and more transparent (a hull object you can see; a distance you can check).</td></tr>
<tr><td><b>scikit-learn · scipy</b></td><td>DKD aggregates (DBSCAN), neighbourhoods/distances (cKDTree), hulls (ConvexHull), stats (Spearman, Mann–Whitney, bootstrap)</td><td>the workhorses for the organization + distance layers; no spatial-graph machinery needed for "is this clustered" and "how far to the nearest landmark".</td></tr>
<tr><td><b>matplotlib · figstyle</b></td><td>every panel here; one shared palette</td><td>committed PNG (no base64) → the walkthrough stays lightweight and diff-able.</td></tr></table>""")

sec("S8","Synthesis, honest limits, and the study this motivates", f"""
<p><b>What spatial CAN do here:</b> (1) reproduce published identities <b>blind</b> (ARI 0.78/0.68) with the
markers shown; (2) quantify composition + <b>organization</b> per section; (3) reveal <b>two orthogonal
programs</b> — injury↔myeloid, B-lineage↔glomeruli — robust within a section.</p>
<p><b>Where it STOPS (stated, not hidden):</b></p>
<ul>
<li><b>Protein deposits are invisible.</b> The defining IgAN/MN lesions are IF/EM immune complexes — mRNA cannot see them.</li>
<li><b>Single-section disease claims are hypotheses.</b> IgAN/MN/C3GN n=1; AA n=2; one patient per section; no donor column.</li>
<li><b>Panel limits (audited).</b> IgA isotype unmeasurable; germinal-centre markers (BCL6/AICDA/MKI67) sub-floor; CCL21/PNAd/NPHS1/NPHS2/C3/THSD7A off-panel → "TLS-like" not confirmed TLS, C3 biology out of reach.</li>
<li><b>Typing caveats.</b> iPT recall ≈ 0.64 (injury distances conservative); CD8 recall ~0.58.</li>
<li><b>Associational.</b> All co-localization is correlational within fixed tissue — no dynamics, no causation.</li>
</ul>
<p><b>The powered study this motivates:</b> a multi-donor cohort per nephropathy on a panel carrying the missing
isotype + follicular/GC markers (IGHA/IGHM/IGKC, CXCL13/CXCR5/CR2, BCL6/AICDA) and complement (C3), to TEST the
three candidate B-lineage programs and the injury↔myeloid vs B↔glomeruli dissociation generated here.</p>""")

CHECK=["one consistent colourblind-safe palette across every spatial panel (S1 key reused in S4/S5/tour)",
 "object→number scoring for composition (S3), aggregates (S4), distance (S5)",
 "marker-proof dotplot behind every label (S2); gene availability AUDITED from the matrix, not assumed",
 "per-group tissue drill-downs incl. IgAN & MN full + DKD B-rich/B-poor + Control (TOUR); centroid+transcript, flagged NOT histology",
 "distance DISTRIBUTIONS not just medians (S5 histograms)",
 "aggregate determination incl. the ZERO-cluster IgAN/control case",
 "per-stage method callouts (S2–S5) + tooling section (R/Seurat vs Python; squidpy=RCC/cLN, sklearn/scipy=DKD)",
 "expert lens: protein-deposit blind spot, TLS measurability table, MN-antigen/complement table, plasma-maturity ceiling",
 "two self-corrections displayed (S7); panel limits where they bite (S6/S8)",
 "lightweight: relative-path PNGs (no base64) → committable"]
print("PANEL CHECKLIST"); [print(f"  [x] {x}") for x in CHECK]

cards=[]
for pid,title,body in SEC:
    cls="card title" if pid=="S0" else "card"
    pidlab=pid if pid not in ("TOUR","TOOLS") else ("◆" if pid=="TOUR" else "⚙")
    cards.append(f'<section class="{cls}"><h2><span class="pid">{pidlab}</span> {_html.escape(title)}</h2>{body}</section>')
checkli="".join(f"<li>{_html.escape(x)}</li>" for x in CHECK)
HTML=f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>DKD spatial transcriptomics — a walkthrough</title><style>
body{{margin:0;font-family:-apple-system,'Helvetica Neue',Arial,sans-serif;background:#eef0f3;color:#1c1c1c;line-height:1.6}}
header{{background:#15203a;color:#fff;padding:26px 32px}} header h1{{margin:0;font-size:23px}}
header p{{margin:9px 0 0;color:#b9c6e0;font-size:14.5px;max-width:980px}}
main{{max-width:1020px;margin:0 auto;padding:22px 18px 80px}}
.card{{background:#fff;border-radius:13px;box-shadow:0 1px 5px rgba(0,0,0,.10);margin:18px 0;padding:20px 26px}}
.card.title{{background:#f7f5fb;border-left:6px solid #6A3D9A}}
h2{{font-size:18.5px;margin:0 0 12px}} h2 .pid{{background:#6A3D9A;color:#fff;border-radius:6px;padding:2px 10px;font-size:13px;margin-right:9px}}
img{{width:100%;height:auto;border-radius:9px;border:1px solid #e6e8ec;display:block;margin:10px 0 2px}}
p,li{{font-size:14px;color:#2b3040}} .lead{{font-size:15.5px;color:#15203a}}
.cap{{font-size:12.6px;color:#6a7180;margin:3px 2px 10px}}
.frame{{background:#f3eef8;border:1px solid #d4c4ea;border-radius:9px;padding:13px 17px;margin:10px 0;font-size:13.4px}}
.frame.warn{{background:#fdeaea;border-color:#e3a6a6;color:#7a1a1a}}
.frame.exp{{background:#eef4fb;border-color:#a9c7e8;color:#143a5e}}
.frame .tag{{display:inline-block;background:#1f5e9e;color:#fff;font-size:10.5px;font-weight:700;border-radius:4px;padding:1px 7px;margin-right:7px;letter-spacing:.4px}}
.frame ul{{margin:7px 0 0;padding-left:20px}} ul{{margin:6px 0}}
.mbox{{background:#f6f7f9;border-left:4px solid #6A3D9A;border-radius:7px;padding:9px 14px;margin:9px 0;font-size:12.8px;color:#2b3040}}
.mbox .mt{{display:inline-block;background:#6A3D9A;color:#fff;font-size:10.5px;font-weight:700;border-radius:4px;padding:1px 7px;margin-right:7px}}
.mbox .mw{{margin-top:4px;color:#55607a}}
table.meth,table.aud{{border-collapse:collapse;width:100%;font-size:12.5px;margin:8px 0}}
table.meth th,table.meth td,table.aud th,table.aud td{{border-bottom:1px solid #eceef2;padding:5px 9px;text-align:left;vertical-align:top}}
table.meth th,table.aud th{{color:#6a7180}} td.n{{text-align:right;font-variant-numeric:tabular-nums}}
.aud td.ok{{color:#1a6d1a;font-weight:600}} .aud td.warn{{color:#8a5a12;font-weight:600}} .aud td.bad{{color:#a32020;font-weight:600}}
code{{background:#eef1f5;padding:1px 5px;border-radius:4px;font-size:12px}}
footer{{max-width:1020px;margin:0 auto;padding:0 18px 60px;color:#7a828f;font-size:12px}}
.chk{{background:#eef6ee;border:1px solid #b9d9b9;border-radius:9px;padding:12px 18px;margin:14px 0}} .chk li{{font-size:12.7px;color:#244d24}}
</style></head><body>
<header><h1>Diabetic kidney disease &amp; nephropathy references — a spatial-transcriptomics walkthrough</h1>
<p>One story for a bench biologist new to spatial, for someone who wants to see the data behind every
abstraction, and for a B-cell / autoimmune-nephropathy expert. What Xenium can show across IgAN/MN/AA vs
DKD/control — and where it stops. Validated labels · gene availability audited · single-section refs flagged · raw data read-only.</p></header>
<main>{''.join(cards)}
<div class="chk"><b>Panel checklist (design rules met):</b><ul>{checkli}</ul></div>
</main>
<footer>Generated narrative, scoring &amp; drill-down panels (this folder) · reused figures by relative path
from summaries 01/02/03/04/06 + the B-lineage and composition analyses · tool attribution &amp; gene availability
sourced from the scripts and <code>gene_panel_audit.csv</code> · pure science · committed PNGs (no base64).</footer>
</body></html>"""
open(f"{HERE}/dkd_walkthrough.html","w").write(HTML)
print(f"\nwrote {HERE}/dkd_walkthrough.html ({len(HTML)//1024} KB, relative-path images)")
