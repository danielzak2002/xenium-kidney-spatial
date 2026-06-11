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
    ax.set_title("Readout B is a PANEL LIMITATION, not a result:\nANGPT2 / CXCL9 / CXCL10 are sub-ambient on DKD Xenium → not measurable cross-cohort",fontsize=11)
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
def build_html():
 PANELS=[
 ("P0",None,"Does the RCC immunoregulatory B-aggregate bias generalize across kidney disease?",
   "Three Xenium cohorts (two ccRCC, one DKD), harmonized on a reliably-detected gene set, asked whether the Treg-favoring B-aggregate bias reproduces."),
 ("P1",f"{OUT}/P1_cohorts.png","Meet the cohorts","Three Xenium kidney datasets — all imaged on Xenium, integrated to ask one question."),
 ("P2",f"{OUT}/P2_funnel.png","The gene funnel (honest reduction)","123 name-intersection → 104 measured on DKD Xenium → 68 reliably above ambient in all three. The cost is stated up front; most drops are DKD-driven sparsity."),
 ("P3a",f"{PF}/INTrel_umap.png","Integration is sound — embedding","Harmony UMAP on the 68-gene set: the three cohorts MIX (left), lineages RESOLVE (middle), and the two RCC cohorts merge as a batch baseline (right)."),
 ("P3b",f"{PF}/INTrel_dotplot.png","Integration is sound — typing","Per-cohort marker dot-plot: each lineage is defined by the right markers in EVERY cohort — typing is consistent, not blended."),
 ("P4",f"{OUT}/P4_lineage.png","Lineage availability (what reliability cost)","Typing stays marker-faithful, but NK collapses to KLRD1-only and the Treg gate loses IL2RA — stated, not hidden."),
 ("P5",f"{OUT}/P5_spatial.png","The B-aggregate object","One representative section per cohort, region-cropped: DBSCAN B-aggregates (blue, hulls) with Treg (red) and cytotoxic (green) cells in and around them — the spatial substrate of Readout A."),
 ("P6",f"{OUT}/P6_pillars.png","Readout A — replication vs separation (the core)","Two complementary pillars, NOT the same claim. LEFT (integrated, reliable set): the Treg-favoring direction REPLICATES across two independent RCC cohorts, but the direct cohort−DKD difference CIs INCLUDE zero (underpowered at this depth — CIs shown, not hidden). RIGHT (native depth): the magnitude/separation (Δ≈+2.6 vs +0.24, non-overlapping CIs)."),
 ("P7",f"{OUT}/P7_readoutB.png","Readout B — a panel limitation","ANGPT2/CXCL9/CXCL10 are sub-ambient on the DKD Xenium panel → the vascular/inflammatory axis is not measurable cross-cohort; the within-RCC signal is null. Reported as a limitation, not a result."),
 ("P8",f"{OUT}/P8_synthesis.png","Synthesis","What generalizes (conserved scaffold + replicated Treg direction) vs what is established only at native depth (the magnitude/separation)."),
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
    P1(); P2(); P4(); P5(); P6(); P7(); P8()
    build_html()   # assemble last so it embeds fresh figures
    print("== walkthrough done ==")
