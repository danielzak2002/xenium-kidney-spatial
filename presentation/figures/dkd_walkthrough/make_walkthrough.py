#!/usr/bin/env python
"""DKD spatial-transcriptomics WALKTHROUGH for a skeptical bench biologist (slides/IHC, new to spatial).

ONE story, layer by layer. Reuses committed figures from summaries 01/02/03/04/06/07 + the new
composition_by_group analysis (referenced by RELATIVE PATH, not base64 -> lightweight & committable).
Generates ONLY narrative + worked-example scoring panels from the validated reannotation labels.

Honest framing is load-bearing: DKD/control is the powered backbone; IgAN/MN/C3GN are SINGLE sections
(reconnaissance, hypothesis-generating, nothing tested); two headline axes are panel-limited; the
credibility comes from blind validation + two displayed self-corrections. Read-only raw; PNG-only.
"""
import os, sys, html as _html
import numpy as np, pandas as pd
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree, ConvexHull
from matplotlib.patches import Polygon
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0,"/Users/danie/ClaudeCode/pilot_analyses/xenium/presentation/figures")
import figstyle; figstyle.apply()

REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
RE=f"{REPO}/analysis/dkd_xenium_reannotation"; DIS=f"{REPO}/analysis/dkd_xenium_disease"
EPI=f"{REPO}/analysis/dkd_epi_endo_stress"
HERE=f"{REPO}/presentation/figures/dkd_walkthrough"; os.makedirs(HERE,exist_ok=True)
def rel(p): return os.path.relpath(p,HERE)            # relative path for committable <img>
def png(name): return f"{HERE}/{name}.png"
def savep(fig,name): fig.savefig(png(name),dpi=170,bbox_inches="tight"); plt.close(fig); print("  [panel]",name)

# coarse-lineage + immune palettes (consistent with the analyses)
LCOL={"Epithelial":"#cfcfcf","Stroma":"#e377c2","Endothelial":"#8c564b","Immune":"#6A3D9A"}
ICOL={"B":"#1f77b4","Plasma":"#ff7f0e","Myeloid":"#9467bd","CD4 T":"#2ca02c","CD8 T":"#17becf"}
INJC="#8c510a"; GLOMC="#E31A1C"

c=pd.read_parquet(f"{RE}/cells.parquet").reset_index(drop=True)
c["sample"]=c.orig_ident.astype(str)
samp=c["sample"].values; lab=c.my_label.values; lin=c.my_lineage.values
XY=c[["spatial_x","spatial_y"]].values
isB=lab=="B"; isP=lab=="Plasma"; Blin=isB|isP; isMye=lab=="Myeloid"
GLOM=np.isin(lab,["Podo","MC","EC_glom","PEC"]); INJ=np.isin(lab,["iPT","iTAL"])

def crop_on_aggregate(s,pad=230):
    """center a crop on the largest B-lineage DBSCAN cluster of sample s."""
    m=samp==s; pts=XY[Blin&m]
    if len(pts)>=20:
        dl=DBSCAN(eps=50,min_samples=20).fit(pts).labels_
        if (dl!=-1).any():
            k=pd.Series(dl[dl!=-1]).value_counts().idxmax(); ctr=pts[dl==k].mean(0)
        else: ctr=pts.mean(0)
    else: ctr=XY[m].mean(0)
    return ctr, (ctr[0]-pad,ctr[0]+pad,ctr[1]-pad,ctr[1]+pad)

def win_mask(s,box):
    x0,x1,y0,y1=box; m=samp==s
    return m & (XY[:,0]>=x0)&(XY[:,0]<=x1)&(XY[:,1]>=y0)&(XY[:,1]<=y1)

# ============================================================
# S1 — slide -> cells -> typed (segmentation + typing bridge), real crop of 1006
# ============================================================
ctr,box=crop_on_aggregate("1006"); w=win_mask("1006",box)
fig,ax=plt.subplots(1,3,figsize=(15,5.2))
ax[0].scatter(XY[w,0],XY[w,1],s=7,c="#b9b9b9",linewidths=0,rasterized=True)
ax[0].set_title("1 · Segmentation\n(cell boundaries drawn on the image)",fontsize=12)
for nm,col in LCOL.items():
    mm=w&(lin==nm); ax[1].scatter(XY[mm,0],XY[mm,1],s=7,c=col,linewidths=0,rasterized=True,label=nm)
ax[1].set_title("2 · Typing\n(identity called from ~5K-gene expression)",fontsize=12); ax[1].legend(fontsize=8,markerscale=1.6,loc="upper right")
ax[2].scatter(XY[w,0],XY[w,1],s=6,c="#e3e3e3",linewidths=0,rasterized=True)
for nm,col in ICOL.items():
    mm=w&(lab==nm); ax[2].scatter(XY[mm,0],XY[mm,1],s=12,c=col,linewidths=0,rasterized=True,label=nm)
ax[2].set_title("3 · The cells we follow\n(immune subtypes, in tissue context)",fontsize=12); ax[2].legend(fontsize=8,markerscale=1.4,loc="upper right")
for a in ax: a.set_aspect("equal"); a.axis("off")
fig.suptitle("Xenium keeps the architecture: every dot is one cell, in place — not a dissociated suspension   (sample 1006, ~460 µm region)",fontsize=13)
savep(fig,"wt_pipeline")

# ============================================================
# S3 — composition object->number worked example (count -> fraction), one Control + 1006
# ============================================================
comp=pd.read_csv(f"{DIS}/composition_by_group.csv")
def lin_fracs(s):
    d=comp[(comp["sample"]==s)&(comp.resolution=="coarse_lineage")].set_index("cell_type").fraction
    return {k:float(d.get(k,0)) for k in ["Epithelial","Immune","Stroma","Endothelial"]}
fig,ax=plt.subplots(1,3,figsize=(15,4.6))
ex="1006"; m=samp==ex
counts={k:int((lin==k).__and__(m).sum()) for k in LCOL}; tot=int(m.sum())
ax[0].axis("off")
txt=("OBJECT -> NUMBER  (composition)\n\nsample 1006\n"
     f"total cells scored = {tot:,}\n\n"+"\n".join(f"{k:12s} {counts[k]:>7,}  ->  {counts[k]/tot*100:5.1f}%" for k in LCOL)+
     "\n\nfraction = (cells of type) / (all cells)\nsums to 100% -> compositional (closure)")
ax[0].text(0.0,0.98,txt,va="top",ha="left",family="monospace",fontsize=11.5)
# stacked bar for 1006
bottom=0
for k in ["Epithelial","Stroma","Endothelial","Immune"]:
    v=counts[k]/tot*100; ax[1].bar(0,v,bottom=bottom,color=LCOL[k],width=0.6,label=k); bottom+=v
ax[1].set_xlim(-0.8,0.8); ax[1].set_xticks([]); ax[1].set_ylabel("% of cells"); ax[1].set_title("1006 (DKD) composition",fontsize=12); ax[1].legend(fontsize=8)
# control contrast
cf=lin_fracs("HK3626"); bottom=0
for k in ["Epithelial","Stroma","Endothelial","Immune"]:
    v=cf[k]*100; ax[2].bar(0,v,bottom=bottom,color=LCOL[k],width=0.6); bottom+=v
ax[2].set_xlim(-0.8,0.8); ax[2].set_xticks([]); ax[2].set_ylabel("% of cells"); ax[2].set_title("HK3626 (Control) — epithelial-dominant, near-bare immune",fontsize=11)
fig.suptitle("Composition is just counting, made transparent: count each type, divide by all cells -> the fractions in S3's dotplots",fontsize=13)
savep(fig,"wt_scoring_composition")

# ============================================================
# S4 — aggregate object->number incl the ZERO-CLUSTER case (B-rich 1006 vs IgAN 1003)
# ============================================================
fig,ax=plt.subplots(1,4,figsize=(19,5.2))
# (1) cells -> B-lineage highlighted
ctr,box=crop_on_aggregate("1006",pad=260); w=win_mask("1006",box)
ax[0].scatter(XY[w,0],XY[w,1],s=6,c="#e0e0e0",linewidths=0,rasterized=True)
bm=w&Blin; ax[0].scatter(XY[bm,0],XY[bm,1],s=16,c="#1f77b4",linewidths=0,rasterized=True)
ax[0].set_title("1 · B-lineage cells\n(blue) among all cells",fontsize=11)
# (2) local density (neighbors within 50um)
pts=XY[bm]; tree=cKDTree(pts); dens=np.array([len(tree.query_ball_point(p,r=50))-1 for p in pts])
sc=ax[1].scatter(pts[:,0],pts[:,1],s=22,c=dens,cmap="viridis",linewidths=0,rasterized=True)
ax[1].set_title("2 · Local density\n(B neighbors within 50 µm)",fontsize=11); plt.colorbar(sc,ax=ax[1],fraction=0.046,pad=0.04)
# (3) DBSCAN core/border/noise + hull
dl=DBSCAN(eps=50,min_samples=20).fit(pts).labels_
ax[2].scatter(pts[dl==-1,0],pts[dl==-1,1],s=14,c="#cccccc",linewidths=0,label="noise")
for k in [kk for kk in set(dl) if kk!=-1]:
    P=pts[dl==k]; ax[2].scatter(P[:,0],P[:,1],s=18,c="#08519c",linewidths=0)
    if len(P)>=4:
        h=ConvexHull(P); ax[2].add_patch(Polygon(P[h.vertices],closed=True,fill=False,edgecolor="#08306b",lw=1.6))
ax[2].set_title("3 · DBSCAN cluster + hull\n(eps=50 µm, minPts=20)",fontsize=11)
# (4) burden number + threshold + zero-cluster IgAN
ax[3].axis("off")
sub=pd.read_csv(f"{DIS}/per_sample_substrate.csv"); sub["sample"]=sub.orig_ident.astype(str); sub=sub.set_index("sample")
b1006=sub.loc["1006","agg_cells_per10k"]; b1003=sub.loc["1003","agg_cells_per10k"]
txt=("OBJECT -> NUMBER  (organization)\n\n"
     "hull = one B-aggregate (the scored object)\n"
     "burden = clustered B-lineage cells per 10k\n\n"
     f"1006 (DKD)   burden = {b1006:.0f}/10k\n"
     f"threshold (natural gap) = 75/10k\n"
     f"   {b1006:.0f} ≥ 75  ->  B-RICH\n\n"
     "SAME ALGORITHM, IgAN 1003:\n"
     f"   burden = {b1003:.0f}/10k  ->  ZERO clusters\n"
     "   (B cells present, but not organized)")
ax[3].text(0.0,0.98,txt,va="top",ha="left",family="monospace",fontsize=11.5)
# inset: IgAN 1003 B-lineage, no hull
m3=samp=="1003"; p3=XY[Blin&m3]
iax=ax[3].inset_axes([0.46,0.0,0.54,0.40])
iax.scatter(XY[m3,0],XY[m3,1],s=1,c="#eeeeee",linewidths=0,rasterized=True)
iax.scatter(p3[:,0],p3[:,1],s=5,c="#1f77b4",linewidths=0,rasterized=True)
iax.set_xticks([]); iax.set_yticks([]); iax.set_title("IgAN 1003 — 0 aggregates",fontsize=8)
for a in ax[:3]: a.set_aspect("equal"); a.axis("off")
fig.suptitle("How a B-aggregate becomes a number — and the honest negative control: the SAME unbiased algorithm finds ZERO clusters in IgAN/controls",fontsize=12.5)
savep(fig,"wt_scoring_aggregate")

# ============================================================
# S5 — distance object->number (B->glomerulus vs myeloid->injury), 1006 crop
# ============================================================
coup=pd.read_csv(f"{DIS}/bcell_damage_coupling.csv").set_index("sample")
ctr,box=crop_on_aggregate("1006",pad=300); w=win_mask("1006",box)
fig,ax=plt.subplots(1,2,figsize=(14,6))
a=ax[0]
a.scatter(XY[w,0],XY[w,1],s=5,c="#ededed",linewidths=0,rasterized=True)
gm=w&GLOM; im=w&INJ; bm=w&isB; mm=w&isMye
a.scatter(XY[gm,0],XY[gm,1],s=12,c=GLOMC,linewidths=0,label="glomerulus")
a.scatter(XY[im,0],XY[im,1],s=12,c=INJC,linewidths=0,label="injured tubule")
a.scatter(XY[bm,0],XY[bm,1],s=22,c="#1f77b4",linewidths=0,label="B cell")
a.scatter(XY[mm,0],XY[mm,1],s=22,c="#9467bd",linewidths=0,label="myeloid")
# nearest-neighbor lines for a few B->glom and myeloid->injury
if gm.sum() and bm.sum():
    tg=cKDTree(XY[gm]); bsel=np.where(bm)[0][:14]
    for i in bsel:
        d,j=tg.query(XY[i]); g=XY[gm][j]; a.plot([XY[i,0],g[0]],[XY[i,1],g[1]],c="#1f77b4",lw=0.7,alpha=0.6)
if im.sum() and mm.sum():
    ti=cKDTree(XY[im]); msel=np.where(mm)[0][:14]
    for i in msel:
        d,j=ti.query(XY[i]); g=XY[im][j]; a.plot([XY[i,0],g[0]],[XY[i,1],g[1]],c="#9467bd",lw=0.7,alpha=0.6)
a.set_aspect("equal"); a.axis("off"); a.legend(fontsize=9,markerscale=1.4,loc="upper right")
a.set_title("Each cell -> distance to nearest landmark (1006 region)",fontsize=12)
# bar of the per-sample medians (from committed CSV)
a=ax[1]
vals=[coup.loc["1006","B_d_glom"],coup.loc["1006","B_d_inj"],coup.loc["1006","Mye_d_inj"]]
labs=["B -> glomerulus","B -> injured tubule","myeloid -> injured tubule"]
cols=["#1f77b4","#1f77b4","#9467bd"]; al=[1,0.45,1]
for i,(v,cc,aa) in enumerate(zip(vals,cols,al)): a.barh(i,v,color=cc,alpha=aa)
a.set_yticks(range(3)); a.set_yticklabels(labs,fontsize=11); a.invert_yaxis()
a.set_xlabel("median distance (µm)")
for i,v in enumerate(vals): a.text(v+1.5,i,f"{v:.0f} µm",va="center",fontsize=11)
a.set_title("Object -> number: B sits near GLOMERULI; myeloid sits near INJURY\n(1006 — two orthogonal spatial logics)",fontsize=11.5)
fig.suptitle("Two spatial programs, measured the same simple way — distance to the nearest landmark cell",fontsize=13)
savep(fig,"wt_scoring_distance")

print("panels generated.\n")

# ============================================================
# HTML — narrative, relative-path images (committable, no base64)
# ============================================================
def img(path,cap=None):
    s=f'<img src="{rel(path)}" loading="lazy">'
    return s+(f'<p class="cap">{cap}</p>' if cap else "")

R1=f"{RE}/figures"; D=f"{DIS}/figures"
SEC=[]
def sec(pid,title,html): SEC.append((pid,title,html))

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
axis of IgAN) is <b>not measurable</b> on this panel; the follicular-TLS marker CXCL13 is <b>presence-only</b>.</li>
<li><b>The credibility comes from two things a slide reader can audit:</b> our cell calls are
<b>validated blind</b> against the published atlas, and we <b>caught and reversed two of our own
over-claims</b> (shown in S7, not hidden).</li>
</ul></div>
<p class="cap">Plain-language bridges used throughout — <b>segmentation</b> = drawing each cell's boundary on
the image; <b>typing</b> = calling identity from expression; <b>aggregate / DBSCAN</b> = an unbiased way to
find cell clusters (shown beside the image so you can see it's a real lymphoid aggregate); <b>ambient /
spillover</b> = transcripts misassigned across cell boundaries, the failure mode we actively control for.</p>""")

sec("S1","The data, for a slide person", f"""
<p>Xenium is <b>in-situ</b>: the tissue is never dissociated. Each transcript is imaged at its location,
assigned to a segmented cell, and ~5,000 genes are measured per cell <b>with the architecture intact</b> —
the spatial equivalent of a multiplexed IHC you can re-stain 5,000 ways. That is the whole reason a
glomerulus stays a glomerulus and an immune aggregate stays an aggregate.</p>
{img(png('wt_pipeline'),"From image to cells we can follow: (1) segmentation draws boundaries; (2) typing calls identity from the ~5K-gene profile; (3) the immune subtypes we trace downstream, still in tissue context. One ~460 µm region of sample 1006.")}
<p class="cap">Cohort: 16 sections — 8 DKD, 3 control, plus IgAN / MN / AA-amyloid (×2) / C3GN as references.
Large image and per-transcript files stay on disk; we work from the cell-by-gene matrix + coordinates.</p>""")

sec("S2","Trust the calls: blind re-annotation, validated against the atlas", f"""
<p>Before any biology: are the cell-type calls real? We re-integrated <b>within this atlas</b> (Harmony on
sample) and re-annotated the whole object <b>blind</b> to the authors' labels, then compared. <b>This is
NOT the separate 68-gene three-cohort integration</b> — that coarser cross-platform exercise is kept apart;
here we use the full ~5K panel.</p>
{img(R1+'/concordance_matrices.png',"Our blind labels vs the authors' published annotation. Segment-level ARI 0.78, immune ARI 0.68 — strong agreement, achieved without looking at their calls.")}
{img(R1+'/umap_yours_vs_theirs.png',"Same cells, our typing vs theirs, on the integrated embedding — the populations land in the same places.")}
<p><b>Takeaway a bench biologist can hold onto:</b> "their published calls and our independent, blind calls
agree" — so everything downstream rests on validated identities, not on our say-so.</p>""")

sec("S3","Top-level composition, by disease group", f"""
<p>Start where IHC starts: what is the tissue made of? Every dot is <b>one section</b>, grouped by disease;
a median bar is drawn <b>only where n&gt;1</b>, so the single IgAN/MN/C3GN dots can't masquerade as a group.</p>
{img(png('wt_scoring_composition'),"Object → number for composition: count cells of each type, divide by all cells. Transparent, and the basis of every dot below.")}
{img(D+'/composition_by_group_coarse.png',"Coarse lineage. Control is epithelial-dominant (74%) with a near-bare immune compartment (2.8%); every disease loses epithelium and gains immune + stroma. The CLR row (closure-aware) tracks the raw ordering.")}
{img(D+'/composition_by_group_immune.png',"Immune drill-down. Single-section standouts: MN most immune/plasma-skewed, C3GN most T-skewed, AA highest B-fraction. NK/DC are not separately typed in the validated labels — not invented.")}
<p class="cap"><b>Closure noted:</b> fractions sum to one, so they are compositional; we read them descriptively
and show a CLR sensitivity view rather than testing them.</p>""")

sec("S4","Organization, not just composition — finding aggregates", f"""
<p>Composition says <i>how much</i> immune; it doesn't say whether the immune cells are <b>organized</b>. To
ask that without cherry-picking, we run <b>DBSCAN</b> — an unbiased clustering that either finds a dense cell
cluster or doesn't. Here is the object→number, shown beside the image so you can see the hull is a real
lymphoid aggregate — and, crucially, the <b>negative control</b>: the same algorithm on IgAN finds nothing.</p>
{img(png('wt_scoring_aggregate'),"B-lineage cells → local density → DBSCAN core/hull (the scored object) → burden per 10k → threshold. Sample 1006 scores 202/10k (B-rich); the SAME algorithm on IgAN 1003 returns ZERO clusters — B cells are present but not organized.")}
{img(D+'/dkd_subgroup.png',"Within DKD, aggregate burden has a clean ~3× gap → a B-rich subgroup = {1006, HK2695}. This reproduces the authors' B-cell-rich subgroup: 8/8 concordant with their B-predominant immune niche.")}
{img(D+'/b_lineage_gallery_16.png',"Every section, grouped: B (blue) + Plasma (orange) with aggregate hulls. B-rich DKD show compact follicular-looking aggregates; IgAN and controls have none under the identical algorithm.")}""")

sec("S5","Two spatial programs (the centerpiece — robust within a section)", f"""
<p>This is where spatial earns its keep, and it holds at n=1 because it is measured <b>within each section</b>
(every cell relative to its own neighbors). Two <b>orthogonal</b> logics:</p>
{img(png('wt_scoring_distance'),"Distance object→number: in 1006, B cells sit ~34 µm from glomeruli but ~94 µm from injured tubules, while myeloid sit ~30 µm from injury. B-lineage organizes near GLOMERULI; myeloid tracks INJURY.")}
{img(EPI+'/coloc_fig_A_specificity.png',"Across sections, tubular injury co-localizes with MYELOID infiltration (ρ = 0.82 [0.46, 0.95]); the B-lineage association does not survive de-circularization (partial 0.13, CI spans 0).")}
{img(EPI+'/coloc_fig_C2_immune_injury.png',"The near-injury immune compartment is myeloid-skewed in 12/16 sections; B-lineage is depleted near injury — consistent with the distance picture.")}
{img(D+'/bcell_fig_C_localization_coupling.png',"The same B-vs-myeloid distance contrast across all samples (from the B-lineage analysis): B-to-glomerulus ≪ B-to-injury; myeloid hugs injury everywhere.")}
<p><b>One sentence:</b> injury recruits myeloid; B-lineage organizes around glomeruli/aggregates — two
different spatial programs in the same kidney.</p>""")

sec("S6","B-lineage programs across nephropathies (reconnaissance)", f"""
<p>Layering mechanism onto the single-section references gives <b>three apparent B-lineage programs</b>.
These are anecdotes to test, not results — and the panel limits are stated exactly where they bite.</p>
{img(D+'/bcell_fig_A_split_isotype.png',"Splitting B vs Plasma: MN is the one plasma-skewed glomerular disease with IgG-dominant plasma. LIMIT: IgA (the defining IgAN isotype) is below panel detection in plasma — the IgA axis is NOT testable here; only the IgG arm is.")}
{img(D+'/bcell_fig_B_tls_state.png',"TLS organization + cell state: only DKD B-rich (HK2695) carries substantial follicular CXCL13. LIMIT: CXCL13 is below the quantitative floor → reported presence-only; CCL19 is the quantitative TLS marker.")}
{img(D+'/bcell_fig_D_glom_crops.png',"Region crops — IgAN / MN / DKD B-rich. DKD B-rich = organized, B-skewed, follicular-flagged; MN = plasma-skewed, IgG-dominant, mixed/T-zone; IgAN = diffuse, peri-vascular, no aggregate, mature plasma.")}
<div class="frame warn"><b>n = 1 per non-DKD condition. Descriptive. Hypothesis-generating. NOT tested.</b>
The IgAN read rests on architecture + plasma maturity, NOT on isotype (which is unmeasurable on this panel).</div>""")

sec("S7","How we avoided fooling ourselves (two self-corrections)", f"""
<p><b>Ambient / spillover</b> in one line: transcripts occasionally get assigned to the wrong neighboring
cell across a segmentation boundary, creating a faint smear of signal where it isn't really produced. If you
don't control for it, you "discover" things that are just spillover. We did — and it cost us two claims.</p>
{img(D+'/baff_ambient_control.png',"Self-correction 1: we initially read a peri-aggregate STROMAL BAFF niche. Under an ambient stress-test it did not survive — non-producer epithelium and endothelium rise just like stroma near aggregates → a local spillover field, not stromal production. The localized-niche claim is RETRACTED.")}
{img(D+'/baff_myeloid_anchor.png',"What survives (the canonical BAFF result, 03/04): a tissue-wide MYELOID producer, cell-intrinsic and reproducible (16/16), ~24× the epithelial floor per-transcriptome — but with NO aggregate-specific niche. APRIL is a NO-GO (sub-floor).")}
{img(EPI+'/coloc_fig_B1_disease.png',"Self-correction 2: an earlier 'tubular injury sits near B-aggregates' reading. De-circularized, the injury association is MYELOID, not B-lineage (S5). The B-specific claim is RETRACTED; the myeloid one is stated in its own right.")}
<p><b>Why this matters to a skeptic:</b> the analysis visibly disproves two of its own attractive stories.
That same discipline is what makes the surviving claims — validated typing, the two spatial programs —
trustworthy.</p>""")

sec("S8","Synthesis, honest limits, and the study this motivates", f"""
<p><b>What spatial transcriptomics CAN do here:</b> (1) reproduce published cell identities <b>blind</b>
(ARI 0.78/0.68); (2) quantify composition and <b>organization</b> (DBSCAN aggregates) per section; (3) reveal
<b>two orthogonal spatial programs</b> — injury↔myeloid, B-lineage↔glomeruli — robust within a section.</p>
<p><b>Where it STOPS (stated, not hidden):</b></p>
<ul>
<li><b>Single-section disease claims are hypotheses.</b> IgAN/MN/C3GN n=1; AA n=2; one patient per section; no donor column → no biological replication for the references.</li>
<li><b>Panel limits.</b> IgA isotype unmeasurable (IGHA1 sub-floor; IGHM/IGKC/IGLC absent) → the defining IgAN axis can't be tested; CXCL13 presence-only → follicular-TLS is flagged, not quantified.</li>
<li><b>Typing caveats.</b> iPT (injured PT) recall ≈ 0.64 → injury distances are conservative; CD8 recall ~0.58.</li>
<li><b>Associational.</b> All co-localization is correlational within fixed tissue — no dynamics, no causation.</li>
</ul>
<p><b>The powered study this motivates:</b> a multi-donor cohort per nephropathy (IgAN, MN, AA, C3GN) on a
panel carrying the missing isotype + follicular markers (IGHA/IGHM/IGKC, CXCL13/CXCR5/CR2), to TEST the three
candidate B-lineage programs and the injury↔myeloid vs B↔glomeruli dissociation that this reconnaissance generated.</p>
<table class="meth"><tr><th>step</th><th>tool</th><th>what it does</th></tr>
<tr><td>typing / integration</td><td>scanpy · harmonypy</td><td>within-atlas Harmony(sample) + blind re-annotation (full ~5K panel)</td></tr>
<tr><td>validation</td><td>scikit-learn</td><td>ARI / confusion vs the authors' published labels</td></tr>
<tr><td>aggregates</td><td>scikit-learn DBSCAN · scipy</td><td>eps=50 µm, minPts=20; cKDTree neighborhoods; ConvexHull</td></tr>
<tr><td>co-localization</td><td>scipy · numpy</td><td>Spearman + bootstrap CIs; de-circularized partials; within-section permutation</td></tr>
<tr><td>composition / state</td><td>pandas · numpy</td><td>per-section fractions; CLR sensitivity; producer-conditioned, ambient-gated detection</td></tr>
<tr><td>figures</td><td>matplotlib · figstyle</td><td>slide-grade panels; region crops; committed PNG</td></tr></table>""")

# panel checklist (design rule)
CHECK=["object→number scoring shown for composition (S3), aggregates (S4), distance (S5)",
 "individual replicates grouped by condition (Control/DKD/IgAN/MN/AA/C3GN) in S3 & S4",
 "all spatial zooms are region-crops (S1, S4, S5, S6)",
 "aggregate determination incl. the ZERO-cluster IgAN case (S4)",
 "two self-corrections DISPLAYED — BAFF niche + damage→B (S7)",
 "panel limits displayed where they bite — IgA untestable, CXCL13 presence-only (S6, S8)",
 "one canonical BAFF (03/04 myeloid anchor; not 02's earlier panel)",
 "numbers sourced from committed CSVs/REPORTs",
 "lightweight: images referenced by RELATIVE PATH (no base64) → committable"]
print("PANEL CHECKLIST"); [print(f"  [x] {x}") for x in CHECK]

cards=[]
for pid,title,body in SEC:
    cls="card title" if pid=="S0" else "card"
    cards.append(f'<section class="{cls}"><h2><span class="pid">{pid}</span> {_html.escape(title)}</h2>{body}</section>')
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
.frame ul{{margin:7px 0 0;padding-left:20px}} ul{{margin:6px 0}}
table.meth{{border-collapse:collapse;width:100%;font-size:12.7px;margin:10px 0}}
table.meth th,table.meth td{{border-bottom:1px solid #eceef2;padding:5px 9px;text-align:left}} table.meth th{{color:#6a7180}}
footer{{max-width:1020px;margin:0 auto;padding:0 18px 60px;color:#7a828f;font-size:12px}}
.chk{{background:#eef6ee;border:1px solid #b9d9b9;border-radius:9px;padding:12px 18px;margin:14px 0}} .chk li{{font-size:12.7px;color:#244d24}}
</style></head><body>
<header><h1>Diabetic kidney disease & nephropathy references — a spatial-transcriptomics walkthrough</h1>
<p>One story for a bench biologist new to spatial: what Xenium can show about how IgAN, MN and AA differ from
DKD and controls — and where the method stops. Validated labels · descriptive · single-section references flagged · raw data read-only.</p></header>
<main>{''.join(cards)}
<div class="chk"><b>Panel checklist (design rules met):</b><ul>{checkli}</ul></div>
</main>
<footer>Generated narrative + scoring panels (this folder) · reused figures referenced by relative path from
summaries 01/02/03/04/06 + the B-lineage and composition analyses · pure science · figures are committed PNGs
(no base64) · numbers sourced from the committed CSVs/REPORTs.</footer>
</body></html>"""
open(f"{HERE}/dkd_walkthrough.html","w").write(HTML)
print(f"\nwrote {HERE}/dkd_walkthrough.html ({len(HTML)//1024} KB, relative-path images)")
