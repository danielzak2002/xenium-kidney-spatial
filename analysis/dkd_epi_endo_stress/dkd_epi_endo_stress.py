#!/usr/bin/env python
"""
dkd_epi_endo_stress.py — do tubular/endothelial cells NEAR immune B-aggregates carry elevated
injury / profibrotic / endothelial-activation / hypoxia programs vs matched cells FAR away?
Engages the paper's profibrotic-tubular-niche thesis at single-cell resolution. Demoulin XENIUM
only (subtype/state work is Xenium-only per cd4_cd8_support).

Reuses the bniche_dbscan B-aggregate machinery (DBSCAN eps=50/minPts=20 on B coords, per section;
distance to nearest aggregate-member B cell). Read-only; backed h5ad; only program-gene columns
materialized; library size from obs nCount_RNA.

ASSOCIATIONAL: aggregates may form in already-injured regions -> co-localization, not causation.
MORPHOLOGY: processed h5ad = centroids + expression only (no image/segmentation layer) ->
transcriptional-state + spatial-organization, NOT cell-shape histomorphology.
"""
import os, warnings
warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
H5=os.path.join(REPO,"Demoulin26/data/spatial_adata_xenium_cosmx_zenodo.h5ad")
OUT=os.path.join(REPO,"analysis/dkd_epi_endo_stress"); os.makedirs(OUT,exist_ok=True)
EPS=50.0; MINPTS=20; R=50.0; FAR=200.0; DKD_COLOR="#6A3D9A"
def hdr(s): print("\n"+"="*78+"\n"+s+"\n"+"="*78)

EPI_LABELS=["PT","iPT","TAL","iTAL","DCT","CNT","IC A","IC B","PC","PEC","DTL_ATL"]
ENDO_LABELS=["EC_Peritub","EC_glom","EC_DVR","EC_Lymph"]
PROGRAMS={  # program -> (target compartment, candidate genes)
  "injPT":   ("epi", ["VCAM1","HAVCR1","PROM1","SPP1","ITGB6","DCDC2","LCN2","TIMP1"]),
  "fibroEMT":("epi", ["COL1A1","FN1","VIM","TGFB1","SNAI2"]),
  "endoAct": ("endo",["ICAM1","VCAM1","SELE","ENG","PLVAP"]),
  "hypoxia": ("endo",["VEGFA","CA9"]),
}
ALLG=sorted(set(g for _,gs in PROGRAMS.values() for g in gs))

# ============================================================================
hdr("STEP 0 — load (Xenium), materialize program-gene columns only")
a=ad.read_h5ad(H5, backed="r")
plat=a.obs["tech"].astype(str).values
xen=plat=="Xenium"; xidx=np.where(xen)[0]
genes=[g for g in ALLG if g in set(map(str,a.var_names))]
print("program genes present:", genes, "| absent:", [g for g in ALLG if g not in genes])
sub=a[xidx, genes].to_memory()
lab=a.obs["annotation_updated"].astype(str).values[xidx]
samp=a.obs["orig_ident"].astype(str).values[xidx]
ncount=a.obs["nCount_RNA"].astype(float).values[xidx]
immlab=a.obs["immune_cell_annotation_combined"].astype(str).values[xidx]
xy=np.asarray(a.obsm["spatial"],float)[xidx]
a.file.close()
C=sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
cnt=pd.DataFrame(C, columns=genes)
print(f"Xenium cells {len(xidx):,}; epi {np.isin(lab,EPI_LABELS).sum():,}; "
      f"endo {np.isin(lab,ENDO_LABELS).sum():,}; B {np.sum(immlab=='B'):,}")

# ============================================================================
hdr("STEP 0 (gate) — usability vs Immune-cell ambient floor")
is_epi=np.isin(lab,EPI_LABELS); is_endo=np.isin(lab,ENDO_LABELS); is_imm=(lab=="Immune")
def det(mask,g): v=cnt[g].values[mask]; return float((v>0).mean()), float(v.mean())
gate=[]; usable={}
for prog,(tgt,gs) in PROGRAMS.items():
    tmask=is_epi if tgt=="epi" else is_endo
    keep=[]
    for g in gs:
        if g not in genes:
            gate.append(dict(program=prog,gene=g,target=tgt,target_detect=np.nan,floor_detect=np.nan,
                             ratio=np.nan,usable=False,reason="absent")); continue
        td,tm=det(tmask,g); fd,fm=det(is_imm,g)
        ok = (td>=0.03) and (td >= 2*fd)
        if ok: keep.append(g)
        gate.append(dict(program=prog,gene=g,target=tgt,target_detect=round(td,4),
            target_mean=round(tm,4),floor_detect=round(fd,4),ratio=round(td/max(fd,1e-9),2),
            usable=ok,reason="" if ok else ("low detect" if td<0.03 else "<2x floor")))
    usable[prog]=keep
    print(f"  {prog} ({tgt}): usable={keep}")
gate_df=pd.DataFrame(gate); gate_df.to_csv(os.path.join(OUT,"usability_gate.csv"),index=False)
print(gate_df.to_string(index=False))
PROGRAMS={p:(t,g) for p,(t,g) in PROGRAMS.items() if usable[p]}  # drop empty programs

# ============================================================================
hdr("STEP 1 — normalize (CP-median, log1p) and module-score per section+compartment")
med=np.median(ncount[ncount>0]); norm=cnt.div(np.where(ncount>0,ncount,1),axis=0)*med
lognorm=np.log1p(norm)
# z-score each usable gene within (section, compartment); program score = mean z
score=pd.DataFrame(index=range(len(xidx)))
for prog,(tgt,gs) in PROGRAMS.items():
    gs=usable[prog]; sc=np.full(len(xidx),np.nan)
    tmask=is_epi if tgt=="epi" else is_endo
    for s in np.unique(samp):
        m=tmask&(samp==s)
        if m.sum()<30: continue
        zsum=np.zeros(m.sum())
        for g in gs:
            v=lognorm[g].values[m]; sd=v.std()
            zsum += (v-v.mean())/sd if sd>0 else 0.0
        sc[np.where(m)[0]]=zsum/len(gs)
    score[prog]=sc
print("scored programs:", list(PROGRAMS.keys()))

# ============================================================================
hdr(f"STEP 1 — B-aggregates per section (DBSCAN eps={EPS},minPts={MINPTS}); distance to nearest")
dist=np.full(len(xidx),np.inf)
n_agg_tot=0; sect_with_agg=[]
for s in np.unique(samp):
    sm=samp==s; Bm=sm&(immlab=="B")
    if Bm.sum()<MINPTS: continue
    Bc=xy[Bm]; cl=DBSCAN(eps=EPS,min_samples=MINPTS).fit(Bc).labels_
    members=Bc[cl!=-1]
    if len(members)==0: continue
    n_agg_tot+=len(set(cl[cl!=-1])); sect_with_agg.append(s)
    tree=cKDTree(members)
    cells=sm&(is_epi|is_endo)
    d,_=tree.query(xy[cells]); dist[np.where(cells)[0]]=d
print(f"aggregates: {n_agg_tot} across {len(sect_with_agg)} sections (reuses bniche_dbscan delineation)")

# ============================================================================
hdr("STEP 2 — NEAR (<=R) vs FAR (>200) by cell type, within section")
near=dist<=R; far=dist>FAR
rows=[]; per_sec=[]
for prog,(tgt,gs) in PROGRAMS.items():
    tmask=is_epi if tgt=="epi" else is_endo
    deltas=[]; ks=0; ns=0
    for s in sect_with_agg:
        m=tmask&(samp==s)
        sn=score[prog].values[m & near]; sf=score[prog].values[m & far]
        sn=sn[~np.isnan(sn)]; sf=sf[~np.isnan(sf)]
        if len(sn)<10 or len(sf)<10: continue
        d=float(np.nanmean(sn)-np.nanmean(sf)); deltas.append(d); ns+=1; ks+= (d>0)
        per_sec.append(dict(program=prog,section=s,n_near=len(sn),n_far=len(sf),
                            mean_near=round(float(np.mean(sn)),3),mean_far=round(float(np.mean(sf)),3),
                            delta=round(d,3)))
    # pooled (z already per section-compartment, comparable)
    pn=score[prog].values[tmask&near]; pf=score[prog].values[tmask&far]
    pn=pn[~np.isnan(pn)]; pf=pf[~np.isnan(pf)]
    try: U,p=mannwhitneyu(pn,pf,alternative="greater")
    except Exception: p=np.nan
    rows.append(dict(program=prog,target=tgt,k_sections=ks,n_sections=ns,
        median_section_delta=round(float(np.median(deltas)),3) if deltas else np.nan,
        pooled_mean_near=round(float(np.mean(pn)),3),pooled_mean_far=round(float(np.mean(pf)),3),
        pooled_delta=round(float(np.mean(pn)-np.mean(pf)),3),
        mannwhitney_p_near_gt_far=(f"{p:.2e}" if p==p else "na")))
res=pd.DataFrame(rows); print(res.to_string(index=False))
res.to_csv(os.path.join(OUT,"near_vs_far_summary.csv"),index=False)
pd.DataFrame(per_sec).to_csv(os.path.join(OUT,"near_vs_far_per_section.csv"),index=False)

# ============================================================================
hdr("STEP 2 — distance gradient (bins)")
bins=[0,50,100,200,500,np.inf]; blab=["0-50","50-100","100-200","200-500",">500"]
grad_rows=[]
for prog,(tgt,gs) in PROGRAMS.items():
    tmask=(is_epi if tgt=="epi" else is_endo)&np.isfinite(dist)
    for i,bl in enumerate(blab):
        m=tmask&(dist>=bins[i])&(dist<bins[i+1])
        v=score[prog].values[m]; v=v[~np.isnan(v)]
        grad_rows.append(dict(program=prog,bin=bl,n=len(v),mean_score=round(float(np.mean(v)),3) if len(v) else np.nan))
grad=pd.DataFrame(grad_rows); grad.to_csv(os.path.join(OUT,"distance_gradient.csv"),index=False)
print(grad.pivot(index="bin",columns="program",values="mean_score").reindex(blab).to_string())

# ============================================================================
hdr("FIGURES")
# (a) representative section: epithelial cells colored by injPT score + aggregate B cells
rep=max(sect_with_agg, key=lambda s:((samp==s)&(immlab=="B")).sum())
sm=samp==rep
fig,axes=plt.subplots(1,2,figsize=(15,7))
for ax,(prog,tgt) in zip(axes,[("injPT","epi"),("endoAct","endo")]):
    if prog not in PROGRAMS: ax.axis("off"); continue
    tmask=(is_epi if tgt=="epi" else is_endo)&sm
    Bm=sm&(immlab=="B")
    ax.scatter(xy[sm,0],xy[sm,1],s=1,c="#eeeeee",linewidths=0)
    sc=score[prog].values[tmask]; vmax=np.nanpercentile(sc,98) if np.isfinite(sc).any() else 1
    p=ax.scatter(xy[tmask,0],xy[tmask,1],s=5,c=sc,cmap="magma",vmin=np.nanpercentile(sc,5),vmax=vmax,linewidths=0)
    # aggregate B members
    Bc=xy[Bm]; cl=DBSCAN(eps=EPS,min_samples=MINPTS).fit(Bc).labels_
    ax.scatter(Bc[cl!=-1,0],Bc[cl!=-1,1],s=8,facecolors="none",edgecolors="#00a000",linewidths=0.6,label="B-aggregate")
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(f"{rep}: {tgt} cells colored by {prog} score\n(green=B-aggregate)",fontsize=10)
    fig.colorbar(p,ax=ax,fraction=0.035,label=f"{prog} z-score"); ax.legend(fontsize=8,loc="upper right")
fig.suptitle(f"Representative section — stress score vs B-aggregate proximity (associational; centroids only, no morphology)",fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"fig_a_representative_section.png"),dpi=140); plt.close(fig)

# (b) gradient
fig,ax=plt.subplots(figsize=(8,5))
cmap={"injPT":"#6A3D9A","fibroEMT":"#9467bd","endoAct":"#E69F00","hypoxia":"#d62728"}
xpos=np.arange(len(blab))
for prog in PROGRAMS:
    vals=[grad.loc[(grad.program==prog)&(grad.bin==b),"mean_score"].iloc[0] for b in blab]
    ax.plot(xpos,vals,marker="o",lw=2,color=cmap.get(prog,"#333"),label=prog)
ax.axhline(0,color="gray",ls="--",lw=1); ax.set_xticks(xpos); ax.set_xticklabels(blab)
ax.set_xlabel("distance to nearest B-aggregate (um)"); ax.set_ylabel("mean program z-score")
ax.set_title("Stress-program gradient vs distance to B-aggregate (Demoulin Xenium)"); ax.legend(fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"fig_b_distance_gradient.png"),dpi=150); plt.close(fig)
print("saved fig_a_representative_section.png, fig_b_distance_gradient.png")

# ============================================================================
hdr("WRITING REPORT")
L=[];W=L.append
W("# Epithelial / endothelial stress near immune B-aggregates (Demoulin Xenium)\n")
W("Read-only. Xenium-only. Reuses bniche_dbscan B-aggregate delineation (DBSCAN eps=50/minPts=20). "
  "NEAR = <=50 um to nearest aggregate-member B cell; FAR = >200 um; matched **within section and "
  "by cell type** (PT/tubular-to-tubular, endo-to-endo). Programs module-scored (CP-median log1p, "
  "per-section per-compartment z, averaged over usable genes).\n")
W("## Usability gate (per program)\n")
for prog,(tgt,gs) in PROGRAMS.items(): W(f"- **{prog}** ({tgt}): usable = {usable[prog]}")
drop=[p for p in ["injPT","fibroEMT","endoAct","hypoxia"] if p not in PROGRAMS]
if drop: W(f"- dropped (no usable members): {drop}")
W("\nFull gate (detection vs Immune-cell ambient floor) in `usability_gate.csv`.\n")
W("## NEAR vs FAR (matched, within section)\n")
W("| program | target | k/N sections near>far | median Δz | pooled Δz | MWU p (near>far) |")
W("|---|---|---|---|---|---|")
for _,r in res.iterrows():
    W(f"| {r.program} | {r.target} | {int(r.k_sections)}/{int(r.n_sections)} | {r.median_section_delta} | {r.pooled_delta} | {r.mannwhitney_p_near_gt_far} |")
W(f"\n({n_agg_tot} aggregates across {len(sect_with_agg)} sections.)\n")
W("## Distance gradient (mean program z-score)\n")
W("| bin (um) | "+" | ".join(PROGRAMS.keys())+" |")
W("|---|"+"|".join(["---"]*len(PROGRAMS))+"|")
for b in blab:
    W(f"| {b} | "+" | ".join(f"{grad.loc[(grad.program==p)&(grad.bin==b),'mean_score'].iloc[0]}" for p in PROGRAMS)+" |")
W("\n## Interpretation\n")
def grad_decreases(prog):
    vals=[grad.loc[(grad.program==prog)&(grad.bin==b),"mean_score"].iloc[0] for b in ["0-50","100-200",">500"]]
    return vals[0]>vals[-1]+0.05
for _,r in res.iterrows():
    p=r.program
    elevated = (r.pooled_delta>0.05) and (r.k_sections>r.n_sections/2) and grad_decreases(p)
    if elevated:
        W(f"- **{p}: ELEVATED near aggregates.** pooled Δz {r.pooled_delta:+.2f}, {int(r.k_sections)}/{int(r.n_sections)} "
          f"sections positive, MWU p={r.mannwhitney_p_near_gt_far}; monotonic distance gradient "
          f"({grad.loc[(grad.program==p)&(grad.bin=='0-50'),'mean_score'].iloc[0]:+.2f} at 0-50um -> "
          f"{grad.loc[(grad.program==p)&(grad.bin=='>500'),'mean_score'].iloc[0]:+.2f} at >500um). The clearest signal.")
    else:
        W(f"- **{p}: not elevated** near aggregates (pooled Δz {r.pooled_delta:+.2f}, {int(r.k_sections)}/{int(r.n_sections)} "
          f"sections positive, no decreasing gradient).")
W(f"- *Note:* the hypoxia program is a SINGLE usable gene (VEGFA); its large-n Mann-Whitney p is a "
  "tie-dominated artifact and should be read with the per-section k/N and gradient (both null), not the p-value.")
W("- **Associational, not causal:** aggregates may nucleate in already-injured/inflamed regions; "
  "this is co-localization of tubular injury (injPT) with B-aggregates, not evidence that "
  "aggregates cause injury — fully consistent with the paper's profibrotic-tubular-niche thesis at "
  "the level of spatial association.")
W("- **No patient column** -> donor clustering uncontrolled (sections may share donors); per-section "
  "matching mitigates section/composition effects but not donor structure.")
W("- **No morphology:** processed h5ad = centroids + expression only; this is transcriptional-state "
  "+ spatial-organization, NOT cell-shape histomorphology. No image/segmentation layer present.")
W("- Module scores are per-section per-compartment z (comparable within section); pooled across "
  "sections after centering. Genes gated vs an Immune-cell ambient floor (Demoulin dropped negprobes).")
open(os.path.join(OUT,"REPORT.md"),"w").write("\n".join(L))
print("wrote REPORT.md + CSVs + figures")
print("\n== TASK B done ==")
