#!/usr/bin/env python
"""B-lineage MECHANISTIC signatures across nephropathies -- EXTENDS summary 02 (does NOT re-derive it).

ANECDOTAL / DESCRIPTIVE. n=1 per non-DKD condition -> NO statistics, NO p-values, hypothesis-
generating only, suggestive-of-different-disease-process. Validated reannotation labels; raw read-only.
Every Ig/TLS/state readout is conditioned on the correct producer cell and gated above ambient floor.

REUSED FROM 02 (pulled, not recomputed): per_sample_substrate.csv (B-lineage burden, aggregate
counts), aggregate_composition.csv (in/around-aggregate hull composition). The per-participant gallery
and subgroup-burden figure also live in 02 and are NOT re-emitted here.

NEW LAYERS ONLY (per condition, per-sample, never pooled):
 1. B:Plasma ratio  -- split 02's combined B-lineage into B vs Plasma (plasma-skew = antibody-mediated
    MN/IgAN; B-skew = TLS-like DKD).
 2. Ig isotype WITHIN plasma -- IGHG1 (IgG) quantitative; IGHA1 (IgA) reported but FAILS the >=3%
    detection floor in plasma (panel cannot usably resolve IgA). IGHM/IGKC/IGLC NOT on the Xenium panel.
 3. TLS-organization markers in the aggregate region -- CCL19 (T-zone/FRC) quantitative; CXCL13/CXCR5/
    CCR7 presence-only flags (detectable above ambient but <1% prevalence); CR2 dropped (ambient-sat.).
 4. Structural localization -- peri-glomerular / peri-tubular / peri-vascular / interstitial; glomerular
    axis emphasised for the glomerular diseases (MN/IgAN).
 5. Damage x B coupling -- per B / Plasma / Myeloid cell, median distance to nearest injured-tubule vs
    nearest glomerulus. Does B/plasma track injury (myeloid-like) or glomeruli (antibody-mediated)?
 6. B / plasma state -- B naive/memory/activation (TCL1A/BANK1/CD27/CD83); plasma maturity
    (MZB1/TNFRSF17/PRDM1).
"""
import os, warnings
warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree, ConvexHull
from matplotlib.patches import Polygon
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
RE=f"{REPO}/analysis/dkd_xenium_reannotation"; DIS=f"{REPO}/analysis/dkd_xenium_disease"
OUT=DIS; FIG=f"{DIS}/figures"; os.makedirs(FIG,exist_ok=True)
EPS,MINPTS,REGION=50.0,20,50.0          # same aggregate definition as 02/substrate.py
BRICH={"1006","HK2695"}
MIN_PLASMA=15                            # suppress within-plasma readouts below this n
# --- gated usable gene set (see ambient-gate preflight; conditioned on producer) ---
IGG="IGHG1"; IGA="IGHA1"                 # IGHG1 PASS (62% plasma,163x); IGHA1 FAIL (0.97% plasma) -> flag
TLS_QUANT=["CCL19"]                      # PASS (4.5% immune, 52x)
TLS_FLAG=["CXCL13","CXCR5","CCR7"]       # detectable>>ambient but <1% prevalence -> presence-only
BSTATE={"TCL1A":"naive","BANK1":"Bcore","CD27":"memory","CD83":"activ"}
PSTATE=["MZB1","TNFRSF17","PRDM1"]
GENES=[IGG,IGA]+TLS_QUANT+TLS_FLAG+list(BSTATE)+PSTATE
def hdr(s): print("\n"+"="*80+"\n"+s+"\n"+"="*80)
def save(fig,n): fig.savefig(f"{FIG}/{n}.png",dpi=170,bbox_inches="tight"); plt.close(fig); print("  [fig]",n)

# ---------------- load (validated labels + gated genes) ----------------
c=pd.read_parquet(f"{RE}/cells.parquet").reset_index(drop=True)
g=pd.read_parquet(f"{DIS}/bcell_genes_percell.parquet").reset_index(drop=True)
assert (c.cell_id.astype(str).values==g.cell_id.astype(str).values).all(), "row alignment"
for col in GENES: c[col]=g[col].values
samp=c.orig_ident.astype(str).values; lab=c.my_label.values
XY=c[["spatial_x","spatial_y"]].values
isB=lab=="B"; isP=lab=="Plasma"; Blin=isB|isP; isMye=lab=="Myeloid"
GLOM=np.isin(lab,["Podo","MC","EC_glom","PEC"]); VESS=np.isin(lab,["EC_Peritub","EC_DVR","VSMC"])
TUB=np.isin(lab,["PT","iPT","TAL","iTAL","DCT","CNT","PC","IC A"]); INJ=np.isin(lab,["iPT","iTAL"])
def condlabel(s):
    cc=c.loc[samp==s,"Condition"].iloc[0]
    if cc=="DKD": return "DKD B-rich" if s in BRICH else "DKD B-poor"
    return {"IgA":"IgAN"}.get(cc,cc)
ORDER=["1006","HK2695","1001","1008","1010","1011","1012","1013","HK3626","HK2753","HK3106",
       "1003","1005","1004","1009","1007"]

# ---------------- pulled-from-02 substrate (burden context; not recomputed) ----------------
SUB=pd.read_csv(f"{DIS}/per_sample_substrate.csv")
COMP=pd.read_csv(f"{DIS}/aggregate_composition.csv")
sub_idx=SUB.set_index("orig_ident"); comp_idx=COMP.set_index("orig_ident")

def det_in(genecol,cellmask,s):
    m=cellmask&(samp==s)
    return (np.nan if m.sum()==0 else float((c[genecol].values[m]>0).mean()))

content=[]; iso=[]; tls=[]; loc=[]; coup=[]; state=[]; members={}
for s in ORDER:
    sm=samp==s; n=int(sm.sum()); nB=int((isB&sm).sum()); nP=int((isP&sm).sum()); cl_=condlabel(s)
    blin_frac=float(sub_idx.loc[s,"Blin_frac"]); nagg=int(sub_idx.loc[s,"n_agg"])
    # 1. B:Plasma split (NEW; 02 had combined B-lineage)
    content.append(dict(sample=s,condition=cl_,n_cells=n,Blin_frac_x02=round(blin_frac,4),
        B_frac=round(nB/n,4),Plasma_frac=round(nP/n,4),B_pct=round(nB/n*100,3),Plasma_pct=round(nP/n*100,3),
        B_to_Plasma=round(nB/max(nP,1),2),n_B=nB,n_Plasma=nP,n_agg_x02=nagg))
    # 2. Ig isotype within plasma (conditioned on plasma)
    iga=det_in(IGA,isP,s); igg=det_in(IGG,isP,s)
    iso.append(dict(sample=s,condition=cl_,n_plasma=nP,
        pct_plasma_IgG=round((igg or 0)*100,1),pct_plasma_IgA_subgate=round((iga or 0)*100,2),
        note=("n_plasma<%d"%MIN_PLASMA if nP<MIN_PLASMA else "")))
    # 3 + crops: DBSCAN aggregate membership (same def as 02) -> region cells
    Bc=XY[Blin&sm]; mem=None
    if Bc.shape[0]>=MINPTS:
        dbl=DBSCAN(eps=EPS,min_samples=MINPTS).fit(Bc).labels_
        if (dbl!=-1).any(): mem=Bc[dbl!=-1]; members[s]=[Bc[dbl==k] for k in set(dbl) if k!=-1]
    rec=dict(sample=s,condition=cl_)
    if mem is not None:
        reg=np.unique(np.concatenate([np.asarray(t,int) for t in cKDTree(XY[sm]).query_ball_point(mem,r=REGION)]))
        gi=np.where(sm)[0][reg]; nreg=len(gi)
        rec["n_region"]=nreg
        for tg in TLS_QUANT: rec[f"{tg}_pct"]=round(float((c[tg].values[gi]>0).mean())*100,2)
        for tg in TLS_FLAG: rec[f"{tg}_ncells"]=int((c[tg].values[gi]>0).sum())  # presence-only flag
    else:
        rec["n_region"]=0
        for tg in TLS_QUANT: rec[f"{tg}_pct"]=np.nan
        for tg in TLS_FLAG: rec[f"{tg}_ncells"]=0
    tls.append(rec)
    # 4. structural localization (positional; nearest structure <=30um else interstitial)
    trees={k:(cKDTree(XY[m&sm]) if (m&sm).sum()>0 else None) for k,m in [("glom",GLOM),("tub",TUB),("vasc",VESS)]}
    def localize(cellmask):
        idx=np.where(cellmask&sm)[0]
        if len(idx)==0: return dict(periglom=np.nan,peritub=np.nan,perivasc=np.nan,interstitial=np.nan)
        D=np.vstack([(trees[k].query(XY[idx])[0] if trees[k] else np.full(len(idx),np.inf)) for k in ("glom","tub","vasc")]).T
        nearest=np.argmin(D,axis=1); mind=D.min(1); out=np.where(mind>30,3,nearest)
        return dict(periglom=round((out==0).mean(),3),peritub=round((out==1).mean(),3),
                    perivasc=round((out==2).mean(),3),interstitial=round((out==3).mean(),3))
    loc.append(dict(sample=s,condition=cl_,**{f"B_{k}":v for k,v in localize(isB).items()},
                                          **{f"P_{k}":v for k,v in localize(isP).items()}))
    # 5. damage x B coupling -- median dist to injured-tubule vs glomerulus (B, Plasma, Myeloid)
    tinj=cKDTree(XY[INJ&sm]) if (INJ&sm).sum()>0 else None
    tglo=cKDTree(XY[GLOM&sm]) if (GLOM&sm).sum()>0 else None
    def meddist(cellmask):
        idx=np.where(cellmask&sm)[0]
        if len(idx)==0: return (np.nan,np.nan)
        di=tinj.query(XY[idx])[0] if tinj else np.full(len(idx),np.nan)
        dg=tglo.query(XY[idx])[0] if tglo else np.full(len(idx),np.nan)
        return (round(float(np.median(di)),1),round(float(np.median(dg)),1))
    bi,bg=meddist(isB); pi,pg=meddist(isP); mi,mg=meddist(isMye)
    coup.append(dict(sample=s,condition=cl_,B_d_inj=bi,B_d_glom=bg,P_d_inj=pi,P_d_glom=pg,
                     Mye_d_inj=mi,Mye_d_glom=mg))
    # 6. B / plasma state (identity-conditioned, gated)
    st=dict(sample=s,condition=cl_,n_B=nB,n_Plasma=nP)
    for gn,nm in BSTATE.items(): st[f"B_{nm}"]=round((det_in(gn,isB,s) or 0)*100,1)
    for gn in PSTATE: st[f"P_{gn}"]=round((det_in(gn,isP,s) or 0)*100,1)
    state.append(st)

CON=pd.DataFrame(content); CON.to_csv(f"{OUT}/bcell_content.csv",index=False)
ISO=pd.DataFrame(iso);     ISO.to_csv(f"{OUT}/bcell_isotype.csv",index=False)
TLSd=pd.DataFrame(tls);    TLSd.to_csv(f"{OUT}/bcell_tls.csv",index=False)
LOC=pd.DataFrame(loc);     LOC.to_csv(f"{OUT}/bcell_localization.csv",index=False)
COUP=pd.DataFrame(coup);   COUP.to_csv(f"{OUT}/bcell_damage_coupling.csv",index=False)
ST=pd.DataFrame(state);    ST.to_csv(f"{OUT}/bcell_state.csv",index=False)

hdr("1. B:Plasma split  (NEW; 02 had combined B-lineage)"); print(CON.to_string(index=False))
hdr("2. Ig isotype within plasma  (IGHG1 gated PASS; IGHA1 sub-gate, panel cannot resolve IgA)"); print(ISO.to_string(index=False))
hdr("3. TLS markers in aggregate region  (CCL19 quantitative; CXCL13/CXCR5/CCR7 presence-only)"); print(TLSd.to_string(index=False))
hdr("4. B/plasma localization vs structures"); print(LOC.to_string(index=False))
hdr("5. damage x B coupling  (median um to injured-tubule vs glomerulus)"); print(COUP.to_string(index=False))
hdr("6. B/plasma state"); print(ST.to_string(index=False))
hdr("(pulled from 02) aggregate hull composition -- focus rows");
print(COMP[COMP.orig_ident.isin(["1006","HK2695","1003","1005","1004","1007"])].to_string(index=False))

# ============================ FIGURES (NEW LAYERS ONLY) ============================
hdr("FIGURES")
CC={"DKD B-rich":"#6A3D9A","DKD B-poor":"#c9b8de","Control":"#888888","IgAN":"#1f77b4",
    "MN":"#E31A1C","AA amyloid":"#FF7F00","C3GN":"#33A02C"}
xl=[f"{r.sample}\n{r.condition}" for r in CON.itertuples()]; x=np.arange(len(CON))
SUP="ANECDOTAL  ·  n=1 per non-DKD  ·  descriptive  ·  hypothesis-generating  ·  NOT statistically tested"

# FIG A -- B:Plasma split + IgG-within-plasma
fig,axes=plt.subplots(1,2,figsize=(16,5))
ax=axes[0]; ax.bar(x-0.2,CON.B_pct,0.4,label="B %",color="#1f77b4"); ax.bar(x+0.2,CON.Plasma_pct,0.4,label="Plasma %",color="#ff7f0e")
ax.set_xticks(x); ax.set_xticklabels(xl,fontsize=7,rotation=90); ax.set_ylabel("% of section cells"); ax.legend(fontsize=9)
ax.set_title("1. B vs Plasma split  (plasma-skew = antibody-mediated; B-skew = TLS-like)",fontsize=10)
ax=axes[1]
for xi,r in zip(x,ISO.itertuples()):
    if r.n_plasma>=MIN_PLASMA: ax.bar(xi,r.pct_plasma_IgG,0.55,color="#9467bd")
    else: ax.text(xi,2,"n<%d"%MIN_PLASMA,rotation=90,fontsize=6,ha="center",color="#999")
ax.set_xticks(x); ax.set_xticklabels(xl,fontsize=7,rotation=90); ax.set_ylabel("% of plasma IGHG1+ (IgG)")
ax.set_title("2. IgG within plasma  (IgA/IGHA1 below panel detection — cannot resolve IgA axis)",fontsize=10)
fig.suptitle(SUP,fontsize=12,color="#a33"); fig.tight_layout(); save(fig,"bcell_fig_A_split_isotype")

# FIG B -- TLS in aggregate + B/plasma state
fig,axes=plt.subplots(1,2,figsize=(16,5))
ax=axes[0]; tt=TLSd[TLSd.n_region>0]
xt=np.arange(len(tt))
ax.bar(xt,tt["CCL19_pct"],0.5,color="#17becf",label="CCL19 % (T-zone/FRC)")
for i,r in enumerate(tt.itertuples()):
    flags=[f for f,tg in zip(["X","x","r"],TLS_FLAG) if getattr(r,f"{tg}_ncells")>0]
    if flags: ax.text(xt[i],tt["CCL19_pct"].fillna(0).values[i]+0.4,"+".join(flags),fontsize=6,ha="center",color="#444",rotation=90)
ax.set_xticks(xt); ax.set_xticklabels([f"{r.sample}\n{r.condition}" for r in tt.itertuples()],fontsize=7,rotation=90)
ax.set_ylabel("% of aggregate-region cells CCL19+"); ax.legend(fontsize=8)
ax.set_title("3. TLS-organization in aggregate (CCL19; X=CXCL13 x=CXCR5 r=CCR7 present)",fontsize=10)
ax=axes[1]; sm_=ST.copy()
for col,cc,nm in [("B_Bcore","#1f77b4","B BANK1"),("B_memory","#2ca02c","B CD27"),("B_activ","#ff7f0e","B CD83"),
                  ("P_MZB1","#9467bd","Pl MZB1"),("P_TNFRSF17","#8c564b","Pl BCMA")]:
    ax.plot(x,sm_[col],"o-",ms=4,lw=1,label=nm)
ax.set_xticks(x); ax.set_xticklabels(xl,fontsize=7,rotation=90); ax.set_ylabel("% of B / plasma +"); ax.legend(fontsize=7,ncol=2)
ax.set_title("6. B activation / plasma maturity state",fontsize=10)
fig.suptitle(SUP,fontsize=12,color="#a33"); fig.tight_layout(); save(fig,"bcell_fig_B_tls_state")

# FIG C -- glomerular-axis localization + damage coupling
fig,axes=plt.subplots(1,2,figsize=(16,5))
ax=axes[0]; bottom=np.zeros(len(LOC))
for col,cc,nm in [("B_periglom","#E31A1C","peri-glom"),("B_peritub","#c9b8de","peri-tubular"),
                  ("B_perivasc","#8c564b","peri-vascular"),("B_interstitial","#cccccc","interstitial")]:
    ax.bar(x,LOC[col]*100,bottom=bottom,color=cc,label=nm); bottom+=(LOC[col].fillna(0).values*100)
ax.set_xticks(x); ax.set_xticklabels(xl,fontsize=7,rotation=90); ax.set_ylabel("% of B cells"); ax.legend(fontsize=8)
ax.set_title("4. B-cell localization (glomerular axis for MN/IgAN)",fontsize=10)
ax=axes[1]
ax.plot(x,COUP.B_d_inj,"o-",color="#1f77b4",label="B → injured tubule")
ax.plot(x,COUP.B_d_glom,"o--",color="#1f77b4",alpha=0.45,label="B → glomerulus")
ax.plot(x,COUP.Mye_d_inj,"^-",color="#9467bd",label="myeloid → injured tubule")
ax.plot(x,COUP.P_d_glom,"s--",color="#ff7f0e",alpha=0.6,label="plasma → glomerulus")
ax.set_xticks(x); ax.set_xticklabels(xl,fontsize=7,rotation=90); ax.set_ylabel("median distance (µm)"); ax.legend(fontsize=7)
ax.set_title("5. damage×B coupling (B tracks injury [myeloid-like] or glomeruli?)",fontsize=10)
fig.suptitle(SUP,fontsize=12,color="#a33"); fig.tight_layout(); save(fig,"bcell_fig_C_localization_coupling")

# FIG D -- focused glomerular-axis crops (IgAN 1003, MN 1005, DKD B-rich 1006). NOT the 16-panel gallery (that's in 02).
LC={"B":"#1f77b4","Plasma":"#ff7f0e","Myeloid":"#2ca02c"}
fig,axes=plt.subplots(1,3,figsize=(17,6))
for ax,(dis,s) in zip(axes,[("IgAN",'1003'),("MN",'1005'),("DKD B-rich",'1006')]):
    sm=samp==s; ax.scatter(XY[sm,0],XY[sm,1],s=1,c="#ececec",linewidths=0,rasterized=True)
    gm=sm&GLOM; ax.scatter(XY[gm,0],XY[gm,1],s=3,c="#f4a0a0",linewidths=0,rasterized=True,label="glomerulus")
    for nm,cc in LC.items():
        m=sm&(lab==nm); ax.scatter(XY[m,0],XY[m,1],s=(7 if nm!="Myeloid" else 4),c=cc,linewidths=0,rasterized=True,label=nm)
    for P in members.get(s,[]):
        if len(P)>=4:
            try: h=ConvexHull(P); ax.add_patch(Polygon(P[h.vertices],closed=True,fill=False,edgecolor="#08306b",lw=1.2))
            except Exception: pass
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(f"{dis} — {s}",fontsize=12)
axes[0].legend(fontsize=8,loc="upper right",markerscale=2)
fig.suptitle("Glomerular-axis crops — B(blue)/Plasma(orange)/Myeloid(green) + glomeruli(pink) + B-aggregate hulls   ·   ANECDOTAL n=1 (per-participant gallery is in 02)",fontsize=11)
fig.tight_layout(); save(fig,"bcell_fig_D_glom_crops")
print("\n== bcell_nephropathy_anecdotal done ==")
