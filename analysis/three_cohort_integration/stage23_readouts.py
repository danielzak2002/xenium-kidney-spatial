#!/usr/bin/env python
"""STAGE 2+3 — cross-cohort readouts, computed PER COHORT on native cells/coords + harmonized labels.
A: immunoregulatory aggregate differential (Treg vs cytotoxic), two RCC cohorts vs DKD.
B: endothelial/inflammatory stress near vs far from B-aggregates (cross-context axis).
Reads cells_labeled.parquet (coords+labels+stress counts only). Reuses bniche DBSCAN/diff code."""
import os, warnings; warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import sys; sys.path.insert(0,"/Users/danie/ClaudeCode/pilot_analyses/xenium/presentation/figures")
import figstyle as fs
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"; OUT=f"{REPO}/analysis/three_cohort_integration"
EPS=50.0; MINPTS=20; R=50.0; FAR=200.0; rng=np.random.default_rng(0)
COH=["RCC_big","RCC_figshare","DKD"]; COL={"RCC_big":"#1F78B4","RCC_figshare":"#6BAED6","DKD":"#6A3D9A"}
df=pd.read_parquet(f"{OUT}/cells_labeled.parquet")
print("cells:",len(df),"| per cohort:",dict(df.groupby('cohort').size()))

# ---- aggregates per (cohort, section); store per-agg composition + member coords ----
def cohort_aggregates(cd):
    aggs=[]   # per-aggregate dicts; members per section for readout B
    sect_members={}
    for s,g in cd.groupby("sample"):
        xy=g[["x","y"]].values; isB=g.is_B.values
        if isB.sum()<MINPTS: continue
        cl=DBSCAN(eps=EPS,min_samples=MINPTS).fit(xy[isB]).labels_
        members=np.where(isB)[0][cl!=-1]
        if len(members)==0: continue
        sect_members[s]=xy[members]
        tree=cKDTree(xy)
        bgT=g.is_Treg.mean(); bgC=g.is_cyto.mean(); treg=g.is_Treg.values; cyto=g.is_cyto.values
        for c in [c for c in np.unique(cl) if c!=-1]:
            mem=np.where(isB)[0][cl==c]
            reg=np.unique(np.concatenate([np.asarray(t,int) for t in tree.query_ball_point(xy[mem],r=R)]))
            aggs.append(dict(section=s,n=len(reg),Treg_in=int(treg[reg].sum()),cyto_in=int(cyto[reg].sum()),
                             bgT=float(bgT),bgC=float(bgC)))
    return pd.DataFrame(aggs),sect_members

# ============================================================================
# READOUT A — differential
# ============================================================================
def pooled_dlog2(a):
    tin=a.Treg_in.values; texp=(a.n*a.bgT).values; ein=a.cyto_in.values; eexp=(a.n*a.bgC).values
    return float(np.log2((tin.sum()+1e-9)/(texp.sum()+1e-9))-np.log2((ein.sum()+1e-9)/(eexp.sum()+1e-9)))
def boot(a,nb=5000,seed=0):
    r=np.random.default_rng(seed); n=len(a); bs=[]
    tin=a.Treg_in.values;texp=(a.n*a.bgT).values;ein=a.cyto_in.values;eexp=(a.n*a.bgC).values
    for _ in range(nb):
        i=r.integers(0,n,n)
        bs.append(np.log2((tin[i].sum()+1e-9)/(texp[i].sum()+1e-9))-np.log2((ein[i].sum()+1e-9)/(eexp[i].sum()+1e-9)))
    return np.percentile(bs,2.5),np.percentile(bs,97.5)
print("\n=== READOUT A: Treg-vs-cytotoxic differential per cohort ===")
AGG={}; rowsA=[]
for c in COH:
    a,mem=cohort_aggregates(df[df.cohort==c]); AGG[c]=(a,mem)
    if len(a)==0: print(f"  {c}: no aggregates"); continue
    d=pooled_dlog2(a); lo,hi=boot(a)
    rowsA.append(dict(cohort=c,n_aggregates=len(a),n_sections=a.section.nunique(),
                      delta_log2=round(d,3),ci_lo=round(lo,3),ci_hi=round(hi,3),fold=round(2**d,2)))
    print(f"  {c}: {len(a)} aggs / {a.section.nunique()} sec | Δlog2 {d:+.2f} [{lo:+.2f},{hi:+.2f}] (~{2**d:.1f}x)")
RA=pd.DataFrame(rowsA); RA.to_csv(f"{OUT}/readoutA_differential.csv",index=False)

fig,ax=plt.subplots(figsize=(9,4.4))
ypos={c:i for i,c in enumerate(COH[::-1])}
for _,r in RA.iterrows():
    y=ypos[r.cohort]; ax.errorbar(r.delta_log2,y,xerr=[[r.delta_log2-r.ci_lo],[r.ci_hi-r.delta_log2]],
        fmt="o",ms=13,color=COL[r.cohort],capsize=7,lw=2.5,markeredgecolor="black",zorder=3)
    ax.text(r.delta_log2,y+0.18,f"Δ {r.delta_log2:+.2f} (~{r.fold:.0f}×)",ha="center",fontsize=11,color=COL[r.cohort],fontweight="bold")
fs.zeroline(ax,0,"v"); ax.set_yticks(list(ypos.values())); ax.set_yticklabels(list(ypos.keys()))
ax.set_ylim(-0.5,len(COH)-0.5); ax.set_xlabel("Δlog₂ = log₂(Treg) − log₂(cytotoxic)  per B-aggregate (count-pooled, bootstrap 95% CI)")
ax.set_title("Readout A · immunoregulatory bias replicated across TWO RCC cohorts vs DKD")
fs.save_fig(fig,"READOUT_A")

# ============================================================================
# READOUT B — endothelial / inflammatory stress near aggregates
# ============================================================================
print("\n=== READOUT B: endothelial/inflammatory stress near B-aggregates ===")
STRESS=["ANGPT2","CXCL9","CXCL10","HLA-DRA"]; ENDOID=["PECAM1","VWF"]
ACT=["ANGPT2"]; INFL=["CXCL9","CXCL10","HLA-DRA"]
# usability gate per cohort (detection in endothelial cells >= 2%)
gate=[]; usable={}
for c in COH:
    e=df[(df.cohort==c)&df.is_endo]
    usable[c]={}
    for g in STRESS+ENDOID:
        col=f"cnt_{g}"; dr=float((e[col]>0).mean()) if col in df else 0.0
        usable[c][g]= dr>=0.02
        gate.append(dict(cohort=c,gene=g,endo_detect=round(dr,4),usable=dr>=0.02))
GATE=pd.DataFrame(gate); GATE.to_csv(f"{OUT}/readoutB_usability.csv",index=False)
print(GATE.pivot(index="gene",columns="cohort",values="endo_detect").to_string())

def module_score(sub,genes,cohort):
    gs=[g for g in genes if usable[cohort].get(g)]
    if not gs: return np.full(len(sub),np.nan)
    z=np.zeros(len(sub))
    for g in gs:
        v=np.log1p(sub[f"cnt_{g}"].values.astype(float)); sd=v.std()
        z+=(v-v.mean())/sd if sd>0 else 0
    return z/len(gs)

rowsB=[]; gradB=[]
for c in COH:
    cd=df[df.cohort==c]; _,mem=AGG[c]
    endo=cd[cd.is_endo].copy()
    if len(endo)<50 or not mem: continue
    # distance of each endo cell to nearest aggregate member, per section
    endo["dist"]=np.inf
    for s,M in mem.items():
        sel=endo["sample"].values==s
        if sel.sum()==0: continue
        d,_=cKDTree(M).query(endo.loc[sel,["x","y"]].values); endo.loc[sel,"dist"]=d
    for mod,genes in [("endo_activation",ACT),("inflammatory",INFL)]:
        sc=module_score(endo,genes,c); endo[f"m_{mod}"]=sc
        near=endo.dist<=R; far=endo.dist>FAR
        sn=sc[near.values]; sf=sc[far.values]; sn=sn[~np.isnan(sn)]; sf=sf[~np.isnan(sf)]
        if len(sn)<10 or len(sf)<10:
            rowsB.append(dict(cohort=c,module=mod,n_near=len(sn),n_far=len(sf),delta=np.nan)); continue
        # per-section matched delta
        deltas=[]
        for s in endo["sample"].unique():
            m=endo["sample"].values==s; nn=sc[(m&near.values)]; ff=sc[(m&far.values)]
            nn=nn[~np.isnan(nn)]; ff=ff[~np.isnan(ff)]
            if len(nn)>=10 and len(ff)>=10: deltas.append(np.nanmean(nn)-np.nanmean(ff))
        rowsB.append(dict(cohort=c,module=mod,n_near=len(sn),n_far=len(sf),
            delta=round(float(np.nanmean(sn)-np.nanmean(sf)),3),
            k_sec_pos=int(np.sum(np.array(deltas)>0)),n_sec=len(deltas)))
        # gradient
        for lo,hi,lab in [(0,50,"0-50"),(50,100,"50-100"),(100,200,"100-200"),(200,500,"200-500"),(500,1e9,">500")]:
            mm=(endo.dist>=lo)&(endo.dist<hi); v=sc[mm.values]; v=v[~np.isnan(v)]
            gradB.append(dict(cohort=c,module=mod,bin=lab,mean=round(float(np.mean(v)),3) if len(v) else np.nan,n=len(v)))
RB=pd.DataFrame(rowsB); RB.to_csv(f"{OUT}/readoutB_near_far.csv",index=False)
GR=pd.DataFrame(gradB); GR.to_csv(f"{OUT}/readoutB_gradient.csv",index=False)
print("\nnear vs far endo stress (delta = near-far):"); print(RB.to_string(index=False))

fig,axes=plt.subplots(1,2,figsize=(13,5))
mods=["endo_activation","inflammatory"]; x=np.arange(len(COH)); w=0.36
for mi,mod in enumerate(mods):
    vals=[RB[(RB.cohort==c)&(RB.module==mod)].delta.iloc[0] if len(RB[(RB.cohort==c)&(RB.module==mod)]) else np.nan for c in COH]
    axes[0].bar(x+(mi-0.5)*w,[v if v==v else 0 for v in vals],w,label=mod,color=["#E69F00","#CC79A7"][mi])
fs.zeroline(axes[0],0,"h"); axes[0].set_xticks(x); axes[0].set_xticklabels(COH,rotation=15); axes[0].set_ylabel("near − far (module z)")
axes[0].set_title("endothelial stress near vs far B-aggregates"); axes[0].legend(frameon=False)
blab=["0-50","50-100","100-200","200-500",">500"]; xb=np.arange(len(blab))
for c in COH:
    sub=GR[(GR.cohort==c)&(GR.module=="inflammatory")]
    if len(sub): axes[1].plot(xb,[sub[sub.bin==b]["mean"].iloc[0] if len(sub[sub.bin==b]) else np.nan for b in blab],marker="o",lw=2,color=COL[c],label=c)
fs.zeroline(axes[1],0,"h"); axes[1].set_xticks(xb); axes[1].set_xticklabels(blab); axes[1].set_xlabel("distance to B-aggregate (µm)")
axes[1].set_ylabel("inflammatory module z"); axes[1].set_title("inflammatory gradient"); axes[1].legend(frameon=False,fontsize=9)
fig.suptitle("Readout B · endothelial/inflammatory stress vs B-aggregate proximity (cross-context axis)",fontsize=14)
fs.save_fig(fig,"READOUT_B")
print("== stage23 done ==")
