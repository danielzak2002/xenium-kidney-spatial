#!/usr/bin/env python
"""STAGE 2b — Readout A on the reliable-set labels + DIRECT cohort-difference bootstrap test."""
import os, warnings; warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import sys; sys.path.insert(0,"/Users/danie/ClaudeCode/pilot_analyses/xenium/presentation/figures")
import figstyle as fs
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"; OUT=f"{REPO}/analysis/three_cohort_integration"
EPS=50.0; MINPTS=20; R=50.0; COH=["RCC_big","RCC_figshare","DKD"]
COL={"RCC_big":"#1F78B4","RCC_figshare":"#6BAED6","DKD":"#6A3D9A","RCC_pooled":"#08519c"}
df=pd.read_parquet(f"{OUT}/cells_labeled_reliable.parquet")

def aggregates(cd):
    rows=[]
    for s,g in cd.groupby("sample"):
        xy=g[["x","y"]].values; isB=g.is_B.values
        if isB.sum()<MINPTS: continue
        cl=DBSCAN(eps=EPS,min_samples=MINPTS).fit(xy[isB]).labels_
        mem0=np.where(isB)[0][cl!=-1]
        if len(mem0)==0: continue
        tree=cKDTree(xy); bgT=g.is_Treg.mean(); bgC=g.is_cyto.mean(); tr=g.is_Treg.values; cy=g.is_cyto.values
        for c in [c for c in np.unique(cl) if c!=-1]:
            mem=np.where(isB)[0][cl==c]
            reg=np.unique(np.concatenate([np.asarray(t,int) for t in tree.query_ball_point(xy[mem],r=R)]))
            rows.append((s,len(reg),int(tr[reg].sum()),int(cy[reg].sum()),float(bgT),float(bgC)))
    return pd.DataFrame(rows,columns=["section","n","Treg_in","cyto_in","bgT","bgC"])

def arrs(a): return a.Treg_in.values,(a.n*a.bgT).values,a.cyto_in.values,(a.n*a.bgC).values
def dl(tin,texp,ein,eexp): return np.log2((tin.sum()+1e-9)/(texp.sum()+1e-9))-np.log2((ein.sum()+1e-9)/(eexp.sum()+1e-9))

print("=== READOUT A (reliable-set labels) ===")
AG={c:aggregates(df[df.cohort==c]) for c in COH}
AG["RCC_pooled"]=pd.concat([AG["RCC_big"],AG["RCC_figshare"]],ignore_index=True)
rowsA=[]
for c in ["RCC_big","RCC_figshare","RCC_pooled","DKD"]:
    a=AG[c]; tin,texp,ein,eexp=arrs(a); d=dl(tin,texp,ein,eexp)
    r=np.random.default_rng(0); n=len(a); bs=[dl(*[v[r.integers(0,n,n)] for v in (tin,texp,ein,eexp)]) for _ in range(5000)]
    lo,hi=np.percentile(bs,[2.5,97.5])
    rowsA.append(dict(cohort=c,n_agg=n,n_sec=a.section.nunique(),delta=round(d,3),lo=round(lo,3),hi=round(hi,3),fold=round(2**d,2)))
    print(f"  {c:13s}: {n} aggs/{a.section.nunique()} sec | Δ {d:+.2f} [{lo:+.2f},{hi:+.2f}] (~{2**d:.1f}x)")
RA=pd.DataFrame(rowsA); RA.to_csv(f"{OUT}/readoutA_reliable.csv",index=False)

# ---- DIRECT cohort-difference bootstrap (resample aggregates within cohort) ----
def diff_boot(ca,cb,nb=5000,seed=1):
    a=AG[ca]; b=AG[cb]; ta,xa,ea,fa=arrs(a); tb,xb,eb,fb=arrs(b)
    r=np.random.default_rng(seed); na=len(a); nb_=len(b); pt=dl(ta,xa,ea,fa)-dl(tb,xb,eb,fb); bs=[]
    for _ in range(nb):
        ia=r.integers(0,na,na); ib=r.integers(0,nb_,nb_)
        bs.append(dl(ta[ia],xa[ia],ea[ia],fa[ia])-dl(tb[ib],xb[ib],eb[ib],fb[ib]))
    lo,hi=np.percentile(bs,[2.5,97.5]); return pt,lo,hi
print("\n=== DIRECT cohort difference (Δ_cohortA − Δ_DKD), bootstrap 95% CI ===")
rowsD=[]
for ca in ["RCC_pooled","RCC_big","RCC_figshare"]:
    pt,lo,hi=diff_boot(ca,"DKD"); excl = (lo>0)or(hi<0)
    rowsD.append(dict(contrast=f"{ca} − DKD",diff=round(pt,3),lo=round(lo,3),hi=round(hi,3),excludes_zero=excl))
    print(f"  {ca:13s} − DKD: {pt:+.2f} [{lo:+.2f},{hi:+.2f}]  excludes 0: {excl}")
RD=pd.DataFrame(rowsD); RD.to_csv(f"{OUT}/readoutA_difference.csv",index=False)

# ---- figure: per-cohort Δ + difference CIs ----
fig,axes=plt.subplots(1,2,figsize=(14,4.6))
ax=axes[0]; order=["RCC_big","RCC_figshare","RCC_pooled","DKD"]; yp={c:i for i,c in enumerate(order[::-1])}
for _,r in RA.iterrows():
    y=yp[r.cohort]; ax.errorbar(r.delta,y,xerr=[[r.delta-r.lo],[r.hi-r.delta]],fmt="o",ms=12,color=COL[r.cohort],capsize=6,lw=2.4,markeredgecolor="black",zorder=3)
    ax.text(r.delta,y+0.2,f"{r.delta:+.2f}",ha="center",fontsize=10,color=COL[r.cohort],fontweight="bold")
fs.zeroline(ax,0,"v"); ax.set_yticks(list(yp.values())); ax.set_yticklabels(list(yp.keys())); ax.set_ylim(-0.5,len(order)-0.5)
ax.set_xlabel("Δlog₂ (Treg − cytotoxic) per aggregate"); ax.set_title("Per-cohort differential (reliable-set labels)")
ax=axes[1]; yp2={r.contrast:i for i,r in RD[::-1].reset_index(drop=True).iterrows()}
for _,r in RD.iterrows():
    y=yp2[r.contrast]; c="#08519c" if "pooled" in r.contrast else "#3182bd"
    ax.errorbar(r["diff"],y,xerr=[[r["diff"]-r.lo],[r.hi-r["diff"]]],fmt="s",ms=11,color=c,capsize=6,lw=2.4,markeredgecolor="black",zorder=3)
    ax.text(r["diff"],y+0.18,f"{r['diff']:+.2f} {'(excl 0)' if r.excludes_zero else '(incl 0)'}",ha="center",fontsize=9,fontweight="bold")
fs.zeroline(ax,0,"v"); ax.set_yticks(list(yp2.values())); ax.set_yticklabels(list(yp2.keys())); ax.set_ylim(-0.5,len(RD)-0.5)
ax.set_xlabel("difference in Δlog₂  (cohort − DKD), bootstrap 95% CI"); ax.set_title("DIRECT cohort-difference test")
fig.suptitle("Readout A · reliable-set replication + direct cohort-difference test",fontsize=14)
fs.save_fig(fig,"READOUT_A_reliable")
print("== stage2b done ==")
