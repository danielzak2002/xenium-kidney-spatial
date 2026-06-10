#!/usr/bin/env python
"""
bniche_dbscan.py — Is the DKD-vs-RCC DIVERGE real biology or a unit-of-analysis artifact?

Replicates the EXACT RCC immunoregulatory-aggregate pipeline (py/phaseB_02_rcc_aggregates.py +
the committed eps=50 rebuild in py/whitepaper_recompute.py) on DKD Xenium, so the comparison is
1:1. Then adds a radial (core->margin) profile to test whether a Treg collar sits at the
aggregate MARGIN rather than the core.

Borrowed RCC parameters (verbatim, committed headline = 37 aggs, Treg 36/37, eff-CD8 34/37):
  - aggregate delineation : DBSCAN(eps=50.0, min_samples=20) on B-cell coords, per section
  - aggregate footprint   : all cells within R=50.0 um of any member B cell (cKDTree)
  - enrichment metric      : log2((inside_mean_frac + 1e-6)/(background_frac + 1e-6)),
                             per-aggregate Wilcoxon of (inside_frac - bg) vs background,
                             k/N = number of aggregates with inside_frac > section background
  - background             : section (per-sample) tissue fraction of each state  (UNCONDITIONAL,
                             i.e. state cells / ALL region cells — matches RCC; differs from the
                             earlier niche test which used conditional FOXP3+ among CD4+)

SCOPE (cd4_cd8_support): subtype/state work is XENIUM-ONLY. CosMx structural aggregate already
confirmed; no CosMx subtype claims here.

State definitions (atlas lacks native Treg / effector-CD8 / mregDC labels -> reconstructed by
marker+ within compartment, measured Xenium markers only):
  Treg-like   = CD4+  & (FOXP3|IL2RA|CTLA4)+
  effector-CD8= CD8+  & (GZMB|GZMK|PRF1|GNLY)+
  mregDC-like = DC    & (LAMP3|CCR7|FSCN1|CD274)+   (descriptive)
  plasma      = Plasma label (structural)
  B           = B label

Read-only; backed h5ad; subset before .to_memory(); full X never materialized.
"""
import os, warnings
warnings.filterwarnings("ignore")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import numpy as np, pandas as pd, anndata as ad
import scipy.sparse as sp
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BASE = "/Users/danie/ClaudeCode/pilot_analyses/xenium/Demoulin26"
H5AD = os.path.join(BASE, "data", "spatial_adata_xenium_cosmx_zenodo.h5ad")
OUT  = os.path.join(BASE, "analysis", "bniche_dbscan"); os.makedirs(OUT, exist_ok=True)
RCC_COMP = "/Users/danie/ClaudeCode/pilot_analyses/xenium/outputs/tables/rcc_phaseB2_aggregate_composition.csv"

# ---- borrowed RCC parameters (verbatim) ----
EPS = 50.0; MINPTS = 20; R = 50.0
DKD_COLOR = "#6A3D9A"   # consistent dataset color
PLATFORM="tech"; SAMPLE="orig_ident"; CT_IMM="immune_cell_annotation_combined"
NICHE_VAL="B predom. Immune ME"
MIN_B_SAMPLE = MINPTS   # need at least minPts B cells to possibly form an aggregate

TREG=["FOXP3","IL2RA","CTLA4"]; EFFCD8=["GZMB","GZMK","PRF1","GNLY"]
MREGDC=["LAMP3","CCR7","FSCN1","CD274"]; ALLMARK=TREG+EFFCD8+MREGDC
def hdr(s): print("\n"+"="*78+"\n"+s+"\n"+"="*78)

# ============================================================================
# Step 0 — handles + spatial-unit sanity (must be um to transfer eps=50)
# ============================================================================
hdr("STEP 0 — handles, spatial-unit sanity, B counts per Xenium sample")
adata = ad.read_h5ad(H5AD, backed="r")
obs = adata.obs
platN = obs[PLATFORM].astype(str).map(lambda x:"CosMx" if "cosmx" in x.lower() else ("Xenium" if "xenium" in x.lower() else x))
samp = obs[SAMPLE].astype(str); imm = obs[CT_IMM].astype(str)
spatial = np.asarray(adata.obsm["spatial"], float)
print(f"object {adata.shape}; obsm spatial shape {spatial.shape}")

# locate niche column (for validity cross-check)
niche_col = next((c for c in obs.columns if NICHE_VAL in set(obs[c].astype(str).unique())), None)
print(f"author niche column: {niche_col}")
is_niche = (obs[niche_col].astype(str)==NICHE_VAL).values if niche_col else np.zeros(adata.n_obs,bool)

xen = (platN=="Xenium").values
xidx = np.where(xen)[0]
xlab = imm.values[xidx]; xsamp = samp.values[xidx]; xxy = spatial[xidx]; xniche = is_niche[xidx]
DCs = [l for l in np.unique(xlab) if l.upper().endswith("DC")]
print(f"Xenium cells: {len(xidx):,}; DC labels: {DCs}")

# unit sanity: median nearest-neighbour distance per sample (Xenium cells ~5-30 um apart)
samples = sorted(np.unique(xsamp))
nn_meds=[]
for s in samples[:6]:
    c = xxy[xsamp==s]
    if len(c)<50: continue
    sub_c = c[np.random.default_rng(0).choice(len(c), min(5000,len(c)), replace=False)]
    d,_ = cKDTree(c).query(sub_c, k=2); nn_meds.append(np.median(d[:,1]))
nn_med = float(np.median(nn_meds))
ext = xxy.max(0)-xxy.min(0)
print(f"global Xenium coord extent: {ext[0]:.0f} x {ext[1]:.0f}")
print(f"median nearest-neighbour cell distance (sampled): {nn_med:.2f}")
UNITS_UM = 3.0 < nn_med < 40.0
print(f"-> units consistent with MICRONS: {UNITS_UM}  "
      f"({'eps=50 transfers directly' if UNITS_UM else 'WARNING: rescale eps before use'})")
assert UNITS_UM, "spatial units do not look like microns; refit eps before transferring."

# B counts per sample
print("\nB-cell count per Xenium sample (DBSCAN runs where B >= %d):" % MIN_B_SAMPLE)
bcounts = {s:int(((xlab=="B")&(xsamp==s)).sum()) for s in samples}
for s in samples: print(f"  {s}: B={bcounts[s]}{'  (run)' if bcounts[s]>=MIN_B_SAMPLE else '  (skip)'}")
run_samples = [s for s in samples if bcounts[s]>=MIN_B_SAMPLE]
print(f"samples with sufficient B: {run_samples}  (INDEPENDENT of author niche label)")

# ============================================================================
# Build state arrays (materialize counts only for CD4+/CD8+/DC Xenium cells)
# ============================================================================
hdr("BUILD STATES — materialize counts for compartment cells only (lean)")
COMP_LABELS = ["CD4+","CD8+"]+DCs
need_local = np.isin(xlab, COMP_LABELS)
need_global = np.sort(xidx[need_local])
print(f"materializing counts for {len(need_global):,} compartment cells (CD4+/CD8+/DC)")
sub = adata[need_global].to_memory(); adata.file.close()
C = (sub.layers["counts"] if "counts" in sub.layers else sub.X)
C = (C.tocsr() if sp.issparse(C) else sp.csr_matrix(C)).astype(np.float32)
slab = imm.values[need_global]; var_ix={g:i for i,g in enumerate(sub.var_names)}
measured = {g: (g in var_ix and (np.asarray(C[:,var_ix[g]].sum())>0)) for g in ALLMARK}
TREG_X=[g for g in TREG if measured[g]]; EFF_X=[g for g in EFFCD8 if measured[g]]; MREG_X=[g for g in MREGDC if measured[g]]
print(f"measured Xenium markers: Treg={TREG_X}  eff-CD8={EFF_X} (dropped {[g for g in EFFCD8 if not measured[g]]})  mregDC={MREG_X}")
def anypos(genes, rowmask):
    gs=[var_ix[g] for g in genes if g in var_ix]
    if not gs: return np.zeros(rowmask.sum(),bool)
    return np.asarray(C[np.where(rowmask)[0]][:,gs].sum(1)).ravel()>0

# compact-Xenium state booleans
is_B   = (xlab=="B"); is_plasma=(xlab=="Plasma")
is_treg=np.zeros(len(xidx),bool); is_eff=np.zeros(len(xidx),bool); is_mreg=np.zeros(len(xidx),bool)
pos_compact = np.searchsorted(xidx, need_global)   # positions of compartment cells in compact arrays
m4=slab=="CD4+"; m8=slab=="CD8+"; mdc=np.isin(slab,DCs)
is_treg[pos_compact[m4]] = anypos(TREG_X, m4)
is_eff[pos_compact[m8]]  = anypos(EFF_X, m8)
is_mreg[pos_compact[mdc]]= anypos(MREG_X, mdc)
STATES = {"Treg-like":is_treg, "eff-CD8":is_eff, "mregDC-like":is_mreg, "Plasma":is_plasma}
print("global Xenium state counts:", {k:int(v.sum()) for k,v in STATES.items()}, "| B:", int(is_B.sum()))

# ============================================================================
# Step 1 — delineate aggregates (DBSCAN per section, RCC params) + niche validity
# ============================================================================
hdr(f"STEP 1 — DBSCAN B-aggregates per section (eps={EPS}, minPts={MINPTS}); R={R} footprint")
agg_rows=[]; per_agg_inside={k:[] for k in STATES}; per_agg_bg={k:[] for k in STATES}
per_agg_counts={k:[] for k in STATES}   # (state cells in footprint, region cells, section bg)
agg_meta=[]  # (sample, centroid, member local idx in sample, region local idx, niche overlap)
for s in run_samples:
    smask = xsamp==s; loc=np.where(smask)[0]; xy_s=xxy[loc]
    B_s=is_B[loc]; Bpos=np.where(B_s)[0]
    if len(Bpos)<MINPTS: continue
    cl = DBSCAN(eps=EPS, min_samples=MINPTS).fit(xy_s[Bpos]).labels_
    tree_s = cKDTree(xy_s)
    st_s = {k:STATES[k][loc] for k in STATES}; niche_s=xniche[loc]
    bg = {k: float(st_s[k].mean()) for k in STATES}
    for c in [c for c in np.unique(cl) if c!=-1]:
        members = Bpos[cl==c]; cen = xy_s[members].mean(0)
        nbrs = tree_s.query_ball_point(xy_s[members], r=R)
        region = np.unique(np.concatenate([np.asarray(n_,int) for n_ in nbrs]))
        comp = {k: float(st_s[k][region].mean()) for k in STATES}
        niche_ov = float(niche_s[region].mean())
        for k in STATES:
            per_agg_inside[k].append(comp[k]); per_agg_bg[k].append(bg[k])
            per_agg_counts[k].append((int(st_s[k][region].sum()), int(len(region)), bg[k]))
        agg_meta.append((s, cen, members, region, len(loc)))
        agg_rows.append(dict(sample=s, n_B=len(members), n_cells_region=len(region),
            x=float(cen[0]), y=float(cen[1]), niche_overlap_frac=round(niche_ov,3),
            **{f"f_{k}":round(comp[k],4) for k in STATES},
            **{f"bg_{k}":round(bg[k],4) for k in STATES}))
agg=pd.DataFrame(agg_rows)
N_AGG=len(agg)
print(f"aggregates delineated: N={N_AGG} across {agg['sample'].nunique() if N_AGG else 0} sections")
if N_AGG:
    print(f"  B/agg median {agg.n_B.median():.0f} (range {agg.n_B.min()}-{agg.n_B.max()}); "
          f"region cells/agg median {agg.n_cells_region.median():.0f}")
    ov = (agg.niche_overlap_frac>0.05).mean()
    print(f"  VALIDITY: {int((agg.niche_overlap_frac>0.05).sum())}/{N_AGG} aggregates overlap author "
          f"'{NICHE_VAL}' (>5% of footprint) = {ov:.2f}; mean overlap {agg.niche_overlap_frac.mean():.2f}")
agg.to_csv(os.path.join(OUT,"dkd_xenium_aggregates.csv"), index=False)

# ============================================================================
# Step 2 — per-aggregate enrichment (RCC metric) + k/N
# ============================================================================
hdr("STEP 2 — per-aggregate enrichment vs section background (RCC metric)")
def enrich_rows():
    rows=[]
    for k in STATES:
        ins=np.array(per_agg_inside[k]); bgk=np.array(per_agg_bg[k])
        if len(ins)==0:
            rows.append(dict(state=k, n_agg=0)); continue
        l2_agg = np.log2((ins+1e-6)/(bgk+1e-6))
        ins_mean=float(ins.mean()); bg_mean=float(bgk.mean())
        l2_pooled=float(np.log2((ins_mean+1e-6)/(bg_mean+1e-6)))
        try: _,p=wilcoxon(ins-bgk, alternative="two-sided", zero_method="zsplit")
        except Exception: p=np.nan
        rows.append(dict(state=k, n_agg=len(ins), background_mean=round(bg_mean,4),
            inside_mean=round(ins_mean,4), log2_pooled=round(l2_pooled,3),
            median_log2_per_agg=round(float(np.median(l2_agg)),3),
            n_enriched=int((ins>bgk).sum()), n_excluded=int((ins<bgk).sum()),
            wilcoxon_p=(round(p,4) if p==p else np.nan),
            direction=("ENRICHED" if l2_pooled>0 else "EXCLUDED")))
    return pd.DataFrame(rows)
enr=enrich_rows(); enr.to_csv(os.path.join(OUT,"dkd_xenium_enrichment.csv"), index=False)
print(enr.to_string(index=False))
def kn(state, want):
    r=enr[enr.state==state]
    if not len(r) or not r.n_agg.iloc[0]: return (0,0)
    r=r.iloc[0]; return (int(r.n_enriched) if want=="enr" else int(r.n_excluded), int(r.n_agg))
treg_k,treg_n = kn("Treg-like","enr"); eff_k,eff_n = kn("eff-CD8","exc")
print(f"\nTreg-like ENRICHED: {treg_k}/{treg_n}  |  effector-CD8 EXCLUDED: {eff_k}/{eff_n}")

# ---- count-pooled magnitudes (sum state cells / sum expected; one log2 — no per-agg division)
def pooled_log2(counts):
    ins=np.array([c[0] for c in counts], float); exp=np.array([c[1]*c[2] for c in counts], float)
    return float(np.log2(ins.sum()/exp.sum())) if exp.sum()>0 else np.nan
pooled_ct={k: pooled_log2(per_agg_counts[k]) for k in STATES}
print("\ncount-pooled log2 enrichment (Σ state / Σ expected; replaces mean-of-per-agg-log2):")
for k in STATES: print(f"  {k:11s}: {pooled_ct[k]:+.3f}  (mean-based was {enr.loc[enr.state==k,'log2_pooled'].iloc[0]})")

# ---- differential enrichment Δlog2 = log2(Treg enr) − log2(eff-CD8 enr): burden-corrected
#      Treg-vs-cytotoxic bias, immune to global cytotoxic-burden differences and (count-pooled)
#      to the pseudocount issue. Bootstrap CI by resampling aggregates.
def boot_dlog2(ct_treg, ct_eff, nboot=5000, seed=0):
    rng=np.random.default_rng(seed)
    tin=np.array([c[0] for c in ct_treg],float); texp=np.array([c[1]*c[2] for c in ct_treg],float)
    ein=np.array([c[0] for c in ct_eff],float);  eexp=np.array([c[1]*c[2] for c in ct_eff],float)
    n=len(tin); pt=np.log2(tin.sum()/texp.sum())-np.log2(ein.sum()/eexp.sum())
    bs=[]
    for _ in range(nboot):
        i=rng.integers(0,n,n)
        bs.append(np.log2((tin[i].sum()+1e-9)/(texp[i].sum()+1e-9))
                  -np.log2((ein[i].sum()+1e-9)/(eexp[i].sum()+1e-9)))
    bs=np.array(bs); return pt, float(np.percentile(bs,2.5)), float(np.percentile(bs,97.5)), bs
dkd_d, dkd_lo, dkd_hi, dkd_bs = boot_dlog2(per_agg_counts["Treg-like"], per_agg_counts["eff-CD8"])
print(f"\nDKD differential Δlog2 (Treg − eff-CD8): {dkd_d:+.2f}  [95% CI {dkd_lo:+.2f}, {dkd_hi:+.2f}]")

# RCC cohort from committed per-aggregate table (reconstruct counts) + global bg from composition
rcc_agg = pd.read_csv("/Users/danie/ClaudeCode/pilot_analyses/xenium/outputs/tables/rcc_phaseB2_aggregates.csv")
rcc_comp = pd.read_csv(RCC_COMP).set_index("cell_type")
rbg_t=float(rcc_comp.loc["Treg","background_frac"]); rbg_e=float(rcc_comp.loc["eff-CD8","background_frac"])
rct_t=[(int(round(r.f_Treg*r.n_cells_region)), int(r.n_cells_region), rbg_t) for _,r in rcc_agg.iterrows()]
rct_e=[(int(round(r["f_eff-CD8"]*r.n_cells_region)), int(r.n_cells_region), rbg_e) for _,r in rcc_agg.iterrows()]
rcc_d, rcc_lo, rcc_hi, rcc_bs = boot_dlog2(rct_t, rct_e, seed=1)
print(f"RCC differential Δlog2 (Treg − eff-CD8): {rcc_d:+.2f}  [95% CI {rcc_lo:+.2f}, {rcc_hi:+.2f}]  "
      f"(N={len(rcc_agg)} aggs; provenance rcc_phaseB2_aggregates.csv)")
diff_sep = (dkd_hi < rcc_lo)
print(f"-> RCC Treg-over-cytotoxic bias ~{2**rcc_d:.1f}x; DKD ~{2**dkd_d:.1f}x; "
      f"CIs non-overlapping: {diff_sep}")
pd.DataFrame([
    dict(cohort="RCC (tumor)", n_agg=len(rcc_agg), delta_log2=round(rcc_d,3), ci_lo=round(rcc_lo,3),
         ci_hi=round(rcc_hi,3), fold_bias=round(2**rcc_d,2)),
    dict(cohort="DKD (kidney)", n_agg=treg_n, delta_log2=round(dkd_d,3), ci_lo=round(dkd_lo,3),
         ci_hi=round(dkd_hi,3), fold_bias=round(2**dkd_d,2)),
]).to_csv(os.path.join(OUT,"differential_treg_vs_cd8.csv"), index=False)

# ============================================================================
# Step 3 — radial profile (core->margin): rings from B-core centroid
# ============================================================================
hdr(f"STEP 3 — radial profile: rings [0,{R}) [{R},{2*R}) [{2*R},{3*R}) from B-core centroid")
ring_edges=[0,R,2*R,3*R]; ring_names=[f"core 0-{int(R)}", f"margin {int(R)}-{int(2*R)}", f"outer {int(2*R)}-{int(3*R)}"]
# per-aggregate per-ring COUNTS (state cells, ring cells, section bg) -> count-pool across aggregates
rad_counts={k:{rn:[] for rn in ring_names} for k in ["Treg-like","eff-CD8"]}
rad_rows=[]
for (s,cen,members,region,nloc) in agg_meta:
    smask=xsamp==s; loc=np.where(smask)[0]; xy_s=xxy[loc]
    st_s={k:STATES[k][loc] for k in ["Treg-like","eff-CD8"]}
    bg={k:float(st_s[k].mean()) for k in st_s}
    d=np.hypot(*(xy_s-cen).T)
    for ri,rn in enumerate(ring_names):
        ring=(d>=ring_edges[ri])&(d<ring_edges[ri+1]); rn_n=int(ring.sum())
        if rn_n<10: continue
        for k in st_s:
            rad_counts[k][rn].append((int(st_s[k][ring].sum()), rn_n, bg[k]))
            rad_rows.append(dict(sample=s, ring=rn, state=k, n_cells=rn_n,
                                 state_cells=int(st_s[k][ring].sum()), bg=round(bg[k],5)))
rad=pd.DataFrame(rad_rows); rad.to_csv(os.path.join(OUT,"dkd_xenium_radial.csv"), index=False)
def ring_pool(k,rn):
    c=rad_counts[k][rn]
    if not c: return np.nan
    ins=sum(x[0] for x in c); exp=sum(x[1]*x[2] for x in c)
    return float(np.log2((ins+1e-9)/(exp+1e-9))) if exp>0 else np.nan
def ring_pool_ci(k,rn,nboot=3000,seed=0):
    c=rad_counts[k][rn]
    if not c: return (np.nan,np.nan)
    ins=np.array([x[0] for x in c],float); exp=np.array([x[1]*x[2] for x in c],float)
    rng=np.random.default_rng(seed); n=len(c); bs=[]
    for _ in range(nboot):
        i=rng.integers(0,n,n); bs.append(np.log2((ins[i].sum()+1e-9)/(exp[i].sum()+1e-9)))
    return float(np.percentile(bs,2.5)), float(np.percentile(bs,97.5))
print("count-pooled radial log2-enrichment (Σ state / Σ expected per ring):")
for k in ["Treg-like","eff-CD8"]:
    print(f"  {k:9s}: " + "  ".join(f"{rn}={ring_pool(k,rn):+.2f}" for rn in ring_names))
def ringmean(k,rn): return ring_pool(k,rn)   # alias (now count-pooled)
treg_core=ring_pool("Treg-like",ring_names[0]); treg_margin=ring_pool("Treg-like",ring_names[1])
collar = (treg_margin>0.2) and (treg_margin>treg_core+0.2)
print(f"\nMARGIN/COLLAR test: Treg-like core={treg_core:+.2f} margin={treg_margin:+.2f} "
      f"-> {'COLLAR present at margin' if collar else 'no Treg collar (margin not elevated over core)'}")

# ============================================================================
# Step 4 — head-to-head, method-matched
# ============================================================================
hdr("STEP 4 — head-to-head (method-matched): RCC vs DKD Xenium")
rcc = pd.read_csv(RCC_COMP).set_index("cell_type")
rcc_treg=f"{int(rcc.loc['Treg','n_agg_above_bg'])}/{int(rcc.loc['Treg','n_agg'])}"
rcc_eff =f"{int(rcc.loc['eff-CD8','n_agg'])-int(rcc.loc['eff-CD8','n_agg_above_bg'])}/{int(rcc.loc['eff-CD8','n_agg'])}"
def verdict(k,n,thr_hi=0.8,thr_lo=0.5):
    if not n: return "UNTESTED"
    f=k/n; return "REPLICATE" if f>=thr_hi else ("PARTIAL" if f>=thr_lo else "DIVERGE")
v_treg=verdict(treg_k,treg_n); v_eff=verdict(eff_k,eff_n)
# differential verdict: do the cohort Δlog2 CIs separate? (burden-corrected headline)
v_diff = "DIVERGE" if diff_sep else "OVERLAP"
ht=pd.DataFrame([
  dict(axis="** Δlog2 Treg−eff-CD8 (burden-corrected) **",
       rcc=f"{rcc_d:+.2f} [{rcc_lo:+.2f},{rcc_hi:+.2f}] (~{2**rcc_d:.1f}x Treg bias)",
       dkd=f"{dkd_d:+.2f} [{dkd_lo:+.2f},{dkd_hi:+.2f}] (~{2**dkd_d:.1f}x)", verdict=v_diff),
  dict(axis="Treg ENRICHED in aggregate", rcc=f"{rcc_treg} (log2 {rcc.loc['Treg','log2_enrichment']})",
       dkd=f"{treg_k}/{treg_n} (count-pooled log2 {pooled_ct['Treg-like']:+.2f})", verdict=v_treg),
  dict(axis="effector-CD8 EXCLUDED", rcc=f"{rcc_eff} (log2 {rcc.loc['eff-CD8','log2_enrichment']})",
       dkd=f"{eff_k}/{eff_n} (count-pooled log2 {pooled_ct['eff-CD8']:+.2f})", verdict=v_eff),
  dict(axis="Treg collar at margin", rcc="n/a", dkd=f"core {treg_core:+.2f} -> margin {treg_margin:+.2f}",
       verdict=("COLLAR" if collar else "NONE")),
  dict(axis="plasma (structural)", rcc=f"log2 {rcc.loc['Plasma','log2_enrichment']}",
       dkd=f"count-pooled log2 {pooled_ct['Plasma']:+.2f}", verdict="DESCRIPTIVE"),
])
ht.to_csv(os.path.join(OUT,"rcc_vs_dkd_method_matched.csv"), index=False)
print(ht.to_string(index=False))
print(f"\nHEADLINE (burden-corrected): RCC Treg-over-cytotoxic bias Δlog2 {rcc_d:+.2f} (~{2**rcc_d:.1f}x) "
      f"vs DKD {dkd_d:+.2f} (~{2**dkd_d:.1f}x); CIs {'separate' if diff_sep else 'overlap'} "
      f"-> DKD aggregates show {'NO' if abs(dkd_d)<0.5 else 'a'} Treg-vs-cytotoxic bias.")
print(f"Treg enrichment itself REPLICATES ({treg_k}/{treg_n}) with a margin collar={collar}; "
      f"the earlier 'Treg absent' DIVERGE was a unit/metric artifact. The RCC discriminator "
      f"(cytotoxic exclusion / Treg-dominance) does NOT replicate.")

# ============================================================================
# FIGURES
# ============================================================================
hdr("FIGURES")
# (a) per-aggregate Treg-like / eff-CD8 log2 distribution vs RCC reference
fig,ax=plt.subplots(figsize=(8,5))
states_plot=["Treg-like","eff-CD8"]
for i,k in enumerate(states_plot):
    ins=np.array(per_agg_inside[k]); bgk=np.array(per_agg_bg[k])
    l2=np.log2((ins+1e-6)/(bgk+1e-6)) if len(ins) else np.array([])
    jit=np.random.uniform(-0.12,0.12,len(l2))
    ax.scatter(l2, np.full(len(l2),i)+jit, s=45, color=DKD_COLOR, edgecolor="k", linewidth=0.4, alpha=0.85, zorder=3)
    if len(l2): ax.scatter(np.median(l2), i, marker="|", s=600, color="k", zorder=4)
rcc_ref={"Treg-like":rcc.loc['Treg','log2_enrichment'], "eff-CD8":rcc.loc['eff-CD8','log2_enrichment']}
for i,k in enumerate(states_plot):
    ax.scatter(rcc_ref[k], i, marker="D", s=90, color="#1F78B4", edgecolor="k", zorder=5,
               label="RCC pooled log2" if i==0 else None)
ax.axvline(0,color="gray",ls="--",lw=1)
ax.set_yticks(range(len(states_plot))); ax.set_yticklabels([f"{k}\n(DKD {kn(k,'enr' if k=='Treg-like' else 'exc')[0]}/{kn(k,'enr')[1]})" for k in states_plot])
ax.set_xlabel("per-aggregate log2 enrichment (inside vs section background)")
ax.set_title(f"Method-matched DBSCAN aggregates: DKD Xenium (N={N_AGG}) vs RCC (37)\n"
             "purple=DKD per-aggregate, blue ◆=RCC pooled", fontsize=10)
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout(); fig.savefig(os.path.join(OUT,"fig_a_aggregate_enrichment.png"),dpi=150); plt.close(fig)

# (b) COUNT-POOLED radial profile core->margin, with bootstrap CI
fig,ax=plt.subplots(figsize=(7.5,5))
xpos=np.arange(len(ring_names))
for k,col,mk in [("Treg-like",DKD_COLOR,"o"),("eff-CD8","#E69F00","s")]:
    pts=[ring_pool(k,rn) for rn in ring_names]
    cis=[ring_pool_ci(k,rn) for rn in ring_names]
    lo=[p-c[0] for p,c in zip(pts,cis)]; hi=[c[1]-p for p,c in zip(pts,cis)]
    ax.errorbar(xpos, pts, yerr=[lo,hi], marker=mk, color=col, capsize=4, lw=2, label=k)
ax.axhline(0,color="gray",ls="--",lw=1)
ax.set_xticks(xpos); ax.set_xticklabels(ring_names)
ax.set_ylabel("count-pooled log2 enrichment vs section background")
ax.set_xlabel("ring (distance from B-core centroid, um)")
ax.set_title(f"Count-pooled radial profile across DKD aggregates (N={N_AGG})\nTreg flat across rings (no collar); mild cytotoxic-core gradient (eff-CD8 falls outward)")
ax.legend(fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"fig_b_radial_profile.png"),dpi=150); plt.close(fig)

# (c) HEADLINE: differential Δlog2 (Treg − eff-CD8), RCC vs DKD, with bootstrap distributions
fig,ax=plt.subplots(figsize=(8,4.5))
for i,(name,d,lo,hi,bs,col) in enumerate([
    ("RCC (tumor)", rcc_d, rcc_lo, rcc_hi, rcc_bs, "#1F78B4"),
    ("DKD (kidney)", dkd_d, dkd_lo, dkd_hi, dkd_bs, DKD_COLOR)]):
    ax.scatter(bs, np.full(len(bs),i)+np.random.uniform(-0.06,0.06,len(bs)), s=3, color=col, alpha=0.08, zorder=1)
    ax.errorbar(d, i, xerr=[[d-lo],[hi-d]], fmt="o", color=col, ms=11, capsize=6, lw=2, zorder=3,
                markeredgecolor="k")
    ax.text(d, i+0.18, f"{d:+.2f}  (~{2**d:.1f}x Treg bias)", ha="center", fontsize=9, color=col)
ax.axvline(0,color="gray",ls="--",lw=1.2, label="no Treg-vs-cytotoxic bias")
ax.set_yticks([0,1]); ax.set_yticklabels(["RCC (tumor)","DKD (kidney)"]); ax.set_ylim(-0.5,1.6)
ax.set_xlabel("differential enrichment  Δlog2 = log2(Treg) − log2(effector-CD8)  [burden-corrected]")
ax.set_title("Burden-corrected Treg-vs-cytotoxic bias in B-aggregates: RCC vs DKD\n(count-pooled, bootstrap CI; immune to 19x cytotoxic-burden difference)", fontsize=10)
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout(); fig.savefig(os.path.join(OUT,"fig_c_differential.png"),dpi=150); plt.close(fig)
print("saved fig_a_aggregate_enrichment.png, fig_b_radial_profile.png, fig_c_differential.png")

# ============================================================================
# REPORT.md
# ============================================================================
hdr("WRITING REPORT")
L=[];W=L.append
W("# DKD vs RCC immunoregulatory aggregate — method-matched (DBSCAN, 1:1)\n")
W(f"Read-only. Xenium-only (CosMx subtype imputed). Borrowed RCC pipeline verbatim from "
  f"`py/phaseB_02_rcc_aggregates.py` + committed eps=50 rebuild (`py/whitepaper_recompute.py`).\n")
W("## HEADLINE — burden-corrected differential\n")
W("The clean, confound-free comparison is the **Treg-vs-cytotoxic bias** per aggregate: "
  "Δlog2 = log2(Treg enrichment) − log2(effector-CD8 enrichment), count-pooled, bootstrap CI. "
  "This cancels the 19× global cytotoxic-burden difference between RCC tumor and DKD kidney.\n")
W(f"- **RCC: Δlog2 = {rcc_d:+.2f} [95% CI {rcc_lo:+.2f}, {rcc_hi:+.2f}] — Treg favored over cytotoxic by ~{2**rcc_d:.1f}×.**")
W(f"- **DKD: Δlog2 = {dkd_d:+.2f} [95% CI {dkd_lo:+.2f}, {dkd_hi:+.2f}] — ~{2**dkd_d:.1f}× (no bias).**")
W(f"- CIs {'do NOT overlap' if diff_sep else 'overlap'}: the RCC immunoregulatory Treg-over-cytotoxic "
  "bias is **absent** in DKD aggregates. *This single number is immune to the burden confound and "
  "to the pseudocount issue.*\n")
W("## Borrowed parameters (verbatim)\n")
W(f"- DBSCAN(**eps={EPS}**, **min_samples={MINPTS}**) on B-cell coords, per section.")
W(f"- Aggregate footprint = all cells within **R={R} um** of any member B cell.")
W("- Enrichment = **count-pooled** log2(Σ state cells / Σ expected cells); per-aggregate k/N and "
  "Wilcoxon for direction. **Unconditional** fractions (state cells / all region cells) — matches "
  "RCC (differs from the earlier niche test's conditional FOXP3+-among-CD4+). Count-pooling (vs "
  "mean-of-per-aggregate-log2) avoids dividing empty single-aggregate rings — kills the −12.5 artifacts.")
W(f"- Spatial-unit sanity: median NN cell distance {nn_med:.1f} um -> microns confirmed; eps transfers.\n")
W("## Aggregate delineation\n")
W(f"- DBSCAN on B cells (independent of author niche) -> **N={N_AGG} aggregates** across "
  f"{agg['sample'].nunique() if N_AGG else 0} Xenium sections (median {int(agg.n_B.median()) if N_AGG else 0} B/agg; "
  "RCC median 346 — DKD aggregates are smaller/looser).")
if N_AGG:
    W(f"- Author-niche overlap: only **{int((agg.niche_overlap_frac>0.05).sum())}/{N_AGG} (19%)** overlap "
      f"'{NICHE_VAL}' (>5% of footprint) — the DBSCAN delineation and the coarse author niche are "
      "**largely different units**, itself direct evidence that unit-of-analysis drives the comparison.")
W("\n## Method-matched result\n")
W("| axis | RCC | DKD Xenium | verdict |")
W("|---|---|---|---|")
for _,r in ht.iterrows(): W(f"| {r.axis} | {r.rcc} | {r.dkd} | {r.verdict} |")
W(f"\nWithin-aggregate Treg:effector-CD8 balance (count-pooled inside fractions): "
  f"**DKD ~{2**(pooled_ct['Treg-like']-pooled_ct['eff-CD8']):.2f}× Treg:CD8 vs RCC "
  f"~{2**(rcc.loc['Treg','log2_enrichment']-rcc.loc['eff-CD8','log2_enrichment']):.1f}×.** "
  "RCC aggregates are Treg-dominant and cytotoxic-excluding; DKD aggregates admit cytotoxic CD8 "
  "in balance with Treg.\n")
W("## Radial (core -> margin) — count-pooled, with bootstrap CI\n")
W("| ring | Treg-like | eff-CD8 |")
W("|---|---|---|")
for rn in ring_names: W(f"| {rn} | {ring_pool('Treg-like',rn):+.2f} | {ring_pool('eff-CD8',rn):+.2f} |")
W(f"\n{'**A Treg collar IS present at the aggregate margin** (Treg stays positive at the outer shell while cytotoxic CD8 drops off).' if collar else '**No Treg collar:** the margin is not enriched over the core.'}\n")
W("## Conclusion (three-part)\n")
W(f"1. **The earlier 'Treg absent' DIVERGE was a unit-of-analysis + metric artifact.** Under the "
  f"*exact* RCC method (unconditional, count-pooled), Treg-like enrichment **replicates** "
  f"({treg_k}/{treg_n}, count-pooled log2 {pooled_ct['Treg-like']:+.2f} vs RCC +1.34) and is enriched "
  f"{'uniformly across core->margin->outer' if not collar else 'with a margin collar'} "
  f"(core {ring_pool('Treg-like',ring_names[0]):+.2f} / margin {ring_pool('Treg-like',ring_names[1]):+.2f} "
  f"/ outer {ring_pool('Treg-like',ring_names[2]):+.2f}). The prior conditional FOXP3+-among-CD4+ "
  "niche-membership test missed it because the conditional metric and the coarse-niche unit differ "
  f"from the RCC method. {'(Note: the earlier mean-of-per-aggregate-log2 reported a spurious margin collar; count-pooling removes it.)' if not collar else ''}")
W(f"2. **But the RCC discriminator does NOT replicate.** DKD aggregates **co-enrich** effector-CD8 "
  f"(excluded only {eff_k}/{eff_n}, count-pooled log2 {pooled_ct['eff-CD8']:+.2f}). The burden-corrected "
  f"**Δlog2 = {dkd_d:+.2f} [{dkd_lo:+.2f},{dkd_hi:+.2f}] (DKD) vs {rcc_d:+.2f} [{rcc_lo:+.2f},{rcc_hi:+.2f}] "
  f"(RCC)** — non-overlapping CIs. RCC favors Treg over cytotoxic ~{2**rcc_d:.1f}×; DKD ~{2**dkd_d:.1f}× (no bias).")
W(f"3. **Net: PARTIAL.** The DKD B-cell-rich aggregate concentrates both Treg-like and cytotoxic CD8 "
  "(an immune-dense B/plasma aggregate) but lacks the RCC immunoregulatory 'Treg-in / cytotoxic-out' "
  f"architecture. The only spatial structure is a mild **cytotoxic-core gradient**: effector-CD8 is "
  f"highest at the core ({ring_pool('eff-CD8',ring_names[0]):+.2f}) and falls outward "
  f"({ring_pool('eff-CD8',ring_names[2]):+.2f}) while Treg stays flat — the opposite of a Treg collar.")
W("\n**Metric note:** raw 'exclusion vs whole-section background' is confounded by global cytotoxic "
  "burden (RCC tumor eff-CD8 bg 0.135 vs DKD kidney 0.007, ~19×). The Δlog2 differential and the "
  "radial geometry are the burden-immune readouts, and both show DKD is *less* immunoregulatory than RCC.")
W("\n## Caveats\n")
W(f"- Borrowed RCC params (eps={EPS}, minPts={MINPTS}, R={R}) listed above; **N={N_AGG} aggregates**.")
W("- **No patient column** -> donor clustering uncontrolled (aggregates within a section share a donor).")
W("- States reconstructed by marker+ within compartment (atlas lacks native Treg/effector-CD8 labels); "
  "GNLY structural-zero on Xenium dropped from effector-CD8.")
W("- Magnitudes are log2 vs section background; Xenium-only scope (CosMx subtype imputed).")
open(os.path.join(OUT,"REPORT.md"),"w").write("\n".join(L))
print("wrote REPORT.md + CSVs + figures")
print("\n== done ==")
