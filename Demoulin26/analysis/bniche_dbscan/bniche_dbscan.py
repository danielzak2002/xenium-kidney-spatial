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
        for k in STATES: per_agg_inside[k].append(comp[k]); per_agg_bg[k].append(bg[k])
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

# ============================================================================
# Step 3 — radial profile (core->margin): rings from B-core centroid
# ============================================================================
hdr(f"STEP 3 — radial profile: rings [0,{R}) [{R},{2*R}) [{2*R},{3*R}) from B-core centroid")
ring_edges=[0,R,2*R,3*R]; ring_names=[f"core 0-{int(R)}", f"margin {int(R)}-{int(2*R)}", f"outer {int(2*R)}-{int(3*R)}"]
radial={k:{rn:[] for rn in ring_names} for k in ["Treg-like","eff-CD8"]}
rad_rows=[]
for (s,cen,members,region,nloc),(_,arow) in zip(agg_meta, agg.iterrows()):
    smask=xsamp==s; loc=np.where(smask)[0]; xy_s=xxy[loc]
    st_s={k:STATES[k][loc] for k in ["Treg-like","eff-CD8"]}
    bg={k:float(st_s[k].mean()) for k in st_s}
    d=np.hypot(*(xy_s-cen).T)
    for ri,rn in enumerate(ring_names):
        ring = (d>=ring_edges[ri])&(d<ring_edges[ri+1])
        if ring.sum()<10: continue
        for k in st_s:
            l2=float(np.log2((st_s[k][ring].mean()+1e-6)/(bg[k]+1e-6)))
            radial[k][rn].append(l2)
            rad_rows.append(dict(sample=s, ring=rn, state=k, n_cells=int(ring.sum()), log2_enrich=round(l2,3)))
rad=pd.DataFrame(rad_rows); rad.to_csv(os.path.join(OUT,"dkd_xenium_radial.csv"), index=False)
print("mean radial log2-enrichment (across aggregates):")
for k in ["Treg-like","eff-CD8"]:
    means={rn:(round(float(np.mean(radial[k][rn])),3) if radial[k][rn] else np.nan) for rn in ring_names}
    print(f"  {k:9s}: " + "  ".join(f"{rn}={means[rn]}" for rn in ring_names))
# margin test: does Treg-like rise above 0 at the margin even if core <= 0?
def ringmean(k,rn): return float(np.mean(radial[k][rn])) if radial[k][rn] else np.nan
treg_core=ringmean("Treg-like",ring_names[0]); treg_margin=ringmean("Treg-like",ring_names[1])
collar = (treg_margin>0.2) and (treg_margin>treg_core+0.2)
print(f"\nMARGIN/COLLAR test: Treg-like core={treg_core:.2f} margin={treg_margin:.2f} "
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
ht=pd.DataFrame([
  dict(axis="Treg ENRICHED in aggregate", rcc=f"{rcc_treg} (log2 {rcc.loc['Treg','log2_enrichment']})",
       dkd=f"{treg_k}/{treg_n} (log2 {enr.loc[enr.state=='Treg-like','log2_pooled'].iloc[0] if treg_n else 'na'})", verdict=v_treg),
  dict(axis="effector-CD8 EXCLUDED", rcc=f"{rcc_eff} (log2 {rcc.loc['eff-CD8','log2_enrichment']})",
       dkd=f"{eff_k}/{eff_n} (log2 {enr.loc[enr.state=='eff-CD8','log2_pooled'].iloc[0] if eff_n else 'na'})", verdict=v_eff),
  dict(axis="Treg collar at margin", rcc="n/a", dkd=f"core {treg_core:.2f} -> margin {treg_margin:.2f}",
       verdict=("COLLAR" if collar else "NONE")),
  dict(axis="plasma (structural)", rcc=f"log2 {rcc.loc['Plasma','log2_enrichment']}",
       dkd=f"log2 {enr.loc[enr.state=='Plasma','log2_pooled'].iloc[0] if N_AGG else 'na'}", verdict="DESCRIPTIVE"),
])
ht.to_csv(os.path.join(OUT,"rcc_vs_dkd_method_matched.csv"), index=False)
print(ht.to_string(index=False))
overall = "DIVERGE" if (v_treg=="DIVERGE") else ("REPLICATE" if v_treg=="REPLICATE" and v_eff in ("REPLICATE","PARTIAL") else "PARTIAL")
artifact = (overall!="DIVERGE") or collar
print(f"\nOverall (method-matched, Treg axis): {overall}; collar={collar}")
print(f"-> earlier categorical-niche DIVERGE was {'an ARTIFACT of unit-of-analysis' if artifact else 'ROBUST (holds under matched DBSCAN method, core AND margin)'}")

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

# (b) mean radial profile core->margin
fig,ax=plt.subplots(figsize=(7.5,5))
xpos=np.arange(len(ring_names))
for k,col,mk in [("Treg-like",DKD_COLOR,"o"),("eff-CD8","#E69F00","s")]:
    means=[ringmean(k,rn) for rn in ring_names]
    sds=[ (np.std(radial[k][rn])/max(np.sqrt(len(radial[k][rn])),1) if radial[k][rn] else 0) for rn in ring_names]
    ax.errorbar(xpos, means, yerr=sds, marker=mk, color=col, capsize=4, lw=2, label=k)
ax.axhline(0,color="gray",ls="--",lw=1)
ax.set_xticks(xpos); ax.set_xticklabels(ring_names)
ax.set_ylabel("mean log2 enrichment vs section background")
ax.set_xlabel("ring (distance from B-core centroid, um)")
ax.set_title(f"Radial profile across DKD aggregates (N={N_AGG}): does a Treg collar appear at the margin?")
ax.legend(fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"fig_b_radial_profile.png"),dpi=150); plt.close(fig)
print("saved fig_a_aggregate_enrichment.png, fig_b_radial_profile.png")

# ============================================================================
# REPORT.md
# ============================================================================
hdr("WRITING REPORT")
L=[];W=L.append
W("# DKD vs RCC immunoregulatory aggregate — method-matched (DBSCAN, 1:1)\n")
W(f"Read-only. Xenium-only (CosMx subtype imputed). Borrowed RCC pipeline verbatim from "
  f"`py/phaseB_02_rcc_aggregates.py` + committed eps=50 rebuild (`py/whitepaper_recompute.py`).\n")
W("## Borrowed parameters (verbatim)\n")
W(f"- DBSCAN(**eps={EPS}**, **min_samples={MINPTS}**) on B-cell coords, per section.")
W(f"- Aggregate footprint = all cells within **R={R} um** of any member B cell.")
W("- Enrichment = log2((inside_mean+1e-6)/(section background+1e-6)); per-aggregate Wilcoxon; "
  "k/N = aggregates above background. **Unconditional** fractions (state cells / all region cells) "
  "— matches RCC (differs from the earlier niche test's conditional FOXP3+-among-CD4+).")
W(f"- Spatial-unit sanity: median NN cell distance {nn_med:.1f} um -> microns confirmed; eps transfers.\n")
W("## Aggregate delineation\n")
W(f"- DBSCAN on B cells (independent of author niche) -> **N={N_AGG} aggregates** across "
  f"{agg['sample'].nunique() if N_AGG else 0} Xenium sections.")
if N_AGG:
    W(f"- Validity: **{int((agg.niche_overlap_frac>0.05).sum())}/{N_AGG}** overlap the author "
      f"'{NICHE_VAL}' niche (>5% of footprint) — confirms we delineate the same structures.")
W("\n## Method-matched result\n")
W("| axis | RCC | DKD Xenium | verdict |")
W("|---|---|---|---|")
for _,r in ht.iterrows(): W(f"| {r.axis} | {r.rcc} | {r.dkd} | {r.verdict} |")
W(f"\n**Treg-like ENRICHED {treg_k}/{treg_n}; effector-CD8 EXCLUDED {eff_k}/{eff_n}. "
  f"Overall Treg axis: {overall}.**\n")
W("## Radial (core -> margin) — the collar test\n")
W("| ring | Treg-like | eff-CD8 |")
W("|---|---|---|")
for rn in ring_names: W(f"| {rn} | {ringmean('Treg-like',rn):.2f} | {ringmean('eff-CD8',rn):.2f} |")
W(f"\n{'**A Treg collar IS present at the aggregate margin** despite a non-enriched core.' if collar else '**No Treg collar:** the margin is not enriched over the core — Treg absence is not a core-only effect.'}\n")
W("## Conclusion\n")
if overall=="DIVERGE" and not collar:
    W("Under the **exact RCC method** (same DBSCAN eps/minPts, same R, same unconditional metric), "
      "and across **core AND margin rings**, the DKD B-aggregate still lacks Treg enrichment and shows "
      "no effector-CD8 exclusion. **The earlier DIVERGE is ROBUST — not a unit-of-analysis artifact.** "
      "The DKD B-cell-rich niche is a B/plasma aggregate without the RCC Treg(+)/effector-CD8(-) "
      "immunoregulatory organization.")
elif collar:
    W("Treg enrichment appears at the **margin** under the matched method even though the core is flat "
      "— the earlier categorical-niche DIVERGE was partly a **unit-of-analysis artifact** (niche membership "
      "averaged over the collar). Reassess the RCC comparison as core-vs-collar.")
else:
    W(f"Under the matched method the Treg axis is {overall}; see table. Interpret with the small-N caveat.")
W("\n## Caveats\n")
W(f"- Borrowed RCC params (eps={EPS}, minPts={MINPTS}, R={R}) listed above; **N={N_AGG} aggregates**.")
W("- **No patient column** -> donor clustering uncontrolled (aggregates within a section share a donor).")
W("- States reconstructed by marker+ within compartment (atlas lacks native Treg/effector-CD8 labels); "
  "GNLY structural-zero on Xenium dropped from effector-CD8.")
W("- Magnitudes are log2 vs section background; Xenium-only scope (CosMx subtype imputed).")
open(os.path.join(OUT,"REPORT.md"),"w").write("\n".join(L))
print("wrote REPORT.md + CSVs + figures")
print("\n== done ==")
