# Peri-aggregate stromal BAFF — ambient/spillover stress-test

The disease analysis reported peri-aggregate **stromal BAFF enriched ~2.1–2.5×** near B-lineage
aggregates and floated a TLS-survival-niche reading. Before relying on it, this stress-tests that
signal against the obvious confound: **B-aggregates are dense regions, and a local ambient/spillover
field would raise BAFF in *any* nearby cell regardless of who produces it.** Reuses the
`dkd_xenium_disease` object (validated reannotation labels + per-cell BAFF/TNFSF13B + nCount).
Read-only. Scripts: `baff_ambient_control.py`; tables `baff_ambient_step{1,2,3}.csv`; figure
`figures/baff_ambient_control.png`.

> No Xenium negative-control-probe / blank-codeword features survive in this deposited object (only
> the real gene *NEGR1*). So the **local non-producer ambient floor = epithelial cells in the same
> window**, and **per-cell total counts (nCount) = the transcript-density proxy**. Aggregate-bearing
> sections tested: the two B-rich DKD (1006, HK2695) + the two non-DKD with aggregates (C3GN 1007,
> AA amyloid 1009). BAFF is rare (0.3–2% detection) so counts are small — pooled where needed.

## STEP 1 — local-ambient control (stromal vs non-producer epithelium, same window)

| zone | cell type | BAFF det % | n | positive |
|---|---|--:|--:|--:|
| near (≤50 µm) | Stromal | **1.96** | 1,326 | 26 |
| near (≤50 µm) | Epithelial (non-producer) | 0.62 | 812 | 5 |
| far (>200 µm) | Stromal | 0.80 | 48,369 | 389 |
| far (>200 µm) | Epithelial (non-producer) | 0.29 | 116,772 | 335 |

Stromal-near (1.96%) exceeds the local epithelial ambient (0.62%) by **3.2×** — which *looks* like
production. **But** the same stroma>epithelium ratio holds in the FAR zone (0.80% vs 0.29% = 2.8×):
this is a **tissue-wide producer-baseline difference (stroma makes more BAFF than epithelium
everywhere), not a peri-aggregate effect.**

## STEP 2 — near vs far within stromal (per section)

| section | cond | stromal near % | stromal far % | near/far |
|---|---|--:|--:|--:|
| 1006 | DKD | 1.88 | 1.19 | 1.58× |
| HK2695 | DKD | 1.75 | 0.73 | 2.40× |
| 1007 | C3GN | 2.47 | 0.67 | 3.69× |
| 1009 | AA amyloid | 2.03 | 0.74 | 2.75× |
| **POOLED** | – | **1.96** | **0.80** | **2.44×** |

Stromal BAFF *does* rise near aggregates (~2.4× pooled, consistent across all four sections). **In
isolation this looks like localization** — which is why STEP 3 is necessary.

## STEP 3 — density confound / specificity (decisive)

Near/far BAFF ratio per cell type, with the transcript-density (nCount) near/far ratio as reference:

| cell type | BAFF near % | BAFF far % | **BAFF near/far** | nCount near/far |
|---|--:|--:|--:|--:|
| Stromal | 1.96 | 0.80 | **2.44×** | 0.98 |
| Epithelial (non-producer) | 0.62 | 0.29 | **2.15×** | 0.87 |
| Endothelial | 1.06 | 0.37 | **2.86×** | 0.98 |
| Myeloid | 3.99 | 3.50 | 1.14× | 0.96 |
| T | 1.18 | 1.19 | 1.00× | 1.09 |
| B-lineage | 1.27 | 1.39 | 0.91× | 1.12 |

**The near-aggregate BAFF rise is NOT stromal-specific.** Non-producer **epithelium rises 2.15×** and
**endothelium 2.86×** — essentially the same as **stromal 2.44×**. Transcript density (nCount) is
**flat (~0.98×)**, so it is not a trivial "more counts per cell" artefact; rather, BAFF specifically
is elevated in the local *ambient pool* around aggregates and is picked up equally by every bystander
parenchymal/stromal cell. The genuine tissue-wide producer, **myeloid, is flat (1.14×)** — it does
not concentrate at aggregates; T and B-lineage are flat too.

## VERDICT — ambient/spillover, NOT localized stromal production

The peri-aggregate stromal BAFF signal **does not survive the controls.** It passes the naive STEP 1/2
look (stromal > local epithelium; stromal near > far) but **fails STEP 3**: the elevation is shared
by non-producer epithelium and endothelium at the same magnitude, with flat transcript density — the
signature of a **local ambient/spillover field** around the (dense, BAFF-bearing) aggregate, not
stromal-specific production.

**What to drop:** the *localized stromal BAFF production* claim and the "stromal TLS survival-niche"
reading built on it.

**What survives (unchanged by this test):**
- **Myeloid BAFF is the robust producer signal** — tissue-wide (near/far 1.14×; 3.5–4.0% detection,
  7.9× the epithelial floor in the disease analysis). Not aggregate-localized, but real and dominant.
- **Stroma > epithelium BAFF baseline** is real but **tissue-wide**, not a niche.
- **B-aggregates do sit in a locally BAFF-elevated field** — but that field is ambient/spillover, so
  it cannot be attributed to peri-aggregate stromal production.

**Caveat / open item:** the receptor-side observation (BAFF-R on B cells elevated ~1.6× in aggregates)
was **not** re-tested here and could be subject to the *same* spillover confound; it should get an
identical non-producer/near-far control before any "B cells receive BAFF in the niche" claim is made.

**Net:** keep myeloid tissue-wide BAFF; **retract the localized peri-aggregate stromal-BAFF niche
claim.** Honest negative result — the conditioned signal was real at the producer-baseline level but
the spatial-localization interpretation was ambient/spillover.
