# Methods

Draft methods section for the AD-vs-PD whole-genome association study across AMP-AD and AMP-PD.
Journal style: what was done and why, in the order it was done. Operational instructions live in
`README.md`; project state in `HANDOFF.md`; the decision history, including hypotheses that were
tested and rejected, in `PROJECT_LOG.md`.

---

## 1. Cohorts and genotype callsets

Whole-genome sequence data were assembled from four independently generated callsets spanning five
clinical cohorts across two consortia.

| Callset | Consortium | Cohorts | Genomes as supplied | Build | Source format |
|---|---|---|---|---|---|
| `wgs_harm` | AMP-AD | ROSMAP, MayoRNAseq, MSBB | 1,894 | GRCh37 → GRCh38 | 24 GATK scatter-interval VCFs |
| `divco_hs` | AMP-AD | Diverse Cohorts | 1,019 | GRCh38 | single merged VCF |
| `wb_dwgs` | AMP-PD | AMP-PD (whole blood) | 10,418 | GRCh38 | PLINK2 pfiles |
| `br_dsnwgs` | AMP-PD | AMP-PD (postmortem brain) | 97 | GRCh38 | single joint-call VCF |

A callset is not a cohort. `wgs_harm` is a joint call over three AMP-AD studies, and a donor
enrolled in both ROSMAP and Diverse Cohorts appears in two callsets as two separately sequenced
genomes. Each genome therefore carries both `source_callset` (where it was sequenced and
joint-called) and `source_dataset` (which study phenotyped the donor); deduplicating on donor
identity would discard real data, and the two fields are kept distinct throughout.

`wgs_harm` was lifted from GRCh37 to GRCh38 with UCSC `liftOver`; the other three were natively
called on GRCh38. This asymmetry is the origin of the technical artifact addressed in §6.

`br_dsnwgs` is a 97-donor joint call and therefore emits no records at sites monomorphic within its
own donors. It contributes to approximately half the final association set. This is expected
behaviour of joint calling at small *n*, not a defect, but it has consequences reported in §10.

## 2. Phenotype definition

**AD was defined neuropathologically throughout.** Clinical dementia diagnosis codes were never
used to define an AD case. Criteria follow the AMP-AD diagnosis-criteria specification (Synapse
`syn51757663`; harmonized data dictionary `syn73713784`).

Two labels were derived per donor. The primary label `pheno` ∈ {AD, PD, control, other}; a secondary
label `dx_detailed` ∈ {PD, AD, MCI, DLB, PSP, control, other} for secondary contrasts. Where `pheno`
is AD, `dx_detailed` is AD by construction, so the AD set is identical under both labels and a donor
cannot be AD in one analysis and MCI in another.

| Cohort | Rule for `pheno` |
|---|---|
| ROSMAP, MSBB | Braak stage + CERAD score. AD = Braak ≥ IV and CERAD moderate/frequent; control = Braak ≤ III and CERAD sparse/none; otherwise other |
| MayoRNAseq | Braak stage + Thal phase. AD = Braak ≥ IV and Thal ≥ 2; control = Braak ≤ III and Thal < 2; otherwise other — then the control-purity screen below |
| Diverse Cohorts | pre-adjudicated `ADoutcome` (`mayoDx` for the Mayo contribution group) |
| AMP-PD | curated `case_control_other_latest`: Case → PD, Control → control, Other excluded |

Missing staging inputs yield a null label, not `other`: `other` means classified as neither case nor
control, whereas null means not classifiable. The distinction matters because nulls are excluded from
every arm rather than pooled into a heterogeneous comparison group.

**Control-purity screen (MayoRNAseq only).** Braak and Thal are both AD-specific axes, so a
progressive supranuclear palsy brain — a 4R tauopathy with little amyloid — scores as a control on
both and would contaminate the control arm. Any neuropathological control whose independent
`diagnosis` field named a disease was demoted to `other`. The screen is one-directional: it never
rescues a null and never overrides an AD call. It was applied to Mayo alone because Mayo is the only
cohort carrying an independent neuropathological diagnosis alongside the staging axes.

**Diverse Cohorts retained its pre-adjudicated call** rather than being recomputed from Braak and
Thal, because the harmonized Thal column is a lossy subset of the plaque data the original
adjudication rested on; recomputing would have nulled donors who have Braak but no Thal.

Braak stage alone was never used for any cohort, as it misclassifies amyloid-negative age-related
tauopathy as AD.

For `dx_detailed`, AMP-PD used `diagnosis_latest` (→ PD / DLB / PSP / control, else other); ROSMAP
used the clinical consensus code `dcfdx_lv` (1 → control, 2–3 → MCI, 4–6 → other); Mayo mapped
`diagnosis` = progressive supranuclear palsy → PSP and otherwise retained `pheno`; MSBB and Diverse
Cohorts retained `pheno`, having no independent specific diagnosis. ROSMAP codes 4 and 5 denote
clinical AD dementia and map to `other`, not AD, because AD is defined neuropathologically here.
`dcfdx_lv` was used in preference to `cogdx` because it is substantially more complete.

**Sample-to-donor resolution.** Genotype samples were resolved to donors per callset: `wgs_harm` via
the ROSMAP WGS QC manifest (`WGS_id` → `projid` → `individualID`), Mayo identity, and the MSBB
biospecimen table (`specimenID` → `individualID`); `divco_hs` via direct identifier, tissue-suffix
stripping, or biospecimen lookup; `wb_dwgs` by identity; `br_dsnwgs` via the AMP-PD sample inventory.
All three `wgs_harm` rules were evaluated against every sample rather than short-circuiting, so a
sample claimed by two studies surfaces as a collision rather than being settled silently by rule
order. Resolution was complete for `divco_hs` (1,019/1,019), `wb_dwgs` (10,418/10,418) and
`br_dsnwgs` (97/97); `wgs_harm` resolved 1,886/1,894 (99.6%), with 8 samples of unknown sex.

Sex was taken per-cohort and never reconciled across sources: each genome takes the sex its own
source study recorded, and cross-source disagreements were emitted unmodified so that the
genotype-based sex check in §3 could adjudicate them. Only unknown sex was left unwritten.

## 3. Per-callset processing

Each callset was processed independently before merging.

**Genotype QC and ancestry (GenoTools v1.3.6).** Per callset: apply the sex-update file, filter to
biallelic PASS SNPs, convert to PLINK bed, project onto a pruned reference panel for ancestry
assignment, and run sample-level QC. Variant identifiers were set to `chr:pos:REF:ALT`.

Ancestry projection assigns each sample to one of 11 strata (EUR, AJ, AAC, AFR, AMR, CAH, CAS, EAS,
FIN, MDE, SAS).

**Variant order is a precondition of ancestry projection.** GenoTools aligns the study genotype
matrix to the reference panel by column position rather than by variant name
(`ancestry.py:247`); the reordering that would enforce name-matching is unreachable during model
training. Where the two orders disagree, every shared column is standardized by another variant's
mean and standard deviation and projected through another variant's loading. Counts still match, no
exception is raised, and the run exits successfully. `br_dsnwgs` was affected — its pgen was ordered
`1,10,11,…,19,2,20,…` from alphabetical concatenation of per-chromosome files, against a numerically
ordered panel — and all 97 samples were assigned to a single stratum with 0.97 nominal model
accuracy. Adding `--sort-vars` at import resolved it. A read-only checker (`scripts/diag_order.py`)
reports rank inversions between a callset and the panel and was run against the other three
callsets, each returning zero inversions.

**Normalization.** All four callsets were rewritten to a common convention — `chr`-prefixed
chromosome codes and `chr:pos:REF:ALT` variant identifiers — so that merge keys are comparable.

**Merge.** The four normalized callsets were merged as a union with PLINK 1.9 (PLINK 2's
non-concatenating `--pmerge-list` is unimplemented in the available build), yielding
**172,497,055 variants × 13,334 samples**.

The §1 table lists 13,428 genomes as supplied, and the four normalized filesets entering the merge
still sum to 13,428 — so no sample was dropped upstream of it. The 94-genome difference is **94
donors carried in both `wgs_harm` and `divco_hs` under the same sample ID**, which PLINK treats as
one individual and fuses into a single merged sample (genotype conflicts set to missing). They are
the only such collisions: every other callset pair is disjoint, and `br_dsnwgs` shares no ID with
any of the other three. 87 of the 94 survive sample exclusion and carry a multi-valued
`source_callset` of `divco_hs|wgs_harm` in the analysis grain; folding them into either parent
callset moves AJ PC1 eta² by at most 0.002 (0.984 → 0.983 / 0.982).

## 4. Relatedness and sample exclusion

Relatedness was estimated with KING-robust on a common variant set built across the full merged
cohort, then computed within ancestry strata. Kinship thresholds: 0.354 for duplicate/monozygotic
pairs, 0.0884 for up to second-degree relatives.

The call-rate threshold defining the common variant set is load-bearing and set to
`--geno 0.005`. It must sit below the smallest callset's share of the cohort (`br_dsnwgs` is
97/13,334 = 0.0073). At a conventional `--geno 0.05`, every variant absent from `br_dsnwgs` was
retained in the common set and all 97 `br_dsnwgs` samples consequently read as ~50% missing — a
value subsequently consumed as a tie-break signal when resolving duplicate pairs. This threshold
requires revisiting if any callset below ~0.5% of the cohort is added.

Samples were excluded for duplication, second-degree-or-closer relatedness, or QC failure:
**12,495 of 13,334 retained** (839 excluded: 499 relatives, 319 duplicates, 21 sex-check failures).

No `--mind` sample-level missingness filter was applied. `br_dsnwgs` is ~50% missing on the common
set for the structural reason given in §1, so a missingness filter would remove all 95 retained
`br_dsnwgs` samples rather than any genuinely low-quality ones.

## 5. Per-ancestry variant QC and principal components

Variant QC and PCA were run within ancestry stratum, since call rate, allele frequency and
Hardy-Weinberg expectation are all ancestry-specific.

Thresholds (fixed across strata): `--geno 0.05`, `--maf 0.01`, `--hwe 1e-6 keep-fewhet`. The
`keep-fewhet` modifier removes heterozygote-excess genotyping artifacts while retaining
heterozygote-deficient sites, which can reflect real population structure.

Principal components were computed on an LD-pruned subset (`--indep-pairwise 1000kb 1 0.1`),
excluding long-range LD and inversion regions via a fixed BED file, retaining 10 PCs. **The
long-range-LD exclusion is applied at the pruning step only**, so it shapes the covariates without
removing those regions from the association set — the MHC is a genuine AD locus and masking it from
association would delete signal the study expects to find.

Six of eleven strata yielded PCs: EUR (n=10,135), AJ (1,518), AAC (226), AFR (191), AMR (185), CAH
(106). CAS, EAS, FIN, MDE and SAS failed LD pruning at their sample sizes and carry no PCs; no
contrast in those strata reached the viability threshold in §7.

## 6. Callset-concordance variant filter

### 6.1 Motivation

Before filtering, ancestry-stratified principal components separated samples by callset rather than
by genetic ancestry. Quantified as η² of callset on each PC, the worst-affected components were
AJ PC1 = 0.984 and EUR PC2 = 0.757. Removing `br_dsnwgs` entirely changed neither value
(0.984 / 0.758), establishing the structure as a genuine `wgs_harm` ↔ `wb_dwgs` effect rather than an
artifact of the smallest callset.

Because PCs enter every association test as covariates, callset structure in the PCs propagates to
every result, and because disease is entangled with callset by construction — AMP-AD supplies the AD
cases, AMP-PD the PD cases — that structure is not separable from the phenotype after the fact.

### 6.2 Design: disease held constant within every comparison

A between-callset allele-frequency difference has two possible causes: technical breakage, or a real
disease effect. A genome-wide filter comparing callsets over all samples cannot distinguish them and
would delete true signal — `rs429358` differs sharply between AD-source and PD-source samples for an
entirely real reason.

Every comparison was therefore made **within a (stratum × diagnosis) cell**. Inside a cell disease
is constant, so a true disease effect cannot generate a between-callset difference and only
technical breakage can. This is what licenses deletion rather than annotation. Two cells carry most
of the weight: controls (`wb_dwgs` vs `wgs_harm` vs `divco_hs`, the axis driving PC1) and AD cases
(`wgs_harm` vs `divco_hs`, lifted versus natively called, disease matched by construction). Cells
below 100 samples per callset were not evaluated. Flags were unioned across cells.

This is not selection on the outcome: the cells are internally disease-constant, and the resulting
list is applied to case-versus-case contrasts.

### 6.3 Two channels

**Frequency (effect-based).** Each callset pair within a cell was compared using `plink --assoc`
— the 1-degree-of-freedom allelic chi-square on the 2×2 allele-count table — with callset membership
as the phenotype. A variant was flagged only if it cleared **both** |ΔAF| > 0.05 **and** z > 5.0
(read from the *P* column as *P* < erfc(z/√2)). Requiring both is a noise calibration: a significance
threshold alone is sample-size dependent, so a 120-sample cell would throw far more flags than a
3,000-sample cell and the per-pair rates would not be comparable. The absolute |ΔAF| floor is
deliberate rather than relative — at MAF ≈ 20% significance already requires |ΔAF| ≈ 0.14 so the floor
never binds, while at MAF ≈ 2% it is reachable at ≈ 0.035, a large relative discordance that is
nonetheless retained, because a filter licensed to delete should err toward keeping.

**Hardy-Weinberg (mechanism-based).** A variant mismapped so that reads pool from two near-identical
genomic locations shows inflated heterozygosity, breaking HWE *in the callset carrying the error*.
The §5 QC computes HWE pooled across callsets, so a deviation confined to `wgs_harm` (1,540 of the
10,135 retained EUR samples) is averaged against a clean majority before the test sees it. This channel
therefore tests HWE **within each callset**, at *P* < 1e-4, in **controls only** — case ascertainment
genuinely distorts HWE at a real disease locus, so testing over an AD-case-heavy callset would flag
APOE for an entirely real reason.

**Excess-over-chance gate.** A (stratum × callset) cell contributes exclusions only when its
rejection count exceeds the number expected by chance at the applied threshold. The ratio of
expected to observed rejections is the Benjamini-Hochberg FDR estimate for that cell, so at or below
1.0× there is nothing to attribute and no multiplier is applied. Measured rates:

| Stratum | Callset | Controls | Rejections | Expected by chance | Ratio |
|---|---|---|---|---|---|
| EUR | `wgs_harm` | 328 | 1,680 | 377 | **4.5×** |
| EUR | `wb_dwgs` | 3,064 | 132 | 377 | 0.35× |
| AJ | `wb_dwgs` | 638 | 97 | 400 | 0.24× |

`wgs_harm` exceeds expectation 4.5-fold from the smallest of the three samples — real heterozygote
excess, consistent with mismapping in the lifted callset. The other two rows sit at or below chance;
underdispersion is expected at these control counts, since the HWE *P* distribution is discrete and
conservative at *n* in the hundreds, which is a second reason a below-chance cell should not
contribute. Gating removed 228 variants from the list.

### 6.4 Result

The final exclusion list contains **4,187 variants**, applied to the association set — not only to
the PCA input, since a variant mismapped badly enough to bend PC1 also produces a spurious
association in the test itself, where no PC adjustment reaches it.

PCA was recomputed on the filtered set within the same job that built the list, so the before/after
comparison is internal to one run and not assembled from two.

Each cell below is the **maximum η² over PC1–PC10** for that stratum, and the worst-affected
component is named alongside it. A maximum over ten values is upward-biased relative to any single
pre-specified component; it is used here as a monitoring statistic, chosen so that no axis can carry
callset structure unreported.

| Stratum | *n* | max η² over PC1–10, unfiltered | filtered | Reduction |
|---|---|---|---|---|
| EUR | 10,135 | 0.757 (PC2) | **0.036** (PC6) | **95.3%** |
| AJ | 1,518 | 0.984 (PC1) | 0.962 (PC1) | 2.2% |
| AAC | 226 | 0.877 (PC2) | 0.876 (PC2) | 0.1% |
| AFR | 191 | 0.708 (PC1) | 0.713 (PC1) | −0.6% |
| CAH | 106 | 0.309 (PC4) | 0.199 (PC4) | 35.6% |
| AMR | 185 | 0.203 (PC8) | 0.170 (PC8) | 16.0% |

Measured on the live 4,187-variant list; the per-stratum and per-PC tables are versioned at
`results/pca/af_filter_effect_gated4187{,_per_pc}.csv`, regenerated by
`review/plot_af_filter_effect.py` from the two manifests step 6 writes in one job.

**In EUR the effect is not a shifted peak but a collapsed profile.** Before filtering, four
components carried appreciable callset structure (PC2 0.757, PC3 0.155, PC1 0.080, PC7 0.038);
afterwards **all ten sit at or below 0.036, and PC2 itself falls to 0.011**. The before and after
entries in the table are therefore different components — worst axis before versus worst axis after,
which is the comparison the claim requires — and the reduction understates the change.

**The 4,187 excluded variants are 0.06% of EUR's 7,538,809 post-QC variants.** That fraction
carried essentially all callset structure in the principal components. The flagged variants are ~7×
enriched for genotype discordance between duplicate sample pairs, where the true genotype is
identical by construction — independent confirmation that they are technical.

**η² is not interpretable to two decimals at the small strata.** AAC, AFR, AMR and CAH have *n* of
106–226 across four or five callset groups, and the estimator is upward-biased and high-variance at
those counts: under a filter that only removes variants, CAH PC6 rose 0.098 → 0.123 and AFR PC6
0.005 → 0.050. Those movements are sampling noise, not filter effects. The EUR and AJ conclusions
rest on *n* = 10,135 and 1,518 and do not depend on this.

**The filter does not resolve AJ**, and the scope of the EUR result should not be extended to it.
No AJ cell was ever evaluated (the largest same-diagnosis callset pair is AJ/control at
`wgs_harm` = 44, below the 100-sample floor), so the list carries no AJ-derived flags and
EUR-derived flags do not transfer. AJ additionally shows no HWE excess (0.24× chance). The leading
hypothesis is that AJ's split reflects real sub-continental structure, in which case PC1 belongs in
the model and the callset-phenotype collinearity is a limitation rather than an artifact. This was
not tested. AAC (0.877) and AFR (0.708) are each driven by a single `br_dsnwgs` sample at 39.7σ and
19.6σ respectively; excluding that one point reduces them to 0.016 and 0.121.

## 7. Association testing

Association was tested per ancestry stratum per contrast with `plink2 --glm`, logistic with
Firth fallback for quasi-separated variants. Covariates were sex and PC1–PC10, variance-standardized.
No mixed model was used: second-degree-and-closer relatives were removed in §4, so samples are
independent, and per-stratum analysis is ancestry-stratified by construction.

Analysis was restricted to common variants (`--maf 0.05`) computed on the samples in each contrast
rather than on the pooled stratum. At these arm sizes, variants below MAF 0.05 are unstable and
batch-artifact-prone; rare and low-frequency variation belongs in a separate burden analysis.

**Differential-missingness filter (per contrast, pre-association).** This is the primary defence
against the cross-cohort confound. A variant called well in one sequencing program and poorly in the
other produces a clean, highly significant, entirely artifactual association; because disease is
entangled with program here, this is the default failure mode rather than an edge case. Per-arm
variant missingness was computed for each contrast and variants were removed when they failed
**both** a 2×2 chi-square at *P* < 1e-4 **and** an absolute between-arm missingness difference
≥ 0.02. Both are required because a chi-square alone scales with *n*: at EUR's sample size a
negligible difference clears any threshold, while at CAH's it may not, and the filter would then
mean something different in every stratum. The 0.02 floor is a meaningful fraction of the 0.05
missingness ceiling already imposed in §5.

Contrasts were enumerated from the harmonized labels, requiring ≥ 20 samples per arm to run and
flagged as adequately powered at ≥ 100 per arm. Arms may be source-restricted (`@amppd`, `@ampad`)
to yield within-program comparisons. Phenotype, covariate and contrast-manifest files were written
by a single upstream definition and read directly by the association step, so the definition of who
is a case exists once.

**44 contrasts ran; 17 met the ≥100-per-arm viability threshold** (EUR 14, AJ 3). Across those 17,
genomic inflation was well controlled: **λ_GC 1.0175–1.0549**, the maximum in `PD_vs_AD`, the
primary and most confounded contrast. No viable contrast required inflation correction.

Among the 27 non-viable contrasts λ ranges more widely, as expected at arm sizes of 20–100: it
reaches 1.3247 in AFR `control_amppd_vs_control_ampad` (22 vs 37), and three AJ contrasts return
λ = 0.0000. The latter is a non-convergence signature rather than a deflated result — those three
(`AD_vs_DLB`, `PD_vs_AD`, `control_amppd_vs_control_ampad`) are the AJ contrasts whose arms sit on
opposite callsets, and every variant returns `ERRCODE=UNFINISHED` with *P* ≈ 1, so the median
chi-square is ~0. AJ PC1 is 0.962 η² on callset (§6.4) and enters as a covariate, so a covariate
separates the arms almost perfectly. All three are below the viability floor and none is reported.

## 8. Control-versus-control artifact scan

A control-versus-control association (AMP-PD controls versus AMP-AD controls, within stratum) has no
true disease signal by construction, so variants reaching significance there are candidate
cohort/batch artifacts. Every genome-wide hit in every primary contrast was annotated with its
control-versus-control *P* value.

**This scan annotates and never subtracts, and the distinction from §6 is a difference in
entitlement rather than a difference of caution.** §6 may delete because disease is held constant
inside each of its cells. This scan may not, because the two control definitions differ across
programs: AMP-AD controls are screened as cognitively normal, AMP-PD controls for absence of PD and
not for absence of AD. A true AD locus is therefore *expected* to separate the two control arms, and
subtracting on that basis would remove real signal.

This was confirmed empirically rather than argued. All four flagged genome-wide hits in EUR
`AD_ampad_vs_control_ampad` — the within-cohort AD contrast, and the cleanest in the study — are
APOE, including both causal variants: `chr19:44908684:T:C` (rs429358, ε4) at *P* = 3.55e-15 and
`chr19:44908822:C:T` (rs7412, ε2). Destructive filtering would have deleted ε4 at 3.55e-15 from the
study's cleanest AD contrast. Two of the three genome-wide control-versus-control hits are
`HLA-DQB1` (chr6:32,661,554 and 32,661,570), ~30 kb from `HLA-DRB1`/`DRB5` — the same screening
asymmetry at the study's other real finding. The third, `chr3:106666502`, is intergenic and 443 kb
from any gene, and is the one plausible genuine artifact.

The scan also corroborates §6: the high-signal contrasts came through with zero control-flagged
variants (`AD_vs_PSP` 2,319 hits / 0 flagged; `PD_vs_control` 2,523 / 0;
`PD_amppd_vs_control_amppd` 2,508 / 0).

The flagging threshold is *P* < 1e-5, deliberately wider than genome-wide significance. All four
APOE control-versus-control *P* values fall between 1.3e-07 and 6.2e-06 — above 5e-8 — so APOE is
not among the three genome-wide control-versus-control hits, and at a 5e-8 flagging threshold the AD
contrasts would lose nothing. The wider net was retained because flagging is annotation: it is the
threshold, not the mechanism, that places APOE in the flagged set.

## 9. Software and environment

PLINK 2.00a6LM; PLINK 1.9 (multi-cohort merge, and the allelic chi-square in §6.3); bcftools;
UCSC `liftOver`; GenoTools v1.3.6; Python 3.11. Gene symbols and coordinates were resolved from a
UCSC refFlat table on GRCh38, verified at load time against APOE's known GRCh38 extent. Analyses ran
on the NIH Biowulf cluster under SLURM.

## 10. Limitations

**No age covariate.** AMP-AD records age at death and AMP-PD age at baseline or analysis. These are
not the same variable, and combining them into one column would silently mix two quantities. Age is
the dominant confounder for both AD and PD, and AD-by-neuropathology donors are older at collection
than the AMP-PD arms. The primary contrast is precisely the one where the two age variables are
least comparable; the within-cohort contrasts are the only ones in which an age term would carry the
same meaning on both arms. Resolving this requires a harmonized age variable that no current
phenotype source emits.

**Disease is confounded with sequencing program by construction.** AMP-AD supplies the AD cases and
AMP-PD the PD cases, so the primary AD-versus-PD contrast cannot separate disease from program by
design. The within-program contrasts (`AD@ampad` vs `control@ampad`, `PD@amppd` vs
`control@amppd`) are the interpretive backbone, and §6, §7's differential-missingness filter, and §8
are three independent defences against the confound rather than one.

**Within-program is not the same as within-callset.** The confound tag measures AMP-PD share and
therefore pools `wb_dwgs` with `br_dsnwgs`. A contrast can score zero program difference and still
have one arm carrying a callset the other lacks entirely. EUR `PD_vs_DLB` is such a case:
`br_dsnwgs`'s 95 retained samples are 71 PD / 21 control / 3 unlabelled and contribute no DLB at
all. It scores Δ = 0.0 and is tagged `within_cohort`, yet the differential-missingness filter
removed **353,068** variants from it — the largest removal of any within-program contrast, the same
order of magnitude as the primary cross-program contrast `PD_vs_AD` (387,207, Δ = 100), and five to
six times every other cross-program contrast (`AD_vs_DLB` 65,245; `AD_vs_control` 63,393;
`control_amppd_vs_control_ampad` 62,344). Roughly 55 `br_dsnwgs` genomes on one arm of a
2,595-sample comparison generated as much technical asymmetry as the entire AMP-AD/AMP-PD split.
`PD_amppd_vs_control_amppd`, which has `br_dsnwgs` on both arms, removed only 5,814, confirming the
mechanism. The results are protected — that contrast returned λ = 1.0402 — but the summary label
understates its technical asymmetry, so per-contrast excluded-variant counts are reported alongside
it.

**`br_dsnwgs` sparsity.** A 97-donor joint call emits nothing at sites monomorphic within its own
donors, so `br_dsnwgs` contributes to roughly half the association set and sits entirely on the
AMP-PD side of the primary contrast.

**AJ callset structure is unresolved.** η² of callset on AJ PC1 remains 0.962 after filtering. AJ
cannot field the primary contrast (AJ/AD is 97, below the 100-per-arm floor), and its three viable
contrasts are all within-cohort, where callset-phenotype collinearity cannot bias the comparison
because both arms share a callset. The sub-continental-structure hypothesis was not tested.

**Five ancestry strata carry no principal components.** CAS, EAS, FIN, MDE and SAS could not be LD
pruned at their sample sizes. No contrast in those strata was adequately powered, so none reached
association testing, but they are absent from the analysis rather than tested and null.

**Definition-dependence of the AD arm.** The neuropathological rule returns null where staging
inputs are missing, and those donors fall out of both arms. The alternative diagnosis columns
available per cohort are not comparable to one another: only ROSMAP's `dcfdx_lv` offers a genuine
clinical-versus-neuropathological contrast, Mayo's `diagnosis` is a second pathology call, MSBB's
CDR measures dementia severity rather than etiology, and Diverse Cohorts' `reag` is a third
pathology instrument. Mayo's Thal phase is ~44% complete against 98.5% for its `diagnosis` column,
which is why the rule nulls 347 of 620 Mayo donors. This is an argument for a sensitivity arm, not
for redefining the primary phenotype.
