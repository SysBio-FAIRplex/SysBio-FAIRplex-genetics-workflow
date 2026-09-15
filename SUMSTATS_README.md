# AD vs PD WGS GWAS — summary statistics

Case-control association results from whole-genome sequence data across AMP-AD and AMP-PD,
**with AD defined by neuropathology, never by clinical dementia code**. GRCh38.

Aggregate summary statistics only — no genotypes, no sample identifiers.

## Contents

- **44 files**, `gwas_<ANCESTRY>_<CONTRAST>.filtered.tsv` — every contrast that ran, i.e. every
  contrast with ≥20 samples in both arms. Smallest arm in the release is 21.
- **`gwas_summary.csv`** — the table below in machine-readable form, same numbers.

Six ancestry strata are represented. CAS, EAS, FIN, MDE and SAS could not be LD-pruned at their
sample sizes, carry no principal components, and produced no contrast.

## All 44 contrasts

`†` = below the ≥100-per-arm threshold we treat as adequately powered. Released, but underpowered.
**Program** is the AMP-AD/AMP-PD split: `within` = both arms from one program, `cross` = one arm
each, `partial` = mixed. **Diff-miss removed** = variants dropped by the per-contrast
differential-missingness filter. **GW hits** = variants at P < 5e-8.

| Ancestry | Contrast | Case arm | *n* | Control arm | *n* | Program | λ_GC | Diff-miss removed | GW hits |
|---|---|---|---:|---|---:|---|---:|---:|---:|
| EUR | `PD_amppd_vs_control_amppd` | `PD@amppd` | 2,595 | `control@amppd` | 3,077 | within | 1.0376 | 5,814 | 2,508 |
| EUR | `PD_vs_control` | `PD` | 2,595 | `control` | 3,484 | within | 1.0379 | 6,331 | 2,523 |
| EUR | `PD_vs_DLB` | `PD` | 2,595 | `DLB` | 2,471 | within | 1.0402 | 353,068 | 318 |
| EUR | `AD_vs_DLB` | `AD` | 820 | `DLB` | 2,471 | cross | 1.0277 | 65,245 | 33 |
| EUR | `AD_vs_control` | `AD` | 820 | `control` | 3,484 | cross | 1.0414 | 63,393 | 34 |
| EUR | `PD_vs_AD` | `PD` | 2,595 | `AD` | 820 | cross | 1.0549 | 387,207 | 49 |
| EUR | `AD_ampad_vs_control_ampad` | `AD@ampad` | 820 | `control@ampad` | 407 | within | 1.0175 | 1 | 23 |
| EUR | `control_amppd_vs_control_ampad` | `control@amppd` | 3,077 | `control@ampad` | 407 | cross | 1.0319 | 62,344 | 3 |
| EUR | `AD_vs_other` | `AD` | 820 | `other` | 395 | partial | 1.0293 | 88 | 0 |
| EUR | `PD_vs_other` | `PD` | 2,595 | `other` | 395 | partial | 1.0198 | 60,071 | 0 |
| EUR | `AD_vs_MCI` | `AD` | 820 | `MCI` | 159 | within | 1.0272 | 14 | 0 |
| EUR | `PD_vs_MCI` | `PD` | 2,595 | `MCI` | 159 | cross | 1.0277 | 58,107 | 3 |
| EUR | `AD_vs_PSP` | `AD` | 820 | `PSP` | 135 | partial | 1.0323 | 56,330 | 2,319 |
| EUR | `PD_vs_PSP` | `PD` | 2,595 | `PSP` | 135 | partial | 1.0189 | 356 | 1 |
| AJ | `PD_amppd_vs_control_amppd` | `PD@amppd` | 543 | `control@amppd` | 642 | within | 1.0252 | 61 | 82 |
| AJ | `PD_vs_control` | `PD` | 543 | `control` | 688 | within | 1.0256 | 72 | 59 |
| AJ | `PD_vs_DLB` | `PD` | 543 | `DLB` | 114 | within | 1.0195 | 16 | 255 |
| AJ | `AD_vs_DLB` † | `AD` | 97 | `DLB` | 114 | cross | **0.0000** ⚠ | 175 | 0 |
| AJ | `AD_vs_control` † | `AD` | 97 | `control` | 688 | cross | 1.0948 | 593 | 0 |
| AJ | `PD_vs_AD` † | `PD` | 543 | `AD` | 97 | cross | **0.0000** ⚠ | 577 | 0 |
| AJ | `AD_vs_other` † | `AD` | 97 | `other` | 65 | partial | 1.1113 | 4 | 0 |
| AJ | `PD_vs_other` † | `PD` | 543 | `other` | 65 | partial | 1.0344 | 101,761 | 0 |
| AJ | `AD_ampad_vs_control_ampad` † | `AD@ampad` | 97 | `control@ampad` | 46 | within | 1.1094 | 0 | 0 |
| AJ | `control_amppd_vs_control_ampad` † | `control@amppd` | 642 | `control@ampad` | 46 | cross | **0.0000** ⚠ | 645 | 0 |
| AAC | `AD_vs_control` † | `AD` | 79 | `control` | 81 | partial | 1.1023 | 0 | 1 |
| AAC | `AD_vs_other` † | `AD` | 79 | `other` | 39 | within | 1.1261 | 0 | 0 |
| AAC | `AD_ampad_vs_control_ampad` † | `AD@ampad` | 79 | `control@ampad` | 34 | within | 1.1562 | 0 | 0 |
| AAC | `control_amppd_vs_control_ampad` † | `control@amppd` | 47 | `control@ampad` | 34 | cross | 1.2264 | 0 | 0 |
| AAC | `PD_amppd_vs_control_amppd` † | `PD@amppd` | 25 | `control@amppd` | 47 | within | 1.2110 | 0 | 0 |
| AAC | `PD_vs_AD` † | `PD` | 25 | `AD` | 79 | cross | 1.1344 | 117 | 0 |
| AAC | `PD_vs_control` † | `PD` | 25 | `control` | 81 | partial | 1.1545 | 29 | 0 |
| AAC | `PD_vs_other` † | `PD` | 25 | `other` | 39 | cross | 1.2484 | 0 | 0 |
| AFR | `AD_vs_control` † | `AD` | 81 | `control` | 59 | partial | 1.0988 | 0 | 0 |
| AFR | `AD_ampad_vs_control_ampad` † | `AD@ampad` | 81 | `control@ampad` | 37 | within | 1.1221 | 0 | 0 |
| AFR | `AD_vs_other` † | `AD` | 81 | `other` | 36 | within | 1.1212 | 3 | 0 |
| AFR | `control_amppd_vs_control_ampad` † | `control@amppd` | 22 | `control@ampad` | 37 | cross | 1.3247 | 0 | 0 |
| AMR | `PD_vs_AD` † | `PD` | 46 | `AD` | 75 | cross | 1.1182 | 2 | 0 |
| AMR | `AD_vs_control` † | `AD` | 75 | `control` | 40 | partial | 1.1368 | 0 | 0 |
| AMR | `PD_vs_control` † | `PD` | 46 | `control` | 40 | partial | 1.1303 | 0 | 0 |
| AMR | `PD_amppd_vs_control_amppd` † | `PD@amppd` | 46 | `control@amppd` | 24 | within | 1.2497 | 0 | 0 |
| AMR | `AD_vs_other` † | `AD` | 75 | `other` | 23 | within | 1.1767 | 7 | 0 |
| AMR | `PD_vs_other` † | `PD` | 46 | `other` | 23 | cross | 1.2597 | 1 | 0 |
| CAH | `AD_vs_control` † | `AD` | 50 | `control` | 22 | partial | 1.2223 | 0 | 0 |
| CAH | `AD_vs_other` † | `AD` | 50 | `other` | 21 | within | 1.2383 | 0 | 0 |

### Reading the table

- **λ_GC = 0.0000 ⚠ means the model did not converge**, not that the result is deflated. In all
  three cases every variant carries `ERRCODE=UNFINISHED` with P ≈ 1. These are AJ's cross-callset
  contrasts, where a covariate (AJ PC1, η² 0.962 on callset) separates the arms almost perfectly.
  **Do not use these three files.**
- **λ_GC above ~1.1 is small-sample noise**, not confounding — every such contrast has an arm under
  100. Nothing here required genomic-control correction.
- **`within` does not mean within-*callset*.** The tag measures the AMP-AD/AMP-PD program split and
  pools AMP-PD's two callsets. EUR `PD_vs_DLB` scores `within` yet lost 353,068 variants, because
  the postmortem-brain callset contributes PD but no DLB. **Use `Diff-miss removed` as the real
  technical-asymmetry signal, not the tag.**
- **`other`** is a heterogeneous residual arm (classified as neither case nor control), not a clean
  comparison group. **`control`** pools both programs' controls, which are differently screened —
  see Limitations.

### Suggested starting points

| use | contrast |
|---|---|
| cleanest AD signal | `gwas_EUR_AD_ampad_vs_control_ampad` — both arms within AMP-AD, λ 1.0175 |
| cleanest PD signal | `gwas_EUR_PD_amppd_vs_control_amppd` — both arms within AMP-PD, λ 1.0376 |
| primary contrast | `gwas_EUR_PD_vs_AD` — λ 1.0549, confounded with program by design |
| negative control | `gwas_EUR_control_amppd_vs_control_ampad` — controls vs controls, no true disease signal. Check any locus of interest against it |

## Columns

Raw PLINK 2 `--glm` output. **Read by header name, never by position.**

| column | meaning |
|---|---|
| `#CHROM`, `POS` | GRCh38 coordinates |
| `ID` | `chr:pos:REF:ALT` |
| `REF`, `ALT` | reference and alternate allele |
| `PROVISIONAL_REF?` | PLINK internal flag — **not a frequency**. Column 5 has been misread as one |
| `A1` | **the effect allele.** `BETA` is per copy of `A1` |
| `OMITTED` | the allele `A1` was tested against |
| `A1_FREQ` | frequency of `A1` |
| `FIRTH?` | `Y` if the Firth fallback was used |
| `TEST` | always `ADD` |
| `OBS_CT` | samples contributing to this variant |
| `BETA`, `SE`, `Z_STAT`, `P` | log-odds, standard error, z, p |
| `ERRCODE` | `.` is a clean fit. **Anything else did not converge — drop those rows.** A non-converged fit still carries a numeric `P` and survives a p-value filter |

## Source data — where every input came from

Verbatim from the acquisition record (`wgs_core.ipynb` §0). All pulls run on the NIH Helix
transfer node. Both AMP-AD and AMP-PD sources are controlled-access.

### Genotypes

| Callset | *n* | Source |
|---|---:|---|
| `wgs_harm` | 1,894 | Synapse folder **`syn11707420`** — 24 per-chromosome GRCh37 VCFs (1–22, X, Y).<br>`synapse get -r syn11707420` |
| `divco_hs` | 1,019 | Synapse **`syn68260951`** (`merged.deduped.vcf.gz`) + **`syn68260952`** (`.tbi`), single merged GRCh38 VCF.<br>`synapse get syn68260951` |
| `wb_dwgs` | 10,418 | `gs://amp-pd-genomics/releases/2023_v4release_1027/wgs-WB-DWGS/plink/pfiles/all_chrs_merged.{pgen,psam,pvar,log}` — AMP-PD ships PLINK 2 pfiles already on GRCh38.<br>`gcloud storage cp … --billing-project <project>` |
| `br_dsnwgs` | 97 | `gs://amp-pd-receipt-2026/mssm/20260626-transfer-mssm-jvcf/*` — single 97-donor joint-call GRCh38 VCF.<br>`gcloud storage cp … --billing-project <project>` |

### Clinical metadata

| Directory | Pulled from | Files it holds |
|---|---|---|
| `WGS_Harmonization/metadata` | Synapse `syn73713766`, `syn73713767`, `syn73713768`, `syn12178037`, `syn21893059` | `MayoRNAseq_individual_metadata_harmonized.csv`, `MSBB_individual_metadata_harmonized.csv`, `MSBB_biospecimen_metadata.csv`, `ROSMAP_clinical_harmonized.csv`, `WGS_sample_QC_info.csv` |
| `DivCo_HS/metadata` | Synapse `syn51757644`, `syn51757645`, `syn51757646` | `AMP-AD_DiverseCohorts_individual_metadata.csv`, `AMP-AD_DiverseCohorts_biospecimen_metadata.csv`, `AMP-AD_DiverseCohorts_assay_WGS_metadata.csv`, and two further files (see Gaps) |
| `amp-pd-genomics/metadata` | `gs://amp-pd-data/releases/2023_v4release_1027/` | `clinical/Demographics.csv`, `amp_pd_participants.csv`, `amp_pd_case_control.csv`, `wgs_BR-DSNWGS_sample_inventory.csv` |

Two Synapse documents define the phenotype rules rather than supplying data: the AMP-AD
diagnosis-criteria specification **`syn51757663`** and the harmonized data dictionary
**`syn73713784`**.

### Gaps in the acquisition record

Stated rather than papered over. None affects the released statistics; all three are
reproducibility gaps.

- **The `WGS_Harmonization` Synapse IDs are recorded as a set, not mapped one-to-one to filenames.**
  The five IDs and the five files are both listed above; the pairing is not written down.
- **`DivCo_HS/metadata` holds five files and only three Synapse IDs are recorded.** The other two
  have no recorded provenance.
- **`wgs_WB-DWGS_samples.csv` is on disk but is not in the acquisition record's AMP-PD file list.**

The GenoTools ancestry reference panel (`ref_panel_gp2_prune_rm_underperform_pos_update`, a pruned
GP2 panel) and the GRCh37→GRCh38 liftover chain are third-party reference data; their acquisition
is documented separately from this release.

## QC pipeline

Every step, in order, with the parameters actually used.

**0 — Input callsets.** 13,428 genomes across four independently generated callsets:
`wgs_harm` (AMP-AD: ROSMAP, MayoRNAseq, MSBB; 1,894; GRCh37, lifted to GRCh38 with UCSC `liftOver`),
`divco_hs` (AMP-AD Diverse Cohorts; 1,019; native GRCh38),
`wb_dwgs` (AMP-PD whole blood; 10,418; native GRCh38),
`br_dsnwgs` (AMP-PD postmortem brain; 97; native GRCh38).
Sources and acquisition commands above.

**0b — VCF → pgen.** `wgs_harm`: per-chromosome `plink2` filter → `--pmerge-list` → `liftOver`
b37→hg38 on the merged file → `--chr 1-22,X,Y`. `divco_hs` and `br_dsnwgs`: `bcftools` filter →
`plink2 --make-pgen`. `wb_dwgs`: none, it ships as pfiles.

**1 — Phenotype definition (neuropathological).**
- ROSMAP, MSBB: AD = Braak ≥ IV **and** CERAD moderate/frequent; control = Braak ≤ III **and** CERAD sparse/none.
- MayoRNAseq: AD = Braak ≥ IV **and** Thal ≥ 2; control = Braak ≤ III **and** Thal < 2. Any neuropathological control whose independent `diagnosis` field named a disease was demoted to `other` (control-purity screen; one-directional — never rescues a null, never overrides an AD call).
- Diverse Cohorts: pre-adjudicated `ADoutcome` (`mayoDx` for the Mayo contribution group).
- AMP-PD: curated `case_control_other_latest` → PD / control; `Other` excluded.
- Missing staging inputs yield **null, not `other`** — those donors are excluded from every arm rather than pooled.
- Braak alone was never used for any cohort. Clinical dementia codes never define an AD case.
- `dx_detailed` (PD / AD / MCI / DLB / PSP / control / other) from `diagnosis_latest` (AMP-PD), `dcfdx_lv` (ROSMAP), `diagnosis` (Mayo).

**2 — Per-callset QC and ancestry** (GenoTools v1.3.6, per callset independently):
- Apply per-callset sex-update file derived from harmonized clinical metadata.
- Filter to biallelic PASS SNPs. (`divco_hs` has no VQSR — FILTER is all `.` — so no PASS filter was applied to it.)
- Convert to PLINK bed; `--sort-vars` at import.
- Ancestry projection onto a pruned GP2 reference panel → one of 11 strata (`--ancestry --ref_panel --ref_labels --all_sample`).
- GenoTools sample-level QC.
- **94 genomes dropped here**, which is why the callset table sums to 13,428 and the merge is 13,334.

**3 — Normalization.** `plink2 --fa <GRCh38> --ref-from-fa force --output-chr chrM --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 1000 missing`, so merge keys are comparable across callsets.

**4 — Merge.** Union merge with PLINK 1.9 → **172,497,055 variants × 13,334 samples**.

**5 — Relatedness.** KING-robust, computed within ancestry stratum on a common variant set built at `--geno 0.005`. Thresholds: **0.354** duplicate/MZ, **0.0884** second-degree-or-closer.
- The `--geno 0.005` value is load-bearing: it must sit below the smallest callset's share of the cohort (97/13,334 = 0.0073). At a conventional 0.05, every variant absent from the 97-donor callset stays in the common set and those samples read as ~50% missing.

**6 — Sample exclusion.** 839 removed → **12,495 retained**: 499 relatives, 319 duplicates, 21 sex-check failures. **No `--mind`** — the 97-donor joint call is structurally ~50% missing on the common set, so a missingness filter would remove those samples rather than any genuinely low-quality ones.

**7 — Per-stratum variant QC and PCA.** Within each ancestry stratum, since call rate, allele frequency and HWE are all ancestry-specific:
- `--geno 0.05 --maf 0.01 --hwe 1e-6 keep-fewhet` (`keep-fewhet` removes heterozygote-excess artifacts while retaining heterozygote-deficient sites).
- 10 PCs on an LD-pruned subset, `--indep-pairwise 1000kb 1 0.1`.
- Long-range-LD and inversion regions excluded **from pruning only**, so the MHC stays in the association set.
- Strata with PCs: EUR 10,135 · AJ 1,518 · AAC 226 · AFR 191 · AMR 185 · CAH 106.

**8 — Callset-concordance variant filter — 4,187 variants deleted.** Comparisons made **inside each (stratum × diagnosis) cell**, so disease is constant and a between-callset difference can only be technical:
- *Frequency channel:* PLINK 1.9 `--assoc` (1-df allelic chi-square) with callset as the phenotype. Flagged only if **both** |ΔAF| > 0.05 **and** z > 5.0.
- *HWE channel:* HWE tested **within each callset**, P < 1e-4, **controls only** (case ascertainment genuinely distorts HWE at a real disease locus). Requires ≥50 controls in the cell.
- Cells with fewer than 100 samples per callset were not evaluated. Flags unioned across cells.
- *Excess-over-chance gate:* a (stratum × callset) cell contributes only when its rejection count exceeds chance at the applied threshold. This removed 228 variants from the list (4,415 → 4,187).
- Effect: EUR's callset η² on the worst principal component falls from 0.757 to 0.036.
- **These variants are absent from every file here** — missing, not non-significant.

**9 — Association testing.** `plink2 --glm hide-covar firth-fallback cols=+a1freq,+beta`, logistic with Firth fallback:
- Covariates: sex and PC1–PC10, variance-standardized.
- `--maf 0.05`, computed on the samples in each contrast rather than pooled.
- No mixed model — second-degree-and-closer relatives already removed, analysis ancestry-stratified.
- *Differential-missingness filter, per contrast, pre-association:* variants removed when they failed **both** a 2×2 chi-square at P < 1e-4 **and** an absolute between-arm missingness difference ≥ 0.02. Counts per contrast in the table above.
- Contrasts require ≥20 per arm to run; ≥100 per arm is flagged as adequately powered.
- **Released rows are `ADD` only with |BETA| ≤ 5**, which drops non-converged Firth fits carrying implausible effect sizes.

**10 — Control-versus-control scan (annotate, never subtract).** Variants significant between the two control arms are candidate cohort artifacts, but the arms are differentially screened — AMP-AD controls assessed as cognitively normal, AMP-PD controls for absence of PD and *not* AD — so a true AD locus is **expected** to separate them. Subtracting on that basis would have deleted APOE-ε4 (P = 3.55e-15) from the cleanest AD contrast. Judge per variant against `gwas_EUR_control_amppd_vs_control_ampad`.

## Limitations

- **No age covariate.** AMP-AD records age at death, AMP-PD age at baseline or analysis — not the same variable, and combining them would silently mix two quantities. Age is the dominant confounder for both diseases, and `PD_vs_AD` is exactly where the two are least comparable.
- **Disease is confounded with sequencing program by construction.** AMP-AD supplies the AD cases, AMP-PD the PD cases. The within-program contrasts are the interpretive backbone.
- **Control arms are not equivalent across programs.** AMP-AD controls are screened as cognitively normal; AMP-PD controls for absence of PD, not of AD. Any `*_vs_control` contrast pooling both inherits this.
- **Postmortem-brain callset sparsity.** A 97-donor joint call emits nothing at sites monomorphic within its own donors, so it contributes to roughly half the association set and sits entirely on the AMP-PD side of the primary contrast.
- **AJ callset structure is unresolved.** η² of callset on AJ PC1 is 0.962 after filtering. AJ's three well-powered contrasts are all within-cohort, where callset-phenotype collinearity cannot bias the comparison; its cross-callset contrasts are the three that failed to converge.
- **The 27 contrasts marked †** have an arm under 100. They are released so you can pick and choose, not because they are powered.

## Access

Released to **`gs://sysbio-gwas/results`**; the GWAS browser reads from there.

Derived from AMP-AD (Synapse) and AMP-PD controlled-access data. Redistribution is governed by
those data use agreements — confirm before widening bucket access.
