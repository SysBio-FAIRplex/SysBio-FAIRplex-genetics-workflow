# demo_sample_check — which pseudobulk donors have genotypes?

A side-quest, self-contained in this directory. The question: of the RNA-seq pseudobulk
donors under `/data/CARD/sysbio/data/beta1/pseudobulk`, how many are also in a genotype
callset — first for the old `AD_versus_PD_GWAS.sh` demo, then for the GWAS running in this
repo.

Nothing here writes outside `demo_sample_check/`. The only things it reads from the parent
repo are `data/*/pgen/*.psam` and two files from `clinical_core_out/`.

## Run it

```bash
./scripts/01_pull_pseudobulk.sh          # rsync from biowulf; only needed to refresh
python3 scripts/02_pseudobulk_ids.py     # 179 sample files -> donor-level table
python3 scripts/03_overlap.py            # donors x callsets
python3 scripts/04_clinical_breakdown.py # sex / pheno / ancestry of the overlap
```

Steps 02–04 are pure local reads and take a couple of seconds. Outputs land in `out/`.

To check a new callset — including the finished GWAS — add one line to `CALLSETS` in
`03_overlap.py` and rerun 03 and 04. That is the whole of "plan step 6"; no logic changes.
Ad-hoc alternative: `python3 scripts/03_overlap.py mylabel=/path/to.psam`.

## The pipeline

| step | reads | writes |
|---|---|---|
| 01 | biowulf `pseudobulk/` | `pseudobulk/` (179 files, tree preserved) |
| 02 | `pseudobulk/` | `out/pseudobulk_donors.tsv`, `out/pseudobulk_donors_unique.txt` |
| 03 | 02 + `clinical_core_out/genome_crosswalk.csv` + callset psams | `out/overlap.tsv`, `out/unmatched_donors.txt` |
| 04 | 03 + `clinical_core_out/analysis_grain.csv` | `out/clinical_breakdown.tsv` |

The directory tree under `pseudobulk/` is **preserved, never flattened**. Two of the six
cohorts (CMD, rasle) put the cell type in the *directory* and name every file
`samples.tsv`/`samples.csv`; the other four put it in the *filename*. Flattening would both
collide and destroy the cohort/tissue/celltype labels step 02 attributes each donor to.

## Where the IDs come from

This is the whole difficulty. It is not a set intersection.

**Pseudobulk side** — six cohorts, six donor columns. Step 02 tries an ordered candidate
list and takes the first present, so a new cohort needs no code change.

| cohort | donor column | example | granularity |
|---|---|---|---|
| CMD | `donor_id`, or `sample_id` in 3 files | `<alnum>`, `<n>-<n>`, `HL<6d>` | **cell-level** |
| PD | `participant_id` | `PM-<site>_<n>_<n>` | donor |
| diversecohorts | `individualID` | `R<7d>` | donor |
| mitrosmap | `individualID` (+`projid`) | `R<7d>` | donor |
| rasle | `donor` | `<n>-<n>` | donor |
| rosmap | `individualID` (+`projid`) | `R<7d>` | donor |

The CMD files are **cell-level**: `ID` is a 10x barcode and
`hypothalamus/c1_4_oligo_mature/samples.tsv` holds 144,247 rows for 11 donors. Taking `ID`
would inflate the denominator ~1000×. Every cohort is deduplicated to distinct donors per
file: 179 files → 14,410 (cohort, tissue, celltype, donor) rows → **1,098 unique donors**.

**Genotype side** — a `.psam` IID is not a person. Each group named samples after the
specimen, the assay, or a sequencing manifest entry, so one donor is `MAP<8d>` in
AMP-AD WGS_Harmonization, `R<7d>` in every ROSMAP-derived pseudobulk, and
`<individualID>_DLPFC_WGS` in DivCo_HS.

**The bridge** is `clinical_core_out/genome_crosswalk.csv`, which `clinical_core.py`
already builds and the rest of this project already depends on — 13,428 IIDs, each with the
`individual_id` it belongs to and the `rule` that established it:

| rule | n | example |
|---|---|---|
| `wb_dwgs:psam_identity` | 10,418 | `BF-<4d>` → `BF-<4d>` |
| `rosmap:wgs_qc_manifest` | 1,196 | `MAP<8d>` → `R<7d>` |
| `divco:individualID_direct` | 620 | `<individualID>` → `<individualID>` |
| `mayo:identity` | 349 | `<individualID>` → `<individualID>` |
| `msbb:biospecimen_specimenID` | 349 | `<specimenID>` → `AMPAD_MSSM_<10d>` |
| `divco:specimenID_suffix_strip` | 293 | `<individualID>_DLPFC_WGS` → `<individualID>` |
| `divco:biospecimen_lookup` | 106 | `<specimenID>-D` → `<individualID>` |
| `br_dsnwgs:sample_inventory` | 97 | `PM-MS_<5d>-BLM0-PVC-DWGS` → `PM-MS_<5d>` |

Step 03 **re-derives none of this**. An earlier version did, with hand-written regexes, and
recovered 317 donors where the crosswalk finds 716. Most of the gap is ROSMAP, where
`MAP<n>` → `R<m>` comes from a manifest and is not obtainable by string surgery at all. If
you are tempted to add a normalization rule to `03_overlap.py`, add it to the crosswalk in
`clinical_core.py` instead, so the whole project benefits and there is one answer.

## Results

**716 / 1,098** donors reach a current callset; **520** reach the demo.

```
cohort           donors   demo_PD_EUR  demo_AD_EUR   BR-DSNWGS   WGS_Harm   DivCo_HS
CMD                 174             0            0           0          0          0
PD                   97             0            0          94          0          0
diversecohorts      167             0           69           0         73         88
mitrosmap            48             0           44           0         45         12
rasle               186             0            0           0          0          0
rosmap              450             0          429           0        439        122
ALL unique         1098             0          520          94        535        216
callset size                     8607         1605          97       1894       1019
```

The demo's reach is a **strict subset** of the current one — 0 donors reach the demo but not
the current callsets. The current run adds **196**:

| source of the gain | n | why the demo missed them |
|---|---|---|
| postmortem AMP-PD (BR-DSNWGS) | 94 | demo's PD arm was WB-DWGS, the living cohorts |
| reachable only via DivCo_HS | 87 | a callset the demo never used |
| in WGS_Harm but not the demo psam | 15 | dropped before the demo's EUR split |

Of those 196: 133 EUR, 40 AFR, 23 not in the grain at all.

Per-cohort rows can sum to more than the unique total: 13 donors are in both `rosmap` and
`diversecohorts`. Step 03 counts such a donor in each cohort it belongs to; step 04 instead
groups them into an explicit combined bucket (`diversecohorts,rosmap`) so its rows stay
mutually exclusive and sum to the total.

After the pipeline's own QC (`analysis_grain.csv`), **655** of the 716 remain; 61 were
dropped by ancestry/relatedness/call-rate filtering. Of those 655: 411 female / 244 male;
278 AD, 175 control, 137 other, 62 PD, 3 blank; 507 EUR / 148 AFR.

### The demo's PD arm reaches zero, and that is structural

`FILTERED.AMP_PD_EUR.psam` is 8,607 samples from AMP-PD's **living** cohorts (`LB`, `PD`,
`PP`, `HB`, `LC`, `SY`, `SU`, `BF`). The only AMP-PD pseudobulk is the **postmortem** arm,
which lives in BR-DSNWGS. The two cannot intersect. All 520 of the demo's reach is its AD
arm.

`FILTERED.AMP_AD_EUR.psam` is the **pre-QC** EUR pool, not the demo's actual input — the
demo ran on `FILTERED.AMP_AD_EUR_callrate_related_het_geno_haplotype` and additionally
`--remove`d `FILTERED.AMP_AD_EUR.related`. So 520 is the demo's *ceiling*, and its real N
is lower. The full QC chain is on biowulf in
`/data/CARD/AD/AMP_AD/jointGenotypingROSMAPMayoRNAseqMSBB/genotools/` as six psams
(`_callrate`, `_callrate_related`, `_callrate_related_het`, `_callrate_related_het_geno`,
`…_haplotype`) if the per-stage attrition is ever worth measuring.

Every demo IID resolved through the crosswalk (`not in xwalk` = 0), so 520 is not an
artifact of unmapped IDs.

### The 382 unmatched are mostly not a mapping failure

**CMD (174) and rasle (186) have no genotypes anywhere in this project, by design.**

- **rasle** is AMP **RA/SLE** — a different Accelerating Medicines Partnership program
  (Rheumatoid Arthritis and Lupus), not AMP-PD. 186 donors, kidney, 1,517 SLE vs 258
  control, IDs `200-NNNN`. Verified: AMP-PD's nine constituent studies are LBD, PPMI,
  PDBP, HBS, LCC, Steady, Sure, BioFIND and Post Mortem Brain; zero `NNN-NNNN` values
  appear in `study_participant_id` (the field an outside study's native ID lands in) or
  anywhere else across the four AMP-PD metadata files or the 13,428-row crosswalk.
- **CMD** is external reference tissue atlases — `dataset=Gaulton_FNIH`,
  `dataset_role=Reference`; kidney is KPMP (`dataSource=KPMP`, `<n>-<n>`), liver `HL######`,
  hypothalamus `<alnum>` (Tadross).

That leaves **22 genuine misses** worth chasing: PD 3, rosmap 7, diversecohorts 10,
mitrosmap 2. Listed in `out/unmatched_donors.txt`.

## Three things found along the way that are worth acting on

**1. AMP-PD postmortem: 97 and 97, but not the same 97.** The BR-DSNWGS pgen has 97
samples and the PD pseudobulk has 97 donors — the overlap is **94**. Three donors have
pseudobulk but never made the pgen — all three are present in the WGS inventory, so they
were lost between inventory and delivery; three others have genotypes but no pseudobulk.
(Both trios are recoverable from `out/`, which stays out of git; the IDs are
controlled-access and are deliberately not reproduced here.)
Comparing counts alone reads as a perfect 97/97. After QC, 84 remain.

**2. Ancestry is unstable — across callsets, and across genotools runs.** Two independent
symptoms of the same thing:

*Within the current run*: 178 individuals carry more than one genotype IID (163
`divco_hs`+`wgs_harm`, 15 within `wgs_harm`). Of those, **58 get different ancestry labels
depending on the callset**, always AFR vs EUR. `sex`, `pheno` and `dx_detailed` never
conflict — only ancestry.

*Between the demo and now*: the demo's psam is **EUR-only**, yet of the 520 donors it
reaches, the current pipeline calls **90 AFR** and is itself split on another **34**
(`AFR|EUR`) — the same conflict set. Only 358 are unambiguously EUR both times. So roughly
a quarter of the demo's "EUR" cohort would not be classified EUR today.

Since the GWAS is ancestry-stratified, this decides which stratum a donor lands in, and it
is not reproducible between runs. This is a main-pipeline concern, not a side-quest one.
`04_clinical_breakdown.py` reports conflicts rather than silently taking the first row.

**3. Three overlapping AMP-PD donors have a blank phenotype** in the grain — genotyped,
pseudobulked, and unlabelled. Worth resolving before relying on the case/control split.
