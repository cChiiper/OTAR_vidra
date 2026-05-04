# VIDRA_5

## Project Status & Provenance

The goal of this pipeline is to run analyses using the most recent available datasets (OTAR 26.03 and AZ v5). The scripts have been developed iteratively and largely incorporate AI-assisted code generation (notably using GPT-5.4). While functional and already producing meaningful results, most parts of the codebase have not yet undergone full manual refactoring and review.

`VIDRA_5` builds upon the original implementation by [Luca Stefanucci]((https://github.com/TrynkaLab/VIDRA/tree/main)), with substantial adaptations inspired by the fork by [Daniel Considine](https://github.com/Daniel-Considine/VIDRA) (referred here as `VIDRA_2`). The current naming (`VIDRA_5`) is provisional and reflects ongoing development rather than a formal versioning scheme.

At present, the codebase remains somewhat heterogeneous, and some legacy components persist. Future updates will focus on:

- streamlining the pipeline structure  
- removing deprecated steps  
- improving code clarity and consistency (this includes this README file)
- refining the analysis framework  
  - including coding variants irrespective of colocalisation status  
  - improving uncertainty handling in the Stan model

The pipeline is, however, stable in its current form and suitable for running end-to-end analyses locally or on google cloud.

The existing GitHub Actions workflow configuration is no longer actively maintained and may be updated or removed in future revisions.

Pipeline paths and execution logic are managed via a `Makefile`. This design choice reflects the modular nature of the workflow, where different stages are built and executed independently, enabling flexible development and testing from scratch.

---

## Data Availability

All input files are publicly available **except** for the AstraZeneca-to-EFO mapping file:

- `az_phewas_to_efo.tsv`

Access to this file may be requested via Open Targets.

## Workflow

`VIDRA_5` is a hybrid VIDRA pipeline.

The expensive preparation and annotation stages are launched explicitly, in the
same spirit as `VIDRA_2`. The Bayesian stage is the only Nextflow workflow and
runs one task per gene.

```text
Step 1: prepare analysis input
  local:  Docker + local Spark
  gcloud: Dataproc Serverless Spark
   ↓
Step 2: annotate coding/ClinVar variants
  local:  Docker + local VEP resources
  gcloud: standalone GCE annotation VM
   ↓
Step 3: Bayesian analysis
  local:  Nextflow local executor
  gcloud: Nextflow + Google Batch, one task per gene
```

Step 1 and Step 2 are outside Nextflow on purpose. Nextflow only sees existing
`vidra_analysis_ready/as_gene=*` partitions and the shared annotation parquet.

## Main Outputs

- `vidra_analysis_ready/`: Step 1 analysis-ready rows, partitioned by `as_gene`.
- `vidra_analysis_ready_manifest/`: Step 1 variant manifest.
- `GWAS_coding_variants_from_cs.parquet`: Step 1 handoff to Step 2 (this is legacy, will be removed).
- `clinvar_variants.parquet`: Step 1 handoff to Step 2 (this is legacy, will be removed).
- `variant_annotations/annotations.parquet`: Step 2 annotation output.
- `vidra_results/as_gene=<gene>/`: Step 3 Bayesian output partitions.
- `pipeline_info/`: Nextflow reports and gcloud annotation logs where applicable.

## Configuration Files

- `conf/make.defaults.mk` contains Make-only launch defaults for local and
  gcloud Step 1/2 commands, image names, paths, VM sizes, and helper targets.
- `conf/make.local.mk` is optional. The Makefile includes it if present; use it
  for machine-specific overrides instead of editing committed defaults. Not present on the repo, may be removed in the future.
- `conf/*.config` are Nextflow configs for Step 3 only.

The Step 1 CLI can run without `lead_variant_effect` and without raw AZ burden,
but the current Make defaults enable both for local and gcloud runs. To run the
release-only fallback path, override these variables to empty in or on the command line.

## Setup Targets

Build all local images:

```bash
make docker_build
```

Push all gcloud images after building:

```bash
make docker_push_gcloud
```

Download the Open Targets release folders needed by local Step 1:

```bash
make local_sync_ot_release_inputs
```

Convert the local AZ rare-variant input from `.xz` to Spark-readable `.bz2`.
This is a one-time helper; if the `.bz2` already exists it does nothing:

```bash
make local_convert_az_variants_bz2
```

Convert the local AZ PheWAS binary and quantitative `.csv.xz` inputs into
parquet datasets using the paths configured in `conf/make.defaults.mk` and
optionally overridden in `conf/make.local.mk`:

```bash
make local_convert_az_phewas_parquet
```

Upload the VEP cache archive and small gcloud test gene list:

```bash
make gcloud_upload_vep_cache
make gcloud_upload_testdata
```

Run Python tests and Nextflow config checks:

```bash
make test
make config_test
```

## Running The Pipeline

Local smoke run:

```bash
make local_step1
make local_step2
make nextflow_local_step3
```

Gcloud Step 1 and Step 2 use the Make variables in `conf/make.defaults.mk`.
Those defaults currently point to the production output root
`gs://vidra-2-0/nextflow/results/gcloud_prod` and the full protein-coding gene
list. A production-style gcloud run is:

```bash
make gcloud_step1
make gcloud_step2
make nextflow_gcloud_prod_step3
```

For a gcloud smoke test, override the Make variables used by Step 1/2 so they
write to the same root that the `gcloud_test` Nextflow profile reads:

```bash
GCLOUD_OUTDIR=gs://vidra-2-0/nextflow/results/gcloud_test \
GCLOUD_GENES=gs://vidra-2-0/nextflow/testdata/genes_ensembl_small.csv \
make gcloud_step1

GCLOUD_OUTDIR=gs://vidra-2-0/nextflow/results/gcloud_test \
make gcloud_step2

make nextflow_gcloud_test_step3
```

Copy gcloud production Bayesian results that are not already present locally:

```bash
make sync_gcloud_prod_results
```

Step 3 commands use `-resume`. Nextflow therefore skips cached per-gene tasks
when the inputs, script, container, and command are unchanged.

## Step 1: Prepare Analysis Input

Script: `tools/prepare_analysis_input.py`

Step 1 runs as one Spark job and computes the intermediate source tables in
memory. The production launchers do not write or reread transient parquet
sections. Persisted outputs are only the public handoff outputs listed above.

Step 1 reads these input families:

- Open Targets colocalisation.
- Open Targets release `credible_set`, `study`, and `variant` in fallback mode.
- Optional `lead_variant_effect` in lead mode, replacing the release
  `credible_set + study + variant` path for GWAS/QTL/coding statistics.
- Open Targets `evidence_gene_burden`, always used as the release burden
  baseline.
- Optional raw AZ PheWAS burden parquet, unioned with release burden when both
  binary and quantitative folders are provided.
- Open Targets `evidence_eva` for ClinVar/EVA evidence.
- AZ rare-variant CSV plus AZ phenotype and gene mapping files.

### Step 1 Optional Modes

`lead_variant_effect` mode is enabled by `--lead_variant_effect_dir`. It changes
GWAS/QTL/coding statistics by using rescaled lead-variant effects. It also
filters lead rows to `p <= 5e-8`, removes rows with
`majorLdPopulationAf.alleleFrequency == 0`, and keeps resolved beta in
`[-3, 3]`.

Raw AZ burden mode is enabled only when both `--az_burden_binary_dir` and
`--az_burden_quantitative_dir` are provided. It also requires
`--target_data_dir` for HGNC to Ensembl mapping. If the raw burden paths are not
provided, Step 1 logs that it is using release burden only.

AZ rare-variant matching always uses the release `variant.variantId` universe as
a diagnostic match index, even when `lead_variant_effect` is enabled. The match
summary is logged, but unmatched AZ rows are not dropped by that diagnostic.

### Step 1 Analysis Schema

- `variant`, `as_gene`, `as_disease`: downstream analysis key.
- `GsourceLab`: source label, where `0` is common coloc, `1` is AZ rare variant,
  `2` is ClinVar/EVA, and `3` is coding GWAS.
- `GqtlLab`: QTL label, where `0` is eQTL, `1` is pQTL, and `2` is not a QTL
  branch.
- `yc`, `ycse`: GWAS-like effect and SE. Common and coding rows use beta/SE. AZ
  rare variants use `log(odds ratio)` and an SE proxy from the OR confidence
  interval. ClinVar rows use `yc=0` because ClinVar is not an effect-size source.
- `xc`, `xcse`: QTL effect and SE. Only common coloc rows carry real QTL values.
- `bO`, `bOse`: gene/disease burden effect and SE joined onto every variant row
  for the same `(as_gene, as_disease)`.
- `as_clinicalSignificance`: ClinVar pathogenicity-like score on `[0, 1]`.

### Step 1 Deduplication

- Common coloc rows are deduplicated by
  `(variant, as_gene, as_disease, GsourceLab, GqtlLab)`, keeping lowest GWAS
  p-value, then lowest QTL p-value, then stable locus IDs.
- Coding GWAS rows use the same source-aware key, keeping lowest GWAS p-value.
- ClinVar rows use the same source-aware key, keeping highest EVA `score`, then
  highest VIDRA clinical-significance rank.
- AZ rare-variant rows use the same source-aware key, keeping lowest `p-value`,
  then highest odds ratio.
- Burden rows are deduplicated by `(targetId, disease)` after release and
  optional raw AZ burden are unioned. Lowest p-value wins; exact ties prefer
  release burden.
- Rows from different source classes are intentionally not collapsed together,
  even when they share variant, gene, and disease.

## Step 2: Variant Annotation

Script: `tools/annotate_variants_cli.py`

Step 2 annotates the distinct coding-GWAS and ClinVar variants needed by Step 3.
It reads:

- `GWAS_coding_variants_from_cs.parquet`
- `clinvar_variants.parquet`

It writes:

- `variant_annotations/annotations.parquet`

Local Step 2 mounts `LOCAL_VEP_RESOURCES_DIR`. Gcloud Step 2 creates a GCE VM,
installs Docker and the Google Cloud CLI, pulls the annotation image, runs the
annotation container, writes progress/status/log files under
`pipeline_info/annotation_vm/<timestamp>/`, and deletes the VM when
`ANNOTATION_VM_AUTO_DELETE=true`.

Required plugins are controlled by Make variables:

- local default: `LOCAL_REQUIRED_VEP_PLUGINS=AlphaMissense`
- gcloud default: `GCLOUD_REQUIRED_VEP_PLUGINS=AlphaMissense,CADD,REVEL`

Plugins that are not required are warning-only and fall back to neutral
defaults. Local `SAVE_INTERMEDIATE_FILES=true` keeps the temporary VEP JSON side
files under `results/.../intermediate/annotate_variants`; otherwise they are
ephemeral.

## Step 3: Bayesian Analysis

Script: `tools/run_bayesian_analysis.py`

Workflow: `main.nf`

Step 3 is the only Nextflow workflow. It requires:

- `bucket_uri`
- `analysis_ready_dir`
- `annotations_dir`

It optionally accepts `genes`. If `genes` is absent, Step 3 runs all discovered
`as_gene=*` partitions. If `genes` is present, it runs the intersection of that
list and discovered partitions.

The shared annotation file is always:

```text
${bucket_uri}/${annotations_dir}/annotations.parquet
```

Each gene partition becomes one `RUN_BAYESIAN_ANALYSIS_PER_GENE` task. Local
profiles run these tasks locally. Gcloud profiles submit them to Google Batch,
with one Batch task per gene.

## Reading `vidra_results`

The final Step 3 output stores evidence labels in `source` and `qtl`.

`source` values:

- `0`: common GWAS / coloc
- `1`: AZ rare variant
- `2`: ClinVar / EVA
- `3`: coding GWAS

`qtl` values:

- `0`: eQTL-supported common variant
- `1`: pQTL-supported common variant
- `2`: not a QTL branch / not applicable

Common source/qtl combinations:

- `(0, 0)`: common GWAS + eQTL
- `(0, 1)`: common GWAS + pQTL
- `(1, 2)`: AZ rare variant
- `(2, 2)`: ClinVar / EVA
- `(3, 2)`: coding GWAS

In `vidra_results`, `source` and `qtl` can be scalar strings for single-variant
models, for example `"0"`, or sorted list strings for multi-variant models, for
example `"[0, 2]"`. List values mean the same `(gene, disease)` model contains
multiple evidence types.

## Bayesian Model Notes

The multi-variant model is `stan_models/VIDRA.stan`. The single-variant model is
`stan_models/VIDRA_single_variant.stan`.

### Burden Evidence Is Scalar

Step 1 joins the selected burden `bO` and `bOse` onto every variant row for the
same `(gene, disease)` pair. That does not mean there are multiple independent
burden tests. Step 3 collapses repeated burden values to one scalar observation
before calling Stan.

This avoids counting the same burden evidence once per variant. If Step 3 sees
conflicting non-zero burden values inside one `(gene, disease)` model, it emits
an error row for that model instead of silently choosing one value.

### Shared Intercept With ClinVar Scale Adaptation

The multi-variant model keeps one shared intercept, `intercept_random`, on the
log-OR scale. Burden evidence informs it once, AZ rare variants use it directly,
and coding GWAS also uses it directly. With no burden evidence, the intercept is
still free through a loose `normal(0, 10)` prior, so coding effects are not
forced through zero.

```stan
intercept_random ~ normal(0, 10);
if (has_burden == 1) {
  bO ~ normal(intercept_random, bOse);
}
```

ClinVar is different because `as_clinicalSignificance` is bounded on `[0, 1]`,
not on the log-OR scale. The same shared intercept is adapted with a sigmoid
inside the ClinVar branch:

```stan
disease_prior[n] ~ student_t(
  nu,
  inv_logit(intercept_random) + slope_random[5] * protein_prior[n],
  0.3
)
```

There is no extra ClinVar-specific prior on `slope_random[5]`; VIDRA_5 uses the
same generic `slope_random ~ normal(0, 5)` prior as the updated VIDRA_2 model.

### Source Branches

- `GsourceLab == 0`: common GWAS/QTL evidence; eQTL and pQTL use latent QTL
  effect `xcest` and source-specific slopes.
- `GsourceLab == 3`: coding GWAS; uses `intercept_random`, coding slope, and
  protein prior.
- `GsourceLab == 1`: AZ rare variants; uses `intercept_random`, AZ slope, and
  protein prior.
- `GsourceLab == 2`: ClinVar; uses `inv_logit(intercept_random)` as baseline and
  adds the ClinVar slope and protein prior on the `[0, 1]` disease-prior scale.

The global `slope` remains a hierarchical summary across source-specific
`slope_random[...]` terms, but only sources present in a gene/disease block are
tied back to the global slope.

### Retained Step 3 Filters And Scale Choices

Both the current VIDRA_2 snapshot and VIDRA_5 drop single-variant
coding-GWAS-only groups in Step 3:

```python
len(group["variant"]) == 1 and GsourceLab == 3
```

The coding-GWAS and AZ branches still use the legacy Student-t scale expression
based on `abs(yOR[n] / protein_prior[n])`. If future runs show numerical
instability for variants with near-zero `protein_prior`, revisit that scale term
separately.

`VIDRA_single_variant.stan` is unchanged relative to the current VIDRA_2
snapshot and does not consume `bO` or `bOse`.

## Remaining VIDRA_2 vs VIDRA_5 Analysis Differences

This section tracks analysis-relevant differences that remain after the VIDRA_5
refactor. It excludes cloud launchers, Docker, Nextflow, and paths unless they
change rows, columns, model inputs, or Bayesian behavior.

Aligned behavior:

- Public Step 1 output column meanings are aligned.
- Source labels and QTL labels are aligned.
- Multi-variant Stan source branches are aligned, except VIDRA_5 passes burden
  as a validated scalar instead of a repeated vector.
- `VIDRA_single_variant.stan` is identical to the current VIDRA_2 snapshot.
- Both versions still drop single-variant coding-GWAS-only Step 3 groups.

Remaining differences (needs to be updated):

| Area | VIDRA_2 behavior | VIDRA_5 behavior | Expected impact |
| --- | --- | --- | --- |
| Burden source policy | Uses raw AZ PheWAS burden as the burden source. | Uses release `evidence_gene_burden`; current Make defaults also union optional raw AZ burden. | Result-changing when release and raw AZ differ or when raw AZ is disabled. |
| Burden dedup after mapping | Chooses lowest raw `pValue`, then after obsolete-EFO remap deduplicates merged `(gene, disease)` rows by largest `abs(bO)`. | Keeps all phenotype-to-EFO mappings, then keeps one `(gene, disease)` row by best p-value; exact ties prefer release burden. | Intentional. VIDRA_5 consistently favors statistical significance rather than switching to effect magnitude after disease remap. |
| Burden Stan input | Passes repeated `bO/bOse` vectors, but updated Stan uses `bO[1]/bOse[1]`. | Collapses to scalar `bO/bOse` before Stan and validates repeated values. | Equivalent when repeated values match; safer when they do not. |
| ClinVar duplicates | Takes first clinical-significance label per row, with older Step 3 dedup able to select by row order. | Takes first label per row, then Step 1 keeps highest EVA `score`, then highest VIDRA clinical-significance rank. | Result-changing only for duplicated ClinVar rows. |
| AZ rare-variant duplicates | Uses broad `dropDuplicates`; older Step 3 dedup can later prefer higher absolute `yc`. | Step 1 keeps lowest `p-value`, then highest odds ratio, then stable tie-breakers. | Result-changing for duplicated AZ rows; p-value is now primary. |
| Step 3 duplicate handling | Performs legacy dedup inside Step 3. | Legacy Step 3 dedup is commented out because Step 1 performs source-aware dedup. | Intentional. Duplicate policy is now explicit in Step 1. |
| Study-to-disease mapping | Uses the study phenotype mapping available in VIDRA_2, including its EFO remap path. | Prefers `traitFromSourceMappedIds`, with `diseaseIds` fallback; same rule in release and lead mode. | Result-changing if study mapping fields disagree. |
| Obsolete EFO remapping | Builds and applies a global obsolete-to-current EFO map. | Relies mostly on already mapped release IDs and curated AZ mappings; no separate global remap pass. | Potential edge case if inputs contain obsolete EFO IDs. |
| Optional `lead_variant_effect` mode | Not present. | Optional analysis mode using rescaled lead statistics, `p <= 5e-8`, non-zero allele frequency, and beta in `[-3, 3]`. | Deliberately result-changing; not VIDRA_2-equivalent. |
| AZ release variant check | Keeps normalized mapped AZ rows. | Logs whether normalized AZ IDs exist in release `variant`, but does not drop unmatched rows. | Diagnostic only. |

For the closest VIDRA_2-like run, disable `lead_variant_effect` and enable raw
AZ burden. Even then, VIDRA_5 still keeps its explicit ClinVar/AZ dedup rules,
p-value-based burden conflict resolution, and scalar burden validation before
Stan.

## Local VEP Data

Local Step 2 expects VEP resources under `LOCAL_VEP_RESOURCES_DIR`. Set that
path in `conf/make.local.mk`, for example:

```make
LOCAL_VEP_RESOURCES_DIR := /path/to/data/VEP
```

Expected minimum layout:

```text
data/VEP/
├── homo_sapiens/
│   └── 115_GRCh38/
└── alphamissense/
    ├── AlphaMissense_hg38.tsv.gz
    └── AlphaMissense_hg38.tsv.gz.tbi
```

Download and unpack the Homo sapiens VEP 115 cache:

```bash
curl -O https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz
tar xzf homo_sapiens_vep_115_GRCh38.tar.gz
```

Download AlphaMissense and generate its tabix index:

```bash
curl -L -o AlphaMissense_hg38.tsv.gz \
  https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg38.tsv.gz
```

On macOS, install `tabix` with:

```bash
brew install htslib
```

## Acknowledgements
This work was carried out at the Wellcome Genome Campus in collaboration with Open Targets, and during an exchange funded by an EMBO Scientific Exchange Grant.
