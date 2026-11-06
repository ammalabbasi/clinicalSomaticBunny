# SomaticBunny

<img src="workflow_logo/SomaticBunny_Logo.png" alt="SomaticBunny Logo" width="300" align="left"/>

**SomaticBunny** is an ensemble pipeline for somatic variant calling that integrates multiple tools to provide comprehensive variant detection in tumor-normal paired samples.

<br clear="left"/>

## Clinical SomaticBunny at a glance

This branch adapts SomaticBunny to clinical panel sequencing — vendor-processed BAMs or FASTQs, with or without a matched normal. Purple steps are new on this branch, blue steps are original SomaticBunny reused or adapted, and grey steps are original steps not used for clinical data. Every validation number in this README comes from simulated reads on public GRCh38; no patient data was used.

![Clinical SomaticBunny workflow](docs/clinical_somaticbunny_workflow.svg)

A slide-ready PNG and PDF of the same figure live in [`docs/`](docs/).


## Table of Contents
- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
  - [Basic Usage](#basic-usage)
  - [Choosing Where Jobs Run](#choosing-where-jobs-run)
  - [Custom Reference Genome](#custom-reference-genome)
  - [Clinical BAMs: panels, tumour-only samples, vendor references](#clinical-bams-panels-tumour-only-samples-vendor-references)
  - [Advanced Options](#advanced-options)
- [Output](#output)
- [Tool Versions](#tool-versions)
- [Benchmark](#benchmark)

## Overview

The SomaticBunny pipeline integrates multiple variant calling tools to improve detection accuracy for:
- Single Nucleotide Variants (SNVs)
- Insertions and Deletions (INDELs)
- Structural Variants (SVs)
- Copy Number Variants (CNVs)

This pipeline is designed to work with paired tumor-normal whole-genome or whole-exome sequencing data.

<p align="center">
  <img src="https://github.com/AlexandrovLab/SomaticBunny/blob/main/workflow_logo/SomaticBunny.png" alt="SomaticBunny Pipeline Workflow"/>
</p>


## Prerequisites

- Conda/Mamba for environment management
- A SLURM cluster is recommended for WGS data (see [Choosing Where Jobs Run](#choosing-where-jobs-run))

## Installation

1. Clone this repository

2. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) as a conda environment

3. Make sure you have the following required files in your working directory:
   - `main.nf`
   - `nextflow.config`
   - `conf/base.config`
   - `sample.csv`

4. Prepare your `sample.csv` file according to your analysis requirements:

   <details>
   <summary><b>For starting from fastq (click to expand)</b></summary>
   
   ```csv
   patient,sample,status,fastq_1,fastq_2
   UCSD101,A,normal,/PATH/TO/ERR5285401_1.fastq.gz,/PATH/TO/ERR5285401_2.fastq.gz
   UCSD101,A,tumor,/PATH/TO/ERR5285402_1.fastq.gz,/PATH/TO/ERR5285402_2.fastq.gz
   UCSD101,B,tumor,/PATH/TO/ERR5285403_1.fastq.gz,/PATH/TO/ERR5285403_2.fastq.gz
   ```
   </details>
   <details>
   <summary><b>For starting from bams (click to expand)</b></summary>
   
   ```csv
   patient,sample,status,bam,bai,sex
   UCSD101,A,normal,/PATH/TO/UCSD101_A_normal.bam,/PATH/TO/UCSD101_A_normal.bam.bai,XY
   UCSD101,A,tumor,/PATH/TO/UCSD101_A_tumor.bam,/PATH/TO/UCSD101_A_tumor.bam.bai,XY
   UCSD101,B,tumor,/PATH/TO/UCSD101_B_tumor.bam,/PATH/TO/UCSD101_B_tumor.bam.bai,XY
   ```
   </details>

   <details>
   <summary><b>For ASCAT analysis (click to expand)</b></summary>
   
   ```csv
   patient,sample,status,fastq_1,fastq_2,sex
   UCSD101,A,normal,/PATH/TO/ERR5285401_1.fastq.gz,/PATH/TO/ERR5285401_2.fastq.gz,XY
   UCSD101,A,tumor,/PATH/TO/ERR5285402_1.fastq.gz,/PATH/TO/ERR5285402_2.fastq.gz,XY
   UCSD101,B,tumor,/PATH/TO/ERR5285403_1.fastq.gz,/PATH/TO/ERR5285403_2.fastq.gz,XY
   ```
   </details>

  > [!NOTE]
  > **One FASTQ pair per sample.** `sample.csv` takes exactly one `fastq_1`/`fastq_2` pair per
  > sample. If a sample was sequenced on several lanes or runs, concatenate the lanes into one
  > pair before running the pipeline. Each BAM therefore carries a single read group:
  > `ID`, `SM` and `LB` are `<patient>_<sample>_<status>`, `PL` is `ILLUMINA`, and `PU` is the
  > flowcell and lane of the first read (or the sample name if the read names are not in the
  > Illumina format). Because merged lanes share one read group, base quality score
  > recalibration models one error profile per sample rather than per lane.
  > bwa-mem2 runs with `-Y` (soft-clipped supplementary alignments) and a fixed batch size
  > `-K 160000000`, so alignments do not depend on the number of CPUs assigned.



## Configuration

For memory management, modify the following parameters if needed:

- Set `$params.mkdup_temp_dir` in `main.nf` to your ideal location (Default: `$projectDir/mkdup_tmp`)
- Set `$workDir` in `nextflow.config` to your ideal location (Default: `./work`)

```nextflow
// In main.nf
params.mkdup_temp_dir = "/path/to/ideal/mkdup_tmp/folder/"

// In nextflow.config
workDir = "/path/to/ideal/intermediate/work/folder/"
```

## Usage

### Basic Usage

1. Request an interactive compute node (TSCC example; on another cluster, use your own partition and account):

```bash
srun -N 1 -n 1 -c 16 --mem 250G -t 80:00:00 -p platinum -q hcp-ddp302 -A ddp302 --pty bash
```

2. Activate your Nextflow environment:

```bash
conda activate env_nf
```

3. Set TSCC temporary directory and Nextflow tmp directory, and memory:

```bash
export TMPDIR=/path/to/restricted/folder/

export NXF_OPTS="-Djava.io.tmpdir=${TMPDIR} -Xms4g -Xmx16g"
```

4. Run the pipeline, choosing where jobs run with `-profile` (see [Choosing Where Jobs Run](#choosing-where-jobs-run)). This example runs on TSCC:

```bash
nextflow run main.nf -profile tscc_slurm,conda --type [genome|exome] --genome [GRCh38|GRCh37|mm39|RN7] --first_step [mapping|markdup|recalibration|variant_calling] --tool [ascat,manta,...] --database_path [/path/to/SomaticBunny_database] --ref [/path/to/ref/] --ref_fai [/path/to/ref_fai] --ref_dict [/path/to/ref_dict] --bed [/path/to/bed] --bed_tbi [/path/to/bed_tbi] -profile [conda, ...]
```

#### Pipeline Steps and Tools

The SomaticBunny pipeline workflow is divided into multiple steps that can be run individually or in sequence using the `--first_step` parameter.

##### Available first_step

| First_step | Description | Required Input | 
|------|-------------| ---------------|
| `mapping` | This will start from performing alignment of FASTQ files to the reference genome | fastq_1 and fastq_2 |
| `markdup` | This will start from mark duplicate reads using the aligned BAM files | raw bam and bai |
| `recalibration` | This will start from performing base quality score recalibration | mark duplicated bam and bai |
| `variant_calling` | This will start from executing variant calling tools | recalibrated bam and bai |

##### Available reference genomes
| Genomes | 
|---------|
| GRCh38  | 
| GRCh37  | 
| mm39  | 
| RN7  | 

> [!NOTE]
> **`--type genome` vs `--type exome`.** The SNV/INDEL callers run on the same regions for both
> data types: the callable-regions BED (`--bed`, by default the primary chromosomes of the
> selected genome). `--type exome` does not restrict calling to exome targets. It only switches
> Strelka and Manta to `--exome` mode, MuSE to `-E`, base quality score recalibration to the
> exome interval list, and mosdepth and ASCAT to their exome resources.

##### Available Tools

Use the `--tool` parameter to specify which variant callers to run. You can select multiple tools by separating them with commas.

| Category | Available Tools | Description |
|----------|----------------|-------------|
| SV Callers | `manta` | Structural variant and indel caller |
| | `delly` | Integrated structural variant detection |
| CNV Callers | `ascat` | Allele-specific copy number analysis (requires sex information) |
| | `cnvkit` | Copy number variation detection from targeted DNA sequencing |

Example of running multiple tools:
```bash
# Run Manta for SVs and ASCAT for CNVs
nextflow run main.nf -profile tscc_slurm,conda --genome GRCh38 --type genome --first_step mapping --tool manta,ascat
```

### Choosing Where Jobs Run

Use `-profile` (one dash, it is a Nextflow option) to choose the scheduler, combined with `conda` for the software. Profiles are separated by a comma, with no space.

| Profile | Runs on |
|---------|---------|
| `tscc_slurm,conda` | UCSD TSCC: each job is submitted to SLURM with the lab's partition, QoS and account |
| `slurm,conda` | Any SLURM cluster: each job is submitted to SLURM with the cluster's defaults, or with your own options (see below) |
| `conda` | The current machine only (no scheduler) |

> [!WARNING]
> Without `tscc_slurm` or `slurm`, all jobs run on the machine where Nextflow is started, which is not practical for WGS data.

#### Running on other SLURM clusters

If your cluster needs a partition, account or QoS, put them in a small config file of your own and pass it with `-c`. You don't need to edit SomaticBunny's `nextflow.config`, so your settings are kept when you update the pipeline.

1. Create `my_cluster.config` (any name works):

```groovy
// my_cluster.config
process.clusterOptions = '-p <partition> -A <account>'
```

2. Run with the `slurm` profile and your config:

```bash
nextflow run main.nf -profile slurm,conda -c my_cluster.config --type genome ...
```

`process.clusterOptions` is passed to every job as extra `sbatch` options, so anything `sbatch` accepts works, e.g. `'--partition=compute --account=mylab --qos=normal'`.

If your nodes are smaller than the default requests (up to 16 CPUs and 250 GB for some steps), cap what any job can request:

```bash
nextflow run main.nf -profile slurm,conda -c my_cluster.config --max_cpus 32 --max_memory 192.GB --max_time 72.h ...
```

#### Software

SomaticBunny installs its tools with Conda (`-profile conda`, or `-profile mamba` for faster environment creation). Container engines such as Docker and Singularity are not supported.

### Custom Reference Genome

When starting the pipeline from BAM files (`--first_step markdup`, `recalibration`, or `variant_calling`), you **must** provide your own reference genome files. This is because BAM files are already aligned to a specific reference, and the pipeline needs the matching reference for downstream analysis.

#### Required parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--ref` | Reference FASTA file | `/path/to/hg38_genome.fa` |
| `--ref_fai` | FASTA index file | `/path/to/hg38_genome.fa.fai` |
| `--ref_dict` | Sequence dictionary | `/path/to/hg38_genome.dict` |
| `--bed` | Callable regions BED file (bgzipped) | `/path/to/callable_regions.bed.gz` |
| `--bed_tbi` | Tabix index for the BED file | `/path/to/callable_regions.bed.gz.tbi` |

All five files must be provided together. The pipeline will exit with an error if any are missing.

#### Example

```bash
nextflow run main.nf \
  -profile tscc_slurm,conda \
  --type genome \
  --genome GRCh38 \
  --first_step variant_calling \
  --sample sample.csv \
  --database_path /path/to/SomaticBunny_database \
  --ref /path/to/hg38_genome.fa \
  --ref_fai /path/to/hg38_genome.fa.fai \
  --ref_dict /path/to/hg38_genome.dict \
  --bed /path/to/hg38_genome_callable_regions.bed.gz \
  --bed_tbi /path/to/hg38_genome_callable_regions.bed.gz.tbi
```

When a custom reference is provided, the pipeline will automatically generate the scattered interval lists required by GATK tools (Mutect2, BaseRecalibrator, etc.) from your reference. When starting from `mapping`, these parameters are optional and will fall back to the pipeline's built-in defaults for the selected `--genome`.

> [!IMPORTANT]
> **Your custom reference genome must use the `chr` prefix for contig names** (e.g., `chr1`, `chr2`, ..., `chrX`, `chrY`). The pipeline's built-in database files (dbSNP, gnomAD, Panel of Normals, known indels, etc.) use `chr`-prefixed contig names. A reference without the `chr` prefix (e.g., `1`, `2`, ..., `X`, `Y`) will cause tools such as Mutect2, MuSE, Strelka, SAGE, and BQSR to fail due to contig name mismatches.
>
> Different GRCh38 builds that include additional or fewer contigs (e.g., alt contigs, decoys, HLA) are supported as long as the standard chromosomes use the `chr` prefix.
>
> BAMs aligned to a no-`chr` reference (b37 style: `1`, `2`, `MT`) can still be used: pass a `chr`-named FASTA with the same sequences as `--ref`, and `--contig_rename auto` (the default) renames the BAM's contigs on the fly. See below.

### Clinical BAMs: panels, tumour-only samples, vendor references

Vendor deliveries (e.g. Tempus, CARIS) are usually processed BAMs from targeted panels, often without a matched normal, and often aligned to b37. Start them at `--first_step variant_calling`. Clinical FASTQs start at `--first_step mapping` (sheet columns `patient,sample,status,fastq_1,fastq_2`, one row per sample, so concatenate lanes first); panel mode and tumour-only patients work there too.

| Option | Default | What it does |
|--------|---------|--------------|
| `--type panel` | | Targeted mode like `exome` (Strelka/Manta `--exome`, MuSE `-E`). Requires `--targets`. Refuses `ascat`/`delly`, which need genome-wide data. A Conpair failure is retried, then skipped for that pair instead of stopping the run |
| `--targets <bed>` | | Capture BED (plain or `.gz`), same contig names as `--ref`. Padded, merged and used for: Mutect2/GetPileupSummaries shards (instead of whole chromosomes), Strelka `--callRegions`, mosdepth and CNVkit targets, BQSR intervals. Also accepted with `--type exome` |
| `--target_padding` | `100` | bp added on each side of every target |
| `--contig_rename` | `auto` | Rename BAM contigs between b37 and `chr` style (`1`↔`chr1`, `MT`↔`chrM`) when the other name exists in `--ref` with the same length. `none` turns this off. Either way the run stops early, with a readable message, if the BAM's primary chromosomes do not match `--ref` |
| `--cnvkit_reference` | `auto` | `pooled`: CNVkit reference built from the run's normals. `flat`: no normals needed; every bin gets the same expected coverage, with GC and repeat corrections but no per-bait capture correction, so large CNAs stay visible and small ones are unreliable. `auto` uses pooled when the sheet has normals and flat when it has none |
| `--outdir <dir>` | `<pipeline>/RESULTS` | Where `RESULTS`, `REPORT` and `pipeline.log` go. Use a protected location for patient data, so results are not written inside the code checkout |

**Tumour-only samples.** A patient with a tumour row and no normal row is no longer dropped. From any entry point (`mapping`, `markdup`, `recalibration` or `variant_calling`) with `--genome GRCh38` or `GRCh37`:

- Mutect2 runs in tumour-only mode (panel of normals + gnomAD germline resource), and SAGE runs without `-reference`.
- CNVkit calls copy number against the pooled normals of the run (at least one paired normal is needed).
- Strelka2, MuSE, Conpair, ASCAT, Delly and Manta need a matched normal, so they are skipped for these patients.
- `POST/tumor_only/<patient>_<sample>/` holds the Mutect2 ∩ SAGE consensus. Sites with population AF > 0.001 (Mutect2 `POPAF` < 3) are soft-filtered as `germline_popaf`, and `*.PASS.vcf.gz` keeps the rest.
- `REPORT/tumor_only_patients.txt` lists them.

Tumour-only calls still contain rare germline variants that are absent from gnomAD and the panel of normals, so they are not equivalent to paired calls. Do not pool them for TMB or signatures without checking the false-positive rate on your assay. For example, run some paired tumours again as tumour-only and compare.

On mm39 and RN7, which have no panel of normals or SAGE resources, a sheet with tumour-only patients stops with an error.

<!-- >>> mutect2_pon -->
**Mutect2 panel of normals from your own normals.** The database PON (`MuTect2.PON.5210` on GRCh38) is generic WGS/WES. Recurrent artefacts specific to one assay (capture, chemistry, vendor pipeline) are not in it, so they show up as low-VAF "somatic" calls, most of all in tumour-only samples. A PON built from normals sequenced on the same assay removes them.

| Option | Default | What it does |
|--------|---------|--------------|
| `--mutect2_cohort_pon build` | off | Build a PON from every normal BAM in the sheet (paired normals and normal-only rows), then call this run's tumours (paired and tumour-only) against it instead of the database PON. Each normal goes through Mutect2 tumour-only with `--max-mnp-distance 0` over the same interval shards as the tumours (target shards with `--targets`). Then GenomicsDBImport over the padded targets (or the shards) and CreateSomaticPanelOfNormals with the gnomAD germline resource. A site enters the PON when at least 2 normals carry it |
| `--mutect2_cohort_pon_name` | `cohort` | Output name: `Mutect2_PON/<name>.pon.vcf.gz(.tbi)` and `<name>.pon.normals.txt` |
| `--mutect2_pon <vcf.gz>` | database PON | Use a prebuilt PON (for example one published by an earlier `build` run) for all Mutect2 calling. Its `.tbi` (or `.idx` for `.vcf`) must be next to it. Cannot be combined with `build` |

Fewer than 2 normals stops the run, and fewer than 10 gives a warning: with only a few normals the PON catches only the most recurrent artefacts. Build one PON per assay, from normals sequenced on that assay. For a cohort that is mostly tumour-only, add the assay's normals as normal-only rows, or build the PON once and reuse it with `--mutect2_pon`. The PON log line and `pipeline.log` say which PON was used.
<!-- <<< mutect2_pon -->

**Sample-sheet checks.** The run stops before any job is submitted if:

- a status is not `tumor`/`normal`;
- a `patient,sample,status` row is duplicated;
- a patient has more than one normal;
- a row has the wrong number of fields;
- an ID contains `tumor`/`normal` as an underscore-separated word.

Example, a Tempus-style panel cohort aligned to b37:

```bash
nextflow run main.nf -profile tscc_slurm,conda \
  --type panel --targets panel_targets.chr.bed --tool cnvkit \
  --genome GRCh37 --first_step variant_calling --sample samples.csv \
  --database_path /path/to/SomaticBunny_database \
  --ref human_g1k_v37.chr.fasta --ref_fai human_g1k_v37.chr.fasta.fai --ref_dict human_g1k_v37.chr.dict \
  --bed callable_regions.bed.gz --bed_tbi callable_regions.bed.gz.tbi \
  --outdir /protected/path/RESULTS
```

Run one assay per run: the target BED and the CNVkit reference are per assay.

<!-- >>> tmb -->
#### Tumour mutational burden (`--tool tmb`)

TMB = eligible somatic variants / callable Mb, per tumour, for pairs and tumour-only patients from every entry point (`--genome GRCh38|GRCh37`).

- **Callable territory**: bases of `--tmb_regions` (default: the cleaned `--targets`; pass a coding-exon BED to count coding territory only) where the tumour has depth ≥ `--tmb_min_depth` and, for pairs, the matched normal has depth ≥ `--tmb_min_normal_depth`. Depth comes from `mosdepth --quantize` with reads of MAPQ ≥ `--tmb_min_mapq` that are not duplicate, secondary, supplementary or QC-failed.
- **Eligible variants**: PASS SNVs, MNVs and indels (one per ALT allele) whose VCF `POS` (the anchor base for indels) is callable and whose tumour VAF, recomputed from the tumour BAM (MAPQ ≥ `--tmb_min_mapq`, base quality ≥ `--tmb_min_baseq`), is ≥ `--tmb_min_vaf`.
  - Pairs count the POST consensus of `--tmb_scenario` (default 4). All four standard VCFs (scenarios 3–6) are written by `MERGE_HEADER`, but scenario 3 does not gate on PASS. Scenario 4 is the 2-of-4 SNV / 2-of-3 indel consensus with the PASS gate. Scenarios 5 and 6 also prefilter each caller on PASS and run the DKFZ bias filter, and 6 adds TLOD ≥ 10 / SomaticEVS ≥ 13. All of them keep the same records as `POST/standard_vcf/*_scenario_N_standard_PASSed.vcf`.
  - Tumour-only patients count `POST/tumor_only/*_tumor_only_consensus.PASS.vcf.gz`. Optionally, `--tmb_germline_calls <dir>` holding PureCN `<patient>_<sample>_variants.csv` files drops variants with `POSTERIOR.SOMATIC` < `--tmb_min_somatic_posterior` (0.5). Variants that PureCN did not assess are kept.
  - `--tmb_coding_only` counts only coding consequences, and only when the VCF carries a VEP `CSQ` or SnpEff `ANN` annotation. The standard consensus VCFs are not annotated, so for now the option is ignored with a warning.

| Option | Default |
|--------|---------|
| `--tmb_regions` | cleaned `--targets` |
| `--tmb_min_depth` / `--tmb_min_normal_depth` | `50` / `20` |
| `--tmb_min_vaf` | `0.05` |
| `--tmb_min_mapq` / `--tmb_min_baseq` | `20` / `20` |
| `--tmb_scenario` | `4` (4, 5 or 6) |
| `--tmb_coding_only` | `false` |
| `--tmb_germline_calls`, `--tmb_min_somatic_posterior` | none, `0.5` |

Outputs go to `TMB/`:

- `<patient>_<sample>.tmb.tsv`: `n_variants`, `n_snv`, `n_indel`, `callable_bp`, `callable_Mb`, `TMB`, the thresholds, `mode`, `source`, exclusion counts and `warning`.
- `<patient>_<sample>.tmb.variants.tsv`: every consensus record with its tumour VAF and the reason it was counted or excluded.
- `<patient>_<sample>.tmb.callable.bed`: the callable territory.
- `tmb_cohort.tsv`: one row per tumour.

When callable territory is under 1 Mb, the `warning` column and the run log say so. Panel TMB from a few hundred kb is noisy, and a 1-Mb panel with 5 mutations has a Poisson 95% interval of about 1.6–11.7 mut/Mb. Panel TMB is not comparable to WES TMB unless it is calibrated on the assay. Tumour-only TMB is inflated by rare germline variants (see above), so it should not be pooled with paired TMB.
<!-- <<< tmb -->
<!-- >>> msi -->
#### Microsatellite instability (`--tool msi`)

[MSIsensor-pro](https://github.com/xjtu-omics/msisensor-pro) 1.3.0 (bioconda, `yml/msisensor_pro.yml`) scores MSI as the percentage of evaluated microsatellites that are unstable.

- **Sites** (`MSI_SCAN`, once per run): `msisensor-pro scan` over the padded `--targets` (panel/exome) or the whole `--ref` (genome). Published to `MSI/sites/msi_sites.tsv`. `--msi_sites <file>` reuses a list and skips the scan.
- **Paired** patients: `msisensor-pro msi` compares the tumour with its own normal.
- **Tumour-only** patients: `msisensor-pro pro` compares the tumour with a baseline of normals from the same assay. The baseline is built from every normal in the sheet, paired normals plus normal-only pool rows, and published to `MSI/baseline/msi_baseline.tsv` (normals listed in `msi_baseline.normals.txt`). Pass it to a later run of the same assay with `--msi_baseline <file>`. If the sheet has tumour-only patients but no normals and no `--msi_baseline` is given, the run stops before any job is submitted. With fewer than 20 normals, the run warns: the authors recommend at least 20 normals from the same sequencer.

| Option | Default | What it does |
|--------|---------|--------------|
| `--msi_threshold_paired` | `10` | MSI-H when score (% unstable sites) ≥ this, paired. MSK-IMPACT's cut-off for MSIsensor on a panel (Middha et al. 2017, JCO Precis Oncol). WES studies often use 3.5 |
| `--msi_threshold_tumor_only` | `20` | MSI-H cut-off for tumour-only scores. MSIsensor-pro publishes **no** cut-off and says to calibrate one per panel against known MSI status, so treat 20 as a placeholder |
| `--msi_min_sites` | `50` | Fewer sites evaluated → `qc = LOW_SITES` (the MSIsensor-pro paper shows it works with ≥ 50 sites) |
| `--msi_min_coverage` | `20` | Read depth a site needs to be evaluated (`-c`; MSIsensor-pro suggests 20 for WES, 15 for WGS) |
| `--msi_baseline_min_normals` | `10` | Baseline sites need coverage in this many normals (`baseline -s`). Capped at the number of normals |

Outputs: `MSI/paired/<patient>_<sample>/` and `MSI/tumor_only/<patient>_<sample>/` hold the raw MSIsensor-pro files and `<patient>_<sample>.msi_status.tsv`. `MSI/msi_summary.tsv` has one row per tumour: `patient, sample, mode, sites_evaluated, unstable_sites, msi_score_pct, threshold_pct, status (MSI-H|MSS), qc (PASS|LOW_SITES)`. Paired and tumour-only scores are on different scales, so do not compare them with one cut-off. Check both thresholds on your assay with samples of known MSI status (IHC/PCR) before using the status.
<!-- <<< msi -->
<!-- >>> umi -->
#### UMI consensus instead of MarkDuplicates (`--umi`, FASTQ entry only)

Panels sequenced with unique molecular identifiers (UMIs) can collapse PCR duplicates into one consensus read per molecule. Errors present in only one read of a family are voted out, so low-VAF calls are less often PCR or sequencing errors. With `--umi` and `--first_step mapping`, each sample goes through [fgbio](https://github.com/fulcrumgenomics/fgbio) 3.1.2 (`yml/umi.yml`):

`FastqToBam` (UMI → `RX`) → `bwa-mem2` + `ZipperBams` → `GroupReadsByUmi` → `CallMolecularConsensusReads` (or `CallDuplexConsensusReads`) → re-align consensus (`bwa-mem2` + `ZipperBams`) → `FilterConsensusReads` → BQSR → callers.

MarkDuplicates is skipped because the consensus reads are already one per molecule. Read groups keep `ID:<patient>` and `SM:<patient>_<sample>_<status>`. Paired and tumour-only patients both work. `--umi` is refused at the other entry points, because vendor BAMs are already deduplicated.

How the UMI is delivered differs between vendors. Give exactly one source:

| UMI source | Option | Example |
|---|---|---|
| Inline, at the start of the reads | `--umi_read_structure '<R1> <R2>'` ([fgbio read structures](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures)) | `'3M2S+T 3M2S+T'` = 3 bp UMI, 2 bp skipped, rest template, on both reads |
| Separate UMI/index read | `--umi_fastq_column` + sheet column `fastq_umi` | default structures `'+T +T +M'`; override with three structures |
| In the read name (`@...:ACGTACGT`, BCL Convert style) | `--umi_from_header` | `--umi_read_structure` may still trim, e.g. `'5S+T 5S+T'`, but cannot contain `M` |

| Option | Default | What it does |
|---|---|---|
| `--umi_group_strategy` | `adjacency` (`paired` with `--umi_duplex`) | `GroupReadsByUmi --strategy`: `identity`, `edit`, `adjacency` or `paired` |
| `--umi_edits` | `1` | UMI mismatches allowed within a family |
| `--umi_min_map_q` | `1` | Reads below this MAPQ are not grouped |
| `--umi_min_reads` | `1` | Reads needed for a consensus (call and filter). `1` keeps singletons. Duplex takes up to three values, e.g. `'2 1 1'` |
| `--umi_duplex` | off | Duplex consensus. Needs a UMI on both ends of the molecule (`M` in both read structures, or two UMIs in the UMI read or read name) |
| `--umi_min_base_quality` | `30` | `FilterConsensusReads`: consensus bases below this become `N`. fgbio's own best practice uses 40–45 with `--umi_min_reads ≥ 2`. At `1`, a high value masks every singleton |
| `--umi_max_base_error_rate` / `--umi_max_no_call_fraction` | `0.2` / `0.2` | `FilterConsensusReads` base masking and read removal |

Outputs go to `UMI/` (`<patient>_<sample>_<status>.*`):

- `fgbio_group_family_size_histogram.txt` and `fgbio_group_metrics.txt` (MultiQC's fgbio module reads the histogram)
- `umi_consensus_filter_counts.tsv`: consensus reads in, and kept after filtering
- `fgbio_duplex.*` (duplex only: CollectDuplexSeqMetrics family sizes, duplex yield, UMI counts)
- `*_umi_consensus.bam` (only with `--publish markdup_bam`)

`--cleanup` does not remove the UMI intermediates. Throughput depends on family sizes. Check the histogram before you trust `--umi_min_reads`.

Synthetic check (simulated reads from a public GRCh38 window; panel of 40 × 250 bp; family sizes 1–8; per-read errors at 15 hotspots, 4 % of reads each; five SNVs at 2–3 % VAF per tumour). The same reads were run with and without `--umi`:

- **Families:** fgbio recovered 120,010 families from 120,000 simulated molecules.
- **Paired, Strelka2:** error-hotspot false positives fell from 12/15 to 0/15. Low-VAF SNVs found went from 3/5 to 4/5.
- **Paired, POST consensus:** hotspot false positives fell from 6 to 0. Low-VAF SNVs found went from 0/5 to 4/5.
- **Tumour-only, Mutect2:** low-VAF SNVs found went from 0/5 to 5/5.
<!-- <<< umi -->
<!-- >>> vep -->
#### Variant annotation and MAF (`--tool vep`)

Ensembl VEP (offline cache) and vcf2maf run on the consensus calls after POST, from any entry point (`mapping`, `markdup`, `recalibration`, `variant_calling`), GRCh38 or GRCh37. Environment: `yml/vep.yml` (ensembl-vep 112.0, vcf2maf 1.6.22).

| Option | Default | What it does |
|--------|---------|--------------|
| `--vep_cache <dir>` | required | VEP cache root holding `homo_sapiens/<version>_<assembly>/` (indexed cache). On TSCC: `/tscc/projects/ps-lalexandrov/shared/vep_ref` (has `112_GRCh38`) |
| `--vep_cache_version` | `112` | Must match the ensembl-vep release in `yml/vep.yml` |
| `--vep_assembly` | `--genome` | `GRCh38` or `GRCh37` |
| `--vep_fasta` | `--ref` | FASTA for HGVS and vcf2maf; same contig names as the VCFs |
| `--vep_scenario` | `scenario_5` | Which paired consensus to annotate: `POST/standard_vcf/<patient>_<sample>_<scenario>_standard_PASSed.vcf` |

What is annotated:

- Paired: the standard consensus VCF of `--vep_scenario`. The default `scenario_5` keeps calls that PASS in each contributing caller and the DKFZ bias filter, without scenario 6's extra TLOD/SEVS cut-offs, which mostly remove low-VAF calls. Scenarios 3 and 4 do not require a per-caller PASS. These standard VCFs have no sample columns, so the MAF depth and count columns (`t_depth`, `t_alt_count`, `n_*`) are empty for paired rows.
- Tumour-only: `POST/tumor_only/<patient>_<sample>/*_tumor_only_consensus.PASS.vcf.gz`. Depths come from Mutect2's tumour column. `Matched_Norm_Sample_Barcode` is `NORMAL`, vcf2maf's default for tumour-only VCFs.

VEP runs with the flags vcf2maf uses, plus `--mane`. The picked transcript (`PICK`) prefers MANE Select, then the Ensembl canonical transcript. vcf2maf reports the canonical transcript of the most severely affected gene, which in Ensembl 112 is the MANE Select transcript for most protein-coding genes. MANE IDs are kept as `MANE_SELECT`/`MANE_PLUS_CLINICAL` columns. Multi-allelic records are split, so each ALT allele gets one MAF row. Each task fails if the VEP VCF or the MAF has a different number of records than its input. Sample barcodes follow the caller VCF columns: `<patient>_<sample>_tumor` and `<patient>_<normal sample>_normal`.

Outputs:

- `VEP/<patient>_<sample>/<patient>_<sample>.<paired|tumor_only>.vep.vcf.gz` (+ `.tbi`)
- `MAF/<patient>_<sample>/<patient>_<sample>.<paired|tumor_only>.maf`
- `MAF/cohort.maf`: all samples, with a final `Analysis_Mode` column (`paired` or `tumor_only`)

A cache for another release or assembly: download the indexed cache `homo_sapiens_vep_<release>_<assembly>.tar.gz` from `https://ftp.ensembl.org/pub/release-<release>/variation/indexed_vep_cache/`, extract it under `--vep_cache`, and pin the same ensembl-vep release in `yml/vep.yml`.
<!-- <<< vep -->
<!-- >>> qc -->
#### Sample identity (somalier) and the MultiQC report

Conpair checks identity only for tumour–normal pairs, so tumour-only cohorts get no identity check from it. `--tool somalier` runs `somalier extract` on every BAM of the run (tumours, normals and normal-only pool rows), from either entry point, then one `somalier relate` across them all. A checker then compares the result with the sample sheet.

| Option | Default | What it does |
|--------|---------|--------------|
| `--somalier_sites` | `release` | `release`: the public somalier sites (GRCh38 `sites.hg38.vcf.gz`, GRCh37 `sites.GRCh37.vcf.gz`), read from `<database>/Databases/<genome>/somalier/` or `--somalier_release_sites <file>`. `targets`: built once per run from common gnomAD SNPs (`--somalier_min_af` 0.2 to `--somalier_max_af` 0.8, from the Mutect2 germline resource) inside the padded `--targets`. Or a path to your own sites VCF |
| `--somalier_swap_max_relatedness` | `0.6` | Two samples of the same patient below this are flagged `swap` |
| `--somalier_dup_min_relatedness` | `0.7` | Samples of different patients at or above this are flagged `duplicate` (same person twice, or a mislabel) |
| `--somalier_min_sites` | `50` | Pairs compared at fewer sites are not judged. Same-patient pairs below this are flagged `insufficient_sites` |
| `--somalier_min_x_sites` | `10` | Informative chrX sites needed before sex is inferred |
| `--multiqc` | `true` | One MultiQC report at the end of the run. `--multiqc false` skips it |

The release sites are mostly exonic and miss most of a small panel. On a few-hundred-kb panel, use `--somalier_sites targets` and check the site count in the `SOMALIER_SITES` log line.

Outputs:

- `somalier/`: `somalier.pairs.tsv`, `somalier.samples.tsv`, `somalier.html`, `identity_sex.tsv` (inferred sex per sample), `extract/*.somalier`. Digests are named `<patient>_<sample>_<status>`, not by the BAM's `@RG SM`.
- `REPORT/identity_flags.tsv` has one row per problem: `swap`, `duplicate`, `insufficient_sites` or `sex_mismatch`. It always has a header, and the log warns when it has rows.
- Sex is checked only when the sheet has an optional `sex` column (`male`/`female`). Normals are checked against the X heterozygosity and X/Y depth. Tumours are flagged only for evidence a deletion cannot create (Y coverage in a declared female, or X heterozygosity in a declared male), because tumours often lose chrY or one chrX.

Thresholds are looser than for germline–germline checks, because copy-number changes and LOH in the tumour lower its relatedness to its own normal. Unrelated people score around 0, first-degree relatives around 0.5, and the same person around 1.0. At moderate purity, LOH sites have intermediate allele balance. Somalier does not genotype them, so they reduce the number of sites compared rather than the relatedness. In the synthetic test, a tumour at 70% purity with LOH across 4 chromosomes still scored 1.00, at 345 shared sites against 397 to 399 for the other pairs. Very pure tumours with extensive LOH turn heterozygous sites into homozygous calls, and that does lower relatedness. For such tumours, lower `--somalier_swap_max_relatedness`. A flag is a reason to check, not proof of a swap.

`MultiQC/multiqc_report.html` combines, where present, FastQC (`mapping` entry), Picard MarkDuplicates metrics (`mapping`/`markdup`, now always kept in `MKDUP/`), GATK BaseRecalibrator tables (`mapping`/`markdup`/`recalibration`), mosdepth, somalier, and `bcftools stats` on the consensus VCFs (`MultiQC/bcftools_stats/`). It runs once, after all of these are finished. The report shows sample IDs, so keep it with the results.
<!-- <<< qc -->
<!-- >>> purecn -->
#### PureCN: purity, ploidy and allele-specific copy number (`--tool purecn`)

[PureCN](https://bioconductor.org/packages/PureCN/) 2.16 (env `yml/purecn.yml`) estimates tumour purity and ploidy and calls integer, allele-specific copy number from targeted data, with or without a matched normal (the open analogue of Tempus's CONA). Needs `--type panel` or `--type exome`, `--genome GRCh38` (`hg38`) or `GRCh37` (`hg19`, chr-named reference). Steps, all in `Modules/PURECN/`:

1. `PURECN_INTERVALS` (once): `IntervalFile.R` on the capture BED (`--targets`, or the genome's exome BED with `--type exome` and no `--targets`): GC content and gene symbols.
2. `PURECN_COVERAGE` (every BAM, normal-only rows included): `Coverage.R`, GC-normalised unless `--purecn_gc_norm false`.
3. `PURECN_NORMALDB` (once): `NormalDB.R` from all normals of the run (at least 2), published to `PureCN/normaldb/` with the intervals it was built on.
4. `PURECN_SITES` (once) + `PURECN_SNPS` (per tumour): allele counts (`bcftools mpileup`) at the gnomAD SNPs with AF ≥ `--purecn_min_pop_af` inside the padded targets (from the Mutect2 germline resource), plus that tumour's Mutect2 PASS SNVs. The Mutect2 calls are only read, never changed. Sites get `INFO/DB` and `INFO/POP_AF`; in paired runs the Mutect2 PASS sites also get `SOMATIC`, and the matched normal is the second sample.
5. `PURECN` (per tumour): `PureCN.R --fun-segmentation PSCBS` (joint log-ratio and B-allele frequency).

| Option | Default | What it does |
|--------|---------|--------------|
| `--purecn_normaldb <dir>` | | Use a prebuilt normal DB (`normalDB*.rds` + `intervals.txt`, i.e. `PureCN/normaldb/` from an earlier run of the **same assay**) instead of building one. Needed for a tumour-only cohort with no normal rows |
| `--purecn_matched_normal_coverage` | `false` | Paired tumours: normalise coverage against the matched normal instead of the normal DB. PureCN recommends the DB (it warns that `--normal` overrides its denoising). Without a DB (sheet with one normal), the matched normal is used |
| `--purecn_min_purity`, `--purecn_max_purity` | `0.15`, `0.95` | Purity search range |
| `--purecn_purity_window` | `0.15` | With a sheet column `purity` (pathology estimate, `0.55` or `55`): PureCN still searches the full range, then the most likely of its solutions within purity ± window is reported. If that is not PureCN's first choice, the outputs are regenerated for it through PureCN's curation path (`<id>.csv` then says `Curated` and why). `<id>_purity_check.tsv` records the decision and every candidate solution; if no solution falls within the window, PureCN's first choice is kept and the file says so. The sheet purity does not narrow the search: on the synthetic test, a 0.40–0.70 search around 0.55 (truth 0.60) produced a doubled-ploidy solution (0.70 / 5.2) |
| `--purecn_mappability <bigWig>` | | Mappability track for `IntervalFile.R` (recommended for real assays) |
| `--purecn_off_target` | `false` | Also use off-target bins |
| `--purecn_gc_norm` | `true` | GC loess normalisation of each BAM's coverage (`Coverage.R`). Keep it on for real data, which has GC bias. `false` skips it; build and use a normal DB with the same setting |
| `--purecn_min_pop_af` | `0.001` | gnomAD AF threshold for the SNP sites and the `DB` flag |
| `--purecn_model`, `--purecn_sex`, `--purecn_args` | `betabin`, `?`, | Passed to `PureCN.R` (`--model`, `--sex`, extra arguments verbatim) |

The run stops at startup if `purecn` is selected with neither normal rows nor `--purecn_normaldb`, if the sheet has one normal and tumour-only patients, if `--type genome`, or if a `purity` value is not a fraction or percentage.

Outputs in `PureCN/<patient>_<sample>/`: `<id>.csv` (purity, ploidy, flags), `<id>_dnacopy.seg` (segments with integer `C`), `<id>_genes.csv` (gene-level copy number, amplification/deletion calls), `<id>_loh.csv` (allele-specific segments, LOH), `<id>_variants.csv` and `<id>.vcf.gz` (somatic/germline posteriors, cellular fraction, multiplicity), `<id>.pdf`, `<id>_chromosomes.pdf`, `<id>_local_optima.pdf`, `<id>.rds` (for manual curation with `PureCN.R --rds`), and `<id>_purecn_input.vcf.gz` (the allele counts). `PureCN/coverage/` keeps every BAM's coverage, so a bigger normal DB can be built later. Paired and tumour-only results of the same tumour are not interchangeable: tumour-only somatic/germline calls rely on population frequencies and the model.

Tested on synthetic data only (`synthetic_test/purecn/`): six 2 Mb regions with planted gains, an amplification, losses and copy-neutral LOH at purity 0.6. The test runs with `--purecn_gc_norm false` only because the simulated coverage has no GC signal: there the loess fit has nothing to model and raised the per-target log-ratio SD from 0.13 to 1.4, so PureCN found no solution.
<!-- <<< purecn -->
<!-- >>> pcgr -->
#### Clinical interpretation with PCGR (`--tool pcgr`)

[PCGR](https://sigven.github.io/pcgr/) 2.3.2 annotates each tumour's consensus SNVs/indels (VEP 115 + CIViC, CGI, ClinVar, COSMIC/TCGA, …), assigns AMP/ASCO/CAP tiers and writes an HTML report, an Excel workbook and TSVs. It runs on pairs and tumour-only patients from every entry point.

**Install once** (conda, no root; the PCGR lock files are in `yml/pcgr/`):

```bash
D=/path/to/somaticbunny_resources/pcgr
bash Modules/PCGR/install_pcgr.sh $D envs             # $D/conda/pcgr + $D/conda/pcgrr (keep them side by side)
bash Modules/PCGR/install_pcgr.sh $D refdata grch38   # PCGR bundle 20260620, 7.7 GB download (7.3 GB unpacked)
bash Modules/PCGR/install_pcgr.sh $D vep GRCh38       # Ensembl VEP 115 cache, 25.3 GB download
# GRCh37: same with `refdata grch37` and `vep GRCh37`
```

**Run:** add `--tool pcgr --pcgr_env $D/conda/pcgr --pcgr_refdata $D/refdata/20260620 --pcgr_vep_dir $D/vep_cache`.

| Option | Default | What it does |
|--------|---------|--------------|
| `--pcgr_scenario` | `5` | Which paired consensus feeds PCGR: `POST/standard_vcf/<id>/<id>_scenario_N_standard_PASSed.vcf`, N = 3 or 4 (2-of-4 SNV / 2-of-3 indel over all caller records), 5 (only records that PASS in each caller, plus the DKFZ bias filter) or 6 (5 with Mutect2 TLOD ≥ 10 and Strelka2 SomaticEVS ≥ 13) |
| `--pcgr_assay` | from `--type` | `TARGETED` for panel, `WES` for exome, `WGS` for genome |
| `--pcgr_effective_target_size_mb` | merged `--targets` size | TMB denominator. Without `--targets`, 34 Mb (PCGR's coding exome); capped at 34 |
| sheet column `tumor_site` | `--pcgr_tumor_site 0` | Optional PCGR primary-site code per tumour row (0 any, 6 breast, 9 colon/rectum, 15 lung, 23 prostate, …; full list in `pcgr --help`). Focuses biomarker matching and signature priors |
| `--pcgr_msi`, `--pcgr_signatures` | on | MSI and signature fitting; see caveats |
| `--pcgr_tumor_dp_min/_af_min`, `--pcgr_tmb_dp_min/_af_min` | unset | PCGR report/TMB depth and VAF filters |
| `--pcgr_gnomad_popmax_af`, `--pcgr_tonly_germline_filters` | `0.001`, on | Tumour-only germline filters (gnomAD popmax AF; likely hom/het germline) |
| `--pcgr_pon_vcf` | unset | Tumour-only: extra panel-of-normals VCF (records flagged `INFO/PANEL_OF_NORMALS`), excluded with `--exclude_pon` |
| `--pcgr_cna_from_purecn` | off | Hook: use PureCN allele-specific segments as PCGR `--input_cna` (converter `Modules/PCGR/purecn_to_pcgr_cna.py`; wire `PURECN` where main.nf marks the hook) |
| `--pcgr_extra_args` | | Passed to `pcgr` verbatim (e.g. `--vep_gencode_basic`) |

**Input VCF.** PCGR reads tumour depth and allele fraction only from INFO. The paired standard VCF has no FORMAT columns, so `Modules/PCGR/pcgr_prepare_vcf.py` keeps CHROM/POS/REF/ALT of the consensus and writes `TDP`/`TAF` (tumour) and `NDP`/`NAF` (matched normal) from the callers' FORMAT read counts: Mutect2 AD/DP, else SAGE AD/DP, else Strelka2 tier-1 counts (`TVAF_SRC` says which). Tumour-only input is `POST/tumor_only/<id>/<id>_tumor_only_consensus.PASS.vcf.gz`, whose Mutect2 FORMAT AD/DP give `TDP`/`TAF`. Only chromosomes 1–22, X, Y are kept. PCGR runs with `--tumor_dp_tag TDP --tumor_af_tag TAF` (and `--control_dp_tag NDP --control_af_tag NAF` for pairs) and `--tumor_only` for tumour-only patients.

**Caveats for panels.** PCGR predicts MSI only for WES/WGS tumour–normal data, so it is off for `--type panel` and for tumour-only patients (use the MSIsensor-pro step). PCGR TMB on a panel smaller than 1 Mb is not robust, and the clean target BED also counts intronic and UTR baits, so set `--pcgr_effective_target_size_mb` to the panel's coding footprint when you know it. Signature fitting needs at least `--pcgr_min_mutations_signatures` SNVs (200; minimum 100), which most panels do not reach, so it is skipped. Tumour-only TMB and signatures still include rare germline variants. OncoKB annotation is off: it needs a licensed API token (`--pcgr_extra_args "--oncokb_api_token …"`).

**Output:** `PCGR/<patient>_<sample>/`: `<id>.pcgr.grch38.html` (report), `.xlsx`, `.snv_indel_ann.tsv.gz`, `.tmb.tsv`, `.pass.vcf.gz` and the rest of PCGR's files, plus the prepared input (`<id>.pcgr_input.vcf.gz`) and a one-line count summary. **These reports contain patient-level clinical results. Keep them on the cluster in a protected `--outdir`. Do not email, upload or commit them.**
<!-- <<< pcgr -->

#### Consensus output: final VCFs with commands in the header, and "no consensus mutations"

The paired consensus is upstream's `POST_SCENARIOS` (`consensus_scenarios.sh`, all six scenarios in one task). It writes valid sites-only VCFs and continues when a tumour has no candidate SNVs. On top of it, the clinical branch adds `POST_FINAL_VCF` for scenarios 3–6:

- `POST/<patient>_<sample>/scenario_N/<id>_final.vcf.gz` (+ `.tbi`): upstream's `<id>_final.vcf`, sorted and indexed. Its header records the exact commands that produced it:
  - `##SomaticBunny_consensusCommand`: the `consensus_scenarios.sh` call, with its arguments.
  - `##SomaticBunny_combineCommand`
  - `##SomaticBunny_bcftoolsAnnotateCommand`, `##SomaticBunny_bcftoolsSortCommand`, `##SomaticBunny_bcftoolsIndexCommand` and `##SomaticBunny_bcftoolsValidateCommand`.
  - bcftools' own `##bcftools_*` lines.

  The run fails if this file does not parse. Upstream's plain `<id>_final.vcf` is unchanged.
- `REPORT/consensus_status.tsv`: one row per tumour and scenario (`n_snv`, `n_indel`, `status`). `NO_CONSENSUS_MUTATIONS` marks a tumour with neither SNVs nor indels, and the log warns. The run continues, and AlleleCounter (VAF) is skipped for those tumours because it crashes on an empty VCF.

VEP (`--vep_scenario`), PCGR (`--pcgr_scenario`), TMB and the MultiQC bcftools stats read these scenario VCFs.

### Advanced Options

#### Publishing Intermediate Files

By default, intermediate files are not saved to reduce disk usage. Use the `--publish` parameter to save specific file types:

```bash
nextflow run main.nf -profile tscc_slurm,conda --type genome --genome GRCh38 --first_step mapping --publish raw_bam,markdup_bam
```

**Available options:**
- `raw_bam`: Publish raw BAM files from BWA-MEM
- `markdup_bam`: Publish mark duplicates BAM files
- `recal_bam`: Publish recalibrated BAM files

#### Resuming Failed Runs

If your pipeline terminates with an error or the interactive node is killed, resume with:

```bash
nextflow run main.nf -profile tscc_slurm,conda --type genome --genome GRCh38 --first_step mapping --tool ascat,manta -resume
```

#### Email Notifications

Receive completion notification:

```bash
nextflow run main.nf -profile tscc_slurm,conda --type genome --genome GRCh38 --first_step mapping --tool ascat,manta -N your_email@example.com
```

## Output

All results and reports are stored in the **RESULTS** folder with the following structure:

```
RESULTS/
├── ASCAT (optional)
├── BAM (optional)
├── CNVkit (optional)
├── Conpair
├── Delly (optional)
├── FASTQC
├── MANTA (optional)
├── MKDUP
├── mosdepth
├── MuSE
├── Mutect2
├── POST
│   ├── <patient>_<sample>/scenario1 ... scenario_6   consensus SNV and INDEL calls per scenario
│   └── VAF                                           allele counts and VAFs for the scenario 4-6 consensus calls
├── RECALIBRATE (optional)
├── REPORT
├── SAGE
└── STRELKA
```

### Variant allele fractions (`RESULTS/POST/VAF`)

For each tumor sample and each of scenarios 4, 5 and 6, the pipeline writes
`<patient>_<sample>_<scenario>_VAF.tsv`: one row per consensus variant, 21 fixed
tab-separated columns, header without `#`, so it loads directly into pandas or R.

| Column | Content |
|--------|---------|
| `sample` | `<patient>_<sample>_<scenario>` |
| `CHROM`, `POS`, `ID`, `REF`, `ALT` | as in the consensus VCF |
| `TYPE` | `SNV`, `INDEL`, `MNP`, `MULTIALLELIC` or `INVALID` |
| `FILTER` | consensus FILTER (`PASS`, `OnePassCaller`, `NoPassCaller`) |
| `CALLERS`, `PASS_CALLERS`, `CALLER_FILTERS` | callers that reported the variant, those that passed, and each caller's filter value |
| `TUMOR_REF`, `TUMOR_ALT`, `TUMOR_DEPTH`, `TUMOR_VAF` | tumor counts and `ALT / DEPTH` |
| `NORMAL_REF`, `NORMAL_ALT`, `NORMAL_DEPTH`, `NORMAL_VAF` | the same for the matched normal |
| `STATUS` | `OK`, `NO_COVERAGE` (depth 0 in tumor or normal), or why counts are `NA` (`MNP`, `MULTIALLELIC`, `INVALID`) |
| `NORMALIZED` | `.`, or `POS:REF>ALT` of the left-aligned indel that was counted when it differs from the VCF record |

SNV counts come from alleleCounter (base quality >= 20, mapping quality >= 35); depth is
alleleCounter's `Good_depth`. INDEL counts are taken from the BAMs with the same thresholds:
a read supports `ALT` only if it carries exactly the called event after left-normalisation,
supports `REF` only if it is aligned without any indel across the event, and is otherwise
excluded; depth is `REF + ALT`. The raw alleleCounter outputs
(`*_tumor_alleleCounter.txt`, `*_normal_alleleCounter.txt`) are kept next to the table.

## Tool Versions

| Category | Tool | Version |
|----------|------|---------|
| QC | FastQC | v0.12.1 |
| Data Processing | Picard | v3.4.0 |
| | samtools | v1.24 |
| | bwa-mem2 | v2.3 |
| | Conpair | v0.2 |
| | GATK4 | v4.6.2.0 |
| | mosdepth | v0.3.14 |
| SNV/INDEL Callers | Strelka2 | v2.9.10 |
| | Mutect2 | v4.6.2.0 |
| | SAGE | v4.2 |
| | MuSE2 | v2.1.2 |
| Structural Variants | Delly | v2.6.0 |
| | Manta | v1.6.0 |
| Copy Number Variants | CNVkit | v0.9.13 |
| | ASCAT | v3.2.0 |

## Benchmark

Code to reproduce the benchmark against nf-core/sarek is in [`benchmark/`](benchmark).
