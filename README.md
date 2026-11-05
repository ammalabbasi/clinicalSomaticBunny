# SomaticBunny

<img src="workflow_logo/SomaticBunny_Logo.png" alt="SomaticBunny Logo" width="300" align="left"/>

**SomaticBunny** is an ensemble pipeline for somatic variant calling that integrates multiple tools to provide comprehensive variant detection in tumor-normal paired samples.

<br clear="left"/>


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
├── RECALIBRATE (optional)
├── REPORT
├── SAGE
└── STRELKA
```

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
