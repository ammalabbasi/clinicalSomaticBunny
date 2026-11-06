# Running the pipeline

[&larr; back to the main README](../README.md)

## Basic Usage

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

### Pipeline Steps and Tools

The SomaticBunny pipeline workflow is divided into multiple steps that can be run individually or in sequence using the `--first_step` parameter.

#### Available first_step

| First_step | Description | Required Input | 
|------|-------------| ---------------|
| `mapping` | This will start from performing alignment of FASTQ files to the reference genome | fastq_1 and fastq_2 |
| `markdup` | This will start from mark duplicate reads using the aligned BAM files | raw bam and bai |
| `recalibration` | This will start from performing base quality score recalibration | mark duplicated bam and bai |
| `variant_calling` | This will start from executing variant calling tools | recalibrated bam and bai |

#### Available reference genomes
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

#### Available Tools

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

## Custom Reference Genome

When starting the pipeline from BAM files (`--first_step markdup`, `recalibration`, or `variant_calling`), you **must** provide your own reference genome files. This is because BAM files are already aligned to a specific reference, and the pipeline needs the matching reference for downstream analysis.

### Required parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--ref` | Reference FASTA file | `/path/to/hg38_genome.fa` |
| `--ref_fai` | FASTA index file | `/path/to/hg38_genome.fa.fai` |
| `--ref_dict` | Sequence dictionary | `/path/to/hg38_genome.dict` |
| `--bed` | Callable regions BED file (bgzipped) | `/path/to/callable_regions.bed.gz` |
| `--bed_tbi` | Tabix index for the BED file | `/path/to/callable_regions.bed.gz.tbi` |

All five files must be provided together. The pipeline will exit with an error if any are missing.

### Example

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

## Advanced Options

### Publishing Intermediate Files

By default, intermediate files are not saved to reduce disk usage. Use the `--publish` parameter to save specific file types:

```bash
nextflow run main.nf -profile tscc_slurm,conda --type genome --genome GRCh38 --first_step mapping --publish raw_bam,markdup_bam
```

**Available options:**
- `raw_bam`: Publish raw BAM files from BWA-MEM
- `markdup_bam`: Publish mark duplicates BAM files
- `recal_bam`: Publish recalibrated BAM files

### Resuming Failed Runs

If your pipeline terminates with an error or the interactive node is killed, resume with:

```bash
nextflow run main.nf -profile tscc_slurm,conda --type genome --genome GRCh38 --first_step mapping --tool ascat,manta -resume
```

### Email Notifications

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
