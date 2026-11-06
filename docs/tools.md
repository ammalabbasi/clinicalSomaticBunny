# Tool packages and versions

[&larr; back to the main README](../README.md)

Every step runs in its own pinned conda environment (`yml/*.yml`); the exact builds are pinned so a fresh
install solves to the same binaries. PCGR uses a fully locked environment (`yml/pcgr/*.lock`). SAGE and GATK ship
with the reference database.

| Stage | Tool | Version | Used for | Env file |
|---|---|---|---|---|
| Alignment (FASTQ entry) | BWA-MEM2 | 2.3 | read mapping | `yml/bwamem2.yml` |
| | samtools | 1.24 | BAM handling throughout | several envs |
| UMI consensus | fgbio | 3.1.2 | `--umi`: GroupReadsByUmi, consensus, filter | `yml/umi.yml` |
| Read QC | FastQC | 0.12.1 | FASTQ QC | `yml/fastqc_env.yml` |
| Duplicates / BQSR | Picard | 3.4.0 | MarkDuplicates, header rename | `yml/picard.yml` |
| | GATK4 | 4.6.2.0 | BQSR, Mutect2, PON, contamination | database bundle |
| Coverage | mosdepth | 0.3.14 | depth over targets; TMB callable bases | `yml/mosdepth_env.yml` |
| Pair QC | Conpair | 0.2 | tumour-normal concordance/contamination | `yml/conpair.yml` |
| SNV/indel callers | Mutect2 (GATK4) | 4.6.2.0 | paired + tumour-only calling | database bundle |
| | Strelka2 | 2.9.10 | paired calling on panel call regions | `yml/strelka_env.yml` |
| | MuSE | 2.1.2 | paired calling | `yml/muse.yml` |
| | SAGE | 4.2 | paired + tumour-only calling (Hartwig) | jar in database |
| Consensus | bcftools | 1.9 / 1.22 | scenario voting; tumour-only consensus; final VCFs | `yml/POST.yml`, `yml/bcftools.yml` |
| | DKFZ bias filter | 1.2.3a | scenario 6 strand-bias filter | `yml/POST.yml` |
| VAF | alleleCounter | 4.3.0 | allele counts at consensus sites | `yml/alleleCounter.yml` |
| Copy number | CNVkit | 0.9.13 | log-ratio segments, pooled/flat reference | `yml/cnvkit.yml` |
| | PureCN | 2.16.1 | purity, ploidy, integer CN, LOH; per-assay normal DB | `yml/purecn.yml` |
| Biomarkers | MSIsensor-pro | 1.3.0 | MSI paired + tumour-only vs baseline | `yml/msisensor_pro.yml` |
| | TMB module | - | mutations per callable Mb | `yml/tmb.yml` |
| | gene_status module | python 3.12 / pysam 0.22.1 / mosdepth 0.3.14 | per-gene status: coverage + consensus variants + CN | `yml/gene_status.yml` |
| Annotation | Ensembl VEP | 112.0 | offline annotation (pinned cache) | `yml/vep.yml` |
| | vcf2maf | 1.6.22 | per-sample MAF + cohort MAF | `yml/vep.yml` |
| Identity QC | somalier | 0.3.7 | duplicates, swaps, sex (works tumour-only) | `yml/somalier.yml` |
| Run report | MultiQC | 1.34 | one report per run (+ bcftools stats 1.25) | `yml/multiqc.yml` |
| Clinical report | PCGR | 2.3.2 | AMP/ASCO/CAP-tiered report per tumour | `yml/pcgr/` (locked) |
| WGS/WES only | Manta | 1.6.0 | SVs (off for panels) | `yml/manta.yml` |
| | Delly | 2.6.0 | SVs (off for panels) | `yml/delly.yml` |
| | ASCAT | 3.2.0 | allele-specific CN (off for panels) | `yml/ascat.yml` |

A note on pinning: loosely specified conda environments resolve differently over time - on this pipeline an
unpinned respec once produced a Strelka build that segfaulted and a consensus env whose solver ran for hours.
Keep the pins.
