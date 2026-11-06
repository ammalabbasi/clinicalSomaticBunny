nextflow.enable.dsl=2

// Per-gene somatic status for one tumour (--tool gene_status), paired or tumour-only.
// Covered genes = the resolved genes BED (main.nf: --genes_bed, else the 4th column of --targets) intersected
// with the padded clean targets; a gene with no targeted bp is omitted entirely (NOT_COVERED). Depth comes from
// `mosdepth --quantize` on the tumour BAM (reads not duplicate/secondary/supplementary/QC-failed, MAPQ >=
// --gene_status_min_mapq); variants from the same consensus VCFs the other tools read (annotated by VEP when
// --tool vep is on); copy number from CNVkit *_calls.cns and, when PureCN ran, its *_loh.csv (integer CN + LOH).
//
// <patient>_<sample>.gene_status.tsv — one row per covered gene; columns:
//   sample            <patient>_<sample>
//   mode              paired | tumor_only
//   source            consensus the variants came from (scenario_N, tumor_only_consensus, *_vep)
//   gene              name from the genes BED
//   chrom,start,end   the gene's full BED interval (0-based half-open)
//   targeted_bp       bp of gene ∩ padded targets
//   callable_bp       targeted bp with tumour depth >= --gene_status_min_depth
//   callable_frac     callable_bp / targeted_bp
//   n_snv,n_indel     consensus records whose POS falls in the gene's targeted intervals (one per ALT)
//   worst_consequence VEP worst consequence over the gene's variants ('.' without --tool vep or variants)
//   vep_symbol        VEP SYMBOLs seen on the gene's variants ('.' without annotation)
//   cnvkit_cn         CNVkit integer cn of the *_calls.cns segment with the largest gene overlap ('.' if none)
//   purecn_cn         PureCN integer C over the gene (largest overlap in *_loh.csv; '.' if PureCN absent)
//   purecn_minor_cn   PureCN minor allele copy number M ('.' if absent)
//   cn_source         purecn | cnvkit | none — which CN feeds the status
//   status            NOT_CALLABLE, or '+'-joined labels of MUTATED, AMP(cn>=5), GAIN(cn>2), DEL(cn<2),
//                     DEEP_DEL(cn=0), LOH (PureCN M=0, C>0); WT when none apply
process GENE_STATUS {
    conda "${params.gene_status_env}"
    scratch true
    label 'process_low'
    publishDir("${params.gene_status_dir}", mode: 'copy', pattern: "${patient_sample}.gene_status.tsv")
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(patient_sample), val(mode), val(source),
          path(vcf, stageAs: 'input/vcf/*'),
          path(tumor_bam, stageAs: 'input/tumor/*'), path(tumor_bai, stageAs: 'input/tumor/*'),
          path(cnvkit_cns, stageAs: 'input/cnvkit/*'),
          path(purecn_loh, stageAs: 'input/purecn/*')
    path(genes_bed, stageAs: 'input/genes/*')
    path(targets_bed, stageAs: 'input/targets/*')

    output:
    tuple val(patient_sample), path("${patient_sample}.gene_status.tsv"), emit: gene_status_tsv

    script:
    def cnvkit_arg = cnvkit_cns ? "--cnvkit_cns ${cnvkit_cns}" : ""
    def purecn_arg = purecn_loh ? "--purecn_loh ${purecn_loh}" : ""
    """
    set -euo pipefail
    # mosdepth finds the index next to the BAM; vendor indexes may be named X.bai instead of X.bam.bai
    [[ -e "${tumor_bam}.bai" ]] || ln -s "\$(basename ${tumor_bai})" "${tumor_bam}.bai"

    MOSDEPTH_Q0=LOW MOSDEPTH_Q1=CALLABLE \\
        mosdepth -t ${task.cpus} -n -F 3844 --mapq ${params.gene_status_min_mapq} \\
        --quantize 0:${params.gene_status_min_depth}: gsq_tumor ${tumor_bam}

    python3 -E ${projectDir}/Modules/GENE_STATUS/gene_status.py \\
        --sample ${patient_sample} --mode ${mode} --source ${source} \\
        --vcf ${vcf} --genes_bed ${genes_bed} --targets_bed ${targets_bed} \\
        --quantized gsq_tumor.quantized.bed.gz \\
        --min_depth ${params.gene_status_min_depth} --min_callable ${params.gene_status_min_callable} \\
        ${cnvkit_arg} ${purecn_arg} \\
        --out ${patient_sample}.gene_status.tsv
    """
}

// Cohort matrix: genes (rows, BED order) x tumours (columns), the status strings of the per-sample TSVs.
// A gene not covered for a sample (e.g. mixed assays) holds '.'.
process GENE_STATUS_MATRIX {
    label 'process_low'
    publishDir("${params.gene_status_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(status_files, stageAs: 'status/*')

    output:
    path("gene_status_matrix.tsv"), emit: matrix

    script:
    """
    set -euo pipefail
    python3 -E ${projectDir}/Modules/GENE_STATUS/gene_status_matrix.py --out gene_status_matrix.tsv status/*
    """
}
