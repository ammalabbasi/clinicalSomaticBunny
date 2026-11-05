nextflow.enable.dsl=2

// bcftools stats on each consensus VCF (paired standard VCFs from MERGE_HEADER, tumour-only *.PASS.vcf.gz), for the
// MultiQC "Bcftools" section. Counts only; the VCFs themselves are not copied.
process BCFTOOLS_STATS {
    conda "${params.multiqc_env}"
    scratch true
    label 'process_low'
    publishDir("${params.multiqc_dir}/bcftools_stats", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(vcf, stageAs: 'input/*')

    output:
    path("*.bcftools_stats.txt"), emit: stats

    script:
    // vcf is staged as input/<name>; keep the output in the task root, named after the file only
    def base = vcf.toString().tokenize('/')[-1].replaceAll(/\.vcf(\.gz)?$/, '')
    """
    bcftools stats ${vcf} > ${base}.bcftools_stats.txt
    """
}
