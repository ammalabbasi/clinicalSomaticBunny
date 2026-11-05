nextflow.enable.dsl=2

// bcftools stats on each consensus VCF (paired standard VCFs from MERGE_HEADER, tumour-only *.PASS.vcf.gz), for the
// MultiQC "Bcftools" section. Counts only; the VCFs themselves are not copied.
// A VCF that bcftools cannot parse is skipped with a warning (no stats file) instead of stopping the run: this is a
// side report. Seen on synthetic data: MERGE_HEADER standard VCFs with an undeclared INFO/caller_filter and an INFO
// key "maxGermlineVAF,LowEVS". Signals (exit 128-254) still fail the task.
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
    path("*.bcftools_stats.txt"), optional: true, emit: stats

    script:
    // vcf is staged as input/<name>; keep the output in the task root, named after the file only
    def base = vcf.toString().tokenize('/')[-1].replaceAll(/\.vcf(\.gz)?$/, '')
    """
    rc=0
    bcftools stats ${vcf} > ${base}.bcftools_stats.txt || rc=\$?
    if [ \$rc -ge 128 ] && [ \$rc -lt 255 ]; then exit \$rc; fi
    if [ \$rc -ne 0 ]; then
        echo "WARNING: bcftools stats could not parse ${base} (exit \$rc); skipped in MultiQC" >&2
        rm -f ${base}.bcftools_stats.txt
    fi
    """
}
