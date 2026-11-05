nextflow.enable.dsl=2

// One MultiQC report over every QC file the run produced (FastQC, Picard MarkDuplicates metrics, GATK
// BaseRecalibrator tables, mosdepth, somalier, bcftools stats). main.nf mixes and collects the published files, so
// this runs once, after everything else. Each input is staged in its own numbered directory so files with the same
// name from different samples cannot collide.
process MULTIQC {
    conda "${params.multiqc_env}"
    scratch true
    label 'process_low'
    publishDir("${params.multiqc_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(qc_files, stageAs: 'mqc_in/f*/*')

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_report_data"), emit: data

    script:
    def cfg = params.multiqc_config ? "--config ${params.multiqc_config}" : ""
    """
    export MPLCONFIGDIR=\$PWD/.mpl HOME=\$PWD
    multiqc -f ${cfg} --cl-config "no_version_check: true" \\
        --title "SomaticBunny QC" -n multiqc_report.html -o . mqc_in/
    """
}
