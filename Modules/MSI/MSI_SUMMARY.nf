nextflow.enable.dsl=2

// Cohort table: one row per tumour (paired and tumour-only), from the per-sample *.msi_status.tsv files.
process MSI_SUMMARY {
    label 'process_low'
    publishDir("${params.msi_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(status_files, stageAs: 'status/*')

    output:
    path("msi_summary.tsv"), emit: summary

    script:
    """
    set -euo pipefail
    first=\$(ls status/* | head -1)
    head -1 "\$first" > msi_summary.tsv
    for f in status/*; do tail -n +2 "\$f"; done | sort -t\$'\\t' -k3,3 -k1,1 -k2,2 >> msi_summary.tsv
    """
}
