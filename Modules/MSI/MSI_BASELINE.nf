nextflow.enable.dsl=2

// Tumour-only MSI baseline (msisensor-pro baseline) from the per-normal `pro` outputs of this run.
// Published to MSI/baseline/ so later runs of the same assay can pass it as --msi_baseline.
process MSI_BASELINE {
    conda "${params.msi_env}"
    scratch true
    label 'process_low'
    publishDir("${params.msi_dir}/baseline", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(normal_all, stageAs: 'normals/*')
    path(sites, stageAs: 'input/*')

    output:
    path("msi_baseline.tsv"), emit: baseline
    path("msi_baseline.normals.txt"), emit: normals

    script:
    def files = normal_all instanceof List ? normal_all : [normal_all]
    def n = files.size()
    def support = Math.min(params.msi_baseline_min_normals as int, n)
    """
    set -euo pipefail
    : > configure.txt
    for f in ${files.join(' ')}; do
        printf '%s\\t%s\\n' "\$(basename "\$f" _all)" "\$(readlink -f "\$f")" >> configure.txt
    done
    cut -f1 configure.txt > msi_baseline.normals.txt
    echo "MSI_BASELINE: ${n} normal(s), sites need support from >= ${support}" >&2
    msisensor-pro baseline -d ${sites} -i configure.txt -o msi_baseline.tsv -s ${support}
    if [ ! -s msi_baseline.tsv ]; then
        echo "ERROR: msisensor-pro baseline wrote no baseline" >&2
        exit 1
    fi
    """
}
