nextflow.enable.dsl=2

// Flat CNVkit reference: every bin gets the same expected coverage (log2 0), with GC and repeat corrections
// from the FASTA but no per-bait capture correction. For cohorts with no normal BAMs (--cnvkit_reference flat,
// or auto with no normal rows). Noisier than a pooled-normal reference; large events remain visible.
process CNVkit_flatref {
    conda "${params.cnvkit_env}"
    scratch true
    label 'process_medium'
    publishDir("${params.cnvkit_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(trigger)

    output:
    path("reference.flat.cnn"), emit: CNVkit_ref_cnn

    script:
    def regions = params.targeted ? "--targets \"${params.mosdepth_bed}\"" : "--method wgs"
    """
    export PYTHONNOUSERSITE=1
    unset PYTHONPATH

    cnvkit.py batch \\
        -n \\
        ${regions} \\
        --fasta "${params.ref}" \\
        --output-reference reference.flat.cnn \\
        -p ${task.cpus}
    test -s reference.flat.cnn
    """
}
