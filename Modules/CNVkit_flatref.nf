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
    // Targeted: on-target bins only (--method amplicon). With hybrid mode a flat reference keeps every antitarget
    // bin at full weight with no per-bin expectation; empty or noisy off-target bins are then centred to 0 and
    // dilute whole-chromosome segments (synthetic test: CN3 gain read as log2 +0.13). CNVkit.nf must use the
    // same method with this reference.
    def regions = params.targeted ? "--method amplicon --targets \"${params.mosdepth_bed}\"" : "--method wgs"
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
