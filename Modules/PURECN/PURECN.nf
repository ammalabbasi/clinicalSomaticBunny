nextflow.enable.dsl=2

// PureCN.R per tumour: purity/ploidy, allele-specific segmentation (PSCBS: joint log-ratio + BAF), integer copy number,
// gene-level calls, LOH and somatic/germline posteriors for every variant in the input VCF.
//   coverage normalisation: the normal DB (tangent normal) when there is one. A paired tumour's matched-normal coverage
//   is used instead only with --purecn_matched_normal_coverage, or when the run has no normal DB (one normal only):
//   PureCN.R warns that --normal overrides the DB's denoising, so the DB is the default.
//   purity search: [--purecn_min_purity, --purecn_max_purity], or the sheet's `purity` +- --purecn_purity_window.
process PURECN {
    tag "${map.patient}_${map.tumor_meta.sample}"
    scratch true
    label 'process_medium'
    conda "${params.purecn_env}"
    publishDir("${params.purecn_dir}/${map.patient}_${map.tumor_meta.sample}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(map), path(vcf, stageAs: 'input/vcf/*'), path(vcf_tbi, stageAs: 'input/vcf/*'),
          path(tumor_cov, stageAs: 'input/t/*'), path(normal_cov, stageAs: 'input/n/*')
    path(normaldb, stageAs: 'input/normaldb')
    path(intervals, stageAs: 'input/intervals/*')

    output:
    tuple val(map), path("${map.patient}_${map.tumor_meta.sample}.csv"), emit: purity_ploidy
    tuple val(map), path("${map.patient}_${map.tumor_meta.sample}_dnacopy.seg"), emit: segments
    path("${map.patient}_${map.tumor_meta.sample}{.rds,.log,*.pdf,_genes.csv,_loh.csv,_variants.csv,_amplification_pvalues.csv,.vcf.gz,.vcf.gz.tbi}"), emit: all

    script:
    def id = "${map.patient}_${map.tumor_meta.sample}"
    def use_db = normaldb ? true : false
    def use_normal = normal_cov && (params.purecn_matched_normal_coverage || !use_db)
    def min_p = params.purecn_min_purity as double
    def max_p = params.purecn_max_purity as double
    def path_purity = params.purecn_purity_map ? params.purecn_purity_map[id.toString()] : null
    if (path_purity != null) {
        def w = params.purecn_purity_window as double
        min_p = Math.max(0.1d, (path_purity as double) - w)
        max_p = Math.min(0.99d, (path_purity as double) + w)
    }
    def purity_args = String.format(Locale.ROOT, "--min-purity %.2f --max-purity %.2f", min_p, max_p)
    """
    set -euo pipefail
    unset R_LIBS R_LIBS_USER R_LIBS_SITE
    X=\$(Rscript --vanilla -e 'cat(system.file("extdata", package = "PureCN"))')
    DB_ARGS=""
    if [ "${use_db}" = "true" ]; then
        db=\$(ls input/normaldb/normalDB*.rds | head -1)
        DB_ARGS="--normaldb \$db"
    fi
    echo "PureCN ${id}: normal DB ${use_db}, matched-normal coverage ${use_normal}, ${purity_args}${path_purity != null ? " (sheet purity ${path_purity})" : ""}"

    Rscript --vanilla \$X/PureCN.R \\
        --out ./${id} \\
        --sampleid ${id} \\
        --tumor ${tumor_cov} \\
        ${use_normal ? "--normal ${normal_cov}" : ""} \\
        \$DB_ARGS \\
        --vcf ${vcf} \\
        --intervals ${intervals} \\
        --genome ${params.purecn_genome} \\
        --fun-segmentation PSCBS \\
        --model ${params.purecn_model} \\
        --sex '${params.purecn_sex}' \\
        ${purity_args} \\
        --seed 123 \\
        --out-vcf \\
        --cores ${task.cpus} \\
        ${params.purecn_args} \\
        --force

    test -s ${id}.csv
    """
}
