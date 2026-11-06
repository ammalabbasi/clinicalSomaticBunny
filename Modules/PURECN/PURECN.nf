nextflow.enable.dsl=2

// PureCN.R per tumour: purity/ploidy, allele-specific segmentation (PSCBS: joint log-ratio + BAF), integer copy number,
// gene-level calls, LOH and somatic/germline posteriors for every variant in the input VCF.
//   coverage normalisation: the normal DB (tangent normal) when there is one. A paired tumour's matched-normal coverage
//   is used instead only with --purecn_matched_normal_coverage, or when the run has no normal DB (one normal only):
//   PureCN.R warns that --normal overrides the DB's denoising, so the DB is the default.
//   purity search: always [--purecn_min_purity, --purecn_max_purity]. A sheet `purity` does NOT bound the search:
//   bounding it (0.40-0.70 around 0.55) sent PureCN to a doubled-ploidy solution on the synthetic test although the
//   true purity (0.60) was inside the bounds. Instead resources/purecn_select.R picks, among PureCN's ranked solutions,
//   the most likely one within sheet purity +- --purecn_purity_window. If that is not rank 1, it writes PureCN's
//   curation file and the outputs are regenerated from the .rds (PureCN.R --rds, PureCN's manual-curation path; the
//   CSV then says Curated=TRUE and why). <id>_purity_check.tsv records the decision and all candidate solutions.
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
    // >>> gene_status
    tuple val(map), path("${map.patient}_${map.tumor_meta.sample}_loh.csv"), optional: true, emit: loh
    // <<< gene_status
    path("${map.patient}_${map.tumor_meta.sample}{.rds,.log,*.pdf,_genes.csv,_loh.csv,_variants.csv,_amplification_pvalues.csv,.vcf.gz,.vcf.gz.tbi}"), emit: all
    path("${map.patient}_${map.tumor_meta.sample}_purity_check.tsv"), optional: true, emit: purity_check

    script:
    def id = "${map.patient}_${map.tumor_meta.sample}"
    def use_db = normaldb ? true : false
    def use_normal = normal_cov && (params.purecn_matched_normal_coverage || !use_db)
    def path_purity = params.purecn_purity_map ? params.purecn_purity_map[id.toString()] : null
    def purity_args = String.format(Locale.ROOT, "--min-purity %.2f --max-purity %.2f",
        params.purecn_min_purity as double, params.purecn_max_purity as double)
    def sheet = path_purity != null ? String.format(Locale.ROOT, "%.4f", path_purity as double) : ""
    def window = String.format(Locale.ROOT, "%.4f", params.purecn_purity_window as double)
    """
    set -euo pipefail
    unset R_LIBS R_LIBS_USER R_LIBS_SITE
    X=\$(Rscript --vanilla -e 'cat(system.file("extdata", package = "PureCN"))')
    DB_ARGS=""
    if [ "${use_db}" = "true" ]; then
        db=\$(ls input/normaldb/normalDB*.rds | head -1)
        DB_ARGS="--normaldb \$db"
    fi
    echo "PureCN ${id}: normal DB ${use_db}, matched-normal coverage ${use_normal}, ${purity_args}${sheet ? ", sheet purity ${sheet} +- ${window} (selection)" : ""}"

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

    if [ -n "${sheet}" ]; then
        rc=0
        Rscript --vanilla ${moduleDir}/resources/purecn_select.R ./${id} ${sheet} ${window} || rc=\$?
        if [ \$rc -eq 3 ]; then
            # regenerate every output for the selected solution (reads ./${id}.csv through readCurationFile)
            Rscript --vanilla \$X/PureCN.R \\
                --rds ./${id}.rds \\
                --out ./${id} \\
                --sampleid ${id} \\
                --intervals ${intervals} \\
                --genome ${params.purecn_genome} \\
                --out-vcf
        elif [ \$rc -ne 0 ]; then
            exit \$rc
        fi
        test -s ${id}_purity_check.tsv
    fi
    """
}
