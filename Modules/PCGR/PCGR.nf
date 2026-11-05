nextflow.enable.dsl=2

// PCGR 2.3.2 (Personal Cancer Genome Reporter): VEP + vcfanno annotation, AMP/ASCO/CAP tiers, TMB, MSI (WES/WGS
// tumour-normal only), mutational signatures, HTML/Excel/TSV report.
//
// Input VCF: SomaticBunny consensus, rewritten by pcgr_prepare_vcf.py into a clean single-sample VCF whose INFO
// carries TDP/TAF (tumour) and NDP/NAF (matched normal, pairs only) taken from caller FORMAT fields, because PCGR
// reads depth/allele fraction from INFO only.
//   mode 'paired'     : POST/standard_vcf/<id>/<id>_scenario_<pcgr_scenario>_standard_PASSed.vcf + Mutect2/SAGE/Strelka2 VCFs
//   mode 'tumor_only' : POST/tumor_only/<id>/<id>_tumor_only_consensus.PASS.vcf.gz (Mutect2 records, FORMAT AD/DP)
// Optional CNA: PureCN <sample>_loh.csv (converted to PCGR's nMajor/nMinor TSV); NO_FILE when absent.
//
// The two conda envs (pcgr, pcgrr) are prebuilt from the PCGR lock files (yml/pcgr/, Modules/PCGR/install_pcgr.sh) and
// must be siblings: pcgr finds the R env as $(dirname $CONDA_PREFIX)/pcgrr.
// Outputs are patient-level clinical interpretations: keep ${params.pcgr_dir} on the cluster.
// assay / target size come in as inputs: a module sees params as they were at include time, before main.nf resolves them.
process PCGR {
    tag "${patient_sample}"
    conda "${params.pcgr_env}"
    scratch true
    label 'process_medium'
    errorStrategy 'retry'
    maxRetries 3
    publishDir("${params.pcgr_dir}", mode: 'copy')

    input:
    tuple val(patient_sample), val(mode), val(tumor_site), val(assay), val(target_size_mb),
          path(sites_vcf, stageAs: 'input/sites/*'),
          path(mutect2_vcf, stageAs: 'input/mutect2/*'),
          path(sage_vcf, stageAs: 'input/sage/*'),
          path(strelka_snv, stageAs: 'input/strelka_snv/*'),
          path(strelka_indel, stageAs: 'input/strelka_indel/*'),
          path(cna, stageAs: 'input/cna/*')

    output:
    tuple val(patient_sample), path("${patient_sample}"), emit: report_dir

    script:
    def assembly = params.genome == 'GRCh37' ? 'grch37' : 'grch38'
    def paired = mode == 'paired'
    def callerArgs = paired ?
        "--mutect2 ${mutect2_vcf} --sage ${sage_vcf} --strelka_snv ${strelka_snv} --strelka_indel ${strelka_indel}" : ''
    def tagArgs = paired ?
        '--tumor_dp_tag TDP --tumor_af_tag TAF --control_dp_tag NDP --control_af_tag NAF' :
        '--tumor_dp_tag TDP --tumor_af_tag TAF'
    def tonlyArgs = paired ? '' :
        ("--tumor_only --gnomad_popmax_af_tolerated ${params.pcgr_gnomad_popmax_af}" +
         (params.pcgr_tonly_germline_filters ? ' --exclude_likely_hom_germline --exclude_likely_het_germline' : '') +
         (params.pcgr_pon_vcf ? " --pon_vcf ${params.pcgr_pon_vcf} --exclude_pon" : ''))
    // PCGR switches MSI off itself for TARGETED or tumour-only; only ask where it can run.
    def msiArg = (paired && assay != 'TARGETED' && params.pcgr_msi) ? '--estimate_msi' : ''
    def sigArg = params.pcgr_signatures ? "--estimate_signatures --min_mutations_signatures ${params.pcgr_min_mutations_signatures}" : ''
    def filterArgs = [
        params.pcgr_tumor_dp_min ? "--tumor_dp_min ${params.pcgr_tumor_dp_min}" : '',
        params.pcgr_tumor_af_min ? "--tumor_af_min ${params.pcgr_tumor_af_min}" : '',
        params.pcgr_tmb_dp_min ? "--tmb_dp_min ${params.pcgr_tmb_dp_min}" : '',
        params.pcgr_tmb_af_min ? "--tmb_af_min ${params.pcgr_tmb_af_min}" : ''
    ].findAll { it }.join(' ')
    """
    set -euo pipefail
    export XDG_CACHE_HOME=\$PWD/.cache   # quarto (HTML report) needs a writable cache
    # Keep the pinned envs isolated: a user-site pandas 1.x (~/.local/lib/python3.12) shadows the env's pandas and
    # pcgr dies with "No module named 'pytz'"; user R libraries could shadow pcgrr's the same way.
    export PYTHONNOUSERSITE=1 R_LIBS_USER=\$PWD/.no_user_rlib R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
    unset PYTHONPATH PERL5LIB R_LIBS
    mkdir -p "\$XDG_CACHE_HOME" ${patient_sample}

    python3 ${projectDir}/Modules/PCGR/pcgr_prepare_vcf.py --mode ${mode} \\
        --sites ${sites_vcf} --fai ${params.ref_fai} --out ${patient_sample}.pcgr_input.vcf ${callerArgs} \\
        | tee ${patient_sample}/${patient_sample}.pcgr_input_summary.txt
    bgzip -f ${patient_sample}.pcgr_input.vcf
    tabix -f -p vcf ${patient_sample}.pcgr_input.vcf.gz
    cp ${patient_sample}.pcgr_input.vcf.gz ${patient_sample}.pcgr_input.vcf.gz.tbi ${patient_sample}/

    cna_arg=""
    if [[ -s "${cna}" ]]; then
        python3 ${projectDir}/Modules/PCGR/purecn_to_pcgr_cna.py ${cna} ${patient_sample}/${patient_sample}.pcgr_cna.tsv \\
            | tee -a ${patient_sample}/${patient_sample}.pcgr_input_summary.txt
        cna_arg="--input_cna ${patient_sample}/${patient_sample}.pcgr_cna.tsv"
    fi

    n=\$(bgzip -dc ${patient_sample}.pcgr_input.vcf.gz | grep -vc '^#' || true)
    if [[ "\$n" -eq 0 && -z "\$cna_arg" ]]; then
        echo "No PASS SNV/indel on primary chromosomes and no CNA input: PCGR not run." \\
            > ${patient_sample}/${patient_sample}.pcgr_skipped.txt
        exit 0
    fi
    vcf_arg="--input_vcf ${patient_sample}.pcgr_input.vcf.gz"
    [[ "\$n" -eq 0 ]] && vcf_arg=""

    pcgr --version
    pcgr \\
        \$vcf_arg \$cna_arg \\
        --vep_dir ${params.pcgr_vep_dir} \\
        --refdata_dir ${params.pcgr_refdata} \\
        --output_dir ${patient_sample} \\
        --genome_assembly ${assembly} \\
        --sample_id ${patient_sample} \\
        --assay ${assay} \\
        --effective_target_size_mb ${target_size_mb} \\
        --tumor_site ${tumor_site} \\
        ${tagArgs} ${tonlyArgs} ${filterArgs} \\
        --estimate_tmb ${msiArg} ${sigArg} \\
        --retained_info_tags TVAF_SRC,CALLERS \\
        --vep_n_forks ${task.cpus} --vcfanno_n_proc ${task.cpus} \\
        ${params.pcgr_extra_args} \\
        --force_overwrite
    """
}
