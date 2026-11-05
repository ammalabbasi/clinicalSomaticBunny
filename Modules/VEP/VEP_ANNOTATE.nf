nextflow.enable.dsl=2

// Ensembl VEP, offline cache, on one consensus VCF:
//   paired     : POST/standard_vcf/<patient_sample>_<vep_scenario>_standard_PASSed.vcf (sites-only; cleaned first)
//   tumor_only : POST/tumor_only/<patient_sample>/<patient_sample>_tumor_only_consensus.PASS.vcf.gz
// Flags are the ones vcf2maf 1.6.22 passes to VEP (so vcf2maf --inhibit-vep reads this output as if it had run
// VEP itself), plus --mane and a pick order that prefers the MANE Select transcript, then Ensembl canonical.
// Multi-allelic records are split first, so the MAF has one row per ALT allele.
process VEP_ANNOTATE {
    conda "${params.vep_env}"
    scratch true
    label 'process_low'
    publishDir("${params.vep_dir}/${patient_sample}", mode: 'copy', pattern: "*.vep.vcf.gz*")
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(patient_sample), val(mode), val(tumor_id), val(normal_id), path(vcf, stageAs: 'input/*')

    output:
    tuple val(patient_sample), val(mode), val(tumor_id), val(normal_id), path("${patient_sample}.${mode}.vep.vcf.gz"), emit: vcf
    path("${patient_sample}.${mode}.vep.vcf.gz.tbi"), emit: tbi

    script:
    def assembly = params.vep_assembly ?: params.genome
    def fasta = params.vep_fasta ?: params.ref
    def prefix = "${patient_sample}.${mode}"
    def prep = (mode == 'paired') ?
        "awk -f ${moduleDir}/standard_vcf_to_sites.awk ${vcf} > raw.vcf" :
        "bcftools view -Ov -o raw.vcf ${vcf}"
    """
    set -euo pipefail
    ${prep}
    bcftools norm -m -any -Ou raw.vcf | bcftools sort -Ov -o sites.vcf -
    n_in=\$(grep -vc '^#' sites.vcf || true)
    echo "${prefix}: \$n_in records to annotate" >&2

    if [[ "\$n_in" == "0" ]]; then
        # Nothing called: keep an empty, valid VCF so the MAF step still writes a header-only MAF.
        cp sites.vcf ${prefix}.vep.vcf
    else
    vep --species homo_sapiens --assembly ${assembly} \\
        --offline --dir ${params.vep_cache} --cache_version ${params.vep_cache_version} \\
        --fasta ${fasta} --format vcf --vcf \\
        --input_file sites.vcf --output_file ${prefix}.vep.vcf --force_overwrite \\
        --no_progress --no_stats --buffer_size 5000 --fork ${task.cpus} \\
        --sift b --polyphen b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype \\
        --canonical --mane --protein --biotype --tsl --variant_class --shift_hgvs 1 \\
        --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 \\
        --flag_pick_allele --pick_order mane_select,canonical,tsl,biotype,rank,ccds,length \\
        --af --af_1kg --af_gnomade --af_gnomadg --max_af --pubmed --regulatory
    fi

    # VEP writes every input record; anything else means records were lost.
    n_out=\$(grep -vc '^#' ${prefix}.vep.vcf || true)
    if [[ "\$n_in" != "\$n_out" ]]; then
        echo "ERROR: ${prefix}: \$n_in records in, \$n_out annotated" >&2
        exit 1
    fi
    bgzip -f ${prefix}.vep.vcf
    tabix -p vcf ${prefix}.vep.vcf.gz
    """
}
