nextflow.enable.dsl=2

// vcf2maf on a VEP-annotated consensus VCF (VEP already ran in VEP_ANNOTATE, so --inhibit-vep).
// Sample IDs follow the caller VCF columns: tumour <patient>_<sample>_tumor, normal <patient>_<normal>_normal.
//   paired     : the standard consensus VCF is sites-only, so t_/n_ depth and count columns stay empty.
//   tumor_only : Mutect2's tumour column is read for depths; Matched_Norm_Sample_Barcode is vcf2maf's
//                default "NORMAL" and the normal columns stay empty (vcf2maf's tumour-only convention).
// One MAF row per VCF record; the task fails if the counts differ.
process VCF2MAF {
    conda "${params.vep_env}"
    scratch true
    label 'process_low'
    publishDir("${params.vep_maf_dir}/${patient_sample}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(patient_sample), val(mode), val(tumor_id), val(normal_id), path(vep_vcf, stageAs: 'input/*')

    output:
    tuple val(patient_sample), val(mode), path("${patient_sample}.${mode}.maf"), emit: maf

    script:
    def assembly = params.vep_assembly ?: params.genome
    def fasta = params.vep_fasta ?: params.ref
    def prefix = "${patient_sample}.${mode}"
    """
    set -euo pipefail
    bgzip -dc ${vep_vcf} > ${prefix}.vep.vcf
    n_vcf=\$(grep -vc '^#' ${prefix}.vep.vcf || true)

    # Tumour column: the expected caller name, else the only sample column (tumour-only VCFs have one).
    vcf_tumor="${tumor_id}"
    samples=\$(bcftools query -l ${prefix}.vep.vcf | tr '\\n' ' ')
    if [[ -n "\$samples" ]] && ! grep -qw -- "${tumor_id}" <<< "\$samples"; then
        set -- \$samples
        if [[ \$# -eq 1 ]]; then
            vcf_tumor=\$1
            echo "WARNING: ${prefix}: no column ${tumor_id}; using the only sample column \$vcf_tumor" >&2
        else
            echo "ERROR: ${prefix}: no column ${tumor_id} among: \$samples" >&2
            exit 1
        fi
    fi

    mkdir -p tmp
    vcf2maf.pl --inhibit-vep \\
        --input-vcf ${prefix}.vep.vcf --output-maf ${prefix}.maf --tmp-dir tmp \\
        --tumor-id ${tumor_id} --vcf-tumor-id "\$vcf_tumor" \\
        --normal-id ${normal_id} --vcf-normal-id ${normal_id} \\
        --ref-fasta ${fasta} --ncbi-build ${assembly} --species homo_sapiens \\
        --cache-version ${params.vep_cache_version} \\
        --retain-ann MANE_SELECT,MANE_PLUS_CLINICAL

    n_maf=\$(grep -v '^#' ${prefix}.maf | tail -n +2 | wc -l)
    echo "${prefix}: \$n_vcf VCF records, \$n_maf MAF rows" >&2
    if [[ "\$n_vcf" != "\$n_maf" ]]; then
        echo "ERROR: ${prefix}: MAF has \$n_maf rows for \$n_vcf VCF records" >&2
        exit 1
    fi
    """
}
