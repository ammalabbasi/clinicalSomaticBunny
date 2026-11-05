nextflow.enable.dsl=2

process CNVkit {
    conda "${params.cnvkit_env}"
    scratch true
    label 'process_high'
    publishDir("${params.cnvkit_dir}", mode: 'copy')

    input:
    val map
    path reference_cnn

    output:
    path("*.bed"), emit: CNVkit_bed
    path("*.cn*"), emit: CNVkit_cn_files
    path("*.pdf"), emit: CNVkit_pdf
    path("*.png"), emit: CNVkit_png

    script:
    def output_prefix = "${map.patient}_${map.sample}"
    // The method must match the reference: a targeted flat reference is built on-target only (amplicon).
    def method_option = !params.targeted ? "--method wgs" : (params.cnvkit_ref_mode == 'flat' ? "--method amplicon" : "")

    """
    export PYTHONNOUSERSITE=1
    unset PYTHONPATH

    cnvkit.py batch \
        "${map.tumor}" \
        -r "${reference_cnn}" \
        ${method_option} \
        -p ${task.cpus} \
        --scatter \
        --diagram

    # Find the segment file CNVkit wrote: <base>.cns
    shopt -s nullglob
    segs=()
    for f in *.cns; do
        b="\${f%.cns}"
        if [ -e "\$b.call.cns" ] && [ -e "\$b.bintest.cns" ]; then
            segs+=("\$f")
        fi
    done
    if [ "\${#segs[@]}" -ne 1 ]; then
        echo "ERROR: expected 1 segment file from cnvkit.py batch, found \${#segs[@]}: \${segs[*]:-none}" >&2
        ls -l >&2
        exit 1
    fi
    test -s "\${segs[0]}"

    cnvkit.py call \
        "\${segs[0]}" \
        -o "${output_prefix}_calls.cns"

    cnvkit.py export bed \
        "${output_prefix}_calls.cns" \
        -o "${output_prefix}_calls.bed"
    """
}