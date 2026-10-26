nextflow.enable.dsl=2

// Pad, sort and merge the cleaned --targets BED (main.nf cleans and checks it at startup); bgzip + tabix the
// result for Strelka --callRegions.
process PREPARE_TARGETS {
    conda "${params.samtools_env}"
    label 'process_low'
    publishDir("${params.intervals_dir}", mode: 'copy')

    input:
    // Staged under input/: writing a work file with the input's name would truncate the original through
    // the staged symlink.
    path(targets, stageAs: 'input/*')
    path(ref_fai)

    output:
    tuple path("targets.padded.bed.gz"), path("targets.padded.bed.gz.tbi"), emit: call_regions
    path("targets.padded.bed"), emit: padded_bed

    script:
    def cat_cmd = targets.name.endsWith('.gz') ? "zcat" : "cat"
    """
    ${cat_cmd} ${targets} | grep -vE '^(track|browser|#)' | awk -v OFS='\\t' 'NF>=3{print \$1,\$2,\$3}' > work.targets.bed
    if [ ! -s work.targets.bed ]; then
        echo "ERROR: --targets has no BED intervals" >&2
        exit 1
    fi

    # Every target contig must exist in the reference (chr1 vs 1 is the usual mismatch).
    awk 'NR==FNR{ok[\$1]=1; next} !(\$1 in ok){print \$1}' ${ref_fai} work.targets.bed | sort -u > missing_contigs.txt
    if [ -s missing_contigs.txt ]; then
        echo "ERROR: --targets uses contigs that are not in --ref: \$(head -5 missing_contigs.txt | tr '\\n' ' ')" >&2
        echo "       (chr prefix mismatch? the BED must use the same contig names as the reference)" >&2
        exit 1
    fi

    awk -v OFS='\\t' -v pad=${params.target_padding} 'NR==FNR{len[\$1]=\$2; next}
        { s=\$2-pad; if (s<0) s=0; e=\$3+pad; if (e>len[\$1]) e=len[\$1]; print \$1, s, e }' ${ref_fai} work.targets.bed \\
      | sort -k1,1 -k2,2n \\
      | awk -v OFS='\\t' 'c==\$1 && \$2<=e { if (\$3>e) e=\$3; next }
                         { if (c!="") print c, s, e; c=\$1; s=\$2; e=\$3 }
                         END { if (c!="") print c, s, e }' > targets.padded.bed
    bgzip -c targets.padded.bed > targets.padded.bed.gz
    tabix -p bed targets.padded.bed.gz
    echo "targets: \$(wc -l < work.targets.bed) input intervals -> \$(wc -l < targets.padded.bed) padded/merged, \$(awk '{s+=\$3-\$2}END{print s}' targets.padded.bed) bp"
    """
}
