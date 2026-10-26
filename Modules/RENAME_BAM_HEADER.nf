process RENAME_BAM_HEADER {
    tag "${meta.patient}_${meta.sample}"
    label 'process_low'
    conda "${params.rename_bam_env}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.patient}_${meta.sample}_${meta.status}.bam"), path("${meta.patient}_${meta.sample}_${meta.status}.bam.bai"), emit: renamed_bam

    script:
    def new_name = "${meta.patient}_${meta.sample}_${meta.status}"
    def output_bam = "${new_name}.bam"
    """
    # Contig names: vendor BAMs are often b37 style (1, 2, MT) while the pipeline database is chr-named.
    # With --contig_rename auto, rename a contig only when the other style exists in --ref with the SAME length.
    samtools view -H ${bam} > orig.hdr
    awk -F'\\t' -v OFS='\\t' -v mode=${params.contig_rename} '
        NR==FNR { len[\$1]=\$2; next }
        /^@SQ/ && mode=="auto" {
            for (i=2; i<=NF; i++) { if (\$i ~ /^SN:/) si=i; if (\$i ~ /^LN:/) l=substr(\$i,4) }
            n=substr(\$si,4)
            if (!(n in len)) {
                if (n=="MT") c="chrM"; else if (n=="chrM") c="MT"; else if (n ~ /^chr/) c=substr(n,4); else c="chr" n
                if ((c in len) && len[c]==l) \$si="SN:" c
            }
        }
        { print }' ${params.ref_fai} orig.hdr > new.hdr

    # Fail fast, with a readable message, if the BAM does not match --ref: primary chromosomes must exist with
    # the same lengths, and contigs shared with --ref must come in the same order (GATK checks the order).
    awk -F'\\t' 'NR==FNR { idx[\$1]=FNR; if (\$1 ~ /^(chr)?([0-9]+|X)\$/) want[\$1]=\$2; next }
        /^@SQ/ { for (i=2; i<=NF; i++) { if (\$i ~ /^SN:/) n=substr(\$i,4); if (\$i ~ /^LN:/) l=substr(\$i,4) }; have[n]=l
                 if (n in idx) { if (idx[n] < last) order_bad=1; last=idx[n] } else extra++ }
        END { for (c in want) if (!(c in have) || have[c]!=want[c]) { bad++; msg=msg " " c "(ref " want[c] ", bam " (c in have ? have[c] : "absent") ")" }
              if (bad) { print "ERROR: BAM contigs do not match --ref:" msg > "/dev/stderr"; exit 1 }
              if (order_bad) { print "ERROR: BAM contig order differs from --ref; use a --ref FASTA in the same contig order as the BAM (or picard ReorderSam)" > "/dev/stderr"; exit 1 }
              if (extra) print "WARNING: " extra " BAM contig(s) are not in --ref; reads on them are ignored or may upset GATK" > "/dev/stderr" }' ${params.ref_fai} new.hdr

    IN=${bam}
    if ! cmp -s orig.hdr new.hdr; then
        echo "renamed contigs to match --ref" >&2
        samtools reheader new.hdr ${bam} > reheadered.bam
        IN=reheadered.bam
    fi

    picard -Xmx${task.memory.toGiga()}g AddOrReplaceReadGroups \\
        I=\$IN \\
        O=${output_bam} \\
        RGID=${new_name} \\
        RGSM=${meta.patient}_${meta.sample}_${meta.status} \\
        RGPL=ILLUMINA \\
        RGLB=${meta.sample} \\
        RGPU=${new_name} \\
        SORT_ORDER=coordinate \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT \\
        USE_JDK_DEFLATER=true \\
        USE_JDK_INFLATER=true

    mv ${new_name}.bai ${output_bam}.bai
    rm -f reheadered.bam
    """
}
