
process RENAME_BAM_HEADER {
    tag "${meta.patient}_${meta.sample}"
    label 'process_low'
    conda "${params.samtools_env}"

    input:
    // Staged under input/ so a BAM already named <patient>_<sample>_<status>.bam cannot collide with the output.
    tuple val(meta), path(bam, stageAs: 'input/*'), path(bai, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.patient}_${meta.sample}_${meta.status}.bam"), path("${meta.patient}_${meta.sample}_${meta.status}.bam.bai"), emit: renamed_bam

    script:
    def new_name   = "${meta.patient}_${meta.sample}_${meta.status}"
    def output_bam = "${new_name}.bam"
    """
    samtools view -H --no-PG ${bam} > orig.hdr

    # Clinical: contig names. Vendor BAMs are often b37 style (1, 2, MT) while the pipeline database is chr-named.
    # With --contig_rename auto, rename a contig only when the other style exists in --ref with the SAME length.
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
        { print }' ${params.ref_fai} orig.hdr > contig.hdr

    # Clinical: fail fast, with a readable message, if the BAM does not match --ref: primary chromosomes must exist
    # with the same lengths, and contigs shared with --ref must come in the same order (GATK checks the order).
    awk -F'\\t' 'NR==FNR { idx[\$1]=FNR; if (\$1 ~ /^(chr)?([0-9]+|X)\$/) want[\$1]=\$2; next }
        /^@SQ/ { for (i=2; i<=NF; i++) { if (\$i ~ /^SN:/) n=substr(\$i,4); if (\$i ~ /^LN:/) l=substr(\$i,4) }; have[n]=l
                 if (n in idx) { if (idx[n] < last) order_bad=1; last=idx[n] } else extra++ }
        END { for (c in want) if (!(c in have) || have[c]!=want[c]) { bad++; msg=msg " " c "(ref " want[c] ", bam " (c in have ? have[c] : "absent") ")" }
              if (bad) { print "ERROR: BAM contigs do not match --ref:" msg > "/dev/stderr"; exit 1 }
              if (order_bad) { print "ERROR: BAM contig order differs from --ref; use a --ref FASTA in the same contig order as the BAM (or picard ReorderSam)" > "/dev/stderr"; exit 1 }
              if (extra) print "WARNING: " extra " BAM contig(s) are not in --ref; reads on them are ignored or may upset GATK" > "/dev/stderr" }' ${params.ref_fai} contig.hdr
    if ! cmp -s orig.hdr contig.hdr; then echo "renamed contigs to match --ref" >&2; fi

    n_rg=\$(grep -c '^@RG' contig.hdr || true)

    if [ "\$n_rg" -gt 0 ]; then
        # Upstream: keep the existing read groups (ID, LB, PU, PL per lane/library); only SM becomes the pipeline sample name
        awk -F'\\t' -v sm='${new_name}' 'BEGIN{OFS="\\t"}
            /^@RG/ { has=0; for (i=2;i<=NF;i++) if (\$i ~ /^SM:/) { \$i="SM:" sm; has=1 }
                     if (!has) \$0 = \$0 OFS "SM:" sm }
            { print }' contig.hdr > new_header.sam
        samtools reheader new_header.sam ${bam} > ${output_bam}
    else
        # Upstream: no read group in the input: add a single one (after the contig rename, if any)
        IN=${bam}
        if ! cmp -s orig.hdr contig.hdr; then
            samtools reheader contig.hdr ${bam} > reheadered.bam
            IN=reheadered.bam
        fi
        samtools addreplacerg -@ ${task.cpus} -m overwrite_all \\
            -r ID:${new_name} -r SM:${new_name} -r LB:${new_name} -r PL:ILLUMINA -r PU:${new_name} \\
            -o ${output_bam} \$IN
        rm -f reheadered.bam
    fi

    samtools index ${output_bam}

    # Guard against silent output corruption (a bad node once wrote a truncated BAM with exit 0,
    # and -resume kept reusing it): quickcheck + require the output to hold as many reads as the input.
    samtools quickcheck ${output_bam}
    n_in=\$(samtools view -c ${bam})
    n_out=\$(samtools view -c ${output_bam})
    [ "\$n_in" -eq "\$n_out" ] || { echo "ERROR: ${output_bam} has \$n_out reads, input had \$n_in" >&2; exit 1; }
    """
}
