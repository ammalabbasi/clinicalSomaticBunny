nextflow.enable.dsl=2

// Microsatellite list for MSIsensor-pro, once per reference/assay. Targeted runs scan only the padded targets
// (extracted with samtools faidx, coordinates shifted back to the reference); genome runs scan all of --ref.
process MSI_SCAN {
    conda "${params.msi_env}"
    scratch true
    label 'process_medium'
    publishDir("${params.msi_dir}/sites", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    // Padded targets BED, or the string 'genome' for a whole-reference scan
    val(regions)

    output:
    path("msi_sites.tsv"), emit: sites

    script:
    if (regions.toString() != 'genome') {
        """
        set -euo pipefail
        # 1-based faidx regions from the 0-based BED
        awk -v OFS='' 'NF>=3 { print \$1, ":", \$2 + 1, "-", \$3 }' ${regions} > regions.txt
        samtools faidx -r regions.txt ${params.ref} > targets.fa
        msisensor-pro scan -d targets.fa -o scan.raw.tsv

        # Back to reference coordinates: sequence "chr:S-E", position p -> p + S - 1. Drop sites within 10 bp of a
        # region edge (the scan sees them truncated, with the wrong flanks). Keep the reference's contig order.
        awk -F'\\t' -v OFS='\\t' 'FNR==NR { rank[\$1] = FNR; next }
            FNR==1 { print > "header.tsv"; next }
            { n = split(\$1, a, ":"); split(a[n], b, "-"); chr = substr(\$1, 1, length(\$1) - length(a[n]) - 1)
              p = \$2; len = b[2] - b[1] + 1; span = \$3 * \$4
              if (p < 10 || p + span > len - 10) next
              \$1 = chr; \$2 = p + b[1] - 1; print rank[chr], \$0 }' ${params.ref_fai} scan.raw.tsv \\
          | sort -k1,1n -k3,3n | cut -f2- > body.tsv
        cat header.tsv body.tsv > msi_sites.tsv
        echo "MSI_SCAN: \$(wc -l < regions.txt) regions -> \$(wc -l < body.tsv) microsatellites" >&2
        if [ ! -s body.tsv ]; then
            echo "ERROR: no microsatellites found in the targets" >&2
            exit 1
        fi
        """
    } else {
        """
        set -euo pipefail
        msisensor-pro scan -d ${params.ref} -o msi_sites.tsv
        """
    }
}
