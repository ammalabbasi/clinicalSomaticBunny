nextflow.enable.dsl=2

// Concatenate the per-sample MAFs into MAF/cohort.maf, adding Analysis_Mode (paired | tumor_only) as the last
// column. Per-sample MAFs are named <patient_sample>.<mode>.maf. Stops if the column headers differ.
process VEP_COHORT_MAF {
    conda "${params.vep_env}"
    scratch true
    label 'process_low'
    publishDir("${params.vep_maf_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    path(mafs, stageAs: 'input/*')

    output:
    path("cohort.maf"), emit: maf

    script:
    """
    set -euo pipefail
    first=1
    header=""
    for f in \$(ls input/*.maf | sort); do
        case "\$f" in
            *.tumor_only.maf) mode=tumor_only ;;
            *.paired.maf)     mode=paired ;;
            *) echo "ERROR: cannot tell the mode of \$f" >&2; exit 1 ;;
        esac
        # awk, not grep | head: under pipefail, head closing the pipe early kills grep with SIGPIPE (exit 141)
        h=\$(awk '!/^#/ { print; exit }' "\$f")
        if [[ \$first -eq 1 ]]; then
            awk '/^#/ { print; exit }' "\$f" > cohort.maf
            printf '%s\\tAnalysis_Mode\\n' "\$h" >> cohort.maf
            header="\$h"
            first=0
        elif [[ "\$h" != "\$header" ]]; then
            echo "ERROR: MAF columns of \$f differ from the first MAF" >&2
            exit 1
        fi
        grep -v '^#' "\$f" | tail -n +2 | awk -v m="\$mode" 'BEGIN{OFS="\\t"} {print \$0, m}' >> cohort.maf
    done
    """
}
