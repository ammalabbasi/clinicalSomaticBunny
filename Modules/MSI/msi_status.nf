nextflow.enable.dsl=2

// Shell snippet shared by MSI_PAIRED and MSI_TUMOR_ONLY: turn the MSIsensor-pro summary
// (header "Total_Number_of_Sites  Number_of_Unstable_Sites  %", then one row) into one status row.
// status: MSI-H when score >= threshold, else MSS. qc: LOW_SITES when fewer than --msi_min_sites were evaluated
// (the status is then not reliable), else PASS.
def msiStatusCmd(String summary, String out, String patient, String sample, String mode, threshold) {
    return """
    awk -F'\\t' -v OFS='\\t' -v pt='${patient}' -v sm='${sample}' -v mode='${mode}' \\
        -v thr='${threshold}' -v minsites='${params.msi_min_sites}' '
        NR == 2 { total = \$1 + 0; unstable = \$2 + 0; score = \$3 + 0; found = 1 }
        END {
            if (!found) { print "ERROR: empty MSIsensor-pro summary ${summary}" > "/dev/stderr"; exit 1 }
            status = (score >= thr + 0) ? "MSI-H" : "MSS"
            qc = (total < minsites + 0) ? "LOW_SITES" : "PASS"
            print "patient", "sample", "mode", "sites_evaluated", "unstable_sites", "msi_score_pct", "threshold_pct", "status", "qc"
            print pt, sm, mode, total, unstable, score, thr, status, qc
        }' ${summary} > ${out}
    cat ${out} >&2
    """
}
