# Turn a MERGE_HEADER "standard" consensus VCF (POST/standard_vcf/*_standard_PASSed.vcf) into a clean
# sites-only VCF that VEP and vcf2maf can read.
#
# Why: in those files the SNV rows have 8 columns with QUAL/FILTER "." and the per-caller INFO strings joined
# with "&", while the indel rows have 9 columns (QUAL holds the caller list, e.g. "mt,sa,st", FILTER is PASS,
# INFO is "." and a 9th column holds the "|"-joined caller INFO). The header declares 8 columns. Every record
# in the file is a consensus call (the file is the scenario's PASSed set), so FILTER is written as PASS.
# The original standard VCF is left untouched; only the caller list is carried over (INFO/callers).
BEGIN { FS = OFS = "\t"; header_done = 0 }
/^##fileformat=/ { if (!seen_ff++) print; next }
/^##contig=/     { print; next }
/^##/            { next }
/^#CHROM/ {
    if (!header_done) {
        print "##INFO=<ID=callers,Number=.,Type=String,Description=\"Callers that called the variant (from the standard consensus VCF)\">"
        print "##FILTER=<ID=PASS,Description=\"Consensus call (SomaticBunny POST scenario)\">"
        print "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
        header_done = 1
    }
    next
}
NF >= 5 {
    filter = (NF >= 7 && $7 != "" && $7 != ".") ? $7 : "PASS"
    info = "."
    if ($6 != "" && $6 != "." && $6 !~ /^[0-9.eE+-]+$/) info = "callers=" $6
    print $1, $2, ($3 == "" ? "." : $3), $4, $5, ".", filter, info
}
