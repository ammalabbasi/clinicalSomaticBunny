# Clinical SomaticBunny: write a valid, sites-only VCF from a MERGE_HEADER "standard" consensus VCF.
# The original file is not modified; this produces an additional *.fixed.vcf.gz (sorted with bcftools afterwards).
#
# What the original rows look like (from scenarioN.sh):
#   SNV rows   (8 columns): QUAL ".", FILTER ".", INFO = <caller INFO, joined with "&"> + ";callers=<list>;caller_filter=<FILTER>"
#              The FILTER copied into caller_filter may itself contain ";" (e.g. "maxGermlineVAF;LowEVS"), which a
#              VCF parser reads as separate, undeclared INFO keys.
#   Indel rows (9 columns): QUAL = <caller list>, FILTER = PASS, INFO = ".", column 9 = ".|<INFO 1>|<INFO 2>|..."
#   Header lines may also appear after #CHROM (the first caller VCF's header, concatenated in).
# What the fixed rows look like (8 columns, every record is a consensus call):
#   QUAL ".", FILTER PASS, INFO = callers=<list>;caller_filter=<list>;SB_CALLER_INFO=<original caller INFO, %-encoded>
# Usage: awk -v provfile=<file with ## provenance lines> -f fix_standard_vcf.awk <standard.vcf> > fixed.unsorted.vcf
function enc(s) {
    gsub(/%/, "%25", s); gsub(/;/, "%3B", s); gsub(/=/, "%3D", s); gsub(/,/, "%2C", s)
    gsub(/ /, "%20", s); gsub(/\t/, "%09", s)
    return s
}
function listval(s) {                       # caller/filter lists: ";" -> "," ; drop characters invalid in VCF
    gsub(/;/, ",", s); gsub(/[ =]/, "_", s)
    return (s == "") ? "." : s
}
BEGIN { FS = OFS = "\t"; done = 0; n_snv = 0; n_indel = 0; n_skip = 0 }
!done && /^##fileformat=/ { if (!ff++) print "##fileformat=VCFv4.2"; next }
!done && /^##INFO=<ID=(callers|caller_filter|caller_filtering_status|SB_CALLER_INFO),/ { next }
!done && /^##/ { print; next }
!done && /^#CHROM/ {
    print "##INFO=<ID=callers,Number=.,Type=String,Description=\"Callers that called the variant (mt=Mutect2, st=Strelka2, ms=MuSE, sa=SAGE)\">"
    print "##INFO=<ID=caller_filter,Number=.,Type=String,Description=\"Filter status carried by the consensus record\">"
    print "##INFO=<ID=SB_CALLER_INFO,Number=1,Type=String,Description=\"Original per-caller INFO from the standard VCF, percent-encoded (%3B ; %3D = %2C , %25 % %20 space); SNVs joined by &, indels by |\">"
    print "##FILTER=<ID=PASS,Description=\"Consensus call (SomaticBunny POST scenario)\">"
    if (provfile != "") while ((getline line < provfile) > 0) print line
    print "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
    done = 1
    next
}
/^#/ { next }                                # stray header lines after #CHROM
NF < 8 { n_skip++; next }
NF == 8 {                                    # SNV row
    info = $8; callers = "."; cf = "."; rest = info
    p = index(info, "caller_filter=")
    if (p) {
        cf = substr(info, p + 14); rest = substr(info, 1, p - 1)
        sub(/;$/, "", rest)
    }
    q = index(rest, "callers=")
    if (q && (q == 1 || substr(rest, q - 1, 1) == ";")) {
        callers = substr(rest, q + 8); rest = substr(rest, 1, q - 1)
        sub(/;$/, "", rest)
    }
    out = "callers=" listval(callers) ";caller_filter=" listval(cf)
    if (rest != "" && rest != ".") out = out ";SB_CALLER_INFO=" enc(rest)
    print $1, $2, ($3 == "" ? "." : $3), $4, $5, ".", "PASS", out
    n_snv++
    next
}
{                                            # indel row (9+ columns)
    rest = $9; sub(/^\.\|/, "", rest)
    if ($8 != "." && $8 != "") rest = $8 "|" rest
    out = "callers=" listval($6) ";caller_filter=" listval($7)
    if (rest != "" && rest != ".") out = out ";SB_CALLER_INFO=" enc(rest)
    print $1, $2, ($3 == "" ? "." : $3), $4, $5, ".", "PASS", out
    n_indel++
}
END { printf "fix_standard_vcf: %d SNV rows, %d indel rows, %d short rows skipped\n", n_snv, n_indel, n_skip > "/dev/stderr" }
