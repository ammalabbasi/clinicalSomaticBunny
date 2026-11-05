# Turn `bcftools mpileup -a FORMAT/AD` output into ref/alt allele counts at GIVEN alleles, for PureCN.
#   usage: awk -f purecn_counts.awk alleles.tsv mpileup.vcf
#   alleles.tsv: CHROM POS REF,ALT (the gnomAD SNP sites + Mutect2 PASS SNVs to report)
# For each mpileup record at a listed site whose REF matches, report AD = (ref reads, reads of the listed ALT) and
# DP = their sum, whatever other alleles mpileup saw. Sites without a single alt read in any sample are dropped.
# A genotyper (bcftools call -C alleles) is NOT used on purpose: it sets ALT to "." and drops the alt count when it
# calls 0/0 in every sample, which is what happens to a low-VAF somatic SNV (seen on synthetic data).
# GT is a plain allele-fraction label (0/0 no alt read, 1/1 alt fraction >= 0.9, else 0/1); PureCN only uses the
# counts, and finds the tumour column by --sampleid.
BEGIN { FS = OFS = "\t" }
NR == FNR {
    split($3, given, ",")
    want_ref[$1 ":" $2] = given[1]
    want_alt[$1 ":" $2] = given[2]
    next
}
/^##/ {
    if ($0 ~ /^##(INFO|FORMAT)=/) next
    print
    next
}
/^#CHROM/ {
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype label from the alt allele fraction (no genotyping)\">"
    print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Reads supporting REF and ALT (bcftools mpileup, given alleles)\">"
    print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"REF + ALT reads\">"
    print
    next
}
{
    k = $1 ":" $2
    if (!(k in want_alt) || toupper($4) != want_ref[k]) next
    n = split($5, alts, ",")
    idx = 0
    for (i = 1; i <= n; i++) if (alts[i] == want_alt[k]) idx = i
    nf = split($9, fmt, ":")
    iad = 0
    for (i = 1; i <= nf; i++) if (fmt[i] == "AD") iad = i
    if (!iad) next
    out = ""
    any = 0
    for (s = 10; s <= NF; s++) {
        split($s, v, ":")
        split(v[iad], ad, ",")
        r = ad[1] + 0
        a = idx ? ad[idx + 1] + 0 : 0
        d = r + a
        gt = (a == 0) ? "0/0" : ((a >= 0.9 * d) ? "1/1" : "0/1")
        out = out OFS gt ":" r "," a ":" d
        if (a > 0) any = 1
    }
    if (!any) next
    print $1, $2, ".", want_ref[k], want_alt[k], ".", ".", ".", "GT:AD:DP" out
}
