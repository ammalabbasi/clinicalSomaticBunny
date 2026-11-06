#!/usr/bin/env python3
"""Per-gene somatic status for one tumour sample (see Modules/GENE_STATUS/GENE_STATUS.nf for the column key).

Covered genes = genes BED (4 columns) ∩ targets BED; a gene with no targeted bp is omitted (NOT_COVERED).
Per covered gene:
  - callable_bp / callable_frac from the tumour's `mosdepth --quantize` output (CALLABLE = depth >= --min_depth),
  - n_snv / n_indel from the consensus VCF, assigned by the record's POS (anchor base) falling in the gene's
    targeted intervals; one count per ALT allele, duplicates (CHROM,POS,REF,ALT) counted once,
  - worst VEP consequence + SYMBOLs when the VCF carries a CSQ annotation (--tool vep),
  - copy number: CNVkit *_calls.cns `cn` and PureCN *_loh.csv `C`/`M`, each from the segment with the largest
    overlap with the gene; the status uses PureCN when present, else CNVkit (column cn_source),
  - status: NOT_CALLABLE when callable_frac < --min_callable; else '+'-joined MUTATED / AMP (cn>=5) /
    GAIN (cn>2) / DEL (cn<2) / DEEP_DEL (cn=0) / LOH (PureCN M=0 and C>0); WT when none apply.
"""
import argparse
import csv
import gzip
import re
import sys

# VEP consequence terms, most severe first (Ensembl severity order; unknown terms rank after these).
SEVERITY = [
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained",
    "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification",
    "feature_elongation", "feature_truncation", "inframe_insertion", "inframe_deletion",
    "missense_variant", "protein_altering_variant", "splice_donor_5th_base_variant",
    "splice_region_variant", "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant", "start_retained_variant", "stop_retained_variant",
    "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant",
    "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
    "NMD_transcript_variant", "non_coding_transcript_variant", "coding_transcript_variant",
    "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification",
    "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification",
    "regulatory_region_variant", "intergenic_variant", "sequence_variant",
]
SEV_RANK = {t: i for i, t in enumerate(SEVERITY)}


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def merge(ivs):
    out = []
    for s, e in sorted(ivs):
        if e <= s:
            continue
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def intersect_lists(a, b):
    """Intersection of two merged interval lists on one chromosome."""
    i = j = 0
    res = []
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            res.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return res


def total_bp(ivs):
    return sum(e - s for s, e in ivs)


def contains(ivs, pos0):
    lo, hi = 0, len(ivs)
    while lo < hi:
        mid = (lo + hi) // 2
        if ivs[mid][1] <= pos0:
            lo = mid + 1
        else:
            hi = mid
    return lo < len(ivs) and ivs[lo][0] <= pos0 < ivs[lo][1]


def read_bed(path, keep_label=None):
    """{chrom: merged [(start, end), ...]}. keep_label filters on column 4 (mosdepth quantize labels)."""
    iv = {}
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            if keep_label is not None and (len(f) < 4 or f[3] != keep_label):
                continue
            iv.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return {c: merge(v) for c, v in iv.items()}


def read_genes(path):
    """4-column BED -> ordered {(gene, chrom): merged [(start, end), ...]} keeping file order of first appearance."""
    genes = {}
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or not f[3].strip():
                continue
            genes.setdefault((f[3].strip(), f[0]), []).append((int(f[1]), int(f[2])))
    return {k: merge(v) for k, v in genes.items()}


def csq_parser(header_lines):
    """(symbol_idx, consequence_idx, pick_idx) of the CSQ Format, or None when the VCF has no CSQ annotation."""
    for h in header_lines:
        if h.startswith("##INFO=<ID=CSQ,"):
            m = re.search(r"Format: ([^\">]+)", h)
            if m:
                fields = m.group(1).strip().split("|")
                if "Consequence" in fields:
                    return (fields.index("SYMBOL") if "SYMBOL" in fields else None,
                            fields.index("Consequence"),
                            fields.index("PICK") if "PICK" in fields else None)
    return None


def info_get(info, key):
    """Value of INFO key up to the next ';' (consensus INFO may join callers' INFO with '&' or '|')."""
    m = re.search(r"(?:^|[;&|])" + re.escape(key) + r"=([^;]*)", info)
    return m.group(1) if m else ""


def pick_overlap(segments, gene_ivs):
    """segments: [(start, end, value)]; return the value with the largest total overlap with gene_ivs, or None."""
    weights = {}
    for s, e, v in segments:
        bp = total_bp(intersect_lists(gene_ivs, [(s, e)]))
        if bp:
            weights[v] = weights.get(v, 0) + bp
    return max(weights, key=weights.get) if weights else None


def read_cnvkit(path):
    """*_calls.cns -> {chrom: [(start, end, cn)]} (integer cn from `cnvkit.py call`)."""
    segs = {}
    with open_text(path) as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        cols = {c.lower(): c for c in rows.fieldnames or []}
        if any(n not in cols for n in ("chromosome", "start", "end", "cn")):
            sys.exit(f"--cnvkit_cns {path}: needs columns chromosome,start,end,cn (CNVkit *_calls.cns)")
        for r in rows:
            try:
                segs.setdefault(r[cols["chromosome"]], []).append(
                    (int(r[cols["start"]]), int(r[cols["end"]]), int(round(float(r[cols["cn"]])))))
            except (ValueError, TypeError):
                continue
    return segs


def read_purecn(path):
    """*_loh.csv -> {chrom: [(start, end, (C, M))]} (PureCN integer total / minor copy number segments)."""
    segs = {}
    with open_text(path) as fh:
        rows = csv.DictReader(fh)
        cols = {c.lower(): c for c in rows.fieldnames or []}
        if any(n not in cols for n in ("chr", "start", "end", "c")):
            sys.exit(f"--purecn_loh {path}: needs columns chr,start,end,C (PureCN *_loh.csv)")
        mcol = cols.get("m")
        for r in rows:
            try:
                c = int(round(float(r[cols["c"]])))
            except (ValueError, TypeError):
                continue
            m = None
            if mcol and r[mcol] not in (None, "", "NA"):
                try:
                    m = int(round(float(r[mcol])))
                except ValueError:
                    m = None
            try:
                segs.setdefault(r[cols["chr"]], []).append(
                    (int(float(r[cols["start"]])), int(float(r[cols["end"]])), (c, m)))
            except (ValueError, TypeError):
                continue
    return segs


def cn_label(c, amp_cn=5):
    if c is None:
        return None
    if c == 0:
        return "DEEP_DEL"
    if c < 2:
        return "DEL"
    if c >= amp_cn:
        return "AMP"
    if c > 2:
        return "GAIN"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--mode", required=True, choices=["paired", "tumor_only"])
    ap.add_argument("--source", required=True, help="label of the consensus the VCF came from")
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--genes_bed", required=True, help="4-column BED: chrom, start, end, gene name")
    ap.add_argument("--targets_bed", required=True, help="padded clean target BED")
    ap.add_argument("--quantized", required=True, help="mosdepth quantized.bed.gz of the tumour")
    ap.add_argument("--min_depth", type=int, required=True)
    ap.add_argument("--min_callable", type=float, required=True)
    ap.add_argument("--cnvkit_cns", help="CNVkit *_calls.cns of this tumour")
    ap.add_argument("--purecn_loh", help="PureCN *_loh.csv of this tumour")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    genes = read_genes(a.genes_bed)
    if not genes:
        sys.exit(f"--genes_bed {a.genes_bed}: no 4-column gene intervals")
    targets = read_bed(a.targets_bed)
    callable_iv = read_bed(a.quantized, keep_label="CALLABLE")

    # Covered genes and their targeted intervals (NOT_COVERED genes are dropped here)
    covered = {}   # (gene, chrom) -> dict
    for (gene, chrom), ivs in genes.items():
        cov = intersect_lists(ivs, targets.get(chrom, []))
        if not cov:
            continue
        cal = intersect_lists(cov, callable_iv.get(chrom, []))
        covered[(gene, chrom)] = {
            "gene": gene, "chrom": chrom, "start": ivs[0][0], "end": ivs[-1][1],
            "targeted": cov, "targeted_bp": total_bp(cov), "callable_bp": total_bp(cal),
            "n_snv": 0, "n_indel": 0, "worst": None, "symbols": set(),
        }
    if not covered:
        print(f"WARNING [{a.sample}]: no gene of {a.genes_bed} overlaps the targets", file=sys.stderr)

    by_chrom = {}
    for key, g in covered.items():
        by_chrom.setdefault(g["chrom"], []).append(key)

    # Consensus variants, assigned to genes by the anchor base
    header, records = [], []
    with open_text(a.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
            elif line.strip():
                records.append(line.rstrip("\n").split("\t"))
    csq = csq_parser(header)
    seen = set()
    for f in records:
        if len(f) < 5:
            continue
        chrom, pos, ref, alts = f[0], int(f[1]), f[3].upper(), f[4].upper()
        filt = f[6] if len(f) > 6 else "."
        info = f[7] if len(f) > 7 else ""
        if filt not in (".", "PASS", "") and "PASS" not in filt.split(";"):
            continue
        for alt in alts.split(","):
            key = (chrom, pos, ref, alt)
            if key in seen:
                continue
            seen.add(key)
            if alt in ("*", ".", "") or alt.startswith("<") or not re.fullmatch(r"[ACGTN]+", ref + alt):
                continue
            is_indel = len(ref) != len(alt)
            pos0 = pos - 1
            hits = [k for k in by_chrom.get(chrom, []) if contains(covered[k]["targeted"], pos0)]
            if not hits:
                continue
            terms, symbols = set(), set()
            if csq:
                sym_i, con_i, pick_i = csq
                entries = [e.split("|") for e in info_get(info, "CSQ").split(",") if e]
                use = [e for e in entries if pick_i is not None and len(e) > pick_i and e[pick_i] == "1"] or entries
                for e in use:
                    if len(e) > con_i and e[con_i]:
                        terms.update(e[con_i].split("&"))
                    if sym_i is not None and len(e) > sym_i and e[sym_i]:
                        symbols.add(e[sym_i])
            for k in hits:
                g = covered[k]
                g["n_indel" if is_indel else "n_snv"] += 1
                g["symbols"].update(symbols)
                for t in terms:
                    if g["worst"] is None or SEV_RANK.get(t, len(SEVERITY)) < SEV_RANK.get(g["worst"], len(SEVERITY)):
                        g["worst"] = t

    cnvkit = read_cnvkit(a.cnvkit_cns) if a.cnvkit_cns else None
    purecn = read_purecn(a.purecn_loh) if a.purecn_loh else None

    cols = ["sample", "mode", "source", "gene", "chrom", "start", "end", "targeted_bp", "callable_bp",
            "callable_frac", "n_snv", "n_indel", "worst_consequence", "vep_symbol", "cnvkit_cn",
            "purecn_cn", "purecn_minor_cn", "cn_source", "status"]
    order = sorted(covered.values(), key=lambda g: (g["chrom"], g["start"], g["gene"]))
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for g in order:
            frac = g["callable_bp"] / g["targeted_bp"]
            ck = pick_overlap(cnvkit.get(g["chrom"], []), g["targeted"]) if cnvkit else None
            pc = pick_overlap(purecn.get(g["chrom"], []), g["targeted"]) if purecn else None
            p_c, p_m = pc if pc else (None, None)
            cn_source = "purecn" if p_c is not None else ("cnvkit" if ck is not None else "none")
            cn = p_c if p_c is not None else ck
            if frac < a.min_callable:
                status = "NOT_CALLABLE"
            else:
                parts = []
                if g["n_snv"] + g["n_indel"] > 0:
                    parts.append("MUTATED")
                lab = cn_label(cn)
                if lab:
                    parts.append(lab)
                if p_m == 0 and p_c is not None and p_c > 0:
                    parts.append("LOH")
                status = "+".join(parts) if parts else "WT"
            row = [a.sample, a.mode, a.source, g["gene"], g["chrom"], g["start"], g["end"],
                   g["targeted_bp"], g["callable_bp"], f"{frac:.4f}", g["n_snv"], g["n_indel"],
                   g["worst"] or ".", ",".join(sorted(g["symbols"])) or ".",
                   ck if ck is not None else ".", p_c if p_c is not None else ".",
                   p_m if p_m is not None else ".", cn_source, status]
            fh.write("\t".join(str(v) for v in row) + "\n")


if __name__ == "__main__":
    main()
