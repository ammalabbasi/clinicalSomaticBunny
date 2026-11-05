#!/usr/bin/env python3
"""Tumour mutational burden (mutations per callable Mb) for one tumour sample.

Callable territory = merged --regions BED  ∩  tumour depth >= --min_depth  [∩ normal depth >= --min_normal_depth]
(depth from `mosdepth --quantize`, see Modules/TMB/TMB.nf). A variant counts when
  - it is a record of the consensus VCF (the caller decides PASS; any non-PASS FILTER value is dropped here too),
  - it is an SNV, MNV or indel (one count per ALT allele; duplicate CHROM/POS/REF/ALT counted once),
  - its anchor base (VCF POS) lies in callable territory,
  - its tumour VAF, recomputed from the tumour BAM (MAPQ/BQ filtered), is >= --min_vaf,
  - optionally: it is coding (CSQ/ANN annotation), and PureCN does not call it germline.

Writes <out>.tsv (one row), <out>.variants.tsv (every record and its decision) and <out>.callable.bed.
"""
import argparse
import csv
import gzip
import re
import sys

import pysam

CODING_TERMS = {
    "missense_variant", "synonymous_variant", "stop_gained", "stop_lost", "start_lost", "stop_retained_variant",
    "frameshift_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant",
    "coding_sequence_variant", "incomplete_terminal_codon_variant", "initiator_codon_variant",
    "start_retained_variant", "disruptive_inframe_insertion", "disruptive_inframe_deletion",
    "conservative_inframe_insertion", "conservative_inframe_deletion",
}


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_bed(path, keep_label=None):
    """Intervals {chrom: [(start, end), ...]} (0-based half-open), merged. keep_label filters column 4."""
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


def intersect(a, b):
    out = {}
    for c in set(a) & set(b):
        x, y, i, j, res = a[c], b[c], 0, 0, []
        while i < len(x) and j < len(y):
            s, e = max(x[i][0], y[j][0]), min(x[i][1], y[j][1])
            if s < e:
                res.append((s, e))
            if x[i][1] < y[j][1]:
                i += 1
            else:
                j += 1
        if res:
            out[c] = res
    return out


def total_bp(iv):
    return sum(e - s for v in iv.values() for s, e in v)


def contains(iv, chrom, pos0):
    v = iv.get(chrom)
    if not v:
        return False
    lo, hi = 0, len(v)
    while lo < hi:
        mid = (lo + hi) // 2
        if v[mid][1] <= pos0:
            lo = mid + 1
        else:
            hi = mid
    return lo < len(v) and v[lo][0] <= pos0 < v[lo][1]


def read_ok(r, min_mapq):
    return not (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail
                or r.mapping_quality < min_mapq)


def snv_vaf(bam, chrom, pos0, ref, alt, min_mapq, min_bq):
    """VAF of an SNV/MNV: alt reads / (ref + alt reads) over the first differing base."""
    k = next((i for i in range(len(ref)) if ref[i] != alt[i]), 0)
    p, rb, ab = pos0 + k, ref[k].upper(), alt[k].upper()
    n_ref = n_alt = 0
    for col in bam.pileup(chrom, p, p + 1, truncate=True, stepper="nofilter", ignore_overlaps=False,
                          min_base_quality=0, max_depth=10_000_000):
        for pr in col.pileups:
            r = pr.alignment
            if pr.is_del or pr.is_refskip or pr.query_position is None or not read_ok(r, min_mapq):
                continue
            if r.query_qualities is not None and r.query_qualities[pr.query_position] < min_bq:
                continue
            b = r.query_sequence[pr.query_position].upper()
            n_ref += b == rb
            n_alt += b == ab
    tot = n_ref + n_alt
    return n_alt, tot, (n_alt / tot if tot else 0.0)


def indel_vaf(bam, chrom, pos0, ref, alt, min_mapq, min_bq):
    """VAF of a left-anchored indel by walking each read's alignment at the anchor (as INDEL_VAF.py)."""
    is_del = len(ref) > len(alt)
    dl = len(ref) - len(alt)
    ins = alt[1:] if not is_del else ""
    n_ref = n_alt = 0
    for r in bam.fetch(chrom, max(0, pos0 - 1), pos0 + max(len(ref), len(alt)) + 1):
        if not read_ok(r, min_mapq) or r.query_sequence is None:
            continue
        pairs = r.get_aligned_pairs(matches_only=False)
        idx = next((i for i, (q, rp) in enumerate(pairs) if rp == pos0), None)
        if idx is None or pairs[idx][0] is None:
            continue
        if r.query_qualities is not None and r.query_qualities[pairs[idx][0]] < min_bq:
            continue
        if is_del:
            nxt = pairs[idx + 1: idx + 1 + dl]
            if len(nxt) == dl and all(q is None for q, _ in nxt):
                n_alt += 1
            elif len(pairs[idx + 1: idx + 2 + dl]) == dl + 1:
                n_ref += 1
        else:
            nxt = pairs[idx + 1: idx + 1 + len(ins)]
            if len(nxt) == len(ins) and all(rp is None for _, rp in nxt):
                qs = [q for q, _ in nxt if q is not None]
                if len(qs) == len(ins) and r.query_sequence[qs[0]: qs[-1] + 1].upper() == ins.upper():
                    n_alt += 1
            elif len(pairs[idx + 1: idx + 2]) > 0:
                n_ref += 1
    tot = n_ref + n_alt
    return n_alt, tot, (n_alt / tot if tot else 0.0)


def annotation_parser(header_lines):
    """Return f(info_str) -> set of consequence terms, or None when the VCF carries no CSQ/ANN annotation."""
    for h in header_lines:
        if h.startswith("##INFO=<ID=CSQ,"):
            m = re.search(r"Format: ([^\">]+)", h)
            if m:
                fields = m.group(1).strip().split("|")
                if "Consequence" in fields:
                    k = fields.index("Consequence")
                    return "CSQ", lambda info: {t for e in info_get(info, "CSQ").split(",") if e
                                                for t in (e.split("|")[k] if len(e.split("|")) > k else "").split("&")}
    for h in header_lines:
        if h.startswith("##INFO=<ID=ANN,"):
            return "ANN", lambda info: {t for e in info_get(info, "ANN").split(",") if e
                                        for t in (e.split("|")[1] if len(e.split("|")) > 1 else "").split("&")}
    return None, None


def info_get(info, key):
    """Value of INFO key up to the next ';'. The consensus INFO joins callers' INFO with '&' or '|', so the key may
    also follow one of those; the value itself keeps its '|' and '&' (CSQ/ANN field and term separators)."""
    m = re.search(r"(?:^|[;&|])" + re.escape(key) + r"=([^;]*)", info)
    return m.group(1) if m else ""


def read_germline(path, min_post):
    """PureCN *_variants.csv -> {(chr, pos, ref, alt): posterior_somatic}."""
    post = {}
    with open(path) as fh:
        rows = csv.DictReader(fh)
        cols = {c.lower(): c for c in rows.fieldnames or []}
        need = ["chr", "start", "ref", "alt", "posterior.somatic"]
        if any(n not in cols for n in need):
            sys.exit(f"--germline_calls {path}: needs columns chr,start,REF,ALT,POSTERIOR.SOMATIC (PureCN *_variants.csv)")
        for r in rows:
            try:
                post[(r[cols["chr"]], int(float(r[cols["start"]])), r[cols["ref"]], r[cols["alt"]])] = \
                    float(r[cols["posterior.somatic"]])
            except ValueError:
                continue
    return post


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--mode", required=True, choices=["paired", "tumor_only"])
    ap.add_argument("--source", required=True, help="label of the consensus the VCF came from")
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--tumor_bam", required=True)
    ap.add_argument("--tumor_bai", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--tumor_quantized", required=True)
    ap.add_argument("--normal_quantized")
    ap.add_argument("--min_depth", type=int, required=True)
    ap.add_argument("--min_normal_depth", type=int, required=True)
    ap.add_argument("--min_vaf", type=float, required=True)
    ap.add_argument("--min_mapq", type=int, default=20)
    ap.add_argument("--min_baseq", type=int, default=20)
    ap.add_argument("--coding_only", action="store_true")
    ap.add_argument("--germline_calls")
    ap.add_argument("--min_somatic_posterior", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    regions = read_bed(a.regions)
    callable_iv = intersect(regions, read_bed(a.tumor_quantized, keep_label="CALLABLE"))
    if a.mode == "paired":
        if not a.normal_quantized:
            sys.exit("paired mode needs --normal_quantized")
        callable_iv = intersect(callable_iv, read_bed(a.normal_quantized, keep_label="CALLABLE"))
    callable_bp = total_bp(callable_iv)
    with open(a.out + ".callable.bed", "w") as fh:
        for c in sorted(callable_iv):
            for s, e in callable_iv[c]:
                fh.write(f"{c}\t{s}\t{e}\n")

    header, records = [], []
    with open_text(a.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
            elif line.strip():
                records.append(line.rstrip("\n").split("\t"))

    notes = []
    coding_status = "off"
    ann_kind, ann_terms = annotation_parser(header)
    if a.coding_only:
        if ann_kind:
            coding_status = f"on ({ann_kind})"
        else:
            coding_status = "ignored (no CSQ/ANN annotation)"
            print(f"WARNING [{a.sample}]: --tmb_coding_only ignored: the VCF has no CSQ/ANN annotation; "
                  "counting all variants in callable territory", file=sys.stderr)
            notes.append("coding_only_ignored_no_annotation")

    germ = None
    germ_status = "off"
    if a.germline_calls:
        if a.mode == "tumor_only":
            germ = read_germline(a.germline_calls, a.min_somatic_posterior)
            germ_status = f"on (POSTERIOR.SOMATIC >= {a.min_somatic_posterior})"
        else:
            germ_status = "not used (paired)"

    bam = pysam.AlignmentFile(a.tumor_bam, "rb", index_filename=a.tumor_bai)
    seen = set()
    n_snv = n_indel = 0
    counts = {}
    with open(a.out + ".variants.tsv", "w") as vt:
        vt.write("chrom\tpos\tref\talt\ttype\ttumor_alt_reads\ttumor_depth\ttumor_vaf\tcallable\tdecision\n")
        for f in records:
            if len(f) < 5:
                continue
            chrom, pos, ref, alts = f[0], int(f[1]), f[3].upper(), f[4].upper()
            filt = f[6] if len(f) > 6 else "."
            info = f[7] if len(f) > 7 else ""
            for alt in alts.split(","):
                key = (chrom, pos, ref, alt)
                if key in seen:
                    continue
                seen.add(key)
                if alt in ("*", ".", "") or alt.startswith("<") or not re.fullmatch(r"[ACGTN]+", ref + alt):
                    vtype, decision = "other", "not_snv_or_indel"
                elif len(ref) == len(alt):
                    vtype = "snv" if len(ref) == 1 else "mnv"
                    decision = None
                elif ref[0] == alt[0]:
                    vtype, decision = "indel", None
                else:
                    vtype, decision = "complex", "not_snv_or_indel"
                pos0 = pos - 1
                inside = contains(callable_iv, chrom, pos0)
                alt_n = depth = 0
                vaf = 0.0
                if decision is None and filt not in (".", "PASS", "") and "PASS" not in filt.split(";"):
                    decision = "not_pass"
                if decision is None and not inside:
                    decision = "outside_callable"
                if decision is None:
                    fn = indel_vaf if vtype == "indel" else snv_vaf
                    alt_n, depth, vaf = fn(bam, chrom, pos0, ref, alt, a.min_mapq, a.min_baseq)
                    if vaf < a.min_vaf:
                        decision = "low_vaf"
                if decision is None and coding_status.startswith("on"):
                    if not (ann_terms(info) & CODING_TERMS):
                        decision = "non_coding"
                if decision is None and germ is not None:
                    p = germ.get(key)
                    if p is not None and p < a.min_somatic_posterior:
                        decision = "germline_purecn"
                if decision is None:
                    decision = "counted"
                    if vtype == "indel":
                        n_indel += 1
                    else:
                        n_snv += 1
                counts[decision] = counts.get(decision, 0) + 1
                vt.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{vtype}\t{alt_n}\t{depth}\t{vaf:.4f}\t"
                         f"{int(inside)}\t{decision}\n")
    bam.close()

    n = n_snv + n_indel
    mb = callable_bp / 1e6
    tmb = n / mb if mb > 0 else float("nan")
    warning = ""
    if mb < 1:
        warning = f"callable_Mb<1 ({mb:.4f} Mb): TMB is unreliable on this little territory"
        print(f"WARNING [{a.sample}]: {warning}", file=sys.stderr)
    if notes:
        warning = ";".join(([warning] if warning else []) + notes)
    excluded = ",".join(f"{k}={v}" for k, v in sorted(counts.items()) if k != "counted") or "none"
    cols = [
        ("sample", a.sample), ("mode", a.mode), ("source", a.source),
        ("n_variants", n), ("n_snv", n_snv), ("n_indel", n_indel),
        ("callable_bp", callable_bp), ("callable_Mb", f"{mb:.6f}"),
        ("TMB", f"{tmb:.4f}" if mb > 0 else "NA"),
        ("regions_bp", total_bp(regions)),
        ("min_depth", a.min_depth), ("min_normal_depth", a.min_normal_depth if a.mode == "paired" else "NA"),
        ("min_vaf", a.min_vaf), ("min_mapq", a.min_mapq), ("min_baseq", a.min_baseq),
        ("coding_only", coding_status), ("germline_filter", germ_status),
        ("excluded", excluded), ("warning", warning or "none"),
    ]
    with open(a.out + ".tsv", "w") as fh:
        fh.write("\t".join(k for k, _ in cols) + "\n")
        fh.write("\t".join(str(v) for _, v in cols) + "\n")


if __name__ == "__main__":
    main()
