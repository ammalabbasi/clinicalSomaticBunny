#!/usr/bin/env python3
"""Write a clean single-sample VCF for PCGR from a SomaticBunny consensus VCF.

PCGR reads tumour (and control) depth and allele fraction only from INFO tags. The SomaticBunny consensus VCFs do
not carry them in a usable form:

  * paired: POST/standard_vcf/<id>/<id>_scenario_N_standard_PASSed.vcf has 8 columns (no FORMAT/sample columns);
    its INFO is the callers' INFO strings joined with '&' or '|', so it is not parseable as VCF INFO.
  * tumour-only: POST/tumor_only/<id>/<id>_tumor_only_consensus.PASS.vcf.gz keeps the Mutect2 record, with
    FORMAT AD/DP for the tumour sample.

So only CHROM/POS/REF/ALT are taken from the consensus, and per-sample read counts come from FORMAT fields:
  paired:      Mutect2 (FORMAT AD, DP) -> SAGE (AD, DP) -> Strelka2 (SNV: tier-1 <ALT>U/<REF>U; indel: TIR/TAR)
  tumour-only: the consensus record itself (Mutect2 FORMAT AD, DP)
The first caller that has the variant wins (INFO/TVAF_SRC records which). For a 2-of-4 SNV / 2-of-3 indel consensus
every site is in at least one of these callers; anything left without counts is reported on stderr and written
without TDP/TAF (PCGR then treats depth/AF as missing).

Output INFO tags: TDP (Integer), TAF (Float), NDP, NAF (paired only), TVAF_SRC, CALLERS. FILTER is PASS.
Only primary chromosomes (1-22, X, Y, with or without chr) are kept; PCGR does not annotate decoys/alts.

Usage:
  pcgr_prepare_vcf.py --mode paired --sites X.vcf --fai ref.fai --out out.vcf \
      --mutect2 m.vcf --sage s.vcf.gz --strelka_snv s.vcf --strelka_indel i.vcf.gz
  pcgr_prepare_vcf.py --mode tumor_only --sites X.vcf.gz --fai ref.fai --out out.vcf
Prints one summary line of counts (no variant records) to stdout.
"""
import argparse
import gzip
import re
import sys

PRIMARY = re.compile(r'^(chr)?([0-9]+|X|Y)$')


def vopen(path):
    with open(path, 'rb') as fh:
        magic = fh.read(2)
    return gzip.open(path, 'rt') if magic == b'\x1f\x8b' else open(path, 'rt')


def split_records(ref, alt):
    """Yield (ref, alt) biallelic records; MNPs are split into SNVs like the POST scenario scripts do."""
    if len(ref) == len(alt) and len(ref) > 1:
        for i in range(len(ref)):
            if ref[i] != alt[i]:
                yield i, ref[i], alt[i]
    else:
        yield 0, ref, alt


def to_int(x):
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


class CallerCounts:
    """Per-site tumour/normal (depth, alt count) from one caller VCF."""

    def __init__(self, path, kind):
        self.kind = kind
        self.sites = {}
        if path and path != 'NONE':
            self._load(path)

    def _columns(self, header_cols, meta):
        samples = header_cols[9:]
        if not samples:
            return None, None
        if self.kind == 'mutect2':
            t = meta.get('tumor_sample')
            n = meta.get('normal_sample')
            ti = samples.index(t) if t in samples else None
            ni = samples.index(n) if n in samples else None
            if ti is None:
                ti = 0 if len(samples) == 1 else None
            return ti, ni
        if self.kind == 'sage':
            ti = next((i for i, s in enumerate(samples) if s.endswith('_tumor')), None)
            ni = next((i for i, s in enumerate(samples) if s.endswith('_normal')), None)
            if ti is None and len(samples) == 1:
                ti = 0
            return ti, ni
        if self.kind in ('strelka_snv', 'strelka_indel'):
            ti = samples.index('TUMOR') if 'TUMOR' in samples else None
            ni = samples.index('NORMAL') if 'NORMAL' in samples else None
            return ti, ni
        return None, None

    def _load(self, path):
        meta = {}
        ti = ni = None
        with vopen(path) as fh:
            for line in fh:
                if line.startswith('##'):
                    m = re.match(r'^##(tumor_sample|normal_sample)=(.+)$', line.rstrip('\n'))
                    if m:
                        meta[m.group(1)] = m.group(2)
                    continue
                if line.startswith('#CHROM'):
                    ti, ni = self._columns(line.rstrip('\n').split('\t'), meta)
                    if ti is None:
                        sys.stderr.write(f'WARNING: {self.kind}: tumour column not found in {path}; not used\n')
                        return
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) < 10:
                    continue
                keys = f[8].split(':')
                tum = dict(zip(keys, f[9 + ti].split(':')))
                nor = dict(zip(keys, f[9 + ni].split(':'))) if ni is not None and len(f) > 9 + ni else None
                for k, alt in enumerate(f[4].split(',')):
                    t = self._counts(tum, f[3], alt, k)
                    n = self._counts(nor, f[3], alt, k) if nor is not None else None
                    for off, r1, a1 in split_records(f[3], alt):
                        self.sites.setdefault((f[0], int(f[1]) + off, r1, a1), (t, n))

    def _counts(self, s, ref, alt, k):
        """(depth, alt_reads, informative_reads) or None."""
        if self.kind in ('mutect2', 'sage'):
            ad = [to_int(x) for x in s.get('AD', '').split(',')]
            if len(ad) < k + 2 or any(x is None for x in ad):
                return None
            dp = to_int(s.get('DP'))
            depth = sum(ad) if dp is None or dp < sum(ad) else dp
            return depth, ad[k + 1], sum(ad)
        if self.kind == 'strelka_snv':
            if len(ref) != 1 or len(alt) != 1:
                return None
            r = to_int(s.get(ref + 'U', '').split(',')[0])
            a = to_int(s.get(alt + 'U', '').split(',')[0])
            dp = to_int(s.get('DP'))
            if r is None or a is None:
                return None
            return (dp if dp is not None and dp >= r + a else r + a), a, r + a
        if self.kind == 'strelka_indel':
            r = to_int(s.get('TAR', '').split(',')[0])
            a = to_int(s.get('TIR', '').split(',')[0])
            dp = to_int(s.get('DP'))
            if r is None or a is None:
                return None
            return (dp if dp is not None and dp >= r + a else r + a), a, r + a
        return None

    def get(self, key):
        return self.sites.get(key)


def callers_of(fields):
    """Caller list from the consensus line: INFO 'callers=' (SNVs) or the QUAL column (indels)."""
    m = re.search(r'callers=([A-Za-z,]+)', '\t'.join(fields[7:]))
    if m:
        return m.group(1)
    if len(fields) > 5 and re.fullmatch(r'[a-z]{2}(,[a-z]{2})*', fields[5]):
        return fields[5]
    return None


def af(counts):
    _, a, informative = counts
    return f'{a / informative:.4f}' if informative else '0'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mode', choices=['paired', 'tumor_only'], required=True)
    p.add_argument('--sites', required=True)
    p.add_argument('--fai', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--mutect2')
    p.add_argument('--sage')
    p.add_argument('--strelka_snv')
    p.add_argument('--strelka_indel')
    a = p.parse_args()

    order = {}
    with open(a.fai) as fh:
        for i, line in enumerate(fh):
            c = line.split('\t')
            if PRIMARY.match(c[0]):
                order[c[0]] = (i, c[1])

    if a.mode == 'paired':
        sources = [('mutect2', CallerCounts(a.mutect2, 'mutect2')),
                   ('sage', CallerCounts(a.sage, 'sage')),
                   ('strelka', CallerCounts(a.strelka_snv, 'strelka_snv')),
                   ('strelka', CallerCounts(a.strelka_indel, 'strelka_indel'))]
    else:
        sources = [('consensus', CallerCounts(a.sites, 'mutect2'))]

    records = {}
    n_in = n_nonprimary = n_filtered = n_nocounts = 0
    src_n = {}
    with vopen(a.sites) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 5:
                continue
            n_in += 1
            if a.mode == 'tumor_only' and len(f) > 6 and f[6] not in ('PASS', '.'):
                n_filtered += 1
                continue
            if f[0] not in order:
                n_nonprimary += 1
                continue
            callers = callers_of(f)
            for alt in f[4].split(','):
                for off, r1, a1 in split_records(f[3], alt):
                    key = (f[0], int(f[1]) + off, r1, a1)
                    if key in records:
                        continue
                    info = []
                    for name, src in sources:
                        hit = src.get(key)
                        if hit and hit[0]:
                            t, n = hit
                            info += [f'TDP={t[0]}', f'TAF={af(t)}']
                            if a.mode == 'paired' and n:
                                info += [f'NDP={n[0]}', f'NAF={af(n)}']
                            info.append(f'TVAF_SRC={name}')
                            src_n[name] = src_n.get(name, 0) + 1
                            break
                    else:
                        n_nocounts += 1
                    if callers:
                        info.append(f'CALLERS={callers}')
                    records[key] = ';'.join(info) if info else '.'

    with open(a.out, 'w') as out:
        out.write('##fileformat=VCFv4.2\n')
        out.write('##source=SomaticBunny pcgr_prepare_vcf.py (' + a.mode + ')\n')
        for c, (i, ln) in sorted(order.items(), key=lambda x: x[1][0]):
            out.write(f'##contig=<ID={c},length={ln}>\n')
        out.write('##INFO=<ID=TDP,Number=1,Type=Integer,Description="Read depth at the site in the tumour (FORMAT DP/AD of TVAF_SRC)">\n')
        out.write('##INFO=<ID=TAF,Number=1,Type=Float,Description="Alt allele fraction in the tumour (alt reads / informative reads)">\n')
        if a.mode == 'paired':
            out.write('##INFO=<ID=NDP,Number=1,Type=Integer,Description="Read depth at the site in the matched normal">\n')
            out.write('##INFO=<ID=NAF,Number=1,Type=Float,Description="Alt allele fraction in the matched normal">\n')
        out.write('##INFO=<ID=TVAF_SRC,Number=1,Type=String,Description="Caller whose FORMAT fields gave TDP/TAF">\n')
        out.write('##INFO=<ID=CALLERS,Number=.,Type=String,Description="SomaticBunny consensus callers (mt=Mutect2, ms=MuSE, st=Strelka2, sa=SAGE)">\n')
        out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for (c, pos, r, alt), info in sorted(records.items(), key=lambda x: (order[x[0][0]][0], x[0][1], x[0][2], x[0][3])):
            out.write(f'{c}\t{pos}\t.\t{r}\t{alt}\t.\tPASS\t{info}\n')

    srcs = ','.join(f'{k}:{v}' for k, v in sorted(src_n.items())) or 'none'
    print(f'pcgr_prepare_vcf mode={a.mode} input_records={n_in} written={len(records)} '
          f'non_primary_dropped={n_nonprimary} non_pass_dropped={n_filtered} '
          f'without_counts={n_nocounts} counts_from={srcs}')


if __name__ == '__main__':
    main()
