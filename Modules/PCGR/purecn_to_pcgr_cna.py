#!/usr/bin/env python3
"""Convert PureCN allele-specific segments (<sample>_loh.csv) to the PCGR --input_cna TSV.

PureCN loh.csv columns used: chr, start, end, C (total copy number), M (minor copy number, NA when the segment has
no informative heterozygous SNPs). PCGR wants Chromosome, Start, End (1-based), nMajor, nMinor as integers.

  nMinor = M, nMajor = C - M
  M is NA -> nMinor = floor(C/2), nMajor = C - nMinor: total copy number (what PCGR uses for gain/amplification/
             deletion) is kept, and no LOH is claimed for a segment where PureCN could not measure it.
Segments with C = NA are dropped. Prints one summary line of counts.

Usage: purecn_to_pcgr_cna.py <loh.csv> <out.tsv>
"""
import csv
import sys


def num(x):
    try:
        return int(round(float(x)))
    except (TypeError, ValueError):
        return None


def main(src, dst):
    kept = dropped = minor_imputed = 0
    with open(src, newline='') as fh, open(dst, 'w') as out:
        rows = csv.DictReader(fh)
        missing = {'chr', 'start', 'end', 'C', 'M'} - set(rows.fieldnames or [])
        if missing:
            sys.exit(f'ERROR: {src} is not a PureCN loh.csv (missing columns: {",".join(sorted(missing))})')
        out.write('Chromosome\tStart\tEnd\tnMajor\tnMinor\n')
        for r in rows:
            c, m, s, e = num(r['C']), num(r['M']), num(r['start']), num(r['end'])
            if c is None or s is None or e is None:
                dropped += 1
                continue
            if m is None or m > c or m < 0:
                m = c // 2
                minor_imputed += 1
            out.write(f"{r['chr']}\t{max(s, 1)}\t{e}\t{c - m}\t{m}\n")
            kept += 1
    print(f'purecn_to_pcgr_cna segments={kept} dropped_no_C={dropped} minor_cn_imputed={minor_imputed}')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
